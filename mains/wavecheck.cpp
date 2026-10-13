// Waveform check on a SOLARIS run that saved the ADC input trace for every
// event (SOL "OneTrace" blocks: analog probe 0 only, no firmware trapezoid
// and no digital trigger probe). Every mapped channel is processed with the
// filter timing it was acquired with, read from the DPP-PHA settings file next
// to the data, and the trapezoid is recomputed offline from the raw trace with
// the same rise, flat top and pole-zero (TrapezoidFilter, the Jordanov-Knoll
// recursion the firmware runs). A second, alternative trapezoid (--altrise,
// --altflat, same pole-zero) is run on the same traces for comparison.
//
// Time zero of every record is the DAQ trigger, ChPreTriggerT into the
// record. The pulse itself is located from the data (10 % of its leading
// edge), because the time filter fires some way up a slow edge and the two
// need not coincide. Measured per channel:
//   * the input pulse: where it starts relative to the trigger, polarity,
//     height, 10-90 % rise time, and the preamp decay time from the tail
//     after isolated pulses;
//   * the trapezoid: where its flat top is, whether the firmware's nominal
//     read point (trigger + rise + peaking position) sits on it, which
//     sample the reported energy actually correlates best with, and the tail
//     after the trapezoid against the pole-zero prediction of
//     doc/polezero.tex, Eq. (c);
//   * the beam-peak resolution of the reported energy against offline
//     estimators: the set trapezoid at the sample the firmware reads, its
//     maximum, and the alternative trapezoid's maximum and flat-top centre.
// Two passes over the files: the first collects per-hit features, the
// correlation scan and the times to the previous and next hit on each
// channel; the second fills the spectra of the time-isolated hits, averages
// the traces of the beam-like ones and keeps a few single events per group.
//
//   wavecheck --settings=<XSetting.dat> --out=<dir> [--nfiles=N]
//             [--maxtraces=N] [--altrise=ns] [--altflat=ns] [--minenergy=ADC]
//             <sol files...>
#include "BinaryUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "SolWave.hpp"
#include "TrapezoidFilter.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace {

using namespace solwave;

const Double_t kIsolatedBefore_us = 60.0;
const Int_t kExamplesPerGroup = 3;
// Candidate read samples for the energy correlation scan: every 4 samples
// (32 ns) from this far before the trigger to this far after it.
const Int_t kScanFrom = -500, kScanTo = 1500, kScanStep = 4;
const Int_t kScanN = (kScanTo - kScanFrom) / kScanStep;
// Energy estimators compared on the time-isolated hits.
enum Estimator {
  kFirmware = 0, // energy reported by the DPP-PHA
  kSetRead = 1,  // set trapezoid at the sample the firmware correlates with
  kSetMax = 2,   // set trapezoid maximum
  kAltMax = 3,   // alternative trapezoid maximum
  kAltRead = 4,  // alternative trapezoid at pulse start + rise + flat / 2
  kNEst = 5
};
const char *kEstName[kNEst] = {"firmware", "set@read", "set max", "alt max",
                               "alt@flat"};

struct Hit {
  Int_t channel, group;
  Double_t t_ns, energy;
  Int_t iStart;          // sample where the input crosses 10 % of its height
  Double_t inPre, inAmp; // input baseline and signed pulse height, ADC
  Double_t noise;        // RMS of the baseline before the pulse, ADC
  Double_t rise_ns;      // 10-90 % rise time of the input
  Double_t trapPre;      // set trapezoid just before the pulse
  Double_t trapNominal;  // set trapezoid at the nominal firmware read window
  Double_t trapMax;      // its maximum
  Double_t tTrapMax_us;  // when that maximum occurs, after the trigger
  Double_t altMax;       // alternative trapezoid maximum
  Double_t altRead;      // alternative trapezoid at start + rise + flat / 2
  Double_t dtPrev_us, dtNext_us;
  Bool_t isolated, beamlike;
};

Bool_t HitBefore(const Hit &a, const Hit &b) { return a.t_ns < b.t_ns; }

// The trapezoid a filter produces from this record, with the baseline taken
// out and the pulse made positive first.
void OfflineTrapezoid(const Int_t *in, Int_t len, Double_t pre, Double_t pol,
                      TrapezoidFilter &filter, std::vector<Double_t> &trap) {
  filter.Reset();
  trap.resize(len);
  for (Int_t i = 0; i < len; i++)
    trap[i] = filter.Process(pol * (Double_t(in[i]) - pre));
}

Double_t MaxFrom(const std::vector<Double_t> &v, Int_t from, Double_t base,
                 Int_t &iMax) {
  Double_t mx = -1.0e30;
  iMax = from;
  for (Int_t i = from; i < Int_t(v.size()); i++)
    if (v[i] - base > mx) {
      mx = v[i] - base;
      iMax = i;
    }
  return mx;
}

// Sample index the pulse start is moved to in the start-aligned average.
const Int_t kAlignRef = 250;

struct ShapeSum {
  std::vector<Double_t> input, trap, alt, aligned, alignedN;
  Long64_t n = 0;
  void Init(Int_t len) {
    input.assign(len, 0.0);
    trap.assign(len, 0.0);
    alt.assign(len, 0.0);
    aligned.assign(len, 0.0);
    alignedN.assign(len, 0.0);
  }
  // The aligned sum divided by how many traces reached each sample.
  std::vector<Double_t> AlignedMean() const {
    std::vector<Double_t> m(aligned.size(), 0.0);
    for (size_t i = 0; i < aligned.size(); i++)
      m[i] = alignedN[i] > 0.0 ? aligned[i] / alignedN[i] : 0.0;
    return m;
  }
};

struct Example {
  Int_t channel;
  Double_t energy, inPre, inAmp, trapPre;
  std::vector<Double_t> input, trap, alt;
};

// Running sums for the correlation of the reported energy with the set
// trapezoid at each candidate sample.
struct ScanSum {
  std::vector<Double_t> sx, sxx, sxy;
  Double_t sy = 0.0, syy = 0.0;
  Long64_t n = 0;
  void Init() {
    sx.assign(kScanN, 0.0);
    sxx.assign(kScanN, 0.0);
    sxy.assign(kScanN, 0.0);
  }
};

TGraph *TraceGraph(const std::vector<Double_t> &v, Int_t pre, Double_t scale,
                   const char *name) {
  TGraph *g = new TGraph(Int_t(v.size()));
  for (Int_t i = 0; i < Int_t(v.size()); i++)
    g->SetPoint(i, (Double_t(i) - Double_t(pre)) * kSample_ns / 1000.0,
                v[i] * scale);
  g->SetName(name);
  return g;
}

// Tail amplitude c of doc/polezero.tex Eq. (c), continuum form, for a
// trapezoid of rise k and rise + flat l on a preamp tau_pre with the
// pole-zero at tau_pz. Positive c is an undershoot.
Double_t TailAmplitude(Double_t tauPre, Double_t tauPz, Double_t kUs,
                       Double_t lUs) {
  if (tauPre <= 0.0 || tauPz <= 0.0 || kUs <= 0.0)
    return 0.0;
  return (tauPz - tauPre) * tauPre / (kUs * tauPz) *
         (1.0 - TMath::Exp(kUs / tauPre)) * (1.0 - TMath::Exp(lUs / tauPre));
}

struct PeakFit {
  Double_t mean = 0.0, sigma = 0.0;
  Bool_t ok = kFALSE;
};

// Gaussian on the main peak above a floor: locate the mode, fit +-12 %.
PeakFit FitPeak(TH1D *h, Double_t floor) {
  PeakFit p;
  if (!h || h->GetEntries() < 50)
    return p;
  const Int_t bLo = h->FindBin(floor);
  h->GetXaxis()->SetRange(bLo, h->GetNbinsX());
  const Double_t mode = h->GetBinCenter(h->GetMaximumBin());
  h->GetXaxis()->SetRange(1, h->GetNbinsX());
  if (mode <= floor)
    return p;
  TF1 *f = new TF1(PlottingUtils::GetRandomName(), "gaus", 0.88 * mode,
                   1.12 * mode);
  f->SetParameters(h->GetMaximum(), mode, 0.03 * mode);
  h->Fit(f, "RQN0");
  p.mean = f->GetParameter(1);
  p.sigma = TMath::Abs(f->GetParameter(2));
  p.ok = p.mean > floor && p.sigma > 0.0 && p.sigma < 0.5 * p.mean;
  delete f;
  return p;
}

struct ChResult {
  Long64_t nHits = 0, nIso = 0, nBeam = 0;
  Double_t mode = 0.0, eLo = 0.0, eHi = 0.0;
  Double_t polarity = 0.0;
  Double_t startMean_us = 0.0, startRms_us = 0.0; // pulse start vs trigger
  Double_t riseMedian_ns = 0.0;
  Double_t noise = 0.0, noiseFrac = 0.0; // baseline RMS, ADC and % of height
  Double_t tauPre_us = 0.0, tauPreErr_us = 0.0;
  Double_t flatFrom_us = 0.0, flatTo_us = 0.0; // 99 % of max, vs trigger
  Double_t nominalOverMax = 0.0; // trapezoid at the nominal read / maximum
  Double_t bestRead_us = 0.0;    // sample the energy correlates best with
  Int_t iBest = 0;
  Double_t bestCorr = 0.0, nominalCorr = 0.0;
  Double_t bestOverMax = 0.0;    // trapezoid at that sample / maximum
  Double_t tailMeasured = 0.0;   // -tail amplitude / height, i.e. c measured
  Double_t tailPredicted = 0.0;  // Eq. (c) with the fitted tau_pre
  PeakFit peak[kNEst];
};

Double_t Median(std::vector<Double_t> &v) {
  if (v.empty())
    return 0.0;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}

Double_t Correlation(const ScanSum &s, Int_t k) {
  const Double_t n = Double_t(s.n);
  if (s.n < 20)
    return 0.0;
  const Double_t cov = s.sxy[k] / n - (s.sx[k] / n) * (s.sy / n);
  const Double_t vx = s.sxx[k] / n - (s.sx[k] / n) * (s.sx[k] / n);
  const Double_t vy = s.syy / n - (s.sy / n) * (s.sy / n);
  if (vx <= 0.0 || vy <= 0.0)
    return 0.0;
  return cov / TMath::Sqrt(vx * vy);
}

void DrawMark(Double_t x, Double_t yLo, Double_t yHi, Int_t colour,
              Int_t style = 2) {
  TLine *l = new TLine(x, yLo, x, yHi);
  l->SetLineStyle(style);
  l->SetLineColor(colour);
  l->SetLineWidth(2);
  l->Draw();
}

} // namespace

int main(int argc, char **argv) {
  TString out = "output/wavecheck";
  TString settingsPath;
  Long64_t maxFiles = 1000, maxTraces = -1;
  Double_t altRise_ns = 2000.0, altFlat_ns = 3000.0, minEnergy = 1000.0;
  std::vector<TString> inputs;
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a.BeginsWith("--out="))
      out = a(6, a.Length() - 6);
    else if (a.BeginsWith("--settings="))
      settingsPath = a(11, a.Length() - 11);
    else if (a.BeginsWith("--nfiles="))
      maxFiles = TString(a(9, a.Length() - 9)).Atoll();
    else if (a.BeginsWith("--maxtraces="))
      maxTraces = TString(a(12, a.Length() - 12)).Atoll();
    else if (a.BeginsWith("--altrise="))
      altRise_ns = TString(a(10, a.Length() - 10)).Atof();
    else if (a.BeginsWith("--altflat="))
      altFlat_ns = TString(a(10, a.Length() - 10)).Atof();
    else if (a.BeginsWith("--minenergy="))
      minEnergy = TString(a(12, a.Length() - 12)).Atof();
    else
      inputs.push_back(a);
  }
  if (inputs.empty() || settingsPath.Length() == 0) {
    std::cerr << "usage: wavecheck --settings=XSetting.dat --out=dir "
                 "[--nfiles=N] [--maxtraces=N] [--altrise=ns] [--altflat=ns] "
                 "[--minenergy=ADC] <sol files>"
              << std::endl;
    return 1;
  }
  if (Long64_t(inputs.size()) > maxFiles)
    inputs.resize(size_t(maxFiles));
  gSystem->mkdir(out, kTRUE);
  const TString absOut = gSystem->IsAbsoluteFileName(out)
                             ? out
                             : TString(gSystem->WorkingDirectory()) + "/" + out;
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, absOut, absOut, kFALSE);
  std::ofstream summary((absOut + "/summary.txt").Data());
  const Int_t altRise = Int_t(TMath::Nint(altRise_ns / kSample_ns));
  const Int_t altFlat = Int_t(TMath::Nint(altFlat_ns / kSample_ns));

  std::vector<ChSettings> cfg = ReadSettings(settingsPath);
  std::vector<TrapezoidFilter *> filters(kMaxCh, nullptr), altFilters(kMaxCh, nullptr);
  const TString head0 = Form("settings %s; alternative trapezoid rise %.0f ns "
                             "flat %.0f ns (same pole-zero)",
                             settingsPath.Data(), altRise_ns, altFlat_ns);
  std::cout << head0 << std::endl;
  summary << head0 << std::endl;
  for (Int_t c = 0; c < kMaxCh; c++) {
    if (GroupOf(c) < 0)
      continue;
    const ChSettings &s = cfg[c];
    const TString line =
        Form("  ch %2d  strip %2d  %-8s  %s  pre %4d  rise %4d  flat %4d  "
             "peak +%3d avg %2d  pz %5.1f us  baseline %s  probe1 %s",
             c, StripOf(c), kGroupName[GroupOf(c)],
             s.enabled ? "on " : "OFF", s.pre, s.rise, s.flat, s.peakDelay,
             s.avg, s.pz_us, s.baseline.Data(), s.probe1.Data());
    std::cout << line << std::endl;
    summary << line << std::endl;
    if (s.rise > 0 && s.pz_us > 0.0) {
      filters[c] = new TrapezoidFilter(s.rise, s.flat, kSample_ns,
                                       s.pz_us * 1000.0);
      altFilters[c] = new TrapezoidFilter(altRise, altFlat, kSample_ns,
                                          s.pz_us * 1000.0);
    }
  }

  // --- pass 1: features and the correlation scan --------------------------------
  std::vector<Hit> hits;
  std::vector<Int_t> traceLen(kMaxCh, 0);
  std::vector<ScanSum> scan(kMaxCh);
  for (Int_t c = 0; c < kMaxCh; c++)
    scan[c].Init();
  std::vector<Double_t> trap, alt;
  Long64_t nBlocks = 0, nTraces = 0, nUnmapped = 0, nNoPulse = 0;
  for (size_t f = 0; f < inputs.size(); f++) {
    SOLReader reader;
    if (!reader.Open(inputs[f].Data())) {
      std::cerr << "cannot open " << inputs[f] << std::endl;
      continue;
    }
    Long64_t nThisFile = 0;
    while (reader.ReadEvent()) {
      nBlocks++;
      const SOLData &ev = reader.GetCurrentEvent();
      if (!ev.hasTraces())
        continue;
      nTraces++;
      if (maxTraces > 0 && nThisFile >= maxTraces)
        break;
      const Int_t ch = Int_t(ev.channel);
      const Int_t group = GroupOf(ch);
      if (group < 0 || !filters[ch]) {
        nUnmapped++;
        continue;
      }
      const ChSettings &s = cfg[ch];
      const Int_t len = Int_t(ev.getSamples());
      const Int_t *in = ev.getAnalog0();
      const Int_t nominalEnd = s.pre + s.rise + s.peakDelay + s.avg;
      if (!in || len < nominalEnd + 64 || s.pre + kScanFrom < 0 ||
          s.pre + kScanTo >= len)
        continue;
      if (traceLen[ch] == 0)
        traceLen[ch] = len;
      if (len != traceLen[ch])
        continue;
      nThisFile++;
      Hit h;
      h.channel = ch;
      h.group = group;
      h.t_ns = Double_t(ev.timestamp);
      h.energy = Double_t(ev.energy);
      Double_t pol = 1.0, amp = 0.0;
      Int_t i90 = 0, iExt = 0;
      if (!LocatePulse(in, len, h.inPre, pol, amp, h.noise, h.iStart, i90,
                       iExt)) {
        nNoPulse++;
        continue;
      }
      h.inAmp = pol * amp;
      h.rise_ns = Double_t(i90 - h.iStart) * kSample_ns;
      OfflineTrapezoid(in, len, h.inPre, pol, *filters[ch], trap);
      OfflineTrapezoid(in, len, h.inPre, pol, *altFilters[ch], alt);
      h.trapPre = MeanD(trap, TMath::Max(0, h.iStart - 60), 56);
      const Int_t peakStart = s.pre + s.rise + s.peakDelay - s.avg / 2;
      h.trapNominal = MeanD(trap, peakStart, s.avg) - h.trapPre;
      Int_t iMx = 0;
      h.trapMax = MaxFrom(trap, h.iStart, h.trapPre, iMx);
      h.tTrapMax_us = Double_t(iMx - s.pre) * kSample_ns / 1000.0;
      const Double_t altPre = MeanD(alt, TMath::Max(0, h.iStart - 60), 56);
      Int_t iAlt = 0;
      h.altMax = MaxFrom(alt, h.iStart, altPre, iAlt);
      const Int_t iAltRead = h.iStart + altRise + altFlat / 2;
      h.altRead = iAltRead + 8 < len ? MeanD(alt, iAltRead - 8, 16) - altPre
                                     : 0.0;
      h.dtPrev_us = -1.0;
      h.dtNext_us = -1.0;
      h.isolated = kFALSE;
      h.beamlike = kFALSE;
      hits.push_back(h);
      if (h.energy >= minEnergy) {
        ScanSum &ss = scan[ch];
        for (Int_t k = 0; k < kScanN; k++) {
          const Double_t x = trap[s.pre + kScanFrom + k * kScanStep] - h.trapPre;
          ss.sx[k] += x;
          ss.sxx[k] += x * x;
          ss.sxy[k] += x * h.energy;
        }
        ss.sy += h.energy;
        ss.syy += h.energy * h.energy;
        ss.n++;
      }
    }
  }
  std::cout << inputs.size() << " files, " << nBlocks << " blocks, " << nTraces
            << " with traces, " << nUnmapped << " on unmapped channels, "
            << nNoPulse << " without a usable pulse, " << hits.size()
            << " kept" << std::endl;
  if (hits.empty()) {
    summary.close();
    return 1;
  }

  std::stable_sort(hits.begin(), hits.end(), HitBefore);
  std::vector<Int_t> lastIdx(kMaxCh, -1);
  std::vector<TH1D *> spec(kMaxCh, nullptr);
  for (size_t i = 0; i < hits.size(); i++) {
    Hit &h = hits[i];
    if (lastIdx[h.channel] >= 0) {
      Hit &p = hits[size_t(lastIdx[h.channel])];
      h.dtPrev_us = (h.t_ns - p.t_ns) / 1000.0;
      p.dtNext_us = h.dtPrev_us;
    }
    lastIdx[h.channel] = Int_t(i);
    if (!spec[h.channel]) {
      spec[h.channel] = new TH1D(Form("spec_ch%d", h.channel),
                                 ";Energy [ADC];Hits", 1024, 0.0, 65536.0);
      spec[h.channel]->SetDirectory(nullptr);
    }
    spec[h.channel]->Fill(h.energy);
  }
  std::vector<ChResult> res(kMaxCh);
  for (Int_t c = 0; c < kMaxCh; c++) {
    if (!spec[c])
      continue;
    spec[c]->GetXaxis()->SetRange(5, 1024); // above 256 ADC
    res[c].mode = spec[c]->GetBinCenter(spec[c]->GetMaximumBin());
    spec[c]->GetXaxis()->SetRange(1, 1024);
    res[c].nHits = Long64_t(spec[c]->GetEntries());
    // Beam-like window: around the peak when the channel has one above the
    // floor, otherwise anything above the floor (quiet channels).
    if (res[c].mode >= minEnergy) {
      res[c].eLo = 0.7 * res[c].mode;
      res[c].eHi = 1.3 * res[c].mode;
    } else {
      res[c].eLo = minEnergy;
      res[c].eHi = 1.0e9;
    }
    Int_t kBest = 0;
    for (Int_t k = 0; k < kScanN; k++) {
      const Double_t corr = Correlation(scan[c], k);
      if (corr > res[c].bestCorr) {
        res[c].bestCorr = corr;
        kBest = k;
      }
    }
    res[c].iBest = cfg[c].pre + kScanFrom + kBest * kScanStep;
    res[c].bestRead_us = Double_t(res[c].iBest - cfg[c].pre) * kSample_ns / 1000.0;
    const Int_t kNominal =
        (cfg[c].rise + cfg[c].peakDelay - kScanFrom + kScanStep / 2) / kScanStep;
    res[c].nominalCorr = (kNominal >= 0 && kNominal < kScanN)
                             ? Correlation(scan[c], kNominal)
                             : 0.0;
  }
  // Isolated: nothing on the channel in the previous 60 us or until the end
  // of the record. Beam-like: isolated and inside the energy window.
  std::set<std::pair<Int_t, Long64_t> > isolatedKeys;
  std::vector<std::vector<Double_t> > rises(kMaxCh), polar(kMaxCh), starts(kMaxCh);
  std::vector<std::vector<Double_t> > noises(kMaxCh), noiseFrac(kMaxCh);
  for (size_t i = 0; i < hits.size(); i++) {
    Hit &h = hits[i];
    const ChResult &r = res[h.channel];
    const Double_t recordAfter_us =
        Double_t(traceLen[h.channel] - cfg[h.channel].pre) * kSample_ns / 1000.0;
    h.isolated = h.dtPrev_us > kIsolatedBefore_us &&
                 h.dtNext_us > recordAfter_us + 1.0;
    h.beamlike = h.isolated && h.energy >= r.eLo && h.energy < r.eHi;
    if (!h.isolated)
      continue;
    isolatedKeys.insert(std::make_pair(h.channel, Long64_t(h.t_ns)));
    res[h.channel].nIso++;
    if (!h.beamlike)
      continue;
    res[h.channel].nBeam++;
    rises[h.channel].push_back(h.rise_ns);
    noises[h.channel].push_back(h.noise);
    if (h.inAmp != 0.0)
      noiseFrac[h.channel].push_back(100.0 * h.noise / TMath::Abs(h.inAmp));
    polar[h.channel].push_back(h.inAmp >= 0.0 ? 1.0 : -1.0);
    starts[h.channel].push_back(Double_t(h.iStart - cfg[h.channel].pre) *
                                kSample_ns / 1000.0);
  }
  for (Int_t c = 0; c < kMaxCh; c++) {
    res[c].riseMedian_ns = Median(rises[c]);
    res[c].noise = Median(noises[c]);
    res[c].noiseFrac = Median(noiseFrac[c]);
    Double_t p = 0.0, s1 = 0.0, s2 = 0.0;
    for (size_t i = 0; i < polar[c].size(); i++) {
      p += polar[c][i];
      s1 += starts[c][i];
      s2 += starts[c][i] * starts[c][i];
    }
    const Double_t n = Double_t(polar[c].size());
    res[c].polarity = n > 0.0 ? p / n : 0.0;
    res[c].startMean_us = n > 0.0 ? s1 / n : 0.0;
    res[c].startRms_us =
        n > 1.0 ? TMath::Sqrt(TMath::Max(0.0, s2 / n - (s1 / n) * (s1 / n)))
                : 0.0;
  }

  // --- pass 2: spectra of isolated hits, shapes of beam-like ones ---------------
  std::vector<ShapeSum> shapeCh(kMaxCh);
  std::vector<ShapeSum> shapeGroup(kNGroups);
  std::vector<std::vector<Example> > examples(kNGroups);
  std::vector<Int_t> groupPre(kNGroups, 0), groupLen(kNGroups, 0);
  std::vector<std::vector<TH1D *> > est(kMaxCh, std::vector<TH1D *>(kNEst, nullptr));
  for (Int_t c = 0; c < kMaxCh; c++) {
    if (traceLen[c] > 0) {
      shapeCh[c].Init(traceLen[c]);
      for (Int_t e = 0; e < kNEst; e++) {
        est[c][e] = new TH1D(Form("est_ch%d_%d", c, e), ";Energy [ADC];Hits",
                             2048, 0.0, 32768.0);
        est[c][e]->SetDirectory(nullptr);
      }
    }
    const Int_t g = GroupOf(c);
    if (g >= 0 && traceLen[c] > 0 && groupLen[g] == 0) {
      groupLen[g] = traceLen[c];
      groupPre[g] = cfg[c].pre;
      shapeGroup[g].Init(traceLen[c]);
    }
  }
  for (size_t f = 0; f < inputs.size(); f++) {
    SOLReader reader;
    if (!reader.Open(inputs[f].Data()))
      continue;
    while (reader.ReadEvent()) {
      const SOLData &ev = reader.GetCurrentEvent();
      if (!ev.hasTraces())
        continue;
      const Int_t ch = Int_t(ev.channel);
      const Int_t g = GroupOf(ch);
      if (g < 0 || !filters[ch])
        continue;
      if (!isolatedKeys.count(std::make_pair(ch, Long64_t(ev.timestamp))))
        continue;
      const Int_t len = Int_t(ev.getSamples());
      if (len != traceLen[ch])
        continue;
      const ChSettings &s = cfg[ch];
      const ChResult &r = res[ch];
      const Int_t *in = ev.getAnalog0();
      const Double_t e = Double_t(ev.energy);
      Double_t inPre = 0.0, pol = 1.0, amp = 0.0, noise = 0.0;
      Int_t iStart = 0, i90 = 0, iExt = 0;
      if (!LocatePulse(in, len, inPre, pol, amp, noise, iStart, i90, iExt))
        continue;
      OfflineTrapezoid(in, len, inPre, pol, *filters[ch], trap);
      OfflineTrapezoid(in, len, inPre, pol, *altFilters[ch], alt);
      const Double_t trPre = MeanD(trap, TMath::Max(0, iStart - 60), 56);
      const Double_t altPre = MeanD(alt, TMath::Max(0, iStart - 60), 56);
      Int_t iMx = 0, iAlt = 0;
      const Double_t trMax = MaxFrom(trap, iStart, trPre, iMx);
      const Double_t altMax = MaxFrom(alt, iStart, altPre, iAlt);
      const Int_t iAltRead = iStart + altRise + altFlat / 2;
      est[ch][kFirmware]->Fill(e);
      est[ch][kSetRead]->Fill(trap[r.iBest] - trPre);
      est[ch][kSetMax]->Fill(trMax);
      est[ch][kAltMax]->Fill(altMax);
      if (iAltRead + 8 < len)
        est[ch][kAltRead]->Fill(MeanD(alt, iAltRead - 8, 16) - altPre);
      if (e < r.eLo || e >= r.eHi || e <= 0.0)
        continue;
      ShapeSum &sc = shapeCh[ch];
      ShapeSum &sg = shapeGroup[g];
      const Bool_t sameLen = (len == groupLen[g] && s.pre == groupPre[g]);
      for (Int_t i = 0; i < len; i++) {
        const Double_t vi = pol * (Double_t(in[i]) - inPre) / amp;
        const Double_t vt = (trap[i] - trPre) / e;
        const Double_t va = (alt[i] - altPre) / e;
        sc.input[i] += vi;
        sc.trap[i] += vt;
        sc.alt[i] += va;
        if (sameLen) {
          sg.input[i] += vi;
          sg.trap[i] += vt;
          sg.alt[i] += va;
        }
        // Start-aligned copy: sample iStart lands on kAlignRef, so the
        // trigger-walk jitter does not smear the leading edge.
        const Int_t j = i + iStart - kAlignRef;
        if (j >= 0 && j < len) {
          const Double_t vj = pol * (Double_t(in[j]) - inPre) / amp;
          sc.aligned[i] += vj;
          sc.alignedN[i] += 1.0;
          if (sameLen) {
            sg.aligned[i] += vj;
            sg.alignedN[i] += 1.0;
          }
        }
      }
      sc.n++;
      if (sameLen)
        sg.n++;
      if (Int_t(examples[g].size()) < kExamplesPerGroup) {
        Example ex;
        ex.channel = ch;
        ex.energy = e;
        ex.inPre = inPre;
        ex.inAmp = pol * amp;
        ex.trapPre = trPre;
        ex.input.resize(len);
        ex.trap = trap;
        ex.alt = alt;
        for (Int_t i = 0; i < len; i++)
          ex.input[i] = Double_t(in[i]);
        examples[g].push_back(ex);
      }
    }
  }

  // --- per-channel results ------------------------------------------------------
  TFile *fout = IO::OpenForWriting("wavecheck.root");
  TTree *tree = new TTree("hits", "per-trace features");
  Hit b;
  tree->Branch("channel", &b.channel, "channel/I");
  tree->Branch("group", &b.group, "group/I");
  tree->Branch("t_ns", &b.t_ns, "t_ns/D");
  tree->Branch("energy", &b.energy, "energy/D");
  tree->Branch("iStart", &b.iStart, "iStart/I");
  tree->Branch("inPre", &b.inPre, "inPre/D");
  tree->Branch("inAmp", &b.inAmp, "inAmp/D");
  tree->Branch("noise", &b.noise, "noise/D");
  tree->Branch("rise_ns", &b.rise_ns, "rise_ns/D");
  tree->Branch("trapPre", &b.trapPre, "trapPre/D");
  tree->Branch("trapNominal", &b.trapNominal, "trapNominal/D");
  tree->Branch("trapMax", &b.trapMax, "trapMax/D");
  tree->Branch("tTrapMax_us", &b.tTrapMax_us, "tTrapMax_us/D");
  tree->Branch("altMax", &b.altMax, "altMax/D");
  tree->Branch("altRead", &b.altRead, "altRead/D");
  tree->Branch("dtPrev_us", &b.dtPrev_us, "dtPrev_us/D");
  tree->Branch("dtNext_us", &b.dtNext_us, "dtNext_us/D");
  tree->Branch("isolated", &b.isolated, "isolated/O");
  tree->Branch("beamlike", &b.beamlike, "beamlike/O");
  for (size_t i = 0; i < hits.size(); i++) {
    b = hits[i];
    tree->Fill();
  }
  tree->Write();

  const TString header =
      Form("%-3s %-5s %-8s %7s %6s %6s %13s %4s %11s %7s %6s %6s %9s %11s %8s "
           "%8s %8s %8s %7s %8s %8s",
           "ch", "strip", "group", "hits", "iso", "beam", "E window", "pol",
           "start us", "rise ns", "noise", "n/h %", "tau us", "flat us", "nom/max", "nom corr",
           "best us", "best/max", "corr", "tail c", "pred c");
  std::cout << header << std::endl;
  summary << header << std::endl;
  std::vector<TGraph *> gIn(kMaxCh, nullptr), gTr(kMaxCh, nullptr);
  for (Int_t c = 0; c < kMaxCh; c++) {
    const Int_t g = GroupOf(c);
    if (g < 0 || shapeCh[c].n < 20)
      continue;
    const ChSettings &s = cfg[c];
    ChResult &r = res[c];
    const Double_t nBeam = Double_t(shapeCh[c].n);
    gIn[c] = TraceGraph(shapeCh[c].input, s.pre, 1.0 / nBeam,
                        Form("input_ch%d", c));
    gTr[c] = TraceGraph(shapeCh[c].trap, s.pre, 1.0 / nBeam,
                        Form("trap_ch%d", c));
    const Double_t tEnd =
        (Double_t(traceLen[c] - s.pre) - 8.0) * kSample_ns / 1000.0;
    // Preamp decay: input tail from 1.5 us past the pulse maximum.
    Double_t inMax = -1.0e30, tInMax = 0.0;
    for (Int_t i = 0; i < gIn[c]->GetN(); i++) {
      Double_t x, y;
      gIn[c]->GetPoint(i, x, y);
      if (x >= r.startMean_us && x < r.startMean_us + 10.0 && y > inMax) {
        inMax = y;
        tInMax = x;
      }
    }
    TF1 *fDecay = new TF1(Form("decay_ch%d", c), "[0] * exp(-x / [1]) + [2]",
                          tInMax + 1.5, tEnd);
    fDecay->SetParameters(1.0, 15.0, 0.0);
    fDecay->SetParLimits(1, 0.5, 500.0);
    gIn[c]->Fit(fDecay, "RQ0");
    r.tauPre_us = fDecay->GetParameter(1);
    r.tauPreErr_us = fDecay->GetParError(1);
    // Trapezoid: its maximum, the stretch within 1 % of it, the level at the
    // nominal read window and at the best-correlated sample.
    Double_t trMax = -1.0e30;
    for (Int_t i = 0; i < gTr[c]->GetN(); i++) {
      Double_t x, y;
      gTr[c]->GetPoint(i, x, y);
      if (x >= r.startMean_us && x < r.startMean_us + 12.0)
        trMax = TMath::Max(trMax, y);
    }
    r.flatFrom_us = -99.0;
    r.flatTo_us = -99.0;
    for (Int_t i = 0; i < gTr[c]->GetN(); i++) {
      Double_t x, y;
      gTr[c]->GetPoint(i, x, y);
      if (x < r.startMean_us)
        continue;
      if (y >= 0.99 * trMax) {
        if (r.flatFrom_us < -50.0)
          r.flatFrom_us = x;
        r.flatTo_us = x;
      }
    }
    Double_t peakSum = 0.0;
    const Int_t peakStart = s.pre + s.rise + s.peakDelay - s.avg / 2;
    for (Int_t i = 0; i < s.avg; i++)
      peakSum += shapeCh[c].trap[peakStart + i] / nBeam;
    r.nominalOverMax = trMax > 0.0 ? (peakSum / Double_t(s.avg)) / trMax : 0.0;
    r.bestOverMax =
        trMax > 0.0 ? (shapeCh[c].trap[r.iBest] / nBeam) / trMax : 0.0;
    // Tail after the trapezoid: from its end (rise + flat + rise past the
    // pulse start, plus the input rise time) plus a margin.
    const Double_t trapEnd_us =
        r.startMean_us + Double_t(2 * s.rise + s.flat) * kSample_ns / 1000.0 +
        r.riseMedian_ns / 1000.0 + 1.0;
    TF1 *fTail = new TF1(Form("tail_ch%d", c), "[0] * exp(-x / [1]) + [2]",
                         trapEnd_us, tEnd);
    fTail->SetParameters(0.0, r.tauPre_us, 0.0);
    fTail->FixParameter(1, r.tauPre_us);
    gTr[c]->Fit(fTail, "RQ0");
    r.tailMeasured = trMax > 0.0 ? -fTail->GetParameter(0) / trMax : 0.0;
    r.tailPredicted =
        TailAmplitude(r.tauPre_us, s.pz_us, Double_t(s.rise) * kSample_ns / 1000.0,
                      Double_t(s.rise + s.flat) * kSample_ns / 1000.0);
    for (Int_t e = 0; e < kNEst; e++)
      r.peak[e] = FitPeak(est[c][e], minEnergy);
    const TString window = r.eHi > 1.0e8 ? Form(">= %.0f", r.eLo)
                                         : Form("%.0f-%.0f", r.eLo, r.eHi);
    const TString line =
        Form("%-3d %-5d %-8s %7lld %6lld %6lld %13s %4.1f %+5.2f+-%4.2f %7.0f "
             "%6.1f %6.2f %4.1f+-%3.1f %+5.2f..%+5.2f %8.3f %8.3f %+8.2f %8.3f "
             "%7.4f %+8.3f %+8.3f",
             c, StripOf(c), kGroupTag[g], r.nHits, r.nIso, r.nBeam,
             window.Data(), r.polarity, r.startMean_us, r.startRms_us,
             r.riseMedian_ns, r.noise, r.noiseFrac, r.tauPre_us, r.tauPreErr_us, r.flatFrom_us,
             r.flatTo_us, r.nominalOverMax, r.nominalCorr, r.bestRead_us,
             r.bestOverMax, r.bestCorr, r.tailMeasured, r.tailPredicted);
    std::cout << line << std::endl;
    summary << line << std::endl;
    gIn[c]->Write();
    gTr[c]->Write();
  }
  // Resolution table: sigma / mean of the main peak per estimator.
  const TString header2 =
      Form("%-3s %-5s %-8s %6s | %-18s %-18s %-18s %-18s %-18s", "ch", "strip",
           "group", "iso", kEstName[0], kEstName[1], kEstName[2], kEstName[3],
           kEstName[4]);
  std::cout << std::endl
            << "beam peak: mean and sigma/mean per energy estimator, on the "
               "time-isolated hits"
            << std::endl
            << header2 << std::endl;
  summary << std::endl << header2 << std::endl;
  std::vector<std::vector<Double_t> > resGroup(kNGroups,
                                               std::vector<Double_t>());
  std::vector<std::vector<Double_t> > resGroupAlt(kNGroups,
                                                  std::vector<Double_t>());
  for (Int_t c = 0; c < kMaxCh; c++) {
    const Int_t g = GroupOf(c);
    if (g < 0 || !est[c][0] || res[c].nIso < 50)
      continue;
    TString line = Form("%-3d %-5d %-8s %6lld |", c, StripOf(c), kGroupTag[g],
                        res[c].nIso);
    for (Int_t e = 0; e < kNEst; e++) {
      const PeakFit &p = res[c].peak[e];
      if (p.ok)
        line += Form(" %7.0f %5.2f%%     ", p.mean, 100.0 * p.sigma / p.mean);
      else
        line += Form(" %-18s", "   -");
    }
    std::cout << line << std::endl;
    summary << line << std::endl;
    if (res[c].peak[kFirmware].ok && res[c].peak[kAltMax].ok) {
      resGroup[g].push_back(100.0 * res[c].peak[kFirmware].sigma /
                            res[c].peak[kFirmware].mean);
      resGroupAlt[g].push_back(100.0 * res[c].peak[kAltMax].sigma /
                               res[c].peak[kAltMax].mean);
    }
  }
  for (Int_t g = 0; g < kNGroups; g++) {
    if (resGroup[g].empty())
      continue;
    const TString line =
        Form("  %-8s median sigma/mean: firmware %.2f %%, alt max %.2f %% "
             "(%zu channels)",
             kGroupTag[g], Median(resGroup[g]), Median(resGroupAlt[g]),
             resGroup[g].size());
    std::cout << line << std::endl;
    summary << line << std::endl;
  }

  // --- plots -----------------------------------------------------------------
  // Averaged input, set trapezoid and alternative trapezoid per group.
  for (Int_t g = 0; g < kNGroups; g++) {
    if (shapeGroup[g].n < 20)
      continue;
    Int_t chRef = -1, nRef = 0;
    Double_t startMean = 0.0, bestRead = 0.0;
    for (Int_t c = 0; c < kMaxCh; c++) {
      if (GroupOf(c) != g || traceLen[c] != groupLen[g] || shapeCh[c].n < 20)
        continue;
      if (chRef < 0)
        chRef = c;
      startMean += res[c].startMean_us;
      bestRead += res[c].bestRead_us;
      nRef++;
    }
    if (chRef < 0)
      continue;
    startMean /= Double_t(nRef);
    bestRead /= Double_t(nRef);
    const ChSettings &s = cfg[chRef];
    const Double_t n = Double_t(shapeGroup[g].n);
    TGraph *in = TraceGraph(shapeGroup[g].input, s.pre, 1.0 / n,
                            Form("input_%s", kGroupTag[g]));
    TGraph *tr = TraceGraph(shapeGroup[g].trap, s.pre, 1.0 / n,
                            Form("trap_%s", kGroupTag[g]));
    TGraph *al = TraceGraph(shapeGroup[g].alt, s.pre, 1.0 / n,
                            Form("alt_%s", kGroupTag[g]));
    TGraph *ali = TraceGraph(shapeGroup[g].AlignedMean(), kAlignRef, 1.0,
                             Form("input_aligned_%s", kGroupTag[g]));
    in->Write();
    tr->Write();
    al->Write();
    ali->Write();
    const Double_t tRise = Double_t(s.rise) * kSample_ns / 1000.0;
    const Double_t tFlat = Double_t(s.flat) * kSample_ns / 1000.0;
    const Double_t tPeak = Double_t(s.rise + s.peakDelay) * kSample_ns / 1000.0;
    for (Int_t zoom = 0; zoom < 2; zoom++) {
      TCanvas *c1 = PlottingUtils::GetConfiguredCanvas(kFALSE);
      TGraph *inD = (TGraph *)in->Clone(Form("in_%s_%d", kGroupTag[g], zoom));
      TGraph *trD = (TGraph *)tr->Clone(Form("tr_%s_%d", kGroupTag[g], zoom));
      TGraph *alD = (TGraph *)al->Clone(Form("al_%s_%d", kGroupTag[g], zoom));
      PlottingUtils::ConfigureGraph(
          inD, kGray + 2,
          Form("%s;t #minus DAQ trigger [#mus];input / height,  trapezoid / "
               "energy",
               kGroupName[g]));
      inD->SetLineWidth(2);
      const Double_t yLo = zoom ? -0.3 : -0.5, yHi = 2.1;
      if (zoom)
        inD->GetXaxis()->SetLimits(startMean - 1.0, startMean + 10.0);
      inD->GetYaxis()->SetRangeUser(yLo, yHi);
      inD->Draw("AL");
      PlottingUtils::ConfigureGraph(alD, kGray + 1);
      alD->SetLineWidth(3);
      alD->SetLineStyle(2);
      alD->Draw("L SAME");
      PlottingUtils::ConfigureGraph(trD, kGroupColour[g]);
      trD->SetLineWidth(3);
      trD->Draw("L SAME");
      DrawMark(startMean, yLo, 1.25, kGray + 1, 3);
      DrawMark(0.0, yLo, 1.25, kBlack);
      DrawMark(tPeak, yLo, 1.25, kRed + 2);
      DrawMark(bestRead, yLo, 1.25, kGreen + 2);
      TLegend *leg = PlottingUtils::AddLegend(0.36, 0.89, 0.60, 0.87);
      leg->SetTextSize(19.0);
      leg->AddEntry(inD, Form("ADC input, %lld beam-like isolated hits",
                              shapeGroup[g].n),
                    "l");
      leg->AddEntry(trD,
                    Form("set trapezoid: rise %.3f, flat %.1f, pz %.0f #mus",
                         tRise, tFlat, s.pz_us),
                    "l");
      leg->AddEntry(alD,
                    Form("alternative: rise %.1f, flat %.1f, same pz",
                         altRise_ns / 1000.0, altFlat_ns / 1000.0),
                    "l");
      leg->AddEntry((TObject *)nullptr,
                    Form("dotted grey: pulse starts at %+.2f #mus;  black: "
                         "DAQ trigger",
                         startMean),
                    "");
      leg->AddEntry((TObject *)nullptr,
                    Form("red: nominal read, trigger + %.2f #mus", tPeak), "");
      leg->AddEntry((TObject *)nullptr,
                    Form("green: energy correlates best at %+.2f #mus", bestRead),
                    "");
      leg->Draw();
      PlottingUtils::SaveFigure(
          c1, Form(zoom ? "shapes_zoom_%s" : "shapes_%s", kGroupTag[g]), "",
          PlotSaveOptions::kLINEAR);
      delete c1;
    }
  }
  // All groups overlaid: inputs (zoom) and set trapezoids (full).
  for (Int_t which = 0; which < 2; which++) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.55, 0.89, 0.55, 0.87);
    leg->SetTextSize(19.0);
    Bool_t first = kTRUE;
    for (Int_t g = 0; g < kNGroups; g++) {
      if (shapeGroup[g].n < 20)
        continue;
      Int_t chRef = -1;
      for (Int_t cc = 0; cc < kMaxCh && chRef < 0; cc++)
        if (GroupOf(cc) == g && traceLen[cc] == groupLen[g])
          chRef = cc;
      const Double_t n = Double_t(shapeGroup[g].n);
      TGraph *gr = TraceGraph(which == 0 ? shapeGroup[g].input
                                         : shapeGroup[g].trap,
                              cfg[chRef].pre, 1.0 / n,
                              Form("all_%d_%s", which, kGroupTag[g]));
      PlottingUtils::ConfigureGraph(
          gr, kGroupColour[g],
          which == 0 ? ";t #minus DAQ trigger [#mus];ADC input / height"
                     : ";t #minus DAQ trigger [#mus];set trapezoid / energy");
      gr->SetLineWidth(3);
      if (first) {
        if (which == 0) {
          gr->GetXaxis()->SetLimits(-5.0, 8.0);
          gr->GetYaxis()->SetRangeUser(-0.2, 1.6);
        } else {
          gr->GetYaxis()->SetRangeUser(-0.5, 2.0);
        }
        gr->Draw("AL");
        first = kFALSE;
      } else {
        gr->Draw("L SAME");
      }
      leg->AddEntry(gr, Form("%s (%lld)", kGroupName[g], shapeGroup[g].n),
                    "l");
    }
    if (first) {
      delete c;
      continue;
    }
    DrawMark(0.0, which == 0 ? -0.2 : -0.5, 1.2, kBlack);
    leg->Draw();
    PlottingUtils::SaveFigure(c, which == 0 ? "inputs_all_groups"
                                            : "trapezoids_all_groups",
                              "", PlotSaveOptions::kLINEAR);
    delete c;
  }
  // Single-event examples per group: raw ADC input, set and alt trapezoids.
  for (Int_t g = 0; g < kNGroups; g++) {
    if (examples[g].empty())
      continue;
    for (Int_t which = 0; which < 3; which++) {
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      TLegend *leg = PlottingUtils::AddLegend(0.45, 0.89, 0.66, 0.87);
      leg->SetTextSize(19.0);
      const Int_t colours[3] = {kBlue + 1, kRed + 1, kGreen + 2};
      Double_t lo = 1.0e30, hi = -1.0e30;
      for (size_t e = 0; e < examples[g].size(); e++) {
        const std::vector<Double_t> &v =
            which == 0 ? examples[g][e].input
                       : (which == 1 ? examples[g][e].trap : examples[g][e].alt);
        for (size_t i = 0; i < v.size(); i++) {
          lo = TMath::Min(lo, v[i]);
          hi = TMath::Max(hi, v[i]);
        }
      }
      const Double_t pad = 0.08 * (hi - lo);
      const char *what = which == 0 ? "ADC input (raw)"
                                    : (which == 1 ? "set trapezoid"
                                                  : "alternative trapezoid");
      for (size_t e = 0; e < examples[g].size(); e++) {
        const Example &ex = examples[g][e];
        TGraph *gr = TraceGraph(
            which == 0 ? ex.input : (which == 1 ? ex.trap : ex.alt),
            cfg[ex.channel].pre, 1.0,
            Form("ex_%s_%d_%zu", kGroupTag[g], which, e));
        PlottingUtils::ConfigureGraph(
            gr, colours[e % 3],
            Form("%s, single events;t #minus DAQ trigger [#mus];%s [ADC]",
                 kGroupName[g], what));
        gr->SetLineWidth(2);
        if (e == 0) {
          gr->GetYaxis()->SetRangeUser(lo - pad, hi + pad);
          gr->Draw("AL");
        } else {
          gr->Draw("L SAME");
        }
        leg->AddEntry(gr,
                      Form("ch %d: firmware energy %.0f, input height %+.0f",
                           ex.channel, ex.energy, ex.inAmp),
                      "l");
      }
      DrawMark(0.0, lo - pad, hi + pad, kBlack);
      leg->Draw();
      const char *suffix = which == 0 ? "input" : (which == 1 ? "trap" : "alt");
      PlottingUtils::SaveFigure(c, Form("example_%s_%s", kGroupTag[g], suffix),
                                "", PlotSaveOptions::kLINEAR);
      delete c;
    }
  }
  // Energy estimators compared, per group, on the channel with most isolated
  // hits: each spectrum in units of its own fitted peak.
  for (Int_t g = 0; g < kNGroups; g++) {
    Int_t chRef = -1;
    for (Int_t c = 0; c < kMaxCh; c++)
      if (GroupOf(c) == g && res[c].peak[kFirmware].ok &&
          (chRef < 0 || res[c].nIso > res[chRef].nIso))
        chRef = c;
    if (chRef < 0)
      continue;
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.14, 0.52, 0.55, 0.87);
    leg->SetTextSize(19.0);
    const Int_t colours[kNEst] = {kBlack, kRed + 1, kOrange + 7, kBlue + 1,
                                  kAzure + 7};
    Bool_t first = kTRUE;
    Double_t ymax = 0.0;
    std::vector<TH1D *> drawn;
    for (Int_t e = 0; e < kNEst; e++) {
      const PeakFit &p = res[chRef].peak[e];
      if (!p.ok)
        continue;
      TH1D *h = new TH1D(Form("cmp_%s_%d", kGroupTag[g], e), "", 160, 0.6, 1.4);
      h->SetDirectory(nullptr);
      TH1D *src = est[chRef][e];
      for (Int_t bin = 1; bin <= src->GetNbinsX(); bin++)
        h->Fill(src->GetBinCenter(bin) / p.mean, src->GetBinContent(bin));
      ymax = TMath::Max(ymax, h->GetMaximum());
      drawn.push_back(h);
      PlottingUtils::ConfigureHistogram(
          h, colours[e],
          Form("%s, ch %d;E / fitted peak;Isolated hits", kGroupName[g],
               chRef));
      h->SetLineWidth(e == kFirmware ? 3 : 2);
      leg->AddEntry(h, Form("%s: peak %.0f, #sigma/#mu %.2f %%", kEstName[e],
                            p.mean, 100.0 * p.sigma / p.mean),
                    "l");
    }
    for (size_t i = 0; i < drawn.size(); i++) {
      if (first) {
        drawn[i]->GetYaxis()->SetRangeUser(0.0, 1.25 * ymax);
        drawn[i]->Draw("HIST");
        first = kFALSE;
      } else {
        drawn[i]->Draw("HIST SAME");
      }
    }
    if (first) {
      delete c;
      continue;
    }
    leg->Draw();
    PlottingUtils::SaveFigure(c, Form("estimators_%s", kGroupTag[g]), "",
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  // Per-channel summaries against strip number.
  const Int_t nQ = 8;
  const char *qName[nQ] = {"tau_pre_us",   "rise_ns",     "nominal_over_max",
                           "best_read_us", "tail_c",      "energy_mode",
                           "pulse_start_us", "resolution"};
  const char *qTitle[nQ] = {
      ";strip;preamp decay time from the input tail [#mus]",
      ";strip;input 10-90 % rise time, median [ns]",
      ";strip;trapezoid at the nominal read point / its maximum",
      ";strip;sample the energy correlates best with, vs trigger [#mus]",
      ";strip;trapezoid tail c (undershoot > 0); open: Eq. (c) prediction",
      ";strip;energy peak [ADC]",
      ";strip;pulse start vs DAQ trigger [#mus]",
      ";strip;beam peak #sigma/#mu [%]: filled firmware, open alt. max"};
  for (Int_t q = 0; q < nQ; q++) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.62, 0.89, 0.62, 0.87);
    leg->SetTextSize(19.0);
    Double_t lo = 1.0e30, hi = -1.0e30;
    std::vector<TGraph *> gs(kNGroups, nullptr), gp(kNGroups, nullptr);
    for (Int_t g = 0; g < kStrip17 + 1; g++) {
      gs[g] = new TGraph();
      gp[g] = new TGraph();
      for (Int_t ch = 0; ch < kMaxCh; ch++) {
        if (GroupOf(ch) != g || shapeCh[ch].n < 20)
          continue;
        const ChResult &r = res[ch];
        Double_t y = 0.0, y2 = 0.0;
        Bool_t has2 = kFALSE;
        switch (q) {
        case 0:
          y = r.tauPre_us;
          break;
        case 1:
          y = r.riseMedian_ns;
          break;
        case 2:
          y = r.nominalOverMax;
          break;
        case 3:
          y = r.bestRead_us;
          break;
        case 4:
          y = r.tailMeasured;
          y2 = r.tailPredicted;
          has2 = kTRUE;
          break;
        case 5:
          y = r.mode;
          break;
        case 6:
          y = r.startMean_us;
          break;
        default:
          if (!r.peak[kFirmware].ok)
            continue;
          y = 100.0 * r.peak[kFirmware].sigma / r.peak[kFirmware].mean;
          if (r.peak[kAltMax].ok) {
            y2 = 100.0 * r.peak[kAltMax].sigma / r.peak[kAltMax].mean;
            has2 = kTRUE;
          }
          break;
        }
        gs[g]->SetPoint(gs[g]->GetN(), Double_t(StripOf(ch)), y);
        lo = TMath::Min(lo, y);
        hi = TMath::Max(hi, y);
        if (has2) {
          gp[g]->SetPoint(gp[g]->GetN(), Double_t(StripOf(ch)), y2);
          lo = TMath::Min(lo, y2);
          hi = TMath::Max(hi, y2);
        }
      }
    }
    if (hi < lo) {
      delete c;
      continue;
    }
    const Double_t pad = 0.15 * TMath::Max(hi - lo, 1.0e-3);
    Bool_t first = kTRUE;
    for (Int_t g = 0; g < kStrip17 + 1; g++) {
      if (gs[g]->GetN() == 0)
        continue;
      PlottingUtils::ConfigureGraph(gs[g], kGroupColour[g], qTitle[q]);
      gs[g]->SetMarkerStyle(20);
      gs[g]->SetMarkerSize(1.6);
      gs[g]->SetMarkerColor(kGroupColour[g]);
      if (first) {
        gs[g]->GetXaxis()->SetLimits(-1.0, 18.0);
        gs[g]->GetYaxis()->SetRangeUser(lo - pad, hi + pad);
        gs[g]->Draw("AP");
        first = kFALSE;
      } else {
        gs[g]->Draw("P SAME");
      }
      leg->AddEntry(gs[g], kGroupName[g], "p");
      if (gp[g]->GetN() > 0) {
        gp[g]->SetMarkerStyle(24);
        gp[g]->SetMarkerSize(1.6);
        gp[g]->SetMarkerColor(kGroupColour[g]);
        gp[g]->Draw("P SAME");
      }
    }
    leg->Draw();
    PlottingUtils::SaveFigure(c, Form("bystrip_%s", qName[q]), "",
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  // Energy spectra per group, all channels of the group summed.
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
    TLegend *leg = PlottingUtils::AddLegend(0.62, 0.89, 0.55, 0.87);
    leg->SetTextSize(19.0);
    Bool_t first = kTRUE;
    Double_t xmax = 0.0;
    for (Int_t ch = 0; ch < kMaxCh; ch++)
      if (GroupOf(ch) >= 0 && spec[ch])
        xmax = TMath::Max(xmax, 1.5 * res[ch].mode);
    for (Int_t g = 0; g < kNGroups; g++) {
      TH1D *h = nullptr;
      for (Int_t ch = 0; ch < kMaxCh; ch++) {
        if (GroupOf(ch) != g || !spec[ch])
          continue;
        if (!h) {
          h = (TH1D *)spec[ch]->Clone(Form("spec_%s", kGroupTag[g]));
          h->SetDirectory(nullptr);
        } else {
          h->Add(spec[ch]);
        }
      }
      if (!h)
        continue;
      PlottingUtils::ConfigureHistogram(h, kGroupColour[g],
                                        ";Energy [ADC];Hits");
      h->GetXaxis()->SetRangeUser(0.0, TMath::Max(xmax, 4096.0));
      h->Draw(first ? "HIST" : "HIST SAME");
      first = kFALSE;
      leg->AddEntry(h, kGroupName[g], "l");
      h->Write();
    }
    leg->Draw();
    PlottingUtils::SaveFigure(c, "spectra_by_group", "",
                              PlotSaveOptions::kLOG);
    delete c;
  }
  fout->Close();
  summary.close();
  std::cout << "wrote " << absOut << std::endl;
  return 0;
}
