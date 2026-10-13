// Optimal (matched, noise-weighted) filter for the pulse height on the raw
// SOL traces, against the firmware trapezoid on the same hits.
//
// Model: a record window is v = A s + n, with s the channel group's pulse
// template (the pulse-start-aligned average written by wavecheck) and n
// stationary noise with power spectrum N(f). The maximum-likelihood height is
//   A = sum_f Re(V_f S_f*) / N_f  /  sum_f |S_f|^2 / N_f,
// evaluated in the time domain as a cross-correlation with the weights
// w = IFFT(S* / N) at the alignment the plain matched filter prefers. The
// f = 0 bin is dropped, which makes the estimate insensitive to the baseline.
// The noise-only resolution is 1 / sqrt(sum_f |S_f|^2 / N_f).
//
// N(f) is the electronic noise, measured on the pulse-free samples before
// each pulse (128-sample stretches, so resolved above ~1 MHz and taken flat
// below). The residual v - A s of a plain matched-filter pass is also
// spectrum-analysed and plotted: it contains the event-to-event variation of
// the pulse shape, which is not noise a linear filter can be optimal against.
//
//   optfilter --settings=<XSetting.dat> --template=<wavecheck.root> --out=<dir>
//             [--nfiles=N] [--maxtraces=N] [--minenergy=ADC] <sol files...>
#include "BinaryUtils.hpp"
#include "IOUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "SolWave.hpp"
#include <TCanvas.h>
#include <TComplex.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVirtualFFT.h>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace {

using namespace solwave;

const Int_t kWindow = 2048;  // samples, 16.4 us
const Int_t kBefore = 125;   // samples of the window before the pulse start
const Int_t kAlignRef = 250; // where wavecheck put the pulse start
const Int_t kShiftScan = 8;  // +- samples searched for the best alignment
const Int_t kNFreq = kWindow / 2 + 1;
const Int_t kPre = 128; // pulse-free stretch used for the noise spectrum
const Int_t kNFreqPre = kPre / 2 + 1;
const Int_t kPreGap = 24; // samples between that stretch and the pulse start

struct Template {
  Bool_t ok = kFALSE;
  std::vector<Double_t> s;      // unit height, kWindow samples
  std::vector<Double_t> re, im; // its spectrum
};

struct ChannelFilter {
  Long64_t nRes = 0, nPre = 0;
  std::vector<Double_t> resPSD; // mean |R_f|^2 of the residuals, kNFreq bins
  std::vector<Double_t> prePSD; // mean |N_f|^2 of pulse-free stretches
  std::vector<Double_t> noise;  // prePSD brought to the kWindow grid
  std::vector<Double_t> w;      // time-domain weights
  Double_t norm = 0.0;          // sum_k w_k s_k
  Double_t sigmaNoise = 0.0;    // noise-only height error, ADC
  TH1D *hFirm = nullptr, *hWhite = nullptr, *hOpt = nullptr;
};

struct PeakFit {
  Double_t mean = 0.0, sigma = 0.0;
  Bool_t ok = kFALSE;
};

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

void ForwardFFT(TVirtualFFT *fft, const std::vector<Double_t> &x, Int_t nFreq,
                std::vector<Double_t> &re, std::vector<Double_t> &im) {
  std::vector<Double_t> buf(x);
  fft->SetPoints(&buf[0]);
  fft->Transform();
  re.resize(nFreq);
  im.resize(nFreq);
  for (Int_t k = 0; k < nFreq; k++)
    fft->GetPointComplex(k, re[k], im[k]);
}

// The window of one record, baseline taken out and the pulse made positive.
Bool_t Window(const Int_t *in, Int_t len, Int_t iStart, Double_t pre,
              Double_t pol, std::vector<Double_t> &v) {
  const Int_t from = iStart - kBefore;
  if (from < 0 || from + kWindow > len)
    return kFALSE;
  v.resize(kWindow);
  for (Int_t k = 0; k < kWindow; k++)
    v[k] = pol * (Double_t(in[from + k]) - pre);
  return kTRUE;
}

Double_t Correlate(const std::vector<Double_t> &w,
                   const std::vector<Double_t> &v, Int_t shift) {
  Double_t s = 0.0;
  for (Int_t k = 0; k < kWindow; k++) {
    const Int_t j = k + shift;
    if (j >= 0 && j < kWindow)
      s += w[k] * v[j];
  }
  return s;
}

// Frequency of bin k of an n-point transform, in MHz.
Double_t FreqMHz(Int_t k, Int_t n) {
  return Double_t(k) / (Double_t(n) * kSample_ns * 1.0e-3);
}

} // namespace

int main(int argc, char **argv) {
  TString out = "output/optfilter", settingsPath, templatePath;
  Long64_t maxFiles = 1000, maxTraces = -1;
  Double_t minEnergy = 500.0;
  std::vector<TString> inputs;
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a.BeginsWith("--out="))
      out = a(6, a.Length() - 6);
    else if (a.BeginsWith("--settings="))
      settingsPath = a(11, a.Length() - 11);
    else if (a.BeginsWith("--template="))
      templatePath = a(11, a.Length() - 11);
    else if (a.BeginsWith("--nfiles="))
      maxFiles = TString(a(9, a.Length() - 9)).Atoll();
    else if (a.BeginsWith("--maxtraces="))
      maxTraces = TString(a(12, a.Length() - 12)).Atoll();
    else if (a.BeginsWith("--minenergy="))
      minEnergy = TString(a(12, a.Length() - 12)).Atof();
    else
      inputs.push_back(a);
  }
  if (inputs.empty() || settingsPath.Length() == 0 ||
      templatePath.Length() == 0) {
    std::cerr << "usage: optfilter --settings=XSetting.dat "
                 "--template=wavecheck.root --out=dir [--nfiles=N] "
                 "[--maxtraces=N] [--minenergy=ADC] <sol files>"
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
  std::vector<ChSettings> cfg = ReadSettings(settingsPath);

  Int_t n = kWindow, nPre = kPre;
  TVirtualFFT *fwd = TVirtualFFT::FFT(1, &n, "R2C ES K");
  TVirtualFFT *bwd = TVirtualFFT::FFT(1, &n, "C2R ES K");
  TVirtualFFT *fwdPre = TVirtualFFT::FFT(1, &nPre, "R2C ES K");
  if (!fwd || !bwd || !fwdPre) {
    std::cerr << "no FFT engine in this ROOT build" << std::endl;
    return 1;
  }

  // Templates per group from the aligned averages, and the time-isolated
  // hit list, both from wavecheck.
  TFile tf(templatePath);
  if (tf.IsZombie()) {
    std::cerr << "cannot open " << templatePath << std::endl;
    return 1;
  }
  Template tmpl[kNGroups];
  for (Int_t g = 0; g < kNGroups; g++) {
    TGraph *gr = (TGraph *)tf.Get(Form("input_aligned_%s", kGroupTag[g]));
    if (!gr || gr->GetN() < kAlignRef - kBefore + kWindow)
      continue;
    tmpl[g].s.resize(kWindow);
    Double_t base = 0.0;
    for (Int_t k = 0; k < 100; k++)
      base += gr->GetY()[kAlignRef - kBefore + k];
    base /= 100.0;
    Double_t mx = 0.0;
    for (Int_t k = 0; k < kWindow; k++) {
      tmpl[g].s[k] = gr->GetY()[kAlignRef - kBefore + k] - base;
      mx = TMath::Max(mx, tmpl[g].s[k]);
    }
    if (mx <= 0.0)
      continue;
    for (Int_t k = 0; k < kWindow; k++)
      tmpl[g].s[k] /= mx;
    ForwardFFT(fwd, tmpl[g].s, kNFreq, tmpl[g].re, tmpl[g].im);
    tmpl[g].ok = kTRUE;
  }
  std::map<std::pair<Int_t, Long64_t>, Bool_t> isolated; // -> beam-like
  {
    TTree *t = (TTree *)tf.Get("hits");
    if (!t) {
      std::cerr << "no hits tree in " << templatePath << std::endl;
      return 1;
    }
    Int_t channel = 0;
    Double_t t_ns = 0.0;
    Bool_t iso = kFALSE, beam = kFALSE;
    t->SetBranchAddress("channel", &channel);
    t->SetBranchAddress("t_ns", &t_ns);
    t->SetBranchAddress("isolated", &iso);
    t->SetBranchAddress("beamlike", &beam);
    for (Long64_t i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i);
      if (iso)
        isolated[std::make_pair(channel, Long64_t(t_ns))] = beam;
    }
  }
  std::cout << isolated.size() << " time-isolated hits listed by wavecheck"
            << std::endl;

  std::vector<ChannelFilter> filt(kMaxCh);
  for (Int_t c = 0; c < kMaxCh; c++) {
    filt[c].resPSD.assign(kNFreq, 0.0);
    filt[c].prePSD.assign(kNFreqPre, 0.0);
  }

  // --- pass 1: noise spectra on beam-like hits -----------------------------------
  std::vector<Double_t> v, r, re, im, q;
  for (size_t f = 0; f < inputs.size(); f++) {
    SOLReader reader;
    if (!reader.Open(inputs[f].Data()))
      continue;
    Long64_t nThisFile = 0;
    while (reader.ReadEvent()) {
      const SOLData &ev = reader.GetCurrentEvent();
      if (!ev.hasTraces())
        continue;
      if (maxTraces > 0 && nThisFile >= maxTraces)
        break;
      nThisFile++;
      const Int_t ch = Int_t(ev.channel);
      const Int_t g = GroupOf(ch);
      if (g < 0 || !tmpl[g].ok)
        continue;
      std::map<std::pair<Int_t, Long64_t>, Bool_t>::const_iterator it =
          isolated.find(std::make_pair(ch, Long64_t(ev.timestamp)));
      if (it == isolated.end() || !it->second)
        continue;
      const Int_t len = Int_t(ev.getSamples());
      const Int_t *in = ev.getAnalog0();
      Double_t pre = 0.0, pol = 1.0, amp = 0.0, noise = 0.0;
      Int_t iStart = 0, i90 = 0, iExt = 0;
      if (!in ||
          !LocatePulse(in, len, pre, pol, amp, noise, iStart, i90, iExt))
        continue;
      ChannelFilter &cf = filt[ch];
      // Electronic noise: the pulse-free stretch before the start.
      const Int_t qFrom = iStart - kPreGap - kPre;
      if (qFrom >= 0) {
        q.resize(kPre);
        Double_t qm = 0.0;
        for (Int_t k = 0; k < kPre; k++) {
          q[k] = Double_t(in[qFrom + k]);
          qm += q[k];
        }
        qm /= Double_t(kPre);
        for (Int_t k = 0; k < kPre; k++)
          q[k] -= qm;
        ForwardFFT(fwdPre, q, kNFreqPre, re, im);
        for (Int_t k = 0; k < kNFreqPre; k++)
          cf.prePSD[k] += re[k] * re[k] + im[k] * im[k];
        cf.nPre++;
      }
      // Residual after a plain matched filter, for the shape-variation plot.
      if (!Window(in, len, iStart, pre, pol, v))
        continue;
      const std::vector<Double_t> &s = tmpl[g].s;
      Double_t sv = 0.0, ss = 0.0;
      for (Int_t k = 0; k < kWindow; k++) {
        sv += s[k] * v[k];
        ss += s[k] * s[k];
      }
      const Double_t a0 = sv / ss;
      r.resize(kWindow);
      for (Int_t k = 0; k < kWindow; k++)
        r[k] = v[k] - a0 * s[k];
      ForwardFFT(fwd, r, kNFreq, re, im);
      for (Int_t k = 0; k < kNFreq; k++)
        cf.resPSD[k] += re[k] * re[k] + im[k] * im[k];
      cf.nRes++;
    }
  }
  // Weights per channel: w = IFFT(S* / N), f = 0 dropped. The kPre-point
  // spectrum is carried to the kWindow grid: per-bin power of stationary
  // noise scales with the transform length, and below its first bin the
  // spectrum is taken flat.
  for (Int_t c = 0; c < kMaxCh; c++) {
    ChannelFilter &cf = filt[c];
    const Int_t g = GroupOf(c);
    if (g < 0 || !tmpl[g].ok || cf.nPre < 100 || cf.nRes < 100)
      continue;
    for (Int_t k = 0; k < kNFreqPre; k++)
      cf.prePSD[k] /= Double_t(cf.nPre);
    for (Int_t k = 0; k < kNFreq; k++)
      cf.resPSD[k] /= Double_t(cf.nRes);
    cf.noise.assign(kNFreq, 0.0);
    const Double_t ratio = Double_t(kWindow) / Double_t(kPre);
    for (Int_t k = 0; k < kNFreq; k++) {
      const Double_t x = Double_t(k) / ratio; // position on the kPre grid
      Double_t p = 0.0;
      if (x < 1.0)
        p = cf.prePSD[1];
      else {
        const Int_t j = TMath::Min(Int_t(x), kNFreqPre - 2);
        const Double_t fr = x - Double_t(j);
        p = cf.prePSD[j] * (1.0 - fr) + cf.prePSD[j + 1] * fr;
      }
      cf.noise[k] = p * ratio;
    }
    Double_t sumSN = 0.0;
    for (Int_t k = 0; k < kNFreq; k++) {
      const Double_t nk = cf.noise[k];
      TComplex wk(0.0, 0.0);
      if (k > 0 && nk > 0.0)
        wk = TComplex(tmpl[g].re[k], -tmpl[g].im[k]) / nk;
      bwd->SetPointComplex(k, wk);
      if (k == 0 || nk <= 0.0)
        continue;
      const Double_t s2 =
          tmpl[g].re[k] * tmpl[g].re[k] + tmpl[g].im[k] * tmpl[g].im[k];
      sumSN += (k == kNFreq - 1 ? 1.0 : 2.0) * s2 / nk;
    }
    bwd->Transform();
    // sum_f V_f W_f = N sum_k v[k] IFFT(W)[-k]: the weights are the inverse
    // transform read backwards (for white noise that is s(t) itself).
    std::vector<Double_t> tmp(kWindow, 0.0);
    bwd->GetPoints(&tmp[0]);
    cf.w.resize(kWindow);
    for (Int_t k = 0; k < kWindow; k++)
      cf.w[k] = tmp[(kWindow - k) % kWindow];
    cf.norm = 0.0;
    for (Int_t k = 0; k < kWindow; k++)
      cf.norm += cf.w[k] * tmpl[g].s[k];
    cf.sigmaNoise = sumSN > 0.0 ? 1.0 / TMath::Sqrt(sumSN) : 0.0;
    cf.hFirm = new TH1D(Form("firm_ch%d", c), ";Energy [ADC];Hits", 8192, 0.0,
                        32768.0);
    cf.hWhite = new TH1D(Form("white_ch%d", c), ";Height [ADC];Hits", 8192,
                         0.0, 32768.0);
    cf.hOpt = new TH1D(Form("opt_ch%d", c), ";Height [ADC];Hits", 8192, 0.0,
                       32768.0);
    cf.hFirm->SetDirectory(nullptr);
    cf.hWhite->SetDirectory(nullptr);
    cf.hOpt->SetDirectory(nullptr);
  }

  // --- pass 2: heights on all time-isolated hits --------------------------------
  for (size_t f = 0; f < inputs.size(); f++) {
    SOLReader reader;
    if (!reader.Open(inputs[f].Data()))
      continue;
    Long64_t nThisFile = 0;
    while (reader.ReadEvent()) {
      const SOLData &ev = reader.GetCurrentEvent();
      if (!ev.hasTraces())
        continue;
      if (maxTraces > 0 && nThisFile >= maxTraces)
        break;
      nThisFile++;
      const Int_t ch = Int_t(ev.channel);
      const Int_t g = GroupOf(ch);
      ChannelFilter &cf = filt[ch];
      if (g < 0 || !cf.hOpt)
        continue;
      if (!isolated.count(std::make_pair(ch, Long64_t(ev.timestamp))))
        continue;
      const Int_t len = Int_t(ev.getSamples());
      const Int_t *in = ev.getAnalog0();
      Double_t pre = 0.0, pol = 1.0, amp = 0.0, noise = 0.0;
      Int_t iStart = 0, i90 = 0, iExt = 0;
      if (!in ||
          !LocatePulse(in, len, pre, pol, amp, noise, iStart, i90, iExt))
        continue;
      if (!Window(in, len, iStart, pre, pol, v))
        continue;
      const std::vector<Double_t> &s = tmpl[g].s;
      // Alignment from the plain matched filter, whose output against
      // shift is a single smooth peak; the optimal filter is read there.
      Double_t ss = 0.0;
      for (Int_t k = 0; k < kWindow; k++)
        ss += s[k] * s[k];
      Double_t bestWhite = -1.0e30;
      Int_t bestShift = 0;
      for (Int_t shift = -kShiftScan; shift <= kShiftScan; shift++) {
        const Double_t a = Correlate(s, v, shift) / ss;
        if (a > bestWhite) {
          bestWhite = a;
          bestShift = shift;
        }
      }
      cf.hFirm->Fill(Double_t(ev.energy));
      cf.hWhite->Fill(bestWhite);
      cf.hOpt->Fill(Correlate(cf.w, v, bestShift) / cf.norm);
    }
  }

  // --- results -----------------------------------------------------------------
  TFile *fout = IO::OpenForWriting("optfilter.root");
  const TString header =
      Form("%-3s %-5s %-8s %6s | %-19s %-19s %-19s %10s", "ch", "strip",
           "group", "n", "firmware", "matched (white)", "optimal (noise)",
           "noise-only");
  std::cout << header << std::endl;
  std::vector<std::vector<Double_t> > resF(kNGroups), resW(kNGroups),
      resO(kNGroups), resN(kNGroups);
  std::vector<PeakFit> pkF(kMaxCh), pkW(kMaxCh), pkO(kMaxCh);
  for (Int_t c = 0; c < kMaxCh; c++) {
    ChannelFilter &cf = filt[c];
    if (!cf.hOpt)
      continue;
    const Int_t g = GroupOf(c);
    pkF[c] = FitPeak(cf.hFirm, minEnergy);
    pkW[c] = FitPeak(cf.hWhite, 0.1 * minEnergy);
    pkO[c] = FitPeak(cf.hOpt, 0.1 * minEnergy);
    TString line = Form("%-3d %-5d %-8s %6lld |", c, StripOf(c), kGroupTag[g],
                        Long64_t(cf.hOpt->GetEntries()));
    const PeakFit *pks[3] = {&pkF[c], &pkW[c], &pkO[c]};
    for (Int_t e = 0; e < 3; e++) {
      if (pks[e]->ok)
        line += Form(" %7.0f %5.2f%%      ", pks[e]->mean,
                     100.0 * pks[e]->sigma / pks[e]->mean);
      else
        line += Form(" %-19s", "   -");
    }
    line += pkO[c].ok ? Form(" %8.3f%%", 100.0 * cf.sigmaNoise / pkO[c].mean)
                      : "        -";
    std::cout << line << std::endl;
    if (pkF[c].ok && pkW[c].ok && pkO[c].ok) {
      resF[g].push_back(100.0 * pkF[c].sigma / pkF[c].mean);
      resW[g].push_back(100.0 * pkW[c].sigma / pkW[c].mean);
      resO[g].push_back(100.0 * pkO[c].sigma / pkO[c].mean);
      resN[g].push_back(100.0 * cf.sigmaNoise / pkO[c].mean);
    }
    cf.hFirm->Write();
    cf.hWhite->Write();
    cf.hOpt->Write();
  }
  for (Int_t g = 0; g < kNGroups; g++) {
    if (resF[g].empty())
      continue;
    std::sort(resF[g].begin(), resF[g].end());
    std::sort(resW[g].begin(), resW[g].end());
    std::sort(resO[g].begin(), resO[g].end());
    std::sort(resN[g].begin(), resN[g].end());
    const size_t m = resF[g].size() / 2;
    std::cout << Form("  %-8s median sigma/mean: firmware %.2f %%, matched "
                      "%.2f %%, optimal %.2f %%, noise-only %.3f %% (%zu "
                      "channels)",
                      kGroupTag[g], resF[g][m], resW[g][m], resO[g][m],
                      resN[g][m], resF[g].size())
              << std::endl;
  }

  // --- plots -----------------------------------------------------------------
  // Spectra per group on the channel with most hits, each in units of its
  // fitted peak.
  for (Int_t g = 0; g < kNGroups; g++) {
    Int_t chRef = -1;
    for (Int_t c = 0; c < kMaxCh; c++)
      if (GroupOf(c) == g && pkF[c].ok && pkO[c].ok &&
          (chRef < 0 ||
           filt[c].hOpt->GetEntries() > filt[chRef].hOpt->GetEntries()))
        chRef = c;
    if (chRef < 0)
      continue;
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.14, 0.52, 0.62, 0.87);
    leg->SetTextSize(19.0);
    TH1D *src[3] = {filt[chRef].hFirm, filt[chRef].hWhite, filt[chRef].hOpt};
    const PeakFit *pks[3] = {&pkF[chRef], &pkW[chRef], &pkO[chRef]};
    const char *names[3] = {"firmware trapezoid", "matched filter (white)",
                            "optimal filter (measured noise)"};
    const Int_t colours[3] = {kBlack, kGray + 2, kRed + 1};
    Double_t ymax = 0.0;
    std::vector<TH1D *> drawn;
    for (Int_t e = 0; e < 3; e++) {
      if (!pks[e]->ok)
        continue;
      TH1D *h =
          new TH1D(Form("cmp_%s_%d", kGroupTag[g], e), "", 100, 0.6, 1.4);
      h->SetDirectory(nullptr);
      for (Int_t bin = 1; bin <= src[e]->GetNbinsX(); bin++)
        h->Fill(src[e]->GetBinCenter(bin) / pks[e]->mean,
                src[e]->GetBinContent(bin));
      ymax = TMath::Max(ymax, h->GetMaximum());
      PlottingUtils::ConfigureHistogram(
          h, colours[e],
          Form("%s, ch %d;E / fitted peak;Isolated hits", kGroupName[g],
               chRef));
      h->SetLineWidth(e == 2 ? 3 : 2);
      leg->AddEntry(h, Form("%s: #sigma/#mu %.2f %%", names[e],
                            100.0 * pks[e]->sigma / pks[e]->mean),
                    "l");
      drawn.push_back(h);
    }
    leg->AddEntry((TObject *)nullptr,
                  Form("electronic-noise limit %.3f %%",
                       100.0 * filt[chRef].sigmaNoise / pkO[chRef].mean),
                  "");
    for (size_t i = 0; i < drawn.size(); i++) {
      if (i == 0) {
        drawn[i]->GetYaxis()->SetRangeUser(0.0, 1.25 * ymax);
        drawn[i]->Draw("HIST");
      } else {
        drawn[i]->Draw("HIST SAME");
      }
    }
    leg->Draw();
    PlottingUtils::SaveFigure(c, Form("optfilter_%s", kGroupTag[g]), "",
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  // Noise spectra and the template spectrum for the best long-L channel,
  // and its weights in time.
  {
    Int_t chRef = -1;
    for (Int_t c = 0; c < kMaxCh; c++)
      if (GroupOf(c) == kLongL && filt[c].hOpt &&
          (chRef < 0 || filt[c].nPre > filt[chRef].nPre))
        chRef = c;
    if (chRef >= 0) {
      const ChannelFilter &cf = filt[chRef];
      const Template &tp = tmpl[kLongL];
      // Amplitude spectral densities per sample: sqrt(PSD / n) has the units
      // of the trace, so the three curves are on one footing.
      TGraph *gn = new TGraph(), *gr = new TGraph(), *gs = new TGraph();
      for (Int_t k = 1; k < kNFreq; k++) {
        gn->SetPoint(gn->GetN(), FreqMHz(k, kWindow),
                     TMath::Sqrt(cf.noise[k] / Double_t(kWindow)));
        gr->SetPoint(gr->GetN(), FreqMHz(k, kWindow),
                     TMath::Sqrt(cf.resPSD[k] / Double_t(kWindow)));
        gs->SetPoint(gs->GetN(), FreqMHz(k, kWindow),
                     pkO[chRef].mean *
                         TMath::Sqrt((tp.re[k] * tp.re[k] + tp.im[k] * tp.im[k]) /
                                     Double_t(kWindow)));
      }
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kTRUE);
      c->SetLogx(1);
      PlottingUtils::ConfigureGraph(
          gs, kBlue + 1,
          Form("long L, ch %d;frequency [MHz];amplitude spectral density "
               "[ADC per sample]",
               chRef));
      gs->SetLineWidth(2);
      gs->GetYaxis()->SetRangeUser(0.3, 3.0e4);
      gs->Draw("AL");
      PlottingUtils::ConfigureGraph(gn, kRed + 1);
      gn->SetLineWidth(2);
      gn->Draw("L SAME");
      PlottingUtils::ConfigureGraph(gr, kGray + 2);
      gr->SetLineWidth(2);
      gr->Draw("L SAME");
      TLegend *leg = PlottingUtils::AddLegend(0.45, 0.89, 0.66, 0.87);
      leg->SetTextSize(19.0);
      leg->AddEntry(gs, "pulse template at the beam-peak height", "l");
      leg->AddEntry(gn, "electronic noise, pulse-free stretches", "l");
      leg->AddEntry(gr, "residual after the matched filter (shape variation)",
                    "l");
      leg->Draw();
      PlottingUtils::SaveFigure(c, "optfilter_spectra", "",
                                PlotSaveOptions::kLOG);
      delete c;
      TGraph *gw = new TGraph();
      Double_t wmax = 0.0;
      for (Int_t k = 0; k < kWindow; k++)
        wmax = TMath::Max(wmax, TMath::Abs(cf.w[k]));
      for (Int_t k = 0; k < kWindow; k++)
        gw->SetPoint(k, Double_t(k - kBefore) * kSample_ns / 1000.0,
                     cf.w[k] / wmax);
      TGraph *gt = new TGraph();
      for (Int_t k = 0; k < kWindow; k++)
        gt->SetPoint(k, Double_t(k - kBefore) * kSample_ns / 1000.0, tp.s[k]);
      TCanvas *c2 = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureGraph(
          gw, kRed + 1,
          Form("long L, ch %d;t #minus pulse start [#mus];weight (unit peak), "
               "template",
               chRef));
      gw->SetLineWidth(2);
      gw->GetYaxis()->SetRangeUser(-1.2, 1.2);
      gw->Draw("AL");
      PlottingUtils::ConfigureGraph(gt, kBlue + 1);
      gt->SetLineWidth(2);
      gt->Draw("L SAME");
      TLegend *leg2 = PlottingUtils::AddLegend(0.55, 0.89, 0.72, 0.87);
      leg2->SetTextSize(19.0);
      leg2->AddEntry(gw, "optimal-filter weights", "l");
      leg2->AddEntry(gt, "template", "l");
      leg2->Draw();
      PlottingUtils::SaveFigure(c2, "optfilter_weights", "",
                                PlotSaveOptions::kLINEAR);
      delete c2;
    }
  }
  fout->Close();
  std::cout << "wrote " << absOut << std::endl;
  return 0;
}
