// wavenoise - electronic noise of the long strip ends, measured on the run-70
// traces (analog probe 0 = ADC input, analog probe 1 = trapezoid output).
//
// Only isolated hits enter (no earlier hit on the channel within 200 us and a
// single trigger in the trace), so the pre-trigger region is baseline. Per
// side it reports
//   * the raw ADC rms per sample around a straight-line fit to the region,
//     the quantity the simulation's noiseSigmaADC stands for;
//   * the trapezoid-output rms around its mean, the noise a baseline-restart
//     comparator on the trapezoid would see;
//   * the rms of 64-sample means of the ADC and of the trapezoid, a check on
//     how much of the noise is low frequency (white ADC noise would drop by 8);
//   * the scale between trace units and the energy word: ADC pulse height /
//     energy and trapezoid height / energy for beam-like hits.
//
//   wavenoise [--nfiles=N] [--maxtraces=N] <run-70 sol files>
#include "BinaryUtils.hpp"
#include <TH1D.h>
#include <TMath.h>
#include <iostream>
#include <vector>

namespace {

const Int_t kPreSamples = 625; // ChPreTriggerT = 5000 ns
const Int_t kRegion = 600;     // samples 0..599, clear of the trigger
const Double_t kIsolated_ns = 200000.0;
const Int_t kNLong = 8;
const Int_t kLongL[kNLong] = {1, 9, 17, 25, 33, 41, 49, 57};
const Int_t kLongR[kNLong] = {4, 12, 20, 28, 36, 44, 52, 60};
const Int_t kMaxCh = 64;

Int_t SideOf(Int_t ch) {
  for (Int_t i = 0; i < kNLong; i++) {
    if (kLongL[i] == ch)
      return 0;
    if (kLongR[i] == ch)
      return 1;
  }
  return -1;
}

// rms of v[0..n) around its least-squares straight line
Double_t DetrendedRms(const Int_t *v, Int_t n) {
  Double_t sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (Int_t i = 0; i < n; i++) {
    sx += i;
    sy += v[i];
    sxx += Double_t(i) * i;
    sxy += Double_t(i) * v[i];
  }
  const Double_t b = (n * sxy - sx * sy) / (n * sxx - sx * sx);
  const Double_t a = (sy - b * sx) / n;
  Double_t s = 0;
  for (Int_t i = 0; i < n; i++) {
    const Double_t d = v[i] - a - b * i;
    s += d * d;
  }
  return TMath::Sqrt(s / (n - 2));
}

Double_t Rms(const Int_t *v, Int_t n) {
  Double_t m = 0, s = 0;
  for (Int_t i = 0; i < n; i++)
    m += v[i];
  m /= n;
  for (Int_t i = 0; i < n; i++)
    s += (v[i] - m) * (v[i] - m);
  return TMath::Sqrt(s / (n - 1));
}

Double_t BlockMeanRms(const Int_t *v, Int_t n, Int_t block) {
  std::vector<Double_t> mu;
  for (Int_t i = 0; i + block <= n; i += block) {
    Double_t s = 0;
    for (Int_t j = 0; j < block; j++)
      s += v[i + j];
    mu.push_back(s / block);
  }
  return TMath::RMS(Int_t(mu.size()), mu.data());
}

} // namespace

int main(int argc, char **argv) {
  Long64_t maxFiles = 1000, maxTraces = 200000;
  std::vector<TString> inputs;
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a.BeginsWith("--nfiles="))
      maxFiles = TString(a(9, a.Length() - 9)).Atoll();
    else if (a.BeginsWith("--maxtraces="))
      maxTraces = TString(a(12, a.Length() - 12)).Atoll();
    else
      inputs.push_back(a);
  }
  if (Long64_t(inputs.size()) > maxFiles)
    inputs.resize(size_t(maxFiles));
  const char *side[2] = {"L", "R"};
  TH1D *hAdc[2], *hTrap[2], *hTrapBlk[2], *hAdcBlk[2], *hScaleIn[2], *hScaleTr[2];
  for (Int_t s = 0; s < 2; s++) {
    hAdc[s] = new TH1D(Form("adc%d", s), "", 4000, 0, 400);
    hTrap[s] = new TH1D(Form("trap%d", s), "", 4000, 0, 400);
    hTrapBlk[s] = new TH1D(Form("trapblk%d", s), "", 4000, 0, 400);
    hAdcBlk[s] = new TH1D(Form("adcblk%d", s), "", 4000, 0, 400);
    hScaleIn[s] = new TH1D(Form("scin%d", s), "", 4000, 0, 40);
    hScaleTr[s] = new TH1D(Form("sctr%d", s), "", 4000, 0, 40);
  }
  std::vector<Double_t> last(kMaxCh, -1.0e18);
  Long64_t used = 0;
  for (size_t f = 0; f < inputs.size() && used < maxTraces; f++) {
    SOLReader reader;
    if (!reader.Open(inputs[f].Data()))
      continue;
    while (reader.ReadEvent() && used < maxTraces) {
      const SOLData &ev = reader.GetCurrentEvent();
      const Int_t ch = Int_t(ev.channel);
      const Int_t s = SideOf(ch);
      if (s < 0)
        continue;
      const Double_t t = Double_t(ev.timestamp);
      const Double_t gap = t - last[ch];
      last[ch] = t;
      if (!ev.hasTraces() || gap < kIsolated_ns)
        continue;
      const Int_t len = Int_t(ev.getSamples());
      const Int_t *in = ev.getAnalog0();
      const Int_t *tr = ev.getAnalog1();
      const UChar_t *trig = ev.getDigital(0);
      if (!in || !tr || !trig || len < kPreSamples)
        continue;
      Int_t nTrig = 0;
      for (Int_t i = 1; i < len; i++)
        if (trig[i] && !trig[i - 1])
          nTrig++;
      if (nTrig != 1)
        continue;
      hAdc[s]->Fill(DetrendedRms(in, kRegion));
      hTrap[s]->Fill(Rms(tr, kRegion));
      hTrapBlk[s]->Fill(BlockMeanRms(tr, kRegion, 64));
      hAdcBlk[s]->Fill(BlockMeanRms(in, kRegion, 64));
      if (ev.energy > 500) {
        Double_t pre = 0, preT = 0, mx = -1e18, mxT = -1e18;
        for (Int_t i = 0; i < kRegion; i++) {
          pre += in[i];
          preT += tr[i];
        }
        pre /= kRegion;
        preT /= kRegion;
        for (Int_t i = kPreSamples; i < len && i < kPreSamples + 1200; i++) {
          mx = TMath::Max(mx, Double_t(in[i]));
          mxT = TMath::Max(mxT, Double_t(tr[i]));
        }
        hScaleIn[s]->Fill((mx - pre) / ev.energy);
        hScaleTr[s]->Fill((mxT - preT) / ev.energy);
      }
      used++;
    }
  }
  std::cout << "isolated single-trigger traces: " << used << std::endl;
  std::cout << "side  n       ADC rms/sample (median, mean)   trapezoid rms (median, mean)"
               "   64-sample-mean trapezoid rms (median)   64-sample-mean ADC rms (median)"
               "   ADC height/energy   trap height/energy"
            << std::endl;
  for (Int_t s = 0; s < 2; s++) {
    Double_t q = 0.5, mA, mT, mB, mAB, mSI, mST;
    hAdcBlk[s]->GetQuantiles(1, &mAB, &q);
    hScaleIn[s]->GetQuantiles(1, &mSI, &q);
    hScaleTr[s]->GetQuantiles(1, &mST, &q);
    hAdc[s]->GetQuantiles(1, &mA, &q);
    hTrap[s]->GetQuantiles(1, &mT, &q);
    hTrapBlk[s]->GetQuantiles(1, &mB, &q);
    std::cout << Form("  %s   %-7.0f %8.2f %8.2f %22.2f %8.2f %24.2f %24.2f %18.3f %20.3f",
                      side[s], hAdc[s]->GetEntries(), mA, hAdc[s]->GetMean(), mT,
                      hTrap[s]->GetMean(), mB, mAB, mSI, mST)
              << std::endl;
  }
  return 0;
}
