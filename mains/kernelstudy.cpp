// How well can the pulse-history correction do, and what limits it?
//
// On simulated data the true amplitude of every pulse is known, so the residual
// a correction leaves can be separated from the intrinsic spread of the beam
// and from the electronic noise. This compares the correction PulseHistory
// computes (12 log-spaced dt bins over 1-100 us) against the same linear model
// with a finer binning, with a longer reach, and with an exponential basis.
//
// Every correction is measured on the first half of the run and evaluated on
// the second half, so a model with more parameters cannot flatter itself by
// fitting the noise of the hits it is scored on.
//
//   kernelstudy <sim file> [<sim file> ...]
#include "Constants.hpp"
#include "InitUtils.hpp"
#include "KernelFit.hpp"
#include "PulseHistory.hpp"
#include "PulseHistoryAdapter.hpp"
#include "SimOutput.hpp"
#include <TMath.h>
#include <TSystem.h>
#include <iostream>
#include <map>
#include <vector>

namespace Constants {
void SetSimPulseHistoryAmpBins(Int_t n); // SimDatasetConfig.cpp
}

namespace {

const Double_t kWindowLo = 0.5;
const Double_t kWindowHi = 1.5;

struct Score {
  Double_t rms = 0.0;      // spread of the corrected pulse height
  Double_t rmsTruth = 0.0; // spread around the gain-scaled true amplitude
  Long64_t n = 0;
};

Score Evaluate(const std::vector<SimHit> &hits,
               const std::vector<Double_t> &energy, Int_t channel,
               Double_t mode, const std::vector<char> *mask = nullptr) {
  Score s;
  const Double_t lo = kWindowLo * mode, hi = kWindowHi * mode;
  std::vector<Double_t> keep, keepTrue;
  Double_t sum = 0.0, sumTrue = 0.0;
  for (size_t i = 0; i < hits.size(); i++) {
    if (hits[i].channel != channel || hits[i].pileup)
      continue;
    if (mask && !(*mask)[i])
      continue;
    if (energy[i] <= lo || energy[i] >= hi)
      continue;
    keep.push_back(energy[i]);
    keepTrue.push_back(hits[i].trueAmp);
    sum += energy[i];
    sumTrue += hits[i].trueAmp;
  }
  s.n = Long64_t(keep.size());
  if (s.n < 10)
    return s;
  const Double_t mean = sum / Double_t(s.n);
  const Double_t gain = sumTrue > 0.0 ? sum / sumTrue : 1.0;
  Double_t v = 0.0, vt = 0.0;
  for (size_t i = 0; i < keep.size(); i++) {
    const Double_t d = keep[i] - mean;
    const Double_t dt = keep[i] - gain * keepTrue[i];
    v += d * d;
    vt += dt * dt;
  }
  s.rms = TMath::Sqrt(v / Double_t(s.n));
  s.rmsTruth = TMath::Sqrt(vt / Double_t(s.n));
  return s;
}

// Hits of events PulseHistory calls beam-like for the channel scored (the
// other long end within PULSE_HISTORY_BEAM_LO..HI of its peak), the
// population its kernels are fitted on; set per channel in main.
const std::vector<char> *gBeamLike = nullptr;

void Report(const TString &label, Int_t npar, const std::vector<SimHit> &hits,
            const std::vector<Double_t> &energy, Int_t channel, Double_t mode) {
  const Score s = Evaluate(hits, energy, channel, mode);
  const Score b = gBeamLike ? Evaluate(hits, energy, channel, mode, gBeamLike)
                            : Score();
  std::cout << Form("    %-44s %5d %9.1f %11.1f %10.1f", label.Data(), npar,
                    s.rms, s.rmsTruth, b.rms)
            << std::endl;
}

// Fit on train, score on test.
void RunVariant(const TString &label, KernelFit &fit,
                const std::vector<SimHit> &train,
                const std::vector<SimHit> &test, Int_t channel, Double_t mode,
                const std::vector<Double_t> *trainWeight,
                const std::vector<Double_t> *testWeight) {
  if (!fit.Measure(train, channel, kWindowLo, kWindowHi, trainWeight))
    return;
  std::vector<Double_t> corrected = fit.Apply(test, channel, testWeight);
  Report(label, fit.NParameters(), test, corrected, channel, mode);
}

std::vector<Double_t> TrueAmplitudes(const std::vector<SimHit> &hits) {
  std::vector<Double_t> v(hits.size(), 0.0);
  for (size_t i = 0; i < hits.size(); i++)
    v[i] = hits[i].trueAmp;
  return v;
}

std::vector<Double_t> RawEnergies(const std::vector<SimHit> &hits) {
  std::vector<Double_t> v(hits.size(), 0.0);
  for (size_t i = 0; i < hits.size(); i++)
    v[i] = hits[i].energy;
  return v;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: kernelstudy <sim file> [...]" << std::endl;
    return 1;
  }
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                gSystem->WorkingDirectory(),
                                gSystem->WorkingDirectory(), kFALSE);
  const char *chName[2] = {"L", "R"};

  for (Int_t f = 1; f < argc; f++) {
    SimOutput::RunHeader header;
    std::vector<SimHit> hits;
    if (!SimOutput::Read(argv[f], header, hits))
      continue;
    const Double_t freeze_us =
        (2.0 * header.trapRise_ns + header.trapFlat_ns) / 1000.0;
    const size_t half = hits.size() / 2;
    std::vector<SimHit> train(hits.begin(), hits.begin() + half);
    std::vector<SimHit> test(hits.begin() + half, hits.end());
    std::cout << std::endl;
    std::cout << argv[f] << ": " << hits.size() << " hits, " << header.rate_kHz
              << " kHz, tau_L " << header.tauL_us << " us, tau_R "
              << header.tauR_us << " us, pole-zero " << header.trapDecay_us
              << " us, baseline freeze " << freeze_us << " us" << std::endl;
    std::cout << "  fitted on the first half, scored on the second"
              << std::endl;

    // PulseHistory itself, measured on the train half and applied to the test
    // half through the adapter.
    std::vector<Int_t> groups = PulseHistory::BuildGroupMap();
    // Once with one kernel linear in the previous amplitude, once with the
    // four amplitude bands the 37Cl configuration uses.
    const Int_t phBands[2] = {1, 4};
    Bool_t phOk[2] = {kFALSE, kFALSE};
    std::vector<Double_t> phEnergy[2];
    for (Int_t k = 0; k < 2; k++) {
      Constants::SetSimPulseHistoryAmpBins(phBands[k]);
      PulseHistoryAdapter::Adapted tr = PulseHistoryAdapter::Build(train);
      PulseHistoryAdapter::Adapted te = PulseHistoryAdapter::Build(test);
      PulseHistory::Result res;
      phEnergy[k] = RawEnergies(test);
      phOk[k] = PulseHistory::Measure(tr.hits, groups, res,
                                      Form("kernelstudy_amp%d", phBands[k]));
      if (phOk[k]) {
        if (gSystem->Getenv("KERNELSTUDY_PH_REPORT"))
          std::cout << PulseHistory::Report(
              res, Form("amp%d", phBands[k]), kTRUE);
        PulseHistory::Apply(te.hits, groups, res);
        phEnergy[k] = PulseHistoryAdapter::CorrectedEnergies(te, test);
      }
    }
    Constants::SetSimPulseHistoryAmpBins(1);

    const std::vector<Double_t> trainTrue = TrueAmplitudes(train);
    const std::vector<Double_t> testTrue = TrueAmplitudes(test);
    for (Int_t ch = 0; ch < 2; ch++) {
      KernelFit probe(KernelFit::kLogBins, 12, 100.0, freeze_us);
      if (!probe.Measure(train, ch, kWindowLo, kWindowHi))
        continue;
      const Double_t mode = probe.Mode();
      // Beam-like for this channel: the other long end of the same event
      // within PULSE_HISTORY_BEAM_LO..HI of its own peak (raw energy).
      KernelFit probeOther(KernelFit::kLogBins, 12, 100.0, freeze_us);
      std::vector<char> beamLike(test.size(), 0);
      if (probeOther.Measure(train, 1 - ch, kWindowLo, kWindowHi)) {
        const Double_t mo = probeOther.Mode();
        const Double_t blo = Constants::cfg.PULSE_HISTORY_BEAM_LO * mo;
        const Double_t bhi = Constants::cfg.PULSE_HISTORY_BEAM_HI * mo;
        std::map<Long64_t, Double_t> otherE;
        for (const SimHit &h : test)
          if (h.channel == 1 - ch)
            otherE[h.eventIndex] = h.energy;
        for (size_t i = 0; i < test.size(); i++) {
          std::map<Long64_t, Double_t>::const_iterator it =
              otherE.find(test[i].eventIndex);
          beamLike[i] = it != otherE.end() && it->second >= blo &&
                        it->second <= bhi;
        }
        gBeamLike = &beamLike;
      }
      std::cout << "  channel " << chName[ch] << ", beam peak " << mode
                << " ADC" << std::endl;
      std::cout << Form("    %-44s %5s %9s %11s %10s", "correction", "npar",
                        "RMS", "RMS vs true", "beam-like")
                << std::endl;
      {
        // Amplitude mix of the channel's hits, in the bands below.
        Double_t band[4] = {0, 0, 0, 0};
        for (const SimHit &h : test)
          if (h.channel == ch)
            band[TMath::Max(0, TMath::Min(3, Int_t(h.energy / mode + 0.5)))] += 1;
        const Double_t tot = band[0] + band[1] + band[2] + band[3];
        std::cout << Form("    amplitude mix: <0.5 %.3f, ~1 %.3f, ~2 %.3f, "
                          ">=2.5 %.3f (x beam peak)",
                          band[0] / tot, band[1] / tot, band[2] / tot,
                          band[3] / tot)
                  << std::endl;
      }
      Report("none", 0, test, RawEnergies(test), ch, mode);
      for (Int_t k = 0; k < 2; k++)
        if (phOk[k])
          Report(Form("PulseHistory, %d bins to %.0f us, %d amp",
                      PulseHistory::kNBins,
                      TMath::Power(10.0, PulseHistory::kLogHi) * 1e6,
                      phBands[k]),
                 1 + phBands[k] * PulseHistory::kNBins,
                 test, phEnergy[k], ch, mode);

      const Int_t nbins[3] = {12, 48, 192};
      for (Int_t k = 0; k < 3; k++) {
        KernelFit fit(KernelFit::kLogBins, nbins[k], 100.0, freeze_us);
        RunVariant(Form("%s to 100 us", fit.Name().Data()), fit, train, test,
                   ch, mode, nullptr, nullptr);
      }
      // The same binned kernel, only reaching further back: this is the
      // smallest change upstream, kLogHi and PULSE_HISTORY_APPLY_MAX_US.
      // 15 bins over 1-316 us keeps the bin width the 12 bins have over
      // 1-100 us.
      const Int_t nbins300[2] = {12, 15};
      for (Int_t k = 0; k < 2; k++) {
        KernelFit fit(KernelFit::kLogBins, nbins300[k], 316.0, freeze_us);
        RunVariant(Form("%s to 316 us", fit.Name().Data()), fit, train, test,
                   ch, mode, nullptr, nullptr);
      }
      // How far back the history has to reach: the trapezoid baseline is an
      // average over a few hundred microseconds, so pulses well beyond the
      // 100 us PulseHistory keeps still move it.
      const Double_t maxDt[3] = {100.0, 300.0, 600.0};
      const Int_t nexp[3] = {5, 7, 8};
      for (Int_t m = 0; m < 3; m++) {
        KernelFit fit(KernelFit::kExponentials, nexp[m], maxDt[m], freeze_us);
        RunVariant(Form("%s + freeze, to %.0f us", fit.Name().Data(),
                        maxDt[m]),
                   fit, train, test, ch, mode, nullptr, nullptr);
      }
      // Same model, predicting from the true amplitude of each past pulse
      // instead of its measured one: the gap is the noise the correction
      // imports through its own predictor.
      KernelFit ideal(KernelFit::kExponentials, 7, 300.0, freeze_us);
      RunVariant("same, noiseless predictor", ideal, train, test, ch, mode,
                 &trainTrue, &testTrue);
      // A kernel per amplitude band of the previous pulse (PulseHistory's
      // PULSE_HISTORY_AMP_BINS): the baseline restart is not linear in the
      // previous amplitude, so a 2x pulse need not move the baseline twice
      // as much as a 1x one.
      for (Int_t nAmp = 2; nAmp <= 4; nAmp += 2) {
        KernelFit fitB(KernelFit::kLogBins, 12, 100.0, freeze_us, nAmp);
        RunVariant(Form("%s to 100 us", fitB.Name().Data()), fitB, train,
                   test, ch, mode, nullptr, nullptr);
        KernelFit fitE(KernelFit::kExponentials, 7, 300.0, freeze_us, nAmp);
        RunVariant(Form("%s + freeze, to 300 us", fitE.Name().Data()), fitE,
                   train, test, ch, mode, nullptr, nullptr);
        KernelFit idealB(KernelFit::kExponentials, 7, 300.0, freeze_us, nAmp);
        RunVariant("same, noiseless predictor", idealB, train, test, ch, mode,
                   &trainTrue, &testTrue);
      }
      gBeamLike = nullptr;
    }
  }
  return 0;
}
