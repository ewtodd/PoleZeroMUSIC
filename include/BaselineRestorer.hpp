#ifndef POLEZERO_BASELINERESTORER_HPP
#define POLEZERO_BASELINERESTORER_HPP

#include <Rtypes.h>
#include <vector>

// DPP-PHA trapezoid baseline: the running mean of the last N trapezoid
// output samples. From each trigger the mean is frozen for the holdoff and the
// window stops taking samples; afterwards it resumes and the pre-trigger
// samples are pushed out one by one as new ones arrive.
//
// Optionally a first-order low-pass is applied on top of the running mean.
// CAEN support says the firmware has none (2026-10-19); it is kept as an
// option for the fits described in the README.
//
// Optionally the restart is conditional (UM5960: "once the trapezoid is
// returned to the baseline, the averaging restarts"): after the holdoff the
// window stays frozen until the trapezoid is back within restartThr of the
// frozen baseline value for restartHold consecutive samples, so a pole-zero
// tail larger than the threshold is never averaged in. The hold keeps the
// trapezoid's falling edge, which sweeps through the band at the end of the
// holdoff, from satisfying the condition. restartThr <= 0 restarts
// unconditionally after the holdoff.
//
// With refreezeSamples > 0 the return is checked only once, at the end of
// each freeze: if the trapezoid is not back within restartThr, the baseline
// is frozen for another refreezeSamples, and so on. This is how CAEN support
// describes EnergyFilterLFLimitation = On (2026-10-19).
class BaselineRestorer {
public:
  BaselineRestorer(Int_t meanSamples, Bool_t lfOn, Double_t lfTauNs,
                   Double_t samplePeriodNs, Double_t restartThr = 0.0,
                   Int_t restartHold = 1, Int_t refreezeSamples = 0);
  // Feed the trapezoid output of sample n. Ignored while frozen.
  void Update(Long64_t n, Double_t trapezoid);
  // Freeze from sample n for holdoffSamples; returns the frozen value.
  Double_t Freeze(Long64_t n, Int_t holdoffSamples);
  Bool_t IsFrozen(Long64_t n) const { return n < frozenUntil_; }
  Double_t Value() const;
  Bool_t Filled() const { return count_ >= Int_t(window_.size()); }
  void Reset();

private:
  // The running mean before the optional low-pass stage.
  Double_t Raw() const;
  std::vector<Double_t> window_;
  Int_t pos_;
  Int_t count_;
  Double_t sum_;
  Long64_t frozenUntil_;
  Bool_t lfOn_;
  Double_t lfAlpha_;
  Double_t lfState_;
  Double_t restartThr_;
  Bool_t awaitingReturn_;
  Int_t restartHold_;
  Int_t inBand_;
  Int_t refreezeSamples_;
};

#endif
