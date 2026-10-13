#ifndef POLEZERO_KERNELFIT_HPP
#define POLEZERO_KERNELFIT_HPP

#include "SimHit.hpp"
#include <Rtypes.h>
#include <TString.h>
#include <vector>

// A pulse-history correction fitted on the simulated hits, with a choice of
// basis, so the binning PulseHistory uses can be compared against alternatives
// on data where the true amplitude of every pulse is known.
//
// The model is PulseHistory's: the deviation of a hit from its channel mean is
// linear in the previous pulses on that channel,
//   shift(t) = c_0 + sum_over_past e_prev * sum_j c_j basis_j(t - t_prev),
// so the fit stays a linear least squares whatever the basis is. Only the
// basis differs:
//   kLogBins       one coefficient per log-spaced dt bin (what PulseHistory
//                  does, with 12 bins from 1 to 100 us);
//   kExponentials  one coefficient per fixed decay constant, which is the
//                  shape a pole-zero mismatch actually produces, plus an
//                  indicator for the dt range in which the trapezoid baseline
//                  is still frozen from the previous trigger.
// With nAmp > 1 every coefficient is split by the amplitude of the previous
// pulse, in PulseHistory's bands around multiples of the beam peak ([0, 0.5),
// [0.5, 1.5), ..., open ended), so a shift that is not linear in the previous
// amplitude can be followed.
class KernelFit {
public:
  enum Basis { kLogBins, kExponentials };

  KernelFit(Basis basis, Int_t nTerms, Double_t maxDt_us, Double_t freeze_us,
            Int_t nAmp = 1);

  // Fit on one channel of a time-ordered hit list. kFALSE if degenerate.
  // weight, when given, replaces the measured energy as the amplitude of each
  // past pulse, so a noiseless predictor can be tested against the real one.
  Bool_t Measure(const std::vector<SimHit> &hits, Int_t channel,
                 Double_t windowLo, Double_t windowHi,
                 const std::vector<Double_t> *weight = nullptr);
  // Corrected energies of that channel's hits, in the input hit order;
  // hits of other channels keep their raw energy.
  std::vector<Double_t> Apply(const std::vector<SimHit> &hits, Int_t channel,
                              const std::vector<Double_t> *weight = nullptr) const;

  Int_t NParameters() const { return Int_t(coeff_.size()); }
  Double_t Mode() const { return mode_; }
  Double_t Mean() const { return mean_; }
  TString Name() const;

private:
  // Basis values at one dt, added into out from index offset (which the
  // caller sized).
  void Evaluate(Double_t dt_us, std::vector<Double_t> &out, Int_t offset) const;
  // Basis functions per amplitude band, and the band of a previous pulse.
  Int_t TermsPerBand() const;
  Int_t BandOf(Double_t e) const;

  Basis basis_;
  Int_t nTerms_;
  Double_t maxDt_us_;
  Double_t freeze_us_;
  Int_t nAmp_;
  Double_t mode_;
  Double_t mean_;
  std::vector<Double_t> tau_us_;
  std::vector<Double_t> coeff_; // [0] is the intercept
};

#endif
