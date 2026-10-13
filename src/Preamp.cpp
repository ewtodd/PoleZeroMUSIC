#include "Preamp.hpp"
#include <TMath.h>

Preamp::Preamp(Double_t tau_ns, Double_t detRise_ns, Double_t samplePeriod_ns)
    : tau_ns_(tau_ns), detRise_ns_(detRise_ns),
      samplePeriod_ns_(samplePeriod_ns), tail_(0.0), current_(-1) {
  decayPerSample_ = TMath::Exp(-samplePeriod_ns_ / tau_ns_);
  rampPlateau_ =
      (tau_ns_ / detRise_ns_) * (1.0 - TMath::Exp(-detRise_ns_ / tau_ns_));
}

Double_t Preamp::ImpulseResponse(Double_t t_ns) const {
  if (t_ns < 0.0)
    return 0.0;
  if (t_ns < detRise_ns_)
    return (tau_ns_ / detRise_ns_) * (1.0 - TMath::Exp(-t_ns / tau_ns_));
  return rampPlateau_ * TMath::Exp(-(t_ns - detRise_ns_) / tau_ns_);
}

void Preamp::AddEvent(Double_t t_ns, Double_t amp) {
  RampEvent ev;
  ev.t_ns = t_ns;
  ev.amp = amp;
  ramping_.push_back(ev);
}

Double_t Preamp::Sample(Long64_t n) {
  if (current_ >= 0) {
    // Consecutive call: decay the tail by one sample.
    tail_ *= decayPerSample_;
  }
  current_ = n;
  const Double_t t = Double_t(n) * samplePeriod_ns_;
  Double_t value = tail_;
  // Events inside their ramp are evaluated explicitly; an event whose ramp
  // has just finished is folded into the tail state at exactly this sample.
  size_t keep = 0;
  for (size_t i = 0; i < ramping_.size(); i++) {
    const RampEvent &ev = ramping_[i];
    const Double_t dt = t - ev.t_ns;
    if (dt < 0.0) {
      ramping_[keep++] = ev;
      continue;
    }
    if (dt < detRise_ns_) {
      value += ev.amp * ImpulseResponse(dt);
      ramping_[keep++] = ev;
      continue;
    }
    const Double_t tailNow = ev.amp * ImpulseResponse(dt);
    tail_ += tailNow;
    value += tailNow;
  }
  ramping_.resize(keep);
  return value;
}

void Preamp::SkipTo(Long64_t n) {
  if (current_ < 0 || n <= current_)
    return;
  const Double_t dt_ns = Double_t(n - current_) * samplePeriod_ns_;
  tail_ *= TMath::Exp(-dt_ns / tau_ns_);
  // Sample(n) will decay by one more sample, so stand one sample before n.
  current_ = n - 1;
  tail_ /= decayPerSample_;
}

void Preamp::Reset() {
  tail_ = 0.0;
  current_ = -1;
  ramping_.clear();
}
