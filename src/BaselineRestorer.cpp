#include "BaselineRestorer.hpp"
#include <TMath.h>

BaselineRestorer::BaselineRestorer(Int_t meanSamples, Bool_t lfOn,
                                   Double_t lfTauNs, Double_t samplePeriodNs,
                                   Double_t restartThr, Int_t restartHold)
    : window_(meanSamples > 0 ? meanSamples : 1, 0.0), pos_(0), count_(0),
      sum_(0.0), frozenUntil_(-1), lfOn_(lfOn && lfTauNs > 0.0),
      lfAlpha_(lfOn_ ? samplePeriodNs / (samplePeriodNs + lfTauNs) : 0.0),
      lfState_(0.0), restartThr_(restartThr), awaitingReturn_(kFALSE),
      restartHold_(restartHold > 0 ? restartHold : 1), inBand_(0) {}

void BaselineRestorer::Update(Long64_t n, Double_t trapezoid) {
  if (n < frozenUntil_)
    return;
  // No baseline exists yet to return to.
  if (awaitingReturn_ && count_ == 0)
    awaitingReturn_ = kFALSE;
  if (awaitingReturn_) {
    if (TMath::Abs(trapezoid - Value()) > restartThr_) {
      inBand_ = 0;
      return;
    }
    if (++inBand_ < restartHold_)
      return;
    awaitingReturn_ = kFALSE;
  }
  const Int_t size = Int_t(window_.size());
  if (count_ >= size)
    sum_ -= window_[pos_];
  else
    count_++;
  window_[pos_] = trapezoid;
  sum_ += trapezoid;
  pos_ = (pos_ + 1) % size;
  // The low-pass starts from the first running mean rather than from 0.
  if (lfOn_)
    lfState_ = count_ == 1 ? Raw() : lfState_ + (Raw() - lfState_) * lfAlpha_;
}

Double_t BaselineRestorer::Raw() const {
  if (count_ == 0)
    return 0.0;
  return sum_ / Double_t(count_);
}

Double_t BaselineRestorer::Freeze(Long64_t n, Int_t holdoffSamples) {
  const Long64_t until = n + holdoffSamples;
  if (until > frozenUntil_)
    frozenUntil_ = until;
  if (restartThr_ > 0.0) {
    awaitingReturn_ = kTRUE;
    inBand_ = 0;
  }
  return Value();
}

Double_t BaselineRestorer::Value() const {
  if (lfOn_)
    return lfState_;
  return Raw();
}

void BaselineRestorer::Reset() {
  window_.assign(window_.size(), 0.0);
  pos_ = 0;
  count_ = 0;
  sum_ = 0.0;
  frozenUntil_ = -1;
  lfState_ = 0.0;
  awaitingReturn_ = kFALSE;
  inBand_ = 0;
}
