#include "TrapezoidFilter.hpp"
#include <TMath.h>

TrapezoidFilter::TrapezoidFilter(Int_t riseSamples, Int_t flatSamples,
                                 Double_t samplePeriod_ns, Double_t tauSet_ns)
    : k_(riseSamples), l_(riseSamples + flatSamples), pos_(0), p_(0.0),
      s_(0.0) {
  M_ = 1.0 / (TMath::Exp(samplePeriod_ns / tauSet_ns) - 1.0);
  norm_ = 1.0 / (Double_t(k_) * M_);
  delay_.assign(k_ + l_ + 1, 0.0);
}

Double_t TrapezoidFilter::Process(Double_t v) {
  const Int_t size = Int_t(delay_.size());
  delay_[pos_] = v;
  const Double_t v_k = delay_[(pos_ - k_ + size) % size];
  const Double_t v_l = delay_[(pos_ - l_ + size) % size];
  const Double_t v_kl = delay_[(pos_ - k_ - l_ + size) % size];
  pos_ = (pos_ + 1) % size;
  const Double_t d = v - v_k - v_l + v_kl;
  p_ += d;
  const Double_t r = p_ + M_ * d;
  s_ += r;
  return s_ * norm_;
}

void TrapezoidFilter::Reset() {
  delay_.assign(delay_.size(), 0.0);
  pos_ = 0;
  p_ = 0.0;
  s_ = 0.0;
}
