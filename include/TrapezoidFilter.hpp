#ifndef POLEZERO_TRAPEZOIDFILTER_HPP
#define POLEZERO_TRAPEZOIDFILTER_HPP

#include <Rtypes.h>
#include <vector>

// Jordanov-Knoll recursive trapezoid with pole-zero cancellation:
//   d[n] = v[n] - v[n-k] - v[n-l] + v[n-k-l]
//   p[n] = p[n-1] + d[n]
//   r[n] = p[n] + M d[n]
//   s[n] = s[n-1] + r[n]
// with k = rise, l = k + flat (samples) and M = 1 / (exp(T_s / tauSet) - 1).
// The output is s / (k M), so an exponential input of amplitude A and decay
// tauSet gives a flat top of height A. The whole chain is a finite-impulse
// response of the last k + l input samples, so after Reset the output is exact
// again once k + l samples have been fed.
class TrapezoidFilter {
public:
  TrapezoidFilter(Int_t riseSamples, Int_t flatSamples,
                  Double_t samplePeriod_ns, Double_t tauSet_ns);
  Double_t Process(Double_t v);
  void Reset();
  Int_t Length() const { return k_ + l_; }
  Double_t M() const { return M_; }

private:
  Int_t k_;
  Int_t l_;
  Double_t M_;
  Double_t norm_;
  std::vector<Double_t> delay_; // circular, size k + l + 1
  Int_t pos_;
  Double_t p_;
  Double_t s_;
};

#endif
