#ifndef POLEZERO_EVENTGENERATOR_HPP
#define POLEZERO_EVENTGENERATOR_HPP

#include "Config.hpp"
#include <TRandom3.h>
#include <vector>

// Poisson beam: exponential dt between particles, one Gaussian energy per
// particle, shared by the L and R channels so the two are directly comparable.
class EventGenerator {
public:
  explicit EventGenerator(const Config &cfg);
  void Generate();

  Long64_t NEvents() const { return Long64_t(time_ns_.size()); }
  Double_t Time_ns(Long64_t i) const { return time_ns_[i]; }
  Double_t Amplitude(Long64_t i) const { return amplitude_[i]; }
  const std::vector<Double_t> &Times_ns() const { return time_ns_; }
  const std::vector<Double_t> &Amplitudes() const { return amplitude_; }

private:
  const Config &cfg_;
  TRandom3 rng_;
  std::vector<Double_t> time_ns_;
  std::vector<Double_t> amplitude_;
};

#endif
