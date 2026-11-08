#ifndef POLEZERO_EVENTGENERATOR_HPP
#define POLEZERO_EVENTGENERATOR_HPP

#include "Config.hpp"
#include <TRandom3.h>
#include <vector>

// Poisson beam: exponential dt between particles, one Gaussian energy per
// particle, shared by the L and R channels so the two are directly comparable.
// A fraction cfg.reactionFraction of the particles react upstream of the two
// long ends, which then see a pair of deposits drawn from cfg.reactionPairs
// instead. Reactions use a random stream of their own, so the beam particles
// are the same, event for event, as in a run without them.
class EventGenerator {
public:
  explicit EventGenerator(const Config &cfg);
  // False if reactions are asked for and their pairs cannot be read.
  Bool_t Generate();

  Long64_t NEvents() const { return Long64_t(time_ns_.size()); }
  Double_t Time_ns(Long64_t i) const { return time_ns_[i]; }
  // The beam particle's amplitude, before any reaction.
  Double_t Amplitude(Long64_t i) const { return amplitude_[i]; }
  // What channel 0 (L) or 1 (R) receives from particle i.
  Double_t ChannelAmplitude(Long64_t i, Int_t channel) const {
    return channel == 0 ? amplitudeL_[i] : amplitudeR_[i];
  }
  Bool_t Reacted(Long64_t i) const { return reacted_[i]; }
  const std::vector<Double_t> &Times_ns() const { return time_ns_; }
  const std::vector<Double_t> &Amplitudes() const { return amplitude_; }

private:
  Bool_t LoadReactionPairs(std::vector<Double_t> &pairL,
                           std::vector<Double_t> &pairR,
                           std::vector<Double_t> &cumulative) const;

  const Config &cfg_;
  TRandom3 rng_;
  TRandom3 reactionRng_;
  std::vector<Double_t> time_ns_;
  std::vector<Double_t> amplitude_;
  std::vector<Double_t> amplitudeL_;
  std::vector<Double_t> amplitudeR_;
  std::vector<Bool_t> reacted_;
};

#endif
