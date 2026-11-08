#include "EventGenerator.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

EventGenerator::EventGenerator(const Config &cfg)
    : cfg_(cfg), rng_(cfg.seed), reactionRng_(cfg.seed + 7919) {}

// "rL rR weight" lines, # for comments: the two long ends' deposits of a
// reacted particle in units of the beam peak, and the pair's weight.
Bool_t EventGenerator::LoadReactionPairs(
    std::vector<Double_t> &pairL, std::vector<Double_t> &pairR,
    std::vector<Double_t> &cumulative) const {
  std::ifstream in(cfg_.reactionPairs.Data());
  if (!in) {
    std::cerr << "Cannot read reaction pairs " << cfg_.reactionPairs
              << std::endl;
    return kFALSE;
  }
  std::string line;
  Double_t total = 0.0;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream fields(line);
    Double_t left = 0.0, right = 0.0, weight = 0.0;
    if (!(fields >> left >> right >> weight) || !(weight > 0.0))
      continue;
    pairL.push_back(left);
    pairR.push_back(right);
    total += weight;
    cumulative.push_back(total);
  }
  if (cumulative.empty()) {
    std::cerr << "No reaction pairs in " << cfg_.reactionPairs << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t EventGenerator::Generate() {
  time_ns_.clear();
  amplitude_.clear();
  amplitudeL_.clear();
  amplitudeR_.clear();
  reacted_.clear();
  time_ns_.reserve(cfg_.nEvents);
  amplitude_.reserve(cfg_.nEvents);
  amplitudeL_.reserve(cfg_.nEvents);
  amplitudeR_.reserve(cfg_.nEvents);
  reacted_.reserve(cfg_.nEvents);
  std::vector<Double_t> pairL, pairR, cumulative;
  const Bool_t reactions = cfg_.reactionFraction > 0.0;
  if (reactions && !LoadReactionPairs(pairL, pairR, cumulative))
    return kFALSE;
  // The first particle arrives well after the start so that the filters and
  // the baseline are settled (10 times the trapezoid length).
  const Double_t settle_ns =
      10.0 * (2.0 * cfg_.trapRise_ns + cfg_.trapFlat_ns);
  Double_t t = settle_ns;
  const Double_t sigma = cfg_.beamPeakADC * cfg_.beamPeakSigmaFrac;
  for (Long64_t i = 0; i < cfg_.nEvents; i++) {
    t += rng_.Exp(cfg_.MeanDt_ns());
    Double_t amp = rng_.Gaus(cfg_.beamPeakADC, sigma);
    if (amp < 0.0)
      amp = 0.0;
    Double_t ampL = amp, ampR = amp;
    Bool_t reacted = kFALSE;
    if (reactions && reactionRng_.Uniform() < cfg_.reactionFraction) {
      const Double_t u = reactionRng_.Uniform(cumulative.back());
      size_t k = size_t(std::upper_bound(cumulative.begin(), cumulative.end(),
                                         u) -
                        cumulative.begin());
      if (k >= cumulative.size())
        k = cumulative.size() - 1;
      ampL = pairL[k] * cfg_.beamPeakADC;
      ampR = pairR[k] * cfg_.beamPeakADC;
      reacted = kTRUE;
    }
    time_ns_.push_back(t);
    amplitude_.push_back(amp);
    amplitudeL_.push_back(ampL);
    amplitudeR_.push_back(ampR);
    reacted_.push_back(reacted);
  }
  return kTRUE;
}
