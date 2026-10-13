#include "EventGenerator.hpp"

EventGenerator::EventGenerator(const Config &cfg) : cfg_(cfg), rng_(cfg.seed) {}

void EventGenerator::Generate() {
  time_ns_.clear();
  amplitude_.clear();
  time_ns_.reserve(cfg_.nEvents);
  amplitude_.reserve(cfg_.nEvents);
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
    time_ns_.push_back(t);
    amplitude_.push_back(amp);
  }
}
