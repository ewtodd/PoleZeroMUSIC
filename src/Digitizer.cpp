#include "Digitizer.hpp"
#include <TMath.h>

Digitizer::Digitizer(Double_t noiseSigmaADC, Double_t offsetADC,
                     Double_t adcMax, UInt_t seed)
    : noiseSigma_(noiseSigmaADC), offset_(offsetADC), adcMax_(adcMax),
      rng_(seed) {}

Int_t Digitizer::Digitize(Double_t analogue) {
  Double_t v = analogue + offset_;
  if (noiseSigma_ > 0.0)
    v += rng_.Gaus(0.0, noiseSigma_);
  Double_t rounded = TMath::Floor(v + 0.5);
  if (rounded < 0.0)
    rounded = 0.0;
  if (rounded > adcMax_)
    rounded = adcMax_;
  return Int_t(rounded);
}
