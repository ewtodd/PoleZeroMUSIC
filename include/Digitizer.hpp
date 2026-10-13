#ifndef POLEZERO_DIGITIZER_HPP
#define POLEZERO_DIGITIZER_HPP

#include <Rtypes.h>
#include <TRandom3.h>

// Gaussian white noise per sample, DC offset, rounding and clamping to the
// ADC range.
class Digitizer {
public:
  Digitizer(Double_t noiseSigmaADC, Double_t offsetADC, Double_t adcMax,
            UInt_t seed);
  Int_t Digitize(Double_t analogue);

private:
  Double_t noiseSigma_;
  Double_t offset_;
  Double_t adcMax_;
  TRandom3 rng_;
};

#endif
