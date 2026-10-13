#ifndef POLEZERO_SOLWAVE_HPP
#define POLEZERO_SOLWAVE_HPP

// Shared pieces for the raw-waveform checks on the 88Sr SOL files: the
// channel map, the DPP-PHA settings dump, and locating a pulse in a record.
#include <Rtypes.h>
#include <TString.h>
#include <vector>

namespace solwave {

const Double_t kSample_ns = 8.0;
const Int_t kMaxCh = 64;

// Channel map of the 88Sr run (screenshot next to data_raw): strip s has its
// left end on 4(s-1)+1 and its right end on 4(s-1); odd strips are long on
// the left, even strips long on the right. Strips 0 and 17, grid and cathode
// are single channels. Everything else on the board is unused.
enum Group {
  kLongL = 0,
  kShortL = 1,
  kLongR = 2,
  kShortR = 3,
  kStrip0 = 4,
  kStrip17 = 5,
  kGrid = 6,
  kCathode = 7,
  kNGroups = 8
};
extern const char *kGroupName[kNGroups];
extern const char *kGroupTag[kNGroups];
extern const Int_t kGroupColour[kNGroups];

Int_t GroupOf(Int_t ch);
Int_t StripOf(Int_t ch);

// The filter timing a channel was acquired with, from the settings dump.
struct ChSettings {
  Bool_t enabled = kFALSE;
  Int_t pre = 0;       // ChPreTriggerT, samples
  Int_t rise = 0;      // EnergyFilterRiseTimeT, samples
  Int_t flat = 0;      // EnergyFilterFlatTopT, samples
  Int_t peakDelay = 0; // PeakingPosition % of the flat top, samples
  Int_t avg = 1;       // EnergyFilterPeakingAvg, samples
  Double_t pz_us = 0.0;
  TString probe1;
  TString baseline;
};

// Lines look like "/ch/12/par/EnergyFilterRiseTimeT   !2! 1234!704".
std::vector<ChSettings> ReadSettings(const TString &path);

Double_t Mean(const Int_t *v, Int_t from, Int_t n);
Double_t MeanD(const std::vector<Double_t> &v, Int_t from, Int_t n);

// Locate the input pulse in a record. Polarity from the sign of the largest
// excursion over the whole record relative to the first samples; the pulse
// start is the 10 % crossing walking back from the extreme, and the
// baseline and its RMS are measured just before that start. Returns kFALSE
// when the pulse begins too early in the record to have a baseline.
Bool_t LocatePulse(const Int_t *in, Int_t len, Double_t &pre, Double_t &pol,
                   Double_t &amp, Double_t &noise, Int_t &iStart, Int_t &i90,
                   Int_t &iExt);

} // namespace solwave

#endif
