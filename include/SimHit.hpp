#ifndef POLEZERO_SIMHIT_HPP
#define POLEZERO_SIMHIT_HPP

#include <Rtypes.h>
#include <vector>

// One accepted trigger on one channel, as written to the "hits" tree.
struct SimHit {
  Int_t channel;         // 0 = L, 1 = R
  ULong64_t timestamp;   // ps, the CoMPASS / RawHit convention
  Double_t energy;       // trapezoid height minus frozen baseline, ADC
  Double_t dtPrev;       // us since the previous accepted hit on this channel
  Double_t dtNext;       // us to the next accepted hit on this channel (-1 if none)
  Double_t trueAmp;      // deposited amplitude, ADC
  Double_t prevTrueAmp;  // of the previous accepted hit on this channel
  Bool_t pileup;         // dtPrev or dtNext < trapRise + trapFlat
  Long64_t eventIndex;   // index into the truth list
};

// Raw waveform, trapezoid output and baseline around a pair of pulses.
struct WaveformSnippet {
  Bool_t found = kFALSE;
  Double_t t0_ns = 0.0;       // time of the first pulse of the pair
  Double_t pairDt_us = 0.0;
  std::vector<Double_t> t_us; // relative to t0
  std::vector<Double_t> raw;
  std::vector<Double_t> trap;
  std::vector<Double_t> baseline;
};

#endif
