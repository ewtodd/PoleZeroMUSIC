#ifndef POLEZERO_SIMOUTPUT_HPP
#define POLEZERO_SIMOUTPUT_HPP

#include "Config.hpp"
#include "EventGenerator.hpp"
#include "SimHit.hpp"
#include <TFile.h>
#include <vector>

// ROOT output of one simulated run: the "hits" tree (time ordered, L and R
// merged), the "truth" tree, the configuration as TParameters and the
// waveform snippets as TGraphs. Reading back is provided for the analysis.
namespace SimOutput {

struct RunHeader {
  Double_t rate_kHz = 0.0;
  Double_t tauL_us = 0.0;
  Double_t tauR_us = 0.0;
  Double_t trapDecay_us = 0.0;
  Double_t trapRise_ns = 0.0;
  Double_t trapFlat_ns = 0.0;
  Double_t beamPeakADC = 0.0;
  Double_t beamPeakSigmaFrac = 0.0;
  Double_t samplePeriod_ns = 0.0;
  Double_t adcOffset = 0.0;
  Int_t baselineMeanSamples = 0;
  Double_t baselineHoldoff_ns = 0.0;
  Long64_t nEvents = 0;
};

// Merge the two channels in time order.
std::vector<SimHit> MergeHits(const std::vector<SimHit> &left,
                              const std::vector<SimHit> &right);

Bool_t Write(const Config &cfg, const EventGenerator &events,
             const std::vector<SimHit> &hits, const WaveformSnippet &left,
             const WaveformSnippet &right);

Bool_t Read(const TString &path, RunHeader &header, std::vector<SimHit> &hits);

} // namespace SimOutput

#endif
