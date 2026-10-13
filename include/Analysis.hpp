#ifndef POLEZERO_ANALYSIS_HPP
#define POLEZERO_ANALYSIS_HPP

#include "PulseHistory.hpp"
#include "SimHit.hpp"
#include "SimOutput.hpp"
#include <TProfile.h>
#include <TString.h>
#include <vector>

// Everything the analysis executable does with one simulated run: spectra,
// energy-vs-dt profiles and their fits, the linearity check, the PulseHistory
// correction (measured on the run itself or applied from another run) and the
// figures. Numbers that go into the summary table are returned in RunSummary.
namespace Analysis {

struct ExpFit {
  Bool_t ok = kFALSE;
  Double_t p0 = 0.0, p1 = 0.0, p1err = 0.0, p2 = 0.0, p2err = 0.0, p3 = 0.0;
  Double_t chi2 = 0.0;
  Int_t ndf = 0;
};

struct LineFit {
  Bool_t ok = kFALSE;
  Double_t slope = 0.0, slopeErr = 0.0, intercept = 0.0, interceptErr = 0.0;
  Double_t chi2 = 0.0;
  Int_t ndf = 0;
};

struct ChannelSummary {
  Long64_t nHits = 0;       // accepted hits
  Long64_t nWindow = 0;     // non-pileup hits inside the beam window
  Double_t meanRaw = 0.0, rmsRaw = 0.0;
  Double_t meanCorr = 0.0, rmsCorr = 0.0;
  Double_t fracDtBelow12 = 0.0;      // over accepted hits
  Double_t fracDtBelow12Truth = 0.0; // over all beam particles
  ExpFit expRaw;            // p0 + p1 exp(-dt / p2) on the raw profile
  LineFit lineRaw;          // straight line on the raw profile
  LineFit lineCorr;         // straight line on the corrected profile
  ExpFit expCorr;           // exponential on the corrected profile
  LineFit linearity57;      // (energy - trueAmp) vs prevTrueAmp, dt 5-7 us
  LineFit linearity812;     // same, dt 8-12 us
};

struct RunSummary {
  TString label;
  SimOutput::RunHeader header;
  ChannelSummary channel[2];
  Bool_t kernelMeasured = kFALSE;
  TString kernelSource;
  Double_t kernelMaxAbs[2] = {0.0, 0.0}; // largest |k| of the L, R kernels
  Double_t phRmsBefore[2] = {0.0, 0.0};  // PulseHistory's own rms numbers
  Double_t phRmsAfter[2] = {0.0, 0.0};
  Double_t phR2[2] = {0.0, 0.0};
};

// Analyse one run. When externalKernel is null the PulseHistory kernel is
// measured on this run (and left in ownKernel for reuse); otherwise the given
// kernel is applied. Figures go to <plots base>/<label>/.
RunSummary AnalyzeRun(const TString &path, const TString &label,
                      PulseHistory::Result *externalKernel,
                      const TString &externalLabel,
                      PulseHistory::Result &ownKernel);

void PrintSummaryTable(const std::vector<RunSummary> &rows);

} // namespace Analysis

#endif
