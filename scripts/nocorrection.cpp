// MUSIC's 37Cl pipeline and reaction search with the pole-zero correction
// switched off: the comparison the paper's reaction figures draw against.
// Everything else is MUSIC's own code and configuration, for the late epoch
// (runs 97-137) that the cross section uses. The results land in
// $MUSIC_RESULTS_DIR, which scripts/nocorrection.sh points away from the
// production results.
// Built and run by scripts/nocorrection.sh; it is not part of the project
// build.
#include "Constants.hpp"
#include "Pipeline.hpp"
#include "StripSumScatter.hpp"
#include <iostream>

int main() {
  // The configuration object behind Constants::cfg is not const.
  DatasetConfig &config = const_cast<DatasetConfig &>(Constants::cfg);
  config.PULSE_HISTORY_CORRECTION = kFALSE;
  for (RunEpoch &epoch : config.EPOCHS)
    epoch.enabled = epoch.name == "late";
  // Built events, calibration and summaries, as main_pipeline.cpp runs them.
  Pipeline::Run();
  // The reaction search's scatters, tags and cache, as
  // main_strip_sum_scatter.cpp runs them, without the simulation overlays
  // and the interactive display that Run() adds.
  Constants::ActivateAnalysisEpoch();
  StripSumScatter scatter;
  if (!scatter.Prepare()) {
    std::cerr << "strip-sum-scatter failed" << std::endl;
    return 1;
  }
  return 0;
}
