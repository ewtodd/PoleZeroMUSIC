// The dataset configuration PulseHistory reads through Constants::cfg. This
// plays the role of analysis/<dataset>/config/Constants.cpp in the MUSIC
// monorepo: a three-channel map (Grid reference, one long L end, one long R
// end) and the pulse-history settings used on the real 37Cl data, except for
// the minimum number of beam-like events, lowered so short test runs can be
// fitted, and the event window, which is tightened to the trigger holdoff
// because the simulated L and R hits carry exactly the seed's timestamp.
#include "Constants.hpp"
#include <mutex>

// Declared in FileSet.hpp; FileSet.cpp is not compiled into this project.
std::mutex g_plot_mutex;
std::mutex g_log_mutex;

namespace Constants {

static DatasetConfig gInstance;

const DatasetConfig &cfg = gInstance;

void InitSimDatasetConfig() {
  gInstance.USE_SOLARIS_DATA = kFALSE;
  gInstance.N_BOARDS = 1;
  gInstance.N_CHANNELS = 3;

  // SAVE_FULL_PLOTS in the pinned tooling (renamed from SAVE_PLOTS); the sim
  // writes few subfiles, so sampling has nothing to skip.
  gInstance.SAVE_FULL_PLOTS = kTRUE;
  gInstance.HAS_CATHODE = kFALSE;
  gInstance.HAS_GRID = kTRUE;
  gInstance.HAS_STRIP0 = kFALSE;
  gInstance.HAS_STRIP17 = kFALSE;

  gInstance.REFERENCE_CHANNEL = "Grid";
  gInstance.REFERENCE_CHANNEL_MIN_ADC = 1000;
  gInstance.REFERENCE_CHANNEL_MAX_ADC = 3000;
  gInstance.EVENT_TIME_WINDOW_US = 0.5;

  gInstance.PULSE_HISTORY_CORRECTION = kTRUE;
  gInstance.PULSE_HISTORY_MIN_EVENTS = 1000;
  // The pinned tooling ships every group off; the simulated long ends get the
  // binned kernel (one coefficient per log dt bin), as in production.
  gInstance.PULSE_HISTORY_GROUPS.long_left.enabled = kTRUE;
  gInstance.PULSE_HISTORY_GROUPS.long_left.kernel = kPulseHistoryBinned;
  gInstance.PULSE_HISTORY_GROUPS.long_right.enabled = kTRUE;
  gInstance.PULSE_HISTORY_GROUPS.long_right.kernel = kPulseHistoryBinned;

  gInstance.STRIP_E_MAX_ADC = 4096.0;

  gInstance.channelMap = {{{0, 0}, "Grid"}, {{0, 1}, "L1"}, {{0, 2}, "R2"}};
}

// Amplitude bands of the kernel (1 here; 4 in the 37Cl configuration), for
// studies that compare the two.
void SetSimPulseHistoryAmpBins(Int_t n) {
  gInstance.PULSE_HISTORY_AMP_BINS = n;
}

struct InitGuard {
  InitGuard() { InitSimDatasetConfig(); }
};
static InitGuard gInit;

} // namespace Constants
