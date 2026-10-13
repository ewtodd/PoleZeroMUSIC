#ifndef POLEZERO_PULSEHISTORYADAPTER_HPP
#define POLEZERO_PULSEHISTORYADAPTER_HPP

#include "BinaryUtils.hpp"
#include "SimHit.hpp"
#include <vector>

// Bridges the simulated hits to the unmodified PulseHistory code. PulseHistory
// works on a time-ordered RawHit stream and seeds every beam-like event from a
// reference channel ("Grid"), so each simulated event is turned into a Grid
// hit at a fixed energy followed by the L and R hits at the same timestamp,
// on the three-channel map defined in SimDatasetConfig.cpp
// (board 0: channel 0 = Grid, 1 = L1 (long L), 2 = R2 (long R)).
namespace PulseHistoryAdapter {

const UShort_t kBoard = 0;
const UShort_t kGridChannel = 0;
const UShort_t kLeftChannel = 1;
const UShort_t kRightChannel = 2;
// Inside the reference window of the config (1000-3000 ADC), like the real
// Grid peak (~1950 ADC).
const UShort_t kGridEnergy = 2000;

struct Adapted {
  std::vector<RawHit> hits;
  // For every RawHit, the index into the SimHit vector it came from, or -1
  // for the synthetic Grid hits.
  std::vector<Long64_t> source;
};

Adapted Build(const std::vector<SimHit> &simHits);

// Corrected energies, in the SimHit order, after PulseHistory::Apply ran on
// adapted.hits.
std::vector<Double_t> CorrectedEnergies(const Adapted &adapted,
                                        const std::vector<SimHit> &simHits);

} // namespace PulseHistoryAdapter

#endif
