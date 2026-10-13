#include "PulseHistoryAdapter.hpp"
#include <TMath.h>

namespace PulseHistoryAdapter {

namespace {

UShort_t ToAdc(Double_t energy) {
  Double_t e = TMath::Floor(energy + 0.5);
  if (e < 0.0)
    e = 0.0;
  if (e > 65535.0)
    e = 65535.0;
  return UShort_t(e);
}

RawHit MakeHit(UShort_t channel, UShort_t energy, ULong64_t timestamp) {
  RawHit h;
  h.board = kBoard;
  h.channel = channel;
  h.energy = energy;
  h.timestamp = timestamp;
  h.flags = 0;
  return h;
}

} // namespace

Adapted Build(const std::vector<SimHit> &simHits) {
  Adapted out;
  out.hits.reserve(simHits.size() * 3 / 2 + 1);
  out.source.reserve(simHits.size() * 3 / 2 + 1);
  Long64_t lastEvent = -1;
  for (size_t i = 0; i < simHits.size(); i++) {
    const SimHit &h = simHits[i];
    if (h.eventIndex != lastEvent) {
      out.hits.push_back(MakeHit(kGridChannel, kGridEnergy, h.timestamp));
      out.source.push_back(-1);
      lastEvent = h.eventIndex;
    }
    const UShort_t ch = h.channel == 0 ? kLeftChannel : kRightChannel;
    out.hits.push_back(MakeHit(ch, ToAdc(h.energy), h.timestamp));
    out.source.push_back(Long64_t(i));
  }
  return out;
}

std::vector<Double_t> CorrectedEnergies(const Adapted &adapted,
                                        const std::vector<SimHit> &simHits) {
  std::vector<Double_t> corrected(simHits.size(), 0.0);
  for (size_t i = 0; i < simHits.size(); i++)
    corrected[i] = simHits[i].energy;
  for (size_t j = 0; j < adapted.hits.size(); j++) {
    const Long64_t src = adapted.source[j];
    if (src < 0)
      continue;
    corrected[src] = Double_t(adapted.hits[j].energy);
  }
  return corrected;
}

} // namespace PulseHistoryAdapter
