// PulseHistory on one 37Cl data file with the 37Cl configuration. Saves what
// the run-97 figures and tables of the paper take (figures --run97): the group
// kernels it fits, the pulse height against the time since the previous pulse
// before and after the correction, the strip spectra and their beam-peak FWHM
// before and after it, and the long groups' spread before correction and
// after it with one and with four amplitude bands. Built and run by
// scripts/run97_figure.sh, which compiles it against the 37Cl configuration
// in ~/MUSIC; it is not part of the project build.
//   run97_pulse_history <out.root> <label> <sol file>
#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include "HitAdapter.hpp"
#include "PulseHistory.hpp"
#include <TFile.h>
#include <TH1D.h>
#include <TMath.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TVectorD.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <mutex>
#include <vector>

// Declared in FileSet.hpp; FileSet.cpp is not linked.
std::mutex g_plot_mutex;
std::mutex g_log_mutex;

namespace {

const Int_t kNStrips = 18;

Int_t IndexOf(const RawHit &hit) {
  return Int_t(hit.board) * Constants::ActiveNChannels() + Int_t(hit.channel);
}

// A channel's index as PulseHistory numbers it, from its name; -1 if absent.
Int_t IndexOfName(const TString &name) {
  const Int_t nch = Constants::ActiveNChannels();
  const std::map<std::pair<Int_t, Int_t>, TString> &map =
      Constants::ActiveChannelMap();
  for (std::map<std::pair<Int_t, Int_t>, TString>::const_iterator it =
           map.begin();
       it != map.end(); ++it)
    if (it->second == name)
      return it->first.first * nch + it->first.second;
  return -1;
}

// The DeltaE of every strip in every event, as MUSIC's events summary has it:
// events are seeded on the reference channel, each end of a strip contributes
// its largest height within the event window, and strips 1-16 sum their two
// ends while strips 0 and 17 are whole. Strips that did not fire are left out.
void StripTotals(const std::vector<RawHit> &hits,
                 const std::vector<Double_t> &energy,
                 std::vector<std::vector<Float_t>> &totals) {
  std::vector<std::vector<Int_t>> ends(kNStrips);
  ends[0].push_back(IndexOfName("Strip0"));
  ends[kNStrips - 1].push_back(IndexOfName("Strip17"));
  for (Int_t s = 1; s < kNStrips - 1; s++) {
    ends[s].push_back(IndexOfName(Form("L%d", s)));
    ends[s].push_back(IndexOfName(Form("R%d", s)));
  }
  Int_t maxIndex = 0;
  for (Int_t s = 0; s < kNStrips; s++)
    for (size_t e = 0; e < ends[s].size(); e++)
      maxIndex = std::max(maxIndex, ends[s][e]);
  // Each end's slot in the per-event array of largest heights.
  std::vector<Int_t> slotOf(maxIndex + 1, -1);
  Int_t nSlots = 0;
  for (Int_t s = 0; s < kNStrips; s++)
    for (size_t e = 0; e < ends[s].size(); e++)
      if (ends[s][e] >= 0)
        slotOf[ends[s][e]] = nSlots++;

  const Int_t reference = IndexOfName(Constants::ActiveReferenceChannel());
  const Double_t window_ps = Constants::ActiveEventTimeWindowUs() * 1.0e6;
  std::vector<Double_t> largest(nSlots, 0.0);
  totals.assign(kNStrips, std::vector<Float_t>());
  for (size_t j = 0; j < hits.size(); j++) {
    if (IndexOf(hits[j]) != reference)
      continue;
    std::fill(largest.begin(), largest.end(), 0.0);
    for (size_t k = j + 1;
         k < hits.size() &&
         Double_t(hits[k].timestamp - hits[j].timestamp) < window_ps;
         k++) {
      const Int_t i = IndexOf(hits[k]);
      if (i < 0 || i > maxIndex || slotOf[i] < 0)
        continue;
      largest[slotOf[i]] = std::max(largest[slotOf[i]], energy[k]);
    }
    for (Int_t s = 0; s < kNStrips; s++) {
      Double_t total = 0.0;
      for (size_t e = 0; e < ends[s].size(); e++)
        if (ends[s][e] >= 0)
          total += largest[slotOf[ends[s][e]]];
      if (total > 0.0)
        totals[s].push_back(Float_t(total));
    }
  }
}

// A strip's beam peak: the most populated bin above a sixteenth of the range
// that holds 99.9 % of its events, then refined to within 20 % of it.
Double_t BeamPeak(std::vector<Float_t> values) {
  if (values.empty())
    return 0.0;
  std::sort(values.begin(), values.end());
  const Double_t top = values[size_t(0.999 * Double_t(values.size() - 1))];
  Double_t peak = 0.0;
  Double_t lo = top / 16.0, hi = top;
  for (Int_t pass = 0; pass < 2; pass++) {
    TH1D histogram("beam_peak", "", 128, lo, hi);
    histogram.SetDirectory(nullptr);
    for (size_t i = 0; i < values.size(); i++)
      histogram.Fill(values[i]);
    peak = histogram.GetBinCenter(histogram.GetMaximumBin());
    lo = 0.8 * peak;
    hi = 1.2 * peak;
  }
  return peak;
}

// A strip's beam-peak FWHM in percent of the peak position, from its
// per-event values rather than the spectra's 0.02 a.u. bins, which are only a
// few per peak width: binned at 0.2 % of the peak over 0.5-1.5 of it,
// smoothed over five bins (1 %), with the half-maximum crossings interpolated
// linearly. 0 when no peak is found.
Double_t FwhmPercent(const std::vector<Float_t> &values, Double_t peak) {
  if (!(peak > 0.0) || values.empty())
    return 0.0;
  const Double_t lo = 0.5, width = 0.002;
  const Int_t nBins = 500, half = 2;
  std::vector<Double_t> counts(nBins, 0.0);
  for (size_t i = 0; i < values.size(); i++) {
    const Int_t b = Int_t(TMath::Floor((values[i] / peak - lo) / width));
    if (b >= 0 && b < nBins)
      counts[b]++;
  }
  std::vector<Double_t> smooth(nBins, 0.0);
  Int_t top = half;
  for (Int_t b = half; b < nBins - half; b++) {
    for (Int_t k = -half; k <= half; k++)
      smooth[b] += counts[b + k] / (2 * half + 1);
    if (smooth[b] > smooth[top])
      top = b;
  }
  const Double_t halfMax = 0.5 * smooth[top];
  Int_t left = top, right = top;
  while (left > half && smooth[left] > halfMax)
    left--;
  while (right < nBins - half - 1 && smooth[right] > halfMax)
    right++;
  if (smooth[left] > halfMax || smooth[right] > halfMax)
    return 0.0;
  // Bin centres in units of the peak; each crossing lies between a bin below
  // half maximum and its neighbour towards the top.
  const Double_t xLeft =
      lo + (left + 0.5) * width +
      width * (halfMax - smooth[left]) / (smooth[left + 1] - smooth[left]);
  const Double_t xRight =
      lo + (right + 0.5) * width -
      width * (halfMax - smooth[right]) / (smooth[right - 1] - smooth[right]);
  return 100.0 * (xRight - xLeft) / (lo + (top + 0.5) * width);
}

// The strip spectra in units of each strip's own beam peak, as the
// calibration anchors them, and each strip's beam-peak FWHM in percent.
void WriteStripSpectra(const std::vector<RawHit> &hits,
                       const std::vector<Double_t> &energy, const char *which) {
  std::vector<std::vector<Float_t>> totals;
  StripTotals(hits, energy, totals);
  TVectorD fwhm(kNStrips);
  for (Int_t s = 0; s < kNStrips; s++) {
    const Double_t peak = BeamPeak(totals[s]);
    TH1D *spectrum = new TH1D(Form("strip_%s_S%d", which, s),
                              ";#DeltaE [a.u.];Counts", 200, 0.0, 4.0);
    if (peak > 0.0)
      for (size_t i = 0; i < totals[s].size(); i++)
        spectrum->Fill(totals[s][i] / peak);
    spectrum->Write();
    fwhm[s] = FwhmPercent(totals[s], peak);
    std::cout << Form("  S%d %s: %zu events, beam peak %.0f ADC, FWHM %.2f %%",
                      s, which, totals[s].size(), peak, fwhm[s])
              << std::endl;
  }
  fwhm.Write(Form("fwhm_%s", which));
}

// The mean trigger rate of the long ends' channels, in kHz: each channel's
// hits over the time between its first and last hit.
Double_t LongRateKHz(const std::vector<RawHit> &hits,
                     const std::vector<Int_t> &groups) {
  const Int_t nIndex = Int_t(groups.size());
  std::vector<Double_t> count(nIndex, 0.0);
  std::vector<ULong64_t> first(nIndex, 0), last(nIndex, 0);
  for (size_t j = 0; j < hits.size(); j++) {
    const Int_t i = IndexOf(hits[j]);
    if (i < 0 || i >= nIndex ||
        (groups[i] != PulseHistory::kLONG_LEFT &&
         groups[i] != PulseHistory::kLONG_RIGHT))
      continue;
    if (count[i] == 0.0)
      first[i] = hits[j].timestamp;
    last[i] = hits[j].timestamp;
    count[i]++;
  }
  Double_t sum = 0.0;
  Int_t channels = 0;
  for (Int_t i = 0; i < nIndex; i++)
    if (count[i] > 1.0 && last[i] > first[i]) {
      sum += count[i] / (Double_t(last[i] - first[i]) * 1.0e-12) / 1.0e3;
      channels++;
    }
  return channels > 0 ? sum / channels : 0.0;
}

// The spread of each corrected group, in percent of each channel's beam peak,
// for every set of energies: hits of the group's channels with no other hit
// on the same channel within the trapezoid's 5 us, and within 0.5-1.5 times
// the channel's beam peak in that set, as figures takes the simulated spread.
// Returns [group][set] for the groups given.
std::vector<std::vector<Double_t>>
GroupSpreads(const std::vector<RawHit> &hits, const std::vector<Int_t> &groups,
             const std::vector<Double_t> &mode,
             const std::vector<const std::vector<Double_t> *> &sets,
             const std::vector<Int_t> &wanted) {
  const ULong64_t isolation_ps = 5000000;
  const Int_t nIndex = Int_t(TMath::Min(groups.size(), mode.size()));
  // Whether each hit's previous and next hits on its channel are far enough;
  // the hits are in time order.
  std::vector<Bool_t> alone(hits.size(), kTRUE);
  std::vector<Long64_t> last(nIndex, -1);
  for (size_t j = 0; j < hits.size(); j++) {
    const Int_t i = IndexOf(hits[j]);
    if (i < 0 || i >= nIndex)
      continue;
    if (last[i] >= 0 &&
        hits[j].timestamp - hits[last[i]].timestamp < isolation_ps) {
      alone[j] = kFALSE;
      alone[last[i]] = kFALSE;
    }
    last[i] = Long64_t(j);
  }
  std::vector<std::vector<Double_t>> spreads(
      wanted.size(), std::vector<Double_t>(sets.size(), 0.0));
  for (size_t w = 0; w < wanted.size(); w++)
    for (size_t k = 0; k < sets.size(); k++) {
      Double_t n = 0.0, sum = 0.0, sum2 = 0.0;
      for (size_t j = 0; j < hits.size(); j++) {
        const Int_t i = IndexOf(hits[j]);
        if (!alone[j] || i < 0 || i >= nIndex || groups[i] != wanted[w] ||
            !(mode[i] > 0.0))
          continue;
        const Double_t x = (*sets[k])[j] / mode[i];
        if (x <= 0.5 || x >= 1.5)
          continue;
        n++;
        sum += x;
        sum2 += x * x;
      }
      const Double_t mean = n > 0.0 ? sum / n : 0.0;
      spreads[w][k] =
          n > 0.0 ? 100.0 * TMath::Sqrt(TMath::Max(0.0, sum2 / n - mean * mean))
                  : 0.0;
    }
  return spreads;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: run97_pulse_history <out.root> <label> <sol file>"
              << std::endl;
    return 1;
  }
  const RunEpoch *late = nullptr;
  for (const RunEpoch &e : Constants::cfg.EPOCHS)
    if (e.name == "late")
      late = &e;
  if (!late) {
    std::cerr << "no late epoch in the 37Cl configuration" << std::endl;
    return 1;
  }
  Constants::SetActiveEpoch(late);
  // The diagnostic profiles are only filled for a file that draws.
  Constants::SetPlotsThisFile(kTRUE);

  std::vector<RawHit> hits;
  SOLReader reader;
  reader.SetSkipTraces(kTRUE);
  if (!reader.Open(argv[3])) {
    std::cerr << "cannot open " << argv[3] << std::endl;
    return 1;
  }
  hits.reserve(60000000);
  while (reader.ReadEvent())
    hits.push_back(SOLHitToRawHit(reader.GetCurrentEvent()));
  reader.Close();

  std::vector<Int_t> groups = PulseHistory::BuildGroupMap();
  PulseHistory::Result res;
  // Measure puts the hits in time order; Apply then corrects them in place.
  if (!PulseHistory::Measure(hits, groups, res, argv[2])) {
    std::cerr << "PulseHistory::Measure failed" << std::endl;
    return 1;
  }
  std::vector<Double_t> raw(hits.size());
  for (size_t j = 0; j < hits.size(); j++)
    raw[j] = hits[j].energy;
  PulseHistory::Apply(hits, groups, res);
  std::vector<Double_t> corrected(hits.size());
  for (size_t j = 0; j < hits.size(); j++)
    corrected[j] = hits[j].energy;
  std::cout << PulseHistory::Report(res, argv[2]);

  // The same correction with four amplitude bands, kernels per band of the
  // previous pulse height (single, double, triple beam and above), for the
  // band table; without diagnostics, which the figures take from the
  // production pass. Measure leaves time-ordered hits in place, so the raw
  // heights restore by index; strip totals use only the hits' channels and
  // times. The configuration object behind Constants::cfg is not const, and
  // PulseHistory reads the band count from it on every Measure.
  const Int_t nBands = 4;
  DatasetConfig &config = const_cast<DatasetConfig &>(Constants::cfg);
  const Int_t productionBands = config.PULSE_HISTORY_AMP_BINS;
  config.PULSE_HISTORY_AMP_BINS = nBands;
  Constants::SetPlotsThisFile(kFALSE);
  for (size_t j = 0; j < hits.size(); j++) {
    hits[j].energy = UShort_t(raw[j]);
    hits[j].flags &= ~PulseHistory::kFlagClamped;
  }
  std::vector<Double_t> correctedBands(hits.size());
  {
    PulseHistory::Result banded;
    if (!PulseHistory::Measure(hits, groups, banded, argv[2])) {
      std::cerr << "PulseHistory::Measure failed with " << nBands << " bands"
                << std::endl;
      return 1;
    }
    PulseHistory::Apply(hits, groups, banded);
  }
  config.PULSE_HISTORY_AMP_BINS = productionBands;
  for (size_t j = 0; j < hits.size(); j++)
    correctedBands[j] = hits[j].energy;

  TFile out(argv[1], "RECREATE");
  // The group kernels, shift per unit previous pulse height in each dt bin;
  // in the order the pipeline log lists them.
  const Int_t nGroups = 3;
  const Int_t group[nGroups] = {PulseHistory::kLONG_LEFT,
                                PulseHistory::kLONG_RIGHT,
                                PulseHistory::kSHORT_LEFT};
  const char *name[nGroups] = {"long_left", "long_right", "short_left"};
  for (Int_t i = 0; i < nGroups; i++) {
    const PulseHistory::Kernel &kernel = res.kernel[group[i]];
    if (!kernel.ok)
      continue;
    TVectorD values(PulseHistory::kNBins);
    for (Int_t b = 0; b < PulseHistory::kNBins; b++)
      values[b] = kernel.k[0][b];
    values.Write(Form("kernel_%s", name[i]));
  }
  // Deviation from the group mean against log10 of the time since the
  // previous pulse in us, per corrected group.
  for (Int_t i = 0; i < nGroups; i++) {
    TProfile *before = res.dtprev_before[group[i]];
    TProfile *after = res.dtprev_after[group[i]];
    if (!before || !after)
      continue;
    before->Write(Form("before_%s", name[i]));
    after->Write(Form("after_%s", name[i]));
  }
  WriteStripSpectra(hits, raw, "raw");
  WriteStripSpectra(hits, corrected, "corrected");
  // Which chain reads each strip's long segment: 1 the right, 0 the left, -1
  // for the unsegmented strips 0 and 17, from PulseHistory's group map.
  TVectorD longRight(kNStrips);
  for (Int_t s = 0; s < kNStrips; s++) {
    const Int_t right = IndexOfName(Form("R%d", s));
    longRight[s] =
        (s == 0 || s == kNStrips - 1 || right < 0)
            ? -1.0
            : (groups[right] == PulseHistory::kLONG_RIGHT ? 1.0 : 0.0);
  }
  longRight.Write("strip_long_right");

  // Spread of the long groups before correction, with the production single
  // band, and with four bands, in percent of each channel's beam peak.
  std::vector<const std::vector<Double_t> *> sets;
  sets.push_back(&raw);
  sets.push_back(&corrected);
  sets.push_back(&correctedBands);
  std::vector<Int_t> wanted;
  wanted.push_back(PulseHistory::kLONG_LEFT);
  wanted.push_back(PulseHistory::kLONG_RIGHT);
  const std::vector<std::vector<Double_t>> spreads =
      GroupSpreads(hits, groups, res.mode, sets, wanted);
  for (size_t w = 0; w < wanted.size(); w++) {
    TVectorD values(Int_t(sets.size()));
    for (size_t k = 0; k < sets.size(); k++)
      values[Int_t(k)] = spreads[w][k];
    values.Write(Form("band_spread_%s", w == 0 ? "long_left" : "long_right"));
    std::cout << Form("  %s spread: raw %.2f, 1 band %.2f, %d bands %.2f %% of "
                      "the beam peak",
                      w == 0 ? "long left" : "long right", spreads[w][0],
                      spreads[w][1], nBands, spreads[w][2])
              << std::endl;
  }
  TParameter<Int_t>("band_count", nBands).Write();
  const Double_t rate = LongRateKHz(hits, groups);
  TParameter<Double_t>("long_rate_kHz", rate).Write();
  std::cout << Form("  long-end rate %.1f kHz per channel", rate) << std::endl;
  out.Close();
  return 0;
}
