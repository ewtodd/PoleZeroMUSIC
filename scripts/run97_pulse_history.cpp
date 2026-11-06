// PulseHistory on one 37Cl data file with the 37Cl configuration, saving the
// pulse height against the time since the previous pulse, before and after
// the correction, for the run-97 figure of doc/kernel_note.tex (drawn by
// figures --run97). Built and run by scripts/run97_figure.sh, which
// compiles it against the 37Cl configuration in ~/MUSIC; it is not part of
// the project build.
//   run97_pulse_history <out.root> <label> <sol file>
#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include "HitAdapter.hpp"
#include "PulseHistory.hpp"
#include <TFile.h>
#include <TProfile.h>
#include <iostream>
#include <mutex>
#include <vector>

// Declared in FileSet.hpp; FileSet.cpp is not linked.
std::mutex g_plot_mutex;
std::mutex g_log_mutex;

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
  if (!PulseHistory::Measure(hits, groups, res, argv[2])) {
    std::cerr << "PulseHistory::Measure failed" << std::endl;
    return 1;
  }
  PulseHistory::Apply(hits, groups, res);
  std::cout << PulseHistory::Report(res, argv[2]);

  // Deviation from the group mean against log10 of the time since the
  // previous pulse in us, per corrected group.
  TFile out(argv[1], "RECREATE");
  const Int_t group[3] = {PulseHistory::kLONG_RIGHT, PulseHistory::kLONG_LEFT,
                          PulseHistory::kSHORT_LEFT};
  const char *name[3] = {"long_right", "long_left", "short_left"};
  for (Int_t i = 0; i < 3; i++) {
    TProfile *before = res.dtprev_before[group[i]];
    TProfile *after = res.dtprev_after[group[i]];
    if (!before || !after)
      continue;
    before->Write(Form("before_%s", name[i]));
    after->Write(Form("after_%s", name[i]));
  }
  out.Close();
  return 0;
}
