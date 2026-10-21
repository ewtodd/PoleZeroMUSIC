// PulseHistory on one 37Cl data file with the 37Cl configuration, saving its
// diagnostic plots (the run-97 energy-vs-gap figure of doc/kernel_note.tex)
// in the chosen format. Built and run by scripts/run97_figure.sh, which
// compiles it against the 37Cl configuration in ~/MUSIC; it is not part of
// the project build.
//   run97_pulse_history <pdf|png> <plots dir> <label> <sol file>
#include "BinaryUtils.hpp"
#include "Constants.hpp"
#include "HitAdapter.hpp"
#include "InitUtils.hpp"
#include "PulseHistory.hpp"
#include <iostream>
#include <mutex>
#include <vector>

// Declared in FileSet.hpp; FileSet.cpp is not linked.
std::mutex g_plot_mutex;
std::mutex g_log_mutex;

int main(int argc, char **argv) {
  if (argc < 5) {
    std::cerr << "usage: run97_pulse_history <pdf|png> <plots dir> <label> "
                 "<sol file>"
              << std::endl;
    return 1;
  }
  const TString format = argv[1];
  InitUtils::SetROOTPreferences(format == "png" ? PlotSaveFormat::kPNG
                                                : PlotSaveFormat::kPDF,
                                argv[2], ".", kFALSE);
  const RunEpoch *late = nullptr;
  for (const RunEpoch &e : Constants::cfg.EPOCHS)
    if (e.name == "late")
      late = &e;
  if (!late) {
    std::cerr << "no late epoch in the 37Cl configuration" << std::endl;
    return 1;
  }
  Constants::SetActiveEpoch(late);
  Constants::SetPlotsThisFile(kTRUE);

  std::vector<RawHit> hits;
  SOLReader reader;
  reader.SetSkipTraces(kTRUE);
  if (!reader.Open(argv[4])) {
    std::cerr << "cannot open " << argv[4] << std::endl;
    return 1;
  }
  hits.reserve(60000000);
  while (reader.ReadEvent())
    hits.push_back(SOLHitToRawHit(reader.GetCurrentEvent()));
  reader.Close();

  std::vector<Int_t> groups = PulseHistory::BuildGroupMap();
  PulseHistory::Result res;
  if (!PulseHistory::Measure(hits, groups, res, argv[3])) {
    std::cerr << "PulseHistory::Measure failed" << std::endl;
    return 1;
  }
  PulseHistory::Apply(hits, groups, res);
  std::cout << PulseHistory::Report(res, argv[3]);
  PulseHistory::SavePlots(res, argv[3]);
  return 0;
}
