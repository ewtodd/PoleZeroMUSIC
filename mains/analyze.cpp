#include "Analysis.hpp"
#include "InitUtils.hpp"
#include "PulseHistory.hpp"
#include <TObjArray.h>
#include <TObjString.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

// analyze --main=<40kHz file>,<10kHz file>,<1kHz file> [--control=<file>]
//         [--signcheck=<file>] [--extra=<label>:<file>,...] [--outdir=output]
// The PulseHistory kernel is measured on the first --main file and applied to
// the others; every other run measures its own.

namespace {

struct RunSpec {
  TString label;
  TString path;
};

std::vector<TString> Split(const TString &s, char sep) {
  std::vector<TString> out;
  TObjArray *parts = s.Tokenize(sep);
  for (Int_t i = 0; i < parts->GetEntries(); i++)
    out.push_back(static_cast<TObjString *>(parts->At(i))->GetString());
  delete parts;
  return out;
}

TString LabelFromPath(const TString &path) {
  TString base = gSystem->BaseName(path);
  if (base.EndsWith(".root"))
    base.Remove(base.Length() - 5);
  if (base.BeginsWith("sim_"))
    base.Remove(0, 4);
  return base;
}

} // namespace

int main(int argc, char **argv) {
  std::vector<RunSpec> mainRuns, standalone;
  TString outdir = "output";
  for (Int_t i = 1; i < argc; i++) {
    TString arg = argv[i];
    if (arg.BeginsWith("--main=")) {
      std::vector<TString> files = Split(arg(7, arg.Length() - 7), ',');
      for (size_t j = 0; j < files.size(); j++)
        mainRuns.push_back({LabelFromPath(files[j]), files[j]});
    } else if (arg.BeginsWith("--control=")) {
      TString f = arg(10, arg.Length() - 10);
      standalone.push_back({LabelFromPath(f), f});
    } else if (arg.BeginsWith("--signcheck=")) {
      TString f = arg(12, arg.Length() - 12);
      standalone.push_back({LabelFromPath(f), f});
    } else if (arg.BeginsWith("--extra=")) {
      std::vector<TString> items = Split(arg(8, arg.Length() - 8), ',');
      for (size_t j = 0; j < items.size(); j++) {
        std::vector<TString> lp = Split(items[j], ':');
        if (lp.size() == 2)
          standalone.push_back({lp[0], lp[1]});
        else
          standalone.push_back({LabelFromPath(items[j]), items[j]});
      }
    } else if (arg.BeginsWith("--outdir=")) {
      outdir = arg(9, arg.Length() - 9);
    } else {
      std::cerr << "Unknown argument " << arg << std::endl;
      return 1;
    }
  }
  if (mainRuns.empty() && standalone.empty()) {
    std::cerr << "Nothing to analyse; see the comment at the top of analyze.cpp"
              << std::endl;
    return 1;
  }
  TString absOut = outdir;
  if (!gSystem->IsAbsoluteFileName(absOut))
    absOut = TString(gSystem->WorkingDirectory()) + "/" + outdir;
  // Plots under outdir; input paths are relative to the working directory.
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, absOut,
                                gSystem->WorkingDirectory(), kFALSE);

  std::vector<Analysis::RunSummary> rows;
  PulseHistory::Result kernel40;
  TString kernelLabel;
  Bool_t haveKernel = kFALSE;
  for (size_t i = 0; i < mainRuns.size(); i++) {
    PulseHistory::Result own;
    if (!haveKernel) {
      rows.push_back(Analysis::AnalyzeRun(mainRuns[i].path, mainRuns[i].label,
                                          nullptr, "", kernel40));
      haveKernel = rows.back().kernelMeasured;
      kernelLabel = mainRuns[i].label;
    } else {
      rows.push_back(Analysis::AnalyzeRun(mainRuns[i].path, mainRuns[i].label,
                                          &kernel40, kernelLabel, own));
    }
  }
  for (size_t i = 0; i < standalone.size(); i++) {
    PulseHistory::Result own;
    rows.push_back(Analysis::AnalyzeRun(standalone[i].path,
                                        standalone[i].label, nullptr, "",
                                        own));
  }
  Analysis::PrintSummaryTable(rows);
  return 0;
}
