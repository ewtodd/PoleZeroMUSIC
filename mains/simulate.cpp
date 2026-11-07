#include "ChannelSimulator.hpp"
#include "Config.hpp"
#include "EventGenerator.hpp"
#include "IOUtils.hpp"
#include "SimOutput.hpp"
#include <TROOT.h>
#include <TSystem.h>
#include <chrono>
#include <iostream>
#include <thread>

namespace {

void RunChannel(ChannelSimulator *sim) { sim->Run(); }

void ReportChannel(const char *name, const ChannelSimulator &sim) {
  std::cout << "  " << name << ": " << sim.Hits().size() << " hits, "
            << sim.Dropped() << " dropped by trigger holdoff, "
            << sim.SamplesGenerated() << " samples generated, "
            << sim.SamplesSkipped() << " skipped analytically, snippet "
            << (sim.Snippet().found ? "saved" : "not found") << std::endl;
}

} // namespace

int main(int argc, char **argv) {
  Config cfg;
  if (!cfg.ParseCommandLine(argc, argv)) {
    std::cerr << "usage: simulate [--name=value ...]" << std::endl;
    return 1;
  }
  gROOT->SetBatch(kTRUE);
  IO::SetRootFilesBaseDir(gSystem->WorkingDirectory());
  cfg.Print(std::cout);

  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  EventGenerator events(cfg);
  events.Generate();
  std::cout << "Generated " << events.NEvents() << " beam particles over "
            << (events.Time_ns(events.NEvents() - 1) - events.Time_ns(0)) *
                   1.0e-9
            << " s" << std::endl;

  ChannelSimulator left(cfg, 0, cfg.TauL_ns(), cfg.noiseSigmaADC_L,
                        cfg.seed + 1001, events);
  ChannelSimulator right(cfg, 1, cfg.TauR_ns(), cfg.noiseSigmaADC_R,
                         cfg.seed + 1002, events);
  std::thread leftThread(RunChannel, &left);
  std::thread rightThread(RunChannel, &right);
  leftThread.join();
  rightThread.join();
  ReportChannel("L", left);
  ReportChannel("R", right);

  std::vector<SimHit> hits = SimOutput::MergeHits(left.Hits(), right.Hits());
  if (!SimOutput::Write(cfg, events, hits, left.Snippet(), right.Snippet()))
    return 1;
  const Double_t seconds =
      std::chrono::duration<Double_t>(std::chrono::steady_clock::now() - t0)
          .count();
  std::cout << "Wrote " << hits.size() << " hits to " << cfg.outputFile
            << " in " << seconds << " s" << std::endl;
  return 0;
}
