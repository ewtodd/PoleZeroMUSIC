// The reaction events MUSIC's 37Cl analysis tagged and its pure-beam sample,
// with and without the pole-zero correction, for the paper's reaction
// figure. Everything comes from MUSIC's own outputs:
//   - the production strip-sum-scatter cache: the tagged events and MUSIC's
//     pure-beam traces (beam_flat), each with its run and grid seed timestamp
//     (cache v46) and its calibrated deposits (TraceEvt);
//   - the events files of scripts/nocorrection.sh, MUSIC's pipeline with the
//     correction off, read with MUSIC's EnergyView (its calibration), for the
//     same events' deposits without the correction, joined on run and seed.
//
// Writes, to <out.root>, the TTree events: run, seed_ts, reac (the tagged
// reaction strip; -1 for the beam), raw_found (the event is complete in the
// uncorrected build too), raw[18] and corrected[18], each strip's calibrated
// deposit.
// Built and run by scripts/reaction_events.sh, which links the local MUSIC
// tooling; it is not part of the project build.
//   reaction_events <production cache> <no-correction results dir> <out.root>
#include "Constants.hpp"
#include "FileSet.hpp"
#include "Normalization.hpp"
#include "StripSumScatter.hpp"
#include <TChain.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <vector>

namespace {

const Int_t kNStrips = 18;

// One event of the cache, with its deposits without the correction once
// found.
struct Event {
  Int_t run = -1;
  ULong64_t seed_ts = 0;
  Int_t reac = -1;
  Int_t rawFound = 0;
  Float_t raw[kNStrips] = {0.0f};
  Float_t corrected[kNStrips] = {0.0f};
};

// The tagged events and the pure beam of a strip-sum-scatter cache.
Bool_t ReadCache(const TString &path, std::vector<Event> &events) {
  TFile in(path, "READ");
  TTree *traces = in.IsZombie() ? nullptr : (TTree *)in.Get("traces");
  if (!traces || !traces->GetBranch("run")) {
    std::cerr << path << ": no strip-sum-scatter cache with run numbers "
              << "(MUSIC cache v46)" << std::endl;
    return kFALSE;
  }
  const Int_t reacMin =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MIN;
  const Int_t reacMax =
      Constants::cfg.STRIP_SUM_SCATTER_CONFIG.REACTION_STRIP_MAX;
  TraceEvt e;
  traces->SetBranchAddress("leftdE", e.leftdE);
  traces->SetBranchAddress("rightdE", e.rightdE);
  traces->SetBranchAddress("strip0dE", &e.strip0dE);
  traces->SetBranchAddress("strip17dE", &e.strip17dE);
  traces->SetBranchAddress("reac_mask", &e.reac_mask);
  traces->SetBranchAddress("beam_flat", &e.beam_flat);
  traces->SetBranchAddress("seed_ts", &e.seed_ts);
  traces->SetBranchAddress("run", &e.run);
  for (Long64_t i = 0; i < traces->GetEntries(); i++) {
    traces->GetEntry(i);
    if (e.reac_mask == 0 && !e.beam_flat)
      continue;
    Event event;
    event.run = e.run;
    event.seed_ts = e.seed_ts;
    // One bit per event, at strip - REACTION_STRIP_MIN.
    for (Int_t k = reacMin; k <= reacMax; k++)
      if (e.reac_mask & (1u << (k - reacMin)))
        event.reac = k;
    for (Int_t s = 0; s < kNStrips; s++)
      event.corrected[s] = Float_t(e.Total(s));
    events.push_back(event);
  }
  return kTRUE;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 4) {
    std::cerr << "usage: reaction_events <production cache> <no-correction "
                 "results dir> <out.root>"
              << std::endl;
    return 1;
  }
  const TString productionCache = argv[1], noCorrection = argv[2],
                outPath = argv[3];
  // The epoch the cross section uses (late, runs 97-137): its file tags and
  // channel map are what the events files carry.
  Constants::ActivateAnalysisEpoch();

  std::vector<Event> events;
  if (!ReadCache(productionCache, events))
    return 1;
  std::map<Int_t, std::map<ULong64_t, size_t>> index;
  Long64_t duplicates = 0;
  for (size_t i = 0; i < events.size(); i++) {
    std::map<ULong64_t, size_t> &run = index[events[i].run];
    if (run.count(events[i].seed_ts))
      duplicates++;
    run[events[i].seed_ts] = i;
  }

  // Their deposits without the correction, from the no-correction events.
  const std::vector<FileSpec> specs = FileSet::BuildRawOrProcessedFileSpecs();
  Long64_t rawDuplicates = 0;
  for (std::map<Int_t, std::map<ULong64_t, size_t>>::const_iterator r =
           index.begin();
       r != index.end(); ++r) {
    TChain chain("events");
    for (size_t f = 0; f < specs.size(); f++)
      if (specs[f].run == r->first) {
        const TString path = noCorrection + "/root_files/" +
                             FileSet::EventsName(specs[f]) + ".root";
        if (!gSystem->AccessPathName(path))
          chain.Add(path);
      }
    EnergyView ev;
    if (chain.GetEntries() == 0 || !ev.Attach(&chain)) {
      std::cerr << "run " << r->first << ": no calibrated events without the "
                << "correction" << std::endl;
      continue;
    }
    // The branches MUSIC's reaction search reads (EnableEventBranches).
    chain.SetBranchStatus("*", 0);
    for (const char *b : {"LeftdE", "RightdE", "Strip0dE", "Strip17dE",
                          "Cathode", "Grid", "SeedTs"})
      chain.SetBranchStatus(b, 1);
    ULong64_t seed = 0;
    chain.SetBranchAddress("SeedTs", &seed);
    Long64_t found = 0;
    for (Long64_t j = 0; j < chain.GetEntries(); j++) {
      chain.GetEntry(j);
      std::map<ULong64_t, size_t>::const_iterator w = r->second.find(seed);
      if (w == r->second.end())
        continue;
      Event &e = events[w->second];
      if (e.rawFound) {
        rawDuplicates++;
        continue;
      }
      ev.Decode();
      for (Int_t s = 0; s < kNStrips; s++)
        e.raw[s] = Float_t(ev.Total(s));
      e.rawFound = 1;
      found++;
    }
    std::cout << Form("run %d: %lld of %zu events without the correction",
                      r->first, found, r->second.size())
              << std::endl;
  }

  TFile out(outPath, "RECREATE");
  Event e;
  TTree tree("events", "MUSIC's tagged reaction events and pure beam, with "
                       "and without the correction");
  tree.Branch("run", &e.run, "run/I");
  tree.Branch("seed_ts", &e.seed_ts, "seed_ts/l");
  tree.Branch("reac", &e.reac, "reac/I");
  tree.Branch("raw_found", &e.rawFound, "raw_found/I");
  tree.Branch("raw", e.raw, Form("raw[%d]/F", kNStrips));
  tree.Branch("corrected", e.corrected, Form("corrected[%d]/F", kNStrips));
  Long64_t nTagged[2] = {0, 0}, nBeam[2] = {0, 0};
  for (size_t i = 0; i < events.size(); i++) {
    e = events[i];
    tree.Fill();
    Long64_t *n = e.reac >= 0 ? nTagged : nBeam;
    n[0]++;
    n[1] += e.rawFound;
  }
  tree.Write();
  out.Close();
  std::cout << Form("found without the correction: %lld of %lld tagged "
                    "events, %lld of %lld beam events",
                    nTagged[1], nTagged[0], nBeam[1], nBeam[0])
            << std::endl;
  if (duplicates > 0 || rawDuplicates > 0)
    std::cout << Form("WARNING: run and seed not unique: %lld in the cache, "
                      "%lld in the uncorrected events",
                      duplicates, rawDuplicates)
              << std::endl;
  return 0;
}
