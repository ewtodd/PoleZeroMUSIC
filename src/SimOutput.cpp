#include "SimOutput.hpp"
#include "IOUtils.hpp"
#include <TGraph.h>
#include <TParameter.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>

namespace SimOutput {

namespace {

Bool_t HitBefore(const SimHit &a, const SimHit &b) {
  if (a.timestamp != b.timestamp)
    return a.timestamp < b.timestamp;
  return a.channel < b.channel;
}

void WriteSnippet(const WaveformSnippet &s, const char *tag) {
  if (!s.found)
    return;
  const Int_t n = Int_t(s.t_us.size());
  TGraph *raw = new TGraph(n, &s.t_us[0], &s.raw[0]);
  raw->SetName(Form("raw_%s", tag));
  raw->SetTitle(";t #minus t_{1} [#mus];ADC");
  raw->Write();
  TGraph *trap = new TGraph(n, &s.t_us[0], &s.trap[0]);
  trap->SetName(Form("trap_%s", tag));
  trap->SetTitle(";t #minus t_{1} [#mus];Trapezoid [ADC]");
  trap->Write();
  TGraph *base = new TGraph(n, &s.t_us[0], &s.baseline[0]);
  base->SetName(Form("baseline_%s", tag));
  base->SetTitle(";t #minus t_{1} [#mus];Baseline [ADC]");
  base->Write();
  TParameter<Double_t>(Form("snippet_t0_ns_%s", tag), s.t0_ns).Write();
  TParameter<Double_t>(Form("snippet_dt_us_%s", tag), s.pairDt_us).Write();
  delete raw;
  delete trap;
  delete base;
}

void WriteParameter(const char *name, Double_t value) {
  TParameter<Double_t>(name, value).Write();
}

Double_t ReadParameter(TFile *f, const char *name) {
  TParameter<Double_t> *p = dynamic_cast<TParameter<Double_t> *>(f->Get(name));
  return p ? p->GetVal() : 0.0;
}

} // namespace

std::vector<SimHit> MergeHits(const std::vector<SimHit> &left,
                              const std::vector<SimHit> &right) {
  std::vector<SimHit> all;
  all.reserve(left.size() + right.size());
  all.insert(all.end(), left.begin(), left.end());
  all.insert(all.end(), right.begin(), right.end());
  std::stable_sort(all.begin(), all.end(), HitBefore);
  return all;
}

Bool_t Write(const Config &cfg, const EventGenerator &events,
             const std::vector<SimHit> &hits, const WaveformSnippet &left,
             const WaveformSnippet &right) {
  TFile *f = IO::OpenForWriting(cfg.outputFile);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << cfg.outputFile << " for writing"
              << std::endl;
    return kFALSE;
  }
  f->cd();

  TTree *tree = new TTree("hits", "Simulated DPP-PHA hits");
  Int_t channel = 0;
  ULong64_t timestamp = 0;
  Double_t energy = 0.0, dtPrev = 0.0, dtNext = 0.0, trueAmp = 0.0,
           prevTrueAmp = 0.0;
  Bool_t pileup = kFALSE;
  Long64_t eventIndex = 0;
  tree->Branch("channel", &channel, "channel/I");
  tree->Branch("timestamp", &timestamp, "timestamp/l");
  tree->Branch("energy", &energy, "energy/D");
  tree->Branch("dtPrev", &dtPrev, "dtPrev/D");
  tree->Branch("dtNext", &dtNext, "dtNext/D");
  tree->Branch("trueAmp", &trueAmp, "trueAmp/D");
  tree->Branch("prevTrueAmp", &prevTrueAmp, "prevTrueAmp/D");
  tree->Branch("pileup", &pileup, "pileup/O");
  tree->Branch("eventIndex", &eventIndex, "eventIndex/L");
  for (size_t i = 0; i < hits.size(); i++) {
    const SimHit &h = hits[i];
    channel = h.channel;
    timestamp = h.timestamp;
    energy = h.energy;
    dtPrev = h.dtPrev;
    dtNext = h.dtNext;
    trueAmp = h.trueAmp;
    prevTrueAmp = h.prevTrueAmp;
    pileup = h.pileup;
    eventIndex = h.eventIndex;
    tree->Fill();
  }
  tree->Write();

  TTree *truth = new TTree("truth", "True beam particles");
  Double_t time_ns = 0.0, amp = 0.0;
  truth->Branch("time_ns", &time_ns, "time_ns/D");
  truth->Branch("amp", &amp, "amp/D");
  for (Long64_t i = 0; i < events.NEvents(); i++) {
    time_ns = events.Time_ns(i);
    amp = events.Amplitude(i);
    truth->Fill();
  }
  truth->Write();

  WriteParameter("rate_kHz", cfg.rate_kHz);
  WriteParameter("tauL_us", cfg.tauL_us);
  WriteParameter("tauR_us", cfg.tauR_us);
  WriteParameter("trapDecay_us", cfg.trapDecay_us);
  WriteParameter("trapRise_ns", cfg.trapRise_ns);
  WriteParameter("trapFlat_ns", cfg.trapFlat_ns);
  WriteParameter("beamPeakADC", cfg.beamPeakADC);
  WriteParameter("beamPeakSigmaFrac", cfg.beamPeakSigmaFrac);
  WriteParameter("samplePeriod_ns", cfg.samplePeriod_ns);
  WriteParameter("adcOffset", cfg.adcOffset);
  WriteParameter("baselineMeanSamples", Double_t(cfg.baselineMeanSamples));
  WriteParameter("baselineHoldoff_ns", cfg.baselineHoldoff_ns);
  WriteParameter("nEvents", Double_t(cfg.nEvents));

  WriteSnippet(left, "L");
  WriteSnippet(right, "R");

  f->Close();
  delete f;
  return kTRUE;
}

Bool_t Read(const TString &path, RunHeader &header,
            std::vector<SimHit> &hits) {
  TFile *f = IO::OpenForReading(path);
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << path << std::endl;
    return kFALSE;
  }
  header.rate_kHz = ReadParameter(f, "rate_kHz");
  header.tauL_us = ReadParameter(f, "tauL_us");
  header.tauR_us = ReadParameter(f, "tauR_us");
  header.trapDecay_us = ReadParameter(f, "trapDecay_us");
  header.trapRise_ns = ReadParameter(f, "trapRise_ns");
  header.trapFlat_ns = ReadParameter(f, "trapFlat_ns");
  header.beamPeakADC = ReadParameter(f, "beamPeakADC");
  header.beamPeakSigmaFrac = ReadParameter(f, "beamPeakSigmaFrac");
  header.samplePeriod_ns = ReadParameter(f, "samplePeriod_ns");
  header.adcOffset = ReadParameter(f, "adcOffset");
  header.baselineMeanSamples = Int_t(ReadParameter(f, "baselineMeanSamples"));
  header.baselineHoldoff_ns = ReadParameter(f, "baselineHoldoff_ns");
  header.nEvents = Long64_t(ReadParameter(f, "nEvents"));

  TTree *tree = dynamic_cast<TTree *>(f->Get("hits"));
  if (!tree) {
    std::cerr << "No hits tree in " << path << std::endl;
    f->Close();
    delete f;
    return kFALSE;
  }
  Int_t channel = 0;
  ULong64_t timestamp = 0;
  Double_t energy = 0.0, dtPrev = 0.0, dtNext = 0.0, trueAmp = 0.0,
           prevTrueAmp = 0.0;
  Bool_t pileup = kFALSE;
  Long64_t eventIndex = 0;
  tree->SetBranchAddress("channel", &channel);
  tree->SetBranchAddress("timestamp", &timestamp);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchAddress("dtPrev", &dtPrev);
  tree->SetBranchAddress("dtNext", &dtNext);
  tree->SetBranchAddress("trueAmp", &trueAmp);
  tree->SetBranchAddress("prevTrueAmp", &prevTrueAmp);
  tree->SetBranchAddress("pileup", &pileup);
  tree->SetBranchAddress("eventIndex", &eventIndex);
  const Long64_t n = tree->GetEntries();
  hits.clear();
  hits.reserve(n);
  for (Long64_t i = 0; i < n; i++) {
    tree->GetEntry(i);
    SimHit h;
    h.channel = channel;
    h.timestamp = timestamp;
    h.energy = energy;
    h.dtPrev = dtPrev;
    h.dtNext = dtNext;
    h.trueAmp = trueAmp;
    h.prevTrueAmp = prevTrueAmp;
    h.pileup = pileup;
    h.eventIndex = eventIndex;
    hits.push_back(h);
  }
  f->Close();
  delete f;
  return kTRUE;
}

} // namespace SimOutput
