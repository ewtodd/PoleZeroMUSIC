#ifndef POLEZERO_CHANNELSIMULATOR_HPP
#define POLEZERO_CHANNELSIMULATOR_HPP

#include "BaselineRestorer.hpp"
#include "Config.hpp"
#include "Digitizer.hpp"
#include "EventGenerator.hpp"
#include "Preamp.hpp"
#include "SimHit.hpp"
#include "TrapezoidFilter.hpp"
#include <vector>

// One readout channel: preamp -> digitizer -> trapezoid -> baseline restorer
// -> trigger and energy extraction, streamed sample by sample on the
// digitizer clock. Idle stretches with no pulse, no pending readout and no
// snippet to record are skipped analytically: the preamp tail is decayed in
// closed form and the filter and baseline are rebuilt from scratch over a
// warm-up of more than their memory length before the next pulse (the
// trapezoid is a finite-impulse-response of k + l samples, and the baseline of
// baselineMeanSamples trapezoid samples), so the skipped output is identical
// in distribution to a continuous run.
class ChannelSimulator {
public:
  ChannelSimulator(const Config &cfg, Int_t channel, Double_t tau_ns,
                   Double_t noiseSigmaADC, UInt_t seed,
                   const EventGenerator &events);
  void Run();

  const std::vector<SimHit> &Hits() const { return hits_; }
  const WaveformSnippet &Snippet() const { return snippet_; }
  Long64_t SamplesGenerated() const { return samplesGenerated_; }
  Long64_t SamplesSkipped() const { return samplesSkipped_; }
  Long64_t Dropped() const { return dropped_; }

private:
  struct Readout {
    Long64_t eventIndex;
    Long64_t startSample;
    Long64_t endSample;
    Double_t sum;
    Double_t frozenBaseline;
    Double_t dtPrev_us;
    Double_t prevAmp;
    Bool_t pileup;
  };

  Long64_t TriggerSample(Long64_t eventIndex) const;
  void ChooseSnippetPair();
  void HandleTriggers(Long64_t n);
  void HandleReadouts(Long64_t n, Double_t trap);
  void RecordSnippet(Long64_t n, Int_t adc, Double_t trap);
  void FillNextTimes();
  Long64_t NextBarrier() const;

  const Config &cfg_;
  Int_t channel_;
  const EventGenerator &events_;
  Preamp preamp_;
  Digitizer digitizer_;
  TrapezoidFilter filter_;
  BaselineRestorer baseline_;

  Int_t warmupSamples_;
  // Trapezoid output for the ADC offset alone (Config::baselineIdeal).
  Double_t idealBaseline_;
  Long64_t resetSample_;
  Long64_t nextEvent_;
  Long64_t lastTriggerSample_;
  Double_t lastTriggerTime_ns_;
  Double_t lastTriggerAmp_;
  std::vector<Readout> pending_;
  std::vector<SimHit> hits_;

  Long64_t snippetPair_;
  Long64_t snippetStart_;
  Long64_t snippetEnd_;
  WaveformSnippet snippet_;

  Long64_t samplesGenerated_;
  Long64_t samplesSkipped_;
  Long64_t dropped_;
};

#endif
