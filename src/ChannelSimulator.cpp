#include "ChannelSimulator.hpp"
#include <TMath.h>
#include <iostream>

ChannelSimulator::ChannelSimulator(const Config &cfg, Int_t channel,
                                   Double_t tau_ns, Double_t noiseSigmaADC,
                                   UInt_t seed, const EventGenerator &events)
    : cfg_(cfg), channel_(channel), events_(events),
      preamp_(tau_ns, cfg.DetRise_ns(), cfg.samplePeriod_ns),
      digitizer_(noiseSigmaADC, cfg.adcOffset, cfg.AdcMax(), seed),
      filter_(cfg.RiseSamples(), cfg.FlatSamples(), cfg.samplePeriod_ns,
              cfg.TrapDecay_ns()),
      baseline_(cfg.BaselineWindowSamples(), cfg.lfFilterOn,
                cfg.LfFilterTauNs(), cfg.samplePeriod_ns,
                cfg.baselineRestartThrADC, cfg.baselineRestartHoldSamples),
      resetSample_(0), nextEvent_(0),
      lastTriggerSample_(-1), lastTriggerTime_ns_(0.0), lastTriggerAmp_(0.0),
      snippetPair_(-1), snippetStart_(-1), snippetEnd_(-1),
      samplesGenerated_(0), samplesSkipped_(0), dropped_(0) {
  // Filter memory, then a full baseline window, then the peaking average and
  // the configured margin: after this many samples every state variable is
  // exact again. With the low-pass stage on, add ~5 time constants so the
  // baseline state has settled too.
  warmupSamples_ = filter_.Length() + cfg.BaselineWindowSamples() +
                   cfg.peakingAvgSamples;
  if (cfg.lfFilterOn && cfg.LfFilterTauNs() > 0.0)
    warmupSamples_ += Int_t(5.0 * cfg.LfFilterTauNs() / cfg.samplePeriod_ns);
  warmupSamples_ += cfg.warmupMarginSamples;
}

Long64_t ChannelSimulator::TriggerSample(Long64_t eventIndex) const {
  return Long64_t(TMath::Ceil(events_.Time_ns(eventIndex) /
                              cfg_.samplePeriod_ns));
}

void ChannelSimulator::ChooseSnippetPair() {
  // First pair with dt within 10 % of snippetDt_us, preceded by a quiet gap
  // so the snippet starts on a settled baseline.
  const Double_t want_ns = cfg_.snippetDt_us * 1000.0;
  const Double_t quiet_ns = 40000.0;
  for (Long64_t i = 1; i + 1 < events_.NEvents(); i++) {
    const Double_t dt = events_.Time_ns(i + 1) - events_.Time_ns(i);
    const Double_t gapBefore = events_.Time_ns(i) - events_.Time_ns(i - 1);
    const Double_t gapAfter =
        i + 2 < events_.NEvents() ? events_.Time_ns(i + 2) - events_.Time_ns(i + 1)
                                  : quiet_ns + 1.0;
    if (TMath::Abs(dt - want_ns) < 0.1 * want_ns && gapBefore > quiet_ns &&
        gapAfter > quiet_ns) {
      snippetPair_ = i;
      const Long64_t t0 = TriggerSample(i);
      snippetStart_ = t0 - Long64_t(10000.0 / cfg_.samplePeriod_ns);
      snippetEnd_ = t0 + Long64_t(60000.0 / cfg_.samplePeriod_ns);
      snippet_.t0_ns = events_.Time_ns(i);
      snippet_.pairDt_us = dt / 1000.0;
      return;
    }
  }
}

Long64_t ChannelSimulator::NextBarrier() const {
  Long64_t barrier = -1;
  if (nextEvent_ < events_.NEvents())
    barrier = TriggerSample(nextEvent_);
  if (snippetPair_ >= 0 && !snippet_.found && snippetStart_ >= 0) {
    if (barrier < 0 || snippetStart_ < barrier)
      barrier = snippetStart_;
  }
  return barrier;
}

void ChannelSimulator::HandleTriggers(Long64_t n) {
  while (nextEvent_ < events_.NEvents() && TriggerSample(nextEvent_) == n) {
    const Long64_t i = nextEvent_;
    nextEvent_++;
    preamp_.AddEvent(events_.Time_ns(i), events_.Amplitude(i));
    if (lastTriggerSample_ >= 0 &&
        n - lastTriggerSample_ < cfg_.TriggerHoldoffSamples()) {
      dropped_++;
      continue;
    }
    Readout r;
    r.eventIndex = i;
    const Int_t avg = TMath::Max(1, cfg_.peakingAvgSamples);
    r.startSample =
        n + cfg_.RiseSamples() + cfg_.PeakingDelaySamples() - avg / 2;
    r.endSample = r.startSample + avg - 1;
    r.sum = 0.0;
    r.frozenBaseline = baseline_.Freeze(n, cfg_.HoldoffSamples());
    if (lastTriggerSample_ >= 0) {
      r.dtPrev_us = (events_.Time_ns(i) - lastTriggerTime_ns_) / 1000.0;
      r.prevAmp = lastTriggerAmp_;
    } else {
      r.dtPrev_us = -1.0;
      r.prevAmp = 0.0;
    }
    r.pileup = kFALSE;
    pending_.push_back(r);
    lastTriggerSample_ = n;
    lastTriggerTime_ns_ = events_.Time_ns(i);
    lastTriggerAmp_ = events_.Amplitude(i);
  }
}

void ChannelSimulator::HandleReadouts(Long64_t n, Double_t trap) {
  size_t keep = 0;
  for (size_t j = 0; j < pending_.size(); j++) {
    Readout &r = pending_[j];
    if (n >= r.startSample && n <= r.endSample)
      r.sum += trap;
    if (n < r.endSample) {
      pending_[keep++] = r;
      continue;
    }
    const Double_t avg = Double_t(r.endSample - r.startSample + 1);
    SimHit h;
    h.channel = channel_;
    h.timestamp = ULong64_t(events_.Time_ns(r.eventIndex) * 1000.0 + 0.5);
    h.energy = r.sum / avg - r.frozenBaseline;
    h.dtPrev = r.dtPrev_us;
    h.dtNext = -1.0;
    h.trueAmp = events_.Amplitude(r.eventIndex);
    h.prevTrueAmp = r.prevAmp;
    h.pileup = r.pileup;
    h.eventIndex = r.eventIndex;
    hits_.push_back(h);
  }
  pending_.resize(keep);
}

void ChannelSimulator::RecordSnippet(Long64_t n, Int_t adc, Double_t trap) {
  if (snippetPair_ < 0 || snippet_.found)
    return;
  if (n < snippetStart_)
    return;
  const Double_t t_us =
      (Double_t(n) * cfg_.samplePeriod_ns - snippet_.t0_ns) / 1000.0;
  snippet_.t_us.push_back(t_us);
  snippet_.raw.push_back(Double_t(adc));
  snippet_.trap.push_back(trap);
  snippet_.baseline.push_back(baseline_.Value());
  if (n >= snippetEnd_)
    snippet_.found = kTRUE;
}

void ChannelSimulator::Run() {
  hits_.clear();
  hits_.reserve(events_.NEvents());
  pending_.clear();
  ChooseSnippetPair();
  preamp_.Reset();
  filter_.Reset();
  baseline_.Reset();
  resetSample_ = 0;
  Long64_t n = 0;
  while (kTRUE) {
    const Bool_t eventsLeft = nextEvent_ < events_.NEvents();
    const Bool_t snippetLeft = snippetPair_ >= 0 && !snippet_.found;
    if (!eventsLeft && pending_.empty() && !snippetLeft)
      break;
    // Skip idle time when nothing is being read out or recorded and the next
    // pulse (or snippet start) is more than two warm-ups away.
    if (pending_.empty() && !(snippetLeft && n >= snippetStart_)) {
      const Long64_t barrier = NextBarrier();
      if (barrier >= 0 && barrier - warmupSamples_ - n > warmupSamples_) {
        const Long64_t target = barrier - warmupSamples_;
        samplesSkipped_ += target - n;
        preamp_.SkipTo(target);
        filter_.Reset();
        baseline_.Reset();
        resetSample_ = target;
        n = target;
      }
    }
    HandleTriggers(n);
    const Double_t analogue = preamp_.Sample(n);
    const Int_t adc = digitizer_.Digitize(analogue);
    const Double_t trap = filter_.Process(Double_t(adc));
    if (n - resetSample_ >= filter_.Length())
      baseline_.Update(n, trap);
    HandleReadouts(n, trap);
    RecordSnippet(n, adc, trap);
    samplesGenerated_++;
    n++;
  }
  FillNextTimes();
}

// dtNext and the pileup flag need the following hit, so they are filled once
// the run is complete. A hit is pileup when the previous or the next trigger
// on the channel is closer than the trapezoid rise + flat top, i.e. when its
// own flat top overlaps another pulse's trapezoid, which is what the DPP-PHA
// pile-up flag marks on both pulses.
void ChannelSimulator::FillNextTimes() {
  const Double_t limit_us = (cfg_.trapRise_ns + cfg_.trapFlat_ns) / 1000.0;
  for (size_t i = 0; i < hits_.size(); i++) {
    SimHit &h = hits_[i];
    if (i + 1 < hits_.size())
      h.dtNext = (events_.Time_ns(hits_[i + 1].eventIndex) -
                  events_.Time_ns(h.eventIndex)) /
                 1000.0;
    const Bool_t prevClose = h.dtPrev >= 0.0 && h.dtPrev < limit_us;
    const Bool_t nextClose = h.dtNext >= 0.0 && h.dtNext < limit_us;
    h.pileup = prevClose || nextClose;
  }
}
