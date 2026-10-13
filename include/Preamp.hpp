#ifndef POLEZERO_PREAMP_HPP
#define POLEZERO_PREAMP_HPP

#include <Rtypes.h>
#include <vector>

// Charge-sensitive preamp with a resistive feedback of decay time tau. The
// detector current is a linear ramp of duration detRise (the drift), so the
// impulse response is the closed form of that ramp convolved with the
// exponential:
//   h(t) = (tau / T_r) (1 - exp(-t / tau))                      0 <= t < T_r
//   h(t) = (tau / T_r) (1 - exp(-T_r / tau)) exp(-(t - T_r) / tau)   t >= T_r
// The waveform is sum_i A_i h(t - t_i). Sample by sample this is evaluated
// exactly with O(1) work per sample: every event whose ramp has finished
// contributes a pure exponential, and the sum of those is one state variable
// decayed by exp(-T_s / tau) each sample; only the (rare) events still inside
// their ramp are evaluated explicitly.
class Preamp {
public:
  Preamp(Double_t tau_ns, Double_t detRise_ns, Double_t samplePeriod_ns);

  // Closed-form impulse response at time t_ns after the event.
  Double_t ImpulseResponse(Double_t t_ns) const;

  // Register an event that starts at t_ns with amplitude amp (ADC). Events
  // must be added in time order and before the sample that contains them.
  void AddEvent(Double_t t_ns, Double_t amp);
  // Waveform value at sample n (t = n * T_s), without offset or noise. Must be
  // called for consecutive n once started; use Skip to jump over idle time.
  Double_t Sample(Long64_t n);
  // Advance the decayed-tail state from the current sample to sample n
  // without evaluating the samples in between (no events may lie in between).
  void SkipTo(Long64_t n);
  void Reset();

  Double_t Tau_ns() const { return tau_ns_; }

private:
  struct RampEvent {
    Double_t t_ns;
    Double_t amp;
  };
  Double_t tau_ns_;
  Double_t detRise_ns_;
  Double_t samplePeriod_ns_;
  Double_t decayPerSample_;
  Double_t rampPlateau_; // (tau / T_r) (1 - exp(-T_r / tau))
  Double_t tail_;        // summed exponential tails at the current sample
  Long64_t current_;     // sample index tail_ refers to
  std::vector<RampEvent> ramping_;
};

#endif
