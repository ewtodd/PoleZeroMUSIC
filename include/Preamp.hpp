#ifndef POLEZERO_PREAMP_HPP
#define POLEZERO_PREAMP_HPP

#include <Rtypes.h>
#include <vector>

// Charge-sensitive preamp with a resistive feedback of decay time tau. The
// detector current is a linear ramp of duration detRise, so the impulse
// response is the closed form of that ramp convolved with the exponential.
// Finished ramps are carried as one exponentially decaying state variable, so
// the sum over all past pulses costs O(1) per sample.
class Preamp {
public:
  Preamp(Double_t tau_ns, Double_t detRise_ns, Double_t samplePeriod_ns);

  // Register an event that starts at t_ns with amplitude amp (ADC). Events
  // must be added in time order and before the sample that contains them.
  void AddEvent(Double_t t_ns, Double_t amp);
  // Waveform value at sample n (t = n * T_s), without offset or noise. Must be
  // called for consecutive n once started.
  Double_t Sample(Long64_t n);
  // Advance the decayed-tail state to sample n without evaluating the samples
  // in between. No event may lie in between.
  void SkipTo(Long64_t n);
  void Reset();

private:
  struct RampEvent {
    Double_t t_ns;
    Double_t amp;
  };
  Double_t ImpulseResponse(Double_t t_ns) const;

  Double_t tau_ns_;
  Double_t detRise_ns_;
  Double_t samplePeriod_ns_;
  Double_t decayPerSample_;
  Double_t rampPlateau_;
  Double_t tail_;
  Long64_t current_;
  std::vector<RampEvent> ramping_;
};

#endif
