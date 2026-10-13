#ifndef POLEZERO_CONFIG_HPP
#define POLEZERO_CONFIG_HPP

#include <Rtypes.h>
#include <TString.h>
#include <iostream>

// Every tunable of the simulation, with defaults taken from the CoMPASS /
// SOLARIS DPP-PHA settings of the 37Cl runs 97-137 (board 66222, file
// /labdata/MUSIC/New-37Cl/data_raw/music_exp1915_097XSetting_66222.dat) and
// from the measured 37Cl beam peak. Any field can be overridden on the command
// line as --name=value (see Config::ParseCommandLine).
struct Config {
  // Digitizer clock: CAEN V2740, 125 MS/s.
  Double_t samplePeriod_ns = 8.0;
  // DPP-PHA energy filter (EnergyFilterRiseTimeT / FlatTopT), long ends.
  Double_t trapRise_ns = 2000.0;
  Double_t trapFlat_ns = 3000.0;
  // Pole-zero setting (EnergyFilterPoleZeroT = 20000 ns on every long end).
  Double_t trapDecay_us = 20.0;
  // EnergyFilterPeakingPosition = 50 % of the flat top.
  Double_t peakingDelay_ns = 1500.0;
  // EnergyFilterPeakingAvg = MediumAVG = 16 samples (SOLARIS_DAQ: OneShot/
  // LowAVG/MediumAVG = 1/4/16); the energy is the mean of this many trapezoid
  // samples centred on the peaking point.
  Int_t peakingAvgSamples = 16;
  // EnergyFilterBaselineAvg = Medium = 1024 samples at 8 ns = 8.2 us (UM7788;
  // the 2740 PHA has no decimation). The long-window stand-in used before,
  // 32768 (262 us), is excluded by the measured baseline memory (README,
  // "Baseline restorer in production"); 16384 is High.
  Int_t baselineMeanSamples = 1024;
  // Clock decimation of the baseline accumulator: the running mean covers
  // baselineMeanSamples samples of a clock divided by this factor, i.e.
  // baselineMeanSamples * baselineDecimation ADC samples. 1 = the 125 MS/s ADC
  // clock (1024 -> 8.2 us); 32 = 1024 samples on a 3.9 MS/s clock (1024 ->
  // 262 us), a reading the measured baseline memory excludes.
  Int_t baselineDecimation = 1;
  // A first-order low-pass on the running-mean baseline (one reading of the
  // DPP-PHA EnergyFilterLFLimitation). The run-97 baseline memory allows at
  // most ~16-32 us; 16.384 us (a 2^11-sample shift) is the best fit together
  // with the restart below. Time constants of 65 us or more are excluded.
  Bool_t lfFilterOn = true;
  Double_t lfFilterTauUs = 16.384;
  // Conditional baseline restart: after the holdoff the baseline stays frozen
  // until the trapezoid is back within this many ADC of the frozen baseline,
  // for baselineRestartHoldSamples consecutive samples (> 1 keeps the
  // trapezoid's falling edge from triggering it). 80 ADC / 2 us fits the run-97
  // baseline memory and the ch-4 overshoot (README); 0 restarts
  // unconditionally after the holdoff, the documented behaviour.
  Double_t baselineRestartThrADC = 80.0;
  Int_t baselineRestartHoldSamples = 250;
  // The DPP-PHA baseline is inhibited for the whole trapezoid
  // (2 * rise + flat = 7000 ns) plus EnergyFilterBaselineGuardT (= 0).
  Double_t baselineHoldoff_ns = 7000.0;
  // Extra slack on top of the filter + baseline memory that the analytic idle
  // skip actually needs to rebuild a trigger exactly. 0 is the exact minimum;
  // raise it only as defensive margin against boundary rounding.
  Int_t warmupMarginSamples = 0;
  // TimeFilterRetriggerGuardT = 0; the real kernels are populated down to
  // dt = 1 us, so the effective re-arm time of the time filter is used.
  Double_t triggerHoldoff_ns = 1000.0;
  // ADC: 16 bit on the V2740 (ADC_Nbit = 16, Energy_Nbit = 16), DCOffset =
  // 50 % of full scale.
  Int_t adcBits = 16;
  Double_t adcOffset = 32768.0;
  // Preamp decay times, read off the hardware: the L chain is a 25 us variant,
  // the R chain a 10 us one. Both long-end groups run a 20 us pole-zero, so L
  // is slightly over-compensated (small positive tail) and R strongly
  // under-compensated (undershoot).
  Double_t tauL_us = 25.0;
  Double_t tauR_us = 10.0;
  // Charge collection (drift) time, a linear ramp.
  Double_t detRise_us = 1.0;
  // Beam peak of the long ends in 37Cl run 100 (L1: 1583, L3: 1644, L5: 1691
  // ADC), and its intrinsic spread.
  Double_t beamPeakADC = 1600.0;
  Double_t beamPeakSigmaFrac = 0.02;
  // White noise per raw sample, in ADC, set so that the trapezoid noise
  // matches the run-70 traces (19 L, 51 R in energy units; the trapezoid
  // passes 0.0905 of white noise): see README. The pulse-height width alone
  // would call for 960.
  Double_t noiseSigmaADC_L = 210.0;
  Double_t noiseSigmaADC_R = 563.0;
  // Beam rate and run length.
  Double_t rate_kHz = 40.0;
  Long64_t nEvents = 200000;
  UInt_t seed = 12345;
  // Output ROOT file of the simulation (or input of the analysis).
  TString outputFile = "output/sim_40kHz.root";
  // Snippet selection: a pair of pulses with dt near this value is saved as
  // raw / trapezoid / baseline TGraphs for figures.
  Double_t snippetDt_us = 5.0;

  // Derived, in samples.
  Int_t RiseSamples() const;
  Int_t FlatSamples() const;
  Int_t PeakingDelaySamples() const;
  Int_t HoldoffSamples() const;
  Int_t TriggerHoldoffSamples() const;
  Int_t BaselineWindowSamples() const {
    return baselineMeanSamples * (baselineDecimation > 0 ? baselineDecimation : 1);
  }
  Double_t LfFilterTauNs() const { return lfFilterTauUs * 1000.0; }
  Double_t AdcMax() const;
  Double_t TrapDecay_ns() const { return trapDecay_us * 1000.0; }
  Double_t TauL_ns() const { return tauL_us * 1000.0; }
  Double_t TauR_ns() const { return tauR_us * 1000.0; }
  Double_t DetRise_ns() const { return detRise_us * 1000.0; }
  Double_t MeanDt_ns() const { return 1.0e6 / rate_kHz; }

  // --name=value for every field above; returns kFALSE on an unknown name.
  Bool_t ParseCommandLine(Int_t argc, char **argv);
  void Print(std::ostream &os) const;
};

#endif
