#include "Config.hpp"
#include <TMath.h>
#include <cstdlib>
#include <cstring>

Int_t Config::RiseSamples() const {
  return Int_t(TMath::Nint(trapRise_ns / samplePeriod_ns));
}
Int_t Config::FlatSamples() const {
  return Int_t(TMath::Nint(trapFlat_ns / samplePeriod_ns));
}
Int_t Config::PeakingDelaySamples() const {
  return Int_t(TMath::Nint(peakingDelay_ns / samplePeriod_ns));
}
Int_t Config::HoldoffSamples() const {
  return Int_t(TMath::Nint(baselineHoldoff_ns / samplePeriod_ns));
}
Int_t Config::TriggerHoldoffSamples() const {
  return Int_t(TMath::Nint(triggerHoldoff_ns / samplePeriod_ns));
}
Double_t Config::AdcMax() const { return TMath::Power(2.0, adcBits) - 1.0; }

namespace {

Bool_t MatchDouble(const TString &name, const TString &value,
                   const char *key, Double_t &field) {
  if (name != key)
    return kFALSE;
  field = value.Atof();
  return kTRUE;
}

Bool_t MatchInt(const TString &name, const TString &value, const char *key,
                Int_t &field) {
  if (name != key)
    return kFALSE;
  field = value.Atoi();
  return kTRUE;
}

} // namespace

Bool_t Config::ParseCommandLine(Int_t argc, char **argv) {
  for (Int_t i = 1; i < argc; i++) {
    TString arg = argv[i];
    if (!arg.BeginsWith("--")) {
      std::cerr << "Unrecognised argument: " << arg << std::endl;
      return kFALSE;
    }
    arg.Remove(0, 2);
    const Int_t eq = arg.Index('=');
    if (eq < 0) {
      std::cerr << "Expected --name=value, got --" << arg << std::endl;
      return kFALSE;
    }
    const TString name = arg(0, eq);
    const TString value = arg(eq + 1, arg.Length() - eq - 1);
    Bool_t ok = MatchDouble(name, value, "samplePeriod_ns", samplePeriod_ns) ||
                MatchDouble(name, value, "trapRise_ns", trapRise_ns) ||
                MatchDouble(name, value, "trapFlat_ns", trapFlat_ns) ||
                MatchDouble(name, value, "trapDecay_us", trapDecay_us) ||
                MatchDouble(name, value, "peakingDelay_ns", peakingDelay_ns) ||
                MatchInt(name, value, "peakingAvgSamples", peakingAvgSamples) ||
                 MatchInt(name, value, "baselineMeanSamples",
                          baselineMeanSamples) ||
                 MatchInt(name, value, "baselineDecimation",
                          baselineDecimation) ||
                 MatchDouble(name, value, "lfFilterTauUs", lfFilterTauUs) ||
                 MatchDouble(name, value, "baselineRestartThrADC",
                             baselineRestartThrADC) ||
                 MatchInt(name, value, "baselineRestartHoldSamples",
                          baselineRestartHoldSamples) ||
                 MatchDouble(name, value, "baselineHoldoff_ns",
                             baselineHoldoff_ns) ||
                MatchInt(name, value, "warmupMarginSamples",
                         warmupMarginSamples) ||
                MatchDouble(name, value, "triggerHoldoff_ns",
                            triggerHoldoff_ns) ||
                MatchInt(name, value, "adcBits", adcBits) ||
                MatchDouble(name, value, "adcOffset", adcOffset) ||
                MatchDouble(name, value, "tauL_us", tauL_us) ||
                MatchDouble(name, value, "tauR_us", tauR_us) ||
                MatchDouble(name, value, "detRise_us", detRise_us) ||
                MatchDouble(name, value, "beamPeakADC", beamPeakADC) ||
                MatchDouble(name, value, "beamPeakSigmaFrac",
                            beamPeakSigmaFrac) ||
                MatchDouble(name, value, "noiseSigmaADC_L", noiseSigmaADC_L) ||
                MatchDouble(name, value, "noiseSigmaADC_R", noiseSigmaADC_R) ||
                MatchDouble(name, value, "rate_kHz", rate_kHz) ||
                MatchDouble(name, value, "snippetDt_us", snippetDt_us);
    if (!ok && name == "nEvents") {
      nEvents = Long64_t(std::atoll(value.Data()));
      ok = kTRUE;
    }
    if (!ok && name == "seed") {
      seed = UInt_t(std::strtoul(value.Data(), nullptr, 10));
      ok = kTRUE;
    }
    if (!ok && name == "outputFile") {
      outputFile = value;
      ok = kTRUE;
    }
    if (!ok && name == "lfFilterOn") {
      lfFilterOn = value == "1" || value == "on" || value == "true" ||
                   value == "True" || value == "kTRUE";
      ok = kTRUE;
    }
    if (!ok && name == "baselineRefreeze") {
      baselineRefreeze = value == "1" || value == "on" || value == "true" ||
                         value == "True" || value == "kTRUE";
      ok = kTRUE;
    }
    if (!ok && name == "baselineIdeal") {
      baselineIdeal = value == "1" || value == "on" || value == "true" ||
                      value == "True" || value == "kTRUE";
      ok = kTRUE;
    }
    if (!ok) {
      std::cerr << "Unknown parameter --" << name << std::endl;
      return kFALSE;
    }
  }
  return kTRUE;
}

void Config::Print(std::ostream &os) const {
  os << "Config:" << std::endl;
  os << "  samplePeriod_ns      = " << samplePeriod_ns << std::endl;
  os << "  trapRise_ns          = " << trapRise_ns << " (" << RiseSamples()
     << " samples)" << std::endl;
  os << "  trapFlat_ns          = " << trapFlat_ns << " (" << FlatSamples()
     << " samples)" << std::endl;
  os << "  trapDecay_us         = " << trapDecay_us << std::endl;
  os << "  peakingDelay_ns      = " << peakingDelay_ns << std::endl;
  os << "  peakingAvgSamples    = " << peakingAvgSamples << std::endl;
  os << "  baselineMeanSamples  = " << baselineMeanSamples << std::endl;
  os << "  baselineDecimation   = " << baselineDecimation << " (window "
     << BaselineWindowSamples() << " samples = "
     << BaselineWindowSamples() * samplePeriod_ns / 1000.0 << " us)" << std::endl;
  os << "  lfFilterOn           = " << (lfFilterOn ? "true" : "false")
     << " (tau " << lfFilterTauUs << " us)" << std::endl;
  os << "  baselineRestartThrADC = " << baselineRestartThrADC
     << " (hold " << baselineRestartHoldSamples << " samples)" << std::endl;
  os << "  baselineRefreeze     = " << (baselineRefreeze ? "true" : "false")
     << std::endl;
  os << "  baselineIdeal        = " << (baselineIdeal ? "true" : "false")
     << std::endl;
  os << "  baselineHoldoff_ns   = " << baselineHoldoff_ns << std::endl;
  os << "  warmupMarginSamples  = " << warmupMarginSamples << std::endl;
  os << "  triggerHoldoff_ns    = " << triggerHoldoff_ns << std::endl;
  os << "  adcBits              = " << adcBits << std::endl;
  os << "  adcOffset            = " << adcOffset << std::endl;
  os << "  tauL_us              = " << tauL_us << std::endl;
  os << "  tauR_us              = " << tauR_us << std::endl;
  os << "  detRise_us           = " << detRise_us << std::endl;
  os << "  beamPeakADC          = " << beamPeakADC << std::endl;
  os << "  beamPeakSigmaFrac    = " << beamPeakSigmaFrac << std::endl;
  os << "  noiseSigmaADC_L      = " << noiseSigmaADC_L << std::endl;
  os << "  noiseSigmaADC_R      = " << noiseSigmaADC_R << std::endl;
  os << "  rate_kHz             = " << rate_kHz << std::endl;
  os << "  nEvents              = " << nEvents << std::endl;
  os << "  seed                 = " << seed << std::endl;
  os << "  outputFile           = " << outputFile << std::endl;
  os << "  snippetDt_us         = " << snippetDt_us << std::endl;
}
