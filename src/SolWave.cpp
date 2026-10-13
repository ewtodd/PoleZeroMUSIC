#include "SolWave.hpp"
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <string>

namespace solwave {

const char *kGroupName[kNGroups] = {"long L",  "short L",  "long R", "short R",
                                    "strip 0", "strip 17", "grid",   "cathode"};
const char *kGroupTag[kNGroups] = {"longL",  "shortL",  "longR", "shortR",
                                   "strip0", "strip17", "grid",  "cathode"};
const Int_t kGroupColour[kNGroups] = {kBlue + 1,   kAzure + 7, kRed + 1,
                                      kOrange + 7, kGreen + 2, kGreen - 6,
                                      kMagenta + 1, kGray + 2};

Int_t GroupOf(Int_t ch) {
  if (ch == 63)
    return kStrip0;
  if (ch == 62)
    return kStrip17;
  if (ch == 58)
    return kGrid;
  if (ch == 59)
    return kCathode;
  if (ch < 0 || ch > 61)
    return -1;
  const Int_t rem = ch % 4;
  if (rem > 1)
    return -1;
  const Int_t strip = ch / 4 + 1;
  const Bool_t odd = (strip % 2) == 1;
  if (rem == 1)
    return odd ? kLongL : kShortL;
  return odd ? kShortR : kLongR;
}

Int_t StripOf(Int_t ch) {
  if (ch == 63)
    return 0;
  if (ch == 62)
    return 17;
  if (ch < 0 || ch > 61 || ch % 4 > 1)
    return -1;
  return ch / 4 + 1;
}

namespace {
TString Trim(const TString &s) {
  TString t = s;
  return t.Strip(TString::kBoth);
}
} // namespace

std::vector<ChSettings> ReadSettings(const TString &path) {
  std::vector<ChSettings> out(kMaxCh);
  std::ifstream in(path.Data());
  if (!in) {
    std::cerr << "cannot read settings " << path << std::endl;
    return out;
  }
  std::string line;
  while (std::getline(in, line)) {
    TString l = line.c_str();
    if (!l.BeginsWith("/ch/"))
      continue;
    const Int_t bang = l.Index('!');
    if (bang < 0)
      continue;
    const TString key = Trim(l(0, bang));
    TString rest = l(bang + 1, l.Length() - bang - 1);
    const Int_t last = rest.Last('!');
    if (last < 0)
      continue;
    const TString value = Trim(rest(last + 1, rest.Length() - last - 1));
    const TString body = key(4, key.Length() - 4); // "12/par/Name"
    const Int_t slash = body.Index('/');
    if (slash < 0 || body.Length() < slash + 6)
      continue;
    const Int_t ch = TString(body(0, slash)).Atoi();
    if (ch < 0 || ch >= kMaxCh)
      continue;
    const TString name = body(slash + 5, body.Length() - slash - 5);
    ChSettings &s = out[ch];
    if (name == "ChEnable")
      s.enabled = (value == "True");
    else if (name == "ChPreTriggerT")
      s.pre = Int_t(TMath::Nint(value.Atof() / kSample_ns));
    else if (name == "EnergyFilterRiseTimeT")
      s.rise = Int_t(TMath::Nint(value.Atof() / kSample_ns));
    else if (name == "EnergyFilterFlatTopT")
      s.flat = Int_t(TMath::Nint(value.Atof() / kSample_ns));
    else if (name == "EnergyFilterPeakingPosition")
      s.peakDelay = Int_t(value.Atof()); // percent for now; fixed below
    else if (name == "EnergyFilterPeakingAvg") {
      if (value == "OneShot")
        s.avg = 1;
      else if (value == "LowAVG")
        s.avg = 4;
      else if (value == "MediumAVG")
        s.avg = 16;
      else if (value == "HighAVG")
        s.avg = 64;
    } else if (name == "EnergyFilterPoleZeroT")
      s.pz_us = value.Atof() / 1000.0;
    else if (name == "WaveAnalogProbe1")
      s.probe1 = value;
    else if (name == "EnergyFilterBaselineAvg")
      s.baseline = value;
  }
  for (Int_t c = 0; c < kMaxCh; c++)
    out[c].peakDelay = Int_t(TMath::Nint(Double_t(out[c].peakDelay) / 100.0 *
                                         Double_t(out[c].flat)));
  return out;
}

Double_t Mean(const Int_t *v, Int_t from, Int_t n) {
  Double_t s = 0.0;
  for (Int_t i = 0; i < n; i++)
    s += Double_t(v[from + i]);
  return n > 0 ? s / Double_t(n) : 0.0;
}

Double_t MeanD(const std::vector<Double_t> &v, Int_t from, Int_t n) {
  Double_t s = 0.0;
  for (Int_t i = 0; i < n; i++)
    s += v[from + i];
  return n > 0 ? s / Double_t(n) : 0.0;
}

Bool_t LocatePulse(const Int_t *in, Int_t len, Double_t &pre, Double_t &pol,
                   Double_t &amp, Double_t &noise, Int_t &iStart, Int_t &i90,
                   Int_t &iExt) {
  const Double_t first = Mean(in, 0, 64);
  Double_t ext = 0.0;
  iExt = 0;
  for (Int_t i = 64; i < len; i++) {
    const Double_t d = Double_t(in[i]) - first;
    if (TMath::Abs(d) > TMath::Abs(ext)) {
      ext = d;
      iExt = i;
    }
  }
  pol = ext >= 0.0 ? 1.0 : -1.0;
  Int_t i10 = iExt;
  while (i10 > 0 && pol * (Double_t(in[i10]) - first) > 0.1 * pol * ext)
    i10--;
  const Int_t bEnd = i10 - 24;
  const Int_t bFrom = TMath::Max(0, bEnd - 200);
  if (bEnd - bFrom < 40)
    return kFALSE;
  pre = Mean(in, bFrom, bEnd - bFrom);
  Double_t ss = 0.0;
  for (Int_t i = bFrom; i < bEnd; i++)
    ss += (Double_t(in[i]) - pre) * (Double_t(in[i]) - pre);
  noise = TMath::Sqrt(ss / Double_t(bEnd - bFrom));
  amp = pol * (Double_t(in[iExt]) - pre);
  if (amp <= 0.0)
    return kFALSE;
  iStart = iExt;
  while (iStart > bEnd && pol * (Double_t(in[iStart]) - pre) > 0.1 * amp)
    iStart--;
  i90 = iStart;
  while (i90 < iExt && pol * (Double_t(in[i90]) - pre) < 0.9 * amp)
    i90++;
  return kTRUE;
}

} // namespace solwave
