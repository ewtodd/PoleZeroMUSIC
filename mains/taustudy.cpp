// Does the energy-vs-dt fit return the preamp decay time it was given?
//
// The four strip groups of MUSIC are two preamps seen through two pole-zero
// settings: the long ends run a 20 us pole-zero and the short ends a 10 us
// one. Simulating both settings with the same two preamps therefore produces
// all four groups, and the extraction can be checked against a known answer.
//
// Two fit models are compared. The one the task specifies,
//   p0 + p1 exp(-dt / p2),
// and the same with a linear term,
//   p0 + p1 exp(-dt / p2) + p3 dt,
// which accounts for the baseline average: a hit with a long gap behind it
// sits in a quiet stretch, where the average over the restorer's window is
// less depressed by pile-up, so the profile is lifted in proportion to the
// gap. Over the fitted range the gap is small against the window and that
// lift is close to linear.
//
// The uncertainty is taken from a bootstrap over the hits rather than from the
// fit's own covariance, and the spread over the start of the fit is added as a
// systematic.
//
//   taustudy <sim file> [<sim file> ...]
#include "InitUtils.hpp"
#include "SimOutput.hpp"
#include <TF1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

namespace {

const Double_t kWindowLo = 0.5;
const Double_t kWindowHi = 1.5;
const Double_t kProfMax_us = 150.0;
const Int_t kProfBins = 150;
const Double_t kFitHi_us = 120.0;
const Int_t kNBoot = 200;

struct Pair {
  Double_t dt;
  Double_t e;
};

struct FitOut {
  Bool_t ok = kFALSE;
  Double_t tau = 0.0, tauErr = 0.0, amp = 0.0, slope = 0.0;
  Double_t chi2 = 0.0;
  Int_t ndf = 0;
};

TProfile *MakeProfile(const char *name) {
  TProfile *p = new TProfile(name, "", kProfBins, 0.0, kProfMax_us);
  p->SetDirectory(nullptr);
  return p;
}

FitOut FitProfile(TProfile *p, Double_t lo, Bool_t withSlope) {
  FitOut r;
  const char *form = withSlope ? "[0] + [1] * exp(-x / [2]) + [3] * x"
                               : "[0] + [1] * exp(-x / [2])";
  TF1 f("taufit", form, lo, kFitHi_us);
  const Double_t tail = p->GetBinContent(p->FindBin(kFitHi_us - 10.0));
  const Double_t head = p->GetBinContent(p->FindBin(lo + 0.5));
  f.SetParameters(tail, head - tail, 12.0);
  if (withSlope)
    f.SetParameter(3, 0.0);
  f.SetParLimits(2, 1.0, 200.0);
  const Int_t status = Int_t(p->Fit(&f, "RQN0"));
  if (status != 0)
    return r;
  r.ok = kTRUE;
  r.amp = f.GetParameter(1);
  r.tau = f.GetParameter(2);
  r.tauErr = f.GetParError(2);
  r.slope = withSlope ? f.GetParameter(3) : 0.0;
  r.chi2 = f.GetChisquare();
  r.ndf = f.GetNDF();
  return r;
}

// Spread of the fitted tau over bootstrap resamples of the hits.
Double_t BootstrapError(const std::vector<Pair> &pairs, Double_t lo,
                        Bool_t withSlope, UInt_t seed) {
  TRandom3 rng(seed);
  const Int_t n = Int_t(pairs.size());
  std::vector<Double_t> taus;
  for (Int_t b = 0; b < kNBoot; b++) {
    TProfile *p = MakeProfile(Form("boot_%u_%d", seed, b));
    for (Int_t i = 0; i < n; i++) {
      const Int_t j = Int_t(rng.Uniform(0.0, Double_t(n)));
      p->Fill(pairs[j].dt, pairs[j].e);
    }
    FitOut f = FitProfile(p, lo, withSlope);
    if (f.ok && f.tau > 1.01 && f.tau < 199.0)
      taus.push_back(f.tau);
    delete p;
  }
  if (taus.size() < 10)
    return 0.0;
  Double_t sum = 0.0;
  for (size_t i = 0; i < taus.size(); i++)
    sum += taus[i];
  const Double_t mean = sum / Double_t(taus.size());
  Double_t v = 0.0;
  for (size_t i = 0; i < taus.size(); i++)
    v += (taus[i] - mean) * (taus[i] - mean);
  return TMath::Sqrt(v / Double_t(taus.size()));
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: taustudy <sim file> [...]" << std::endl;
    return 1;
  }
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG,
                                gSystem->WorkingDirectory(),
                                gSystem->WorkingDirectory(), kFALSE);
  const char *chName[2] = {"L", "R"};
  const Int_t nStart = 5;
  const Double_t start_us[nStart] = {8.0, 10.0, 12.0, 16.0, 20.0};

  std::cout << Form("%-26s %6s %8s %26s %26s", "group", "true", "PZ",
                    "p0 + p1 exp(-dt/p2)",
                    "p0 + p1 exp(-dt/p2) + p3 dt")
            << std::endl;
  std::cout << Form("%-26s %6s %8s %10s %7s %7s %10s %7s %7s", "", "us", "us",
                    "tau", "boot", "range", "tau", "boot", "range")
            << std::endl;

  for (Int_t f = 1; f < argc; f++) {
    SimOutput::RunHeader header;
    std::vector<SimHit> hits;
    if (!SimOutput::Read(argv[f], header, hits))
      continue;
    const Double_t freeze_us =
        (2.0 * header.trapRise_ns + header.trapFlat_ns) / 1000.0;
    for (Int_t ch = 0; ch < 2; ch++) {
      const Double_t tauTrue = ch == 0 ? header.tauL_us : header.tauR_us;
      TH1D spec("taustudy_spec", "", 128, 0.0, 4096.0);
      spec.SetDirectory(nullptr);
      for (size_t i = 0; i < hits.size(); i++)
        if (hits[i].channel == ch)
          spec.Fill(hits[i].energy);
      spec.GetXaxis()->SetRange(9, 128);
      const Double_t mode = spec.GetBinCenter(spec.GetMaximumBin());
      std::vector<Pair> pairs;
      for (size_t i = 0; i < hits.size(); i++) {
        const SimHit &h = hits[i];
        if (h.channel != ch || h.pileup || h.dtPrev < 0.0)
          continue;
        if (h.energy <= kWindowLo * mode || h.energy >= kWindowHi * mode)
          continue;
        Pair p;
        p.dt = h.dtPrev;
        p.e = h.energy;
        pairs.push_back(p);
      }
      if (pairs.size() < 1000)
        continue;

      Double_t tau[2] = {0.0, 0.0}, boot[2] = {0.0, 0.0}, range[2] = {0.0, 0.0};
      for (Int_t m = 0; m < 2; m++) {
        const Bool_t withSlope = m == 1;
        TProfile *p = MakeProfile(Form("prof_%d_%d_%d", f, ch, m));
        for (size_t i = 0; i < pairs.size(); i++)
          p->Fill(pairs[i].dt, pairs[i].e);
        FitOut base = FitProfile(p, freeze_us + 1.0, withSlope);
        tau[m] = base.ok ? base.tau : 0.0;
        // Spread over where the fit starts, as a systematic.
        Double_t lo = 1.0e9, hi = -1.0e9;
        for (Int_t s = 0; s < nStart; s++) {
          FitOut v = FitProfile(p, start_us[s], withSlope);
          if (!v.ok || v.tau <= 1.01 || v.tau >= 199.0)
            continue;
          lo = TMath::Min(lo, v.tau);
          hi = TMath::Max(hi, v.tau);
        }
        range[m] = hi > lo ? 0.5 * (hi - lo) : 0.0;
        delete p;
        boot[m] = BootstrapError(pairs, freeze_us + 1.0, withSlope,
                                 UInt_t(1000 * f + 10 * ch + m));
      }
      const char *group = ch == 0 ? "L" : "R";
      const char *end = header.trapDecay_us > 15.0 ? "long" : "short";
      std::cout << Form("%-26s %6.1f %8.1f %10.2f %7.2f %7.2f %10.2f %7.2f %7.2f",
                        Form("%s %s (%s)", end, group,
                             ch == 0 ? "25 us preamp" : "10 us preamp"),
                        tauTrue, header.trapDecay_us, tau[0], boot[0], range[0],
                        tau[1], boot[1], range[1])
                << std::endl;
    }
  }
  return 0;
}
