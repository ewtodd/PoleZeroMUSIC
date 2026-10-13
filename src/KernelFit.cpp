#include "KernelFit.hpp"
#include <TDecompSVD.h>
#include <TH1D.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <deque>

namespace {

struct Past {
  Double_t t_us;
  Double_t e;
  Int_t band;
};

Double_t TimeUs(const SimHit &h) { return Double_t(h.timestamp) * 1.0e-6; }

} // namespace

KernelFit::KernelFit(Basis basis, Int_t nTerms, Double_t maxDt_us,
                     Double_t freeze_us, Int_t nAmp)
    : basis_(basis), nTerms_(nTerms), maxDt_us_(maxDt_us),
      freeze_us_(freeze_us), nAmp_(TMath::Max(1, nAmp)), mode_(0.0),
      mean_(0.0) {
  if (basis_ == kExponentials) {
    // Geometric ladder from 4 us, which brackets both the preamp decay and
    // the slower pile-up depression of the baseline average.
    for (Int_t j = 0; j < nTerms_; j++)
      tau_us_.push_back(4.0 * TMath::Power(2.0, Double_t(j)));
  }
}

TString KernelFit::Name() const {
  TString n = basis_ == kLogBins ? TString::Format("%d log bins", nTerms_)
                                 : TString::Format("%d exponentials", nTerms_);
  if (nAmp_ > 1)
    n += TString::Format(" x %d amp", nAmp_);
  return n;
}

Int_t KernelFit::TermsPerBand() const {
  return nTerms_ + ((basis_ == kExponentials && freeze_us_ > 0.0) ? 1 : 0);
}

Int_t KernelFit::BandOf(Double_t e) const {
  // PulseHistory::AmpBinOf.
  if (nAmp_ <= 1 || !(mode_ > 0.0))
    return 0;
  const Int_t a = Int_t(e / mode_ + 0.5);
  return TMath::Max(0, TMath::Min(a, nAmp_ - 1));
}

void KernelFit::Evaluate(Double_t dt_us, std::vector<Double_t> &out,
                         Int_t offset) const {
  if (basis_ == kLogBins) {
    // Log-spaced from 1 us to maxDt, the same layout as PulseHistory.
    if (dt_us < 1.0 || dt_us >= maxDt_us_)
      return;
    const Double_t lo = TMath::Log10(1.0);
    const Double_t hi = TMath::Log10(maxDt_us_);
    const Int_t b =
        Int_t((TMath::Log10(dt_us) - lo) / (hi - lo) * Double_t(nTerms_));
    if (b >= 0 && b < nTerms_)
      out[offset + b] += 1.0;
    return;
  }
  if (dt_us >= maxDt_us_)
    return;
  for (Int_t j = 0; j < nTerms_; j++)
    out[offset + j] += TMath::Exp(-dt_us / tau_us_[j]);
  if (freeze_us_ > 0.0 && dt_us < freeze_us_)
    out[offset + nTerms_] += 1.0;
}

Bool_t KernelFit::Measure(const std::vector<SimHit> &hits, Int_t channel,
                          Double_t windowLo, Double_t windowHi,
                          const std::vector<Double_t> *weight) {
  const Int_t np = 1 + nAmp_ * TermsPerBand();
  coeff_.assign(np, 0.0);

  TH1D spec("kernelfit_spec", "", 128, 0.0, 4096.0);
  spec.SetDirectory(nullptr);
  for (size_t i = 0; i < hits.size(); i++)
    if (hits[i].channel == channel)
      spec.Fill(hits[i].energy);
  spec.GetXaxis()->SetRange(9, 128);
  mode_ = spec.GetBinCenter(spec.GetMaximumBin());
  if (!(mode_ > 0.0))
    return kFALSE;
  const Double_t lo = windowLo * mode_, hi = windowHi * mode_;

  Double_t sum = 0.0;
  Long64_t n = 0;
  for (size_t i = 0; i < hits.size(); i++) {
    const SimHit &h = hits[i];
    if (h.channel != channel || h.pileup)
      continue;
    if (h.energy > lo && h.energy < hi) {
      sum += h.energy;
      n++;
    }
  }
  if (n < 100)
    return kFALSE;
  mean_ = sum / Double_t(n);

  TMatrixD A(np, np);
  TVectorD b(np);
  A.Zero();
  b.Zero();
  std::deque<Past> past;
  std::vector<Double_t> x(np, 0.0);
  Long64_t nfit = 0;
  for (size_t i = 0; i < hits.size(); i++) {
    const SimHit &h = hits[i];
    if (h.channel != channel)
      continue;
    const Double_t t = TimeUs(h);
    while (!past.empty() && past.front().t_us < t - maxDt_us_)
      past.pop_front();
    const Bool_t use = !h.pileup && h.energy > lo && h.energy < hi;
    if (use) {
      std::fill(x.begin(), x.end(), 0.0);
      x[0] = 1.0;
      std::vector<Double_t> f(np, 0.0);
      for (std::deque<Past>::const_iterator it = past.begin();
           it != past.end(); ++it) {
        std::fill(f.begin(), f.end(), 0.0);
        Evaluate(t - it->t_us, f, 1 + it->band * TermsPerBand());
        for (Int_t j = 1; j < np; j++)
          x[j] += f[j] * it->e;
      }
      const Double_t y = h.energy - mean_;
      for (Int_t a = 0; a < np; a++) {
        if (x[a] == 0.0)
          continue;
        b[a] += x[a] * y;
        for (Int_t c = 0; c < np; c++)
          A(a, c) += x[a] * x[c];
      }
      nfit++;
    }
    Past p;
    p.t_us = t;
    p.e = weight ? (*weight)[i] : h.energy;
    p.band = BandOf(h.energy);
    past.push_back(p);
  }
  if (nfit < 100)
    return kFALSE;

  // Column scaling then SVD, as PulseHistory does: the far bins are nearly
  // collinear with the intercept and a plain inversion gives up on them.
  TVectorD scale(np);
  for (Int_t a = 0; a < np; a++)
    scale[a] = A(a, a) > 0.0 ? 1.0 / TMath::Sqrt(A(a, a)) : 1.0;
  TMatrixD As(np, np);
  TVectorD bs(np);
  for (Int_t a = 0; a < np; a++) {
    bs[a] = b[a] * scale[a];
    for (Int_t c = 0; c < np; c++)
      As(a, c) = A(a, c) * scale[a] * scale[c];
  }
  TDecompSVD svd(As);
  svd.SetTol(1.0e-10);
  Bool_t ok = kFALSE;
  TVectorD ps = svd.Solve(bs, ok);
  if (!ok)
    return kFALSE;
  for (Int_t a = 0; a < np; a++)
    coeff_[a] = ps[a] * scale[a];
  return kTRUE;
}

std::vector<Double_t> KernelFit::Apply(const std::vector<SimHit> &hits,
                                       Int_t channel,
                                       const std::vector<Double_t> *weight) const {
  const Int_t np = Int_t(coeff_.size());
  std::vector<Double_t> out(hits.size(), 0.0);
  for (size_t i = 0; i < hits.size(); i++)
    out[i] = hits[i].energy;
  if (np == 0)
    return out;
  std::deque<Past> past;
  std::vector<Double_t> f(np, 0.0);
  for (size_t i = 0; i < hits.size(); i++) {
    const SimHit &h = hits[i];
    if (h.channel != channel)
      continue;
    const Double_t t = TimeUs(h);
    while (!past.empty() && past.front().t_us < t - maxDt_us_)
      past.pop_front();
    Double_t shift = coeff_[0];
    for (std::deque<Past>::const_iterator it = past.begin(); it != past.end();
         ++it) {
      std::fill(f.begin(), f.end(), 0.0);
      Evaluate(t - it->t_us, f, 1 + it->band * TermsPerBand());
      for (Int_t j = 1; j < np; j++)
        shift += coeff_[j] * f[j] * it->e;
    }
    out[i] = h.energy - shift;
    Past p;
    p.t_us = t;
    p.e = weight ? (*weight)[i] : h.energy;
    p.band = BandOf(h.energy);
    past.push_back(p);
  }
  return out;
}
