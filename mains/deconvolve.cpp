// Deconvolve the preamp from the averaged long-strip pulse to get the
// detector current, then push that current through a slower preamp and
// compare with strips 0 and 17.
//
// A charge-sensitive preamp with feedback time constant tau turns a current
// i(t) into V(t) = int i(t') exp(-(t - t') / tau) dt', so the current is
//   i(t) = dV/dt + V / tau.
// tau is fixed by requiring that no current is left once the charge is
// collected: over a window well after the pulse, <dV/dt> + <V>/tau = 0. The
// current is then convolved with exp(-t / tau2) for the strip-0/17 preamps.
//
// Input: the pulse-start-aligned group averages written by wavecheck
// (input_aligned_<group>, t in us from the 10 % crossing, unit height).
//
//   deconvolve --in=<wavecheck.root> --out=<dir> [--tau2=100]
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

namespace {

const Double_t kSample_us = 0.008;
const Double_t kZeroFrom_us = 5.0, kZeroTo_us = 24.0; // no-current window
const Int_t kSmooth = 5; // samples, moving average before differentiating

struct Trace {
  std::vector<Double_t> t, v;
  Int_t n = 0;
};

Trace ReadTrace(TFile &f, const char *name) {
  Trace tr;
  TGraph *g = (TGraph *)f.Get(name);
  if (!g)
    return tr;
  tr.n = g->GetN();
  tr.t.resize(tr.n);
  tr.v.resize(tr.n);
  for (Int_t i = 0; i < tr.n; i++)
    g->GetPoint(i, tr.t[i], tr.v[i]);
  // Baseline from the samples before the pulse, then unit height.
  Double_t base = 0.0;
  Int_t nb = 0;
  for (Int_t i = 0; i < tr.n; i++)
    if (tr.t[i] > -1.9 && tr.t[i] < -0.4) {
      base += tr.v[i];
      nb++;
    }
  if (nb > 0)
    base /= Double_t(nb);
  Double_t mx = 0.0;
  for (Int_t i = 0; i < tr.n; i++) {
    tr.v[i] -= base;
    mx = TMath::Max(mx, tr.v[i]);
  }
  if (mx > 0.0)
    for (Int_t i = 0; i < tr.n; i++)
      tr.v[i] /= mx;
  return tr;
}

std::vector<Double_t> Smooth(const std::vector<Double_t> &v, Int_t w) {
  std::vector<Double_t> out(v.size(), 0.0);
  const Int_t n = Int_t(v.size());
  for (Int_t i = 0; i < n; i++) {
    Double_t s = 0.0;
    Int_t c = 0;
    for (Int_t k = -w; k <= w; k++) {
      const Int_t j = i + k;
      if (j >= 0 && j < n) {
        s += v[j];
        c++;
      }
    }
    out[i] = s / Double_t(c);
  }
  return out;
}

// tau from <dV/dt> + <V>/tau = 0 over the no-current window.
Double_t TauFromZeroCurrent(const Trace &tr) {
  Double_t sv = 0.0, sd = 0.0;
  for (Int_t i = 1; i + 1 < tr.n; i++) {
    if (tr.t[i] < kZeroFrom_us || tr.t[i] > kZeroTo_us)
      continue;
    sv += tr.v[i];
    sd += (tr.v[i + 1] - tr.v[i - 1]) / (2.0 * kSample_us);
  }
  return sd != 0.0 ? -sv / sd : 0.0;
}

std::vector<Double_t> Current(const Trace &tr, Double_t tau) {
  const std::vector<Double_t> v = Smooth(tr.v, kSmooth);
  std::vector<Double_t> i(tr.n, 0.0);
  for (Int_t k = 1; k + 1 < tr.n; k++)
    i[k] = (v[k + 1] - v[k - 1]) / (2.0 * kSample_us) + v[k] / tau;
  return i;
}

// V(t) of a preamp with time constant tau2 driven by i(t), unit peak.
std::vector<Double_t> Through(const std::vector<Double_t> &i, Double_t tau2) {
  std::vector<Double_t> v(i.size(), 0.0);
  const Double_t decay = TMath::Exp(-kSample_us / tau2);
  Double_t acc = 0.0, mx = 0.0;
  for (size_t k = 0; k < i.size(); k++) {
    acc = acc * decay + i[k] * kSample_us;
    v[k] = acc;
    mx = TMath::Max(mx, acc);
  }
  if (mx > 0.0)
    for (size_t k = 0; k < v.size(); k++)
      v[k] /= mx;
  return v;
}

TGraph *MakeGraph(const std::vector<Double_t> &t, const std::vector<Double_t> &v,
                  Double_t scale, const char *name) {
  TGraph *g = new TGraph(Int_t(t.size()));
  for (size_t k = 0; k < t.size(); k++)
    g->SetPoint(Int_t(k), t[k], v[k] * scale);
  g->SetName(name);
  return g;
}

Double_t PeakOf(const std::vector<Double_t> &v) {
  Double_t mx = 0.0;
  for (size_t k = 0; k < v.size(); k++)
    mx = TMath::Max(mx, v[k]);
  return mx > 0.0 ? mx : 1.0;
}

// Width at half maximum and time to collect 90 % of the charge, in us.
void CurrentTiming(const Trace &tr, const std::vector<Double_t> &i,
                   Double_t &fwhm, Double_t &t90, Double_t &residual) {
  const Double_t pk = PeakOf(i);
  Double_t tLo = 0.0, tHi = 0.0, total = 0.0;
  Bool_t up = kFALSE;
  for (Int_t k = 0; k < tr.n; k++) {
    if (tr.t[k] < -1.0 || tr.t[k] > 6.0)
      continue;
    total += i[k] * kSample_us;
    if (!up && i[k] >= 0.5 * pk) {
      tLo = tr.t[k];
      up = kTRUE;
    }
    if (up && i[k] >= 0.5 * pk)
      tHi = tr.t[k];
  }
  fwhm = tHi - tLo;
  Double_t acc = 0.0;
  t90 = 0.0;
  for (Int_t k = 0; k < tr.n; k++) {
    if (tr.t[k] < -1.0 || tr.t[k] > 6.0)
      continue;
    acc += i[k] * kSample_us;
    if (acc >= 0.9 * total) {
      t90 = tr.t[k];
      break;
    }
  }
  Double_t res = 0.0;
  Int_t nr = 0;
  for (Int_t k = 0; k < tr.n; k++)
    if (tr.t[k] >= kZeroFrom_us && tr.t[k] <= kZeroTo_us) {
      res += i[k];
      nr++;
    }
  residual = nr > 0 ? res / Double_t(nr) / pk : 0.0;
}

} // namespace

int main(int argc, char **argv) {
  TString in, out = "output/deconvolve";
  Double_t tau2 = 100.0;
  for (Int_t a = 1; a < argc; a++) {
    TString s = argv[a];
    if (s.BeginsWith("--in="))
      in = s(5, s.Length() - 5);
    else if (s.BeginsWith("--out="))
      out = s(6, s.Length() - 6);
    else if (s.BeginsWith("--tau2="))
      tau2 = TString(s(7, s.Length() - 7)).Atof();
  }
  if (in.Length() == 0) {
    std::cerr << "usage: deconvolve --in=wavecheck.root --out=dir [--tau2=us]"
              << std::endl;
    return 1;
  }
  gSystem->mkdir(out, kTRUE);
  const TString absOut = gSystem->IsAbsoluteFileName(out)
                             ? out
                             : TString(gSystem->WorkingDirectory()) + "/" + out;
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, absOut, absOut, kFALSE);

  TFile f(in);
  if (f.IsZombie()) {
    std::cerr << "cannot open " << in << std::endl;
    return 1;
  }
  const Int_t nRef = 2;
  const char *refTag[nRef] = {"longL", "longR"};
  const char *refName[nRef] = {"long L", "long R"};
  const Int_t refColour[nRef] = {kBlue + 1, kRed + 1};
  const Int_t nTest = 2;
  const char *testTag[nTest] = {"strip17", "strip0"};
  const char *testName[nTest] = {"strip 17", "strip 0"};
  const Int_t testColour[nTest] = {kGreen + 2, kMagenta + 1};

  Trace ref[nRef], test[nTest];
  std::vector<Double_t> cur[nRef], curTest[nTest];
  Double_t tau[nRef];
  for (Int_t r = 0; r < nRef; r++) {
    ref[r] = ReadTrace(f, Form("input_aligned_%s", refTag[r]));
    if (ref[r].n == 0) {
      std::cerr << "missing input_aligned_" << refTag[r] << std::endl;
      return 1;
    }
    tau[r] = TauFromZeroCurrent(ref[r]);
    cur[r] = Current(ref[r], tau[r]);
    Double_t fwhm = 0.0, t90 = 0.0, residual = 0.0;
    CurrentTiming(ref[r], cur[r], fwhm, t90, residual);
    std::cout << Form("%s: preamp tau from zero current after collection "
                      "%.2f us; current FWHM %.2f us, 90 %% of charge by "
                      "%.2f us after the 10 %% point; residual current in "
                      "%.0f-%.0f us %.4f of peak",
                      refName[r], tau[r], fwhm, t90, kZeroFrom_us, kZeroTo_us,
                      residual)
              << std::endl;
  }
  for (Int_t s = 0; s < nTest; s++) {
    test[s] = ReadTrace(f, Form("input_aligned_%s", testTag[s]));
    if (test[s].n == 0)
      continue;
    // The current this trace implies if its preamp were a clean tau2.
    curTest[s] = Current(test[s], tau2);
  }

  // Prediction: the long-L current through the slow preamp, compared with
  // each test strip sample by sample.
  std::vector<Double_t> pred = Through(cur[0], tau2);
  for (Int_t s = 0; s < nTest; s++) {
    if (test[s].n == 0)
      continue;
    Double_t maxDiff = 0.0, tMax = 0.0, late = 0.0;
    Int_t nLate = 0;
    const Int_t n = TMath::Min(test[s].n, Int_t(pred.size()));
    for (Int_t k = 0; k < n; k++) {
      if (test[s].t[k] < -0.5 || test[s].t[k] > 26.0)
        continue;
      const Double_t d = test[s].v[k] - pred[k];
      if (TMath::Abs(d) > maxDiff) {
        maxDiff = TMath::Abs(d);
        tMax = test[s].t[k];
      }
      if (test[s].t[k] > 10.0) {
        late += d;
        nLate++;
      }
    }
    std::cout << Form("%s against long-L current through %.0f us: largest "
                      "difference %.3f of height at %.2f us, mean difference "
                      "past 10 us %+.3f",
                      testName[s], tau2, maxDiff, tMax,
                      nLate > 0 ? late / Double_t(nLate) : 0.0)
              << std::endl;
  }

  // --- plots -----------------------------------------------------------------
  // Currents, unit peak.
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.42, 0.89, 0.58, 0.87);
    leg->SetTextSize(19.0);
    Bool_t first = kTRUE;
    for (Int_t r = 0; r < nRef; r++) {
      TGraph *g = MakeGraph(ref[r].t, cur[r], 1.0 / PeakOf(cur[r]),
                            Form("current_%s", refTag[r]));
      PlottingUtils::ConfigureGraph(
          g, refColour[r],
          ";t #minus pulse start [#mus];detector current, unit peak");
      g->SetLineWidth(3);
      if (first) {
        g->GetXaxis()->SetLimits(-1.0, 8.0);
        g->GetYaxis()->SetRangeUser(-0.4, 1.3);
        g->Draw("AL");
        first = kFALSE;
      } else {
        g->Draw("L SAME");
      }
      leg->AddEntry(g, Form("%s: dV/dt + V/#tau, #tau = %.1f #mus", refName[r],
                            tau[r]),
                    "l");
    }
    for (Int_t s = 0; s < nTest; s++) {
      if (test[s].n == 0)
        continue;
      TGraph *g = MakeGraph(test[s].t, curTest[s], 1.0 / PeakOf(curTest[s]),
                            Form("current_%s", testTag[s]));
      PlottingUtils::ConfigureGraph(g, testColour[s]);
      g->SetLineWidth(2);
      g->SetLineStyle(2);
      g->Draw("L SAME");
      leg->AddEntry(g,
                    Form("%s: same, assuming a clean %.0f #mus preamp",
                         testName[s], tau2),
                    "l");
    }
    TLine *zero = new TLine(-1.0, 0.0, 8.0, 0.0);
    zero->SetLineStyle(3);
    zero->Draw();
    leg->Draw();
    PlottingUtils::SaveFigure(c, "deconv_current", "",
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  // Each test strip against the prediction.
  for (Int_t s = 0; s < nTest; s++) {
    if (test[s].n == 0)
      continue;
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.36, 0.89, 0.60, 0.87);
    leg->SetTextSize(19.0);
    TGraph *gm = MakeGraph(test[s].t, test[s].v, 1.0,
                           Form("measured_%s", testTag[s]));
    PlottingUtils::ConfigureGraph(
        gm, testColour[s],
        Form("%s;t #minus pulse start [#mus];pulse, unit height",
             testName[s]));
    gm->SetLineWidth(3);
    gm->GetXaxis()->SetLimits(-2.0, 27.0);
    gm->GetYaxis()->SetRangeUser(-0.1, 1.5);
    gm->Draw("AL");
    TGraph *gp = MakeGraph(ref[0].t, pred, 1.0, Form("predicted_%s", testTag[s]));
    PlottingUtils::ConfigureGraph(gp, kBlack);
    gp->SetLineWidth(2);
    gp->SetLineStyle(2);
    gp->Draw("L SAME");
    TGraph *gl = MakeGraph(ref[0].t, ref[0].v, 1.0, Form("longL_%s", testTag[s]));
    PlottingUtils::ConfigureGraph(gl, kGray + 1);
    gl->SetLineWidth(2);
    gl->Draw("L SAME");
    leg->AddEntry(gm, Form("%s, measured", testName[s]), "l");
    leg->AddEntry(gp,
                  Form("long-L current through a %.0f #mus preamp", tau2),
                  "l");
    leg->AddEntry(gl, Form("long L as measured (%.1f #mus preamp)", tau[0]),
                  "l");
    leg->Draw();
    PlottingUtils::SaveFigure(c, Form("deconv_%s", testTag[s]), "",
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  std::cout << "wrote " << absOut << std::endl;
  return 0;
}
