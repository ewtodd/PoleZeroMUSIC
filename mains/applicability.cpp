// Figures for the applicability conditions of the pulse-history correction.
//
//  1. applicability_plane: the (rate, preamp decay time) plane with the
//     boundaries that the DAQ parameters set: the ADC range (the stacked
//     tails must not clip), trapezoid pile-up (hits within one trapezoid
//     length are rejected, not corrected), and the counting argument behind
//     the analytic baseline term (tau_pre << W). The MUSIC operating points
//     are marked.
//  2. applicability_grid: from scripts/applicability.sh, the corrected and
//     uncorrected pulse-height spread against rate for several pole-zero
//     mismatches, each relative to the matched case at the same rate.
//
//   applicability [--grid=output/applicability/grid.txt] [--out=doc/figs]
//                 [--format=pdf|png]
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace {

// Tail amplitude c of doc/polezero.tex Eq. (c), continuum form.
Double_t TailAmplitude(Double_t tauPre, Double_t tauPz, Double_t kUs,
                       Double_t lUs) {
  return (tauPz - tauPre) * tauPre / (kUs * tauPz) *
         (1.0 - TMath::Exp(kUs / tauPre)) * (1.0 - TMath::Exp(lUs / tauPre));
}

// Largest R tau_pre (pulses stacked on the input) for which the stacked
// tails stay inside the ADC: with the baseline at 10 % of the range and a
// pulse of height a (fraction of the range), mean + 3 sigma of the pedestal
// a (x + 3 sqrt(x / 2)) must stay below 0.9.
Double_t MaxStack(Double_t a) {
  Double_t lo = 0.0, hi = 1.0e4;
  for (Int_t i = 0; i < 60; i++) {
    const Double_t x = 0.5 * (lo + hi);
    if (a * (x + 3.0 * TMath::Sqrt(0.5 * x)) < 0.9)
      lo = x;
    else
      hi = x;
  }
  return lo;
}

struct GridRow {
  Double_t rate = 0.0, pz = 0.0, raw = 0.0, corr = 0.0;
};

std::vector<GridRow> ReadGrid(const TString &path) {
  std::vector<GridRow> rows;
  std::ifstream in(path.Data());
  std::string label;
  GridRow g;
  while (in >> label >> g.rate >> g.pz >> g.raw >> g.corr)
    rows.push_back(g);
  return rows;
}

} // namespace

int main(int argc, char **argv) {
  TString grid = "output/applicability/grid.txt", out = "doc/figs";
  PlotSaveFormat format = PlotSaveFormat::kPDF;
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a.BeginsWith("--grid="))
      grid = a(7, a.Length() - 7);
    else if (a.BeginsWith("--out="))
      out = a(6, a.Length() - 6);
    else if (a == "--format=png")
      format = PlotSaveFormat::kPNG;
    else if (a == "--format=pdf")
      format = PlotSaveFormat::kPDF;
  }
  InitUtils::SetROOTPreferences(format, out, ".", kFALSE);

  // ---- 1. the (rate, tau_pre) plane ---------------------------------------
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    c->SetLogx(1);
    c->SetLogy(1);
    const Double_t rLo = 1.0, rHi = 300.0, tLo = 2.0, tHi = 300.0;
    // Saturation curves tau = x_max(a) / R for three pulse-height fractions.
    const Int_t nA = 3;
    const Double_t aFrac[nA] = {0.04, 0.07, 0.18};
    const char *aName[nA] = {"pulse 4 % of range", "pulse 7 % of range",
                             "pulse 18 % of range"};
    const Int_t aColour[nA] = {kGreen + 2, kBlue + 1, kRed + 1};
    TLegend *leg = PlottingUtils::AddLegend(0.13, 0.56, 0.12, 0.40);
    Bool_t first = kTRUE;
    for (Int_t k = 0; k < nA; k++) {
      const Double_t xMax = MaxStack(aFrac[k]);
      TGraph *g = new TGraph();
      for (Int_t i = 0; i <= 200; i++) {
        const Double_t r = rLo * TMath::Power(rHi / rLo, Double_t(i) / 200.0);
        g->SetPoint(i, r, xMax / (r * 1.0e-3)); // R in kHz, tau in us
      }
      PlottingUtils::ConfigureGraph(
          g, aColour[k],
          ";hit rate per channel [kHz];preamplifier decay time #tau_{pre} "
          "[#mus]");
      if (first) {
        g->GetXaxis()->SetLimits(rLo, rHi);
        g->GetYaxis()->SetRangeUser(tLo, tHi);
        g->Draw("AL");
        first = kFALSE;
      } else {
        g->Draw("L SAME");
      }
      leg->AddEntry(g, Form("clips, %s", aName[k]), "l");
    }
    // Trapezoid pile-up: R (2k + l) T = 0.3.
    const Double_t trapLen[2] = {7.0, 4.4};
    const char *trapName[2] = {"7 #mus trapezoid", "4.4 #mus trapezoid"};
    for (Int_t k = 0; k < 2; k++) {
      const Double_t r = 0.3 / trapLen[k] * 1000.0;
      TLine *l = new TLine(r, tLo, r, tHi);
      l->SetLineStyle(k == 0 ? 2 : 3);
      l->SetLineWidth(PlottingUtils::GetLineWidth());
      l->SetLineColor(kGray + 2);
      l->Draw();
      leg->AddEntry(l, Form("30 %% pile-up, %s", trapName[k]), "l");
    }
    // Counting argument for the analytic baseline term: tau_pre < W / 10.
    TLine *w = new TLine(rLo, 26.2, rHi, 26.2);
    w->SetLineStyle(4);
    w->SetLineWidth(PlottingUtils::GetLineWidth());
    w->SetLineColor(kMagenta + 1);
    w->Draw();
    leg->AddEntry(w, "#tau_{pre} = W / 10", "l");
    // Operating points.
    const Int_t nP = 4;
    const Double_t pR[nP] = {27.0, 27.0, 40.0, 40.0};
    const Double_t pT[nP] = {25.0, 10.0, 10.0, 100.0};
    const char *pName[nP] = {"^{37}Cl L", "^{37}Cl R", "^{88}Sr strips",
                             "^{88}Sr strips 0/17"};
    const Int_t pStyle[nP] = {20, 21, 22, 29};
    const Double_t pDx[nP] = {1.12, 1.12, 1.12, 1.12};
    const Double_t pDy[nP] = {1.10, 0.72, 1.15, 1.10};
    for (Int_t k = 0; k < nP; k++) {
      TGraph *p = new TGraph(1);
      p->SetPoint(0, pR[k], pT[k]);
      p->SetMarkerStyle(pStyle[k]);
      p->SetMarkerSize(2.2);
      p->SetMarkerColor(kBlack);
      p->Draw("P SAME");
      TLatex *t = new TLatex(pR[k] * pDx[k], pT[k] * pDy[k], pName[k]);
      t->SetTextFont(43);
      t->SetTextSize(30);
      t->Draw();
    }
    leg->Draw();
    PlottingUtils::SaveFigure(c, "applicability_plane", "",
                              PlotSaveOptions::kLOG);
    delete c;
  }

  // ---- 2. the simulated grid ------------------------------------------------
  std::vector<GridRow> rows = ReadGrid(grid);
  if (rows.empty()) {
    std::cout << "no grid at " << grid << ", plane figure only" << std::endl;
    return 0;
  }
  std::map<Double_t, Double_t> floorAt; // matched raw RMS per rate
  for (size_t i = 0; i < rows.size(); i++)
    if (TMath::Abs(rows[i].pz - 10.0) < 1.0e-6)
      floorAt[rows[i].rate] = rows[i].raw;
  std::map<Double_t, TGraph *> gRaw, gCorr;
  for (size_t i = 0; i < rows.size(); i++) {
    const GridRow &g = rows[i];
    if (!floorAt.count(g.rate) || floorAt[g.rate] <= 0.0)
      continue;
    if (!gRaw.count(g.pz)) {
      gRaw[g.pz] = new TGraph();
      gCorr[g.pz] = new TGraph();
    }
    gRaw[g.pz]->SetPoint(gRaw[g.pz]->GetN(), g.rate, g.raw / floorAt[g.rate]);
    gCorr[g.pz]->SetPoint(gCorr[g.pz]->GetN(), g.rate,
                          g.corr / floorAt[g.rate]);
    std::cout << Form("rate %5.1f kHz  pz %5.1f us  raw %.1f  corr %.1f  "
                      "floor %.1f  -> %.3f / %.3f",
                      g.rate, g.pz, g.raw, g.corr, floorAt[g.rate],
                      g.raw / floorAt[g.rate], g.corr / floorAt[g.rate])
              << std::endl;
  }
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  c->SetLogx(1);
  TLegend *leg = PlottingUtils::AddLegend(0.14, 0.56, 0.48, 0.87);
  const Int_t colours[6] = {kAzure + 7, kGray + 2, kGreen + 2,
                            kOrange + 7, kRed + 1, kMagenta + 1};
  Int_t idx = 0;
  Bool_t first = kTRUE;
  for (std::map<Double_t, TGraph *>::iterator it = gRaw.begin();
       it != gRaw.end(); ++it, idx++) {
    const Double_t pz = it->first;
    const Double_t cTail = TailAmplitude(10.0, pz, 2.0, 5.0);
    TGraph *r = it->second, *k = gCorr[pz];
    PlottingUtils::ConfigureGraph(
        r, colours[idx % 6],
        ";hit rate per channel [kHz];spread / matched spread");
    r->SetLineStyle(2);
    PlottingUtils::ConfigureGraph(k, colours[idx % 6]);
    k->SetMarkerStyle(20);
    k->SetMarkerSize(1.4);
    k->SetMarkerColor(colours[idx % 6]);
    if (first) {
      r->GetXaxis()->SetLimits(4.0, 100.0);
      r->GetYaxis()->SetRangeUser(0.8, 3.0);
      r->Draw("AL");
      first = kFALSE;
    } else {
      r->Draw("L SAME");
    }
    k->Draw("LP SAME");
    leg->AddEntry(k, Form("#tau_{pz} = %.0f #mus, c = %+.2f", pz, cTail),
                  "lp");
  }
  TGraph *dashed = new TGraph(2);
  dashed->SetPoint(0, 4.0, 1.0);
  dashed->SetPoint(1, 100.0, 1.0);
  dashed->SetLineStyle(3);
  dashed->SetLineWidth(PlottingUtils::GetLineWidth());
  dashed->Draw("L SAME");
  leg->AddEntry((TObject *)nullptr, "dashed: uncorrected", "");
  leg->Draw();
  PlottingUtils::SaveFigure(c, "applicability_grid", "",
                            PlotSaveOptions::kLINEAR);
  delete c;
  return 0;
}
