// plotbaseline - one figure comparing the baseline-restorer variants against
// the run-97 R data, for the CAEN EnergyFilterBaselineAvg / LF-filter question.
//
// The R-channel energy-vs-dt profile (normalised to the 38-60 us level,
// beam-peak units) is drawn for three baseline configurations:
//   A: 1024 samples @ 8 ns  (8.2 us window), no LF filter   -> excess, ruled out
//   B: 1024 samples @ /32   (262 us window), no LF filter  -> deficit, matches
//   C: 1024 samples @ 8 ns (8.2 us window) + LF filter     -> deficit, matches
// plus the measured run-97 R long-ends profile as points.
//
//   plotbaseline --a=<file> --b=<file> --c=<file> [--out=<dir>] [--format=png|pdf]
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLine.h>
#include <TPad.h>
#include <TProfile.h>
#include <TTree.h>
#include <iostream>
#include <vector>

namespace {

const Double_t kBeamPeak = 1600.0;
const Double_t kWindowLo = 0.5;
const Double_t kWindowHi = 1.5;

// R-channel (channel 1) deviation profile, normalised to the 38-60 us level.
TGraph *RProfile(const char *fname, Double_t xlo, Double_t xhi, Int_t nbins) {
  TFile *f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    std::cerr << "cannot open " << fname << std::endl;
    return new TGraph();
  }
  TTree *t = dynamic_cast<TTree *>(f->Get("hits"));
  Double_t energy = 0, dtPrev = 0;
  Int_t channel = -1;
  Bool_t pileup = kFALSE;
  t->SetBranchAddress("energy", &energy);
  t->SetBranchAddress("dtPrev", &dtPrev);
  t->SetBranchAddress("channel", &channel);
  t->SetBranchAddress("pileup", &pileup);
  TProfile *h = new TProfile("pb_h", "", nbins, xlo, xhi);
  TProfile *ref = new TProfile("pb_ref", "", 1, 38.0, 60.0);
  h->SetDirectory(nullptr);
  ref->SetDirectory(nullptr);
  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; i++) {
    t->GetEntry(i);
    if (channel != 1 || pileup || dtPrev < 0)
      continue;
    if (energy < kWindowLo * kBeamPeak || energy > kWindowHi * kBeamPeak)
      continue;
    h->Fill(dtPrev, energy);
    if (dtPrev >= 38.0 && dtPrev < 60.0)
      ref->Fill(dtPrev, energy);
  }
  const Double_t refv = ref->GetBinContent(1);
  TGraph *g = new TGraph();
  for (Int_t b = 1; b <= nbins; b++) {
    if (h->GetBinEntries(b) < 20)
      continue; // skip sparse bins (e.g. inside the trigger holdoff)
    g->SetPoint(g->GetN(), h->GetBinCenter(b),
                h->GetBinContent(b) / refv - 1.0);
  }
  f->Close();
  return g;
}

void PlotVariants(const char *fa, const char *fb, const char *fc) {
  TGraph *gA = RProfile(fa, 4.0, 40.0, 72);
  TGraph *gB = RProfile(fb, 4.0, 40.0, 72);
  TGraph *gC = RProfile(fc, 4.0, 40.0, 72);

  // Measured run-97 R long ends (native bins, from the paper).
  TGraph *gReal = new TGraph(6);
  const Double_t rx[] = {5.5, 7.0, 10.0, 15.0, 22.0, 32.0};
  const Double_t ry[] = {-0.197, -0.184, -0.142, -0.087, -0.037, -0.006};
  for (Int_t i = 0; i < 6; i++)
    gReal->SetPoint(i, rx[i], ry[i]);

  std::vector<TString> labels;
  labels.push_back("8.2 #mu s window, no LF  (ruled out)");
  labels.push_back("262 #mu s window  (1024 samples @ /32)");
  labels.push_back("8.2 #mu s window + LF filter  (#tau = 200 #mu s)");
  labels.push_back("run 97, R long ends  (data)");
  const Int_t legendW = Int_t(PlottingUtils::LegendWidthPx(labels));

  TPad *plot, *legpad;
  TCanvas *c =
      PlottingUtils::GetConfiguredCanvasWithSideLegend(plot, legpad, legendW);

  plot->SetLeftMargin(0.14);
  plot->cd();
  TH1F *frame = plot->DrawFrame(4.0, -0.26, 40.0, 0.10);
  frame->GetXaxis()->SetTitle("time since previous hit  #Delta t  [#mu s]");
  frame->GetYaxis()->SetTitle(
      "pulse height #minus 38-60 #mu s level  (beam-peak units)");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();

  TLine *zero = new TLine(4.0, 0.0, 40.0, 0.0);
  zero->SetLineStyle(2);
  zero->SetLineColor(kGray + 2);
  zero->SetLineWidth(PlottingUtils::GetLineWidth());
  zero->Draw();

  const Width_t lw = PlottingUtils::GetLineWidth();
  gA->SetLineColor(kRed);
  gA->SetLineStyle(2);
  gA->SetLineWidth(lw);
  gB->SetLineColor(kBlue + 1);
  gB->SetLineStyle(1);
  gB->SetLineWidth(lw);
  gC->SetLineColor(kGreen + 2);
  gC->SetLineStyle(7);
  gC->SetLineWidth(lw);
  gReal->SetMarkerStyle(20);
  gReal->SetMarkerColor(kBlack);
  gReal->SetMarkerSize(1.2);

  gB->Draw("L SAME");
  gA->Draw("L SAME");
  gC->Draw("L SAME");
  gReal->Draw("P SAME");

  legpad->cd();
  TLegend *leg = PlottingUtils::AddLegend(0.02, 0.98, 0.50, 0.95);
  leg->AddEntry(gA, labels[0], "L");
  leg->AddEntry(gB, labels[1], "L");
  leg->AddEntry(gC, labels[2], "L");
  leg->AddEntry(gReal, labels[3], "P");
  leg->Draw();

  PlottingUtils::ScaleFigure(plot);
  PlottingUtils::SaveFigure(c, "baseline_variants", "", PlotSaveOptions::kLINEAR);
  delete c;
}

} // namespace

int main(int argc, char **argv) {
  PlotSaveFormat format = PlotSaveFormat::kPNG;
  TString out = "doc/figs";
  TString fa = "/tmp/nix-shell.3kmr0P/opencode/cup/verify/caseA_8us.root";
  TString fb = "/tmp/nix-shell.3kmr0P/opencode/cup/verify/caseB_262us.root";
  TString fc = "/tmp/nix-shell.3kmr0P/opencode/cup/verify/caseC_8us_LF.root";
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a == "--format=png")
      format = PlotSaveFormat::kPNG;
    else if (a == "--format=pdf")
      format = PlotSaveFormat::kPDF;
    else if (a.BeginsWith("--out="))
      out = a(6, a.Length() - 6);
    else if (a.BeginsWith("--a="))
      fa = a(5, a.Length() - 5);
    else if (a.BeginsWith("--b="))
      fb = a(5, a.Length() - 5);
    else if (a.BeginsWith("--c="))
      fc = a(5, a.Length() - 5);
  }
  InitUtils::SetROOTPreferences(format, out, ".", kFALSE);
  PlotVariants(fa.Data(), fb.Data(), fc.Data());
  return 0;
}
