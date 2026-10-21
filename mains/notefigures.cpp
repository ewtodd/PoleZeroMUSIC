// Figures for doc/kernel_note.tex, in the PlottingUtils style.
//
//   note_kernel_run97   the pulse-history kernels fitted by PulseHistory on
//                       run 97 (long R, long L, short L), read from the 37Cl
//                       pipeline log;
//   note_sim_R/L        corrected and uncorrected pulse-height spread against
//                       rate, relative to the spread with no earlier pulses,
//                       for three readings of the baseline restorer, read from
//                       analyze logs (labels <reading>_r<rate in kHz>).
//
//   notefigures --log=<pipeline_fused.log> [--run=run97_c000]
//               --sim=<analysis.log>[,<analysis.log>...] [--out=doc/figs]
//               [--format=png|pdf]|pdf]   (png for the note: PDFs of side-legend
//                                       canvases get a wrong page size)|pdf]
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TBox.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPad.h>
#include <TString.h>
#include <TSystem.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

PlotSaveFormat gFormat = PlotSaveFormat::kPDF;

// PlottingUtils::SaveFigure prints a side-legend canvas, which is much wider
// than tall, to PDF on a rotated portrait page that LaTeX places badly; ROOT's
// "Portrait" option gives an unrotated page of the canvas's own shape.
void Save(TCanvas *c, const TString &name) {
  if (gFormat == PlotSaveFormat::kPNG) {
    PlottingUtils::SaveFigure(c, name, "", PlotSaveOptions::kLINEAR);
    return;
  }
  c->Modified();
  c->Update();
  const TString dir = PlottingUtils::GetPlotsBaseDir();
  gSystem->mkdir(dir, kTRUE);
  c->Print(dir + "/" + name + ".pdf", "pdf Portrait");
}

// PulseHistory's kernel bins: 12 per decade from 1 us, 30 of them.
const Int_t kNBins = 30;
Double_t BinCentre(Int_t i) {
  return TMath::Power(10.0, 0.0 + (i + 0.5) / 12.0);
}

// The three "kernel:" lines of a run's block in the pipeline log, in the
// order the log prints them: long L, long R, short L.
std::vector<std::vector<Double_t>> ReadKernels(const TString &log,
                                               const TString &run) {
  std::vector<std::vector<Double_t>> out;
  std::ifstream in(log.Data());
  std::string line;
  Bool_t inRun = kFALSE;
  while (std::getline(in, line)) {
    TString l(line.c_str());
    if (l.Contains("pulse history ")) {
      if (inRun)
        break;
      inRun = l.Contains(Form("pulse history %s:", run.Data()));
      continue;
    }
    if (!inRun || !l.Contains("kernel:"))
      continue;
    std::vector<Double_t> k;
    TObjArray *tok = l.Tokenize(" ");
    for (Int_t i = 0; i < tok->GetEntries(); i++) {
      TString t = ((TObjString *)tok->At(i))->GetString();
      const Ssiz_t c = t.Index("us:");
      if (c > 0)
        k.push_back(TString(t(c + 3, t.Length() - c - 3)).Atof());
    }
    delete tok;
    if (Int_t(k.size()) == kNBins)
      out.push_back(k);
  }
  return out;
}

struct Spread {
  Double_t raw[2] = {0, 0}, corr[2] = {0, 0}; // [L, R]
};

// "=== <reading>_r<rate>: ..." blocks and their "  L: ... RMS a -> corrected
// mean m RMS b ADC" lines.
void ReadSpreads(const TString &log,
                 std::map<TString, std::map<Double_t, Spread>> &out) {
  std::ifstream in(log.Data());
  std::string line;
  TString reading;
  Double_t rate = -1;
  while (std::getline(in, line)) {
    TString l(line.c_str());
    if (l.BeginsWith("=== ")) {
      TString label = l(4, l.Index(":") - 4);
      const Ssiz_t r = label.Last('r');
      reading = label(0, r - 1);
      rate = TString(label(r + 1, label.Length() - r - 1)).Atof();
      continue;
    }
    const Int_t side = l.BeginsWith("  L: ") ? 0 : l.BeginsWith("  R: ") ? 1 : -1;
    if (side < 0 || rate < 0 || !l.Contains("corrected mean"))
      continue;
    const Ssiz_t a = l.Index(" RMS ");
    const Ssiz_t b = l.Index(" RMS ", a + 5);
    out[reading][rate].raw[side] = TString(l(a + 5, 10)).Atof();
    out[reading][rate].corr[side] = TString(l(b + 5, 10)).Atof();
  }
}

void KernelFigure(const std::vector<std::vector<Double_t>> &k) {
  const char *name[3] = {"long L", "long R", "short L"};
  const Int_t order[3] = {1, 0, 2}; // R first: the largest effect
  std::vector<TString> labels;
  for (Int_t j = 0; j < 3; j++)
    labels.push_back(name[order[j]]);
  labels.push_back("read point inside the earlier trapezoid");
  labels.push_back("baseline still frozen");
  TPad *plot, *legpad;
  TCanvas *c = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legpad, Int_t(PlottingUtils::LegendWidthPx(labels)));
  plot->cd();
  plot->SetLogx();
  plot->SetLeftMargin(0.16);
  TH1F *frame = plot->DrawFrame(1.0, -0.2, 350.0, 0.75);
  frame->GetXaxis()->SetTitle("time since the earlier pulse  #Deltat  [#mus]");
  frame->GetYaxis()->SetTitle("shift per unit earlier pulse height");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();
  // Gap ranges fixed by the filter settings: the earlier trapezoid still
  // covers the read point below t_m... (k+l) - t_m = 3.5 us; the baseline
  // latched before the earlier pulse is still used below the 7 us freeze.
  TBox *inside = new TBox(1.0, -0.2, 3.5, 0.75);
  inside->SetFillColorAlpha(kGray, 0.35);
  inside->Draw();
  TBox *frozen = new TBox(3.5, -0.2, 7.0, 0.75);
  frozen->SetFillColorAlpha(kAzure - 9, 0.35);
  frozen->Draw();
  TLine *zero = new TLine(1.0, 0.0, 350.0, 0.0);
  zero->SetLineStyle(2);
  zero->SetLineColor(kGray + 2);
  zero->SetLineWidth(PlottingUtils::GetLineWidth());
  zero->Draw();
  const std::vector<Int_t> colours = PlottingUtils::GetDefaultColors();
  std::vector<TGraph *> g;
  for (Int_t j = 0; j < 3; j++) {
    const std::vector<Double_t> &v = k[order[j]];
    TGraph *gr = new TGraph();
    for (Int_t i = 0; i < kNBins; i++)
      gr->SetPoint(i, BinCentre(i), v[i]);
    PlottingUtils::ConfigureGraph(gr, colours[j]);
    gr->SetMarkerStyle(20 + j);
    gr->SetMarkerColor(colours[j]);
    gr->SetMarkerSize(1.2);
    gr->Draw("LP SAME");
    g.push_back(gr);
  }
  plot->RedrawAxis();
  legpad->cd();
  TLegend *leg = PlottingUtils::AddLegend(0.02, 0.98, 0.55, 0.95);
  for (Int_t j = 0; j < 3; j++)
    leg->AddEntry(g[j], labels[j], "LP");
  leg->AddEntry(inside, labels[3], "F");
  leg->AddEntry(frozen, labels[4], "F");
  leg->Draw();
  PlottingUtils::ScaleFigure(plot);
  Save(c, "note_kernel_run97");
  delete c;
}

void SpreadFigure(const std::map<TString, std::map<Double_t, Spread>> &s,
                  Int_t side) {
  struct R {
    const char *key, *label;
  };
  const std::vector<R> readings = {
      {"ideal", "baseline never takes in the tail"},
      {"wait", "frozen until within 80 ADC for 2 #mus"},
      {"follow", "averaging restarts after the 7 #mus freeze"}};
  std::vector<TString> labels;
  for (const R &r : readings)
    labels.push_back(r.label);
  labels.push_back("before correction (dashed)");
  labels.push_back("after correction (solid)");
  TPad *plot, *legpad;
  TCanvas *c = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legpad, Int_t(PlottingUtils::LegendWidthPx(labels)));
  plot->cd();
  plot->SetLeftMargin(0.16);
  const Double_t yHi = side == 1 ? 2.8 : 2.3;
  TH1F *frame = plot->DrawFrame(0.0, 0.9, 85.0, yHi);
  frame->GetXaxis()->SetTitle("rate per channel  [kHz]");
  frame->GetYaxis()->SetTitle("spread / spread of isolated pulses");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  TLatex *chain = new TLatex(4.0, yHi - 0.08 * (yHi - 0.9),
                             side ? "long R ends" : "long L ends");
  chain->SetTextFont(42);
  chain->SetTextSize(0.06);
  chain->SetTextAlign(13);
  chain->Draw();
  TLine *one = new TLine(0.0, 1.0, 85.0, 1.0);
  one->SetLineStyle(3);
  one->SetLineColor(kGray + 2);
  one->SetLineWidth(PlottingUtils::GetLineWidth());
  one->Draw();
  const std::vector<Int_t> colours = PlottingUtils::GetDefaultColors();
  std::vector<TGraph *> corr;
  for (size_t j = 0; j < readings.size(); j++) {
    auto it = s.find(readings[j].key);
    if (it == s.end() || !it->second.count(0.5))
      continue;
    const Double_t floor = it->second.at(0.5).corr[side];
    TGraph *gr = new TGraph(), *gc = new TGraph();
    for (const auto &p : it->second) {
      if (p.first < 1.0)
        continue;
      gr->SetPoint(gr->GetN(), p.first, p.second.raw[side] / floor);
      gc->SetPoint(gc->GetN(), p.first, p.second.corr[side] / floor);
      std::cout << Form("  %s %-7s %5.1f kHz: %6.1f -> %6.1f ADC (floor %.1f)",
                        side ? "R" : "L", readings[j].key, p.first,
                        p.second.raw[side], p.second.corr[side], floor)
                << std::endl;
    }
    PlottingUtils::ConfigureGraph(gr, colours[j]);
    gr->SetLineStyle(2);
    gr->SetMarkerStyle(24 + j);
    gr->SetMarkerColor(colours[j]);
    gr->SetMarkerSize(1.2);
    gr->Draw("LP SAME");
    PlottingUtils::ConfigureGraph(gc, colours[j]);
    gc->SetMarkerStyle(20 + j);
    gc->SetMarkerColor(colours[j]);
    gc->SetMarkerSize(1.2);
    gc->Draw("LP SAME");
    corr.push_back(gc);
  }
  legpad->cd();
  TLegend *leg = PlottingUtils::AddLegend(0.02, 0.98, 0.55, 0.95);
  for (size_t j = 0; j < corr.size(); j++)
    leg->AddEntry(corr[j], labels[j], "LP");
  leg->AddEntry((TObject *)nullptr, labels[3], "");
  leg->AddEntry((TObject *)nullptr, labels[4], "");
  leg->Draw();
  PlottingUtils::ScaleFigure(plot);
  Save(c, side ? "note_sim_R" : "note_sim_L");
  delete c;
}

} // namespace

int main(int argc, char **argv) {
  TString out = "doc/figs", log, run = "run97_c000", sims;
  PlotSaveFormat &format = gFormat;
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a.BeginsWith("--out="))
      out = a(6, a.Length() - 6);
    else if (a.BeginsWith("--log="))
      log = a(6, a.Length() - 6);
    else if (a.BeginsWith("--run="))
      run = a(6, a.Length() - 6);
    else if (a.BeginsWith("--sim="))
      sims = a(6, a.Length() - 6);
    else if (a == "--format=png")
      format = PlotSaveFormat::kPNG;
  }
  InitUtils::SetROOTPreferences(format, out, ".", kFALSE);
  if (!log.IsNull()) {
    std::vector<std::vector<Double_t>> k = ReadKernels(log, run);
    if (k.size() < 3)
      std::cerr << "found " << k.size() << " kernels for " << run << " in "
                << log << std::endl;
    else
      KernelFigure(k);
  }
  if (!sims.IsNull()) {
    std::map<TString, std::map<Double_t, Spread>> s;
    TObjArray *files = sims.Tokenize(",");
    for (Int_t i = 0; i < files->GetEntries(); i++)
      ReadSpreads(((TObjString *)files->At(i))->GetString(), s);
    delete files;
    SpreadFigure(s, 1);
    SpreadFigure(s, 0);
  }
  return 0;
}
