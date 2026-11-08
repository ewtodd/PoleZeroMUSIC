#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "Preamp.hpp"
#include "PulseHistory.hpp"
#include "PulseHistoryAdapter.hpp"
#include "SimHit.hpp"
#include "TrapezoidFilter.hpp"
#include <TArrow.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPad.h>
#include <TParameter.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

// The 37Cl runs 97-137 DPP-PHA settings (board 66222).
const Double_t kSamplePeriodNs = 8.0;
const Double_t kTrapRiseNs = 2000.0;
const Double_t kTrapFlatNs = 3000.0;
const Double_t kPeakingDelayNs = 1500.0;
const Double_t kBaselineHoldoffNs = 7000.0;
const Double_t kPoleZeroUs = 20.0;

// PulseHistory's kernel bins: 12 per decade from 1 us, 30 of them.
const Int_t kNBins = 30;

Int_t Samples(Double_t ns) { return Int_t(TMath::Nint(ns / kSamplePeriodNs)); }

Double_t BinCentre(Int_t i) { return TMath::Power(10.0, (i + 0.5) / 12.0); }

// ----------------------------------------------------------- input parsing
std::vector<std::vector<Double_t>> ReadKernels(const TString &log,
                                               const TString &run) {
  std::vector<std::vector<Double_t>> kernels;
  std::ifstream in(log.Data());
  std::string line;
  Bool_t inRun = kFALSE;
  while (std::getline(in, line)) {
    TString text(line.c_str());
    if (text.Contains("pulse history ")) {
      if (inRun)
        break;
      inRun = text.Contains(Form("pulse history %s:", run.Data()));
      continue;
    }
    if (!inRun || !text.Contains("kernel:"))
      continue;
    std::vector<Double_t> values;
    TObjArray *tokens = text.Tokenize(" ");
    for (Int_t i = 0; i < tokens->GetEntries(); i++) {
      TString token = ((TObjString *)tokens->At(i))->GetString();
      const Ssiz_t col = token.Index("us:");
      if (col > 0)
        values.push_back(
            TString(token(col + 3, token.Length() - col - 3)).Atof());
    }
    delete tokens;
    if (Int_t(values.size()) == kNBins)
      kernels.push_back(values);
  }
  return kernels;
}

struct Spread {
  Double_t raw[2];
  Double_t corr[2];
  Spread() {
    raw[0] = 0.0;
    raw[1] = 0.0;
    corr[0] = 0.0;
    corr[1] = 0.0;
  }
};

typedef std::map<Double_t, Spread> RateSpreads;
typedef std::map<TString, RateSpreads> ReadingSpreads;

// The spread table's runs (scripts/spread_inputs.sh): three behaviours of the
// baseline restorer, each at the isolated rate and the three table rates.
struct SpreadReading {
  const char *key;
  const char *label;
};
const Int_t kNSpreadReadings = 3;
const SpreadReading kSpreadReadings[kNSpreadReadings] = {
    {"ideal", "never takes in the tail"},
    {"wait", "frozen until within \\SI{80}{ADC}"},
    {"follow", "restarts after the freeze"}};
const Double_t kIsolatedRate = 0.5;
const Int_t kNSpreadRates = 3;
const Double_t kSpreadRates[kNSpreadRates] = {10.0, 44.0, 80.0};
// The beam window of the former analyze program, in units of the beam peak:
// wide enough for the pole-zero distortions, narrow enough to exclude the
// merged double pulses at twice the peak.
const Double_t kSpreadWindowLo = 0.5;
const Double_t kSpreadWindowHi = 1.5;

// ----------------------------------------------------------- shared drawing
const Int_t kRightColour = kRed + 1;
const Int_t kLeftColour = kBlue + 1;
const Int_t kShortLeftColour = kCyan + 2;
const Int_t kMatchedColour = kGreen + 2;

TH1F *Frame(TVirtualPad *pad, Double_t x1, Double_t y1, Double_t x2,
            Double_t y2, const char *xTitle, const char *yTitle) {
  TH1F *frame = pad->DrawFrame(x1, y1, x2, y2);
  frame->GetXaxis()->SetTitle(xTitle);
  frame->GetYaxis()->SetTitle(yTitle);
  frame->GetXaxis()->SetNdivisions(506);
  return frame;
}

// A legend inside the frame, sized to its labels and anchored on the right
// (or the left) edge.
TLegend *Legend(const std::vector<TString> &labels, Double_t y1, Double_t y2,
                const TString &header = "", Bool_t left = kFALSE) {
  const Double_t width =
      PlottingUtils::LegendWidthPx(labels, header) / gPad->GetWw();
  const Double_t x2 = left ? 0.17 + width : 0.88;
  TLegend *legend = PlottingUtils::AddLegend(x2 - width, x2, y1, y2);
  legend->SetFillStyle(1001);
  if (!header.IsNull())
    legend->SetHeader(header);
  return legend;
}

void Reference(Double_t x1, Double_t x2, Double_t y) {
  TLine *line = new TLine(x1, y, x2, y);
  line->SetLineStyle(3);
  line->SetLineColor(kGray + 2);
  line->SetLineWidth(PlottingUtils::GetLineWidth());
  line->Draw();
}

void TimeMark(Double_t x, Double_t yLo, Double_t yHi) {
  TLine *line = new TLine(x, yLo, x, yHi);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->SetLineWidth(PlottingUtils::GetLineWidth());
  line->Draw();
}

void TimeArrow(Double_t x1, Double_t x2, Double_t y, Double_t textY,
               const char *text) {
  TArrow *arrow = new TArrow(x1, y, x2, y, 0.012, "<|>");
  arrow->SetLineWidth(PlottingUtils::GetLineWidth());
  arrow->SetFillColor(kBlack);
  arrow->Draw();
  TLatex *label = new TLatex(0.5 * (x1 + x2), textY, text);
  label->SetTextFont(43);
  label->SetTextSize(30);
  label->SetTextAlign(21);
  label->Draw();
}

void Note(Double_t x, Double_t y, const char *text, Int_t align = 12,
          Int_t colour = kBlack) {
  TLatex *label = new TLatex(x, y, text);
  label->SetTextFont(43);
  label->SetTextSize(30);
  label->SetTextAlign(align);
  label->SetTextColor(colour);
  label->Draw();
}

void PanelLabel(const char *text) {
  TLatex *label = PlottingUtils::AddText(text, 0.02, 0.97);
  label->SetTextFont(63);
  label->SetTextSize(36);
  label->SetTextAlign(13);
}

// --------------------------------------------- simulated filter outputs
// c of the note's Eq. (eq:c): after (k+l)T a unit pulse leaves
// -c exp(-t/tau_pre).
Double_t TailAmplitude(Double_t tauPre, Double_t tauPz, Double_t riseUs,
                       Double_t flatTopEndUs) {
  return (tauPz - tauPre) * tauPre / (riseUs * tauPz) *
         (1.0 - TMath::Exp(riseUs / tauPre)) *
         (1.0 - TMath::Exp(flatTopEndUs / tauPre));
}

TString TailLabel(Double_t c) {
  if (TMath::Abs(c) < 0.005)
    return "#it{c} = 0";
  if (c < 0.0)
    return Form("#it{c} = #minus#kern[0.15]{%.2f}", -c);
  return Form("#it{c} = +%.2f", c);
}

// Filter output for unit pulses at arrivalsUs, against time from the first
// arrival. A one-sample detector ramp keeps the trapezoid corners exact.
TGraph *StepResponse(Double_t tauPreUs, Double_t tauPzUs,
                     const std::vector<Double_t> &arrivalsUs, Double_t tMinUs,
                     Double_t tMaxUs, Int_t stride = 1) {
  const Long64_t pad = Long64_t(2000.0 / kSamplePeriodNs);
  const Long64_t nSamples =
      pad + Long64_t(tMaxUs * 1000.0 / kSamplePeriodNs) + 1;
  Preamp preamp(tauPreUs * 1000.0, kSamplePeriodNs, kSamplePeriodNs);
  TrapezoidFilter filter(Samples(kTrapRiseNs), Samples(kTrapFlatNs),
                         kSamplePeriodNs, tauPzUs * 1000.0);
  for (Int_t a = 0; a < Int_t(arrivalsUs.size()); a++)
    preamp.AddEvent(Double_t(pad) * kSamplePeriodNs + arrivalsUs[a] * 1000.0,
                    1.0);
  TGraph *graph = new TGraph();
  Int_t point = 0;
  for (Long64_t i = 0; i < nSamples; i++) {
    const Double_t out = filter.Process(preamp.Sample(i));
    const Double_t t = Double_t(i - pad) * kSamplePeriodNs / 1000.0;
    if (t >= tMinUs && (i - pad) % stride == 0)
      graph->SetPoint(point++, t, out);
  }
  return graph;
}

// The raw preamplifier waveform for a single unit pulse.
TGraph *RawPulse(Double_t tauPreUs, Double_t tMinUs, Double_t tMaxUs,
                 Int_t stride = 1) {
  const Long64_t pad = Long64_t(2000.0 / kSamplePeriodNs);
  const Long64_t nSamples =
      pad + Long64_t(tMaxUs * 1000.0 / kSamplePeriodNs) + 1;
  Preamp preamp(tauPreUs * 1000.0, kSamplePeriodNs, kSamplePeriodNs);
  preamp.AddEvent(Double_t(pad) * kSamplePeriodNs, 1.0);
  TGraph *graph = new TGraph();
  Int_t point = 0;
  for (Long64_t i = 0; i < nSamples; i++) {
    const Double_t value = preamp.Sample(i);
    const Double_t t = Double_t(i - pad) * kSamplePeriodNs / 1000.0;
    if (t >= tMinUs && (i - pad) % stride == 0)
      graph->SetPoint(point++, t, value);
  }
  return graph;
}

// (Fig. 1a) The matched trapezoid: flat top of height A, back to zero after
// (k + l)T, with the energy read at t_m.
void TrapezoidFigure() {
  const Double_t tauPre = 10.0;
  const Double_t kUs = kTrapRiseNs / 1000.0;
  const Double_t lUs = (kTrapRiseNs + kTrapFlatNs) / 1000.0;
  const Double_t endUs = kUs + lUs;
  const Double_t tm = kUs + kPeakingDelayNs / 1000.0;
  const Double_t tMin = -1.0, tMax = 20.0;
  const Double_t yLo = -0.1, yHi = 1.95;
  const Double_t markTop = 1.15;

  std::vector<Double_t> one(1, 0.0);
  TGraph *out = StepResponse(tauPre, tauPre, one, tMin, tMax, 8);
  TGraph *in = RawPulse(tauPre, tMin, tMax, 8);

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, tMin, yLo, tMax, yHi, "Time [#mus]", "Amplitude [a.u.]");

  TimeMark(0.0, yLo, markTop);
  TimeMark(kUs, yLo, markTop);
  TimeMark(lUs, yLo, markTop);
  TimeMark(endUs, yLo, markTop);
  TimeMark(tm, yLo, 1.0);

  PlottingUtils::ConfigureGraph(in, kGray + 2);
  in->SetLineStyle(2);
  in->Draw("L SAME");
  PlottingUtils::ConfigureGraph(out, kMatchedColour);
  out->Draw("L SAME");
  canvas->RedrawAxis();

  Note(0.5 * kUs, 1.18, "rise", 21);
  Note(0.5 * (kUs + lUs), 1.18, "flat top", 21);
  Note(0.5 * (lUs + endUs), 1.18, "fall", 21);
  TimeArrow(0.0, kUs, 1.40, 1.46, "#it{kT}");
  TimeArrow(lUs, endUs, 1.40, 1.46, "#it{kT}");
  TimeArrow(0.0, lUs, 1.58, 1.64, "#it{lT}: rise + flat top");
  TimeArrow(0.0, endUs, 1.76, 1.82, "(#it{k}+#it{l})#it{T}: back to zero");

  TGraph *read = new TGraph(1);
  read->SetPoint(0, tm, out->Eval(tm));
  read->SetMarkerStyle(20);
  read->SetMarkerSize(2.0);
  read->SetMarkerColor(kBlack);
  read->Draw("P SAME");
  TLine *leader = new TLine(tm + 0.15, 1.04, endUs + 0.6, 1.10);
  leader->SetLineWidth(PlottingUtils::GetLineWidth());
  leader->Draw();
  Note(endUs + 0.8, 1.10, "energy read at #it{t}_{m}", 12);
  Note(endUs + 2.0, 0.12, "no tail:  #tau_{pz} = #tau_{pre}", 12);

  std::vector<TString> labels;
  labels.push_back("Preamplifier Output #it{v}(#it{t})");
  labels.push_back("Filter Output #hat{#it{s}}(#it{t})");
  TLegend *legend = Legend(labels, 0.70, 0.86);
  legend->AddEntry(in, labels[0], "l");
  legend->AddEntry(out, labels[1], "l");
  legend->Draw();
  PanelLabel("(a)");

  PlottingUtils::SaveFigure(canvas, "trapezoid", "", PlotSaveOptions::kLINEAR);
  delete canvas;
}

// (Fig. 1b) The same preamp through matched, under- and over-compensated
// pole--zero settings; the tails after (k + l)T are the point.
void PoleZeroFigure() {
  const Double_t tauPre = 10.0;
  const Double_t kUs = kTrapRiseNs / 1000.0;
  const Double_t lUs = (kTrapRiseNs + kTrapFlatNs) / 1000.0;
  const Double_t endUs = kUs + lUs;
  const Double_t tauPz[3] = {20.0, 5.0, 10.0};
  const Int_t colour[3] = {kRightColour, kLeftColour, kMatchedColour};
  const Double_t tMin = -1.0, tMax = 45.0;
  const Double_t yLo = -0.45, yHi = 0.45;
  std::vector<Double_t> one(1, 0.0);

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, tMin, yLo, tMax, yHi, "Time [#mus]",
        "Filter Output #hat{#it{s}}(t) [a.u.]");
  Reference(tMin, tMax, 0.0);
  TimeMark(endUs, yLo, yHi);

  TGraph *graphs[3];
  std::vector<TString> labels;
  for (Int_t v = 0; v < 3; v++) {
    graphs[v] = StepResponse(tauPre, tauPz[v], one, tMin, tMax, 6);
    PlottingUtils::ConfigureGraph(graphs[v], colour[v]);
    graphs[v]->Draw("L SAME");
    labels.push_back(
        Form("#tau_{pz} = %.0f #mus%s,  %s", tauPz[v],
             tauPz[v] == tauPre ? " = #tau_{pre}" : "",
             TailLabel(TailAmplitude(tauPre, tauPz[v], kUs, lUs)).Data()));
  }
  canvas->RedrawAxis();

  Note(endUs + 1.55, 0.31, "over-compensated", 10, colour[1]);
  Note(endUs + 4, -0.14, "under-compensated", 13, colour[0]);
  Note(endUs + 0.4, yLo + 0.06, "(#it{k}+#it{l})#it{T}", 10);

  TLegend *legend = Legend(labels, 0.16, 0.33);
  const Int_t order[3] = {2, 0, 1};
  for (Int_t j = 0; j < 3; j++) {
    const Int_t v = order[j];
    legend->AddEntry(graphs[v], labels[v], "l");
  }
  legend->Draw();
  PanelLabel("(b)");

  PlottingUtils::SaveFigure(canvas, "polezero", "", PlotSaveOptions::kLINEAR);
  delete canvas;
}

// (Fig. 2) One pulse on the right chain (tau_pre = 10 us) through the 37Cl
// settings (tau_pz = 20 us).
void SinglePulseFigure() {
  const Double_t tauPre = 10.0;
  const Double_t kUs = kTrapRiseNs / 1000.0;
  const Double_t lUs = (kTrapRiseNs + kTrapFlatNs) / 1000.0;
  const Double_t endUs = kUs + lUs;
  const Double_t tm = kUs + kPeakingDelayNs / 1000.0;
  const Double_t tf = kBaselineHoldoffNs / 1000.0;
  const Double_t tMin = -1.0, tMax = 20.0;
  const Double_t yLo = -0.6, yHi = 2.0;
  const Double_t markTop = 1.15;
  const Double_t tail = TailAmplitude(tauPre, kPoleZeroUs, kUs, lUs);

  std::vector<Double_t> one(1, 0.0);
  TGraph *out = StepResponse(tauPre, kPoleZeroUs, one, tMin, tMax, 8);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, tMin, yLo, tMax, yHi, "Time [#mus]",
        "Filter Output #hat{#it{s}}(t) [a.u.]");
  Reference(tMin, tMax, 0.0);
  TimeMark(0.0, yLo, markTop);
  TimeMark(kUs, yLo, markTop);
  TimeMark(lUs, yLo, markTop);
  TimeMark(endUs, yLo, markTop);
  PlottingUtils::ConfigureGraph(out, kRightColour);
  out->Draw("L SAME");
  canvas->RedrawAxis();

  Note(0.5 * kUs, 1.18, "rise", 21);
  Note(0.5 * (kUs + lUs), 1.18, "flat top", 21);
  Note(0.5 * (lUs + endUs), 1.18, "fall", 21);
  TimeArrow(0.0, kUs, 1.42, 1.48, "#it{kT}");
  TimeArrow(kUs, tm, 1.42, 1.48, "#it{t}_{p}");
  TimeArrow(0.0, tm, 1.70, 1.76, "#it{t}_{m} = #it{kT} + #it{t}_{p}");
  TimeMark(tm, yLo, out->Eval(tm));
  TGraph *read = new TGraph(1);
  read->SetPoint(0, tm, out->Eval(tm));
  read->SetMarkerStyle(20);
  read->SetMarkerSize(2.0);
  read->SetMarkerColor(kBlack);
  read->Draw("P SAME");
  TLine *leader = new TLine(tm + 0.15, out->Eval(tm) + 0.03, endUs + 1.1, 0.95);
  leader->SetLineWidth(PlottingUtils::GetLineWidth());
  leader->Draw();
  Note(endUs + 1.3, 0.95, "pulse height measured at #it{t}_{m}", 12);
  TimeArrow(0.0, tf, -0.50, -0.44, "baseline frozen for #it{t}_{f}");
  Note(endUs + 1.5, -0.27,
       Form("tail:  #hat{#it{s}} = #minus #it{c} "
            "e^{#minus#kern[0.15]{#it{t}}/#tau_{pre}},  %s",
            TailLabel(tail).Data()),
       13, kRightColour);
  PlottingUtils::AddText("#tau_{pre} = 10 #mus,  #tau_{pz} = 20 #mus", 0.87,
                         0.84);

  PlottingUtils::SaveFigure(canvas, "single", "", PlotSaveOptions::kLINEAR);
  delete canvas;
}

// A second pulse 6 us after the first, inside the 7 us baseline freeze. The
// second read point is lower than for an isolated pulse by the first pulse's
// tail.
void PileupPairFigure() {
  const Double_t tauPre = 10.0;
  const Double_t kUs = kTrapRiseNs / 1000.0;
  const Double_t lUs = (kTrapRiseNs + kTrapFlatNs) / 1000.0;
  const Double_t tm = kUs + kPeakingDelayNs / 1000.0;
  const Double_t dt = 6.0;
  const Double_t tMin = -1.0, tMax = 20.0;
  const Double_t yLo = -0.6, yHi = 2.0;
  const Double_t markTop = 1.15;

  std::vector<Double_t> one(1, 0.0);
  std::vector<Double_t> two(2);
  two[0] = 0.0;
  two[1] = dt;
  std::vector<Double_t> second(1, dt);

  TGraph *pair = StepResponse(tauPre, kPoleZeroUs, two, tMin, tMax, 8);
  TGraph *first = StepResponse(tauPre, kPoleZeroUs, one, tMin, tMax, 8);
  TGraph *alone = StepResponse(tauPre, kPoleZeroUs, second, tMin, tMax, 8);

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, tMin, yLo, tMax, yHi, "Time [#mus]",
        "Filter Output #hat{#it{s}}(t) [a.u.]");
  PlottingUtils::ConfigureGraph(pair, kRightColour);
  pair->Draw("L SAME");
  PlottingUtils::ConfigureGraph(first, kGray + 2);
  first->SetLineStyle(2);
  first->Draw("L SAME");
  PlottingUtils::ConfigureGraph(alone, kLeftColour);
  alone->SetLineStyle(2);
  alone->Draw("L SAME");
  canvas->RedrawAxis();

  TimeMark(0.0, yLo, markTop);
  TimeMark(dt, yLo, markTop);
  TimeMark(dt + tm, yLo, markTop);
  Note(0.1, 1.2, "pulse 1", 21);
  Note(dt, 1.2, "pulse 2", 21);
  Note(dt + tm, 1.17, "#it{E}_{2} read", 21);
  TimeArrow(0.0, dt, 1.42, 1.48, "#Delta#it{t}");
  TimeArrow(dt, dt + tm, 1.42, 1.48, "#it{t}_{m}");
  TimeArrow(0.0, dt + tm, 1.72, 1.78, "#Delta#it{t} + #it{t}_{m}");
  TimeArrow(0.0, kBaselineHoldoffNs / 1000.0, -0.50, -0.44,
            "baseline frozen for #it{t}_{f}");

  const Double_t yAlone = alone->Eval(dt + tm);
  const Double_t yPair = pair->Eval(dt + tm);
  const Double_t yTail = first->Eval(dt + tm);
  TLine *loss = new TLine(dt + tm, yPair, dt + tm, yAlone);
  loss->SetLineWidth(3 * PlottingUtils::GetLineWidth());
  loss->SetLineColor(kBlack);
  loss->Draw();
  TLine *tailBar = new TLine(dt + tm, yTail, dt + tm, 0.0);
  tailBar->SetLineWidth(3 * PlottingUtils::GetLineWidth());
  tailBar->SetLineColor(kBlack);
  tailBar->Draw();
  TGraph *reads = new TGraph(3);
  reads->SetPoint(0, dt + tm, yAlone);
  reads->SetPoint(1, dt + tm, yPair);
  reads->SetPoint(2, dt + tm, yTail);
  reads->SetMarkerStyle(20);
  reads->SetMarkerSize(2.0);
  reads->SetMarkerColor(kBlack);
  reads->Draw("P SAME");

  const Double_t xLossText = 14.2;
  TLine *leadLoss =
      new TLine(dt + tm + 0.15, 0.5 * (yAlone + yPair), xLossText - 0.2, 0.70);
  leadLoss->SetLineWidth(PlottingUtils::GetLineWidth());
  leadLoss->Draw();
  Note(xLossText, 0.70, Form("loss = %.3f", yAlone - yPair), 12);
  const Double_t xTailText = 11.5;
  TLine *leadTail =
      new TLine(dt + tm + 0.15, 0.5 * yTail, xTailText - 0.2, -0.45);
  leadTail->SetLineWidth(PlottingUtils::GetLineWidth());
  leadTail->Draw();
  Note(xTailText, -0.45, Form("tail of pulse 1 = %.3f", yTail), 12);

  std::vector<TString> labels;
  labels.push_back("Measured");
  labels.push_back("Pulse 1");
  labels.push_back("Pulse 2");
  TLegend *legend = Legend(labels, 0.62, 0.87);
  legend->AddEntry(pair, labels[0], "l");
  legend->AddEntry(first, labels[1], "l");
  legend->AddEntry(alone, labels[2], "l");
  legend->Draw();

  PlottingUtils::SaveFigure(canvas, "pileup_pair", "",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
  std::cout << Form("pileup: dt = %.1f us, loss = %.3f, tail = %.3f", dt,
                    yAlone - yPair, yTail)
            << std::endl;
}

// What the figures take from a simulated run's parameters.
struct RunInfo {
  Double_t beamPeak = 1600.0;
  Double_t tauL_us = 0.0;
  Double_t tauR_us = 0.0;
  Double_t tauPz_us = 0.0;
};

Double_t RunParameter(TFile &file, const char *name, Double_t fallback) {
  TParameter<Double_t> *parameter = (TParameter<Double_t> *)file.Get(name);
  return parameter ? parameter->GetVal() : fallback;
}

Bool_t ReadHits(const TString &path, std::vector<SimHit> &hits,
                RunInfo &info) {
  TFile in(path, "READ");
  if (in.IsZombie()) {
    std::cerr << "cannot open " << path << std::endl;
    return kFALSE;
  }
  TTree *tree = (TTree *)in.Get("hits");
  if (!tree) {
    std::cerr << "no hits tree in " << path << std::endl;
    return kFALSE;
  }
  info.beamPeak = RunParameter(in, "beamPeakADC", 1600.0);
  info.tauL_us = RunParameter(in, "tauL_us", 0.0);
  info.tauR_us = RunParameter(in, "tauR_us", 0.0);
  info.tauPz_us = RunParameter(in, "trapDecay_us", 0.0);

  Int_t channel = 0;
  ULong64_t timestamp = 0;
  Double_t energy = 0.0, dtPrev = 0.0, dtNext = 0.0, trueAmp = 0.0,
           prevTrueAmp = 0.0;
  Bool_t pileup = kFALSE, reaction = kFALSE;
  Long64_t eventIndex = 0;
  tree->SetBranchAddress("channel", &channel);
  tree->SetBranchAddress("timestamp", &timestamp);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchAddress("dtPrev", &dtPrev);
  tree->SetBranchAddress("dtNext", &dtNext);
  tree->SetBranchAddress("trueAmp", &trueAmp);
  tree->SetBranchAddress("prevTrueAmp", &prevTrueAmp);
  tree->SetBranchAddress("pileup", &pileup);
  tree->SetBranchAddress("eventIndex", &eventIndex);
  // Runs written before reactions were simulated have no reaction branch.
  if (tree->GetBranch("reaction"))
    tree->SetBranchAddress("reaction", &reaction);
  hits.clear();
  hits.reserve(size_t(tree->GetEntries()));
  for (Long64_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    SimHit hit;
    hit.channel = channel;
    hit.timestamp = timestamp;
    hit.energy = energy;
    hit.dtPrev = dtPrev;
    hit.dtNext = dtNext;
    hit.trueAmp = trueAmp;
    hit.prevTrueAmp = prevTrueAmp;
    hit.pileup = pileup;
    hit.eventIndex = eventIndex;
    hit.reaction = reaction;
    hits.push_back(hit);
  }
  in.Close();
  return kTRUE;
}

// Pulse-height spectra of one simulated run. Non-pileup hits enter the
// histograms.
void DrawSpectrum(const std::vector<SimHit> &hits,
                  const std::vector<Double_t> &energies, Double_t beamPeak,
                  const TString &label, const TString &name,
                  const char *panel) {
  TH1D *left = new TH1D("spectrum_L", ";Pulse Height [ADC];Counts", 300, 0.0,
                        2.5 * beamPeak);
  TH1D *right = new TH1D("spectrum_R", ";Pulse Height [ADC];Counts", 300, 0.0,
                         2.5 * beamPeak);
  left->SetDirectory(nullptr);
  right->SetDirectory(nullptr);

  for (Int_t i = 0; i < Int_t(hits.size()); i++) {
    if (hits[i].pileup)
      continue;
    TH1D *hist = hits[i].channel == 0 ? left : right;
    hist->Fill(energies[i]);
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  // Five-digit counts need more room than the style's left margin. Drawing
  // applies the style, so the title offset has to be set afterwards.
  canvas->SetLeftMargin(0.18);
  const Double_t yMax =
      1.45 * TMath::Max(left->GetMaximum(), right->GetMaximum());
  left->SetMaximum(yMax);
  PlottingUtils::ConfigureAndDrawHistogram(left, kLeftColour,
                                           ";Pulse Height [ADC];Counts");
  left->GetYaxis()->SetTitleOffset(1.4);
  PlottingUtils::ConfigureHistogram(right, kRightColour);
  right->Draw("HIST SAME");

  std::vector<TString> labels;
  labels.push_back("#tau_{pz} = #tau_{pre}");
  labels.push_back("#tau_{pz} #neq #tau_{pre}");
  TLegend *legend = Legend(labels, 0.74, 0.87);
  legend->AddEntry(left, labels[0], "l");
  legend->AddEntry(right, labels[1], "l");
  legend->Draw();
  PlottingUtils::AddText(label, 0.88, 0.92);
  PanelLabel(panel);
  PlottingUtils::SaveFigure(canvas, name, "", PlotSaveOptions::kLINEAR);
  delete canvas;
}

// The PulseHistory correction of one run, with the kernel measured on the run
// itself by the unmodified MUSIC code.
Bool_t CorrectRun(const std::vector<SimHit> &hits, const TString &label,
                  std::vector<Double_t> &corrected) {
  PulseHistoryAdapter::Adapted adapted = PulseHistoryAdapter::Build(hits);
  std::vector<Int_t> groups = PulseHistory::BuildGroupMap();
  PulseHistory::Result kernel;
  if (!PulseHistory::Measure(adapted.hits, groups, kernel, label)) {
    std::cerr << "PulseHistory::Measure failed for " << label << std::endl;
    return kFALSE;
  }
  kernel.n_corrected = 0;
  kernel.n_clamped = 0;
  PulseHistory::Apply(adapted.hits, groups, kernel);
  corrected = PulseHistoryAdapter::CorrectedEnergies(adapted, hits);
  return kTRUE;
}

void SpectrumFigure(const TString &path, const TString &label) {
  std::vector<SimHit> hits;
  RunInfo info;
  if (!ReadHits(path, hits, info))
    return;

  const char *panel = label == "1kHz" ? "(a)" : "(b)";
  std::vector<Double_t> raw(hits.size(), 0.0);
  for (Int_t i = 0; i < Int_t(hits.size()); i++)
    raw[i] = hits[i].energy;
  DrawSpectrum(hits, raw, info.beamPeak, label,
               Form("spectrum_%s", label.Data()), panel);

  std::vector<Double_t> corrected;
  if (!CorrectRun(hits, label, corrected))
    return;
  DrawSpectrum(hits, corrected, info.beamPeak, label,
               Form("spectrum_%s_corrected", label.Data()), panel);
}

// Standard deviation about the mean, the former analyze program's "RMS".
Double_t Rms(const std::vector<Double_t> &values) {
  if (values.empty())
    return 0.0;
  Double_t sum = 0.0, sum2 = 0.0;
  for (size_t i = 0; i < values.size(); i++) {
    sum += values[i];
    sum2 += values[i] * values[i];
  }
  const Double_t mean = sum / Double_t(values.size());
  const Double_t variance = sum2 / Double_t(values.size()) - mean * mean;
  return variance > 0.0 ? TMath::Sqrt(variance) : 0.0;
}

// The spread of one run, raw and corrected, per chain: non-pileup pulse
// heights within the beam window, as the former analyze program took them.
Bool_t MeasureSpread(const TString &path, const TString &label, Spread &spread,
                     RunInfo &info) {
  std::vector<SimHit> hits;
  if (!ReadHits(path, hits, info))
    return kFALSE;
  std::vector<Double_t> corrected;
  if (!CorrectRun(hits, label, corrected))
    return kFALSE;
  const Double_t lo = kSpreadWindowLo * info.beamPeak;
  const Double_t hi = kSpreadWindowHi * info.beamPeak;
  std::vector<Double_t> raw[2], corr[2];
  for (size_t i = 0; i < hits.size(); i++) {
    const SimHit &hit = hits[i];
    if (hit.pileup)
      continue;
    if (hit.energy > lo && hit.energy < hi)
      raw[hit.channel].push_back(hit.energy);
    if (corrected[i] > lo && corrected[i] < hi)
      corr[hit.channel].push_back(corrected[i]);
  }
  for (Int_t ch = 0; ch < 2; ch++) {
    spread.raw[ch] = Rms(raw[ch]);
    spread.corr[ch] = Rms(corr[ch]);
  }
  std::cout << Form("spread %s: L %.1f -> %.1f, R %.1f -> %.1f ADC",
                    label.Data(), spread.raw[0], spread.corr[0],
                    spread.raw[1], spread.corr[1])
            << std::endl;
  return kTRUE;
}

// ----------------------------------------------------------- data figures
// (Fig. 3) The kernels PulseHistory fitted on run 97.
void KernelFigure(const std::vector<std::vector<Double_t>> &kernels) {
  const char *name[3] = {"Long Strips, Left", "Long Strips, Right",
                         "Short Strips, Left"};
  const Int_t colour[3] = {kLeftColour, kRightColour, kShortLeftColour};
  const Int_t order[3] = {1, 0, 2};
  const Double_t xLo = 1.0, xHi = 350.0, yLo = -0.2, yHi = 0.75;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetLogx();
  TH1F *frame = Frame(canvas, xLo, yLo, xHi, yHi,
                      "Time Between Pulses #Delta#it{t} [#mus]",
                      "Shift / Previous Pulse Height");
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();

  Reference(xLo, xHi, 0.0);

  std::vector<TGraph *> graphs;
  for (Int_t j = 0; j < 3; j++) {
    const std::vector<Double_t> &values = kernels[order[j]];
    TGraph *graph = new TGraph();
    for (Int_t i = 0; i < kNBins; i++)
      graph->SetPoint(i, BinCentre(i), values[i]);
    PlottingUtils::ConfigureGraph(graph, colour[order[j]]);
    graph->SetMarkerStyle(20 + j);
    graph->SetMarkerColor(colour[order[j]]);
    graph->SetMarkerSize(1.2);
    graph->Draw("LP SAME");
    graphs.push_back(graph);
  }
  canvas->RedrawAxis();

  std::vector<TString> labels;
  for (Int_t j = 0; j < 3; j++)
    labels.push_back(name[order[j]]);
  TLegend *legend = Legend(labels, 0.64, 0.86);
  for (Int_t j = 0; j < 3; j++)
    legend->AddEntry(graphs[j], labels[j], "LP");
  legend->Draw();
  PlottingUtils::SaveFigure(canvas, "kernel_run97", "",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
}

// The simulated rms spread of each baseline reading, before and after
// correction, as a complete booktabs table float with its caption and label
// (tab:spread), for the paper to \input (needs siunitx and amsmath). The
// isolated column is the corrected spread at the isolated rate. The caption's
// decay times are those of the runs, which main checks agree; spreads holds
// every reading at every rate.
void SpreadTable(const ReadingSpreads &spreads, const RunInfo &runs,
                 const TString &path) {
  std::ofstream table(path.Data());
  table << "% Written by mains/figures.cpp; edit that, not this file."
        << std::endl;
  table << "\\begin{table}[htb]" << std::endl;
  table << "\\centering" << std::endl;
  table << Form("\\caption{Simulated rms pulse-height spread, in ADC, of "
                "pulses within 0.5--1.5 times the beam peak with no neighbour "
                "within \\SI{5}{\\micro\\second}, before and after correction, "
                "and for isolated pulses (\\SI{0.5}{\\kilo\\hertz}), for three "
                "behaviours of the baseline restorer. The runs used "
                "$\\tau_{\\text{pre}} = \\SI{%g}{\\micro\\second}$ on the "
                "right segments and \\SI{%g}{\\micro\\second} on the left, "
                "with $\\tau_{\\text{pz}} = \\SI{%g}{\\micro\\second}$.}",
                runs.tauR_us, runs.tauL_us, runs.tauPz_us)
        << std::endl;
  table << "\\label{tab:spread}" << std::endl;
  // Scoped by the table; without them the tabular overflows the text width of
  // the elsarticle preprint.
  table << "\\small" << std::endl;
  table << "\\setlength{\\tabcolsep}{4pt}" << std::endl;
  table << "\\begin{tabular}{@{}lccccccc@{}}" << std::endl;
  table << "\\toprule" << std::endl;
  table << "&";
  for (Int_t i = 0; i < kNSpreadRates; i++)
    table << Form(" & \\multicolumn{2}{c}{\\SI{%.0f}{\\kilo\\hertz}}",
                  kSpreadRates[i]);
  table << " \\\\" << std::endl;
  table << "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}"
        << std::endl;
  table << "Baseline & isolated";
  for (Int_t i = 0; i < kNSpreadRates; i++)
    table << " & before & after";
  table << " \\\\" << std::endl;

  const Int_t nSides = 2;
  const Int_t sides[nSides] = {1, 0};
  const char *sideNames[nSides] = {"Long right segments",
                                   "Long left segments"};
  for (Int_t s = 0; s < nSides; s++) {
    table << "\\midrule" << std::endl;
    table << "\\multicolumn{8}{l}{\\textit{" << sideNames[s] << "}} \\\\"
          << std::endl;
    for (Int_t j = 0; j < kNSpreadReadings; j++) {
      const RateSpreads &byRate = spreads.find(kSpreadReadings[j].key)->second;
      table << kSpreadReadings[j].label
            << Form(" & %.1f",
                    byRate.find(kIsolatedRate)->second.corr[sides[s]]);
      for (Int_t i = 0; i < kNSpreadRates; i++) {
        const Spread &spread = byRate.find(kSpreadRates[i])->second;
        table << Form(" & %.1f & %.1f", spread.raw[sides[s]],
                      spread.corr[sides[s]]);
      }
      table << " \\\\" << std::endl;
    }
  }
  table << "\\bottomrule" << std::endl;
  table << "\\end{tabular}" << std::endl;
  table << "\\end{table}" << std::endl;
  std::cout << "wrote " << path << std::endl;
}

TGraphErrors *ProfileGraph(const TProfile *profile) {
  TGraphErrors *graph = new TGraphErrors();
  for (Int_t bin = 1; bin <= profile->GetNbinsX(); bin++) {
    if (profile->GetBinEntries(bin) <= 0)
      continue;
    const Int_t n = graph->GetN();
    graph->SetPoint(n, TMath::Power(10.0, profile->GetBinCenter(bin)),
                    profile->GetBinContent(bin));
    graph->SetPointError(n, 0.0, profile->GetBinError(bin));
  }
  return graph;
}

// (Fig. 5) Run 97 long-right profile before and after correction.
void DeviationFigure(const TString &file) {
  TFile in(file, "READ");
  const TProfile *beforeProfile = (const TProfile *)in.Get("before_long_right");
  const TProfile *afterProfile = (const TProfile *)in.Get("after_long_right");
  if (!beforeProfile || !afterProfile) {
    std::cerr << "no long-right profiles in " << file << std::endl;
    return;
  }
  TGraphErrors *before = ProfileGraph(beforeProfile);
  TGraphErrors *after = ProfileGraph(afterProfile);
  const Double_t xLo = 1.0, xHi = 350.0, yLo = -150.0, yHi = 300.0;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetLogx();
  TH1F *frame = Frame(canvas, xLo, yLo, xHi, yHi,
                      "Time Between Pulses #Delta#it{t} [#mus]",
                      "Pulse Height #minus Average [ADC]");
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();
  Reference(xLo, xHi, 0.0);

  PlottingUtils::ConfigureGraph(before, kRightColour);
  before->SetLineStyle(2);
  before->SetMarkerStyle(24);
  before->SetMarkerColor(kRightColour);
  before->SetMarkerSize(1.2);
  before->Draw("LP SAME");
  PlottingUtils::ConfigureGraph(after, kRightColour);
  after->SetMarkerStyle(20);
  after->SetMarkerColor(kRightColour);
  after->SetMarkerSize(1.2);
  after->Draw("LP SAME");
  canvas->RedrawAxis();

  std::vector<TString> labels;
  labels.push_back("Before");
  labels.push_back("After");
  TLegend *legend = Legend(labels, 0.70, 0.86, "", kTRUE);
  legend->AddEntry(before, labels[0], "LP");
  legend->AddEntry(after, labels[1], "LP");
  legend->Draw();
  PlottingUtils::SaveFigure(canvas, "run97_deviation_R", "",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
}

Bool_t HaveInput(const TString &file, const char *figure) {
  if (!gSystem->AccessPathName(file))
    return kTRUE;
  std::cout << "skipping " << figure << ": no " << file << std::endl;
  return kFALSE;
}

// PlotSaveFormat::kPNGAndPDF writes both side by side; the PNGs belong in
// their own png/ directory.
void MovePngs(const TString &dir) {
  const TString pngDir = dir + "/png";
  gSystem->mkdir(pngDir, kTRUE);
  void *listing = gSystem->OpenDirectory(dir);
  if (!listing)
    return;
  const char *entry;
  while ((entry = gSystem->GetDirEntry(listing))) {
    const TString name = entry;
    if (name.EndsWith(".png"))
      gSystem->Rename(dir + "/" + name, pngDir + "/" + name);
  }
  gSystem->FreeDirectory(listing);
}

} // namespace

int main(int argc, char **argv) {
  TString out = "doc/figs", run = "run97_c000";
  TString log = TString(gSystem->HomeDirectory()) +
                "/MUSIC/analysis/37Cl/pipeline_fused.log";
  TString spreadDir = "output/spread";
  TString run97 = "output/run97/deviation.root";
  TString long1 = "output/figures/long_1kHz.root";
  TString long40 = "output/figures/long_40kHz.root";
  PlotSaveFormat format = PlotSaveFormat::kPNGAndPDF;
  for (Int_t i = 1; i < argc; i++) {
    TString argument = argv[i];
    if (argument.BeginsWith("--out="))
      out = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--log="))
      log = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--run="))
      run = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--spread="))
      spreadDir = argument(9, argument.Length() - 9);
    else if (argument.BeginsWith("--run97="))
      run97 = argument(8, argument.Length() - 8);
    else if (argument.BeginsWith("--long1="))
      long1 = argument(8, argument.Length() - 8);
    else if (argument.BeginsWith("--long40="))
      long40 = argument(9, argument.Length() - 9);
    else if (argument == "--format=png")
      format = PlotSaveFormat::kPNG;
    else if (argument == "--format=pdf")
      format = PlotSaveFormat::kPDF;
  }
  InitUtils::SetROOTPreferences(format, out, ".", kFALSE);

  TrapezoidFigure();
  PoleZeroFigure();
  SinglePulseFigure();
  PileupPairFigure();
  if (HaveInput(long1, "spectrum_1kHz"))
    SpectrumFigure(long1, "1kHz");
  if (HaveInput(long40, "spectrum_40kHz"))
    SpectrumFigure(long40, "40kHz");

  if (HaveInput(log, "kernel_run97")) {
    std::vector<std::vector<Double_t>> kernels = ReadKernels(log, run);
    if (kernels.size() < 3) {
      std::cerr << "found " << kernels.size() << " kernels for " << run
                << " in " << log << std::endl;
    } else {
      KernelFigure(kernels);
    }
  }

  if (format != PlotSaveFormat::kPNG) {
    std::vector<TString> spreadFiles, spreadLabels;
    std::vector<Double_t> allRates(1, kIsolatedRate);
    for (Int_t i = 0; i < kNSpreadRates; i++)
      allRates.push_back(kSpreadRates[i]);
    Bool_t haveSpreads = kTRUE;
    for (Int_t j = 0; j < kNSpreadReadings && haveSpreads; j++) {
      for (size_t i = 0; i < allRates.size() && haveSpreads; i++) {
        const TString label =
            Form("%s_r%g", kSpreadReadings[j].key, allRates[i]);
        spreadLabels.push_back(label);
        spreadFiles.push_back(spreadDir + "/" + label + ".root");
        haveSpreads = HaveInput(spreadFiles.back(), "spread_table");
      }
    }
    ReadingSpreads spreads;
    RunInfo runs;
    for (Int_t j = 0; j < kNSpreadReadings && haveSpreads; j++) {
      for (size_t i = 0; i < allRates.size() && haveSpreads; i++) {
        const size_t n = j * allRates.size() + i;
        RunInfo info;
        haveSpreads = MeasureSpread(spreadFiles[n], spreadLabels[n],
                                    spreads[kSpreadReadings[j].key][allRates[i]],
                                    info);
        if (n == 0)
          runs = info;
        else if (info.tauL_us != runs.tauL_us ||
                 info.tauR_us != runs.tauR_us ||
                 info.tauPz_us != runs.tauPz_us) {
          std::cout << "skipping spread_table: " << spreadFiles[n]
                    << " has other decay times than " << spreadFiles[0]
                    << std::endl;
          haveSpreads = kFALSE;
        }
      }
    }
    if (haveSpreads)
      SpreadTable(spreads, runs, out + "/spread_table.tex");
  }

  if (HaveInput(run97, "run97_deviation_R"))
    DeviationFigure(run97);
  if (format == PlotSaveFormat::kPNGAndPDF)
    MovePngs(out);
  return 0;
}
