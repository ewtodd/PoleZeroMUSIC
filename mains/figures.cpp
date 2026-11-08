#include "Constants.hpp"
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
#include <TVectorD.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// src/SimDatasetConfig.cpp: the kernel's amplitude bands, read by
// PulseHistory on every Measure.
namespace Constants {
void SetSimPulseHistoryAmpBins(Int_t n);
} // namespace Constants

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
// The group kernels scripts/run97_figure.sh saves (long left, long right,
// short left), each the shift per unit previous pulse height in kNBins bins;
// empty if the file has none.
std::vector<std::vector<Double_t>> ReadKernels(const TString &file) {
  std::vector<std::vector<Double_t>> kernels;
  TFile in(file, "READ");
  const Int_t nGroups = 3;
  const char *names[nGroups] = {"kernel_long_left", "kernel_long_right",
                                "kernel_short_left"};
  for (Int_t i = 0; i < nGroups; i++) {
    const TVectorD *values = (const TVectorD *)in.Get(names[i]);
    if (!values || values->GetNrows() != kNBins)
      return std::vector<std::vector<Double_t>>();
    std::vector<Double_t> kernel(kNBins);
    for (Int_t b = 0; b < kNBins; b++)
      kernel[b] = (*values)[b];
    kernels.push_back(kernel);
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

// The spread table's runs (scripts/spread_inputs.sh): two behaviours of the
// baseline restorer, each at the isolated rate and the three table rates.
// Ideal subtracts the true baseline. Conditional restart is the firmware's
// with EnergyFilterLFLimitation on, as on every channel of the 37Cl runs
// (data_raw/music_exp1915_097XSetting_66222.dat): after the freeze, averaging
// resumes only once the trapezoid is back near the baseline. CAEN does not say
// how near. The simulation uses 80 ADC held for 2 us, its default; at
// tauR = 10 us, run 97 favours an earlier restart, about 40 ADC for 1 us.
struct SpreadReading {
  const char *key;
  const char *label;
};
const Int_t kNSpreadReadings = 2;
const SpreadReading kSpreadReadings[kNSpreadReadings] = {
    {"ideal", "ideal"}, {"wait", "conditional restart"}};
// The conditional restart as scripts/spread_inputs.sh simulates it
// (--baselineRestartThrADC=80 and the default 250-sample hold).
const Double_t kRestartThrADC = 80.0;
const Double_t kRestartHoldUs = 2.0;
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

// A legend inside the frame, sized to its labels and anchored a small inset
// from the frame's right (or left) edge, wherever the pad's margins put it.
TLegend *Legend(const std::vector<TString> &labels, Double_t y1, Double_t y2,
                const TString &header = "", Bool_t left = kFALSE) {
  const Double_t width =
      PlottingUtils::LegendWidthPx(labels, header) / gPad->GetWw();
  const Double_t inset = 0.02;
  const Double_t x2 = left ? gPad->GetLeftMargin() + inset + width
                           : 1.0 - gPad->GetRightMargin() - inset;
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

Bool_t ReadHits(const TString &path, std::vector<SimHit> &hits, RunInfo &info) {
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
// The correction fits one kernel per band of the previous pulse height
// (PULSE_HISTORY_AMP_BINS); production uses one.
Bool_t MeasureSpread(const TString &path, const TString &label, Spread &spread,
                     RunInfo &info, Int_t bands = 1) {
  std::vector<SimHit> hits;
  if (!ReadHits(path, hits, info))
    return kFALSE;
  const Int_t productionBands = Constants::cfg.PULSE_HISTORY_AMP_BINS;
  Constants::SetSimPulseHistoryAmpBins(bands);
  std::vector<Double_t> corrected;
  const Bool_t ok = CorrectRun(hits, label, corrected);
  Constants::SetSimPulseHistoryAmpBins(productionBands);
  if (!ok)
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
                    label.Data(), spread.raw[0], spread.corr[0], spread.raw[1],
                    spread.corr[1])
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
  TH1F *frame = Frame(canvas, xLo, yLo, xHi, yHi, "Time Between Pulses [#mus]",
                      "Shift / Earlier Pulse Height");
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
  table << "\\begin{table}[htbp]" << std::endl;
  table << "\\centering" << std::endl;
  table << Form("\\caption{Simulated rms pulse-height spread, in ADC, of "
                "pulses within 0.5--1.5 times the beam peak with no neighbour "
                "within \\SI{5}{\\micro\\second}, before and after correction, "
                "and for isolated pulses (\\SI{0.5}{\\kilo\\hertz}), for two "
                "behaviours of the baseline restorer: ideal, the true "
                "baseline is subtracted; and conditional restart, the "
                "firmware's with its LF limitation on, as in the $^{37}$Cl "
                "runs, in which averaging resumes after the "
                "\\SI{%g}{\\micro\\second} freeze only once the trapezoid is "
                "back near the baseline, simulated as within \\SI{%g}{ADC} "
                "for \\SI{%g}{\\micro\\second}. The runs used "
                "$\\tau_{\\text{pre}} = \\SI{%g}{\\micro\\second}$ on the "
                "Long Strips, Right and \\SI{%g}{\\micro\\second} on the Long "
                "Strips, Left, with $\\tau_{\\text{pz}} = "
                "\\SI{%g}{\\micro\\second}$.}",
                kBaselineHoldoffNs / 1000.0, kRestartThrADC, kRestartHoldUs,
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
  const char *sideNames[nSides] = {"Long Strips, Right", "Long Strips, Left"};
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

// The runs of scripts/reaction_inputs.sh that the reaction table shows.
struct ReactionRun {
  const char *key;
  Double_t rate_kHz;
};
const Int_t kNReactionRuns = 3;
const ReactionRun kReactionRuns[kNReactionRuns] = {
    {"ideal", 44.0}, {"ideal", 80.0}, {"wait", 80.0}};
// Beam particles, reacted ones, and beam particles whose previous pulse on the
// same end was a reacted one.
const Int_t kNReactionClasses = 3;

// Spread of the height about the chain's gain times the true amplitude, in
// ADC, before (raw) and after (corr) correction: [chain][class].
struct ReactionScore {
  Double_t raw[2][kNReactionClasses];
  Double_t corr[2][kNReactionClasses];
};

Double_t SpreadOf(Double_t n, Double_t sum, Double_t sum2) {
  if (!(n > 0.0))
    return 0.0;
  const Double_t mean = sum / n;
  return TMath::Sqrt(TMath::Max(0.0, sum2 / n - mean * mean));
}

const char *ReadingLabel(const char *key) {
  for (Int_t j = 0; j < kNSpreadReadings; j++)
    if (TString(kSpreadReadings[j].key) == key)
      return kSpreadReadings[j].label;
  return key;
}

// The correction of one reaction run against truth. Only particles with no
// other particle, recorded or not, within the trapezoid's rise and flat top
// enter: the pileup flag knows only recorded triggers, and a particle lost in
// the trigger holdoff merges into the one before. The gains come from the
// beam particles.
Bool_t ScoreReactions(const TString &path, const TString &label,
                      ReactionScore &score, RunInfo &info, Double_t &fraction) {
  std::vector<SimHit> hits;
  if (!ReadHits(path, hits, info))
    return kFALSE;
  std::vector<Double_t> times;
  {
    TFile in(path, "READ");
    TTree *truth = (TTree *)in.Get("truth");
    if (!truth) {
      std::cerr << "no truth tree in " << path << std::endl;
      return kFALSE;
    }
    Double_t time_ns = 0.0;
    truth->SetBranchAddress("time_ns", &time_ns);
    for (Long64_t e = 0; e < truth->GetEntries(); e++) {
      truth->GetEntry(e);
      times.push_back(time_ns);
    }
    fraction = RunParameter(in, "reactionFraction", 0.0);
  }
  std::vector<Double_t> corrected;
  if (!CorrectRun(hits, label, corrected))
    return kFALSE;

  const Double_t isolation_ns = kTrapRiseNs + kTrapFlatNs;
  std::vector<Bool_t> alone(hits.size(), kFALSE);
  Double_t sumRaw[2] = {0.0, 0.0}, sumCorr[2] = {0.0, 0.0},
           sumTrue[2] = {0.0, 0.0};
  for (size_t i = 0; i < hits.size(); i++) {
    const Long64_t e = hits[i].eventIndex;
    alone[i] = (e == 0 || times[e] - times[e - 1] >= isolation_ns) &&
               (e + 1 >= Long64_t(times.size()) ||
                times[e + 1] - times[e] >= isolation_ns);
    if (!alone[i] || hits[i].reaction)
      continue;
    sumRaw[hits[i].channel] += hits[i].energy;
    sumCorr[hits[i].channel] += corrected[i];
    sumTrue[hits[i].channel] += hits[i].trueAmp;
  }
  Double_t n[2][kNReactionClasses] = {{0.0}};
  Double_t sum[2][kNReactionClasses][2] = {{{0.0}}};
  Double_t sum2[2][kNReactionClasses][2] = {{{0.0}}};
  Bool_t lastReacted[2] = {kFALSE, kFALSE};
  for (size_t i = 0; i < hits.size(); i++) {
    const SimHit &hit = hits[i];
    const Int_t ch = hit.channel;
    const Bool_t afterReacted = lastReacted[ch];
    lastReacted[ch] = hit.reaction;
    if (!alone[i] || !(sumTrue[ch] > 0.0))
      continue;
    const Int_t c = hit.reaction ? 1 : (afterReacted ? 2 : 0);
    const Double_t residual[2] = {
        hit.energy - sumRaw[ch] / sumTrue[ch] * hit.trueAmp,
        corrected[i] - sumCorr[ch] / sumTrue[ch] * hit.trueAmp};
    n[ch][c]++;
    for (Int_t k = 0; k < 2; k++) {
      sum[ch][c][k] += residual[k];
      sum2[ch][c][k] += residual[k] * residual[k];
    }
  }
  for (Int_t ch = 0; ch < 2; ch++)
    for (Int_t c = 0; c < kNReactionClasses; c++) {
      score.raw[ch][c] = SpreadOf(n[ch][c], sum[ch][c][0], sum2[ch][c][0]);
      score.corr[ch][c] = SpreadOf(n[ch][c], sum[ch][c][1], sum2[ch][c][1]);
    }
  std::cout << Form("reactions %s: right beam %.1f -> %.1f, reacted %.1f -> "
                    "%.1f (n %.0f); left beam %.1f -> %.1f, reacted %.1f -> "
                    "%.1f (n %.0f) ADC",
                    label.Data(), score.raw[1][0], score.corr[1][0],
                    score.raw[1][1], score.corr[1][1], n[1][1], score.raw[0][0],
                    score.corr[0][0], score.raw[0][1], score.corr[0][1],
                    n[0][1])
            << std::endl;
  return kTRUE;
}

// The correction of reacted particles against truth, as a complete booktabs
// table float with its caption and label (tab:reactions), for the paper to
// \input (needs siunitx and amsmath).
void ReactionTable(const std::vector<ReactionScore> &scores,
                   const RunInfo &runs, Double_t fraction,
                   const TString &path) {
  std::ofstream table(path.Data());
  table << "% Written by mains/figures.cpp; edit that, not this file."
        << std::endl;
  table << "\\begin{table}[htbp]" << std::endl;
  table << "\\centering" << std::endl;
  table << Form("\\caption{Spread, in ADC, of the simulated pulse height "
                "about the gain times the true amplitude, before and after "
                "correction, for beam particles, reacted particles, and beam "
                "particles whose previous pulse on the same strip was a "
                "reacted one. Particles with another particle within "
                "\\SI{5}{\\micro\\second} are left out; the gain is that of "
                "the beam particles. \\SI{%g}{\\percent} of the particles "
                "reacted, raised from about $2\\times10^{-5}$ for statistics; "
                "their deposits in the Long Strips are those of $(\\alpha,n)$ "
                "and $(\\alpha,p)$ reactions in the $^{37}$Cl simulation, "
                "weighted by their TALYS cross sections. The baseline "
                "behaviours are those of Table~\\ref{tab:spread}. The runs "
                "used "
                "$\\tau_{\\text{pre}} = \\SI{%g}{\\micro\\second}$ on the "
                "Long Strips, Right and \\SI{%g}{\\micro\\second} on the Long "
                "Strips, Left, with $\\tau_{\\text{pz}} = "
                "\\SI{%g}{\\micro\\second}$.}",
                100.0 * fraction, runs.tauR_us, runs.tauL_us, runs.tauPz_us)
        << std::endl;
  table << "\\label{tab:reactions}" << std::endl;
  // Scoped by the table; without them the tabular overflows the text width of
  // the elsarticle preprint.
  table << "\\small" << std::endl;
  table << "\\setlength{\\tabcolsep}{4pt}" << std::endl;
  table << "\\begin{tabular}{@{}llcccccc@{}}" << std::endl;
  table << "\\toprule" << std::endl;
  table << "& & \\multicolumn{2}{c}{beam} & \\multicolumn{2}{c}{reacted} & "
           "\\multicolumn{2}{c}{after reacted} \\\\"
        << std::endl;
  table << "\\cmidrule(lr){3-4} \\cmidrule(lr){5-6} \\cmidrule(lr){7-8}"
        << std::endl;
  table << "Baseline & rate & before & after & before & after & before & "
           "after \\\\"
        << std::endl;
  const Int_t nSides = 2;
  const Int_t sides[nSides] = {1, 0};
  const char *sideNames[nSides] = {"Long Strips, Right", "Long Strips, Left"};
  for (Int_t s = 0; s < nSides; s++) {
    table << "\\midrule" << std::endl;
    table << "\\multicolumn{8}{l}{\\textit{" << sideNames[s] << "}} \\\\"
          << std::endl;
    for (Int_t j = 0; j < kNReactionRuns; j++) {
      table << ReadingLabel(kReactionRuns[j].key)
            << Form(" & \\SI{%.0f}{\\kilo\\hertz}", kReactionRuns[j].rate_kHz);
      for (Int_t c = 0; c < kNReactionClasses; c++)
        table << Form(" & %.1f & %.1f", scores[j].raw[sides[s]][c],
                      scores[j].corr[sides[s]][c]);
      table << " \\\\" << std::endl;
    }
  }
  table << "\\bottomrule" << std::endl;
  table << "\\end{tabular}" << std::endl;
  table << "\\end{table}" << std::endl;
  std::cout << "wrote " << path << std::endl;
}

// The band table: the correction with one kernel per band of the previous
// pulse height, bands [0, 0.5), [0.5, 1.5), ..., [kNBands - 1.5, inf) times
// the beam peak, against the production single kernel. The spread table's
// baseline readings at its three rates, and run 97.
const Int_t kNBands = 4;

// What scripts/run97_figure.sh measures for the band table: the long groups'
// spread before correction, with one band and with `bands` bands, in percent
// of each channel's beam peak, [0] left and [1] right as the simulated
// chains; and the long ends' trigger rate.
struct DataBands {
  Double_t spread[2][3];
  Int_t bands = 0;
  Double_t rate_kHz = 0.0;
};

Bool_t ReadDataBands(const TString &file, DataBands &data) {
  TFile in(file, "READ");
  const TVectorD *chain[2] = {
      (const TVectorD *)in.Get("band_spread_long_left"),
      (const TVectorD *)in.Get("band_spread_long_right")};
  const TParameter<Int_t> *bands = (TParameter<Int_t> *)in.Get("band_count");
  const TParameter<Double_t> *rate =
      (TParameter<Double_t> *)in.Get("long_rate_kHz");
  if (!chain[0] || !chain[1] || !bands || !rate || chain[0]->GetNrows() != 3 ||
      chain[1]->GetNrows() != 3) {
    std::cout << "skipping band_table: no band spreads in " << file
              << " (rerun scripts/run97_figure.sh)" << std::endl;
    return kFALSE;
  }
  for (Int_t c = 0; c < 2; c++)
    for (Int_t k = 0; k < 3; k++)
      data.spread[c][k] = (*chain[c])[k];
  data.bands = bands->GetVal();
  data.rate_kHz = rate->GetVal();
  if (data.bands != kNBands) {
    std::cout << "skipping band_table: " << file << " has " << data.bands
              << " bands, the simulation " << kNBands << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

// The band table as a complete booktabs table float with its caption and
// label (tab:bands), for the paper to \input (needs siunitx). single holds
// the spread table's spreads, banded the same runs with kNBands bands; both
// in ADC, shown in percent of the simulated beam peak to sit beside run 97.
void BandTable(const ReadingSpreads &single, const ReadingSpreads &banded,
               Double_t beamPeak, const DataBands &data, const TString &path) {
  std::ofstream table(path.Data());
  table << "% Written by mains/figures.cpp; edit that, not this file."
        << std::endl;
  table << "\\begin{table}[htbp]" << std::endl;
  table << "\\centering" << std::endl;
  // The band count spelled out in the caption's prose.
  const char *const words[] = {"zero", "one", "two",   "three", "four",
                               "five", "six", "seven", "eight", "nine"};
  const TString count =
      kNBands < 10 ? TString(words[kNBands]) : TString::Format("%d", kNBands);
  table << Form("\\caption{Rms pulse-height spread, in percent of the beam "
                "peak, of pulses within 0.5--1.5 times the beam peak with no "
                "other pulse on the same strip within "
                "\\SI{5}{\\micro\\second}, "
                "before correction and after it with one kernel and with one "
                "kernel for each of %s bands of the previous pulse height, "
                "below 0.5, ",
                count.Data());
  for (Int_t b = 1; b < kNBands - 1; b++)
    table << Form("%g--%g, ", b - 0.5, b + 0.5);
  table << Form("and above %g times the beam peak. The simulated runs are "
                "those of Table~\\ref{tab:spread}. Run 97 pools all Long "
                "Strips, Right and all Long Strips, Left, each strip in units "
                "of its own beam peak, and is corrected with the kernels "
                "fitted on it.}",
                kNBands - 1.5)
        << std::endl;
  table << "\\label{tab:bands}" << std::endl;
  // Scoped by the table; without them the tabular overflows the text width of
  // the elsarticle preprint.
  table << "\\small" << std::endl;
  table << "\\setlength{\\tabcolsep}{4pt}" << std::endl;
  // The baseline reading heads each block, so the rows need only the rate.
  table << "\\begin{tabular}{@{}lcccccc@{}}" << std::endl;
  table << "\\toprule" << std::endl;
  table << "& \\multicolumn{3}{c}{Long Strips, Right} & "
           "\\multicolumn{3}{c}{Long Strips, Left} \\\\"
        << std::endl;
  table << "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}" << std::endl;
  table << Form("Rate & before & 1 band & %d bands & before & 1 band & %d "
                "bands \\\\",
                kNBands, kNBands)
        << std::endl;
  const Int_t sides[2] = {1, 0};
  const Double_t percent = 100.0 / beamPeak;
  for (Int_t j = 0; j < kNSpreadReadings; j++) {
    const RateSpreads &one = single.find(kSpreadReadings[j].key)->second;
    const RateSpreads &many = banded.find(kSpreadReadings[j].key)->second;
    table << "\\midrule" << std::endl;
    table << "\\multicolumn{7}{l}{\\textit{Simulation, "
          << kSpreadReadings[j].label << "}} \\\\" << std::endl;
    for (Int_t i = 0; i < kNSpreadRates; i++) {
      const Spread &a = one.find(kSpreadRates[i])->second;
      const Spread &b = many.find(kSpreadRates[i])->second;
      table << Form("\\SI{%.0f}{\\kilo\\hertz}", kSpreadRates[i]);
      for (Int_t s = 0; s < 2; s++)
        table << Form(" & %.2f & %.2f & %.2f", percent * a.raw[sides[s]],
                      percent * a.corr[sides[s]], percent * b.corr[sides[s]]);
      table << " \\\\" << std::endl;
    }
  }
  table << "\\midrule" << std::endl;
  table << "\\multicolumn{7}{l}{\\textit{Run 97}} \\\\" << std::endl;
  table << Form("\\SI{%.0f}{\\kilo\\hertz}", data.rate_kHz);
  for (Int_t s = 0; s < 2; s++)
    table << Form(" & %.2f & %.2f & %.2f", data.spread[sides[s]][0],
                  data.spread[sides[s]][1], data.spread[sides[s]][2]);
  table << " \\\\" << std::endl;
  table << "\\bottomrule" << std::endl;
  table << "\\end{tabular}" << std::endl;
  table << "\\end{table}" << std::endl;
  std::cout << "wrote " << path << std::endl;
}

// Run 97's strip beam-peak FWHM as LaTeX macros for the text, in percent of
// the peak, before and after the correction, for the strips S1-S16 whose long
// segment the right chain reads (pole-zero mismatched) and those whose long
// segment the left chain reads: the mean over the strips, and the narrowest
// and widest strip. The spectra figures show the peaks themselves.
void FwhmMacros(const TString &file, const TString &path) {
  TFile in(file, "READ");
  const TVectorD *fwhm[2] = {(const TVectorD *)in.Get("fwhm_raw"),
                             (const TVectorD *)in.Get("fwhm_corrected")};
  const TVectorD *longRight = (const TVectorD *)in.Get("strip_long_right");
  const Int_t firstStrip = 1, lastStrip = 16;
  if (!fwhm[0] || !fwhm[1] || !longRight || fwhm[0]->GetNrows() <= lastStrip ||
      fwhm[1]->GetNrows() <= lastStrip || longRight->GetNrows() <= lastStrip) {
    std::cout << "skipping run97_fwhm: no strip FWHM in " << file
              << " (rerun scripts/run97_figure.sh)" << std::endl;
    return;
  }
  // [0] long segment on the right chain, [1] on the left; [before, after].
  Double_t n[2] = {0.0, 0.0};
  Double_t mean[2][2] = {{0.0}}, lo[2][2] = {{1.0e9, 1.0e9}, {1.0e9, 1.0e9}},
           hi[2][2] = {{0.0}};
  for (Int_t s = firstStrip; s <= lastStrip; s++) {
    if ((*longRight)[s] < 0.0)
      continue;
    const Int_t g = (*longRight)[s] > 0.5 ? 0 : 1;
    n[g]++;
    for (Int_t k = 0; k < 2; k++) {
      const Double_t value = (*fwhm[k])[s];
      if (!(value > 0.0)) {
        std::cout << "skipping run97_fwhm: no beam peak found in strip " << s
                  << std::endl;
        return;
      }
      mean[g][k] += value;
      lo[g][k] = TMath::Min(lo[g][k], value);
      hi[g][k] = TMath::Max(hi[g][k], value);
    }
  }
  if (!(n[0] > 0.0) || !(n[1] > 0.0)) {
    std::cout << "skipping run97_fwhm: strips of only one kind in " << file
              << std::endl;
    return;
  }
  std::cout << "run 97 strip FWHM [% of peak], before -> after:";
  for (Int_t s = firstStrip; s <= lastStrip; s++)
    std::cout << Form(" S%d%s %.1f->%.1f", s,
                      (*longRight)[s] > 0.5 ? "(R)" : "(L)", (*fwhm[0])[s],
                      (*fwhm[1])[s]);
  std::cout << std::endl;
  std::ofstream macros(path.Data());
  macros << "% Written by mains/figures.cpp; edit that, not this file."
         << std::endl;
  macros << "% Run 97's strip beam-peak FWHM, in percent of the peak, before "
            "and after the"
         << std::endl;
  macros << "% correction, for the strips S1-S16 with a Long Strip, Right "
            "(LongRight) and"
         << std::endl;
  macros << "% those with a Long Strip, Left (LongLeft):" << std::endl;
  macros << "% the mean over the strips, then the narrowest and the widest "
            "strip."
         << std::endl;
  const char *side[2] = {"LongRight", "LongLeft"};
  const char *when[2] = {"Before", "After"};
  for (Int_t g = 0; g < 2; g++)
    for (Int_t k = 0; k < 2; k++) {
      macros << Form("\\newcommand{\\Fwhm%s%s}{%.1f}", side[g], when[k],
                     mean[g][k] / n[g])
             << std::endl;
      macros << Form("\\newcommand{\\Fwhm%s%sMin}{%.1f}", side[g], when[k],
                     lo[g][k])
             << std::endl;
      macros << Form("\\newcommand{\\Fwhm%s%sMax}{%.1f}", side[g], when[k],
                     hi[g][k])
             << std::endl;
    }
  std::cout << "wrote " << path << std::endl;
}

// A profile's bin means against the bin centre, leaving out bins with fewer
// than minEntries entries.
TGraphErrors *ProfileGraph(const TProfile *profile, Double_t minEntries) {
  TGraphErrors *graph = new TGraphErrors();
  for (Int_t bin = 1; bin <= profile->GetNbinsX(); bin++) {
    if (profile->GetBinEntries(bin) < TMath::Max(minEntries, 1.0))
      continue;
    const Int_t n = graph->GetN();
    graph->SetPoint(n, TMath::Power(10.0, profile->GetBinCenter(bin)),
                    profile->GetBinContent(bin));
    graph->SetPointError(n, 0.0, profile->GetBinError(bin));
  }
  return graph;
}

// Run 97's long-right profile before and after correction: the mean, over the
// pulses in each bin of the time since the previous pulse, of the pulse height
// less its channel's average over all beam-like pulses (so the profile
// averages to zero). Bins with fewer than kMinProfilePulses pulses are left
// out; at 44 kHz that ends the profile near 220 us, short of the kernel's
// 316 us reach, which the fit fills from every earlier pulse rather than the
// previous one alone.
const Double_t kMinProfilePulses = 200.0;

void DeviationFigure(const TString &file) {
  TFile in(file, "READ");
  const TProfile *beforeProfile = (const TProfile *)in.Get("before_long_right");
  const TProfile *afterProfile = (const TProfile *)in.Get("after_long_right");
  if (!beforeProfile || !afterProfile) {
    std::cerr << "no long-right profiles in " << file << std::endl;
    return;
  }
  TGraphErrors *before = ProfileGraph(beforeProfile, kMinProfilePulses);
  TGraphErrors *after = ProfileGraph(afterProfile, kMinProfilePulses);
  const Double_t xLo = 1.0, xHi = 350.0, yLo = -150.0, yHi = 300.0;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetLogx();
  // The mean in TLatex angle brackets, since "Mean" does not fit on one line.
  // They come from the Symbol font, which ROOT PDFs name but do not embed: a
  // viewer with no Symbol font of its own leaves them out.
  TH1F *frame = Frame(canvas, xLo, yLo, xHi, yHi, "Time Between Pulses [#mus]",
                      "#LTPulse Height #minus#kern[0.12]{Average}#GT [ADC]");
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

// Run 97's strip spectra before or after the correction, in units of each
// strip's beam peak: the figure MUSIC's events summary draws, for the
// segmented strips 1-16 only, with the side legend sized to its entries.
void StripSpectraFigure(const TString &file, Bool_t corrected) {
  const Int_t firstStrip = 1;
  const Int_t nStrips = 16;
  std::vector<TH1D *> spectra;
  std::vector<TString> labels;
  {
    TFile in(file, "READ");
    for (Int_t s = firstStrip; s < firstStrip + nStrips; s++) {
      const TH1D *stored = (const TH1D *)in.Get(
          Form("strip_%s_S%d", corrected ? "corrected" : "raw", s));
      if (!stored) {
        std::cerr << "no strip spectra in " << file
                  << " (rerun scripts/run97_figure.sh)" << std::endl;
        return;
      }
      TH1D *spectrum = (TH1D *)stored->Clone(PlottingUtils::GetRandomName());
      spectrum->SetDirectory(nullptr);
      spectra.push_back(spectrum);
      labels.push_back(Form("S%d", s));
    }
  }
  TPad *plot = nullptr;
  TPad *legendPad = nullptr;
  const Double_t legendWidthPx = PlottingUtils::LegendWidthPx(labels);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvasWithSideLegend(
      plot, legendPad, Int_t(legendWidthPx), 800, kTRUE);
  const std::vector<Int_t> colours = PlottingUtils::GetDefaultColors();
  Double_t yMax = 0.0;
  for (Int_t s = 0; s < nStrips; s++)
    yMax = TMath::Max(yMax, spectra[s]->GetMaximum());
  spectra[0]->SetMaximum(3.0 * yMax);
  spectra[0]->SetMinimum(0.5);
  PlottingUtils::ConfigureAndDrawHistogram(spectra[0], colours[0],
                                           ";#DeltaE [a.u.];Counts");
  for (Int_t s = 1; s < nStrips; s++) {
    PlottingUtils::ConfigureHistogram(spectra[s],
                                      colours[s % Int_t(colours.size())]);
    spectra[s]->Draw("HIST SAME");
  }
  plot->RedrawAxis();
  TLegend *legend = PlottingUtils::AddSideLegend(plot, legendPad, nStrips);
  // Two corrections to AddSideLegend's box. Sixteen of its 50 px rows are
  // taller than the frame, so the box would fill the canvas and cover the last
  // x-axis label: fit it to the frame's height (both pads span the canvas
  // height, so the frame's margins carry over). And it spans the whole legend
  // pad, from the frame's right edge to the canvas edge, which also holds the
  // plot's right margin: narrow it to the width its entries need, centred
  // between the frame and the edge. ROOT paints a drawn legend from its pad
  // coordinates, which in the axis-less legend pad equal NDC, so both are set.
  const Double_t frameBottom = plot->GetBottomMargin();
  const Double_t frameTop = 1.0 - plot->GetTopMargin();
  const Double_t legendPadPx = legendPad->GetWw() * legendPad->GetAbsWNDC();
  const Double_t boxWidth = legendWidthPx / legendPadPx;
  const Double_t boxLeft = 0.5 * (1.0 - boxWidth);
  const Double_t boxRight = boxLeft + boxWidth;
  legend->SetX1NDC(boxLeft);
  legend->SetX2NDC(boxRight);
  legend->SetY1NDC(frameBottom);
  legend->SetY2NDC(frameTop);
  legend->SetX1(boxLeft);
  legend->SetX2(boxRight);
  legend->SetY1(frameBottom);
  legend->SetY2(frameTop);
  for (Int_t s = 0; s < nStrips; s++)
    legend->AddEntry(spectra[s], labels[s], "l");
  PlottingUtils::AddText(corrected ? "After Correction" : "Before Correction",
                         0.88, 0.92);
  PanelLabel(corrected ? "(b)" : "(a)");
  PlottingUtils::SaveFigure(
      canvas, corrected ? "run97_spectra_corrected" : "run97_spectra", "",
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
  TString out = "doc/figs";
  TString spreadDir = "output/spread";
  TString reactionDir = "output/reactions";
  TString run97 = "output/run97/run97.root";
  TString long1 = "output/figures/long_1kHz.root";
  TString long40 = "output/figures/long_40kHz.root";
  PlotSaveFormat format = PlotSaveFormat::kPNGAndPDF;
  for (Int_t i = 1; i < argc; i++) {
    TString argument = argv[i];
    if (argument.BeginsWith("--out="))
      out = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--spread="))
      spreadDir = argument(9, argument.Length() - 9);
    else if (argument.BeginsWith("--reactions="))
      reactionDir = argument(12, argument.Length() - 12);
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

  if (HaveInput(run97, "kernel_run97")) {
    std::vector<std::vector<Double_t>> kernels = ReadKernels(run97);
    if (kernels.size() < 3)
      std::cerr << "no kernels in " << run97
                << " (rerun scripts/run97_figure.sh)" << std::endl;
    else
      KernelFigure(kernels);
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
        haveSpreads =
            MeasureSpread(spreadFiles[n], spreadLabels[n],
                          spreads[kSpreadReadings[j].key][allRates[i]], info);
        if (n == 0)
          runs = info;
        else if (info.tauL_us != runs.tauL_us || info.tauR_us != runs.tauR_us ||
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

    // The band table: two of the spread table's readings corrected again with
    // kNBands amplitude bands, beside what run 97 shows.
    DataBands dataBands;
    if (haveSpreads && HaveInput(run97, "band_table") &&
        ReadDataBands(run97, dataBands)) {
      ReadingSpreads banded;
      Bool_t haveBands = kTRUE;
      for (Int_t j = 0; j < kNSpreadReadings && haveBands; j++)
        for (Int_t i = 0; i < kNSpreadRates && haveBands; i++) {
          const char *key = kSpreadReadings[j].key;
          const TString label = Form("%s_r%g", key, kSpreadRates[i]);
          RunInfo info;
          haveBands =
              MeasureSpread(spreadDir + "/" + label + ".root",
                            label + Form(" %d bands", kNBands),
                            banded[key][kSpreadRates[i]], info, kNBands);
        }
      if (haveBands)
        BandTable(spreads, banded, runs.beamPeak, dataBands,
                  out + "/band_table.tex");
    }

    Bool_t haveReactions = kTRUE;
    std::vector<TString> reactionFiles;
    for (Int_t j = 0; j < kNReactionRuns && haveReactions; j++) {
      reactionFiles.push_back(reactionDir + Form("/%s_r%g.root",
                                                 kReactionRuns[j].key,
                                                 kReactionRuns[j].rate_kHz));
      haveReactions = HaveInput(reactionFiles.back(), "reaction_table");
    }
    std::vector<ReactionScore> scores(kNReactionRuns);
    RunInfo reactionRuns;
    Double_t fraction = 0.0;
    for (Int_t j = 0; j < kNReactionRuns && haveReactions; j++) {
      RunInfo info;
      Double_t runFraction = 0.0;
      haveReactions =
          ScoreReactions(reactionFiles[j], gSystem->BaseName(reactionFiles[j]),
                         scores[j], info, runFraction);
      if (j == 0) {
        reactionRuns = info;
        fraction = runFraction;
      }
    }
    if (haveReactions)
      ReactionTable(scores, reactionRuns, fraction,
                    out + "/reaction_table.tex");
  }

  if (HaveInput(run97, "run97 figures")) {
    DeviationFigure(run97);
    StripSpectraFigure(run97, kFALSE);
    StripSpectraFigure(run97, kTRUE);
    if (format != PlotSaveFormat::kPNG)
      FwhmMacros(run97, out + "/run97_fwhm.tex");
  }
  if (format == PlotSaveFormat::kPNGAndPDF)
    MovePngs(out);
  return 0;
}
