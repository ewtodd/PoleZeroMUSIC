#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "Preamp.hpp"
#include "TrapezoidFilter.hpp"
#include <TArrow.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPad.h>
#include <TProfile.h>
#include <TString.h>
#include <TSystem.h>
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

// "=== <reading>_r<rate>: ..." blocks and their
// "  L: ... RMS a -> corrected mean m RMS b ADC" lines.
void ReadSpreads(const TString &log, ReadingSpreads &out) {
  std::ifstream in(log.Data());
  std::string line;
  TString reading;
  Double_t rate = -1.0;
  while (std::getline(in, line)) {
    TString text(line.c_str());
    if (text.BeginsWith("=== ")) {
      TString label = text(4, text.Index(":") - 4);
      const Ssiz_t lastR = label.Last('r');
      reading = label(0, lastR - 1);
      rate = TString(label(lastR + 1, label.Length() - lastR - 1)).Atof();
      continue;
    }
    const Int_t side = text.BeginsWith("  L: ")   ? 0
                       : text.BeginsWith("  R: ") ? 1
                                                  : -1;
    if (side < 0 || rate < 0.0 || !text.Contains("corrected mean"))
      continue;
    const Ssiz_t first = text.Index(" RMS ");
    const Ssiz_t second = text.Index(" RMS ", first + 5);
    out[reading][rate].raw[side] = TString(text(first + 5, 10)).Atof();
    out[reading][rate].corr[side] = TString(text(second + 5, 10)).Atof();
  }
}

// ----------------------------------------------------------- shared drawing
const Int_t kRightColour = kRed + 1;
const Int_t kLeftColour = kBlue + 1;
const Int_t kShortLeftColour = kCyan + 2;
const Int_t kMatchedColour = kGreen + 2;
const Int_t kReadingColour[3] = {kOrange + 1, kMagenta + 1, kTeal + 3};

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
        "Filter Output  #hat{#it{s}}(t) [a.u.]");
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

  Note(endUs + 1.55, 0.31, "Over-compensated", 10, colour[1]);
  Note(endUs + 4, -0.14, "Under-compensated", 13, colour[0]);
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
        "Filter Output  #hat{#it{s}}(t)  [a.u.]");
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
  Frame(canvas, tMin, yLo, tMax, yHi,
        "time since pulse 1  #it{t}  [#mus]",
        "filter output  #hat{#it{s}}  [a.u.]");
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
  Note(0.0, 1.17, "pulse 1", 21);
  Note(dt, 1.17, "pulse 2", 21);
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

  const Double_t xText = 14.2;
  TLine *leadLoss =
      new TLine(dt + tm + 0.15, 0.5 * (yAlone + yPair), xText - 0.2, 0.70);
  leadLoss->SetLineWidth(PlottingUtils::GetLineWidth());
  leadLoss->Draw();
  Note(xText, 0.70, Form("loss = %.3f", yAlone - yPair), 12);
  TLine *leadTail =
      new TLine(dt + tm + 0.15, 0.5 * yTail, xText - 0.2, -0.45);
  leadTail->SetLineWidth(PlottingUtils::GetLineWidth());
  leadTail->Draw();
  Note(xText, -0.45, Form("tail of pulse 1 = %.3f", yTail), 12);

  std::vector<TString> labels;
  labels.push_back("both pulses");
  labels.push_back("pulse 1 alone");
  labels.push_back("pulse 2 alone");
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

// ----------------------------------------------------------- data figures
// (Fig. 3) The kernels PulseHistory fitted on run 97.
void KernelFigure(const std::vector<std::vector<Double_t>> &kernels) {
  const char *name[3] = {"long left", "long right", "short left"};
  const Int_t colour[3] = {kLeftColour, kRightColour, kShortLeftColour};
  const Int_t order[3] = {1, 0, 2};
  const Double_t xLo = 1.0, xHi = 350.0, yLo = -0.2, yHi = 0.75;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetLogx();
  TH1F *frame = Frame(canvas, xLo, yLo, xHi, yHi,
                      "time since the earlier pulse  #Delta#it{t}  [#mus]",
                      "shift per unit earlier pulse height");
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();

  TBox *inside = new TBox(xLo, yLo, 3.5, yHi);
  inside->SetFillColorAlpha(kGray, 0.35);
  inside->Draw();
  TBox *frozen = new TBox(3.5, yLo, 7.0, yHi);
  frozen->SetFillColorAlpha(kAzure - 9, 0.35);
  frozen->Draw();
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
  labels.push_back("read point in the earlier trapezoid");
  labels.push_back("baseline still frozen");
  TLegend *legend = Legend(labels, 0.52, 0.86);
  for (Int_t j = 0; j < 3; j++)
    legend->AddEntry(graphs[j], labels[j], "LP");
  legend->AddEntry(inside, labels[3], "F");
  legend->AddEntry(frozen, labels[4], "F");
  legend->Draw();
  PlottingUtils::SaveFigure(canvas, "kernel_run97", "",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
}

// (Fig. 4) Simulated spread against rate for the three baseline readings.
void SpreadFigure(const ReadingSpreads &spreads, Int_t side) {
  struct Reading {
    const char *key;
    const char *label;
  };
  const std::vector<Reading> readings = {
      {"ideal", "never takes in the tail"},
      {"wait", "frozen until within 80 ADC"},
      {"follow", "restarts after the freeze"}};
  const Double_t xHi = 85.0, yLo = 0.9;
  const Double_t yHi = side == 1 ? 3.0 : 2.3;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, 0.0, yLo, xHi, yHi, "rate per channel  [kHz]",
        "rms / rms of isolated pulses");
  Reference(0.0, xHi, 1.0);

  std::vector<TGraph *> corrected;
  std::vector<TString> labels;
  for (Int_t j = 0; j < Int_t(readings.size()); j++) {
    ReadingSpreads::const_iterator reading = spreads.find(readings[j].key);
    if (reading == spreads.end())
      continue;
    RateSpreads::const_iterator floorPoint = reading->second.find(0.5);
    if (floorPoint == reading->second.end())
      continue;
    const Double_t floor = floorPoint->second.corr[side];
    TGraph *before = new TGraph();
    TGraph *after = new TGraph();
    for (RateSpreads::const_iterator point = reading->second.begin();
         point != reading->second.end(); ++point) {
      if (point->first < 1.0)
        continue;
      before->SetPoint(before->GetN(), point->first,
                       point->second.raw[side] / floor);
      after->SetPoint(after->GetN(), point->first,
                      point->second.corr[side] / floor);
      std::cout << Form("  %s %-7s %5.1f kHz: %6.1f -> %6.1f ADC (floor %.1f)",
                        side ? "R" : "L", readings[j].key, point->first,
                        point->second.raw[side], point->second.corr[side],
                        floor)
                << std::endl;
    }
    PlottingUtils::ConfigureGraph(before, kReadingColour[j]);
    before->SetLineStyle(2);
    before->SetMarkerStyle(24 + j);
    before->SetMarkerColor(kReadingColour[j]);
    before->SetMarkerSize(1.2);
    before->Draw("LP SAME");
    PlottingUtils::ConfigureGraph(after, kReadingColour[j]);
    after->SetMarkerStyle(20 + j);
    after->SetMarkerColor(kReadingColour[j]);
    after->SetMarkerSize(1.2);
    after->Draw("LP SAME");
    corrected.push_back(after);
    labels.push_back(readings[j].label);
  }
  canvas->RedrawAxis();

  PlottingUtils::AddText(side ? "long right segments" : "long left segments",
                         0.87, 0.865);
  if (side == 1) {
    TLegend *legend = Legend(labels, 0.645, 0.86, "baseline", kTRUE);
    for (Int_t j = 0; j < Int_t(corrected.size()); j++)
      legend->AddEntry(corrected[j], labels[j], "LP");
    legend->Draw();
  } else {
    std::vector<TString> keys;
    keys.push_back("before correction");
    keys.push_back("after correction");
    TLegend *legend = Legend(keys, 0.70, 0.86, "", kTRUE);
    TGraph *beforeKey = new TGraph();
    PlottingUtils::ConfigureGraph(beforeKey, kGray + 2);
    beforeKey->SetLineStyle(2);
    beforeKey->SetMarkerStyle(24);
    beforeKey->SetMarkerColor(kGray + 2);
    beforeKey->SetMarkerSize(1.2);
    TGraph *afterKey = new TGraph();
    PlottingUtils::ConfigureGraph(afterKey, kGray + 2);
    afterKey->SetMarkerStyle(20);
    afterKey->SetMarkerColor(kGray + 2);
    afterKey->SetMarkerSize(1.2);
    legend->AddEntry(beforeKey, keys[0], "LP");
    legend->AddEntry(afterKey, keys[1], "LP");
    legend->Draw();
  }
  PanelLabel(side ? "(a)" : "(b)");
  PlottingUtils::SaveFigure(canvas, side ? "sim_R" : "sim_L", "",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
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
                      "time since the previous pulse  #Delta#it{t}  [#mus]",
                      "pulse height #minus#kern[0.3]{mean}  [ADC]");
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
  labels.push_back("before correction");
  labels.push_back("after correction");
  TLegend *legend = Legend(labels, 0.70, 0.86, "", kTRUE);
  legend->AddEntry(before, labels[0], "LP");
  legend->AddEntry(after, labels[1], "LP");
  legend->Draw();
  PlottingUtils::AddText("long right segments, run 97", 0.87, 0.865);
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
  TString sims = "output/paper_sim/analysis.log";
  TString run97 = "output/run97/deviation.root";
  PlotSaveFormat format = PlotSaveFormat::kPNGAndPDF;
  for (Int_t i = 1; i < argc; i++) {
    TString argument = argv[i];
    if (argument.BeginsWith("--out="))
      out = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--log="))
      log = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--run="))
      run = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--sim="))
      sims = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--run97="))
      run97 = argument(8, argument.Length() - 8);
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

  if (HaveInput(log, "kernel_run97")) {
    std::vector<std::vector<Double_t>> kernels = ReadKernels(log, run);
    if (kernels.size() < 3) {
      std::cerr << "found " << kernels.size() << " kernels for " << run
                << " in " << log << std::endl;
    } else {
      KernelFigure(kernels);
    }
  }

  ReadingSpreads spreads;
  Bool_t haveSims = kTRUE;
  TObjArray *files = sims.Tokenize(",");
  for (Int_t i = 0; i < files->GetEntries(); i++) {
    const TString file = ((TObjString *)files->At(i))->GetString();
    if (HaveInput(file, "sim_R and sim_L"))
      ReadSpreads(file, spreads);
    else
      haveSims = kFALSE;
  }
  delete files;
  if (haveSims) {
    SpreadFigure(spreads, 1);
    SpreadFigure(spreads, 0);
  }

  if (HaveInput(run97, "run97_deviation_R"))
    DeviationFigure(run97);
  if (format == PlotSaveFormat::kPNGAndPDF)
    MovePngs(out);
  return 0;
}
