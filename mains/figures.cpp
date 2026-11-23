#include "Config.hpp"
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
#include <TDecompChol.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMatrixD.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPad.h>
#include <TParameter.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorD.h>
#include <algorithm>
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

// The rms of a chain's heights, and their mean, before (raw) and after (corr)
// correction, in ADC; [0] the left chain, [1] the right.
struct Spread {
  Double_t raw[2];
  Double_t corr[2];
  Double_t rawMean[2];
  Double_t corrMean[2];
  Spread() {
    for (Int_t ch = 0; ch < 2; ch++) {
      raw[ch] = 0.0;
      corr[ch] = 0.0;
      rawMean[ch] = 0.0;
      corrMean[ch] = 0.0;
    }
  }
  // In percent of the chain's own mean height, which carries its gain.
  Double_t RawPercent(Int_t ch) const {
    return rawMean[ch] > 0.0 ? 100.0 * raw[ch] / rawMean[ch] : 0.0;
  }
  Double_t CorrPercent(Int_t ch) const {
    return corrMean[ch] > 0.0 ? 100.0 * corr[ch] / corrMean[ch] : 0.0;
  }
};

typedef std::map<Double_t, Spread> RateSpreads;
typedef std::map<TString, RateSpreads> ReadingSpreads;

// The simulated runs (scripts/sim_inputs.sh), each at the two rates, with the
// left chain matched to the pole-zero and the right mismatched. The baseline
// restorer is CAEN's description of EnergyFilterLFLimitation = On, as on every
// channel of the 37Cl runs (data_raw/music_exp1915_097XSetting_66222.dat): at
// the end of each freeze the baseline resumes if the trapezoid is back within
// kRestartThrADC of it, and otherwise freezes again. CAEN does not say how
// close, nor what decides the check, so two readings of it; the first is the
// one the figures show.
//   average - the running average of the trapezoid over the last
//             kCheckSamples, kept up during the freeze
//   single  - the last trapezoid sample alone
struct SpreadReading {
  const char *key;
  const char *label;
};
const Int_t kNSpreadReadings = 2;
const SpreadReading kSpreadReadings[kNSpreadReadings] = {
    {"average", "running-average check"}, {"single", "single-sample check"}};
// As scripts/sim_inputs.sh simulates them (--baselineRestartThrADC=20,
// --baselineCheckSamples=1024 for the average).
const Double_t kRestartThrADC = 20.0;
const Int_t kCheckSamples = 1024;
const Int_t kNSpreadRates = 2;
const Double_t kSpreadRates[kNSpreadRates] = {10.0, 50.0};
// The run of each reading's file name, e.g. output/sim/average_r50.root.
TString SimRun(const TString &dir, const char *key, Double_t rate) {
  return dir + Form("/%s_r%g.root", key, rate);
}
// The beam window of the spreads, in units of each chain's beam peak: the
// window PulseHistory selects its beam pulses in (PULSE_HISTORY_OWN_LO/HI), wide
// enough for the pole-zero distortions, below the merged double pulses at
// twice the peak.
const Double_t kSpreadWindowLo = 0.5;
const Double_t kSpreadWindowHi = 1.6;

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

// Text sizes in pixels on the standard 1200x800 canvas: the notes and
// legends of the figures outside the paper (30 px), and the three levels of
// the paper's figures: legends (34 px), in-figure annotations (38 px; in the
// preprint they print at the body text's 12 pt), and axis titles and labels
// (the style's 0.06 of the canvas's 800 px side).
const Double_t kNotePx = 30.0;
const Double_t kLegendPx = 34.0;
const Double_t kAnnotationPx = 38.0;
const Double_t kAxisTextPx = 48.0;

// The text size of the paper's figures' legends: 34 px, or another with
// --legend-px=, to compare. Each figure's room for its legend follows the
// size.
Double_t gLegendPx = kLegendPx;
// The text size of the paper's in-figure annotations, corner notes and panel
// letters: the style's 38 px, or another with --annotation-px=, to compare.
Double_t gAnnotationPx = kAnnotationPx;
Bool_t BigLegend() { return gLegendPx > 0.5 * (kNotePx + kAxisTextPx); }
// A legend's height in NDC on the standard canvas, for rows of gLegendPx text
// with 10 px between them, but never under the style's 50 px rows.
Double_t LegendHeight(Int_t rows) {
  return rows * TMath::Max(gLegendPx + 10.0, 50.0) / 800.0;
}


// A note in the top margin, at the annotations' size: right-aligned with
// the frame's right edge and standing just above its top, so the larger text
// clears the frame.
TLatex *CornerText(const char *text) {
  TLatex *label = PlottingUtils::AddText(text, 1.0 - gPad->GetRightMargin(),
                                         1.0 - gPad->GetTopMargin() + 0.012);
  label->SetTextFont(43);
  label->SetTextSize(gAnnotationPx);
  label->SetTextAlign(31);
  return label;
}

// A legend inside the frame, sized to its labels and anchored a small inset
// from the frame's right (or left) edge, wherever the pad's margins put it;
// with several columns, each column is as wide as the widest label needs.
TLegend *Legend(const std::vector<TString> &labels, Double_t y1, Double_t y2,
                const TString &header = "", Bool_t left = kFALSE,
                Double_t sizePx = kNotePx, Int_t columns = 1) {
  const Double_t width = columns *
                         PlottingUtils::LegendWidthPx(labels, header, sizePx) /
                         gPad->GetWw();
  const Double_t inset = 0.02;
  const Double_t x2 = left ? gPad->GetLeftMargin() + inset + width
                           : 1.0 - gPad->GetRightMargin() - inset;
  TLegend *legend = PlottingUtils::AddLegend(x2 - width, x2, y1, y2);
  legend->SetFillStyle(1001);
  legend->SetTextSize(sizePx);
  if (columns > 1) {
    legend->SetNColumns(columns);
    // A column's symbol takes this share of the column, as in one column.
    legend->SetMargin(0.25);
  }
  if (!header.IsNull())
    legend->SetHeader(header);
  return legend;
}

void Reference(Double_t x1, Double_t x2, Double_t y,
               Width_t width = PlottingUtils::GetLineWidth()) {
  TLine *line = new TLine(x1, y, x2, y);
  line->SetLineStyle(3);
  line->SetLineColor(kGray + 2);
  line->SetLineWidth(width);
  line->Draw();
}

void TimeMark(Double_t x, Double_t yLo, Double_t yHi,
              Width_t width = PlottingUtils::GetLineWidth()) {
  TLine *line = new TLine(x, yLo, x, yHi);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->SetLineWidth(width);
  line->Draw();
}

void TimeArrow(Double_t x1, Double_t x2, Double_t y, Double_t textY,
               const char *text, Double_t sizePx = kNotePx,
               Width_t width = PlottingUtils::GetLineWidth()) {
  TArrow *arrow = new TArrow(x1, y, x2, y, 0.012, "<|>");
  arrow->SetLineWidth(width);
  arrow->SetFillColor(kBlack);
  arrow->Draw();
  TLatex *label = new TLatex(0.5 * (x1 + x2), textY, text);
  label->SetTextFont(43);
  label->SetTextSize(sizePx);
  label->SetTextAlign(21);
  label->Draw();
}

void Note(Double_t x, Double_t y, const char *text, Int_t align = 12,
          Int_t colour = kBlack, Double_t sizePx = kNotePx) {
  TLatex *label = new TLatex(x, y, text);
  label->SetTextFont(43);
  label->SetTextSize(sizePx);
  label->SetTextAlign(align);
  label->SetTextColor(colour);
  label->Draw();
}

void PanelLabel(const char *text, Double_t sizePx = gAnnotationPx) {
  TLatex *label = PlottingUtils::AddText(text, 0.02, 0.97);
  label->SetTextFont(63);
  label->SetTextSize(sizePx);
  label->SetTextAlign(13);
}

// --------------------------------------------- simulated filter outputs
// gamma of the paper's Eq. (eq:gamma): after (k+l)T a unit pulse leaves
// gamma exp(-t/tau_pre), an undershoot (gamma < 0) when tau_pre < tau_pz.
Double_t TailAmplitude(Double_t tauPre, Double_t tauPz, Double_t riseUs,
                       Double_t flatTopEndUs) {
  return (tauPre - tauPz) * tauPre / (riseUs * tauPz) *
         (1.0 - TMath::Exp(riseUs / tauPre)) *
         (1.0 - TMath::Exp(flatTopEndUs / tauPre));
}

TString TailLabel(Double_t gamma) {
  if (TMath::Abs(gamma) < 0.005)
    return "#gamma = 0";
  if (gamma < 0.0)
    return Form("#gamma = #minus#kern[0.15]{%.2f}", -gamma);
  return Form("#gamma = +%.2f", gamma);
}

// Filter output for pulses at arrivalsUs, of the given amplitudes (unit pulses
// when there are none), against time from the first arrival. A one-sample
// detector ramp keeps the trapezoid corners exact.
TGraph *StepResponse(
    Double_t tauPreUs, Double_t tauPzUs,
    const std::vector<Double_t> &arrivalsUs, Double_t tMinUs, Double_t tMaxUs,
    Int_t stride = 1,
    const std::vector<Double_t> &amplitudes = std::vector<Double_t>()) {
  const Long64_t pad = Long64_t(2000.0 / kSamplePeriodNs);
  const Long64_t nSamples =
      pad + Long64_t(tMaxUs * 1000.0 / kSamplePeriodNs) + 1;
  Preamp preamp(tauPreUs * 1000.0, kSamplePeriodNs, kSamplePeriodNs);
  TrapezoidFilter filter(Samples(kTrapRiseNs), Samples(kTrapFlatNs),
                         kSamplePeriodNs, tauPzUs * 1000.0);
  for (Int_t a = 0; a < Int_t(arrivalsUs.size()); a++)
    preamp.AddEvent(Double_t(pad) * kSamplePeriodNs + arrivalsUs[a] * 1000.0,
                    amplitudes.empty() ? 1.0 : amplitudes[a]);
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
// (k + l)T, with the energy read at t_m. Black and grey only; at the
// annotations' size the curves are labelled where they run instead of in a
// legend.
void TrapezoidFigure() {
  const Double_t tauPre = 10.0;
  const Double_t kUs = kTrapRiseNs / 1000.0;
  const Double_t lUs = (kTrapRiseNs + kTrapFlatNs) / 1000.0;
  const Double_t endUs = kUs + lUs;
  const Double_t tm = kUs + kPeakingDelayNs / 1000.0;
  const Double_t tMin = -1.0, tMax = 12.0;
  const Double_t yLo = -0.1, yHi = 2.5;
  const Double_t markTop = 1.08;
  const Width_t lineWidth = PlottingUtils::GetLineWidth();
  const Double_t px = gAnnotationPx;

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
  in->SetLineWidth(lineWidth);
  in->Draw("L SAME");
  PlottingUtils::ConfigureGraph(out, kBlack);
  out->SetLineWidth(lineWidth);
  out->Draw("L SAME");
  canvas->RedrawAxis();

  Note(0.5 * kUs, 1.12, "rise", 21, kBlack, px);
  Note(0.5 * (kUs + lUs), 1.12, "flat top", 21, kBlack, px);
  Note(0.5 * (lUs + endUs), 1.12, "fall", 21, kBlack, px);
  TimeArrow(0.0, kUs, 1.50, 1.56, "#it{kT}", px, lineWidth);
  TimeArrow(lUs, endUs, 1.50, 1.56, "#it{kT}", px, lineWidth);
  TimeArrow(0.0, lUs, 1.86, 1.92, "#it{lT}", px, lineWidth);
  TimeArrow(0.0, endUs, 2.22, 2.28, "(#it{k}+#it{l})#it{T}", px, lineWidth);

  TGraph *read = new TGraph(1);
  read->SetPoint(0, tm, out->Eval(tm));
  read->SetMarkerStyle(20);
  read->SetMarkerSize(2.0);
  read->SetMarkerColor(kBlack);
  read->Draw("P SAME");
  Note(tm + 0.12, 0.12, "#it{t}_{m}", 12, kBlack, px);

  Note(endUs + 0.6, 0.66, "preamplifier", 12, kBlack, px);
  Note(endUs + 0.3, 0.16, "filter output", 12, kBlack, px);
  Note(tMax - 0.5, 1.22, "#tau_{pz} = #tau_{pre}", 32, kBlack, px);
  PanelLabel("(a)", px);

  PlottingUtils::SaveFigure(canvas, "trapezoid", "", PlotSaveOptions::kLINEAR);
  delete canvas;
}

// (Fig. 1b) The same preamp through under- and over-compensated pole-zero
// settings; the tails after (k + l)T are the point. Panel (a) shows the
// matched case.
void PoleZeroFigure() {
  const Double_t tauPre = 10.0;
  const Double_t kUs = kTrapRiseNs / 1000.0;
  const Double_t lUs = (kTrapRiseNs + kTrapFlatNs) / 1000.0;
  const Double_t endUs = kUs + lUs;
  const Int_t nCases = 2;
  const Double_t tauPz[nCases] = {20.0, 5.0};
  const Int_t colour[nCases] = {kRightColour, kLeftColour};
  const Double_t tMin = -1.0, tMax = 45.0;
  const Double_t yLo = BigLegend() ? -0.62 : -0.45, yHi = 0.45;
  const Width_t lineWidth = PlottingUtils::GetLineWidth();
  const Double_t px = gAnnotationPx;
  std::vector<Double_t> one(1, 0.0);

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, tMin, yLo, tMax, yHi, "Time [#mus]", "Amplitude [a.u.]");
  Reference(tMin, tMax, 0.0);
  TimeMark(endUs, yLo, yHi);

  TGraph *graphs[nCases];
  std::vector<TString> labels;
  for (Int_t v = 0; v < nCases; v++) {
    graphs[v] = StepResponse(tauPre, tauPz[v], one, tMin, tMax, 6);
    PlottingUtils::ConfigureGraph(graphs[v], colour[v]);
    graphs[v]->SetLineWidth(lineWidth);
    graphs[v]->Draw("L SAME");
    labels.push_back(
        Form("#tau_{pz} = %.0f #mus,  %s", tauPz[v],
             TailLabel(TailAmplitude(tauPre, tauPz[v], kUs, lUs)).Data()));
  }
  canvas->RedrawAxis();

  Note(endUs + 2.0, 0.37, "over-compensated", 12, colour[1], px);
  Note(endUs + 6.0, -0.15, "under-compensated", 13, colour[0], px);
  Note(endUs - 0.5, yLo + 0.04, "(#it{k}+#it{l})#it{T}", 31, kBlack, px);

  TLegend *legend =
      Legend(labels, 0.18, 0.18 + LegendHeight(2), "", kFALSE, gLegendPx);
  for (Int_t v = 0; v < nCases; v++)
    legend->AddEntry(graphs[v], labels[v], "l");
  legend->Draw();
  TLatex *tauNote = PlottingUtils::AddText("#tau_{pre} = 10 #mus", 0.87, 0.84);
  tauNote->SetTextFont(43);
  tauNote->SetTextSize(px);
  PanelLabel("(b)", px);

  PlottingUtils::SaveFigure(canvas, "polezero", "", PlotSaveOptions::kLINEAR);
  delete canvas;
}

// (Fig. correction a) One pulse on the right chain (tau_pre = 10 us) through
// the 37Cl settings (tau_pz = 20 us), in black, as Fig. 1(a); the caption
// gives the two decay times.
void SinglePulseFigure() {
  const Double_t tauPre = 10.0;
  const Double_t kUs = kTrapRiseNs / 1000.0;
  const Double_t lUs = (kTrapRiseNs + kTrapFlatNs) / 1000.0;
  const Double_t endUs = kUs + lUs;
  const Double_t tm = kUs + kPeakingDelayNs / 1000.0;
  const Double_t tf = kBaselineHoldoffNs / 1000.0;
  const Double_t tMin = -1.0, tMax = 14.0;
  const Double_t yLo = -0.6, yHi = 2.4;
  const Double_t markTop = 1.08;
  const Double_t tail = TailAmplitude(tauPre, kPoleZeroUs, kUs, lUs);
  const Double_t px = gAnnotationPx;
  const Width_t lineWidth = PlottingUtils::GetLineWidth();

  std::vector<Double_t> one(1, 0.0);
  TGraph *out = StepResponse(tauPre, kPoleZeroUs, one, tMin, tMax, 8);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, tMin, yLo, tMax, yHi, "Time [#mus]",
        "Filter Output #hat{#it{s}} [a.u.]");
  Reference(tMin, tMax, 0.0);
  TimeMark(0.0, yLo, markTop);
  TimeMark(kUs, yLo, markTop);
  TimeMark(lUs, yLo, markTop);
  TimeMark(endUs, yLo, markTop);
  PlottingUtils::ConfigureGraph(out, kBlack);
  out->SetLineWidth(lineWidth);
  out->Draw("L SAME");
  canvas->RedrawAxis();

  Note(0.5 * kUs, 1.12, "rise", 21, kBlack, px);
  Note(0.5 * (kUs + lUs), 1.12, "flat top", 21, kBlack, px);
  Note(0.5 * (lUs + endUs), 1.12, "fall", 21, kBlack, px);
  TimeArrow(0.0, tm, 1.50, 1.58, "#it{t}_{m}", px, lineWidth);
  // The freeze's label is wider than its arrow: left-aligned above it.
  TimeArrow(0.0, tf, 1.90, 1.98, "", px, lineWidth);
  Note(0.0, 1.98, "baseline frozen for #it{t}_{f}", 11, kBlack, px);
  TimeMark(tm, yLo, out->Eval(tm));
  TGraph *read = new TGraph(1);
  read->SetPoint(0, tm, out->Eval(tm));
  read->SetMarkerStyle(20);
  read->SetMarkerSize(2.0);
  read->SetMarkerColor(kBlack);
  read->Draw("P SAME");
  // A unit pulse: the height read is the gain g of the paper's Eq. (gain).
  // The leader leaves the read point below the sloping flat top.
  TLine *leader =
      new TLine(tm + 0.2, out->Eval(tm) - 0.04, endUs + 0.2, 0.62);
  leader->SetLineWidth(lineWidth);
  leader->Draw();
  Note(endUs + 0.4, 0.62,
       Form("#it{g} = #hat{#it{s}}_{1}(#it{t}_{m}) = %.2f", out->Eval(tm)), 12,
       kBlack, px);
  Note(endUs + 0.4, -0.38, Form("tail,  %s", TailLabel(tail).Data()), 12,
       kBlack, px);
  PanelLabel("(a)");

  PlottingUtils::SaveFigure(canvas, "single", "", PlotSaveOptions::kLINEAR);
  delete canvas;
}

// (Fig. correction b) What the correction models, on a short train on the
// right chain (tau_pre = 10 us) with an ideal baseline: the height read on the
// last pulse is its height alone plus the tail each earlier pulse leaves at
// that read point, kappa(dt_i) h~_i with h~_i the earlier pulse's height alone
// (the paper's Eq. (kernel)). Kept in its first design and colours, with the
// paper's text sizes: annotations at gAnnotationPx, the legend at gLegendPx.
void CorrectionTrainFigure() {
  const Double_t tauPre = 10.0;
  const Double_t tm = (kTrapRiseNs + kPeakingDelayNs) / 1000.0;
  const Int_t nPulses = 4, last = nPulses - 1;
  const Double_t arrival[nPulses] = {0.0, 8.0, 16.0, 24.0};
  const Double_t amplitude[nPulses] = {1.0, 1.5, 0.8, 1.0};
  const Double_t tRead = arrival[last] + tm;
  const Double_t tMin = -1.0, tMax = 44.0;
  const Double_t yLo = -0.7, yHi = 2.75;
  const Double_t px = gAnnotationPx;

  const std::vector<Double_t> times(arrival, arrival + nPulses);
  const std::vector<Double_t> heights(amplitude, amplitude + nPulses);
  TGraph *train =
      StepResponse(tauPre, kPoleZeroUs, times, tMin, tMax, 8, heights);
  // Each pulse alone, from its arrival: before it, it is the zero line.
  std::vector<TGraph *> alone;
  for (Int_t i = 0; i < nPulses; i++)
    alone.push_back(StepResponse(
        tauPre, kPoleZeroUs, std::vector<Double_t>(1, arrival[i]), arrival[i],
        tMax, 8, std::vector<Double_t>(1, amplitude[i])));

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, tMin, yLo, tMax, yHi, "Time [#mus]",
        "Filter Output #hat{#it{s}} [a.u.]");
  Reference(tMin, tMax, 0.0);
  for (Int_t i = 0; i < last; i++) {
    PlottingUtils::ConfigureGraph(alone[i], kGray + 2);
    alone[i]->SetLineStyle(2);
    alone[i]->Draw("L SAME");
  }
  PlottingUtils::ConfigureGraph(alone[last], kLeftColour);
  alone[last]->SetLineStyle(2);
  alone[last]->Draw("L SAME");
  PlottingUtils::ConfigureGraph(train, kRightColour);
  train->Draw("L SAME");
  canvas->RedrawAxis();

  // Each earlier pulse's corrected height (its height alone), under its flat
  // top, and its time before the last pulse.
  TGraph *earlierReads = new TGraph();
  for (Int_t i = 0; i < last; i++) {
    const Double_t h = alone[i]->Eval(arrival[i] + tm);
    earlierReads->SetPoint(i, arrival[i] + tm, h);
    Note(arrival[i] + tm, h + 0.10, Form("#tilde{#it{h}}_{%d}", i + 1), 21,
         kBlack, px);
    const Double_t y = 1.70 + 0.32 * (last - 1 - i);
    TimeArrow(arrival[i], arrival[last], y, y + 0.07,
              Form("#Delta#it{t}_{%d}", i + 1), px);
  }
  // The read time, from just below the tails, clear of their note.
  TimeMark(tRead, -0.35, alone[last]->Eval(tRead));

  // At the last read point: the height alone, the height read, and the
  // tails that make up the difference.
  const Double_t yAlone = alone[last]->Eval(tRead);
  const Double_t yRead = train->Eval(tRead);
  Double_t tailSum = 0.0;
  TGraph *tails = new TGraph();
  for (Int_t i = 0; i < last; i++) {
    const Double_t tail = alone[i]->Eval(tRead);
    tails->SetPoint(i, tRead, tail);
    tailSum += tail;
  }
  TLine *shift = new TLine(tRead, yRead, tRead, yAlone);
  shift->SetLineWidth(3 * PlottingUtils::GetLineWidth());
  shift->Draw();
  TGraph *reads = new TGraph(2);
  reads->SetPoint(0, tRead, yAlone);
  reads->SetPoint(1, tRead, yRead);
  reads->SetMarkerStyle(20);
  reads->SetMarkerSize(2.0);
  reads->SetMarkerColor(kBlack);
  reads->Draw("P SAME");
  tails->SetMarkerStyle(20);
  tails->SetMarkerSize(1.2);
  tails->SetMarkerColor(kGray + 2);
  tails->Draw("P SAME");

  earlierReads->SetMarkerStyle(20);
  earlierReads->SetMarkerSize(1.2);
  earlierReads->SetMarkerColor(kBlack);
  earlierReads->Draw("P SAME");

  const Double_t xText = tRead + 2.8;
  TLine *leadShift = new TLine(tRead + 0.2, 0.5 * (yAlone + yRead),
                               xText - 0.2, 0.5 * (yAlone + yRead));
  leadShift->SetLineWidth(PlottingUtils::GetLineWidth());
  leadShift->Draw();
  Note(xText, 0.5 * (yAlone + yRead),
       "#sum_{#it{i}} #kappa(#Delta#it{t}_{#it{i}}) #tilde{#it{h}}_{#it{i}}",
       12, kBlack, px);
  // The tails' note right-aligned below the train, its leader from the lowest
  // tail at the read point.
  const Double_t xTailText = tMax - 0.5, yTailText = -0.5;
  TLine *leadTails = new TLine(tRead + 0.2, tails->GetY()[last - 1],
                               tRead + 1.2, yTailText + 0.10);
  leadTails->SetLineWidth(PlottingUtils::GetLineWidth());
  leadTails->Draw();
  Note(xTailText, yTailText, "tails of pulses 1-3", 32, kBlack, px);

  std::vector<TString> labels;
  labels.push_back("Measured");
  labels.push_back("Pulses 1-3 alone");
  labels.push_back("Pulse 4 alone");
  TLegend *legend =
      Legend(labels, 0.87 - LegendHeight(3), 0.87, "", kFALSE, gLegendPx);
  legend->AddEntry(train, labels[0], "l");
  legend->AddEntry(alone[0], labels[1], "l");
  legend->AddEntry(alone[last], labels[2], "l");
  legend->Draw();
  PanelLabel("(b)");

  PlottingUtils::SaveFigure(canvas, "correction_train", "",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
  std::cout << Form("correction train: pulse 4 alone %.3f, read %.3f, "
                    "difference %.3f, sum of tails %.3f",
                    yAlone, yRead, yRead - yAlone, tailSum)
            << std::endl;
}

// What the figures take from a simulated run's parameters.
struct RunInfo {
  Double_t beamPeak = 1600.0;
  Double_t tauL_us = 0.0;
  Double_t tauR_us = 0.0;
  Double_t tauPz_us = 0.0;
  Double_t sigmaFrac = 0.0;
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
  info.sigmaFrac = RunParameter(in, "beamPeakSigmaFrac", 0.0);

  Int_t channel = 0;
  ULong64_t timestamp = 0;
  Double_t energy = 0.0, dtPrev = 0.0, dtNext = 0.0, trueAmp = 0.0,
           prevTrueAmp = 0.0;
  Bool_t pileup = kFALSE;
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
    hits.push_back(hit);
  }
  in.Close();
  return kTRUE;
}

// The pulse-height spectrum of one chain of a simulated run (channel 0 the
// left, 1 the right), from its non-pileup hits. The corrected heights are
// whole ADC counts, as PulseHistory stores them, so the bins are a whole number
// of counts wide with edges on multiples of that width: every bin then holds
// the same number of possible heights, where 13.3-count bins alternately held
// 13 and 14 and combed the corrected spectra.
TH1D *ChainSpectrum(const std::vector<SimHit> &hits,
                    const std::vector<Double_t> &energies, Int_t channel,
                    Double_t beamPeak) {
  const Double_t binWidth =
      TMath::Max(1.0, Double_t(TMath::Nint(2.5 * beamPeak / 300.0)));
  const Int_t nBins = Int_t(TMath::Ceil(2.5 * beamPeak / binWidth));
  TH1D *spectrum =
      new TH1D(PlottingUtils::GetRandomName(), ";Pulse Height [ADC];Counts",
               nBins, 0.0, nBins * binWidth);
  spectrum->SetDirectory(nullptr);
  for (size_t i = 0; i < hits.size(); i++)
    if (!hits[i].pileup && hits[i].channel == channel)
      spectrum->Fill(energies[i]);
  return spectrum;
}

// One panel of the spectra figure: a run's matched (left chain) and
// mismatched (right chain) spectra, before the correction dashed and after it
// solid, in one frame.
void DrawSpectra(const std::vector<SimHit> &hits,
                 const std::vector<Double_t> &raw,
                 const std::vector<Double_t> &corrected, Double_t beamPeak,
                 const TString &rateLabel, const TString &name,
                 const char *panel) {
  // [chain][0 before, 1 after]
  TH1D *spectra[2][2];
  Double_t yMax = 0.0;
  for (Int_t chain = 0; chain < 2; chain++)
    for (Int_t k = 0; k < 2; k++) {
      spectra[chain][k] =
          ChainSpectrum(hits, k == 0 ? raw : corrected, chain, beamPeak);
      yMax = TMath::Max(yMax, spectra[chain][k]->GetMaximum());
    }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  // Five-digit counts need more room than the style's left margin. Drawing
  // applies the style, so the title offset has to be set afterwards.
  canvas->SetLeftMargin(0.18);
  const Int_t colour[2] = {kLeftColour, kRightColour};
  // The dashed spectra first, so the corrected ones are drawn over them. The
  // top half of the frame is left for the legend.
  spectra[0][0]->SetMaximum((BigLegend() ? 2.1 : 1.45) * yMax);
  PlottingUtils::ConfigureAndDrawHistogram(spectra[0][0], colour[0],
                                           ";Pulse Height [ADC];Counts");
  spectra[0][0]->GetYaxis()->SetTitleOffset(1.4);
  for (Int_t k = 0; k < 2; k++)
    for (Int_t chain = 0; chain < 2; chain++) {
      if (k == 0 && chain == 0)
        continue;
      PlottingUtils::ConfigureHistogram(spectra[chain][k], colour[chain]);
      spectra[chain][k]->Draw("HIST SAME");
    }
  for (Int_t chain = 0; chain < 2; chain++)
    spectra[chain][0]->SetLineStyle(2);

  // Colour says which chain, the line style before or after correction.
  TLine *dashed = new TLine();
  dashed->SetLineStyle(2);
  dashed->SetLineWidth(spectra[0][1]->GetLineWidth());
  TLine *solid = new TLine();
  solid->SetLineWidth(spectra[0][1]->GetLineWidth());
  std::vector<TString> labels;
  labels.push_back("#tau_{pz} = #tau_{pre}");
  labels.push_back("#tau_{pz} #neq #tau_{pre}");
  labels.push_back("Before Correction");
  labels.push_back("After Correction");
  TLegend *legend =
      Legend(labels, 0.87 - LegendHeight(4), 0.87, "", kFALSE, gLegendPx);
  legend->AddEntry(spectra[0][1], labels[0], "l");
  legend->AddEntry(spectra[1][1], labels[1], "l");
  legend->AddEntry(dashed, labels[2], "l");
  legend->AddEntry(solid, labels[3], "l");
  legend->Draw();
  CornerText(rateLabel);
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

// (Fig. simkernel) The kernels the unmodified MUSIC code fits on one
// simulated run, matched and mismatched chain, against the analytic kernel of
// the paper's Eq. (kerneltail), gamma exp(-(dt + t_m)/tau_pre) / g, with g the
// height of an isolated unit pulse. That is the kernel of an ideal baseline
// once the earlier pulse's trapezoid has ended, at (k + l)T plus the
// simulation's detector ramp after it, so it is drawn from
// dt = (k + l)T + ramp - t_m on; the firmware baseline of the run adds its own
// share to the fitted kernel. For the matched chain gamma = 0 and the analytic
// kernel is the zero line, so it is not drawn. A bin with no pulse pairs
// (closer than the firmware resolves two pulses) has no coefficient and is
// left out.
void CorrectionKernelFigure(const TString &path, const TString &runLabel) {
  std::vector<SimHit> hits;
  RunInfo info;
  if (!ReadHits(path, hits, info))
    return;
  PulseHistoryAdapter::Adapted adapted = PulseHistoryAdapter::Build(hits);
  std::vector<Int_t> groups = PulseHistory::BuildGroupMap();
  PulseHistory::Result result;
  if (!PulseHistory::Measure(adapted.hits, groups, result,
                             "correction kernel " + runLabel)) {
    std::cerr << "PulseHistory::Measure failed for " << path << std::endl;
    return;
  }
  const Double_t kUs = kTrapRiseNs / 1000.0;
  const Double_t lUs = (kTrapRiseNs + kTrapFlatNs) / 1000.0;
  const Double_t tm = (kTrapRiseNs + kPeakingDelayNs) / 1000.0;
  const Double_t dtTail = kUs + lUs + Config().detRise_us - tm;
  const Double_t xLo = 1.0, xHi = 350.0, yLo = -0.45, yHi = 0.15;
  const Int_t nChains = 2;
  const Int_t group[nChains] = {PulseHistory::kLONG_RIGHT,
                                PulseHistory::kLONG_LEFT};
  const Double_t tau[nChains] = {info.tauR_us, info.tauL_us};
  const Int_t colour[nChains] = {kRightColour, kLeftColour};
  Bool_t matched[nChains];
  TString name[nChains];
  for (Int_t j = 0; j < nChains; j++) {
    matched[j] = TMath::Abs(tau[j] - info.tauPz_us) < 1e-9;
    name[j] = matched[j] ? "#tau_{pz} = #tau_{pre}" : "#tau_{pz} #neq #tau_{pre}";
  }

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetLogx();
  TH1F *frame = Frame(canvas, xLo, yLo, xHi, yHi, "Time Between Pulses [#mus]",
                      "Shift / Earlier Pulse Height");
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();
  Reference(xLo, xHi, 0.0);

  std::vector<TGraph *> fitted, form;
  for (Int_t j = 0; j < nChains; j++) {
    const PulseHistory::Kernel &kernel = result.kernel[group[j]];
    TGraph *points = new TGraph();
    for (Int_t b = 0; b < kNBins && kernel.ok; b++)
      if (kernel.k[0][b] != 0.0)
        points->SetPoint(points->GetN(), BinCentre(b), kernel.k[0][b]);
    points->SetMarkerStyle(j == 0 ? 20 : 21);
    points->SetMarkerSize(1.2);
    points->SetMarkerColor(colour[j]);
    points->SetLineColor(colour[j]);
    fitted.push_back(points);

    const Double_t c = TailAmplitude(tau[j], info.tauPz_us, kUs, lUs);
    TGraph *isolated = StepResponse(tau[j], info.tauPz_us,
                                    std::vector<Double_t>(1, 0.0), 0.0,
                                    tm + 1.0);
    const Double_t gain = isolated->Eval(tm);
    delete isolated;
    TGraph *curve = nullptr;
    if (!matched[j]) {
      curve = new TGraph();
      const Int_t nPoints = 200;
      for (Int_t p = 0; p <= nPoints; p++) {
        const Double_t dt =
            dtTail * TMath::Power(xHi / dtTail, Double_t(p) / nPoints);
        curve->SetPoint(p, dt, c * TMath::Exp(-(dt + tm) / tau[j]) / gain);
      }
      PlottingUtils::ConfigureGraph(curve, colour[j]);
      curve->Draw("L SAME");
    }
    form.push_back(curve);
    std::cout << Form("correction kernel %s: gamma = %.3f, g = %.3f, fitted %s",
                      matched[j] ? "matched" : "mismatched", c, gain,
                      kernel.ok ? "ok" : "FAILED")
              << std::endl;
  }
  for (Int_t j = 0; j < nChains; j++)
    fitted[j]->Draw("P SAME");
  canvas->RedrawAxis();

  std::vector<TString> labels;
  std::vector<TObject *> entries;
  std::vector<const char *> options;
  for (Int_t j = 0; j < nChains; j++) {
    labels.push_back(Form("Fitted, %s", name[j].Data()));
    entries.push_back(fitted[j]);
    options.push_back("p");
    if (form[j]) {
      labels.push_back(Form("Analytic, %s", name[j].Data()));
      entries.push_back(form[j]);
      options.push_back("l");
    }
  }
  TLegend *legend = Legend(labels, 0.20, 0.20 + LegendHeight(labels.size()), "",
                           kFALSE, gLegendPx);
  for (size_t e = 0; e < labels.size(); e++)
    legend->AddEntry(entries[e], labels[e], options[e]);
  legend->Draw();
  CornerText(runLabel);
  PlottingUtils::SaveFigure(canvas, "correction_kernel", "",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
}

// (Fig. rates) One panel of the spectra figure, from one simulated run: the
// spectra before and after the PulseHistory correction, with the kernel
// measured on the run itself.
void SpectrumFigure(const TString &path, const TString &rateLabel,
                    const TString &name, const char *panel) {
  std::vector<SimHit> hits;
  RunInfo info;
  if (!ReadHits(path, hits, info))
    return;
  std::vector<Double_t> raw(hits.size(), 0.0);
  for (size_t i = 0; i < hits.size(); i++)
    raw[i] = hits[i].energy;
  std::vector<Double_t> corrected;
  if (!CorrectRun(hits, name, corrected))
    return;
  DrawSpectra(hits, raw, corrected, info.beamPeak, rateLabel, name, panel);
}

// The mean of the values, 0 for none.
Double_t Mean(const std::vector<Double_t> &values) {
  if (values.empty())
    return 0.0;
  Double_t sum = 0.0;
  for (size_t i = 0; i < values.size(); i++)
    sum += values[i];
  return sum / Double_t(values.size());
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

// A chain's beam peak: the mode of its uncorrected spectrum of non-pileup
// hits, as PulseHistory takes each channel's, in bins of 1 % of the nominal
// peak from 0.3 to 2 of it, smoothed over five bins.
Double_t ChainPeak(const std::vector<SimHit> &hits, Int_t chain,
                   Double_t nominal) {
  const Double_t lo = 0.3, width = 0.01;
  const Int_t nBins = 170, half = 2;
  std::vector<Double_t> counts(nBins, 0.0);
  for (const SimHit &h : hits) {
    if (h.channel != chain || h.pileup)
      continue;
    const Int_t b = Int_t(TMath::Floor((h.energy / nominal - lo) / width));
    if (b >= 0 && b < nBins)
      counts[b]++;
  }
  Int_t top = half;
  Double_t best = -1.0;
  for (Int_t b = half; b < nBins - half; b++) {
    Double_t s = 0.0;
    for (Int_t k = -half; k <= half; k++)
      s += counts[b + k];
    if (s > best) {
      best = s;
      top = b;
    }
  }
  return nominal * (lo + (top + 0.5) * width);
}

// The spread of one run, raw and corrected, per chain: non-pileup pulse
// heights within the beam window, each in
// units of its chain's beam peak.
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
  const Double_t peak[2] = {ChainPeak(hits, 0, info.beamPeak),
                            ChainPeak(hits, 1, info.beamPeak)};
  std::vector<Double_t> raw[2], corr[2];
  for (size_t i = 0; i < hits.size(); i++) {
    const SimHit &hit = hits[i];
    if (hit.pileup)
      continue;
    const Double_t lo = kSpreadWindowLo * peak[hit.channel];
    const Double_t hi = kSpreadWindowHi * peak[hit.channel];
    if (hit.energy > lo && hit.energy < hi)
      raw[hit.channel].push_back(hit.energy);
    if (corrected[i] > lo && corrected[i] < hi)
      corr[hit.channel].push_back(corrected[i]);
  }
  for (Int_t ch = 0; ch < 2; ch++) {
    spread.raw[ch] = Rms(raw[ch]);
    spread.corr[ch] = Rms(corr[ch]);
    spread.rawMean[ch] = Mean(raw[ch]);
    spread.corrMean[ch] = Mean(corr[ch]);
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

// (Fig. run70) What run 70's waveforms show of the baseline restorer: the
// mean change of the firmware baseline from one pulse to the next on the
// Long Strips, Right, against the time between them, with CAEN's description
// of LF limitation at the run's own baseline average for three thresholds
// and for no check. scripts/run70_baseline.sh writes the curves.
void Run70BaselineFigure(const TString &path) {
  TFile in(path, "READ");
  TGraphErrors *data = (TGraphErrors *)in.Get("data");
  const Int_t nCurves = 4;
  const char *curveName[nCurves] = {"model_thr10", "model_thr20", "model_thr80",
                                    "model_nocheck"};
  const char *curveLabel[nCurves] = {"10 ADC Threshold", "20 ADC Threshold",
                                     "80 ADC Threshold", "No check"};
  const Int_t curveColour[nCurves] = {kTeal + 3, kOrange + 1, kMagenta + 1,
                                      kGray + 2};
  const Int_t curveStyle[nCurves] = {7, 1, 2, 9};
  TGraph *curves[nCurves];
  Bool_t haveCurves = data != nullptr;
  for (Int_t c = 0; c < nCurves; c++) {
    curves[c] = (TGraph *)in.Get(curveName[c]);
    haveCurves = haveCurves && curves[c];
  }
  if (!haveCurves) {
    std::cerr << "no curves in " << path
              << " (rerun scripts/run70_baseline.sh)" << std::endl;
    return;
  }
  // Room above the curves for the legend's five rows (more for large legends).
  const Double_t xLo = 6.0, xHi = 450.0, yLo = -15.0,
                 yHi = BigLegend() ? 185.0 : 145.0;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  canvas->SetLogx();
  TH1F *frame = Frame(canvas, xLo, yLo, xHi, yHi, "Time Between Pulses [#mus]",
                      "#LT Baseline Change #GT [ADC]");
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();
  Reference(xLo, xHi, 0.0);

  for (Int_t c = 0; c < nCurves; c++) {
    PlottingUtils::ConfigureGraph(curves[c], curveColour[c]);
    curves[c]->SetLineStyle(curveStyle[c]);
    curves[c]->Draw("L SAME");
  }
  PlottingUtils::ConfigureGraph(data, kRightColour);
  data->SetMarkerStyle(20);
  data->SetMarkerColor(kRightColour);
  data->SetMarkerSize(1.2);
  data->Draw("P SAME");
  canvas->RedrawAxis();

  std::vector<TString> labels;
  labels.push_back("Diagnostic Run");
  for (Int_t c = 0; c < nCurves; c++)
    labels.push_back(curveLabel[c]);
  TLegend *legend =
      Legend(labels, 0.87 - LegendHeight(5), 0.87, "", kTRUE, gLegendPx);
  legend->AddEntry(data, labels[0], "PE");
  for (Int_t c = 0; c < nCurves; c++)
    legend->AddEntry(curves[c], labels[c + 1], "L");
  legend->Draw();
  PlottingUtils::SaveFigure(canvas, "run70_baseline_change", "",
                            PlotSaveOptions::kLINEAR);
  delete canvas;
}

// The simulated rms spread of each baseline reading at each rate, before and
// after correction, in percent of the chain's mean height, for the mismatched and the
// matched chain, as a booktabs tabular for the paper to \input inside its own
// table float, after the caption and label (tab:spread); needs siunitx and
// amsmath. What a caption needs from the runs is written as comments: the
// decay times and amplitude spread are those of the runs, which main checks
// agree; spreads holds every reading at every rate.
void SpreadTable(const ReadingSpreads &spreads, const RunInfo &runs,
                 const TString &path) {
  std::ofstream table(path.Data());
  table << "% Written by mains/figures.cpp; edit that, not this file."
        << std::endl;
  table << Form("%% For the caption: rms of the pulse heights within %g--%g "
                "times the chain's beam peak (the mode of its uncorrected "
                "spectrum), pileup left out, in percent of their mean; beam amplitudes with a %g %% Gaussian spread; after each "
                "trigger the baseline is frozen for %g us; at the end of each "
                "freeze it resumes if the trapezoid is within %g ADC of it and "
                "otherwise freezes again, the trapezoid taken as its running "
                "average over the last %d samples, kept up during the freeze "
                "(running-average check), or its last sample (single-sample "
                "check); tau_pre = %g us (matched), %g us (mismatched), tau_pz "
                "= %g us.",
                kSpreadWindowLo, kSpreadWindowHi,
                100.0 * runs.sigmaFrac, kBaselineHoldoffNs / 1000.0,
                kRestartThrADC, kCheckSamples, runs.tauL_us, runs.tauR_us,
                runs.tauPz_us)
        << std::endl;
  const Int_t nColumns = 1 + 2 * kNSpreadRates;
  // Scoped by the paper's table float, as the other generated tables.
  table << "\\small" << std::endl;
  table << "\\setlength{\\tabcolsep}{4pt}" << std::endl;
  table << "\\begin{tabular}{@{}l";
  for (Int_t i = 0; i < 2 * kNSpreadRates; i++)
    table << "c";
  table << "@{}}" << std::endl;
  table << "\\toprule" << std::endl;
  for (Int_t i = 0; i < kNSpreadRates; i++)
    table << Form(" & \\multicolumn{2}{c}{\\SI{%g}{\\kilo\\hertz}}",
                  kSpreadRates[i]);
  table << " \\\\" << std::endl;
  for (Int_t i = 0; i < kNSpreadRates; i++)
    table << Form("\\cmidrule(lr){%d-%d} ", 2 + 2 * i, 3 + 2 * i);
  table << std::endl;
  table << "Baseline";
  for (Int_t i = 0; i < kNSpreadRates; i++)
    table << " & before & after";
  table << " \\\\" << std::endl;

  const Int_t nSides = 2;
  const Int_t sides[nSides] = {1, 0};
  const char *sideNames[nSides] = {"$\\tau_{\\text{pz}} \\neq \\tau_{\\text{pre}}$",
                                   "$\\tau_{\\text{pz}} = \\tau_{\\text{pre}}$"};
  for (Int_t s = 0; s < nSides; s++) {
    table << "\\midrule" << std::endl;
    table << Form("\\multicolumn{%d}{l}{%s} \\\\", nColumns, sideNames[s])
          << std::endl;
    for (Int_t j = 0; j < kNSpreadReadings; j++) {
      const RateSpreads &byRate = spreads.find(kSpreadReadings[j].key)->second;
      table << kSpreadReadings[j].label;
      for (Int_t i = 0; i < kNSpreadRates; i++) {
        const Spread &spread = byRate.find(kSpreadRates[i])->second;
        table << Form(" & %.2f & %.2f", spread.RawPercent(sides[s]),
                      spread.CorrPercent(sides[s]));
      }
      table << " \\\\" << std::endl;
    }
  }
  table << "\\bottomrule" << std::endl;
  table << "\\end{tabular}" << std::endl;
  std::cout << "wrote " << path << std::endl;
}

// The amplitude-dependent kernels: the correction of run 97 with one kernel
// per band of the previous pulse height, bands [0, 0.5), [0.5, 1.5), ...,
// [bands - 1.5, inf) times the beam peak, against the production single
// kernel.
//
// What scripts/run97_figure.sh measures for the band macros: the long groups'
// spread before correction, with one band and with `bands` bands, rms in
// percent of the mean height as the simulated spreads, [0] left and [1] right
// as the simulated chains; and the long ends' trigger rate.
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
    std::cout << "skipping band_macros: no band spreads in " << file
              << " (rerun scripts/run97_figure.sh)" << std::endl;
    return kFALSE;
  }
  for (Int_t c = 0; c < 2; c++)
    for (Int_t k = 0; k < 3; k++)
      data.spread[c][k] = (*chain[c])[k];
  data.bands = bands->GetVal();
  data.rate_kHz = rate->GetVal();
  return kTRUE;
}

// Run 97's Long Strips, Right as LaTeX macros for the paper (Section 5, the
// abstract and the conclusion), in percent of the mean height: the spread
// before correction, with one kernel and with one per band, pooled over the
// strips; the number of bands and the long ends' rate.
void BandMacros(const DataBands &data, const TString &path) {
  std::ofstream macros(path.Data());
  macros << "% Written by mains/figures.cpp; edit that, not this file."
         << std::endl;
  macros << "% Spreads are rms in percent of the mean height. Run 97: Long "
            "Strips, Right, each strip in units of its own beam peak, "
            "corrected with the kernels fitted on it. Bins of the previous "
            "pulse height, in beam peaks: below 0.5, ";
  for (Int_t b = 1; b < data.bands - 1; b++)
    macros << Form("%g--%g, ", b - 0.5, b + 0.5);
  macros << Form("above %g.", data.bands - 1.5) << std::endl;
  macros << Form("\\newcommand{\\BandCount}{%d}", data.bands) << std::endl;
  macros << Form("\\newcommand{\\BandRunRate}{%.0f}", data.rate_kHz)
         << std::endl;
  macros << Form("\\newcommand{\\BandRunBefore}{%.2f}", data.spread[1][0])
         << std::endl;
  macros << Form("\\newcommand{\\BandRunOne}{%.2f}", data.spread[1][1])
         << std::endl;
  macros << Form("\\newcommand{\\BandRunMany}{%.2f}", data.spread[1][2])
         << std::endl;
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
// out. Drawn to kDeviationMaxUs only: past it a long gap to the previous pulse
// also means a quiet stretch with fewer earlier pulses of any age, which the
// kernel sums over, so this one-neighbour profile no longer shows what the
// correction does there.
const Double_t kMinProfilePulses = 200.0;
const Double_t kDeviationMaxUs = 40.0;

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
  const Double_t xLo = 1.0, xHi = kDeviationMaxUs, yLo = -150.0, yHi = 300.0;
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
// strip's beam peak, for the segmented strips 1-16, as panels (a) and (b) of
// one paper figure. At the axis text size a legend of the sixteen strips does
// not fit the frame, so the strips are coloured by parity, as the paper's text
// compares them: odd strips, whose long segment is on the left, in the left
// chain's colour, and even strips, whose long segment is on the right, in the
// right chain's. The top of the log scale leaves room for the legend.
void StripSpectraFigure(const TString &file, Bool_t corrected) {
  const Int_t firstStrip = 1;
  const Int_t nStrips = 16;
  std::vector<TH1D *> spectra;
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
    }
  }
  // The log axis goes on a pad spanning the canvas, with the canvas's own
  // margins, grid and ticks: saving with the linear option, as the paper's
  // file names need, resets the canvas's own log axis.
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  TPad *plot = new TPad(PlottingUtils::GetRandomName(), "", 0.0, 0.0, 1.0, 1.0);
  plot->SetMargin(canvas->GetLeftMargin(), canvas->GetRightMargin(),
                  canvas->GetBottomMargin(), canvas->GetTopMargin());
  plot->SetGrid(canvas->GetGridx(), canvas->GetGridy());
  plot->SetTicks(canvas->GetTickx(), canvas->GetTicky());
  plot->SetLogy();
  plot->Draw();
  plot->cd();
  Double_t yMax = 0.0;
  for (Int_t s = 0; s < nStrips; s++)
    yMax = TMath::Max(yMax, spectra[s]->GetMaximum());
  spectra[0]->SetMaximum((BigLegend() ? 100.0 : 3.0) * yMax);
  spectra[0]->SetMinimum(0.5);
  const Int_t parityColour[2] = {kRightColour, kLeftColour}; // even, odd
  for (Int_t s = 0; s < nStrips; s++) {
    const Int_t odd = (firstStrip + s) % 2;
    if (s == 0)
      PlottingUtils::ConfigureAndDrawHistogram(spectra[s], parityColour[odd],
                                               ";#DeltaE [a.u.];Counts");
    else {
      PlottingUtils::ConfigureHistogram(spectra[s], parityColour[odd]);
      spectra[s]->Draw("HIST SAME");
    }
  }
  plot->RedrawAxis();
  std::vector<TString> labels;
  labels.push_back("Odd Strips");
  labels.push_back("Even Strips");
  TLegend *legend =
      Legend(labels, 0.87 - LegendHeight(2), 0.87, "", kFALSE, gLegendPx);
  legend->AddEntry(spectra[0], labels[0], "l"); // S1
  legend->AddEntry(spectra[1], labels[1], "l"); // S2
  legend->Draw();
  CornerText(corrected ? "After Correction" : "Before Correction");
  PanelLabel(corrected ? "(b)" : "(a)");
  PlottingUtils::SaveFigure(
      canvas, corrected ? "run97_spectra_corrected" : "run97_spectra", "",
      PlotSaveOptions::kLINEAR);
  delete canvas;
}

// ------------------------------------------------- reaction events (37Cl)
// What scripts/reaction_events.sh collects from MUSIC's own outputs: the
// reaction events the production analysis tagged and its pure-beam sample,
// with and without the correction (the latter from scripts/nocorrection.sh,
// MUSIC's pipeline with the correction off). The figure shows a few events of
// each, of one reaction strip, the one MUSIC's own strip-sum plot shows.
const Int_t kReactionStrip = 3;
const Int_t kShownEvents = 10;
const Int_t kReactionColour = kOrange + 7;
const Int_t kBeamColour = kGray + 2;

// The events tagged in reaction strip `reac`, or the pure beam for reac = -1:
// each strip's calibrated deposit before ([0]) or after ([1]) the
// correction. Events that are not complete without the correction are left
// out of both, so the two hold the same events.
Bool_t ReadReactionTraces(const TString &file, Int_t reac,
                          std::vector<std::vector<Float_t>> traces[2]) {
  TFile in(file, "READ");
  TTree *tree = in.IsZombie() ? nullptr : (TTree *)in.Get("events");
  if (!tree) {
    std::cerr << "no events tree in " << file
              << " (run scripts/reaction_events.sh)" << std::endl;
    return kFALSE;
  }
  const Int_t nStrips = 18;
  Float_t raw[nStrips], corrected[nStrips];
  Int_t strip = 0, rawFound = 0;
  tree->SetBranchAddress("raw", raw);
  tree->SetBranchAddress("corrected", corrected);
  tree->SetBranchAddress("reac", &strip);
  tree->SetBranchAddress("raw_found", &rawFound);
  for (Int_t v = 0; v < 2; v++)
    traces[v].clear();
  for (Long64_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if (!rawFound || strip != reac)
      continue;
    traces[0].push_back(std::vector<Float_t>(raw, raw + nStrips));
    traces[1].push_back(std::vector<Float_t>(corrected, corrected + nStrips));
  }
  return !traces[1].empty();
}

// kShownEvents of a set of traces, spread evenly through it (and so through
// the runs), each as a line over the strips. Returns the last, for the
// legend.
TGraph *DrawTraces(const std::vector<std::vector<Float_t>> &traces,
                   Int_t firstStrip, Int_t lastStrip, Int_t colour) {
  const size_t n = traces.size();
  const size_t shown = TMath::Min(size_t(kShownEvents), n);
  TGraph *event = nullptr;
  for (size_t k = 0; k < shown; k++) {
    const std::vector<Float_t> &trace = traces[k * n / shown];
    event = new TGraph();
    for (Int_t s = firstStrip; s <= lastStrip; s++)
      event->SetPoint(event->GetN(), s, trace[s]);
    event->SetLineColor(colour);
    event->SetLineWidth(PlottingUtils::GetLineWidth());
    event->Draw("L SAME");
  }
  return event;
}

// (Fig. reactions) kShownEvents of the events tagged in reaction strip
// kReactionStrip and kShownEvents of MUSIC's pure-beam events, strip by
// strip, before or after the correction: the same events in both.
void ReactionTracesFigure(const TString &file, Bool_t corrected) {
  std::vector<std::vector<Float_t>> reactions[2], beam[2];
  if (!ReadReactionTraces(file, kReactionStrip, reactions) ||
      !ReadReactionTraces(file, -1, beam))
    return;
  const Int_t v = corrected ? 1 : 0;
  const Int_t firstStrip = 1, lastStrip = 16;
  // Room above the traces for the legend (more for large legends).
  const Double_t yLo = 0.2, yHi = BigLegend() ? 1.8 : 1.6;

  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, firstStrip - 0.5, yLo, lastStrip + 0.5, yHi, "Strip",
        "#DeltaE [a.u.]");
  TGraph *beamLine = DrawTraces(beam[v], firstStrip, lastStrip, kBeamColour);
  TGraph *reactionLine =
      DrawTraces(reactions[v], firstStrip, lastStrip, kReactionColour);
  TimeMark(kReactionStrip, yLo, yHi);
  canvas->RedrawAxis();

  std::vector<TString> labels;
  labels.push_back("Beam");
  labels.push_back(Form("(#alpha,n) in Strip %d", kReactionStrip));
  TLegend *legend =
      Legend(labels, 0.87 - LegendHeight(2), 0.87, "", kFALSE, gLegendPx);
  legend->AddEntry(beamLine, labels[0], "l");
  legend->AddEntry(reactionLine, labels[1], "l");
  legend->Draw();
  CornerText(corrected ? "After Correction" : "Before Correction");
  PanelLabel(corrected ? "(b)" : "(a)");
  PlottingUtils::SaveFigure(
      canvas, corrected ? "reaction_traces_after" : "reaction_traces_before",
      "", PlotSaveOptions::kLINEAR);
  delete canvas;
  std::cout << Form("reaction traces (%s): %d of %zu events in strip %d, %d "
                    "of %zu beam events",
                    corrected ? "after" : "before",
                    Int_t(TMath::Min(size_t(kShownEvents), reactions[v].size())),
                    reactions[v].size(), kReactionStrip,
                    Int_t(TMath::Min(size_t(kShownEvents), beam[v].size())),
                    beam[v].size())
            << std::endl;
}

// The strips summed for the scatter's Delta E: the reaction strip and the
// three after it, where the products' higher deposit shows.
const Int_t kScatterFirstStrip = kReactionStrip;
const Int_t kScatterLastStrip = kReactionStrip + 3;

// (Fig. reaction scatter) Every event tagged in reaction strip
// kReactionStrip and MUSIC's pure-beam events, Delta E summed over strips
// kScatterFirstStrip-kScatterLastStrip against E summed over all strips,
// before or after the correction: the same events in both.
void ReactionScatterFigure(const TString &file, Bool_t corrected) {
  std::vector<std::vector<Float_t>> reactions[2], beam[2];
  if (!ReadReactionTraces(file, kReactionStrip, reactions) ||
      !ReadReactionTraces(file, -1, beam))
    return;
  const Int_t v = corrected ? 1 : 0;
  auto Scatter = [&](const std::vector<std::vector<Float_t>> &events,
                     Int_t colour, Int_t style) {
    TGraph *graph = new TGraph();
    for (const std::vector<Float_t> &event : events) {
      Double_t dE = 0.0, e = 0.0;
      for (size_t s = 0; s < event.size(); s++) {
        e += event[s];
        if (Int_t(s) >= kScatterFirstStrip && Int_t(s) <= kScatterLastStrip)
          dE += event[s];
      }
      graph->SetPoint(graph->GetN(), e, dE);
    }
    graph->SetMarkerStyle(style);
    graph->SetMarkerSize(0.6);
    graph->SetMarkerColor(colour);
    graph->SetLineColor(colour);
    return graph;
  };
  const Double_t xLo = 13.0, xHi = 20.0, yLo = 3.0, yHi = 6.0;
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  Frame(canvas, xLo, yLo, xHi, yHi, "#it{E} [a.u.]",
        Form("#DeltaE_{%d-%d} [a.u.]", kScatterFirstStrip, kScatterLastStrip));
  TGraph *reactionPoints = Scatter(reactions[v], kReactionColour, 20);
  TGraph *beamPoints = Scatter(beam[v], kBeamColour, 21);
  reactionPoints->Draw("P SAME");
  beamPoints->Draw("P SAME");
  canvas->RedrawAxis();

  std::vector<TString> labels;
  labels.push_back("Beam");
  labels.push_back(Form("(#alpha,n) in Strip %d", kReactionStrip));
  TLegend *legend = Legend(labels, 0.73, 0.87);
  legend->AddEntry(beamPoints, labels[0], "p");
  legend->AddEntry(reactionPoints, labels[1], "p");
  legend->Draw();
  PlottingUtils::AddText(corrected ? "After Correction" : "Before Correction",
                         0.88, 0.92);
  PanelLabel(corrected ? "(b)" : "(a)");
  PlottingUtils::SaveFigure(
      canvas, corrected ? "reaction_scatter_after" : "reaction_scatter_before",
      "", PlotSaveOptions::kLINEAR);
  delete canvas;
  std::cout << Form("reaction scatter (%s): %zu events in strip %d, %zu beam "
                    "events",
                    corrected ? "after" : "before", reactions[v].size(),
                    kReactionStrip, beam[v].size())
            << std::endl;
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

// ------------------------------------------------- the first-order check
// Section 3.3: PulseHistory fits and applies the kernel with the measured
// heights of the earlier pulses, a first-order stand-in for the corrected
// heights of Eq. (kernel). What that costs is checked here by redoing the
// binned kernel fit as plain least squares, with an intercept, of each seed's
// height less the seeds' mean on the summed earlier heights in each bin of the
// time since them (PulseHistory's bins: 12 per decade from 1 us, kNBins of
// them), on one chain of a simulated run, three ways:
//   first order - fitted and applied with the measured earlier heights
//   exact       - fitted and applied with the true ones, the chain's gain times
//                 the true amplitude, which only a simulation knows
//   two passes  - fitted and applied with the heights the first order corrected
// Seeds are the non-pileup hits whose event has both chains within
// kSeedLo-kSeedHi of their beam peak, like PulseHistory's beam-like events;
// the same window as the spreads.
const Double_t kSeedLo = kSpreadWindowLo;
const Double_t kSeedHi = kSpreadWindowHi;
const Int_t kNFirstOrder = 3;

// One chain's hits in time order, as indices into hits.
std::vector<size_t> ChainHits(const std::vector<SimHit> &hits, Int_t chain) {
  std::vector<size_t> index;
  for (size_t i = 0; i < hits.size(); i++)
    if (hits[i].channel == chain)
      index.push_back(i);
  std::sort(index.begin(), index.end(), [&](size_t a, size_t b) {
    return hits[a].timestamp < hits[b].timestamp;
  });
  return index;
}

// The kernel bin of the time dt_us since an earlier pulse, -1 outside it.
Int_t KernelBin(Double_t dt_us) {
  if (!(dt_us >= 1.0))
    return -1;
  const Int_t bin = Int_t(12.0 * TMath::Log10(dt_us));
  return bin < kNBins ? bin : -1;
}

// The kernel fitted on one chain: the seeds' measured heights against the
// summed `history` heights of the chain's earlier hits in each bin. Bins no
// pulse pair falls in are left at zero, as PulseHistory leaves them out.
std::vector<Double_t> FitChainKernel(const std::vector<SimHit> &hits,
                                     const std::vector<size_t> &chain,
                                     const std::vector<Bool_t> &seed,
                                     const std::vector<Double_t> &history) {
  const Double_t reachUs = TMath::Power(10.0, kNBins / 12.0);
  const Int_t np = 1 + kNBins;
  Double_t sum = 0.0, n = 0.0;
  for (size_t p = 0; p < chain.size(); p++)
    if (seed[chain[p]]) {
      sum += hits[chain[p]].energy;
      n++;
    }
  const Double_t mean = sum / n;
  TMatrixD A(np, np);
  TVectorD b(np);
  std::vector<Double_t> x(np);
  for (size_t p = 0; p < chain.size(); p++) {
    const size_t j = chain[p];
    if (!seed[j])
      continue;
    std::fill(x.begin(), x.end(), 0.0);
    x[0] = 1.0;
    for (size_t q = p; q-- > 0;) {
      const size_t i = chain[q];
      const Double_t dt = Double_t(hits[j].timestamp - hits[i].timestamp) * 1e-6;
      if (dt >= reachUs)
        break;
      const Int_t bin = KernelBin(dt);
      if (bin >= 0)
        x[1 + bin] += history[i];
    }
    const Double_t y = hits[j].energy - mean;
    for (Int_t a = 0; a < np; a++) {
      if (x[a] == 0.0)
        continue;
      b[a] += x[a] * y;
      for (Int_t c = 0; c < np; c++)
        A(a, c) += x[a] * x[c];
    }
  }
  // Least squares on the bins with pairs: unit diagonal and a negligible
  // ridge, as PulseHistory solves it.
  std::vector<Int_t> keep;
  for (Int_t a = 0; a < np; a++)
    if (A(a, a) > 0.0)
      keep.push_back(a);
  const Int_t nk = Int_t(keep.size());
  TMatrixD As(nk, nk);
  TVectorD bs(nk), scale(nk);
  for (Int_t a = 0; a < nk; a++)
    scale[a] = 1.0 / TMath::Sqrt(A(keep[a], keep[a]));
  for (Int_t a = 0; a < nk; a++) {
    bs[a] = b[keep[a]] * scale[a];
    for (Int_t c = 0; c < nk; c++)
      As(a, c) = A(keep[a], keep[c]) * scale[a] * scale[c];
    As(a, a) += 1e-8;
  }
  TDecompChol chol(As);
  Bool_t ok = chol.Decompose();
  TVectorD p = ok ? chol.Solve(bs, ok) : TVectorD(nk);
  std::vector<Double_t> kernel(kNBins, 0.0);
  for (Int_t a = 0; a < nk && ok; a++)
    if (keep[a] > 0)
      kernel[keep[a] - 1] = p[a] * scale[a];
  if (!ok)
    std::cerr << "first-order check: the kernel fit failed" << std::endl;
  return kernel;
}

// The chain's heights less the kernel times the `history` heights of its
// earlier hits; other hits keep their measured height.
std::vector<Double_t> ApplyChainKernel(const std::vector<SimHit> &hits,
                                       const std::vector<size_t> &chain,
                                       const std::vector<Double_t> &kernel,
                                       const std::vector<Double_t> &history) {
  const Double_t reachUs = TMath::Power(10.0, kNBins / 12.0);
  std::vector<Double_t> corrected(hits.size());
  for (size_t i = 0; i < hits.size(); i++)
    corrected[i] = hits[i].energy;
  for (size_t p = 0; p < chain.size(); p++) {
    const size_t j = chain[p];
    Double_t shift = 0.0;
    for (size_t q = p; q-- > 0;) {
      const size_t i = chain[q];
      const Double_t dt = Double_t(hits[j].timestamp - hits[i].timestamp) * 1e-6;
      if (dt >= reachUs)
        break;
      const Int_t bin = KernelBin(dt);
      if (bin >= 0)
        shift += kernel[bin] * history[i];
    }
    corrected[j] -= shift;
  }
  return corrected;
}

// The spread of one chain as the spread table measures it: non-pileup heights
// within the beam window around the chain's beam peak, rms in percent of their
// mean.
Double_t ChainSpreadPercent(const std::vector<SimHit> &hits,
                            const std::vector<Double_t> &heights, Int_t chain,
                            Double_t peak) {
  std::vector<Double_t> values;
  for (size_t i = 0; i < hits.size(); i++)
    if (hits[i].channel == chain && !hits[i].pileup &&
        heights[i] > kSpreadWindowLo * peak &&
        heights[i] < kSpreadWindowHi * peak)
      values.push_back(heights[i]);
  const Double_t mean = Mean(values);
  return mean > 0.0 ? 100.0 * Rms(values) / mean : 0.0;
}

// The three corrections of the mismatched chain of one run: spreads before
// and after, and each kernel's mean ratio to the analytic one over
// 4.5-20 us, where Eq. (kerneltail) holds.
struct FirstOrderResult {
  Double_t before = 0.0;
  Double_t matched = 0.0;
  Double_t spread[kNFirstOrder] = {0.0, 0.0, 0.0};
  Double_t ratio[kNFirstOrder] = {0.0, 0.0, 0.0};
};

Bool_t FirstOrderCheck(const TString &path, FirstOrderResult &result) {
  std::vector<SimHit> hits;
  RunInfo info;
  if (!ReadHits(path, hits, info))
    return kFALSE;
  // Each chain's gain from hits with nothing on it for 100 us before, by then
  // ten mismatched decay times: measured over true height.
  Double_t gain[2];
  for (Int_t c = 0; c < 2; c++) {
    Double_t se = 0.0, sa = 0.0;
    for (const SimHit &h : hits)
      if (h.channel == c && !h.pileup && h.dtPrev > 100.0) {
        se += h.energy;
        sa += h.trueAmp;
      }
    gain[c] = sa > 0.0 ? se / sa : 0.0;
  }
  // Each chain's beam peak, for the seed and spread windows.
  const Double_t chainPeak[2] = {ChainPeak(hits, 0, info.beamPeak),
                                 ChainPeak(hits, 1, info.beamPeak)};
  // Seeds: the event's hits on both chains inside their window.
  Long64_t nEvents = 0;
  for (const SimHit &h : hits)
    nEvents = TMath::Max(nEvents, h.eventIndex + 1);
  std::vector<Double_t> eventHeight[2] = {
      std::vector<Double_t>(nEvents, -1.0), std::vector<Double_t>(nEvents, -1.0)};
  for (const SimHit &h : hits)
    eventHeight[h.channel][h.eventIndex] = h.energy;
  std::vector<Bool_t> seed(hits.size(), kFALSE);
  for (size_t i = 0; i < hits.size(); i++) {
    if (hits[i].pileup)
      continue;
    Bool_t inside = kTRUE;
    for (Int_t c = 0; c < 2; c++) {
      const Double_t peak = chainPeak[c];
      const Double_t e = eventHeight[c][hits[i].eventIndex];
      inside = inside && e > kSeedLo * peak && e < kSeedHi * peak;
    }
    seed[i] = inside;
  }

  std::vector<Double_t> measured(hits.size()), truth(hits.size());
  for (size_t i = 0; i < hits.size(); i++) {
    measured[i] = hits[i].energy;
    truth[i] = gain[hits[i].channel] * hits[i].trueAmp;
  }
  const Int_t mismatched = 1, matched = 0;
  const std::vector<size_t> chain = ChainHits(hits, mismatched);
  result.before =
      ChainSpreadPercent(hits, measured, mismatched, chainPeak[mismatched]);
  const std::vector<size_t> matchedChain = ChainHits(hits, matched);
  result.matched = ChainSpreadPercent(
      hits,
      ApplyChainKernel(hits, matchedChain,
                       FitChainKernel(hits, matchedChain, seed, measured),
                       measured),
      matched, chainPeak[matched]);

  // The analytic kernel of the mismatched chain, as the kernel figure has it.
  const Double_t kUs = kTrapRiseNs / 1000.0;
  const Double_t lUs = (kTrapRiseNs + kTrapFlatNs) / 1000.0;
  const Double_t tm = (kTrapRiseNs + kPeakingDelayNs) / 1000.0;
  const Double_t tau = info.tauR_us;
  const Double_t gamma = TailAmplitude(tau, info.tauPz_us, kUs, lUs);
  TGraph *isolated = StepResponse(tau, info.tauPz_us,
                                  std::vector<Double_t>(1, 0.0), 0.0, tm + 1.0);
  const Double_t g = isolated->Eval(tm);
  delete isolated;

  std::vector<Double_t> firstOrder;
  for (Int_t v = 0; v < kNFirstOrder; v++) {
    const std::vector<Double_t> &history =
        v == 0 ? measured : (v == 1 ? truth : firstOrder);
    const std::vector<Double_t> kernel =
        FitChainKernel(hits, chain, seed, history);
    const std::vector<Double_t> corrected =
        ApplyChainKernel(hits, chain, kernel, history);
    if (v == 0)
      firstOrder = corrected;
    result.spread[v] =
        ChainSpreadPercent(hits, corrected, mismatched, chainPeak[mismatched]);
    Double_t sum = 0.0, n = 0.0;
    for (Int_t b = 0; b < kNBins; b++) {
      const Double_t dt = BinCentre(b);
      if (dt < 4.5 || dt > 20.0)
        continue;
      sum += kernel[b] / (gamma / g * TMath::Exp(-(dt + tm) / tau));
      n++;
    }
    result.ratio[v] = sum / n;
  }
  std::cout << Form("first-order check %s: gains L %.3f R %.3f; mismatched "
                    "spread before %.2f %%, first order %.2f %%, exact %.2f %%, "
                    "two passes %.2f %% (matched %.2f %%); kernel / analytic at "
                    "4.5-20 us: %.3f, %.3f, %.3f",
                    gSystem->BaseName(path), gain[0], gain[1], result.before,
                    result.spread[0], result.spread[1], result.spread[2],
                    result.matched, result.ratio[0], result.ratio[1],
                    result.ratio[2])
            << std::endl;
  return kTRUE;
}

// The first-order check of the figures' reading at both rates, printed: what
// Section 4.1 states in words, that fitting with the measured rather than the
// true earlier heights moves the rms by a fraction of a percent, and that a
// second pass does not bring the kernel closer to the analytic one.
void FirstOrderReport(const TString &simDir) {
  for (Int_t r = 0; r < kNSpreadRates; r++) {
    FirstOrderResult result;
    FirstOrderCheck(SimRun(simDir, kSpreadReadings[0].key, kSpreadRates[r]),
                    result);
  }
}

int main(int argc, char **argv) {
  TString out = "doc/figs";
  TString simDir = "output/sim";
  TString run97 = "output/run97/run97.root";
  TString reactionEvents = "output/reaction_events/reaction_events.root";
  TString run70 = "output/run70/run70_baseline.root";
  PlotSaveFormat format = PlotSaveFormat::kPNGAndPDF;
  for (Int_t i = 1; i < argc; i++) {
    TString argument = argv[i];
    if (argument.BeginsWith("--out="))
      out = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--sim="))
      simDir = argument(6, argument.Length() - 6);
    else if (argument.BeginsWith("--run97="))
      run97 = argument(8, argument.Length() - 8);
    else if (argument.BeginsWith("--reaction-events="))
      reactionEvents = argument(18, argument.Length() - 18);
    else if (argument.BeginsWith("--run70="))
      run70 = argument(8, argument.Length() - 8);
    else if (argument.BeginsWith("--annotation-px="))
      gAnnotationPx = TString(argument(16, argument.Length() - 16)).Atof();
    else if (argument.BeginsWith("--legend-px="))
      gLegendPx = TString(argument(12, argument.Length() - 12)).Atof();
    else if (argument == "--format=png")
      format = PlotSaveFormat::kPNG;
    else if (argument == "--format=pdf")
      format = PlotSaveFormat::kPDF;
  }
  InitUtils::SetROOTPreferences(format, out, ".", kFALSE);

  TrapezoidFigure();
  PoleZeroFigure();
  SinglePulseFigure();
  CorrectionTrainFigure();
  // The simulated figures show the first baseline reading: the spectra at both
  // rates, the kernels at the higher one.
  const char *shown = kSpreadReadings[0].key;
  const Double_t highRate = kSpreadRates[kNSpreadRates - 1];
  const TString high = SimRun(simDir, shown, highRate);
  if (HaveInput(high, "correction_kernel"))
    CorrectionKernelFigure(high, Form("%g kHz", highRate));
  for (Int_t i = 0; i < kNSpreadRates; i++) {
    const TString path = SimRun(simDir, shown, kSpreadRates[i]);
    const TString name = Form("spectrum_%gkHz", kSpreadRates[i]);
    if (HaveInput(path, name))
      SpectrumFigure(path, Form("%g kHz", kSpreadRates[i]), name,
                     i == 0 ? "(a)" : "(b)");
  }

  if (HaveInput(run97, "kernel_run97")) {
    std::vector<std::vector<Double_t>> kernels = ReadKernels(run97);
    if (kernels.size() < 3)
      std::cerr << "no kernels in " << run97
                << " (rerun scripts/run97_figure.sh)" << std::endl;
    else
      KernelFigure(kernels);
  }

  if (format != PlotSaveFormat::kPNG) {
    ReadingSpreads spreads;
    RunInfo runs;
    Bool_t haveSpreads = kTRUE;
    for (Int_t j = 0; j < kNSpreadReadings && haveSpreads; j++)
      for (Int_t i = 0; i < kNSpreadRates && haveSpreads; i++) {
        const TString path = SimRun(simDir, kSpreadReadings[j].key,
                                    kSpreadRates[i]);
        RunInfo info;
        haveSpreads =
            HaveInput(path, "spread_table") &&
            MeasureSpread(path, gSystem->BaseName(path),
                          spreads[kSpreadReadings[j].key][kSpreadRates[i]],
                          info);
        if (!haveSpreads)
          break;
        if (j == 0 && i == 0)
          runs = info;
        else if (info.tauL_us != runs.tauL_us || info.tauR_us != runs.tauR_us ||
                 info.tauPz_us != runs.tauPz_us ||
                 info.sigmaFrac != runs.sigmaFrac) {
          std::cout << "skipping spread_table: " << path
                    << " has other decay times or amplitude spread than "
                    << SimRun(simDir, kSpreadReadings[0].key, kSpreadRates[0])
                    << std::endl;
          haveSpreads = kFALSE;
        }
      }
    if (haveSpreads)
      SpreadTable(spreads, runs, out + "/spread_table.tex");

    // What the first-order approximation of Section 3.3 costs.
    if (haveSpreads)
      FirstOrderReport(simDir);
  }

  if (HaveInput(run97, "run97 figures")) {
    // Run 97 corrected with one kernel per amplitude band, against one kernel.
    DataBands dataBands;
    if (ReadDataBands(run97, dataBands))
      BandMacros(dataBands, out + "/band_macros.tex");
    DeviationFigure(run97);
    StripSpectraFigure(run97, kFALSE);
    StripSpectraFigure(run97, kTRUE);
  }
  if (HaveInput(reactionEvents, "reaction figures"))
    for (Int_t corrected = 0; corrected < 2; corrected++) {
      ReactionTracesFigure(reactionEvents, corrected);
      ReactionScatterFigure(reactionEvents, corrected);
    }
  if (HaveInput(run70, "run70_baseline_change"))
    Run70BaselineFigure(run70);
  if (format == PlotSaveFormat::kPNGAndPDF)
    MovePngs(out);
  return 0;
}
