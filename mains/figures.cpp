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
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TPad.h>
#include <TParameter.h>
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
  const Double_t yLo = -0.45, yHi = 0.45;
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
    graphs[v]->Draw("L SAME");
    labels.push_back(
        Form("#tau_{pz} = %.0f #mus,  %s", tauPz[v],
             TailLabel(TailAmplitude(tauPre, tauPz[v], kUs, lUs)).Data()));
  }
  canvas->RedrawAxis();

  Note(endUs + 1.55, 0.31, "over-compensated", 10, colour[1]);
  Note(endUs + 4, -0.14, "under-compensated", 13, colour[0]);
  Note(endUs + 0.4, yLo + 0.06, "(#it{k}+#it{l})#it{T}", 10);

  TLegend *legend = Legend(labels, 0.16, 0.27);
  for (Int_t v = 0; v < nCases; v++)
    legend->AddEntry(graphs[v], labels[v], "l");
  legend->Draw();
  PlottingUtils::AddText("#tau_{pre} = 10 #mus", 0.87, 0.84);
  PanelLabel("(b)");

  PlottingUtils::SaveFigure(canvas, "polezero", "", PlotSaveOptions::kLINEAR);
  delete canvas;
}

// (Fig. correction a) One pulse on the right chain (tau_pre = 10 us) through the 37Cl
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
        "Filter Output #hat{#it{s}} [a.u.]");
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
  TimeArrow(0.0, tm, 1.42, 1.48, "#it{t}_{m}");
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
  // A unit pulse: the height read is the gain g of the paper's Eq. (gain).
  Note(endUs + 1.3, 0.77,
       Form("#it{g} = #hat{#it{s}}_{1}(#it{t}_{m}) = %.2f", out->Eval(tm)), 12);
  TimeArrow(0.0, tf, -0.50, -0.44, "baseline frozen for #it{t}_{f}");
  Note(endUs + 1.5, -0.27,
       Form("tail:  #hat{#it{s}} = #gamma "
            "e^{#minus#kern[0.15]{#it{t}}/#tau_{pre}},  %s",
            TailLabel(tail).Data()),
       13, kRightColour);
  PlottingUtils::AddText("#tau_{pre} = 10 #mus,  #tau_{pz} = 20 #mus", 0.87,
                         0.84);
  PanelLabel("(a)");

  PlottingUtils::SaveFigure(canvas, "single", "", PlotSaveOptions::kLINEAR);
  delete canvas;
}

// (Fig. correction b) What the correction models, on a short train on the
// right chain (tau_pre = 10 us) with an ideal baseline: the height read on the
// last pulse is its height alone plus the tail each earlier pulse leaves at
// that read point, kappa(dt_i) h~_i with h~_i the earlier pulse's height alone
// (the paper's Eq. (kernel)).
void CorrectionTrainFigure() {
  const Double_t tauPre = 10.0;
  const Double_t tm = (kTrapRiseNs + kPeakingDelayNs) / 1000.0;
  const Int_t nPulses = 4, last = nPulses - 1;
  const Double_t arrival[nPulses] = {0.0, 8.0, 16.0, 24.0};
  const Double_t amplitude[nPulses] = {1.0, 1.5, 0.8, 1.0};
  const Double_t tRead = arrival[last] + tm;
  const Double_t tMin = -1.0, tMax = 44.0;
  const Double_t yLo = -0.7, yHi = 2.4;

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
    Note(arrival[i] + tm, h + 0.08, Form("#tilde{#it{h}}_{%d}", i + 1), 21);
    const Double_t y = 1.65 + 0.24 * (last - 1 - i);
    TimeArrow(arrival[i], arrival[last], y, y + 0.06,
              Form("#Delta#it{t}_{%d}", i + 1));
  }
  TimeMark(tRead, yLo, alone[last]->Eval(tRead));

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
       12);
  const Double_t xTailText = tRead + 4.5, yTailText = -0.5;
  TLine *leadTails = new TLine(tRead + 0.2, tails->GetY()[last - 1],
                               xTailText - 0.2, yTailText);
  leadTails->SetLineWidth(PlottingUtils::GetLineWidth());
  leadTails->Draw();
  Note(xTailText, yTailText, "tails of pulses 1-3", 12);

  std::vector<TString> labels;
  labels.push_back("Measured");
  labels.push_back("Pulses 1-3 alone");
  labels.push_back("Pulse 4 alone");
  TLegend *legend = Legend(labels, 0.66, 0.87);
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
  // The dashed spectra first, so the corrected ones are drawn over them.
  spectra[0][0]->SetMaximum(1.45 * yMax);
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
  TLegend *legend = Legend(labels, 0.61, 0.87);
  legend->AddEntry(spectra[0][1], labels[0], "l");
  legend->AddEntry(spectra[1][1], labels[1], "l");
  legend->AddEntry(dashed, labels[2], "l");
  legend->AddEntry(solid, labels[3], "l");
  legend->Draw();
  PlottingUtils::AddText(rateLabel, 0.88, 0.92);
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
  // Where the earlier pulse's trapezoid has ended and the analytic kernel
  // begins to hold: (k + l)T + t_r - t_m, t_r the simulation's detector ramp.
  TimeMark(dtTail, yLo, yHi);

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
    labels.push_back(Form("Fitted Kernel (%s)", name[j].Data()));
    entries.push_back(fitted[j]);
    options.push_back("p");
    if (form[j]) {
      labels.push_back(Form("Analytic Kernel (%s)", name[j].Data()));
      entries.push_back(form[j]);
      options.push_back("l");
    }
  }
  TLegend *legend = Legend(labels, 0.17, 0.17 + 0.06 * labels.size());
  for (size_t e = 0; e < labels.size(); e++)
    legend->AddEntry(entries[e], labels[e], options[e]);
  legend->Draw();
  PlottingUtils::AddText(runLabel, 0.88, 0.92);
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
  table << Form("%% For the caption: rms of the pulse heights within 0.5--1.5 "
                "times the beam peak, pileup left out, in percent of the "
                "chain's mean height; beam amplitudes with a %g %% Gaussian spread; after each "
                "trigger the baseline is frozen for %g us; at the end of each "
                "freeze it resumes if the trapezoid is within %g ADC of it and "
                "otherwise freezes again, the trapezoid taken as its running "
                "average over the last %d samples, kept up during the freeze "
                "(running-average check), or its last sample (single-sample "
                "check); tau_pre = %g us (matched), %g us (mismatched), tau_pz "
                "= %g us.",
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

// The amplitude-dependent kernels: the correction with one kernel per band of
// the previous pulse height, bands [0, 0.5), [0.5, 1.5), ..., [kNBands - 1.5,
// inf) times the beam peak, against the production single kernel, on the
// figures' run and on run 97.
const Int_t kNBands = 4;

// What scripts/run97_figure.sh measures for the band macros: the long groups'
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
    std::cout << "skipping band_macros: no band spreads in " << file
              << " (rerun scripts/run97_figure.sh)" << std::endl;
    return kFALSE;
  }
  for (Int_t c = 0; c < 2; c++)
    for (Int_t k = 0; k < 3; k++)
      data.spread[c][k] = (*chain[c])[k];
  data.bands = bands->GetVal();
  data.rate_kHz = rate->GetVal();
  if (data.bands != kNBands) {
    std::cout << "skipping band_macros: " << file << " has " << data.bands
              << " bands, the simulation " << kNBands << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

// The amplitude-dependent kernels of the paper's Section (baseline), as LaTeX
// macros for the text, in percent of the mean height: the mismatched chain's
// spread before correction, with one kernel and with kNBands, on the figures'
// run (the first reading at the higher rate), and the same for run 97's Long
// Strips, Right, pooled over the strips; the number of bins and both rates.
void BandMacros(const Spread &one, const Spread &many, Double_t rate,
                const DataBands &data, const TString &path) {
  std::ofstream macros(path.Data());
  macros << "% Written by mains/figures.cpp; edit that, not this file."
         << std::endl;
  macros << Form("%% Simulated: %s at %g kHz, mismatched chain. Run 97: Long "
                 "Strips, Right, each strip in units of its own beam peak, "
                 "corrected with the kernels fitted on it. Bins of the "
                 "previous pulse height, in beam peaks: below 0.5, ",
                 kSpreadReadings[0].label, rate);
  for (Int_t b = 1; b < kNBands - 1; b++)
    macros << Form("%g--%g, ", b - 0.5, b + 0.5);
  macros << Form("above %g.", kNBands - 1.5) << std::endl;
  macros << Form("\\newcommand{\\BandCount}{%d}", kNBands) << std::endl;
  macros << Form("\\newcommand{\\BandSimRate}{%g}", rate) << std::endl;
  macros << Form("\\newcommand{\\BandSimBefore}{%.2f}", one.RawPercent(1))
         << std::endl;
  macros << Form("\\newcommand{\\BandSimOne}{%.2f}", one.CorrPercent(1))
         << std::endl;
  macros << Form("\\newcommand{\\BandSimMany}{%.2f}", many.CorrPercent(1))
         << std::endl;
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
// strip's beam peak: the figure MUSIC's events summary draws, for the
// segmented strips 1-16 only, with the side legend sized to its entries. The
// paper shows the two as separate figures, so they carry no panel letter.
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
  const Double_t yLo = 0.2, yHi = 1.6;

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
  TLegend *legend = Legend(labels, 0.73, 0.87);
  legend->AddEntry(beamLine, labels[0], "l");
  legend->AddEntry(reactionLine, labels[1], "l");
  legend->Draw();
  PlottingUtils::AddText(corrected ? "After Correction" : "Before Correction",
                         0.88, 0.92);
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
  TString simDir = "output/sim";
  TString run97 = "output/run97/run97.root";
  TString reactionEvents = "output/reaction_events/reaction_events.root";
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

    // The figures' run corrected again with kNBands amplitude bands, beside
    // what run 97 shows.
    DataBands dataBands;
    Spread banded;
    RunInfo info;
    if (haveSpreads && HaveInput(run97, "band_macros") &&
        ReadDataBands(run97, dataBands) &&
        MeasureSpread(high, TString(gSystem->BaseName(high)) +
                                Form(" %d bands", kNBands),
                      banded, info, kNBands))
      BandMacros(spreads[shown][highRate], banded, highRate,
                 dataBands, out + "/band_macros.tex");
  }

  if (HaveInput(run97, "run97 figures")) {
    DeviationFigure(run97);
    StripSpectraFigure(run97, kFALSE);
    StripSpectraFigure(run97, kTRUE);
    if (format != PlotSaveFormat::kPNG)
      FwhmMacros(run97, out + "/run97_fwhm.tex");
  }
  if (HaveInput(reactionEvents, "reaction figures"))
    for (Int_t corrected = 0; corrected < 2; corrected++)
      ReactionTracesFigure(reactionEvents, corrected);
  if (format == PlotSaveFormat::kPNGAndPDF)
    MovePngs(out);
  return 0;
}
