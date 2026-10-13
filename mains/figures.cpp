// Figures for doc/polezero.tex.
//
// Five plots, written as PDF into doc/figs:
//   geometric_series   what the deconvolution stage does to a sampled
//                      exponential, and what the M d[n] term restores;
//   filter_tail        the normalised filter output for the four combinations
//                      of preamplifier and pole-zero setting used by MUSIC;
//   baseline_slope     the fitted linear coefficient against rate, with the
//                      prediction of the baseline-depletion argument;
//   profile_fit        an energy-vs-dt profile with the three- and
//                      four-parameter fits and their residuals;
//   tau_recovery       the decay time recovered from each strip group.
//
//   figures --out=<dir> --long40=<file> --long10=<file> --long1=<file>
//           --long27=<file> --short27=<file>
#include "Config.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include "Preamp.hpp"
#include "SimOutput.hpp"
#include "TrapezoidFilter.hpp"
#include <TArrow.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TProfile.h>
#include <iostream>
#include <vector>

namespace {

const Double_t kWindowLo = 0.5;
const Double_t kWindowHi = 1.5;

// Tail amplitude of the filter output, Eq. (15) of the note. Positive is an
// undershoot.
Double_t TailAmplitude(Double_t tauPre, Double_t tauPz, Double_t riseUs,
                       Double_t flatTopEndUs) {
  return (tauPz - tauPre) * tauPre / (riseUs * tauPz) *
         (1.0 - TMath::Exp(riseUs / tauPre)) *
         (1.0 - TMath::Exp(flatTopEndUs / tauPre));
}

TProfile *BuildProfile(const std::vector<SimHit> &hits, Int_t channel,
                       const TString &name, Double_t maxUs, Int_t bins,
                       Double_t &mode) {
  TH1D spec("fig_spec", "", 128, 0.0, 4096.0);
  spec.SetDirectory(nullptr);
  for (size_t i = 0; i < hits.size(); i++)
    if (hits[i].channel == channel)
      spec.Fill(hits[i].energy);
  spec.GetXaxis()->SetRange(9, 128);
  mode = spec.GetBinCenter(spec.GetMaximumBin());
  TProfile *p = new TProfile(name, "", bins, 0.0, maxUs);
  p->SetDirectory(nullptr);
  for (size_t i = 0; i < hits.size(); i++) {
    const SimHit &h = hits[i];
    if (h.channel != channel || h.pileup || h.dtPrev < 0.0)
      continue;
    if (h.energy <= kWindowLo * mode || h.energy >= kWindowHi * mode)
      continue;
    p->Fill(h.dtPrev, h.energy);
  }
  return p;
}

TF1 *FitProfile(TProfile *p, Double_t lo, Double_t hi, Bool_t withSlope,
                const TString &name) {
  const char *form = withSlope ? "[0] + [1] * exp(-x / [2]) + [3] * x"
                               : "[0] + [1] * exp(-x / [2])";
  TF1 *f = new TF1(name, form, lo, hi);
  const Double_t tail = p->GetBinContent(p->FindBin(hi - 10.0));
  const Double_t head = p->GetBinContent(p->FindBin(lo + 0.5));
  f->SetParameters(tail, head - tail, 12.0);
  if (withSlope)
    f->SetParameter(3, 0.0);
  f->SetParLimits(2, 1.0, 200.0);
  p->Fit(f, "RQN0");
  return f;
}

Double_t MeanAmplitude(const std::vector<SimHit> &hits, Int_t channel) {
  Double_t sum = 0.0;
  Long64_t n = 0;
  for (size_t i = 0; i < hits.size(); i++)
    if (hits[i].channel == channel) {
      sum += hits[i].trueAmp;
      n++;
    }
  return n > 0 ? sum / Double_t(n) : 0.0;
}

// ---------------------------------------------------------------- figure 1
void GeometricSeries() {
  const Int_t n = 26;
  const Double_t a = TMath::Exp(-0.2); // tau = 5 sampling periods
  const Double_t M = a / (1.0 - a);
  TGraph *samples = new TGraph(n);
  TGraph *partial = new TGraph(n);
  TGraph *restored = new TGraph(n);
  TGraph *total = new TGraph(n);
  Double_t running = 0.0;
  for (Int_t i = 0; i < n; i++) {
    const Double_t term = TMath::Power(a, Double_t(i));
    running += term;
    samples->SetPoint(i, Double_t(i), term);
    partial->SetPoint(i, Double_t(i), running);
    restored->SetPoint(i, Double_t(i), M * term);
    total->SetPoint(i, Double_t(i), running + M * term);
  }
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(
      partial, kBlue + 1,
      ";Sample index n;Value in units of the first sample");
  partial->GetYaxis()->SetRangeUser(0.0, 1.35 / (1.0 - a));
  partial->Draw("AL");
  PlottingUtils::ConfigureGraph(samples, kGray + 2);
  samples->SetLineStyle(2);
  samples->Draw("L SAME");
  PlottingUtils::ConfigureGraph(restored, kOrange + 7);
  restored->Draw("L SAME");
  PlottingUtils::ConfigureGraph(total, kRed + 1);
  total->Draw("L SAME");
  TLine *asymptote = new TLine(0.0, 1.0 / (1.0 - a), Double_t(n - 1),
                               1.0 / (1.0 - a));
  asymptote->SetLineStyle(3);
  asymptote->SetLineWidth(PlottingUtils::GetLineWidth());
  asymptote->Draw();
  // The shortfall at one sample, drawn as an arrow-free double tick.
  const Int_t mark = 8;
  Double_t x = 0.0, ylow = 0.0;
  partial->GetPoint(mark, x, ylow);
  TLine *gap = new TLine(x, ylow, x, 1.0 / (1.0 - a));
  gap->SetLineColor(kBlack);
  gap->SetLineWidth(2 * PlottingUtils::GetLineWidth());
  gap->Draw();
  TLegend *leg = PlottingUtils::AddLegend(0.42, 0.89, 0.20, 0.58);
  leg->AddEntry(samples, "a^{n}", "l");
  leg->AddEntry(partial, "G_{n} = #Sigma_{j#leq n} a^{j}", "l");
  leg->AddEntry(restored, "M_{pre} a^{n}", "l");
  leg->AddEntry(total, "G_{n} + M_{pre} a^{n}", "l");
  leg->AddEntry(asymptote, "1/(1#minus a)", "l");
  leg->AddEntry(gap, "area still to come", "l");
  leg->Draw();
  PlottingUtils::SaveFigure(c, "geometric_series", "",
                            PlotSaveOptions::kLINEAR);
  delete c;
}

// ---------------------------------------------------------------- figure 2
void FilterTail(const Config &cfg) {
  const Double_t tauPre[4] = {25.0, 10.0, 25.0, 10.0};
  const Double_t tauPz[4] = {20.0, 20.0, 10.0, 10.0};
  const char *label[4] = {"long L, 25 / 20 #mus", "long R, 10 / 20 #mus",
                          "short L, 25 / 10 #mus", "short R, 10 / 10 #mus"};
  const Int_t colour[4] = {kBlue + 1, kRed + 1, kAzure + 7, kGreen + 2};
  const Long64_t nSamples = Long64_t(45000.0 / cfg.samplePeriod_ns);
  const Long64_t start = Long64_t(2000.0 / cfg.samplePeriod_ns);
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  TLegend *leg = PlottingUtils::AddLegend(0.40, 0.89, 0.60, 0.87);
  TGraph *g[4];
  for (Int_t v = 0; v < 4; v++) {
    Preamp preamp(tauPre[v] * 1000.0, cfg.DetRise_ns(), cfg.samplePeriod_ns);
    TrapezoidFilter filter(cfg.RiseSamples(), cfg.FlatSamples(),
                           cfg.samplePeriod_ns, tauPz[v] * 1000.0);
    preamp.AddEvent(Double_t(start) * cfg.samplePeriod_ns, 1.0);
    g[v] = new TGraph();
    Int_t point = 0;
    for (Long64_t i = 0; i < nSamples; i++) {
      const Double_t out = filter.Process(preamp.Sample(i));
      const Double_t t =
          (Double_t(i - start) * cfg.samplePeriod_ns) / 1000.0;
      if (t > -2.0)
        g[v]->SetPoint(point++, t, out);
    }
    PlottingUtils::ConfigureGraph(
        g[v], colour[v],
        ";t #minus trigger [#mus];Filter output #hat{s} / A");
    if (v == 0) {
      g[v]->GetYaxis()->SetRangeUser(-0.45, 0.45);
      g[v]->Draw("AL");
    } else {
      g[v]->Draw("L SAME");
    }
    leg->AddEntry(g[v],
                  Form("%s: c = %+.2f", label[v],
                       TailAmplitude(tauPre[v], tauPz[v],
                                     cfg.trapRise_ns / 1000.0,
                                     (cfg.trapRise_ns + cfg.trapFlat_ns) /
                                         1000.0)),
                  "l");
  }
  leg->Draw();
  TLine *zero = new TLine(-2.0, 0.0, 43.0, 0.0);
  zero->SetLineStyle(3);
  zero->SetLineWidth(PlottingUtils::GetLineWidth());
  zero->Draw();
  PlottingUtils::SaveFigure(c, "filter_tail", "", PlotSaveOptions::kLINEAR);
  delete c;
}

// ---------------------------------------------------------------- figure 2b
// Filter output for unit pulses on the long-R chain against time from the
// first arrival. A one-sample detector ramp stands in for the step input of
// the derivation, so the trapezoid corners sit exactly at k, l and k + l.
TGraph *StepResponse(const Config &cfg, Double_t tauPreUs, Double_t tauPzUs,
                     const std::vector<Double_t> &arrivalsUs, Double_t tMinUs,
                     Double_t tMaxUs, Int_t stride = 1) {
  const Long64_t pad = Long64_t(2000.0 / cfg.samplePeriod_ns);
  const Long64_t nSamples =
      pad + Long64_t(tMaxUs * 1000.0 / cfg.samplePeriod_ns) + 1;
  Preamp preamp(tauPreUs * 1000.0, cfg.samplePeriod_ns, cfg.samplePeriod_ns);
  TrapezoidFilter filter(cfg.RiseSamples(), cfg.FlatSamples(),
                         cfg.samplePeriod_ns, tauPzUs * 1000.0);
  for (size_t a = 0; a < arrivalsUs.size(); a++)
    preamp.AddEvent(
        Double_t(pad) * cfg.samplePeriod_ns + arrivalsUs[a] * 1000.0, 1.0);
  TGraph *g = new TGraph();
  Int_t point = 0;
  for (Long64_t i = 0; i < nSamples; i++) {
    const Double_t out = filter.Process(preamp.Sample(i));
    const Double_t t = Double_t(i - pad) * cfg.samplePeriod_ns / 1000.0;
    // Every stride-th point: a dashed line drawn through thousands of
    // 8 ns segments comes out solid in ROOT's raster output.
    if (t >= tMinUs && (i - pad) % stride == 0)
      g->SetPoint(point++, t, out);
  }
  return g;
}

// A double-headed arrow between two times with its label above.
void TimeArrow(Double_t x1, Double_t x2, Double_t y, const char *text) {
  TArrow *arrow = new TArrow(x1, y, x2, y, 0.012, "<|>");
  arrow->SetLineWidth(PlottingUtils::GetLineWidth());
  arrow->SetFillColor(kBlack);
  arrow->Draw();
  TLatex *label = new TLatex(0.5 * (x1 + x2), y + 0.04, text);
  label->SetTextFont(43);
  label->SetTextSize(30);
  label->SetTextAlign(21);
  label->Draw();
}

// A dashed vertical at one time.
void TimeMark(Double_t x, Double_t yLo, Double_t yHi) {
  TLine *line = new TLine(x, yLo, x, yHi);
  line->SetLineStyle(2);
  line->SetLineColor(kGray + 2);
  line->SetLineWidth(PlottingUtils::GetLineWidth());
  line->Draw();
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

void PileupTiming(const Config &cfg) {
  const Double_t tauPre = 10.0, tauPz = 20.0;
  const Double_t kUs = cfg.trapRise_ns / 1000.0;
  const Double_t lUs = (cfg.trapRise_ns + cfg.trapFlat_ns) / 1000.0;
  const Double_t tp = cfg.peakingDelay_ns / 1000.0;
  const Double_t tm = kUs + tp;
  const Double_t tf = cfg.baselineHoldoff_ns / 1000.0;
  const Double_t dt = 6.0;
  const Double_t tMin = -1.0, tMax = 20.0;
  const Double_t yLo = -0.6, yHi = 2.0;
  const Double_t c = TailAmplitude(tauPre, tauPz, kUs, lUs);

  // ---- (a) one pulse on its own
  std::vector<Double_t> one(1, 0.0);
  TGraph *single = StepResponse(cfg, tauPre, tauPz, one, tMin, tMax, 8);
  TCanvas *ca = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(
      single, kRed + 1,
      ";t #minus arrival of the pulse [#mus];Filter output #hat{s} / A");
  single->GetXaxis()->SetLimits(tMin, tMax);
  single->GetYaxis()->SetRangeUser(yLo, yHi);
  single->Draw("AL");
  TimeMark(0.0, yLo, 1.15);
  TimeMark(kUs, yLo, 1.15);
  TimeMark(lUs, yLo, 1.15);
  TimeMark(kUs + lUs, yLo, 1.15);
  Note(0.5 * kUs, 1.17, "rise", 21);
  Note(0.5 * (kUs + lUs), 1.17, "flat top", 21);
  Note(0.5 * (lUs + kUs + lUs), 1.17, "fall", 21);
  Note(kUs + lUs + 3.0, 1.17, "tail", 21);
  TimeArrow(0.0, kUs, 1.42, "kT_{clk}");
  TimeArrow(kUs, tm, 1.42, "t_{p}");
  TimeArrow(0.0, tm, 1.72, "t_{m} = kT_{clk} + t_{p}");
  TimeMark(tm, single->Eval(tm), 1.42);
  TGraph *read1 = new TGraph(1);
  read1->SetPoint(0, tm, single->Eval(tm));
  read1->SetMarkerStyle(20);
  read1->SetMarkerSize(2.0);
  read1->SetMarkerColor(kBlack);
  read1->Draw("P SAME");
  TLine *leadRead = new TLine(tm + 0.15, single->Eval(tm) + 0.03, 8.3, 1.0);
  leadRead->SetLineWidth(PlottingUtils::GetLineWidth());
  leadRead->Draw();
  Note(8.5, 1.0, "energy read at t_{m}", 12);
  TimeArrow(0.0, tf, -0.50, "baseline frozen, t_{f}");
  Note(kUs + lUs + 2.5, -0.30, "tail  #minus c A e^{#minus t/#tau_{pre}}", 12,
       kRed + 1);
  PlottingUtils::SaveFigure(ca, "pileup_single", "",
                            PlotSaveOptions::kLINEAR);
  delete ca;

  // ---- (b) a second pulse on the tail of the first
  std::vector<Double_t> two(2);
  two[0] = 0.0;
  two[1] = dt;
  std::vector<Double_t> second(1, dt);
  TGraph *pair = StepResponse(cfg, tauPre, tauPz, two, tMin, tMax, 8);
  TGraph *first = StepResponse(cfg, tauPre, tauPz, one, tMin, tMax, 8);
  TGraph *alone = StepResponse(cfg, tauPre, tauPz, second, tMin, tMax, 8);
  TCanvas *cb = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(
      pair, kRed + 1,
      ";t #minus arrival of pulse 1 [#mus];Filter output #hat{s} / A");
  pair->GetXaxis()->SetLimits(tMin, tMax);
  pair->GetYaxis()->SetRangeUser(yLo, yHi);
  pair->Draw("AL");
  PlottingUtils::ConfigureGraph(first, kGray + 2);
  first->SetLineStyle(2);
  first->Draw("L SAME");
  PlottingUtils::ConfigureGraph(alone, kBlue + 1);
  alone->SetLineStyle(2);
  alone->Draw("L SAME");
  TimeMark(0.0, yLo, 1.15);
  TimeMark(dt, yLo, 1.15);
  TimeMark(dt + tm, yLo, 1.15);
  Note(0.0, 1.17, "pulse 1", 21);
  Note(dt, 1.17, "pulse 2", 21);
  Note(dt + tm, 1.17, "E_{2} read", 21);
  TimeArrow(0.0, dt, 1.42, "#Deltat");
  TimeArrow(dt, dt + tm, 1.42, "t_{m}");
  TimeArrow(0.0, dt + tm, 1.72, "#Deltat + t_{m}");
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
  // Leader lines out to the free space on the right for the two bars.
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
  TLegend *leg = PlottingUtils::AddLegend(0.62, 0.89, 0.62, 0.87);
  leg->AddEntry(pair, "both pulses", "l");
  leg->AddEntry(first, "pulse 1 alone", "l");
  leg->AddEntry(alone, "pulse 2 alone", "l");
  leg->Draw();
  PlottingUtils::SaveFigure(cb, "pileup_pair", "", PlotSaveOptions::kLINEAR);
  delete cb;
  std::cout << Form("pileup timing: t_m = %.1f us, dt = %.1f us, loss %.4f, "
                    "tail of pulse 1 %.4f, c e^(-(dt+t_m)/tau) = %.4f",
                    tm, dt, yAlone - yPair, yTail,
                    c * TMath::Exp(-(dt + tm) / tauPre))
            << std::endl;
}

// ---------------------------------------------------------------- figure 3
void BaselineSlope(const std::vector<TString> &files, Double_t windowUs,
                   Double_t riseUs, Double_t flatEndUs, Double_t freezeUs) {
  const Int_t nch = 2;
  const Int_t colour[nch] = {kBlue + 1, kRed + 1};
  const char *chName[nch] = {"L", "R"};
  TGraphErrors *measured[nch];
  TGraph *predicted[nch];
  for (Int_t ch = 0; ch < nch; ch++) {
    measured[ch] = new TGraphErrors();
    predicted[ch] = new TGraph();
  }
  for (size_t f = 0; f < files.size(); f++) {
    SimOutput::RunHeader header;
    std::vector<SimHit> hits;
    if (!SimOutput::Read(files[f], header, hits))
      continue;
    for (Int_t ch = 0; ch < nch; ch++) {
      Double_t mode = 0.0;
      TProfile *p = BuildProfile(hits, ch, Form("bs_%zu_%d", f, ch), 150.0,
                                 150, mode);
      TF1 *fit = FitProfile(p, freezeUs + 1.0, 130.0, kTRUE,
                            Form("bsfit_%zu_%d", f, ch));
      const Double_t tauPre = ch == 0 ? header.tauL_us : header.tauR_us;
      const Double_t c =
          TailAmplitude(tauPre, header.trapDecay_us, riseUs, flatEndUs);
      const Double_t amp = MeanAmplitude(hits, ch);
      // p3 = - R c A tau exp(-t_freeze / tau) / W, with R per microsecond.
      // The exponential is there because the baseline is frozen for t_freeze
      // after every trigger, so the averager never sees the first and largest
      // part of each tail.
      const Double_t pred = -(header.rate_kHz * 1.0e-3) * c * amp * tauPre *
                            TMath::Exp(-freezeUs / tauPre) / windowUs;
      const Int_t np = measured[ch]->GetN();
      measured[ch]->SetPoint(np, header.rate_kHz, fit->GetParameter(3));
      measured[ch]->SetPointError(np, 0.0, fit->GetParError(3));
      predicted[ch]->SetPoint(np, header.rate_kHz, pred);
      std::cout << Form("  %s  rate %5.1f kHz   p3 fitted %+7.3f +- %.3f   "
                        "predicted %+7.3f  (no freeze factor %+7.3f)",
                        chName[ch], header.rate_kHz, fit->GetParameter(3),
                        fit->GetParError(3), pred,
                        pred / TMath::Exp(-freezeUs / tauPre))
                << std::endl;
      // Diagnostics: the fitted decay time, and p3 when the decay time is
      // held at the true value, to show how much of p3 is fit degeneracy.
      TF1 *fixedTau = new TF1(Form("bsfix_%zu_%d", f, ch),
                              "[0] + [1] * exp(-x / [2]) + [3] * x",
                              freezeUs + 1.0, 130.0);
      fixedTau->SetParameters(fit->GetParameter(0), fit->GetParameter(1),
                              tauPre, fit->GetParameter(3));
      fixedTau->FixParameter(2, tauPre);
      p->Fit(fixedTau, "RQN0");
      std::cout << Form("      p1 %+8.2f  p2 %5.1f +- %4.1f us  |  p2 fixed at "
                        "%4.1f: p3 %+7.3f +- %.3f   chi2/ndf %.2f -> %.2f",
                        fit->GetParameter(1), fit->GetParameter(2),
                        fit->GetParError(2), tauPre,
                        fixedTau->GetParameter(3), fixedTau->GetParError(3),
                        fit->GetChisquare() / fit->GetNDF(),
                        fixedTau->GetChisquare() / fixedTau->GetNDF())
                << std::endl;
      delete p;
    }
  }
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  TH1F *frame = c->DrawFrame(0.0, -0.62, 44.0, 0.62);
  frame->SetTitle(";Beam rate [kHz];p_{3} [ADC / #mus]");
  frame->GetXaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetXaxis()->SetLabelSize(0.06);
  frame->GetYaxis()->SetLabelSize(0.06);
  frame->GetYaxis()->SetTitleOffset(1.2);
  TLegend *leg = PlottingUtils::AddLegend(0.14, 0.58, 0.60, 0.87);
  for (Int_t ch = 0; ch < nch; ch++) {
    PlottingUtils::ConfigureGraph(predicted[ch], colour[ch]);
    predicted[ch]->SetLineStyle(2);
    predicted[ch]->Draw("L SAME");
    PlottingUtils::ConfigureGraph(measured[ch], colour[ch]);
    measured[ch]->Draw("P SAME");
    leg->AddEntry(measured[ch], Form("%s fitted", chName[ch]), "p");
    leg->AddEntry(predicted[ch], Form("%s predicted", chName[ch]), "l");
  }
  leg->Draw();
  TLine *zero = new TLine(0.0, 0.0, 44.0, 0.0);
  zero->SetLineStyle(3);
  zero->SetLineWidth(PlottingUtils::GetLineWidth());
  zero->Draw();
  PlottingUtils::SaveFigure(c, "baseline_slope", "",
                            PlotSaveOptions::kLINEAR);
  delete c;
}

// ---------------------------------------------------------------- figure 4
void ProfileFit(const TString &file) {
  SimOutput::RunHeader header;
  std::vector<SimHit> hits;
  if (!SimOutput::Read(file, header, hits))
    return;
  const Double_t freeze =
      (2.0 * header.trapRise_ns + header.trapFlat_ns) / 1000.0;
  Double_t mode = 0.0;
  TProfile *p = BuildProfile(hits, 1, "pf_R", 150.0, 150, mode);
  TF1 *f3 = FitProfile(p, freeze + 1.0, 130.0, kFALSE, "pf_three");
  TF1 *f4 = FitProfile(p, freeze + 1.0, 130.0, kTRUE, "pf_four");
  std::cout << Form("profile fit: three-parameter p2 = %.2f us, "
                    "four-parameter p2 = %.2f us, true %.1f us",
                    f3->GetParameter(2), f4->GetParameter(2), header.tauR_us)
            << std::endl;

  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  TPad *top = new TPad("pf_top", "", 0.0, 0.30, 1.0, 1.0);
  TPad *bot = new TPad("pf_bot", "", 0.0, 0.0, 1.0, 0.30);
  top->SetBottomMargin(0.03);
  top->SetLeftMargin(0.15);
  bot->SetTopMargin(0.03);
  bot->SetBottomMargin(0.32);
  bot->SetLeftMargin(0.15);
  top->Draw();
  bot->Draw();
  top->cd();
  PlottingUtils::ConfigureHistogram(p, kBlack, ";;Mean pulse height [ADC]");
  p->SetMarkerStyle(20);
  p->SetMarkerSize(0.7);
  p->GetXaxis()->SetLabelSize(0.0);
  p->GetYaxis()->SetRangeUser(1180.0, 1700.0);
  p->GetYaxis()->SetTitleSize(0.075);
  p->GetYaxis()->SetLabelSize(0.07);
  p->GetYaxis()->SetTitleOffset(0.95);
  p->Draw("E1");
  f3->SetLineColor(kOrange + 7);
  f3->SetLineWidth(PlottingUtils::GetLineWidth());
  f3->SetLineStyle(2);
  f3->Draw("SAME");
  f4->SetLineColor(kBlue + 1);
  f4->SetLineWidth(PlottingUtils::GetLineWidth());
  f4->Draw("SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.40, 0.89, 0.12, 0.46);
  leg->AddEntry(p, Form("simulation, #tau_{pre} = %.0f #mus", header.tauR_us),
                "lp");
  leg->AddEntry(f3, Form("three parameters: p_{2} = %.2f #mus",
                         f3->GetParameter(2)),
                "l");
  leg->AddEntry(f4, Form("with p_{3}#Deltat: p_{2} = %.2f #mus",
                         f4->GetParameter(2)),
                "l");
  leg->Draw();
  bot->cd();
  TGraph *res3 = new TGraph();
  TGraph *res4 = new TGraph();
  Int_t n3 = 0, n4 = 0;
  for (Int_t b = 1; b <= p->GetNbinsX(); b++) {
    const Double_t x = p->GetBinCenter(b);
    if (x < freeze + 1.0 || x > 130.0 || p->GetBinEntries(b) < 5)
      continue;
    res3->SetPoint(n3++, x, p->GetBinContent(b) - f3->Eval(x));
    res4->SetPoint(n4++, x, p->GetBinContent(b) - f4->Eval(x));
  }
  TH1F *rframe = bot->DrawFrame(0.0, -13.0, 150.0, 13.0);
  rframe->SetTitle(";#Deltat to previous hit [#mus];data #minus fit");
  rframe->GetXaxis()->SetTitleSize(0.15);
  rframe->GetYaxis()->SetTitleSize(0.15);
  rframe->GetXaxis()->SetLabelSize(0.13);
  rframe->GetYaxis()->SetLabelSize(0.13);
  rframe->GetXaxis()->SetTitleOffset(0.95);
  rframe->GetYaxis()->SetTitleOffset(0.47);
  rframe->GetYaxis()->SetNdivisions(505);
  PlottingUtils::ConfigureGraph(res3, kOrange + 7);
  res3->SetLineStyle(2);
  res3->Draw("L SAME");
  PlottingUtils::ConfigureGraph(res4, kBlue + 1);
  res4->Draw("L SAME");
  TLine *zero = new TLine(0.0, 0.0, 150.0, 0.0);
  zero->SetLineStyle(3);
  zero->SetLineWidth(PlottingUtils::GetLineWidth());
  zero->Draw();
  PlottingUtils::SaveFigure(c, "profile_fit", "", PlotSaveOptions::kLINEAR);
  delete c;
}

// ---------------------------------------------------------------- figure 5
void TauRecovery(const TString &longFile, const TString &shortFile) {
  const TString files[2] = {longFile, shortFile};
  const char *names[4] = {"long L", "long R", "short L", "short R"};
  TGraphErrors *three = new TGraphErrors();
  TGraphErrors *four = new TGraphErrors();
  TGraph *truth = new TGraph();
  Int_t slot = 0;
  for (Int_t f = 0; f < 2; f++) {
    SimOutput::RunHeader header;
    std::vector<SimHit> hits;
    if (!SimOutput::Read(files[f], header, hits))
      continue;
    const Double_t freeze =
        (2.0 * header.trapRise_ns + header.trapFlat_ns) / 1000.0;
    for (Int_t ch = 0; ch < 2; ch++) {
      const Double_t tauPre = ch == 0 ? header.tauL_us : header.tauR_us;
      truth->SetPoint(truth->GetN(), Double_t(slot), tauPre);
      if (TMath::Abs(tauPre - header.trapDecay_us) < 0.01) {
        std::cout << Form("  %-8s matched pole-zero, no tail to fit",
                          names[slot])
                  << std::endl;
        slot++;
        continue;
      }
      Double_t mode = 0.0;
      TProfile *p = BuildProfile(hits, ch, Form("tr_%d_%d", f, ch), 150.0, 150,
                                 mode);
      TF1 *f3 = FitProfile(p, freeze + 1.0, 130.0, kFALSE,
                           Form("trf3_%d_%d", f, ch));
      TF1 *f4 = FitProfile(p, freeze + 1.0, 130.0, kTRUE,
                           Form("trf4_%d_%d", f, ch));
      three->SetPoint(three->GetN(), Double_t(slot) - 0.16,
                      f3->GetParameter(2));
      three->SetPointError(three->GetN() - 1, 0.0, f3->GetParError(2));
      four->SetPoint(four->GetN(), Double_t(slot) + 0.16, f4->GetParameter(2));
      four->SetPointError(four->GetN() - 1, 0.0, f4->GetParError(2));
      std::cout << Form("  %-8s input %4.1f   three %5.2f +- %.2f   "
                        "four %5.2f +- %.2f  us",
                        names[slot], tauPre, f3->GetParameter(2),
                        f3->GetParError(2), f4->GetParameter(2),
                        f4->GetParError(2))
                << std::endl;
      slot++;
      delete p;
    }
  }
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  TH1F *frame = new TH1F("tr_frame", "", 4, -0.5, 3.5);
  frame->SetDirectory(nullptr);
  for (Int_t i = 0; i < 4; i++)
    frame->GetXaxis()->SetBinLabel(i + 1, names[i]);
  frame->SetMinimum(0.0);
  frame->SetMaximum(46.0);
  frame->SetTitle(";;Recovered #tau_{pre} [#mus]");
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetLabelSize(0.06);
  frame->GetYaxis()->SetTitleOffset(1.2);
  frame->GetXaxis()->SetLabelSize(0.075);
  frame->Draw("AXIS");
  truth->SetMarkerStyle(29);
  truth->SetMarkerSize(3.2);
  truth->SetMarkerColor(kGreen + 2);
  truth->Draw("P SAME");
  PlottingUtils::ConfigureGraph(three, kOrange + 7);
  three->SetMarkerStyle(21);
  three->Draw("P SAME");
  PlottingUtils::ConfigureGraph(four, kBlue + 1);
  four->Draw("P SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.16, 0.60, 0.60, 0.84);
  leg->AddEntry(truth, "input #tau_{pre}", "p");
  leg->AddEntry(three, "three parameters", "p");
  leg->AddEntry(four, "with p_{3}#Deltat", "p");
  leg->Draw();
  PlottingUtils::SaveFigure(c, "tau_recovery", "", PlotSaveOptions::kLINEAR);
  delete c;
}

} // namespace

int main(int argc, char **argv) {
  TString out = "doc/figs";
  TString long40, long10, long1, long27, short27;
  PlotSaveFormat format = PlotSaveFormat::kPDF;
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a.BeginsWith("--out="))
      out = a(6, a.Length() - 6);
    else if (a == "--format=png")
      format = PlotSaveFormat::kPNG;
    else if (a == "--format=pdf")
      format = PlotSaveFormat::kPDF;
    else if (a.BeginsWith("--long40="))
      long40 = a(9, a.Length() - 9);
    else if (a.BeginsWith("--long10="))
      long10 = a(9, a.Length() - 9);
    else if (a.BeginsWith("--long1="))
      long1 = a(8, a.Length() - 8);
    else if (a.BeginsWith("--long27="))
      long27 = a(9, a.Length() - 9);
    else if (a.BeginsWith("--short27="))
      short27 = a(10, a.Length() - 10);
  }
  InitUtils::SetROOTPreferences(format, out, ".", kFALSE);
  Config cfg;
  GeometricSeries();
  FilterTail(cfg);
  PileupTiming(cfg);
  std::cout << "baseline slope, fitted against predicted:" << std::endl;
  std::vector<TString> rates;
  rates.push_back(long40);
  rates.push_back(long10);
  rates.push_back(long1);
  BaselineSlope(rates,
                Double_t(cfg.baselineMeanSamples) * cfg.samplePeriod_ns /
                    1000.0,
                cfg.trapRise_ns / 1000.0,
                (cfg.trapRise_ns + cfg.trapFlat_ns) / 1000.0,
                cfg.baselineHoldoff_ns / 1000.0);
  ProfileFit(long40);
  TauRecovery(long27, short27);
  return 0;
}
