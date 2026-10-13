#include "Analysis.hpp"
#include "IOUtils.hpp"
#include "PlottingUtils.hpp"
#include "PulseHistoryAdapter.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TList.h>
#include <TLine.h>
#include <TMath.h>
#include <TProfile.h>
#include <TTree.h>
#include <iomanip>
#include <iostream>

namespace Analysis {

namespace {

const Int_t kColorL = kBlue + 1;
const Int_t kColorR = kRed + 1;
const Int_t kColorCorr = kGreen + 2;
const Double_t kProfileMaxUs = 60.0;
const Int_t kProfileBins = 120;
// The figures keep the 0-60 us range the study asks for, but a 25 us decay
// needs a longer lever arm than that, so the exponential is fitted on a wider
// profile of the same hits.
const Double_t kFitMaxUs = 150.0;
const Int_t kFitBins = 150;
// Beam window for means and RMS, in units of the beam peak: wide enough for
// the pole-zero distortions, narrow enough to exclude the merged double
// pulses at twice the peak (the same idea as PulseHistory's own-chain window
// of 0.5-1.6 x mode).
const Double_t kWindowLo = 0.5;
const Double_t kWindowHi = 1.5;
// Legend text in pixels (font 43); PlottingUtils uses 30, too large for the
// entries here.
const Double_t kLegendTextSize = 22.0;

const char *ChannelName(Int_t c) { return c == 0 ? "L" : "R"; }
Int_t ChannelColor(Int_t c) { return c == 0 ? kColorL : kColorR; }

// The baseline is frozen for the holdoff after each trigger. A window short
// enough to refill inside the plotted range then mixes pre- and post-pulse
// samples until it holds only the latter, and only from there is the profile
// a single exponential in tau, so the fit waits for that. A window longer than
// the range never refills: the restorer is static on this time scale and the
// fit can start as soon as the freeze ends.
Double_t FitLow_us(const SimOutput::RunHeader &h) {
  const Double_t freeze_us = h.baselineHoldoff_ns / 1000.0;
  const Double_t window_us =
      Double_t(h.baselineMeanSamples) * h.samplePeriod_ns / 1000.0;
  if (window_us > kProfileMaxUs)
    return freeze_us + 1.0;
  return freeze_us + window_us + 1.0;
}

void MeanRms(const std::vector<Double_t> &values, Double_t &mean,
             Double_t &rms) {
  mean = 0.0;
  rms = 0.0;
  if (values.empty())
    return;
  Double_t sum = 0.0, sum2 = 0.0;
  for (size_t i = 0; i < values.size(); i++) {
    sum += values[i];
    sum2 += values[i] * values[i];
  }
  mean = sum / Double_t(values.size());
  const Double_t var = sum2 / Double_t(values.size()) - mean * mean;
  rms = var > 0.0 ? TMath::Sqrt(var) : 0.0;
}

// A profile whose exponential amplitude is within 2 sigma of zero, or whose
// time constant ran into a limit, is flat: its time constant means nothing.
Bool_t IsFlat(const ExpFit &e) {
  return !e.ok || TMath::Abs(e.p1) < 2.0 * e.p1err || e.p2 <= 0.51 ||
         e.p2 >= 199.0;
}

TString TauLabel(const ExpFit &e) {
  if (IsFlat(e))
    return "flat";
  return Form("%.1f+-%.1f", e.p2, e.p2err);
}

// p0 + p1 exp(-dt/p2) + p3 dt. The linear term is not cosmetic: it stands for
// the trapezoid baseline, an average over a few hundred microseconds, which a
// hit with a long gap behind it finds less depressed by pile-up, so the
// profile is lifted roughly in proportion to the gap. Without it the fitted
// time constant comes out about 15 percent short on every group that has a
// pole-zero mismatch (bin/taustudy).
ExpFit FitExponential(TProfile *prof, Double_t lo, Double_t hi) {
  ExpFit r;
  TF1 *f = new TF1(PlottingUtils::GetRandomName(),
                   "[0] + [1] * exp(-x / [2]) + [3] * x", lo, hi);
  const Double_t tail = prof->GetBinContent(prof->FindBin(hi - 5.0));
  const Double_t head = prof->GetBinContent(prof->FindBin(lo + 0.5));
  f->SetParameters(tail, head - tail, 10.0, 0.0);
  f->SetParLimits(2, 0.5, 200.0);
  f->SetLineColor(kBlack);
  f->SetLineWidth(2);
  TFitResultPtr res = prof->Fit(f, "RQS");
  if (Int_t(res) == 0) {
    r.ok = kTRUE;
    r.p0 = f->GetParameter(0);
    r.p1 = f->GetParameter(1);
    r.p1err = f->GetParError(1);
    r.p2 = f->GetParameter(2);
    r.p2err = f->GetParError(2);
    r.p3 = f->GetParameter(3);
    r.chi2 = f->GetChisquare();
    r.ndf = f->GetNDF();
  }
  return r;
}

LineFit FitLine(TProfile *prof, Double_t lo, Double_t hi,
                Int_t minEntries = 1) {
  LineFit r;
  // Bins with fewer than minEntries hits carry no usable error estimate (a
  // single-entry bin has error zero and would dominate a chi2 fit), so they
  // are emptied before the fit; ROOT skips bins with zero error.
  for (Int_t b = 1; b <= prof->GetNbinsX(); b++) {
    if (prof->GetBinEntries(b) < minEntries) {
      prof->SetBinContent(b, 0.0);
      prof->SetBinEntries(b, 0.0);
      prof->SetBinError(b, 0.0);
    }
  }
  TF1 *f = new TF1(PlottingUtils::GetRandomName(), "[0] + [1] * x", lo, hi);
  f->SetParameters(prof->GetBinContent(prof->FindBin(0.5 * (lo + hi))), 0.0);
  f->SetLineColor(kBlack);
  f->SetLineStyle(2);
  f->SetLineWidth(2);
  TFitResultPtr res = prof->Fit(f, "RQS+");
  if (Int_t(res) == 0) {
    r.ok = kTRUE;
    r.intercept = f->GetParameter(0);
    r.interceptErr = f->GetParError(0);
    r.slope = f->GetParameter(1);
    r.slopeErr = f->GetParError(1);
    r.chi2 = f->GetChisquare();
    r.ndf = f->GetNDF();
  }
  return r;
}

TH1D *MakeSpectrum(const TString &name, Double_t beamPeak) {
  TH1D *h = new TH1D(name, ";Pulse height [ADC];Hits", 300, 0.0,
                     2.5 * beamPeak);
  h->SetDirectory(nullptr);
  return h;
}

TProfile *MakeProfile(const TString &name) {
  TProfile *p = new TProfile(name, ";#Deltat to previous hit [#mus];Mean pulse height [ADC]",
                             kProfileBins, 0.0, kProfileMaxUs);
  p->SetDirectory(nullptr);
  return p;
}

TProfile *MakeFitProfile(const TString &name) {
  TProfile *p = new TProfile(name, "", kFitBins, 0.0, kFitMaxUs);
  p->SetDirectory(nullptr);
  return p;
}

void StyleProfile(TProfile *p, Int_t color) {
  PlottingUtils::ConfigureHistogram(p, color);
  p->SetMarkerColor(color);
  p->SetMarkerStyle(20);
  p->SetMarkerSize(0.8);
}

void DrawSpectra(TH1D *hL, TH1D *hR, const ChannelSummary *cs,
                 Bool_t corrected, const TString &label,
                 const TString &name) {
  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  const Double_t ymax =
      1.45 * TMath::Max(hL->GetMaximum(), hR->GetMaximum());
  hL->SetMaximum(ymax);
  PlottingUtils::ConfigureAndDrawHistogram(hL, kColorL);
  PlottingUtils::ConfigureHistogram(hR, kColorR);
  hR->Draw("HIST SAME");
  TLegend *leg = PlottingUtils::AddLegend(0.58, 0.89, 0.74, 0.87);
  leg->SetTextSize(kLegendTextSize);
  for (Int_t ch = 0; ch < 2; ch++) {
    const Double_t mean = corrected ? cs[ch].meanCorr : cs[ch].meanRaw;
    const Double_t rms = corrected ? cs[ch].rmsCorr : cs[ch].rmsRaw;
    leg->AddEntry(ch == 0 ? hL : hR,
                  Form("%s: mean %.0f, RMS %.1f", ChannelName(ch), mean, rms),
                  "l");
  }
  leg->Draw();
  PlottingUtils::AddText(label, 0.88, 0.92);
  PlottingUtils::SaveFigure(c, name, label, PlotSaveOptions::kLINEAR);
  delete c;
}

} // namespace

RunSummary AnalyzeRun(const TString &path, const TString &label,
                      PulseHistory::Result *externalKernel,
                      const TString &externalLabel,
                      PulseHistory::Result &ownKernel) {
  RunSummary summary;
  summary.label = label;
  std::vector<SimHit> hits;
  if (!SimOutput::Read(path, summary.header, hits)) {
    std::cerr << "Skipping " << label << std::endl;
    return summary;
  }
  const SimOutput::RunHeader &hdr = summary.header;
  const Double_t peak = hdr.beamPeakADC;
  const Double_t fitLo = FitLow_us(hdr);
  std::cout << std::endl;
  std::cout << "=== " << label << ": " << path << " (" << hits.size()
            << " hits, " << hdr.rate_kHz << " kHz, tauL " << hdr.tauL_us
            << " us, tauR " << hdr.tauR_us << " us, pole-zero "
            << hdr.trapDecay_us << " us; baseline " << hdr.baselineMeanSamples
            << " samples, exponential fits over " << fitLo << "-"
            << kFitMaxUs - 20.0 << " us)" << std::endl;

  // --- PulseHistory correction, on the unmodified MUSIC code -------------
  PulseHistoryAdapter::Adapted adapted = PulseHistoryAdapter::Build(hits);
  std::vector<Int_t> groups = PulseHistory::BuildGroupMap();
  PulseHistory::Result *kernel = externalKernel;
  if (!kernel) {
    if (PulseHistory::Measure(adapted.hits, groups, ownKernel, label)) {
      summary.kernelMeasured = kTRUE;
      summary.kernelSource = label;
      kernel = &ownKernel;
    } else {
      std::cout << "  PulseHistory::Measure failed on " << label << std::endl;
    }
  } else {
    summary.kernelSource = externalLabel;
  }
  std::vector<Double_t> corrected(hits.size(), 0.0);
  for (size_t i = 0; i < hits.size(); i++)
    corrected[i] = hits[i].energy;
  if (kernel) {
    kernel->n_corrected = 0;
    kernel->n_clamped = 0;
    PulseHistory::Apply(adapted.hits, groups, *kernel);
    corrected = PulseHistoryAdapter::CorrectedEnergies(adapted, hits);
    if (!externalKernel) {
      std::cout << PulseHistory::Report(*kernel, label);
      PulseHistory::SavePlots(*kernel, label);
    } else {
      std::cout << "  applied the " << externalLabel << " kernel: "
                << kernel->n_corrected << " hits corrected, "
                << kernel->n_clamped << " clamped" << std::endl;
    }
    const Int_t groupOf[2] = {PulseHistory::kLongLeft,
                              PulseHistory::kLongRight};
    for (Int_t ch = 0; ch < 2; ch++) {
      const PulseHistory::Kernel &K = kernel->kernel[groupOf[ch]];
      summary.phRmsBefore[ch] = K.rms_before;
      summary.phRmsAfter[ch] = K.rms_after;
      summary.phR2[ch] = K.r2;
      // The first four bins (dt < 4.6 us) sit inside the trapezoid overlap of
      // the previous pulse and are populated in every run, control included;
      // the pole-zero effect is judged on the bins beyond.
      for (Int_t b = 4; b < PulseHistory::kNBins; b++)
        summary.kernelMaxAbs[ch] =
            TMath::Max(summary.kernelMaxAbs[ch], TMath::Abs(K.k[0][b]));
    }
  }

  // --- Histograms --------------------------------------------------------
  TH1D *specRaw[2], *specCorr[2];
  TProfile *profRaw[2], *profCorr[2], *profAll[2], *profFit[2], *profFitCorr[2];
  TProfile *lin57[2], *lin812[2];
  std::vector<Double_t> windowRaw[2], windowCorr[2];
  Long64_t nBelow12[2] = {0, 0}, nDt[2] = {0, 0};
  const Double_t sig = peak * hdr.beamPeakSigmaFrac;
  const Double_t linLo = peak - 4.0 * sig, linHi = peak + 4.0 * sig;
  for (Int_t ch = 0; ch < 2; ch++) {
    const char *cn = ChannelName(ch);
    specRaw[ch] = MakeSpectrum(Form("spec_raw_%s_%s", cn, label.Data()), peak);
    specCorr[ch] =
        MakeSpectrum(Form("spec_corr_%s_%s", cn, label.Data()), peak);
    profRaw[ch] = MakeProfile(Form("prof_raw_%s_%s", cn, label.Data()));
    profCorr[ch] = MakeProfile(Form("prof_corr_%s_%s", cn, label.Data()));
    profAll[ch] = MakeProfile(Form("prof_all_%s_%s", cn, label.Data()));
    profFit[ch] = MakeFitProfile(Form("prof_fit_%s_%s", cn, label.Data()));
    profFitCorr[ch] =
        MakeFitProfile(Form("prof_fitcorr_%s_%s", cn, label.Data()));
    lin57[ch] = new TProfile(
        Form("lin57_%s_%s", cn, label.Data()),
        ";Previous pulse true amplitude [ADC];Pulse height #minus true amplitude [ADC]",
        20, linLo, linHi);
    lin812[ch] = new TProfile(
        Form("lin812_%s_%s", cn, label.Data()),
        ";Previous pulse true amplitude [ADC];Pulse height #minus true amplitude [ADC]",
        20, linLo, linHi);
    lin57[ch]->SetDirectory(nullptr);
    lin812[ch]->SetDirectory(nullptr);
  }
  for (size_t i = 0; i < hits.size(); i++) {
    const SimHit &h = hits[i];
    const Int_t ch = h.channel;
    summary.channel[ch].nHits++;
    if (h.dtPrev >= 0.0) {
      nDt[ch]++;
      if (h.dtPrev < 12.0)
        nBelow12[ch]++;
      profAll[ch]->Fill(h.dtPrev, h.energy);
    }
    if (h.pileup)
      continue;
    specRaw[ch]->Fill(h.energy);
    specCorr[ch]->Fill(corrected[i]);
    // Profiles and fits use the beam window, like PulseHistory's beam-like
    // selection, so merged double pulses do not enter.
    const Bool_t inWindow =
        h.energy > kWindowLo * peak && h.energy < kWindowHi * peak;
    if (h.dtPrev >= 0.0 && inWindow) {
      profRaw[ch]->Fill(h.dtPrev, h.energy);
      profCorr[ch]->Fill(h.dtPrev, corrected[i]);
      profFit[ch]->Fill(h.dtPrev, h.energy);
      profFitCorr[ch]->Fill(h.dtPrev, corrected[i]);
      if (h.dtPrev >= 5.0 && h.dtPrev < 7.0)
        lin57[ch]->Fill(h.prevTrueAmp, h.energy - h.trueAmp);
      if (h.dtPrev >= 8.0 && h.dtPrev < 12.0)
        lin812[ch]->Fill(h.prevTrueAmp, h.energy - h.trueAmp);
    }
    if (inWindow)
      windowRaw[ch].push_back(h.energy);
    if (corrected[i] > kWindowLo * peak && corrected[i] < kWindowHi * peak)
      windowCorr[ch].push_back(corrected[i]);
  }
  Double_t truthFrac = 0.0;
  {
    TFile *tf = IO::OpenForReading(path);
    TTree *truth = tf ? dynamic_cast<TTree *>(tf->Get("truth")) : nullptr;
    if (truth) {
      Double_t time_ns = 0.0, prev = -1.0;
      Long64_t below = 0, total = 0;
      truth->SetBranchAddress("time_ns", &time_ns);
      for (Long64_t i = 0; i < truth->GetEntries(); i++) {
        truth->GetEntry(i);
        if (prev >= 0.0) {
          total++;
          if (time_ns - prev < 12000.0)
            below++;
        }
        prev = time_ns;
      }
      truthFrac = total ? Double_t(below) / Double_t(total) : 0.0;
    }
    if (tf) {
      tf->Close();
      delete tf;
    }
  }
  for (Int_t ch = 0; ch < 2; ch++) {
    ChannelSummary &cs = summary.channel[ch];
    cs.nWindow = Long64_t(windowRaw[ch].size());
    MeanRms(windowRaw[ch], cs.meanRaw, cs.rmsRaw);
    MeanRms(windowCorr[ch], cs.meanCorr, cs.rmsCorr);
    cs.fracDtBelow12 = nDt[ch] ? Double_t(nBelow12[ch]) / Double_t(nDt[ch])
                               : 0.0;
    cs.fracDtBelow12Truth = truthFrac;
  }

  // --- Spectra -----------------------------------------------------------
  DrawSpectra(specRaw[0], specRaw[1], summary.channel, kFALSE, label,
              "spectrum_LR_raw");
  if (kernel)
    DrawSpectra(specCorr[0], specCorr[1], summary.channel, kTRUE, label,
                "spectrum_LR_corrected");

  // --- Profiles and fits ---------------------------------------------------
  const Double_t ylo = 0.6 * peak, yhi = 1.4 * peak;
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    StyleProfile(profAll[0], kColorL);
    StyleProfile(profAll[1], kColorR);
    profAll[0]->SetMinimum(ylo);
    profAll[0]->SetMaximum(2.2 * peak);
    profAll[0]->Draw("E1");
    profAll[1]->Draw("E1 SAME");
    TLegend *leg = PlottingUtils::AddLegend(0.50, 0.89, 0.72, 0.87);
    leg->SetTextSize(kLegendTextSize);
    leg->AddEntry(profAll[0], "L, all accepted hits", "lp");
    leg->AddEntry(profAll[1], "R, all accepted hits", "lp");
    leg->Draw();
    PlottingUtils::AddText(label, 0.88, 0.92);
    PlottingUtils::SaveFigure(c, "profile_LR_all_hits", label,
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  for (Int_t ch = 0; ch < 2; ch++) {
    ChannelSummary &cs = summary.channel[ch];
    cs.expRaw = FitExponential(profFit[ch], fitLo, kFitMaxUs - 20.0);
    cs.lineRaw = FitLine(profRaw[ch], fitLo, kProfileMaxUs);
    if (kernel)
      cs.lineCorr = FitLine(profCorr[ch], fitLo, kProfileMaxUs);
    cs.linearity57 = FitLine(lin57[ch], linLo, linHi, 20);
    cs.linearity812 = FitLine(lin812[ch], linLo, linHi, 20);
  }
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    StyleProfile(profRaw[0], kColorL);
    StyleProfile(profRaw[1], kColorR);
    profRaw[0]->SetMinimum(ylo);
    profRaw[0]->SetMaximum(yhi);
    profRaw[0]->Draw("E1");
    profRaw[1]->Draw("E1 SAME");
    // The exponential is fitted on a wider profile than this one, so it is
    // drawn here from its parameters rather than carried on the histogram.
    for (Int_t ch = 0; ch < 2; ch++) {
      const ExpFit &e = summary.channel[ch].expRaw;
      if (IsFlat(e))
        continue;
      TF1 *shown = new TF1(PlottingUtils::GetRandomName(),
                           "[0] + [1] * exp(-x / [2]) + [3] * x", fitLo,
                           kProfileMaxUs);
      shown->SetParameters(e.p0, e.p1, e.p2, e.p3);
      shown->SetLineColor(kBlack);
      shown->SetLineWidth(2);
      shown->Draw("SAME");
    }
    TLegend *leg = PlottingUtils::AddLegend(0.30, 0.89, 0.72, 0.87);
    leg->SetTextSize(kLegendTextSize);
    for (Int_t ch = 0; ch < 2; ch++) {
      const ExpFit &e = summary.channel[ch].expRaw;
      if (IsFlat(e))
        leg->AddEntry(profRaw[ch], Form("%s: flat", ChannelName(ch)), "lp");
      else
        leg->AddEntry(profRaw[ch],
                      Form("%s: #tau_{fit} = %.1f #pm %.1f #mus, A = %.0f",
                           ChannelName(ch), e.p2, e.p2err, e.p1),
                      "lp");
    }
    leg->Draw();
    PlottingUtils::AddText(label, 0.88, 0.92);
    PlottingUtils::SaveFigure(c, "profile_LR_fit", label,
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  if (kernel) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    StyleProfile(profRaw[1], kColorR);
    StyleProfile(profCorr[1], kColorCorr);
    profRaw[1]->SetMinimum(ylo);
    profRaw[1]->SetMaximum(yhi);
    profRaw[1]->Draw("E1");
    profCorr[1]->Draw("E1 SAME");
    TLegend *leg = PlottingUtils::AddLegend(0.40, 0.89, 0.66, 0.87);
    leg->SetTextSize(kLegendTextSize);
    leg->AddEntry(profRaw[1],
                  Form("R raw: RMS %.1f", summary.channel[1].rmsRaw), "lp");
    leg->AddEntry(profCorr[1],
                  Form("R corrected: RMS %.1f", summary.channel[1].rmsCorr),
                  "lp");
    leg->AddEntry((TObject *)nullptr,
                  Form("kernel %s, residual slope %.2f #pm %.2f ADC/#mus",
                       summary.kernelSource.Data(),
                       summary.channel[1].lineCorr.slope,
                       summary.channel[1].lineCorr.slopeErr),
                  "");
    leg->Draw();
    PlottingUtils::AddText(label, 0.88, 0.92);
    PlottingUtils::SaveFigure(c, "profile_R_corrected", label,
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    StyleProfile(lin57[1], kColorR);
    StyleProfile(lin812[1], kColorCorr);
    Double_t lo = TMath::Min(lin57[1]->GetMinimum(), lin812[1]->GetMinimum());
    Double_t hi = TMath::Max(lin57[1]->GetMaximum(), lin812[1]->GetMaximum());
    lin57[1]->SetMinimum(lo - 0.2 * (hi - lo) - 20.0);
    lin57[1]->SetMaximum(hi + 0.4 * (hi - lo) + 20.0);
    lin57[1]->Draw("E1");
    lin812[1]->Draw("E1 SAME");
    TLegend *leg = PlottingUtils::AddLegend(0.30, 0.89, 0.72, 0.87);
    leg->SetTextSize(kLegendTextSize);
    const LineFit &a = summary.channel[1].linearity57;
    const LineFit &b = summary.channel[1].linearity812;
    leg->AddEntry(lin57[1],
                  Form("R, #Deltat 5-7 #mus: slope %.3f #pm %.3f", a.slope,
                       a.slopeErr),
                  "lp");
    leg->AddEntry(lin812[1],
                  Form("R, #Deltat 8-12 #mus: slope %.3f #pm %.3f", b.slope,
                       b.slopeErr),
                  "lp");
    leg->Draw();
    PlottingUtils::AddText(label, 0.88, 0.92);
    PlottingUtils::SaveFigure(c, "linearity_R", label,
                              PlotSaveOptions::kLINEAR);
    delete c;
  }

  // --- Waveform snippets ---------------------------------------------------
  {
    TFile *f = IO::OpenForReading(path);
    for (Int_t ch = 0; ch < 2 && f && !f->IsZombie(); ch++) {
      const char *cn = ChannelName(ch);
      TGraph *raw = dynamic_cast<TGraph *>(f->Get(Form("raw_%s", cn)));
      TGraph *trap = dynamic_cast<TGraph *>(f->Get(Form("trap_%s", cn)));
      TGraph *base = dynamic_cast<TGraph *>(f->Get(Form("baseline_%s", cn)));
      if (!raw || !trap || !base)
        continue;
      // The raw signal is shown without its DC offset and the trapezoid and
      // baseline without the constant the DC offset produces in the filter
      // (its value at the start of the snippet).
      const Double_t trapDc = trap->GetY()[0];
      TGraph *rawS = new TGraph(raw->GetN());
      TGraph *trapS = new TGraph(trap->GetN());
      TGraph *baseS = new TGraph(base->GetN());
      for (Int_t i = 0; i < raw->GetN(); i++) {
        rawS->SetPoint(i, raw->GetX()[i], raw->GetY()[i] - hdr.adcOffset);
        trapS->SetPoint(i, trap->GetX()[i], trap->GetY()[i] - trapDc);
        baseS->SetPoint(i, base->GetX()[i], base->GetY()[i] - trapDc);
      }
      TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
      PlottingUtils::ConfigureGraph(rawS, kGray + 1,
                                    ";t #minus t_{1} [#mus];Signal [ADC]");
      rawS->GetYaxis()->SetRangeUser(-0.6 * peak, 1.4 * peak);
      rawS->Draw("AL");
      PlottingUtils::ConfigureGraph(trapS, ChannelColor(ch));
      trapS->SetLineWidth(3);
      trapS->Draw("L SAME");
      PlottingUtils::ConfigureGraph(baseS, kBlack);
      baseS->SetLineStyle(2);
      baseS->SetLineWidth(3);
      baseS->Draw("L SAME");
      TLegend *leg = PlottingUtils::AddLegend(0.45, 0.89, 0.68, 0.87);
      leg->SetTextSize(kLegendTextSize);
      leg->AddEntry(rawS, "Preamp signal + noise", "l");
      leg->AddEntry(trapS, Form("Trapezoid output, %s", cn), "l");
      leg->AddEntry(baseS, "Baseline restorer", "l");
      leg->Draw();
      PlottingUtils::AddText(Form("%s, %s: #tau %.0f #mus, PZ %.0f #mus",
                                  label.Data(), cn,
                                  ch == 0 ? hdr.tauL_us : hdr.tauR_us,
                                  hdr.trapDecay_us),
                             0.88, 0.92);
      PlottingUtils::SaveFigure(c, Form("waveform_%s", cn), label,
                                PlotSaveOptions::kLINEAR);
      delete c;
      // Tail zoom.
      TCanvas *cz = PlottingUtils::GetConfiguredCanvas(kFALSE);
      trapS->GetYaxis()->SetRangeUser(-0.4 * peak, 0.2 * peak);
      trapS->SetTitle(";t #minus t_{1} [#mus];Trapezoid #minus DC term [ADC]");
      trapS->Draw("AL");
      baseS->Draw("L SAME");
      TLine *zero = new TLine(trapS->GetX()[0], 0.0,
                              trapS->GetX()[trapS->GetN() - 1], 0.0);
      zero->SetLineStyle(3);
      zero->Draw();
      TLegend *legz = PlottingUtils::AddLegend(0.50, 0.89, 0.72, 0.87);
      legz->SetTextSize(kLegendTextSize);
      legz->AddEntry(trapS, Form("Trapezoid output, %s", cn), "l");
      legz->AddEntry(baseS, "Baseline restorer", "l");
      legz->Draw();
      PlottingUtils::AddText(label, 0.88, 0.92);
      PlottingUtils::SaveFigure(cz, Form("trapezoid_tail_%s", cn), label,
                                PlotSaveOptions::kLINEAR);
      delete cz;
    }
    if (f) {
      f->Close();
      delete f;
    }
  }

  // --- Console report -------------------------------------------------------
  for (Int_t ch = 0; ch < 2; ch++) {
    const ChannelSummary &cs = summary.channel[ch];
    std::cout << "  " << ChannelName(ch) << ": " << cs.nHits << " hits, "
              << cs.nWindow << " non-pileup in window; mean " << std::fixed
              << std::setprecision(1) << cs.meanRaw << " RMS " << cs.rmsRaw
              << " -> corrected mean " << cs.meanCorr << " RMS " << cs.rmsCorr
              << " ADC; frac(dt < 12 us) " << std::setprecision(3)
              << cs.fracDtBelow12 << std::endl;
    if (IsFlat(cs.expRaw))
      std::cout << "     raw profile: flat (exponential amplitude "
                << std::setprecision(1) << cs.expRaw.p1 << " +- "
                << cs.expRaw.p1err << " ADC)";
    else
      std::cout << "     raw profile fit " << std::setprecision(1) << "p0 "
                << cs.expRaw.p0 << ", p1 " << cs.expRaw.p1 << " +- "
                << cs.expRaw.p1err << " ADC, tau " << std::setprecision(2)
                << cs.expRaw.p2 << " +- " << cs.expRaw.p2err
                << " us, p3 " << std::setprecision(3) << cs.expRaw.p3
                << " ADC/us (chi2/ndf " << std::setprecision(1)
                << cs.expRaw.chi2 << "/" << cs.expRaw.ndf << ")";
    std::cout << "; linear slope " << std::setprecision(3)
              << cs.lineRaw.slope << " +- " << cs.lineRaw.slopeErr
              << " ADC/us" << std::endl;
    if (kernel)
      std::cout << "     corrected profile: linear slope "
                << std::setprecision(3) << cs.lineCorr.slope << " +- "
                << cs.lineCorr.slopeErr << " ADC/us" << std::endl;
    std::cout << "     linearity (energy - trueAmp) vs prevTrueAmp: dt 5-7 us "
              << "slope " << std::setprecision(4) << cs.linearity57.slope
              << " +- " << cs.linearity57.slopeErr << " (chi2/ndf "
              << std::setprecision(1) << cs.linearity57.chi2 << "/"
              << cs.linearity57.ndf << "); dt 8-12 us slope "
              << std::setprecision(4) << cs.linearity812.slope << " +- "
              << cs.linearity812.slopeErr << " (chi2/ndf "
              << std::setprecision(1) << cs.linearity812.chi2 << "/"
              << cs.linearity812.ndf << ")" << std::endl;
  }
  std::cout.unsetf(std::ios::fixed);
  std::cout << std::setprecision(6);

  for (Int_t ch = 0; ch < 2; ch++) {
    delete specRaw[ch];
    delete specCorr[ch];
    delete profRaw[ch];
    delete profCorr[ch];
    delete profAll[ch];
    delete profFit[ch];
    delete profFitCorr[ch];
    delete lin57[ch];
    delete lin812[ch];
  }
  return summary;
}

void PrintSummaryTable(const std::vector<RunSummary> &rows) {
  std::cout << std::endl;
  std::cout << "Summary (RMS over non-pileup hits within "
            << kWindowLo << "-" << kWindowHi
            << " x beam peak; tau_fit from p0 + p1 exp(-dt/p2) on the "
               "energy-vs-dt profile with a linear baseline term, fitted from the end of the freeze)"
            << std::endl;
  std::cout << Form("%-18s %6s %3s %8s %8s %9s %9s %11s %8s %11s %9s %8s %8s",
                    "run", "kHz", "ch", "N", "mean", "RMS raw", "RMS corr",
                    "kernel", "tau_true", "tau_fit", "p1 [ADC]", "f<12 hit",
                    "f<12 tru")
            << std::endl;
  for (size_t i = 0; i < rows.size(); i++) {
    const RunSummary &r = rows[i];
    for (Int_t ch = 0; ch < 2; ch++) {
      const ChannelSummary &cs = r.channel[ch];
      const Double_t tauTrue = ch == 0 ? r.header.tauL_us : r.header.tauR_us;
      // A profile whose exponential amplitude is within 2 sigma of zero is
      // reported as flat rather than with a meaningless time constant.
      const TString tau = TauLabel(cs.expRaw);
      std::cout << Form("%-18s %6.1f %3s %8lld %8.1f %9.1f %9.1f %11s %8.1f "
                        "%11s %9.0f %8.3f %8.3f",
                        r.label.Data(), r.header.rate_kHz, ChannelName(ch),
                        cs.nWindow, cs.meanRaw, cs.rmsRaw, cs.rmsCorr,
                        r.kernelSource.Length() ? r.kernelSource.Data() : "-",
                        tauTrue, tau.Data(),
                        IsFlat(cs.expRaw) ? 0.0 : cs.expRaw.p1,
                        cs.fracDtBelow12, cs.fracDtBelow12Truth)
                << std::endl;
    }
  }
  std::cout << std::endl;
  std::cout << "PulseHistory kernels (max |k| over the dt bins beyond 4.6 us; rms over its "
               "beam-like selection; residual slope of the corrected "
               "energy-vs-dt profile)"
            << std::endl;
  std::cout << Form("%-18s %3s %11s %8s %9s %9s %6s %22s", "run", "ch",
                    "kernel", "max|k|", "rms bef", "rms aft", "R^2",
                    "corr slope [ADC/us]")
            << std::endl;
  for (size_t i = 0; i < rows.size(); i++) {
    const RunSummary &r = rows[i];
    for (Int_t ch = 0; ch < 2; ch++) {
      const ChannelSummary &cs = r.channel[ch];
      if (r.kernelMeasured)
        std::cout << Form("%-18s %3s %11s %8.3f %9.1f %9.1f %6.3f %12.2f +- %.2f",
                          r.label.Data(), ChannelName(ch),
                          r.kernelSource.Data(), r.kernelMaxAbs[ch],
                          r.phRmsBefore[ch], r.phRmsAfter[ch], r.phR2[ch],
                          cs.lineCorr.slope, cs.lineCorr.slopeErr)
                  << std::endl;
      else
        std::cout << Form("%-18s %3s %11s %8s %9s %9s %6s %12.2f +- %.2f",
                          r.label.Data(), ChannelName(ch),
                          r.kernelSource.Length() ? r.kernelSource.Data()
                                                  : "-",
                          "-", "-", "-", "-", cs.lineCorr.slope,
                          cs.lineCorr.slopeErr)
                  << std::endl;
    }
  }
}

} // namespace Analysis
