// Real-waveform check on SOLARIS run 70 (waveforms saved for every event,
// analog probe 0 = ADC input, analog probe 1 = trapezoid output, digital
// probe 0 = trigger). For the long L and long R ends it measures directly:
//   * the preamp decay time, from the ADC input after isolated pulses, per
//     channel and per side;
//   * the trapezoid tail after isolated pulses (undershoot or overshoot);
//   * the baseline the firmware used for each event, as the trapezoid at the
//     peaking point minus the reported energy, against the time since the
//     previous pulse and against the trapezoid level just before the trigger.
// Two passes over the files: the first collects per-hit features and the
// time to the previous hit on each channel, the second averages the traces
// of the isolated hits.
//
//   realwave --out=output/realwave/run70.root [--nfiles=N] <sol files...>
#include "BinaryUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace {

const Double_t kSample_ns = 8.0;
const Int_t kPreSamples = 625;       // ChPreTriggerT = 5000 ns
const Int_t kRiseSamples = 250;      // 2000 ns
const Int_t kPeakDelaySamples = 187; // 50 % of 3000 ns flat top
const Int_t kPeakAvg = 16;
const Int_t kPreWindow = 64;
const Double_t kIsolated_us = 60.0;
const Int_t kNLong = 8;
const Int_t kLongL[kNLong] = {1, 9, 17, 25, 33, 41, 49, 57};
const Int_t kLongR[kNLong] = {4, 12, 20, 28, 36, 44, 52, 60};
const Int_t kMaxCh = 64;

struct Rec {
  Int_t channel;
  Int_t side; // 0 = L, 1 = R
  Double_t t_ns;
  Double_t energy;
  Double_t inputPre, inputPeak;
  Double_t trapPre, trapPeak;
  Int_t nTrig;
  Double_t dtPrev_us;
  Bool_t isolated;
};

Int_t SideOf(Int_t ch) {
  for (Int_t i = 0; i < kNLong; i++) {
    if (kLongL[i] == ch)
      return 0;
    if (kLongR[i] == ch)
      return 1;
  }
  return -1;
}

Double_t Mean(const Int_t *v, Int_t from, Int_t n) {
  Double_t s = 0.0;
  for (Int_t i = 0; i < n; i++)
    s += Double_t(v[from + i]);
  return s / Double_t(n);
}

Bool_t RecBefore(const Rec &a, const Rec &b) { return a.t_ns < b.t_ns; }

struct ShapeSum {
  std::vector<Double_t> input, trap;
  Long64_t n = 0;
  void Init(Int_t len) {
    input.assign(len, 0.0);
    trap.assign(len, 0.0);
  }
};

TGraph *ShapeGraph(const std::vector<Double_t> &sum, Long64_t n,
                   const char *name, const char *title) {
  TGraph *g = new TGraph(Int_t(sum.size()));
  for (Int_t i = 0; i < Int_t(sum.size()); i++)
    g->SetPoint(i, (Double_t(i) - kPreSamples) * kSample_ns / 1000.0,
                n > 0 ? sum[i] / Double_t(n) : 0.0);
  g->SetName(name);
  g->SetTitle(title);
  return g;
}

struct DecayFit {
  Double_t tau = 0.0, tauErr = 0.0, amp = 0.0, offset = 0.0, chi2 = 0.0;
  Int_t ndf = 0;
};

DecayFit FitDecay(TGraph *g, Double_t lo, Double_t hi, Double_t amp0,
                  Double_t tau0, Bool_t draw) {
  TF1 *f = new TF1(PlottingUtils::GetRandomName(), "[0] * exp(-x / [1]) + [2]",
                   lo, hi);
  f->SetParameters(amp0, tau0, 0.0);
  f->SetParLimits(1, 0.5, 200.0);
  f->SetLineColor(kBlack);
  f->SetLineStyle(2);
  g->Fit(f, draw ? "RQ+" : "RQ0+");
  DecayFit d;
  d.amp = f->GetParameter(0);
  d.tau = f->GetParameter(1);
  d.tauErr = f->GetParError(1);
  d.offset = f->GetParameter(2);
  d.chi2 = f->GetChisquare();
  d.ndf = f->GetNDF();
  return d;
}

} // namespace

int main(int argc, char **argv) {
  TString out = "output/realwave/run70.root";
  Long64_t maxFiles = 1000;
  std::vector<TString> inputs;
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a.BeginsWith("--out="))
      out = a(6, a.Length() - 6);
    else if (a.BeginsWith("--nfiles="))
      maxFiles = TString(a(9, a.Length() - 9)).Atoll();
    else
      inputs.push_back(a);
  }
  if (inputs.empty()) {
    std::cerr << "usage: realwave [--out=file.root] [--nfiles=N] <sol files>"
              << std::endl;
    return 1;
  }
  if (Long64_t(inputs.size()) > maxFiles)
    inputs.resize(size_t(maxFiles));
  TString outDir = gSystem->DirName(out);
  TString absOut = gSystem->IsAbsoluteFileName(outDir)
                       ? outDir
                       : TString(gSystem->WorkingDirectory()) + "/" + outDir;
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, absOut,
                                gSystem->WorkingDirectory(), kFALSE);
  const Int_t peakStart =
      kPreSamples + kRiseSamples + kPeakDelaySamples - kPeakAvg / 2;

  // --- pass 1: features ------------------------------------------------------
  std::vector<Rec> recs;
  Int_t traceLen = 0;
  Long64_t nBlocks = 0, nTraces = 0;
  for (size_t f = 0; f < inputs.size(); f++) {
    SOLReader reader;
    if (!reader.Open(inputs[f].Data())) {
      std::cerr << "cannot open " << inputs[f] << std::endl;
      continue;
    }
    while (reader.ReadEvent()) {
      nBlocks++;
      const SOLData &ev = reader.GetCurrentEvent();
      if (!ev.hasTraces())
        continue;
      nTraces++;
      const Int_t side = SideOf(Int_t(ev.channel));
      if (side < 0)
        continue;
      const Int_t len = Int_t(ev.getSamples());
      const Int_t *in = ev.getAnalog0();
      const Int_t *tr = ev.getAnalog1();
      const UChar_t *trig = ev.getDigital(0);
      if (!in || !tr || !trig || len < peakStart + kPeakAvg)
        continue;
      if (traceLen == 0)
        traceLen = len;
      if (len != traceLen)
        continue;
      Rec r;
      r.channel = Int_t(ev.channel);
      r.side = side;
      r.t_ns = Double_t(ev.timestamp);
      r.energy = Double_t(ev.energy);
      r.inputPre = Mean(in, 0, 250); // -5 to -3 us
      r.trapPre = Mean(tr, kPreSamples - kPreWindow - 8, kPreWindow);
      r.trapPeak = Mean(tr, peakStart, kPeakAvg);
      Double_t mx = -1.0e9;
      for (Int_t i = 0; i < kPreSamples + 2 * kRiseSamples && i < len; i++)
        mx = TMath::Max(mx, Double_t(in[i]));
      r.inputPeak = mx;
      r.nTrig = 0;
      for (Int_t i = 1; i < len; i++)
        if (trig[i] && !trig[i - 1])
          r.nTrig++;
      r.dtPrev_us = -1.0;
      r.isolated = kFALSE;
      recs.push_back(r);
    }
  }
  std::cout << inputs.size() << " files, " << nBlocks << " blocks, " << nTraces
            << " with traces, " << recs.size() << " on long ends, trace length "
            << traceLen << std::endl;
  if (recs.empty())
    return 1;

  std::stable_sort(recs.begin(), recs.end(), RecBefore);
  std::vector<Double_t> last(kMaxCh, -1.0);
  std::vector<TH1D *> spec(kMaxCh, nullptr);
  for (size_t i = 0; i < recs.size(); i++) {
    Rec &r = recs[i];
    if (last[r.channel] >= 0.0)
      r.dtPrev_us = (r.t_ns - last[r.channel]) / 1000.0;
    last[r.channel] = r.t_ns;
    if (!spec[r.channel]) {
      spec[r.channel] = new TH1D(Form("spec_ch%d", r.channel), "", 128, 0.0,
                                 4096.0);
      spec[r.channel]->SetDirectory(nullptr);
    }
    spec[r.channel]->Fill(r.energy);
  }
  std::vector<Double_t> mode(kMaxCh, 0.0);
  for (Int_t c = 0; c < kMaxCh; c++) {
    if (!spec[c])
      continue;
    spec[c]->GetXaxis()->SetRange(9, 128);
    mode[c] = spec[c]->GetBinCenter(spec[c]->GetMaximumBin());
  }
  // Isolated: nothing on the channel in the previous 60 us, a single trigger
  // in the record (nothing in the following 27 us), beam-like energy.
  std::set<std::pair<Int_t, Long64_t> > isolatedKeys;
  Long64_t nIso[2] = {0, 0};
  for (size_t i = 0; i < recs.size(); i++) {
    Rec &r = recs[i];
    r.isolated = r.dtPrev_us > kIsolated_us && r.nTrig == 1 &&
                 r.energy > 0.7 * mode[r.channel] &&
                 r.energy < 1.3 * mode[r.channel];
    if (r.isolated) {
      isolatedKeys.insert(std::make_pair(r.channel, Long64_t(r.t_ns)));
      nIso[r.side]++;
    }
  }
  std::cout << "isolated beam-like hits: L " << nIso[0] << " R " << nIso[1]
            << std::endl;

  // --- pass 2: average shapes of isolated hits, per channel ------------------
  std::vector<ShapeSum> shapeCh(kMaxCh);
  ShapeSum shapeSide[2];
  for (Int_t c = 0; c < kMaxCh; c++)
    shapeCh[c].Init(traceLen);
  shapeSide[0].Init(traceLen);
  shapeSide[1].Init(traceLen);
  for (size_t f = 0; f < inputs.size(); f++) {
    SOLReader reader;
    if (!reader.Open(inputs[f].Data()))
      continue;
    while (reader.ReadEvent()) {
      const SOLData &ev = reader.GetCurrentEvent();
      if (!ev.hasTraces())
        continue;
      const Int_t side = SideOf(Int_t(ev.channel));
      if (side < 0)
        continue;
      if (!isolatedKeys.count(
              std::make_pair(Int_t(ev.channel), Long64_t(ev.timestamp))))
        continue;
      const Int_t len = Int_t(ev.getSamples());
      if (len != traceLen)
        continue;
      const Int_t *in = ev.getAnalog0();
      const Int_t *tr = ev.getAnalog1();
      const Double_t e = Double_t(ev.energy);
      const Double_t inPre = Mean(in, 0, 250);
      const Double_t trPre = Mean(tr, 0, 250);
      ShapeSum &sc = shapeCh[ev.channel];
      ShapeSum &ss = shapeSide[side];
      for (Int_t i = 0; i < len; i++) {
        const Double_t vi = (Double_t(in[i]) - inPre) / e;
        const Double_t vt = (Double_t(tr[i]) - trPre) / e;
        sc.input[i] += vi;
        sc.trap[i] += vt;
        ss.input[i] += vi;
        ss.trap[i] += vt;
      }
      sc.n++;
      ss.n++;
    }
  }

  // --- output --------------------------------------------------------------
  TFile *fout = IO::OpenForWriting(out);
  TTree *tree = new TTree("wave", "Run 70 waveform features, long ends");
  Rec b;
  tree->Branch("channel", &b.channel, "channel/I");
  tree->Branch("side", &b.side, "side/I");
  tree->Branch("t_ns", &b.t_ns, "t_ns/D");
  tree->Branch("energy", &b.energy, "energy/D");
  tree->Branch("inputPre", &b.inputPre, "inputPre/D");
  tree->Branch("inputPeak", &b.inputPeak, "inputPeak/D");
  tree->Branch("trapPre", &b.trapPre, "trapPre/D");
  tree->Branch("trapPeak", &b.trapPeak, "trapPeak/D");
  tree->Branch("nTrig", &b.nTrig, "nTrig/I");
  tree->Branch("dtPrev_us", &b.dtPrev_us, "dtPrev_us/D");
  tree->Branch("isolated", &b.isolated, "isolated/O");
  const char *sideName[2] = {"L", "R"};
  const Int_t sideColor[2] = {kBlue + 1, kRed + 1};
  TProfile *pEnergy[2], *pTrapPre[2], *pBaseUsed[2], *pTrack[2];
  TH1D *hBaseIso[2];
  for (Int_t s = 0; s < 2; s++) {
    pEnergy[s] = new TProfile(Form("pEnergy_%s", sideName[s]),
                              ";#Deltat to previous hit [#mus];Energy [ADC]",
                              60, 0.0, 60.0);
    pTrapPre[s] = new TProfile(
        Form("pTrapPre_%s", sideName[s]),
        ";#Deltat to previous hit [#mus];Trapezoid just before trigger [ADC]",
        60, 0.0, 60.0);
    pBaseUsed[s] = new TProfile(
        Form("pBaseUsed_%s", sideName[s]),
        ";#Deltat to previous hit [#mus];Baseline used by firmware [ADC]", 60,
        0.0, 60.0);
    pTrack[s] = new TProfile(
        Form("pTrack_%s", sideName[s]),
        ";#Deltat to previous hit [#mus];Baseline used #minus trapezoid before trigger [ADC]",
        60, 0.0, 60.0);
    hBaseIso[s] = new TH1D(Form("hBaseIso_%s", sideName[s]),
                           ";Baseline used #minus trapezoid before trigger [ADC];Isolated hits",
                           200, -500.0, 500.0);
  }
  // Reference: trapezoid level just before isolated pulses, per side.
  Double_t ref[2] = {0.0, 0.0};
  Long64_t nref[2] = {0, 0};
  for (size_t i = 0; i < recs.size(); i++)
    if (recs[i].isolated) {
      ref[recs[i].side] += recs[i].trapPre;
      nref[recs[i].side]++;
    }
  for (Int_t s = 0; s < 2; s++)
    ref[s] = nref[s] ? ref[s] / Double_t(nref[s]) : 0.0;
  for (size_t i = 0; i < recs.size(); i++) {
    b = recs[i];
    tree->Fill();
    const Int_t s = b.side;
    if (b.dtPrev_us < 0.0 || b.energy < 0.7 * mode[b.channel] ||
        b.energy > 1.3 * mode[b.channel])
      continue;
    const Double_t baseUsed = b.trapPeak - b.energy;
    pEnergy[s]->Fill(b.dtPrev_us, b.energy);
    pTrapPre[s]->Fill(b.dtPrev_us, b.trapPre - ref[s]);
    pBaseUsed[s]->Fill(b.dtPrev_us, baseUsed - ref[s]);
    pTrack[s]->Fill(b.dtPrev_us, baseUsed - b.trapPre);
    if (b.isolated)
      hBaseIso[s]->Fill(baseUsed - b.trapPre);
  }
  tree->Write();

  // --- shapes, decay fits ----------------------------------------------------
  std::cout << std::endl;
  std::cout << "Preamp decay from the ADC input of isolated pulses, fit "
               "A exp(-t/tau) + C over 0.5-26.5 us after the trigger; "
               "trapezoid tail fit over 9-26.5 us"
            << std::endl;
  std::cout << Form("%-8s %6s %8s %14s %8s %14s %8s", "channel", "n",
                    "mode", "tau_in [us]", "amp_in", "tau_trap [us]",
                    "amp_trap")
            << std::endl;
  for (Int_t s = 0; s < 2; s++) {
    for (Int_t k = 0; k < kNLong; k++) {
      const Int_t c = s == 0 ? kLongL[k] : kLongR[k];
      if (shapeCh[c].n < 50)
        continue;
      TGraph *gi = ShapeGraph(shapeCh[c].input, shapeCh[c].n,
                              Form("input_ch%d", c), "");
      TGraph *gt = ShapeGraph(shapeCh[c].trap, shapeCh[c].n,
                              Form("trap_ch%d", c), "");
      DecayFit di = FitDecay(gi, 0.5, 26.5, 0.8, 15.0, kFALSE);
      DecayFit dt = FitDecay(gt, 9.0, 26.5, 0.1, 12.0, kFALSE);
      gi->Write();
      gt->Write();
      std::cout << Form("%s%-6d %6lld %8.0f %6.2f +- %4.2f %8.3f %6.1f +- %4.1f %8.3f",
                        sideName[s], c, shapeCh[c].n, mode[c], di.tau,
                        di.tauErr, di.amp, dt.tau, dt.tauErr, dt.amp)
                << std::endl;
    }
  }
  TGraph *gIn[2], *gTr[2];
  DecayFit fitIn[2], fitTr[2];
  for (Int_t s = 0; s < 2; s++) {
    gIn[s] = ShapeGraph(shapeSide[s].input, shapeSide[s].n,
                        Form("input_%s", sideName[s]),
                        ";t #minus trigger [#mus];ADC input / energy");
    gTr[s] = ShapeGraph(shapeSide[s].trap, shapeSide[s].n,
                        Form("trap_%s", sideName[s]),
                        ";t #minus trigger [#mus];Trapezoid / energy");
    gIn[s]->Write();
    gTr[s]->Write();
  }
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.35, 0.89, 0.68, 0.87);
    for (Int_t s = 0; s < 2; s++) {
      PlottingUtils::ConfigureGraph(gIn[s], sideColor[s],
                                    ";t #minus trigger [#mus];ADC input / energy");
      gIn[s]->SetLineWidth(3);
      if (s == 0) {
        gIn[s]->GetYaxis()->SetRangeUser(-0.1, 1.2);
        gIn[s]->Draw("AL");
      } else {
        gIn[s]->Draw("L SAME");
      }
      fitIn[s] = FitDecay(gIn[s], 0.5, 26.5, 0.8, 15.0, kTRUE);
      leg->AddEntry(gIn[s],
                    Form("%s input, %lld isolated pulses: #tau = %.1f #pm %.1f #mus",
                         sideName[s], shapeSide[s].n, fitIn[s].tau,
                         fitIn[s].tauErr),
                    "l");
      std::cout << "  " << sideName[s] << " side input: tau = " << fitIn[s].tau
                << " +- " << fitIn[s].tauErr << " us, A " << fitIn[s].amp
                << ", C " << fitIn[s].offset << " (chi2/ndf " << fitIn[s].chi2
                << "/" << fitIn[s].ndf << ")" << std::endl;
    }
    leg->Draw();
    PlottingUtils::AddText("run 70, long ends", 0.88, 0.92);
    PlottingUtils::SaveFigure(c, "input_average_LR", "", PlotSaveOptions::kLINEAR);
    delete c;
  }
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.35, 0.89, 0.68, 0.87);
    for (Int_t s = 0; s < 2; s++) {
      PlottingUtils::ConfigureGraph(gTr[s], sideColor[s],
                                    ";t #minus trigger [#mus];Trapezoid / energy");
      gTr[s]->SetLineWidth(3);
      if (s == 0) {
        gTr[s]->GetYaxis()->SetRangeUser(-0.4, 1.3);
        gTr[s]->Draw("AL");
      } else {
        gTr[s]->Draw("L SAME");
      }
      fitTr[s] = FitDecay(gTr[s], 9.0, 26.5, 0.1, 12.0, kTRUE);
      leg->AddEntry(gTr[s],
                    Form("%s trapezoid, pole-zero %s #mus: tail A = %.3f, #tau = %.1f #pm %.1f #mus",
                         sideName[s], s == 0 ? "20" : "10", fitTr[s].amp,
                         fitTr[s].tau, fitTr[s].tauErr),
                    "l");
      std::cout << "  " << sideName[s] << " side trapezoid tail: A = "
                << fitTr[s].amp << ", tau = " << fitTr[s].tau << " +- "
                << fitTr[s].tauErr << " us, C " << fitTr[s].offset
                << std::endl;
    }
    leg->Draw();
    PlottingUtils::AddText("run 70, long ends", 0.88, 0.92);
    PlottingUtils::SaveFigure(c, "trapezoid_average_LR", "",
                              PlotSaveOptions::kLINEAR);
    delete c;
  }
  {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    for (Int_t s = 0; s < 2; s++) {
      gTr[s]->GetYaxis()->SetRangeUser(-0.15, 0.3);
      gTr[s]->Draw(s == 0 ? "AL" : "L SAME");
    }
    TLegend *leg = PlottingUtils::AddLegend(0.50, 0.89, 0.72, 0.87);
    leg->AddEntry(gTr[0], "L trapezoid (pole-zero 20 #mus)", "l");
    leg->AddEntry(gTr[1], "R trapezoid (pole-zero 10 #mus)", "l");
    leg->Draw();
    PlottingUtils::AddText("run 70, long ends", 0.88, 0.92);
    PlottingUtils::SaveFigure(c, "trapezoid_tail_zoom_LR", "",
                              PlotSaveOptions::kLINEAR);
    delete c;
  }

  // --- profiles --------------------------------------------------------------
  struct ProfFig {
    TProfile **p;
    const char *name;
    Double_t lo, hi;
  };
  ProfFig figs[4] = {{pEnergy, "energy_vs_dt", 800.0, 2400.0},
                     {pTrapPre, "trapezoid_before_trigger_vs_dt", -400.0, 600.0},
                     {pBaseUsed, "baseline_used_vs_dt", -400.0, 600.0},
                     {pTrack, "baseline_minus_trapezoid_vs_dt", -400.0, 600.0}};
  for (Int_t k = 0; k < 4; k++) {
    TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
    TLegend *leg = PlottingUtils::AddLegend(0.45, 0.89, 0.72, 0.87);
    for (Int_t s = 0; s < 2; s++) {
      TProfile *p = figs[k].p[s];
      PlottingUtils::ConfigureHistogram(p, sideColor[s]);
      p->SetMarkerColor(sideColor[s]);
      p->SetMarkerStyle(20);
      p->SetMinimum(figs[k].lo);
      p->SetMaximum(figs[k].hi);
      p->Draw(s == 0 ? "E1" : "E1 SAME");
      p->Write();
      leg->AddEntry(p, Form("%s long ends (baseline %s)", sideName[s],
                            s == 0 ? "Medium, 1024" : "High, 16384"),
                    "lp");
    }
    leg->Draw();
    PlottingUtils::AddText("run 70", 0.88, 0.92);
    PlottingUtils::SaveFigure(c, figs[k].name, "", PlotSaveOptions::kLINEAR);
    delete c;
  }
  for (Int_t s = 0; s < 2; s++) {
    hBaseIso[s]->Write();
    std::cout << "  " << sideName[s]
              << " isolated hits: baseline used - trapezoid before trigger = "
              << hBaseIso[s]->GetMean() << " +- " << hBaseIso[s]->GetRMS()
              << " ADC (n " << hBaseIso[s]->GetEntries() << ")" << std::endl;
  }
  fout->Close();
  delete fout;
  return 0;
}
