// Energy-vs-dt profiles per strip group from a SOL hit stream (no traces
// needed): long L, long R, short L, short R. Energies are normalised to the
// channel's beam peak and the deviation is reported relative to the 38-60 us
// level, which maps the sign and size of the pulse-history effect onto the
// (preamp tau, pole-zero) combination of each group.
//
//   realdt --label=run084 --out=output/realdt <sol file>
#include "BinaryUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TProfile.h>
#include <TSystem.h>
#include <algorithm>
#include <iostream>
#include <vector>

namespace {

const Int_t kMaxCh = 64;
const Int_t kNGroups = 4;
const char *kGroupName[kNGroups] = {"long L", "long R", "short L", "short R"};
const char *kGroupTag[kNGroups] = {"Llong", "Rlong", "Lshort", "Rshort"};
const Int_t kGroupColor[kNGroups] = {kBlue + 1, kRed + 1, kAzure + 7,
                                     kOrange + 7};

// Channel map of the 37Cl runs: R strips on 0, 4, 8, ...; L strips on 1, 5,
// 9, ...; strip n = channel / 4 + 1. Long ends: L odd strips, R even strips.
Int_t GroupOf(Int_t ch) {
  if (ch > 61 || (ch % 4) > 1)
    return -1;
  const Int_t strip = ch / 4 + 1;
  const Bool_t isL = (ch % 4) == 1;
  const Bool_t odd = (strip % 2) == 1;
  if (isL)
    return odd ? 0 : 2;
  return odd ? 3 : 1;
}

struct Hit {
  Int_t ch;
  Double_t t_ns;
  Double_t e;
};

Bool_t HitBefore(const Hit &a, const Hit &b) { return a.t_ns < b.t_ns; }

} // namespace

int main(int argc, char **argv) {
  TString label = "run";
  TString outDir = "output/realdt";
  std::vector<TString> inputs;
  Long64_t maxHits = 30000000;
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a.BeginsWith("--label="))
      label = a(8, a.Length() - 8);
    else if (a.BeginsWith("--out="))
      outDir = a(6, a.Length() - 6);
    else if (a.BeginsWith("--maxhits="))
      maxHits = TString(a(10, a.Length() - 10)).Atoll();
    else
      inputs.push_back(a);
  }
  if (inputs.empty()) {
    std::cerr << "usage: realdt --label=runNNN [--out=dir] <sol files>"
              << std::endl;
    return 1;
  }
  TString absOut = gSystem->IsAbsoluteFileName(outDir)
                       ? outDir
                       : TString(gSystem->WorkingDirectory()) + "/" + outDir;
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, absOut,
                                gSystem->WorkingDirectory(), kFALSE);

  std::vector<Hit> hits;
  hits.reserve(size_t(maxHits));
  for (size_t f = 0; f < inputs.size(); f++) {
    SOLReader reader;
    reader.SetSkipTraces(kTRUE);
    if (!reader.Open(inputs[f].Data())) {
      std::cerr << "cannot open " << inputs[f] << std::endl;
      continue;
    }
    while (reader.ReadEvent() && Long64_t(hits.size()) < maxHits) {
      const SOLData &ev = reader.GetCurrentEvent();
      if (GroupOf(Int_t(ev.channel)) < 0)
        continue;
      Hit h;
      h.ch = Int_t(ev.channel);
      h.t_ns = Double_t(ev.timestamp);
      h.e = Double_t(ev.energy);
      hits.push_back(h);
    }
  }
  std::cout << label << ": " << hits.size() << " strip hits" << std::endl;
  if (hits.empty())
    return 1;
  std::stable_sort(hits.begin(), hits.end(), HitBefore);
  const Double_t span_s = (hits.back().t_ns - hits.front().t_ns) * 1.0e-9;

  // Beam peak per channel: mode of the spectrum above 256 ADC.
  std::vector<TH1D *> spec(kMaxCh, nullptr);
  for (size_t i = 0; i < hits.size(); i++) {
    const Int_t c = hits[i].ch;
    if (!spec[c]) {
      spec[c] = new TH1D(Form("spec_%d", c), "", 128, 0.0, 4096.0);
      spec[c]->SetDirectory(nullptr);
    }
    spec[c]->Fill(hits[i].e);
  }
  std::vector<Double_t> mode(kMaxCh, 0.0);
  std::vector<Long64_t> nch(kMaxCh, 0);
  for (Int_t c = 0; c < kMaxCh; c++) {
    if (!spec[c])
      continue;
    spec[c]->GetXaxis()->SetRange(9, 128);
    mode[c] = spec[c]->GetBinCenter(spec[c]->GetMaximumBin());
    nch[c] = Long64_t(spec[c]->GetEntries());
  }

  TProfile *prof[kNGroups];
  for (Int_t g = 0; g < kNGroups; g++) {
    prof[g] = new TProfile(Form("p_%s_%s", kGroupTag[g], label.Data()),
                           ";#Deltat to previous hit on the channel [#mus];Energy / beam peak",
                           300, 0.0, 150.0);
    prof[g]->SetDirectory(nullptr);
  }
  std::vector<Double_t> last(kMaxCh, -1.0);
  for (size_t i = 0; i < hits.size(); i++) {
    const Hit &h = hits[i];
    const Int_t g = GroupOf(h.ch);
    if (last[h.ch] >= 0.0 && mode[h.ch] > 0.0) {
      const Double_t r = h.e / mode[h.ch];
      if (r > 0.5 && r < 1.6)
        prof[g]->Fill((h.t_ns - last[h.ch]) / 1000.0, r);
    }
    last[h.ch] = h.t_ns;
  }

  // Deviation relative to the 120-150 us level, in dt bands.
  const Int_t nb = 9;
  const Double_t lo[nb] = {5.0, 6.0, 8.0, 12.0, 18.0, 26.0, 38.0, 60.0, 90.0};
  const Double_t hi[nb] = {6.0, 8.0, 12.0, 18.0, 26.0, 38.0, 60.0, 90.0, 120.0};
  std::cout << "  rate proxy: " << hits.size() / span_s / 32.0
            << " hits/s per strip end over " << span_s << " s" << std::endl;
  std::cout << Form("  %-8s %8s %8s", "group", "peak", "n") ;
  for (Int_t b = 0; b < nb; b++)
    std::cout << Form(" %5.0f-%-4.0f", lo[b], hi[b]);
  std::cout << std::endl;
  for (Int_t g = 0; g < kNGroups; g++) {
    Double_t refSum = 0.0, refN = 0.0;
    for (Int_t bin = prof[g]->FindBin(120.0); bin <= prof[g]->FindBin(149.9);
         bin++) {
      refSum += prof[g]->GetBinContent(bin) * prof[g]->GetBinEntries(bin);
      refN += prof[g]->GetBinEntries(bin);
    }
    const Double_t ref = refN > 0.0 ? refSum / refN : 1.0;
    Double_t peakMean = 0.0;
    Int_t peakN = 0;
    for (Int_t c = 0; c < kMaxCh; c++)
      if (GroupOf(c) == g && mode[c] > 0.0) {
        peakMean += mode[c];
        peakN++;
      }
    std::cout << Form("  %-8s %8.0f %8.0f", kGroupName[g],
                      peakN ? peakMean / peakN : 0.0, prof[g]->GetEntries());
    for (Int_t b = 0; b < nb; b++) {
      Double_t s = 0.0, n = 0.0;
      for (Int_t bin = prof[g]->FindBin(lo[b]); bin < prof[g]->FindBin(hi[b]);
           bin++) {
        s += prof[g]->GetBinContent(bin) * prof[g]->GetBinEntries(bin);
        n += prof[g]->GetBinEntries(bin);
      }
      std::cout << Form(" %+10.3f", n > 0.0 ? s / n - ref : 0.0);
    }
    std::cout << std::endl;
  }

  // Decay constant of each group's dt dependence. The residual left by a
  // pole-zero mismatch decays with the preamp's real time constant whatever
  // the pole-zero setting is, and whatever fraction of the amplitude survives
  // the selection, so the long and short ends of one chain must give the same
  // tau even though their pole-zero settings and their amplitudes differ.
  // The start of the fit is scanned as a stability check. On the real data it
  // is not stable, because the profiles are not a single exponential: besides
  // the pole-zero tail they carry a rate correlation, since a hit with a long
  // gap behind it sits in a quiet stretch where the long baseline average is
  // less depressed by pile-up, which lifts the profile at large dt. Comparing
  // the whole shape against the simulation is the meaningful test.
  const Int_t nfit = 5;
  const Double_t fit_lo[nfit] = {8.0, 12.0, 16.0, 20.0, 25.0};
  const Double_t kFitHi_us = 120.0;
  std::cout << "  tau of the dt dependence, p0 + p1 exp(-dt/p2) up to "
            << kFitHi_us << " us, against the start of the fit:" << std::endl;
  std::cout << Form("    %-8s %9s", "group", "p1");
  for (Int_t k = 0; k < nfit; k++)
    std::cout << Form(" %13s", Form("from %.0f us", fit_lo[k]));
  std::cout << std::endl;
  for (Int_t g = 0; g < kNGroups; g++) {
    if (prof[g]->GetEntries() < 1000)
      continue;
    std::cout << Form("    %-8s", kGroupName[g]);
    Double_t amp = 0.0;
    TString cells;
    for (Int_t k = 0; k < nfit; k++) {
      TF1 *f = new TF1(Form("f_%s_%d", kGroupTag[g], k),
                       "[0] + [1] * exp(-x / [2])", fit_lo[k], kFitHi_us);
      f->SetParameters(1.0,
                       prof[g]->GetBinContent(prof[g]->FindBin(fit_lo[k])) - 1.0,
                       15.0);
      f->SetParLimits(2, 1.0, 200.0);
      f->SetLineColor(kGroupColor[g]);
      f->SetLineStyle(2);
      prof[g]->Fit(f, "RQ0");
      if (k == 0)
        amp = f->GetParameter(1);
      cells += Form(" %6.1f+-%-5.1f", f->GetParameter(2), f->GetParError(2));
    }
    std::cout << Form(" %+9.4f", amp) << cells << std::endl;
  }

  TCanvas *c = PlottingUtils::GetConfiguredCanvas(kFALSE);
  TLegend *leg = PlottingUtils::AddLegend(0.55, 0.89, 0.66, 0.87);
  for (Int_t g = 0; g < kNGroups; g++) {
    PlottingUtils::ConfigureHistogram(prof[g], kGroupColor[g]);
    prof[g]->SetMarkerColor(kGroupColor[g]);
    prof[g]->SetMarkerStyle(20);
    prof[g]->SetMarkerSize(0.7);
    prof[g]->SetMinimum(0.7);
    prof[g]->SetMaximum(1.4);
    prof[g]->Draw(g == 0 ? "E1" : "E1 SAME");
    leg->AddEntry(prof[g], kGroupName[g], "lp");
  }
  leg->Draw();
  PlottingUtils::AddText(label, 0.88, 0.92);
  PlottingUtils::SaveFigure(c, Form("energy_vs_dt_%s", label.Data()), "",
                            PlotSaveOptions::kLINEAR);
  TFile *fout = IO::OpenForWriting(outDir + "/" + label + ".root");
  for (Int_t g = 0; g < kNGroups; g++)
    prof[g]->Write();
  fout->Close();
  delete fout;
  delete c;
  return 0;
}
