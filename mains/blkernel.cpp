// blkernel - how much a pulse shifts the energy of later pulses on the same
// channel, against how long before them it came, measured from a hit stream.
// It shows how long the trapezoid baseline remembers a pulse, which tells a
// long averaging window (or a low-pass) from a short one.
//
// Every hit on a pole-zero mismatched channel leaves a trapezoid undershoot
// (R) or overshoot (L) that the baseline averages in for as long as the
// baseline remembers it. Regressing a hit's energy on the summed amplitude of
// the earlier hits on the same channel, in bins of their age, traces that
// memory directly: a boxcar gives a flat kernel with a sharp edge at the
// window length, a low-pass a smooth exponential. Only hits isolated from
// their neighbours (dtPrev > 40 us, dtNext > 10 us) enter as the response, so
// the previous pulse's own pole-zero tail is small. The regression is within
// channel (per-channel intercepts), so beam-peak calibration differences do
// not leak into the kernel. A far bin (default 60-200 ms) is a nuisance term
// that soaks up slow beam-intensity drifts; it is not plotted.
//
// The undershoot on R raises the energy and the overshoot on L lowers it, so
// a real baseline kernel flips sign between the R and L long ends, while a
// rate artefact (gain sag, pile-up) has the same sign on both.
//
//   blkernel --label=run097 [--out=dir] [--maxhits=N] [--group=Rlong|Llong|chN]
//            [--sim=label:file.root ...] <sol files>
#include "BinaryUtils.hpp"
#include "InitUtils.hpp"
#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVectorD.h>
#include <algorithm>
#include <deque>
#include <iostream>
#include <vector>

namespace {

const Int_t kMaxCh = 64;
// Age bin edges [us]. The last bin is the nuisance term.
const std::vector<Double_t> kAgeDefault = {40,   80,   130,  200,   260,   330,   450,
                                           650,  1000, 1500, 2500,  4000,  6000,  10000,
                                           16000, 25000, 40000, 60000, 200000};

struct Hit {
  Int_t ch;
  Double_t t_ns;
  Double_t e;
};

Bool_t HitBefore(const Hit &a, const Hit &b) { return a.t_ns < b.t_ns; }

// Same channel map as realdt: R strips on 0, 4, 8, ...; L strips on 1, 5, 9,
// ...; long ends are L odd strips and R even strips.
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

struct Kernel {
  TGraphErrors *g = nullptr;
  Long64_t nFit = 0;
};

// hits sorted in time; select(ch) picks the channels that enter.
// Only responses with nMin <= (hits on the channel in the preceding 1 ms) <
// nMax enter, which splits the fit by local rate.
template <typename Sel>
Kernel Measure(const std::vector<Hit> &hits, Sel select,
               const std::vector<Double_t> &age, Int_t nMin = 0,
               Int_t nMax = 1 << 30) {
  const Int_t nk = Int_t(age.size()) - 1;
  // Beam peak per channel: mode of the spectrum above 256 ADC.
  std::vector<TH1D *> spec(kMaxCh, nullptr);
  for (const Hit &h : hits) {
    if (!select(h.ch))
      continue;
    if (!spec[h.ch]) {
      spec[h.ch] = new TH1D(Form("bk_spec_%d", h.ch), "", 512, 0.0, 4096.0);
      spec[h.ch]->SetDirectory(nullptr);
    }
    spec[h.ch]->Fill(h.e);
  }
  std::vector<Double_t> mode(kMaxCh, 0.0);
  for (Int_t c = 0; c < kMaxCh; c++)
    if (spec[c]) {
      spec[c]->GetXaxis()->SetRange(33, 512);
      mode[c] = spec[c]->GetBinCenter(spec[c]->GetMaximumBin());
      delete spec[c];
    }

  // Per-channel history (hit times and the running sum of their amplitudes,
  // so a bin's summed amplitude is two binary searches) and per-channel cross
  // products [1, x_k, y].
  const Int_t np = nk + 2;
  std::vector<std::vector<Double_t>> histT(kMaxCh), histC(kMaxCh, std::vector<Double_t>(1, 0.0));
  std::vector<TMatrixDSym> S(kMaxCh, TMatrixDSym(np));
  // Candidate waiting for its dtNext: index into hits, and its regressors.
  std::vector<Long64_t> cand(kMaxCh, -1);
  std::vector<std::vector<Double_t>> candRow(kMaxCh);
  std::vector<Double_t> last(kMaxCh, -1.0);
  // Start of the channel's uninterrupted history; a gap over 0.5 s (a run
  // boundary) restarts it, and only responses with a full history enter.
  std::vector<Double_t> start(kMaxCh, -1.0);
  const Double_t maxAge_ns = age.back() * 1000.0;

  auto commit = [&](Int_t c) {
    const std::vector<Double_t> &r = candRow[c];
    for (Int_t a = 0; a < np; a++)
      for (Int_t b = 0; b <= a; b++)
        S[c](a, b) += r[a] * r[b];
  };

  for (size_t i = 0; i < hits.size(); i++) {
    const Hit &h = hits[i];
    const Int_t c = h.ch;
    if (!select(c) || mode[c] <= 0.0)
      continue;
    const Double_t en = h.e / mode[c];
    // Resolve the pending candidate with this hit as its successor.
    if (cand[c] >= 0 && h.t_ns - hits[cand[c]].t_ns > 10000.0)
      commit(c);
    cand[c] = -1;
    const std::vector<Double_t> &T = histT[c];
    const std::vector<Double_t> &C = histC[c];
    // Number of earlier hits on the channel at or before time x.
    auto upto = [&T](Double_t x) {
      return Long64_t(std::upper_bound(T.begin(), T.end(), x) - T.begin());
    };
    if (last[c] < 0.0 || h.t_ns - last[c] > 5.0e8)
      start[c] = h.t_ns;
    const Bool_t full = last[c] >= 0.0 && h.t_ns - start[c] > maxAge_ns;
    const Int_t recent = Int_t(Long64_t(T.size()) - upto(h.t_ns - 1.0e6));
    if (full && h.t_ns - last[c] > 40000.0 && en > 0.5 && en < 1.6 &&
        recent >= nMin && recent < nMax) {
      std::vector<Double_t> &r = candRow[c];
      r.assign(np, 0.0);
      r[0] = 1.0;
      // Hits aged [age_k, age_k+1) have times in (t - age_k+1, t - age_k].
      for (Int_t k = 0; k < nk; k++)
        r[k + 1] = C[upto(h.t_ns - age[k] * 1000.0)] - C[upto(h.t_ns - age[k + 1] * 1000.0)];
      r[np - 1] = en;
      cand[c] = Long64_t(i);
    }
    histT[c].push_back(h.t_ns);
    histC[c].push_back(histC[c].back() + std::min(std::max(en, 0.0), 3.0));
    last[c] = h.t_ns;
  }

  // Within-channel pooled normal equations: centre each channel's cross
  // products on its own means, then sum.
  TMatrixDSym A(nk);
  TVectorD b(nk);
  Double_t yy = 0.0;
  Long64_t n = 0;
  for (Int_t c = 0; c < kMaxCh; c++) {
    const Double_t m = S[c](0, 0);
    if (m < 50)
      continue;
    n += Long64_t(m);
    for (Int_t i = 0; i < nk; i++) {
      for (Int_t j = 0; j <= i; j++) {
        const Double_t v = S[c](i + 1, j + 1) - S[c](i + 1, 0) * S[c](j + 1, 0) / m;
        A(i, j) += v;
        if (i != j)
          A(j, i) += v;
      }
      b(i) += S[c](np - 1, i + 1) - S[c](np - 1, 0) * S[c](i + 1, 0) / m;
    }
    yy += S[c](np - 1, np - 1) - S[c](np - 1, 0) * S[c](np - 1, 0) / m;
  }
  Kernel k;
  k.nFit = n;
  k.g = new TGraphErrors();
  if (n < 1000)
    return k;
  TMatrixDSym inv(A);
  inv.Invert();
  const TVectorD beta = inv * b;
  const Double_t s2 = (yy - beta * b) / Double_t(n - nk - 1);
  for (Int_t i = 0; i < nk - 1; i++) { // last bin is the nuisance term
    const Int_t p = k.g->GetN();
    k.g->SetPoint(p, TMath::Sqrt(age[i] * age[i + 1]), beta(i));
    k.g->SetPointError(p, 0.0, TMath::Sqrt(s2 * inv(i, i)));
  }
  return k;
}

void Print(const char *name, const Kernel &k) {
  std::cout << Form("  %-24s n=%-9lld", name, k.nFit);
  for (Int_t i = 0; i < k.g->GetN(); i++)
    std::cout << Form(" %+.5f(%.5f)", k.g->GetY()[i], k.g->GetEY()[i]);
  std::cout << std::endl;
}

std::vector<Hit> ReadSim(const TString &file) {
  std::vector<Hit> hits;
  TFile *f = TFile::Open(file);
  if (!f || f->IsZombie())
    return hits;
  TTree *t = dynamic_cast<TTree *>(f->Get("hits"));
  Int_t ch = 0;
  ULong64_t ts = 0;
  Double_t e = 0;
  t->SetBranchAddress("channel", &ch);
  t->SetBranchAddress("timestamp", &ts);
  t->SetBranchAddress("energy", &e);
  for (Long64_t i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    // ps -> ns; the simulated R channel (1) is tagged as R long channel 4,
    // the simulated L channel (0) as L long channel 1.
    if (ch == 1)
      hits.push_back({4, Double_t(ts) * 1.0e-3, e});
    else if (ch == 0)
      hits.push_back({1, Double_t(ts) * 1.0e-3, e});
  }
  f->Close();
  return hits;
}

} // namespace

int main(int argc, char **argv) {
  TString label = "run";
  TString outDir = "output/blkernel";
  TString group = "Rlong";
  Long64_t maxHits = 200000000;
  std::vector<TString> inputs, sims;
  // --slices=a,b: also fit local-rate slices (hits in the preceding 1 ms
  // below a, and at or above b) for R and L, data and every simulation, and
  // write them all to --scoreout for a joint model comparison.
  Int_t sliceLo = -1, sliceHi = -1;
  TString scoreOut;
  for (Int_t i = 1; i < argc; i++) {
    TString a = argv[i];
    if (a.BeginsWith("--label="))
      label = a(8, a.Length() - 8);
    else if (a.BeginsWith("--out="))
      outDir = a(6, a.Length() - 6);
    else if (a.BeginsWith("--group="))
      group = a(8, a.Length() - 8);
    else if (a.BeginsWith("--maxhits="))
      maxHits = TString(a(10, a.Length() - 10)).Atoll();
    else if (a.BeginsWith("--slices=")) {
      TString v = a(9, a.Length() - 9);
      sliceLo = TString(v(0, v.Index(","))).Atoi();
      sliceHi = TString(v(v.Index(",") + 1, v.Length())).Atoi();
    } else if (a.BeginsWith("--scoreout="))
      scoreOut = a(11, a.Length() - 11);
    else if (a.BeginsWith("--sim="))
      sims.push_back(a(6, a.Length() - 6));
    else
      inputs.push_back(a);
  }
  TString absOut = gSystem->IsAbsoluteFileName(outDir)
                       ? outDir
                       : TString(gSystem->WorkingDirectory()) + "/" + outDir;
  InitUtils::SetROOTPreferences(PlotSaveFormat::kPNG, absOut,
                                gSystem->WorkingDirectory(), kFALSE);
  const std::vector<Double_t> &age = kAgeDefault;

  std::vector<Hit> hits;
  // Timestamps restart with every run: a file whose clock starts below what
  // has been read so far is shifted 1 s past it, so runs never interleave and
  // no history (max age 200 ms) reaches across a run boundary.
  Double_t tMax = -1.0, offset = 0.0;
  // --group=chN keeps only that channel in memory.
  const Int_t onlyCh = group.BeginsWith("ch") ? TString(group(2, group.Length() - 2)).Atoi() : -1;
  for (size_t f = 0; f < inputs.size() && Long64_t(hits.size()) < maxHits; f++) {
    SOLReader reader;
    reader.SetSkipTraces(kTRUE);
    if (!reader.Open(inputs[f].Data())) {
      std::cerr << "cannot open " << inputs[f] << std::endl;
      continue;
    }
    Bool_t first = kTRUE;
    while (reader.ReadEvent() && Long64_t(hits.size()) < maxHits) {
      const SOLData &ev = reader.GetCurrentEvent();
      if (GroupOf(Int_t(ev.channel)) < 0 || (onlyCh >= 0 && Int_t(ev.channel) != onlyCh))
        continue;
      const Double_t ts = Double_t(ev.timestamp);
      if (first && ts + offset < tMax)
        offset = tMax + 1.0e9 - ts;
      first = kFALSE;
      hits.push_back({Int_t(ev.channel), ts + offset, Double_t(ev.energy)});
      tMax = TMath::Max(tMax, ts + offset);
    }
  }
  std::stable_sort(hits.begin(), hits.end(), HitBefore);
  std::cout << label << ": " << hits.size() << " strip hits" << std::endl;
  if (!hits.empty()) {
    const Double_t span_s = (hits.back().t_ns - hits.front().t_ns) * 1.0e-9;
    std::vector<Long64_t> nch(kMaxCh, 0);
    for (const Hit &h : hits)
      nch[h.ch]++;
    if (onlyCh >= 0)
      std::cout << Form("  ch %d: %.1f kHz (live span, run gaps removed)", onlyCh,
                        [&] {
                          Double_t live = 0.0;
                          for (size_t i = 1; i < hits.size(); i++) {
                            const Double_t d = hits[i].t_ns - hits[i - 1].t_ns;
                            if (d < 5.0e8)
                              live += d;
                          }
                          return live > 0 ? hits.size() / (live * 1e-9) / 1e3 : 0.0;
                        }())
                << std::endl;
    for (Int_t g = 0; g < 2; g++) {
      Long64_t n = 0;
      Int_t used = 0;
      for (Int_t c = 0; c < kMaxCh; c++)
        if (GroupOf(c) == g && nch[c] > 0) {
          n += nch[c];
          used++;
        }
      std::cout << Form("  %s long: %d channels, %.1f kHz per channel over %.1f s",
                        g ? "R" : "L", used, used ? n / span_s / used / 1e3 : 0.0, span_s)
                << std::endl;
    }
    // Amplitude mix relative to each channel's beam peak, in PulseHistory's
    // bands [0, 0.5), [0.5, 1.5), [1.5, 2.5), >= 2.5: what an amplitude-banded
    // kernel would have to follow.
    std::vector<TH1D *> spec(kMaxCh, nullptr);
    for (const Hit &h : hits) {
      if (GroupOf(h.ch) < 0)
        continue;
      if (!spec[h.ch]) {
        spec[h.ch] = new TH1D(Form("bk_mix_%d", h.ch), "", 512, 0.0, 4096.0);
        spec[h.ch]->SetDirectory(nullptr);
      }
      spec[h.ch]->Fill(h.e);
    }
    std::vector<Double_t> peak(kMaxCh, 0.0);
    for (Int_t c = 0; c < kMaxCh; c++)
      if (spec[c]) {
        spec[c]->GetXaxis()->SetRange(33, 512);
        peak[c] = spec[c]->GetBinCenter(spec[c]->GetMaximumBin());
        delete spec[c];
      }
    for (Int_t g = 0; g < 2; g++) {
      Double_t band[4] = {0, 0, 0, 0};
      for (const Hit &h : hits) {
        if (GroupOf(h.ch) != g || peak[h.ch] <= 0.0)
          continue;
        band[TMath::Min(3, Int_t(h.e / peak[h.ch] + 0.5))] += 1.0;
      }
      const Double_t tot = band[0] + band[1] + band[2] + band[3];
      if (tot > 0)
        std::cout << Form("  %s long amplitude mix: <0.5 %.3f, ~1 %.3f, ~2 %.3f, "
                          ">=2.5 %.3f (x beam peak)",
                          g ? "R" : "L", band[0] / tot, band[1] / tot,
                          band[2] / tot, band[3] / tot)
                  << std::endl;
    }
  }

  // Data kernels: the requested group, plus the opposite-sign control.
  std::vector<std::pair<TString, Kernel>> data;
  auto byGroup = [](Int_t g) { return [g](Int_t c) { return GroupOf(c) == g; }; };
  if (!hits.empty()) {
    if (group.BeginsWith("ch")) {
      const Int_t want = TString(group(2, group.Length() - 2)).Atoi();
      data.push_back({Form("%s ch %d", label.Data(), want),
                      Measure(hits, [want](Int_t c) { return c == want; }, age)});
    } else {
      data.push_back({label + " R long", Measure(hits, byGroup(1), age)});
      data.push_back({label + " L long (sign control)", Measure(hits, byGroup(0), age)});
      // R - L: the baseline kernel adds with opposite signs, a rate artefact
      // with the same sign, so the difference keeps the former only.
      Kernel d;
      d.nFit = data[0].second.nFit;
      d.g = new TGraphErrors();
      for (Int_t i = 0; i < data[0].second.g->GetN(); i++) {
        const TGraphErrors *r = data[0].second.g, *l = data[1].second.g;
        d.g->SetPoint(i, r->GetX()[i], r->GetY()[i] - l->GetY()[i]);
        d.g->SetPointError(i, 0.0, TMath::Hypot(r->GetEY()[i], l->GetEY()[i]));
      }
      data.push_back({label + " R #minus L", d});
    }
  }
  std::vector<std::pair<TString, Kernel>> simk;
  for (const TString &s : sims) {
    const Int_t colon = s.Index(":");
    const TString name = s(0, colon), file = s(colon + 1, s.Length() - colon - 1);
    std::vector<Hit> sh = ReadSim(file);
    std::stable_sort(sh.begin(), sh.end(), HitBefore);
    simk.push_back({"sim: " + name, Measure(sh, [](Int_t c) { return c == 4; }, age)});
  }

  if (!scoreOut.IsNull()) {
    // Every source x {R, L} x {all, low-rate, high-rate}, named
    // <source>__<R|L>__<all|lo|hi>; the source is "data" or the sim name.
    TFile *fs = IO::OpenForWriting(scoreOut);
    auto writeSet = [&](const TString &src, const std::vector<Hit> &hv) {
      const Int_t lims[3][2] = {{0, 1 << 30}, {0, sliceLo}, {sliceHi, 1 << 30}};
      const char *sl[3] = {"all", "lo", "hi"};
      for (Int_t g = 0; g < 2; g++)
        for (Int_t k = 0; k < (sliceLo >= 0 ? 3 : 1); k++) {
          Kernel kk = Measure(hv, [g](Int_t c) { return GroupOf(c) == (g ? 0 : 1); },
                              age, lims[k][0], lims[k][1]);
          fs->cd();
          kk.g->Write(Form("%s__%s__%s", src.Data(), g ? "L" : "R", sl[k]));
          std::cout << Form("  [score] %-22s %s %-3s n=%lld", src.Data(), g ? "L" : "R",
                            sl[k], kk.nFit)
                    << std::endl;
        }
    };
    if (!hits.empty())
      writeSet("data", hits);
    for (const TString &s : sims) {
      const Int_t colon = s.Index(":");
      std::vector<Hit> sh = ReadSim(s(colon + 1, s.Length() - colon - 1));
      std::stable_sort(sh.begin(), sh.end(), HitBefore);
      writeSet(s(0, colon), sh);
    }
    fs->Close();
  }

  std::cout << "  kernel [energy / beam peak per earlier beam-peak hit], age bins (us):";
  for (size_t i = 0; i + 2 < age.size(); i++)
    std::cout << Form(" %.0f-%.0f", age[i], age[i + 1]);
  std::cout << std::endl;
  for (auto &d : data)
    Print(d.first, d.second);
  for (auto &d : simk)
    Print(d.first, d.second);

  // Figure: kernels on a log age axis, data as points, simulation as lines.
  std::vector<TString> labels;
  for (auto &d : data)
    labels.push_back(d.first);
  for (auto &d : simk)
    labels.push_back(d.first);
  const Int_t legendW = Int_t(PlottingUtils::LegendWidthPx(labels));
  TPad *plot, *legpad;
  TCanvas *c = PlottingUtils::GetConfiguredCanvasWithSideLegend(plot, legpad, legendW);
  plot->SetLeftMargin(0.19);
  plot->SetLogx();
  plot->cd();
  Double_t ylo = 0.0, yhi = 0.0;
  for (auto *v : {&data, &simk})
    for (auto &d : *v)
      for (Int_t i = 0; i < d.second.g->GetN(); i++) {
        ylo = std::min(ylo, d.second.g->GetY()[i] - d.second.g->GetEY()[i]);
        yhi = std::max(yhi, d.second.g->GetY()[i] + d.second.g->GetEY()[i]);
      }
  const Double_t pad = 0.1 * (yhi - ylo);
  // The kernels are flat zero well before the far bins; show up to 4 ms.
  const Double_t xhi = std::min(4000.0, age[age.size() - 2]);
  TH1F *frame = plot->DrawFrame(age.front(), ylo - pad, xhi, yhi + pad);
  frame->GetXaxis()->SetTitle("age of earlier hit on the channel [#mus]");
  frame->GetYaxis()->SetTitle("#DeltaE per earlier hit  (beam units)");
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  TLine *zero = new TLine(age.front(), 0.0, xhi, 0.0);
  zero->SetLineStyle(2);
  zero->SetLineColor(kGray + 2);
  zero->Draw();
  const Int_t dataCol[] = {kBlack, kGray + 1, kRed + 1, kBlue + 1};
  const Int_t simCol[] = {kBlue + 1, kGreen + 2, kOrange + 7, kMagenta + 1, kCyan + 2};
  legpad->cd();
  TLegend *leg = PlottingUtils::AddLegend(0.02, 0.98, 0.50, 0.95);
  plot->cd();
  for (size_t i = 0; i < simk.size(); i++) {
    TGraphErrors *g = simk[i].second.g;
    g->SetLineColor(simCol[i % 5]);
    g->SetFillColorAlpha(simCol[i % 5], 0.25);
    g->SetLineWidth(PlottingUtils::GetLineWidth());
    g->Draw("L3 SAME");
    leg->AddEntry(g, simk[i].first, "L");
  }
  for (size_t i = 0; i < data.size(); i++) {
    TGraphErrors *g = data[i].second.g;
    g->SetMarkerStyle(i == 0 ? 20 : 24);
    g->SetMarkerColor(dataCol[i % 4]);
    g->SetLineColor(dataCol[i % 4]);
    g->SetMarkerSize(1.2);
    g->Draw("PE SAME");
    leg->AddEntry(g, data[i].first, "PE");
  }
  legpad->cd();
  leg->Draw();
  PlottingUtils::ScaleFigure(plot);
  PlottingUtils::SaveFigure(c, Form("baseline_kernel_%s", label.Data()), "",
                            PlotSaveOptions::kLINEAR);
  TFile *fout = IO::OpenForWriting(outDir + "/" + label + "_kernel.root");
  for (auto *v : {&data, &simk})
    for (auto &d : *v)
      d.second.g->Write(TString(d.first).ReplaceAll(" ", "_").ReplaceAll(":", ""));
  fout->Close();
  return 0;
}
