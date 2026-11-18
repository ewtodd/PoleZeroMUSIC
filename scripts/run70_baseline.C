// Step 2 of scripts/run70_baseline.sh: what run 70's waveforms say about the
// DPP-PHA baseline restorer with EnergyFilterLFLimitation = On.
//
// In run 70 both long groups read ~22 us preamps; the right group ran a 10 us
// pole-zero, so its trapezoid overshoots after each pulse, and the left group
// a 20 us one, near matched. For every traced hit the firmware's baseline is
// B = trapezoid at the reading point - e / k, the trapezoid taken from the
// EnergyFilter probe (not baseline-subtracted) and k the energy scale. Both are
// found here on isolated hits (no hit on the channel in the previous 300 us):
// e is regressed on the trapezoid's 16-sample mean at each sample of the
// extracted window, less the pre-trigger level, and the sample with the
// smallest residual is the reading point.
//
// The observable is dB = B2 - B1 between a hit and the next one on the same
// channel, which needs neither the pre-trigger level nor a resting reference:
// a baseline frozen between the two gives 0. When the second hit has no record
// of its own (dt < 27 us) its trapezoid is read off the first hit's record.
//
// CAEN's description is fitted to dB(dt): frozen from the first trigger; at
// the end of each 7 us freeze (from the start of the trapezoid's rise) the
// baseline resumes if |trapezoid - B| < thr, else freezes again; once running
// it is the mean of the last W samples. The trapezoid is the measured average
// isolated shape per unit energy times the first hit's energy, its tail
// extended past the record by an exponential fitted at 18-31 us. thr is in
// probe units, close to ADC (k ~ 1).
//
// Writes LaTeX macros for the paper's text and a plot of dB against dt.
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>
#include <TTree.h>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <vector>

namespace {
const int NDS = 160, W0 = 950, NW = 300;
struct Hit {
  int ch, traced;
  unsigned short e;
  unsigned long long t;
  float pre1;
  std::vector<float> ds1;
  std::vector<short> win;
  double B = 0;
};
// The 37Cl map: right long ends on channels 8n + 4, left on 8n + 1.
bool IsRight(int c) { return c % 8 == 4; }
bool IsLeft(int c) { return c % 8 == 1; }
// The trapezoid's 16-sample mean at sample s, relative to the pre-trigger
// level.
double Mean16(const Hit &h, int s) {
  double m = 0;
  for (int q = -8; q < 8; q++)
    m += h.win[s - W0 + q];
  return m / 16;
}
// The 0.2 us trapezoid means, interpolated at sample s.
double Coarse(const Hit &h, double s) {
  const double x = s / 25.0 - 0.5;
  const int i = int(TMath::Floor(x));
  if (i < 0 || i + 1 >= NDS)
    return NAN;
  const double f = x - i;
  return h.ds1[i] * (1 - f) + h.ds1[i + 1] * f;
}
struct Pair {
  double dt, e1, dB;
};
} // namespace

void run70_baseline(const char *hitsPath, const char *macroPath,
                    const char *plotPath) {
  TFile f(hitsPath);
  TTree *t = (TTree *)f.Get("hits");
  int ch, traced;
  unsigned short e;
  unsigned long long ts;
  float pre1, ds1[NDS];
  short win[NW];
  t->SetBranchAddress("ch", &ch);
  t->SetBranchAddress("hasTrace", &traced);
  t->SetBranchAddress("e", &e);
  t->SetBranchAddress("t", &ts);
  t->SetBranchAddress("pre1", &pre1);
  t->SetBranchAddress("ds1", ds1);
  t->SetBranchAddress("win", win);
  std::map<int, std::vector<Hit>> byChannel;
  for (long i = 0; i < t->GetEntries(); i++) {
    t->GetEntry(i);
    if (!IsRight(ch) && !IsLeft(ch))
      continue;
    Hit h;
    h.ch = ch;
    h.traced = traced;
    h.e = e;
    h.t = ts;
    h.pre1 = pre1;
    if (traced) {
      h.ds1.assign(ds1, ds1 + NDS);
      h.win.assign(win, win + NW);
    }
    byChannel[ch].push_back(h);
  }
  for (auto &kv : byChannel)
    std::sort(kv.second.begin(), kv.second.end(),
              [](const Hit &a, const Hit &b) { return a.t < b.t; });

  // [0] right, [1] left.
  const char *groupName[2] = {"right (10 us pole-zero)", "left (20 us pole-zero)"};
  int sRead[2];
  double kE[2];
  for (int g = 0; g < 2; g++) {
    double best = 1e99;
    for (int s = W0 + 8; s < W0 + NW - 8; s += 2) {
      double n = 0, sx = 0, sxx = 0, sxe = 0, se = 0, see = 0;
      for (auto &kv : byChannel) {
        if ((g == 0) != IsRight(kv.first))
          continue;
        const std::vector<Hit> &v = kv.second;
        for (size_t j = 1; j < v.size(); j++) {
          if (!v[j].traced || v[j].t - v[j - 1].t < 300000)
            continue;
          const double x = Mean16(v[j], s);
          n++;
          sx += x;
          sxx += x * x;
          sxe += x * v[j].e;
          se += v[j].e;
          see += double(v[j].e) * v[j].e;
        }
      }
      const double mx = sx / n, me = se / n;
      const double vxx = sxx / n - mx * mx, vxe = sxe / n - mx * me,
                   vee = see / n - me * me;
      const double residual = TMath::Sqrt(TMath::Max(0.0, vee - vxe * vxe / vxx));
      if (residual < best) {
        best = residual;
        sRead[g] = s;
        kE[g] = vxe / vxx;
      }
    }
    printf("%s: reading point sample %d, e = %.4f x trapezoid (residual %.1f)\n",
           groupName[g], sRead[g], kE[g], best);
  }
  for (auto &kv : byChannel) {
    const int g = IsRight(kv.first) ? 0 : 1;
    for (auto &h : kv.second)
      if (h.traced)
        h.B = h.pre1 + Mean16(h, sRead[g]) - h.e / kE[g];
  }

  const int NB = 14;
  const double edges[NB + 1] = {7,  10, 13, 16,  19,  22,  26.5, 35,
                                45, 60, 80, 110, 160, 250, 400};
  const double thrs[] = {10, 20, 40, 60, 80, 100, 130, 160, 200, 250, 300, 400, 1e6};
  const int nThr = sizeof(thrs) / sizeof(thrs[0]);
  const double Ws[] = {8.192, 32.768, 131.072};
  const int nW = 3;
  TCanvas canvas("c", "", 1500, 600);
  canvas.Divide(2, 1);
  std::vector<TObject *> keep;
  std::ofstream macros(macroPath);
  macros << "% Written by scripts/run70_baseline.C; edit that, not this file."
         << std::endl;
  macros << "% Run 70, right long group (10 us pole-zero on ~22 us preamps): the "
            "firmware baseline change between consecutive hits on a channel, "
            "against CAEN's description of LF limitation On. Thresholds in "
            "trapezoid-probe units (~ADC). Chi2 is the best over the baseline "
            "window lengths 8.2, 32.8, 131.1 us."
         << std::endl;
  for (int g = 0; g < 2; g++) {
    std::vector<double> shape(NDS, 0);
    double nShape = 0;
    std::vector<Pair> pairs;
    for (auto &kv : byChannel) {
      if ((g == 0) != IsRight(kv.first))
        continue;
      const std::vector<Hit> &v = kv.second;
      for (size_t j = 1; j + 1 < v.size(); j++) {
        if (!v[j].traced || v[j].t - v[j - 1].t < 300000)
          continue;
        // The isolated shape from hits with nothing for 40 us after them.
        if (v[j + 1].t - v[j].t >= 40000 && v[j].e > 800) {
          for (int q = 0; q < NDS; q++)
            shape[q] += v[j].ds1[q] / v[j].e;
          nShape++;
        }
        // A pair with the next hit, with no third hit within 27 us after it.
        const double dt = (v[j + 1].t - v[j].t) * 1e-3;
        if (j + 2 < v.size() && v[j + 2].t - v[j + 1].t < 27000)
          continue;
        double B2;
        if (v[j + 1].traced)
          B2 = v[j + 1].B;
        else if (dt >= 7 && dt < 26.5) {
          const double P2 = Coarse(v[j], sRead[g] + dt / 0.008);
          if (std::isnan(P2))
            continue;
          B2 = v[j].pre1 + P2 - v[j + 1].e / kE[g];
        } else
          continue;
        if (dt < 7 || dt > 400)
          continue;
        pairs.push_back({dt, double(v[j].e), B2 - v[j].B});
      }
    }
    for (auto &x : shape)
      x /= nShape;
    // The start of the trapezoid's rise in record time: 2 % of its peak.
    const double peak = *std::max_element(shape.begin(), shape.end());
    double tStart = 0;
    for (int q = 0; q < NDS; q++)
      if (shape[q] > 0.02 * peak) {
        tStart = 0.2 * q + 0.1;
        break;
      }
    TGraph graphShape;
    for (int q = 0; q < NDS; q++)
      graphShape.SetPoint(q, 0.2 * q + 0.1, shape[q]);
    TF1 tail("tail", "[0]*exp(-(x-18)/[1])", 18, 31.8);
    tail.SetParameters(shape[90], 20);
    graphShape.Fit(&tail, "QR");
    // Per unit energy, at record time tr.
    auto S = [&](double tr) {
      if (tr < 0.1)
        return shape[0];
      if (tr < 31.8)
        return graphShape.Eval(tr);
      return tail.Eval(tr);
    };
    // The trapezoid where its fall ends, per unit energy: the first 0.2 us
    // step after the steep fall began that drops by under 1 % of the peak.
    int qPeak = int(std::max_element(shape.begin(), shape.end()) - shape.begin());
    int qEnd = NDS - 1;
    bool falling = false;
    for (int q = qPeak + 1; q < NDS; q++) {
      const double drop = shape[q - 1] - shape[q];
      if (!falling && drop > 0.05 * peak)
        falling = true;
      else if (falling && drop < 0.01 * peak) {
        qEnd = q;
        break;
      }
    }
    const double endOfFall = shape[qEnd];
    printf("%s: %zu pairs; trapezoid at the end of its fall (%.1f us after "
           "the rise starts) %.3f of e; tail time constant %.1f us\n",
           groupName[g], pairs.size(), 0.2 * qEnd + 0.1 - tStart, endOfFall,
           tail.GetParameter(1));

    std::vector<double> mean(NB), error(NB), count(NB), centre(NB);
    for (int b = 0; b < NB; b++) {
      double s = 0, s2 = 0, k = 0, sx = 0;
      for (auto &p : pairs)
        if (p.dt >= edges[b] && p.dt < edges[b + 1]) {
          s += p.dB;
          s2 += p.dB * p.dB;
          k++;
          sx += p.dt;
        }
      count[b] = k;
      mean[b] = k ? s / k : 0;
      error[b] = k > 1 ? TMath::Sqrt((s2 / k - mean[b] * mean[b]) / k) : 1e9;
      centre[b] = k ? sx / k : 0;
    }
    auto predict = [&](double thr, double Wus, double dt, double e1) {
      const double step = 0.04; // us
      double tRun = -1;
      for (int m = 1; m < 200; m++) {
        const double tc = 7.0 * m;
        if (tc >= dt)
          break;
        if (TMath::Abs(e1 * S(tStart + tc)) < thr) {
          tRun = tc;
          break;
        }
      }
      if (tRun < 0)
        return 0.0;
      // The mean of the last W: the samples since the restart for the latest
      // part of it, the old baseline (B1, the resting level) for the rest.
      const double used = TMath::Min(dt - tRun, Wus);
      double sumLast = 0;
      for (double x = dt - used; x < dt; x += step)
        sumLast += e1 * S(tStart + x);
      return sumLast * step / Wus;
    };
    auto chi2 = [&](double thr, double W) {
      double chi = 0;
      for (int b = 0; b < NB; b++) {
        if (count[b] < 20)
          continue;
        double s = 0, k = 0;
        for (auto &p : pairs)
          if (p.dt >= edges[b] && p.dt < edges[b + 1]) {
            s += predict(thr, W, p.dt, p.e1);
            k++;
          }
        chi += TMath::Power((mean[b] - s / k) / error[b], 2);
      }
      return chi;
    };
    int nBins = 0;
    for (int b = 0; b < NB; b++)
      if (count[b] >= 20)
        nBins++;
    printf("  chi2 over %d gap bins (7-400 us), threshold x window:\n", nBins);
    std::vector<double> bestOverW(nThr, 1e99);
    double bestChi = 1e99, bestThr = 0, bestW = 0;
    for (int i = 0; i < nThr; i++) {
      printf("  thr %7.0f:", thrs[i]);
      for (int w = 0; w < nW; w++) {
        const double chi = chi2(thrs[i], Ws[w]);
        printf("  W %6.1f us %9.1f", Ws[w], chi);
        bestOverW[i] = TMath::Min(bestOverW[i], chi);
        if (chi < bestChi) {
          bestChi = chi;
          bestThr = thrs[i];
          bestW = Ws[w];
        }
      }
      printf("\n");
    }
    if (g == 0) {
      double maxGap = 0;
      for (int b = 0; b < NB; b++)
        if (count[b] >= 20)
          maxGap = TMath::Max(maxGap, TMath::Abs(mean[b]) + error[b]);
      auto at = [&](double thr) {
        for (int i = 0; i < nThr; i++)
          if (thrs[i] == thr)
            return bestOverW[i];
        return -1.0;
      };
      macros << Form("\\newcommand{\\RunSeventyPairs}{%zu}", pairs.size())
             << std::endl;
      macros << Form("\\newcommand{\\RunSeventyBins}{%d}", nBins) << std::endl;
      macros << Form("\\newcommand{\\RunSeventyOvershoot}{%.0f}",
                     100.0 * endOfFall)
             << std::endl;
      macros << Form("\\newcommand{\\RunSeventyMaxChange}{%.0f}", maxGap)
             << std::endl;
      macros << Form("\\newcommand{\\RunSeventyChiTen}{%.0f}", at(10))
             << std::endl;
      macros << Form("\\newcommand{\\RunSeventyChiTwenty}{%.0f}", at(20))
             << std::endl;
      macros << Form("\\newcommand{\\RunSeventyChiEighty}{%.0f}", at(80))
             << std::endl;
      macros << Form("\\newcommand{\\RunSeventyChiRestart}{%.0f}", at(1e6))
             << std::endl;
    }

    TGraphErrors *data = new TGraphErrors();
    TGraph *model = new TGraph(), *restart = new TGraph(), *restartShort = new TGraph();
    printf("  data against the best model (thr %.0f, W %.1f us):\n", bestThr, bestW);
    for (int b = 0; b < NB; b++) {
      if (count[b] < 20)
        continue;
      double s = 0, s2 = 0, s3 = 0, k = 0;
      for (auto &p : pairs)
        if (p.dt >= edges[b] && p.dt < edges[b + 1]) {
          s += predict(bestThr, bestW, p.dt, p.e1);
          s2 += predict(1e6, bestW, p.dt, p.e1);
          s3 += predict(1e6, 8.192, p.dt, p.e1);
          k++;
        }
      printf("   %6.1f-%-6.1f n %5.0f  dB %7.1f +- %4.1f  model %7.1f  restart "
             "after 7 us: %7.1f (W %.0f us), %7.1f (W 8.2 us)\n",
             edges[b], edges[b + 1], count[b], mean[b], error[b], s / k,
             s2 / k, bestW, s3 / k);
      const int i = data->GetN();
      data->SetPoint(i, centre[b], mean[b]);
      data->SetPointError(i, 0, error[b]);
      model->SetPoint(i, centre[b], s / k);
      restart->SetPoint(i, centre[b], s2 / k);
      restartShort->SetPoint(i, centre[b], s3 / k);
    }
    canvas.cd(g + 1)->SetLogx();
    data->SetTitle(Form("%s: firmware baseline change to the next pulse;time to "
                        "the next pulse [#mus];#DeltaB [probe units]",
                        groupName[g]));
    data->SetMarkerStyle(20);
    data->Draw("AP");
    data->GetYaxis()->SetRangeUser(-20, g == 0 ? 120 : 60);
    model->SetLineColor(kRed);
    model->SetLineWidth(2);
    model->Draw("L SAME");
    restart->SetLineColor(kBlue);
    restart->SetLineStyle(2);
    restart->Draw("L SAME");
    restartShort->SetLineColor(kGreen + 2);
    restartShort->SetLineStyle(3);
    restartShort->Draw("L SAME");
    TLegend *legend = new TLegend(0.45, 0.68, 0.89, 0.89);
    legend->AddEntry(data, "run 70", "p");
    legend->AddEntry(model, Form("CAEN's description, thr %.0f, W %.0f #mus", bestThr, bestW), "l");
    legend->AddEntry(restart, Form("restart after 7 #mus, W %.0f #mus", bestW), "l");
    legend->AddEntry(restartShort, "restart after 7 #mus, W 8.2 #mus", "l");
    legend->Draw();
    keep.insert(keep.end(), {data, model, restart, restartShort, legend});
  }
  canvas.SaveAs(plotPath);
  printf("wrote %s and %s\n", macroPath, plotPath);
}
