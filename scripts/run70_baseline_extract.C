// Step 1 of scripts/run70_baseline.sh: one pass over run 70's converted
// traces (data_raw/run70.root, probe 0 ADCInput, probe 1 EnergyFilter, the
// trapezoid before baseline subtraction) keeping every long-channel hit with
// its time and energy, and for hits with a 4000-sample record a compact copy
// of the trapezoid (trace1) and the input (trace0):
//   25-sample means (0.2 us) of each relative to its first samples, the
//   full-resolution trapezoid around the peaking point, and the 50 % crossing
//   of the trapezoid's rising edge.
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
const int MAXM = 100, MAXL = 5000, NDS = 160, W0 = 950, NW = 300;
static int tr0[MAXM][MAXL], tr1[MAXM][MAXL];
void run70_baseline_extract(const char *in, const char *out, long nEv = -1) {
  TFile f(in);
  TTree *t = (TTree *)f.Get("tree");
  unsigned int multi; unsigned short ch[MAXM], e[MAXM]; unsigned long long et[MAXM]; int len[MAXM];
  t->SetBranchStatus("*", 0);
  for (auto b : {"multi","ch","e","e_t","traceLen","trace0","trace1"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("multi", &multi); t->SetBranchAddress("ch", ch); t->SetBranchAddress("e", e);
  t->SetBranchAddress("e_t", et); t->SetBranchAddress("traceLen", len);
  t->SetBranchAddress("trace0", tr0); t->SetBranchAddress("trace1", tr1);
  TFile o(out, "RECREATE");
  TTree h("hits", "long-channel hits");
  int c, ev, hasTrace, t50; unsigned short en; unsigned long long ts; float pre1, pre0, ds1[NDS], ds0[NDS]; short win[NW];
  h.Branch("ch", &c); h.Branch("ev", &ev); h.Branch("e", &en); h.Branch("t", &ts); h.Branch("hasTrace", &hasTrace);
  h.Branch("pre1", &pre1); h.Branch("pre0", &pre0); h.Branch("t50", &t50);
  h.Branch("ds1", ds1, Form("ds1[%d]/F", NDS)); h.Branch("ds0", ds0, Form("ds0[%d]/F", NDS));
  h.Branch("win", win, Form("win[%d]/S", NW));
  const long n = nEv < 0 ? t->GetEntries() : nEv;
  for (long i = 0; i < n; i++) {
    t->GetEntry(i);
    if (i % 2000 == 0) { printf("event %ld / %ld\n", i, n); fflush(stdout); }
    for (unsigned k = 0; k < multi; k++) {
      int cc = ch[k];
      bool isLong = (cc < 64) && ((cc % 4 == 0) || (cc % 8 == 1));  // R even at 4n, L odd at 8n+1
      if (!isLong) continue;
      c = cc; ev = int(i); en = e[k]; ts = et[k]; hasTrace = (len[k] == 4000);
      pre1 = pre0 = 0; t50 = -1;
      for (int j = 0; j < NDS; j++) ds1[j] = ds0[j] = 0;
      for (int j = 0; j < NW; j++) win[j] = 0;
      if (hasTrace) {
        for (int s = 300; s < 600; s++) pre1 += tr1[k][s]; pre1 /= 300;
        for (int s = 0; s < 25; s++) pre0 += tr0[k][s]; pre0 /= 25;
        for (int j = 0; j < NDS; j++) {
          double a = 0, b = 0;
          for (int s = 25 * j; s < 25 * j + 25; s++) { a += tr1[k][s]; b += tr0[k][s]; }
          ds1[j] = a / 25 - pre1; ds0[j] = b / 25 - pre0;
        }
        for (int j = 0; j < NW; j++) win[j] = short(tr1[k][W0 + j] - pre1);
        int imax = 600; for (int s = 600; s < 1400; s++) if (tr1[k][s] > tr1[k][imax]) imax = s;
        double half = 0.5 * (tr1[k][imax] - pre1) + pre1;
        for (int s = 600; s < imax; s++) if (tr1[k][s] >= half) { t50 = s; break; }
      }
      h.Fill();
    }
  }
  o.cd(); h.Write(); o.Close();
  printf("done\n");
}
