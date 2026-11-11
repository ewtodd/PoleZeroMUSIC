// Long-end deposit pairs of reacted 37Cl + 4He particles, from MUSIC's 37Cl
// simulation, for simulate --reactionPairs. Every (alpha,n) and (alpha,p)
// event reacting in strips 2-8 gives one pair per long-left / long-right
// neighbour pair (odd strip k, even strip k+1) at or after its reaction strip,
// each end in units of the unreacted beam's mean deposit there. A pair's
// weight is its channel's TALYS cross section at the strip's beam energy,
// shared over that strip's events and pairs, so channels and strips enter at
// their relative rates. The file ends with the physical rate of these
// reactions per beam particle, a "# rate" comment line that figures reads for
// the reaction table's caption: the summed cross sections times the He areal
// density of one strip.
//   root -l -b -q 'scripts/reaction_pairs.C("<sim_root_files>", "<pairs.txt>")'
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

const Int_t kNStrips = 16;
const Long64_t kEventsPerFile = 3000;
const Double_t kMassAlpha = 4.0026;
const Double_t kMassCl37 = 36.9659;
// The He areal density of one strip, as MUSIC's CrossSection.cpp computes it
// (ideal gas at the fill temperature over one strip's active length) at the
// 37Cl runs' pressure (CROSS_SECTION_CONFIG.GAS_PRESSURE_TORR).
const Double_t kGasPressureTorr = 400.0;
const Double_t kGasTemperatureK = 293.0;
const Double_t kStripLengthCm = 1.578;
const Double_t kBoltzmann = 1.380649e-23; // J/K
const Double_t kPaPerTorr = 133.322368;
const Double_t kCm2PerMb = 1.0e-27;

// TALYS cross section [mb] of one exit channel ("n", "p"), from the table in
// dir whose centre-of-mass energy is closest to ecm and which lists it.
Double_t TalysSigma(const TString &dir, const TString &channel, Double_t ecm) {
  void *listing = gSystem->OpenDirectory(dir);
  if (!listing)
    return -1.0;
  Double_t best = -1.0, bestDistance = 1.0e9;
  const char *entry;
  while ((entry = gSystem->GetDirEntry(listing))) {
    const TString name = entry;
    const Ssiz_t at = name.Index("Ecm");
    if (!name.EndsWith(".tab") || at < 0)
      continue;
    const Double_t tableEcm = TString(name(at + 3, name.Length())).Atof();
    std::ifstream in((dir + "/" + name).Data());
    std::string line;
    Double_t sigma = 0.0;
    Bool_t listed = kFALSE;
    while (std::getline(in, line)) {
      std::istringstream fields(line);
      std::string key, which;
      Double_t value = 0.0;
      if (!(fields >> key >> which >> value))
        continue;
      if (key == "channel" && which == channel.Data()) {
        sigma += value;
        listed = kTRUE;
      }
    }
    if (listed && TMath::Abs(tableEcm - ecm) < bestDistance) {
      bestDistance = TMath::Abs(tableEcm - ecm);
      best = sigma;
    }
  }
  gSystem->FreeDirectory(listing);
  return best;
}

void reaction_pairs(const char *simDir, const char *outPath) {
  const TString dir = simDir;
  Float_t left[kNStrips], right[kNStrips];
  Double_t beamL[kNStrips] = {0.0}, beamR[kNStrips] = {0.0};
  {
    TFile file(dir + "/traces_37Cl_beam.root");
    TTree *tree = (TTree *)file.Get("events_MeV");
    if (!tree) {
      std::cerr << "no beam traces in " << dir << std::endl;
      return;
    }
    tree->SetBranchAddress("LeftdE", left);
    tree->SetBranchAddress("RightdE", right);
    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
      tree->GetEntry(i);
      for (Int_t k = 0; k < kNStrips; k++) {
        beamL[k] += left[k];
        beamR[k] += right[k];
      }
    }
    for (Int_t k = 0; k < kNStrips; k++) {
      beamL[k] /= Double_t(tree->GetEntries());
      beamR[k] /= Double_t(tree->GetEntries());
    }
  }

  std::ofstream out(outPath);
  out << "# Long-end deposit pairs of reacted particles, from " << dir
      << std::endl;
  out << "# (scripts/reaction_pairs.C): rL rR weight, rL the left end of an "
         "odd strip k, rR the right end of strip k+1, each over the "
         "unreacted beam's mean there; weight in mb"
      << std::endl;
  const Int_t nChannels = 2;
  const char *channels[nChannels] = {"n", "p"};
  Double_t sumW = 0.0, sumL = 0.0, sumR = 0.0, sumLL = 0.0, sumRR = 0.0,
           sumLR = 0.0;
  // Each channel's cross section at each reaction strip, summed: the rate per
  // beam particle per unit areal density.
  Double_t sumSigma = 0.0;
  Long64_t rows = 0;
  for (Int_t c = 0; c < nChannels; c++) {
    for (Int_t s = 2; s <= 8; s++) {
      TFile file(Form("%s/traces_37Cl_a%s_s%d.root", simDir, channels[c], s));
      TTree *tree = (TTree *)file.Get("events_MeV");
      TTree *mc = (TTree *)file.Get("MC");
      if (!tree || !mc) {
        std::cerr << "no traces for (a," << channels[c] << ") in strip " << s
                  << std::endl;
        continue;
      }
      Int_t reactionStrip = -1;
      Float_t beamEnergy = 0.0;
      mc->SetBranchAddress("reaction_strip", &reactionStrip);
      mc->SetBranchAddress("beam_energy_reaction", &beamEnergy);
      std::vector<Long64_t> used;
      Double_t sumElab = 0.0;
      for (Long64_t i = 0;
           i < mc->GetEntries() && Long64_t(used.size()) < kEventsPerFile;
           i++) {
        mc->GetEntry(i);
        if (reactionStrip != s)
          continue;
        used.push_back(i);
        sumElab += beamEnergy;
      }
      if (used.empty())
        continue;
      const Double_t ecm = sumElab / Double_t(used.size()) * kMassAlpha /
                           (kMassCl37 + kMassAlpha);
      const Double_t sigma = TalysSigma(dir + "/talys", channels[c], ecm);
      if (!(sigma > 0.0)) {
        std::cerr << "no TALYS (a," << channels[c] << ") at Ecm " << ecm
                  << std::endl;
        continue;
      }
      sumSigma += sigma;
      Int_t nPairs = 0;
      for (Int_t k = 1; k + 1 <= kNStrips; k += 2)
        if (k >= s)
          nPairs++;
      const Double_t weight = sigma / (Double_t(used.size()) * nPairs);
      tree->SetBranchAddress("LeftdE", left);
      tree->SetBranchAddress("RightdE", right);
      for (size_t e = 0; e < used.size(); e++) {
        tree->GetEntry(used[e]);
        for (Int_t k = 1; k + 1 <= kNStrips; k += 2) {
          if (k < s)
            continue;
          const Double_t rL = left[k - 1] / beamL[k - 1];
          const Double_t rR = right[k] / beamR[k];
          out << Form("%.5f %.5f %.6g", rL, rR, weight) << std::endl;
          sumW += weight;
          sumL += weight * rL;
          sumR += weight * rR;
          sumLL += weight * rL * rL;
          sumRR += weight * rR * rR;
          sumLR += weight * rL * rR;
          rows++;
        }
      }
      std::cout << Form("(a,%s) strip %d: Ecm %.3f MeV, sigma %.3g mb, %zu "
                        "events x %d pairs",
                        channels[c], s, ecm, sigma, used.size(), nPairs)
                << std::endl;
    }
  }
  const Double_t mL = sumL / sumW, mR = sumR / sumW;
  const Double_t sL = TMath::Sqrt(sumLL / sumW - mL * mL);
  const Double_t sR = TMath::Sqrt(sumRR / sumW - mR * mR);
  const Double_t areal = kGasPressureTorr * kPaPerTorr /
                         (kBoltzmann * kGasTemperatureK) * 1.0e-6 *
                         kStripLengthCm;
  const Double_t rate = sumSigma * kCm2PerMb * areal;
  out << Form("# rate %.4g per beam particle: (alpha,n) and (alpha,p) in "
              "strips 2-8, %.4g mb summed, %.4g He atoms/cm^2 per strip",
              rate, sumSigma, areal)
      << std::endl;
  std::cout << Form("wrote %lld pairs to %s: rL %.3f +- %.3f, rR %.3f +- "
                    "%.3f, correlation %+.2f; physical rate %.3g per beam "
                    "particle",
                    rows, outPath, mL, sL, mR, sR,
                    (sumLR / sumW - mL * mR) / (sL * sR), rate)
            << std::endl;
}
