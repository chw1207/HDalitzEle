#include "PhotonSel.h"
#include "Utilities.h"
#include "TMath.h"

ROOT::RVec<int> phoSel::HggPresel(
    const int nPho,
    const float rhoAll,
    const ROOT::RVec<float>& phoSCEta,
    const ROOT::RVec<float>& phoPFChIso,
    const ROOT::RVec<float>& phoPFPhoIso,
    const ROOT::RVec<float>& phoTrkIsoHollowConeDR03,
    const ROOT::RVec<float>& phoR9Full5x5,
    // const ROOT::RVec<float>& phoCalibEt,
    const ROOT::RVec<float>& phoEt,
    const ROOT::RVec<float>& phoSigmaIEtaIEtaFull5x5,
    const ROOT::RVec<float>& phoHoverE
){
    ROOT::RVec<int> v;
    v.reserve(nPho);

    const std::vector<float> bins = {0., 0.9, 1.5, 2, 2.2, 3};
    const std::vector<float> reff = {0.16544, 0.16544, 0.13212, 0.13212, 0.13212};
    for (int phoIdx = 0; phoIdx < nPho; phoIdx++){
        // to mitigate PU effect
        const int effbin = utils::FindBins(bins, fabs(phoSCEta[phoIdx]));
        const float phoEffArea = (effbin != -1) ? reff[effbin] : 0.;
        const float phoPFPhoIso_corr = TMath::Max(phoPFPhoIso[phoIdx] - rhoAll * phoEffArea, (float) 0.);

        const bool isEB = (fabs(phoSCEta[phoIdx]) < 1.4442);
        const bool isEE = (fabs(phoSCEta[phoIdx]) > 1.566 && fabs(phoSCEta[phoIdx]) < 2.5);

        const bool isHR9_EB = isEB && phoR9Full5x5[phoIdx] > 0.85;
        const bool isLR9_EB = isEB && phoR9Full5x5[phoIdx] <= 0.85 && phoR9Full5x5[phoIdx] > 0.5 && phoPFPhoIso_corr < 4 && phoTrkIsoHollowConeDR03[phoIdx] < 6 && phoSigmaIEtaIEtaFull5x5[phoIdx] < 0.015;
        const bool isHR9_EE = isEE && phoR9Full5x5[phoIdx] > 0.9;
        const bool isLR9_EE = isEE && phoR9Full5x5[phoIdx] <= 0.9 && phoR9Full5x5[phoIdx] > 0.8 && phoPFPhoIso_corr < 4 && phoTrkIsoHollowConeDR03[phoIdx] < 6 && phoSigmaIEtaIEtaFull5x5[phoIdx] < 0.035;

        // cuts here mimic the miniAOD photon cuts and the non-category based trigger cuts
        const bool isAOD = phoHoverE[phoIdx] < 0.08 && (phoR9Full5x5[phoIdx] > 0.8 || phoPFChIso[phoIdx] < 20. || (phoPFChIso[phoIdx]/phoEt[phoIdx] < 0.3));
        // const bool base_pT_cut = phoCalibEt[phoIdx] > 25.;

        // Hgg preselection without passElectronVeto
        // const int isHgg = (base_pT_cut && isAOD && (isHR9_EB || isLR9_EB || isHR9_EE || isLR9_EE)) ? 1 : 0;
        const int isHgg = (isAOD && (isHR9_EB || isLR9_EB || isHR9_EE || isLR9_EE)) ? 1 : 0;
        v.emplace_back(isHgg);
    }

    return v;
}