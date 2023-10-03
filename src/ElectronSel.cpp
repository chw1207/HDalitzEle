#include <vector>
#include <cmath>
#include <algorithm>
#include "ElectronSel.h"
#include "Utilities.h"

double eleSel::MVARawCutWP90(float pt, int eleCase){
    switch (eleCase){
        case 0: // EB1_5
            return 2.84704783417 - std::exp(-pt / 3.32529515837) * 9.38050947827;
        case 1: // EB2_5
            return 2.84704783417 - std::exp(-pt / 3.32529515837) * 9.38050947827;
        case 2: // EE_5
            return 2.84704783417 - std::exp(-pt / 3.32529515837) * 9.38050947827;
        case 3: // EB1_10
            return 2.84704783417 - std::exp(-pt / 3.32529515837) * 9.38050947827;
        case 4: // EB2_10
            return 2.84704783417 - std::exp(-pt / 3.32529515837) * 9.38050947827;
        case 5: // EE_10
            return 2.84704783417 - std::exp(-pt / 3.32529515837) * 9.38050947827;
    }
    return 999.;
}


ROOT::RVec<int> eleSel::Fall17V2ID(
    const int nEle,
    const ROOT::RVec<float>& elePt,
    const ROOT::RVec<float>& eleSCEta,
    const ROOT::RVec<float>& eleIDMVAIso
){
    ROOT::RVec<int> v;
    v.clear();

    for (int i = 0; i < nEle; i++){
        int eleCase = -1;
        if (elePt[i] < 10 && abs(eleSCEta[i]) < 0.8)
            eleCase = 0;
        else if (elePt[i] < 10 && abs(eleSCEta[i]) >= 0.8 && abs(eleSCEta[i]) < 1.479)
            eleCase = 1;
        else if (elePt[i] < 10 && abs(eleSCEta[i]) >= 1.479)
            eleCase = 2;
        else if (elePt[i] >= 10 && abs(eleSCEta[i]) < 0.8)
            eleCase = 3;
        else if (elePt[i] >= 10 && abs(eleSCEta[i]) >= 0.8 && abs(eleSCEta[i]) < 1.479)
            eleCase = 4;
        else if (elePt[i] >= 10 && abs(eleSCEta[i]) >= 1.479)
            eleCase = 5;

        double rawcut = MVARawCutWP90(elePt[i], eleCase);

        // The raw cut of Fall17 MVA ID need to convert to the MVA score between -1 and 1!
        // https://github.com/cms-sw/cmssw/blob/master/RecoEgamma/ElectronIdentification/src/ElectronMVAEstimatorRun2.cc#L116
        // https://github.com/cms-sw/cmssw/blob/master/CondFormats/GBRForest/interface/GBRForest.h#L57-L60
        const double mvacut = 2.0 / (1.0 + std::exp(-2.0 * rawcut)) - 1;
        if (eleIDMVAIso[i] > mvacut)
            v.push_back(1);
        else
            v.push_back(0);
    }

    return v;
}


ROOT::RVec<int> eleSel::HggPresel(
    const int nEle,
    const ROOT::RVec<float>& eleSCEta,
    const ROOT::RVec<float>& eleSCPhi,
    const int nPho,
    const float rhoAll,
    const ROOT::RVec<float>& phoSCEta,
    const ROOT::RVec<float>& phoSCPhi,
    const ROOT::RVec<float>& phoPFChIso,
    const ROOT::RVec<float>& phoPFPhoIso,
    const ROOT::RVec<float>& phoTrkIsoHollowConeDR03,
    const ROOT::RVec<float>& phoR9Full5x5,
    // const ROOT::RVec<float>& phoCalibEt,
    const ROOT::RVec<float>& phoEt,
    const ROOT::RVec<float>& phoSigmaIEtaIEtaFull5x5,
    const ROOT::RVec<float>& phoHoverE
){
    ROOT::RVec<int> v(nEle);

    const std::vector<float> bins = {0., 0.9, 1.5, 2, 2.2, 3};
    const std::vector<float> reff = {0.16544, 0.16544, 0.13212, 0.13212, 0.13212};
    for (int i = 0; i < nEle; i++){
        // find the phton matched to the electron
        int phoIdx = -1;
        for (int j = 0; j < nPho; j++){
            if ((phoSCEta[j] == eleSCEta[i]) && (phoSCPhi[j] == eleSCPhi[i])){
                phoIdx = j;
                break;
            }
        }

        int isHgg = 0;
        if (phoIdx != -1){ // prevent nPho < 1
            // to mitigate PU effect
            const int effbin = utils::FindBins(bins, fabs(phoSCEta[phoIdx]));
            const float phoEffArea = (effbin != -1) ? reff[effbin] : 0.;
            const float phoPFPhoIso_corr = std::max(phoPFPhoIso[phoIdx] - rhoAll * phoEffArea, (float) 0.);

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
            isHgg = (isAOD && (isHR9_EB || isLR9_EB || isHR9_EE || isLR9_EE)) ? 1 : 0;
        }
        v[i] = isHgg;
    }

    return v;
}


ROOT::RVec<int> eleSel::LateConvVeto(
    const ROOT::RVec<float>& eleSCEta,
    const ROOT::RVec<float>& eleSCPhi,
    
    const int nConv,
    const ROOT::RVec<float>& convFitPairPX,
    const ROOT::RVec<float>& convFitPairPY,
    const ROOT::RVec<float>& convFitPairPZ,
    const ROOT::RVec<float>& convVtxRadius
){
    ROOT::RVec<int> v;
    for (size_t i = 0; i < eleSCEta.size(); i++){
        int pass_ = 1;
        if (nConv > 0){
            int selected_conversion_index = -1;
            for (int j = 0; j < nConv; j++){
                ROOT::Math::XYZVector RefPairMo(convFitPairPX[j], convFitPairPY[j], convFitPairPZ[j]);

                float dR = ROOT::VecOps::DeltaR(eleSCEta[i], (float) RefPairMo.Eta(), eleSCPhi[i], (float) RefPairMo.Phi());
                float dEta = RefPairMo.Eta() - eleSCEta[i];
                float dPhi = ROOT::VecOps::DeltaPhi(eleSCPhi[i], (float) RefPairMo.Phi());
                if (dR > 0.1)
                    continue;
                if (dEta > 999.)
                    continue;
                if (dPhi > 999.)
                    continue;

                selected_conversion_index = j;
                break;
            }   
            if (selected_conversion_index != -1 && convVtxRadius[selected_conversion_index] > 40) // convVtxRadius stored in ggNtuple is radius^2
                pass_ = 0; // if matched vtx has radius > sqrt{40} cm then veto the electron
        }
        v.emplace_back(pass_);
    }

    return v;
}



