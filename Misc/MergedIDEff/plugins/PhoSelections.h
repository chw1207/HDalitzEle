#include "ROOT/RVec.hxx"
#include "TVector3.h"
#include <cmath>


ROOT::RVec<int> FSRSelec(TLorentzVector mu1, TLorentzVector mu2, ROOT::RVec<TLorentzVector> &phoP4){
    ROOT::RVec<int> v;
    v.clear();
    for (size_t i = 0; i < phoP4.size(); i++){
        int isfsr = 0;
        const bool cut1 = mu1.DeltaR(phoP4[i]) < 0.8 || mu2.DeltaR(phoP4[i]) < 0.8;
        const bool cut2 = mu1.DeltaR(phoP4[i]) > 0.1 && mu2.DeltaR(phoP4[i]) > 0.1;
        const bool cut3 = ((mu1 + mu2).M() + (mu1 + mu2 + phoP4[i]).M()) < 180;
        const bool cut4 = (mu1 + mu2 + phoP4[i]).M() > 80 && (mu1 + mu2 + phoP4[i]).M() < 100;
        if (cut1 && cut2 && cut3 && cut4) isfsr = 1;
        v.push_back(isfsr);
    }

    return v;
}


int GetZPho(TLorentzVector mu1, TLorentzVector mu2, ROOT::RVec<TLorentzVector> phoP4, ROOT::RVec<int> &isGoodPho){
    int finipho = -1;

    float mindiff = 999.;
    for (size_t i = 0; i < phoP4.size(); i++){
        if (isGoodPho[i] == 1){
            float Muu = (mu1 + mu2).M();
	        float Muug = (mu1 + mu2 + phoP4[i]).M();
	        float diff = abs(91.18 - Muug); //close to Z mass

            if (mindiff > diff){
                mindiff = diff;
		        finipho = i;
            }
        }
    }
    return finipho;
}


ROOT::RVec<int> HggPreSelection(float rhoAll, int nPho, ROOT::RVec<float>& phoSCEta, ROOT::RVec<float>& phoPFChIso, ROOT::RVec<float>& phoPFPhoIso, ROOT::RVec<float>& phoTrkIsoHollowConeDR03, ROOT::RVec<float>& phoR9Full5x5, ROOT::RVec<float>& phoCalibEt, ROOT::RVec<float>& phoSigmaIEtaIEtaFull5x5, ROOT::RVec<float>& phoHoverE){
    ROOT::RVec<int> vec;
    vec.clear();

    for (int i = 0; i < nPho; i++){
        int isHgg = 0;

        float corr = 0;
        if (phoSCEta[i] > 0. && phoSCEta[i] <= 1.5) corr = 0.16544;
        else if (phoSCEta[i] > 1.5 && phoSCEta[i] <= 3.) corr = 0.13212;

        float phoPFChIso_corr = phoPFChIso[i] - rhoAll * corr;
        float phoPFPhoIso_corr = phoPFPhoIso[i] - rhoAll * corr;
        float phoTrkIsoHollowConeDR03_corr = phoTrkIsoHollowConeDR03[i] - rhoAll * corr;

        bool isEB = (fabs(phoSCEta[i]) < 1.4442);
        bool isEE = (fabs(phoSCEta[i]) > 1.566 && fabs(phoSCEta[i]) < 2.5);

        bool isAOD = (phoR9Full5x5[i] > 0.8 && phoPFChIso_corr < 20.) || (phoPFChIso_corr/phoCalibEt[i] < 0.3);

        // region
        bool isHR9_EB = phoR9Full5x5[i] > 0.85 && isEB;
        bool isLR9_EB = phoR9Full5x5[i] <= 0.85 && phoR9Full5x5[i] >= 0.5 && isEB;
        bool isHR9_EE = phoR9Full5x5[i] > 0.9 && isEE;
        bool isLR9_EE = phoR9Full5x5[i] <= 0.9 && phoR9Full5x5[i] >= 0.8 && isEE;


        if (isAOD){
            if (isLR9_EB){
                if (phoHoverE[i] < 0.08 && phoSigmaIEtaIEtaFull5x5[i] < 0.015 && phoPFPhoIso_corr < 4. && phoTrkIsoHollowConeDR03_corr < 6.){
                    isHgg = 1;
                }
            }
            if (isHR9_EB){
                if (phoHoverE[i] < 0.08) isHgg = 1;
            }
            if (isLR9_EE){
                if (phoHoverE[i] < 0.08 && phoSigmaIEtaIEtaFull5x5[i] < 0.035 && phoPFPhoIso_corr < 4. && phoTrkIsoHollowConeDR03_corr < 6.){
                    isHgg = 1;
                }
            }
            if (isHR9_EE){
                if (phoHoverE[i] < 0.08) isHgg = 1;
            }
        }
        vec.push_back(isHgg);
    }

    return vec;
}


int GetMatchEle(float phoSCEta, ROOT::RVec<float>& eleSCEta){
    int eleidx = -1;

    for (size_t i = 0; i < eleSCEta.size(); i++){
        if (fabs(phoSCEta - eleSCEta[i]) < 0.00001){
            eleidx = i;
            break;
        }
    }

    return eleidx;
}


int ConvMatch(float phoSCEta, float phoSCPhi, float phoSCE, int nConv, ROOT::RVec<int>& convNTrks, ROOT::RVec<float>& convVtxX, ROOT::RVec<float>& convVtxY, ROOT::RVec<float>& convVtxZ, ROOT::RVec<float>& convFitPairPX, ROOT::RVec<float>& convFitPairPY, ROOT::RVec<float>& convFitPairPZ, ROOT::RVec<float>& convFitProb){
    TVector3 SC;
    SC.SetPtEtaPhi(phoSCE/cosh(phoSCEta), phoSCEta, phoSCPhi);

    float mindR = 999;
    int selected_conversion_index = -1;
    for (int i = 0; i < nConv; i++){
        if (convNTrks[i] != 2) continue;

        float pairPt = sqrt(convFitPairPX[i]*convFitPairPX[i] + convFitPairPY[i]*convFitPairPY[i]);
        if (pairPt < 10.) continue;
        if (convFitProb[i] < 1e-6) continue;

        TVector3 VtxtoSC;
        VtxtoSC.SetXYZ(SC.x() - convVtxX[i], SC.y() - convVtxY[i], SC.z() - convVtxZ[i]);

        TVector3 RefPairMo;
        RefPairMo.SetXYZ(convFitPairPX[i], convFitPairPY[i], convFitPairPZ[i]);

        float dR = VtxtoSC.DeltaR(RefPairMo);
        if(dR < mindR){
            mindR = dR;
            selected_conversion_index = i;
        }
    }

    if (mindR < 0.1)
        return selected_conversion_index;

    return -1;
}