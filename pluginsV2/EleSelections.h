#include "TMath.h"
#include "ROOT/RVec.hxx"
using namespace ROOT::VecOps;
// Function to check if gamma* will pass the Hgg preselection
//* Description:
//*     1) It requires the photon which is matched to the electron by SC to pass H->gg preselection
//*     2) H->gg preselection: https://arxiv.org/pdf/1804.02716.pdf
//*     3) https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/Taggers/python/flashggPreselectedDiPhotons_cfi.py
//*     4) Run1: https://github.com/andreypz/nwu-dalitz-analysis/blob/master/zgamma/egamma.C#L493
ROOT::RVec<int> HggPreSelection(float rhoAll, int nEle, ROOT::RVec<float> eleSCEta, int nPho, ROOT::RVec<float> phoSCEta, ROOT::RVec<float> phoPFChIso, ROOT::RVec<float> phoPFPhoIso, ROOT::RVec<float> phoTrkIsoHollowConeDR03, ROOT::RVec<float> phoR9Full5x5, ROOT::RVec<float> phoCalibEt, ROOT::RVec<float> phoSigmaIEtaIEtaFull5x5, ROOT::RVec<float> phoHoverE){

    ROOT::RVec<int> vec;
    vec.clear();
    for (int ele_ind = 0; ele_ind < nEle; ele_ind++){
        
        int isHgg = 0;
        int pho_ind = -1;
        int ele_passInd = 999;
        for (int i = 0; i < nPho; i++){
            if (fabs(phoSCEta[i] - eleSCEta[ele_ind]) < 0.0001){
                pho_ind = i;
                break;
            }
        }
        
        if (pho_ind != -1){
            float corr = 0; 
            if (phoSCEta[pho_ind] > 0. && phoSCEta[pho_ind] <= 1.5) corr = 0.16544;
            else if (phoSCEta[pho_ind] > 1.5 && phoSCEta[pho_ind] <= 3.) corr = 0.13212;

            float phoPFChIso_corr = phoPFChIso[pho_ind] - rhoAll * corr;
            float phoPFPhoIso_corr = phoPFPhoIso[pho_ind] - rhoAll * corr;
            float phoTrkIsoHollowConeDR03_corr = phoTrkIsoHollowConeDR03[pho_ind] - rhoAll * corr;

            bool isEB = (fabs(phoSCEta[pho_ind]) < 1.4442);
            bool isEE = (fabs(phoSCEta[pho_ind]) > 1.566 && fabs(phoSCEta[pho_ind]) < 2.5);

            bool isAOD = (phoR9Full5x5[pho_ind] > 0.8 && phoPFChIso_corr < 20.) || (phoPFChIso_corr/phoCalibEt[pho_ind] < 0.3);

            // region
            bool isHR9_EB = phoR9Full5x5[pho_ind] > 0.85 && isEB;
            bool isLR9_EB = phoR9Full5x5[pho_ind] <= 0.85 && phoR9Full5x5[pho_ind] >= 0.5 && isEB;
            bool isHR9_EE = phoR9Full5x5[pho_ind] > 0.9 && isEE;
            bool isLR9_EE = phoR9Full5x5[pho_ind] <= 0.9 && phoR9Full5x5[pho_ind] >= 0.8 && isEE;


            if (isAOD){
                if (isLR9_EB){
                    if (phoHoverE[pho_ind] < 0.08 && phoSigmaIEtaIEtaFull5x5[pho_ind] < 0.015 && phoPFPhoIso_corr < 4. && phoTrkIsoHollowConeDR03_corr < 6.){
                        isHgg = 1;
                    }
                }
                if (isHR9_EB){
                    if (phoHoverE[pho_ind] < 0.08) isHgg = 1;
                }
                if (isLR9_EE){
                    if (phoHoverE[pho_ind] < 0.08 && phoSigmaIEtaIEtaFull5x5[pho_ind] < 0.035 && phoPFPhoIso_corr < 4. && phoTrkIsoHollowConeDR03_corr < 6.){
                        isHgg = 1;
                    }
                }
                if (isHR9_EE){
                    if (phoHoverE[pho_ind] < 0.08) isHgg = 1;
                }
            }
        }
        
        vec.push_back(isHgg);
    }
    
    return vec;
}

// Function for official electron MVA ID (Fall17 V2 ID)
//* Description:
//*     1) 90% WP is used
//*     2) threashold values are found in: https://github.com/lsoffi/egm_tnp_analysis/blob/egm_tnp_CleanedCodeForUL_17March2020/etc/config/settings_ele_UL2017.py#L6
//*     3) SFs measurements: https://indico.cern.ch/event/944042/contributions/3966334/attachments/2084133/3501021/EGM_UL18_EleIDSF_Final.pdf
//*     4) SFs root files: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018#SFs_for_Electrons_UL_2018
ROOT::RVec<int> EleFall17V2ID(int nEle, ROOT::RVec<float> eleSCEta, ROOT::RVec<float> eleIDMVAIso){

    const vector<float> MVAcut = {0.913286, 0.805013, 0.358969};

    ROOT::RVec<int> v;
    v.clear();
    for (int i = 0; i < nEle; i++){
        bool mva_cut = false;

        if (fabs(eleSCEta[i]) < 0.8){
            if (eleIDMVAIso[i] > MVAcut[0]) mva_cut = true;
        }

        if (fabs(eleSCEta[i]) > 0.8 && fabs(eleSCEta[i]) < 1.479){
            if (eleIDMVAIso[i] > MVAcut[1]) mva_cut = true;
        }

        if (fabs(eleSCEta[i]) > 1.479){
            if (eleIDMVAIso[i] > MVAcut[2]) mva_cut = true;
        }

        if (mva_cut) v.push_back(1);
        else v.push_back(0);
    }
    return v;
}