using VecF_t = const ROOT::RVec<float>&;
using VecI_t = const ROOT::RVec<int>&;
auto HggPreSelection(float rhoAll, int nEle, VecF_t eleSCEta, int nPho, VecF_t phoSCEta, VecF_t phoPFChIso, VecF_t phoPFPhoIso, VecF_t phoTrkIsoHollowConeDR03, VecF_t phoR9Full5x5, VecF_t phoCalibEt, VecF_t phoSigmaIEtaIEtaFull5x5, VecF_t phoHoverE){

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