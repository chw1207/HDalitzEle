// H->gg preselection: https://arxiv.org/pdf/1804.02716.pdf
// https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/Taggers/python/flashggPreselectedDiPhotons_cfi.py
// require the photon which is matched to the electron by SC to pass H->gg preselection
bool preselection(TreeReader &data, int ele_ind){
    Int_t nPho                          = data.GetInt("nPho");
    Int_t* phoEleVeto                   = data.GetPtrInt("phoEleVeto");
    Float_t* phoSCEta                   = data.GetPtrFloat("phoSCEta");
    Float_t* phoSCPhi                   = data.GetPtrFloat("phoSCPhi");
    Float_t* phoPFChIso                 = data.GetPtrFloat("phoPFChIso");
    Float_t* phoPFPhoIso                = data.GetPtrFloat("phoPFPhoIso");
    Float_t* phoCalibEt                 = data.GetPtrFloat("phoCalibEt");
    Float_t* phoSigmaIEtaIEtaFull5x5    = data.GetPtrFloat("phoSigmaIEtaIEtaFull5x5");
    Float_t* phoHoverE                  = data.GetPtrFloat("phoHoverE");
    Float_t* phoTrkIsoHollowConeDR03    = data.GetPtrFloat("phoTrkIsoHollowConeDR03"); //![FIXED ME]
    Float_t* phoR9Full5x5               = data.GetPtrFloat("phoR9Full5x5");

    Float_t rhoAll                      = data.GetFloat("rhoAll");//![FIXED ME]
    Float_t* eleSCEta                   = data.GetPtrFloat("eleSCEta");
    Float_t* eleSCPhi                   = data.GetPtrFloat("eleSCPhi");

    //Run1: https://github.com/andreypz/nwu-dalitz-analysis/blob/master/zgamma/egamma.C#L493
    int pho_ind = -1;
    for (int i = 0; i < nPho; i++){
        if (fabs(phoSCEta[i] - eleSCEta[ele_ind]) < 0.0001){
            pho_ind = i;
            break;
        }
    }

    bool ispre = false;

    if (pho_ind == -1){
        return ispre; // no matched photon
    }

    // rho correction for isolation variables
    float corr = 0; 
    // float rhoAll = 0.;
    if (phoSCEta[pho_ind] > 0. && phoSCEta[pho_ind] <= 1.5){
        corr = 0.16544;
    }
    else if (phoSCEta[pho_ind] > 1.5 && phoSCEta[pho_ind] <= 3.){
        corr = 0.13212;
    }
    float phoPFChIso_corr = phoPFChIso[pho_ind] - rhoAll * corr;
    float phoPFPhoIso_corr = phoPFPhoIso[pho_ind] - rhoAll * corr;
    float phoTrkIsoHollowConeDR03_corr = phoTrkIsoHollowConeDR03[pho_ind] - rhoAll * corr;

    // acceptance
    bool isEB = (fabs(phoSCEta[pho_ind]) < 1.4442);
    bool isEE = (fabs(phoSCEta[pho_ind]) > 1.566 && fabs(phoSCEta[pho_ind]) < 2.5);
    
    // pass electron conversion safe veto
    bool passVeto = phoEleVeto[pho_ind] == 1;

    // HLT trail Pt threshold
    bool acc_Pt = (phoCalibEt[pho_ind] > 25.);
    
    // To mimic the MiniAOD photon cuts
    bool isAOD = (phoR9Full5x5[pho_ind] > 0.8 && phoPFChIso_corr < 20.) || (phoPFChIso_corr/phoCalibEt[pho_ind] < 0.3);

    // region
    bool isHR9_EB = phoR9Full5x5[pho_ind] > 0.85 && isEB;
    bool isLR9_EB = phoR9Full5x5[pho_ind] <= 0.85 && phoR9Full5x5[pho_ind] >= 0.5 && isEB;
    bool isHR9_EE = phoR9Full5x5[pho_ind] > 0.9 && isEE;
    bool isLR9_EE = phoR9Full5x5[pho_ind] <= 0.9 && phoR9Full5x5[pho_ind] >= 0.8 && isEE; 


    if (isAOD){
        if (isLR9_EB){
            if (phoHoverE[pho_ind] < 0.08 && phoSigmaIEtaIEtaFull5x5[pho_ind] < 0.015 && phoPFPhoIso_corr < 4. && phoTrkIsoHollowConeDR03_corr < 6.){
                // && phoTrkIsoHollowConeDR03_corr < 6.
                ispre = true;
            }
        }
        if (isHR9_EB){
            if (phoHoverE[pho_ind] < 0.08) ispre = true;
        }
        if (isLR9_EE){
            if (phoHoverE[pho_ind] < 0.08 && phoSigmaIEtaIEtaFull5x5[pho_ind] < 0.035 && phoPFPhoIso_corr < 4. && phoTrkIsoHollowConeDR03_corr < 6.){
                // && phoTrkIsoHollowConeDR03_corr < 6.
                ispre = true;
            }
        }
        if (isHR9_EE){
            if (phoHoverE[pho_ind] < 0.08) ispre = true;
        }
    }
    
    return ispre;    
}