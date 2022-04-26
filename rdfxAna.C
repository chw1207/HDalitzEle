#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "./pluginsV2/help.h"
#include "./pluginsV2/puweicalc.h"
#include "./pluginsV2/EleSelections.h"
#include "./pluginsV2/JetSelections.h"
#include "./pluginsV2/PhoR9SSCorrection.h"
#include "./pluginsV2/GenMatching.h"


// TODO:
//! [1] add uncertaities
//*     1) Hgg preselection
//! [2] some SFs need to be added
//*     3) electron conversion veto and hits
//*     4) di-photn and dielectron trigger SFs
//*     5) Hgg preselection


using namespace std;
using namespace ROOT::VecOps;


int make_cat(bool isM2, bool isM1, bool isResolved, bool isHVbf, bool isLVbf, bool isBst, bool isEE, bool isEBHR9, bool isEBLR9){
    int cat = 0;
    if (isM2){
        if (isHVbf == 1) cat = 1;
        else if (isLVbf == 1) cat = 2;
        else if (isBst == 1) cat = 3;
        else if (isEBHR9 == 1) cat = 4;
        else if (isEBLR9 == 1) cat = 5;
        else if (isEE == 1) cat = 6;
        else {
            cout << "[ERROR] No proper category for M2" << endl;
            exit(-1);
        }
    }
    else if (isM1 == 1){
        if (isHVbf == 1) cat = 7;
        else if (isLVbf == 1) cat = 8;
        else if (isBst == 1) cat = 9;
        else if (isEBHR9 == 1) cat = 10;
        else if (isEBLR9 == 1) cat = 11;
        else if (isEE == 1) cat = 12;
        else {
            cout << "[ERROR] No proper category for M1" << endl;
            exit(-1);
        }
    }
    else if (isResolved == 1) cat = 13;
    else{
        cout << "[ERROR] No proper category for Resolved" << endl;
        exit(-1);
    }
    return cat;
}


template <typename T>
auto FindGoodPho(T &df, int year) {
    const float mva_EB = (year == 2016) ? 0.2 : -0.02;
    const float mva_EE = (year == 2016) ? 0.2 : -0.26;

    auto nf = df.Define("isEBPho",      Form("abs(phoSCEta) < 1.4442 && phoIDMVA > %f", mva_EB))
                .Define("isEEPho",      Form("abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5 && phoIDMVA > %f", mva_EE))
                .Define("phoCuts",      "phoEleVeto != 0 && phoCalibEt > 15")
                .Define("isGoodPho",    "(isEBPho || isEEPho) && phoCuts")
                .Define("phoIdx1",      "if (Sum(isGoodPho) > 0) return hp::getIdx(isGoodPho, phoCalibEt)[0]; else return -1;");
    return nf;
}


template <typename T>
auto FindGoodEle(T &df){
    auto nf = df.Define("isEBEle",      "abs(eleSCEta) < 1.4442")
                .Define("isEEEle",      "abs(eleSCEta) > 1.566 && abs(eleSCEta) < 2.5")
                .Define("isHggEle",     "HggPreSelection(rhoAll, nEle, eleSCEta, nPho, phoSCEta, phoPFChIso, phoPFPhoIso, phoTrkIsoHollowConeDR03, phoR9Full5x5, phoCalibEt, phoSigmaIEtaIEtaFull5x5, phoHoverE)")
                .Define("isMVAEle",     "EleFall17V2ID(nEle, eleSCEta, eleIDMVAIso)")
                .Define("eleCuts",      "(isEBEle || isEEEle) && eleCalibPt > 10.")

                .Define("isM2GoodEle",  "eleClass == 0 && isHggEle && eleCuts && eleConvVeto == 1 && eleTrkMissHits < 2 && eleSubTrkMissHits < 2")
                .Define("isM1GoodEle",  "eleClass == 1 && isHggEle && eleCuts && eleConvVeto == 1 && eleTrkMissHits < 2")
                .Define("isReGoodEle",  "eleClass >= 2 && isMVAEle && eleCuts")

                .Define("m2Idx1",       "if (Sum(isM2GoodEle) > 0) return hp::getIdx(isM2GoodEle, eleCalibPt)[0]; else return -1;")
                .Define("m1Idx1",       "if (Sum(isM1GoodEle) > 0) return hp::getIdx(isM1GoodEle, eleCalibPt)[0]; else return -1;")
                .Define("reIdx1",       "if (Sum(isReGoodEle) > 1) return hp::getIdx(isReGoodEle, eleCalibPt)[0]; else return -1;")
                .Define("reIdx2",       "if (Sum(isReGoodEle) > 1) return hp::getIdx(isReGoodEle, eleCalibPt)[1]; else return -1;");
    return nf;
}


template <typename T>
auto FilterGoodEvents(T &df, int year){
    // HLT trigger
    //* Di-photon trigger (for merged)
    //*     1) HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v (2016)
    //*     2) HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v (2017 and 2018)
    //* Di-electron trigger (for resolved)
    //*     1) HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v (2016)
    //*     2) HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v (2017 and 2018)

    string isDiEleTrigger = " ";
    if (year == 2016) isDiEleTrigger = "(HLTEleMuX >> 40 & 1) == 1";
    else isDiEleTrigger = "((HLTEleMuX >> 40 & 1) == 1) || ((HLTEleMuX >> 5 & 1) == 1)";

    auto nf = df.Filter("Sum(isM2GoodEle) > 0 || Sum(isM1GoodEle) ||  Sum(isReGoodEle) > 1", "good ele")
                .Filter("Sum(isGoodPho) > 0", "good pho")

                .Define("isDiPhoHLT",   "((HLTPho >> 14) & 1) == 1")
                .Define("isDiEleHLT",   isDiEleTrigger)
                .Define("isM2",         "Sum(isM2GoodEle) > 0 && isDiPhoHLT && TMath::Max(eleCalibPt[m2Idx1], phoCalibEt[phoIdx1]) > 33. && TMath::Min(eleCalibPt[m2Idx1], phoCalibEt[phoIdx1]) > 25.")
                .Define("isM1",         "isM2 != 1 && Sum(isM1GoodEle) > 0 && isDiPhoHLT && TMath::Max(eleCalibPt[m1Idx1], phoCalibEt[phoIdx1]) > 33. && TMath::Min(eleCalibPt[m1Idx1], phoCalibEt[phoIdx1]) > 25.")
                .Define("isResolved",   "isM2 != 1 && isM1 != 1 && Sum(isReGoodEle) > 1 && isDiEleHLT && TMath::Max(eleCalibPt[reIdx1], eleCalibPt[reIdx2]) > 25. && TMath::Min(eleCalibPt[reIdx1], eleCalibPt[reIdx2]) > 15.")

                .Filter("isM2 || isM1 || isResolved", "HLT/PtCut")

                // kinematic selections: construct physics object
                .Define("eleIdx1",      "if (isM2) return m2Idx1; else if (isM1) return m1Idx1; else return reIdx1;")
                .Define("eleIdx2",      "if (isM2) return m2Idx1; else if (isM1) return m1Idx1; else return reIdx2;")
                .Define("ele1",         "hp::P4Mass(eleCalibPt[eleIdx1], eleEta[eleIdx1], elePhi[eleIdx1], 0.000511)")
                .Define("ele2",         "hp::P4Mass(eleCalibPt[eleIdx2], eleEta[eleIdx2], elePhi[eleIdx2], 0.000511)")
                .Define("pho",          "hp::P4Mass(phoCalibEt[phoIdx1], phoEta[phoIdx1], phoPhi[phoIdx1], 0.)")

                // only for merged electrons
                .Define("gsf1",         "hp::P4Mass(eleTrkPt[eleIdx1], eleTrkEta[eleIdx1], eleTrkPhi[eleIdx1], 0.000511)")
                .Define("gsf2",         "hp::P4Mass(eleSubTrkPt[eleIdx1], eleSubTrkEta[eleIdx1], eleSubTrkPhi[eleIdx1], 0.000511)")

                .Define("diEle",        "if (isResolved) return (ele1 + ele2); else return (gsf1 + gsf2);")
                .Define("H",            "if (isResolved) return (ele1 + ele2 + pho); else return (ele1 + pho);")
                .Define("Hgsf",         "if (isResolved) return H; else return (gsf1 + gsf2 + pho);")

                //! FIXEDME:
                //*     1) 3-body invariant mass cut -> 105 is ued for now
                //*     2) Ratio cut also need to modify once mass cut is changed
                .Define("GsfPtCut",     "if (isM2) return (gsf1.Pt() + gsf2.Pt()) > 44.; else return true;")
                .Define("RatioCut",     "if (isResolved) return (diEle.Pt()/H.M() > 25./105.); else return (ele1.Pt()/H.M() > 25./105.);")
                .Define("diMassCut",    "if (isM1) return true; else return diEle.M() < 50.")

                .Filter("GsfPtCut", "GsfPtCut")
                .Filter("diMassCut", "diMassCut")
                .Filter("RatioCut && pho.Pt()/H.M() > 25./105.", "RatioCut")
                .Filter("ele1.DeltaR(pho) > 1. && ele2.DeltaR(pho) > 1.", "dRCut")
                .Filter("H.M() > 105. && H.M() < 170.", "mHCut");
    return nf;
}


template <typename T>
auto CatGoodEvents(T &df, int year, bool isMC, string sysOption){

    // Conditionally applying a filter and a definition are simple with ROOT::RDF::RNode
    auto ff = ROOT::RDF::RNode(df);
    if (!isMC){
        ff = ff.Define("jetSmearPt",        "jetPt")
               .Define("jetSmearEn",        "jetEn")
               .Define("phoCorrR95x5",      "phoR9Full5x5[phoIdx1]");
    }
    else{
        // systematics uncertainties
        //*     1) UnPhoR9: estimate the effect on the signal yields w/o appling SS correction on R9
        //*     2) UnJECXX/UnJERXX: estimate the effect on the signal yields of jet smearing and correction
        if (sysOption == "UnPhoR9")
            ff = ff.Define("phoCorrR95x5",  "phoR9Full5x5[phoIdx1]");
        else
            ff = ff.Define("phoCorrR95x5",  Form("doR9SSCorrections(phoCalibEt[phoIdx1], phoPhi[phoIdx1], phoSCEta[phoIdx1], phoSCEtaWidth[phoIdx1], phoSCPhiWidth[phoIdx1], phoSigmaIEtaIEtaFull5x5[phoIdx1], phoSigmaIEtaIPhiFull5x5[phoIdx1], phoR9Full5x5[phoIdx1], phoE2x2Full5x5[phoIdx1], phoE5x5Full5x5[phoIdx1], rho, %d)", year));

        if (sysOption == "Nominal" || sysOption == "UnPhoR9")
            ff = ff.Define("jetSmearPt",    "jetP4Smear * jetPt")
                   .Define("jetSmearEn",    "jetP4Smear * jetEn");
        else if (sysOption == "UnJERUp")
            ff = ff.Define("jetSmearPt",    "jetP4SmearUp * jetPt")
                   .Define("jetSmearEn",    "jetP4SmearUp * jetEn");
        else if (sysOption == "UnJERDo")
            ff = ff.Define("jetSmearPt",    "jetP4SmearDo * jetPt")
                   .Define("jetSmearEn",    "jetP4SmearDo * jetEn");
        else if (sysOption == "UnJECUp")
            ff = ff.Define("jetSmearPt",    "jetP4Smear * jetPt + jetJECUnc")
                   .Define("jetSmearEn",    "jetP4Smear * jetEn + jetJECUnc");
        else if (sysOption == "UnJECDo")
            ff = ff.Define("jetSmearPt",    "jetP4Smear * jetPt - jetJECUnc")
                   .Define("jetSmearEn",    "jetP4Smear * jetEn - jetJECUnc");
        else{
            cout << "[ERROR] This systematics options is not allowed! -> " << sysOption << endl;
            exit(-1);
        }
    }

    auto nf = ff.Define("dR_ele1_jet",  "dRVector(ele1, jetSmearPt, jetEta, jetPhi, jetSmearEn)")
                .Define("dR_pho_jet",   "dRVector(pho, jetSmearPt, jetEta, jetPhi, jetSmearEn)")
                .Define("JetID",        Form("CutBasedJet(nJet, jetEta, jetNHF, jetNEF, jetNNP, jetNCH, jetCHF, jetCEF, %d)", year))
                .Define("isGoodJet",    "jetSmearPt > 30. && abs(jetEta) < 4.7 && JetID && dR_ele1_jet > 0.4 && dR_pho_jet > 0.4")
                .Define("vbftag",       "VBFtag(isGoodJet, jetSmearPt, jetEta, jetPhi, jetSmearEn, H)")

                .Define("isHVbf",       "vbftag[0] == 1")
                .Define("isLVbf",       "vbftag[0] == 2")
                .Define("jetIdx1",      "vbftag[1]")
                .Define("jetIdx2",      "vbftag[2]")
                .Define("isBst",        "H.Pt() > 60.")
                .Define("isEBHR9",      "abs(phoSCEta[phoIdx1]) < 1.4442 && phoCorrR95x5 > 0.96")
                .Define("isEBLR9",      "abs(phoSCEta[phoIdx1]) < 1.4442 && phoCorrR95x5 <= 0.96")
                .Define("isEE",         "abs(phoSCEta[phoIdx1]) > 1.566 && abs(phoSCEta[phoIdx1]) < 2.5")

                .Define("category",     "make_cat(isM2, isM1, isResolved, isHVbf, isLVbf, isBst, isEE, isEBHR9, isEBLR9)");
    return nf;
}


template <typename T>
auto CheckGenerator(T &df, bool isMC){
    if (!isMC) return df;

    auto nf = df.Define("GenIdx_reco1",     "GenMatch(ele1, nMC, mcPt, mcEta, mcPhi, mcMass, mcPID, mcMomPID, mcGMomPID, mcStatusFlag)")
                .Define("GenIdx_reco2",     "GenMatch(ele2, nMC, mcPt, mcEta, mcPhi, mcMass, mcPID, mcMomPID, mcGMomPID, mcStatusFlag)")
                .Define("nReco1MatchedGen", "GenIdx_reco1.size()")
                .Define("nReco2MatchedGen", "GenIdx_reco2.size()")
                .Define("isTrueM2",         "isM2 && nReco1MatchedGen > 1") // Reco1Gen2, nGsfMatchToReco >= 2
                .Define("isTrueM1",         "isM1 && nReco1MatchedGen > 1") // Reco1Gen2, nGsfMatchToReco == 1
                .Define("isTrueResolved",   "isResolved && nReco1MatchedGen == 1 && nReco2MatchedGen == 1 && GenIdx_reco1[0] != GenIdx_reco2[0]");
    return nf;
}


template <typename T>
auto DefineFinalVars(T &df){
    auto nf = df.Define("CMS_higgs_mass",           "H.M()")

                .Define("eleCalibPt_lep1",          "eleCalibPt[eleIdx1]")
                .Define("eleEta_lep1",              "eleEta[eleIdx1]")
                .Define("elePhi_lep1",              "elePhi[eleIdx1]")
                .Define("eleSCEn_lep1",             "eleSCEn[eleIdx1]")
                .Define("eleSCEta_lep1",            "eleSCEta[eleIdx1]")
                .Define("eleSCPhi_lep1",            "eleSCPhi[eleIdx1]")
                .Define("eleEn_lep1",               "eleEn[eleIdx1]")
                .Define("eleCalibEn_lep1",          "eleCalibEn[eleIdx1]")
                .Define("eleEcalEn_lep1",           "eleEcalEn[eleIdx1]")
                .Define("eleGsfDeltaR_lep1",        "gsfDeltaR[eleIdx1]")
                .Define("eleGsfPtRatio_lep1",       "gsfPtRatio[eleIdx1]")
                .Define("eleGsfRelPtRatio_lep1",    "gsfRelPtRatio[eleIdx1]")
                .Define("eleScale_up_lep1",         "hp::Mod(eleScale_stat_up[eleIdx1], eleScale_syst_up[eleIdx1], eleScale_gain_up[eleIdx1])")
                .Define("eleScale_dn_lep1",         "hp::Mod(eleScale_stat_dn[eleIdx1], eleScale_syst_dn[eleIdx1], eleScale_gain_dn[eleIdx1])")
                .Define("eleResol_up_lep1",         "hp::Mod(eleResol_rho_up[eleIdx1], eleResol_phi_up[eleIdx1])")
                .Define("eleResol_dn_lep1",         "hp::Mod(eleResol_rho_dn[eleIdx1], eleResol_phi_dn[eleIdx1])")

                .Define("eleCalibPt_lep2",          "eleCalibPt[eleIdx2]")
                .Define("eleEta_lep2",              "eleEta[eleIdx2]")
                .Define("elePhi_lep2",              "elePhi[eleIdx2]")
                .Define("eleSCEn_lep2",             "eleSCEn[eleIdx2]")
                .Define("eleSCEta_lep2",            "eleSCEta[eleIdx2]")
                .Define("eleSCPhi_lep2",            "eleSCPhi[eleIdx2]")
                .Define("eleEn_lep2",               "eleEn[eleIdx2]")
                .Define("eleCalibEn_lep2",          "eleCalibEn[eleIdx2]")
                .Define("eleEcalEn_lep2",           "eleEcalEn[eleIdx2]")
                .Define("eleGsfDeltaR_lep2",        "gsfDeltaR[eleIdx2]")
                .Define("eleGsfPtRatio_lep2",       "gsfPtRatio[eleIdx2]")
                .Define("eleGsfRelPtRatio_lep2",    "gsfRelPtRatio[eleIdx2]")
                .Define("eleScale_up_lep2",         "hp::Mod(eleScale_stat_up[eleIdx2], eleScale_syst_up[eleIdx2], eleScale_gain_up[eleIdx2])")
                .Define("eleScale_dn_lep2",         "hp::Mod(eleScale_stat_dn[eleIdx2], eleScale_syst_dn[eleIdx2], eleScale_gain_dn[eleIdx2])")
                .Define("eleResol_up_lep2",         "hp::Mod(eleResol_rho_up[eleIdx2], eleResol_phi_up[eleIdx2])")
                .Define("eleResol_dn_lep2",         "hp::Mod(eleResol_rho_dn[eleIdx2], eleResol_phi_dn[eleIdx2])")

                .Define("phoCalibEt_lep1",          "phoCalibEt[phoIdx1]")
                .Define("phoR9_lep1",               "phoR9[phoIdx1]")
                .Define("phoCorrR95x5_lep1",        "phoCorrR95x5")
                .Define("phoEta_lep1",              "phoEta[phoIdx1]")
                .Define("phoPhi_lep1",              "phoPhi[phoIdx1]")
                .Define("phoSCEta_lep1",            "phoSCEta[phoIdx1]")
                .Define("phoSCPhi_lep1",            "phoSCPhi[phoIdx1]")
                .Define("phoScale_up_lep1",         "hp::Mod(phoScale_stat_up[phoIdx1], phoScale_syst_up[phoIdx1], phoScale_gain_up[phoIdx1])")
                .Define("phoScale_dn_lep1",         "hp::Mod(phoScale_stat_dn[phoIdx1], phoScale_syst_dn[phoIdx1], phoScale_gain_dn[phoIdx1])")
                .Define("phoResol_up_lep1",         "hp::Mod(phoResol_rho_up[phoIdx1], phoResol_phi_up[phoIdx1])")
                .Define("phoResol_dn_lep1",         "hp::Mod(phoResol_rho_dn[phoIdx1], phoResol_phi_dn[phoIdx1])")

                .Define("dR_ele1g",                 "ele1.DeltaR(pho)")
                .Define("dR_ele2g",                 "ele2.DeltaR(pho)")
                .Define("jet1",                     "hp::P4En(jetSmearPt[jetIdx1], jetEta[jetIdx1], jetPhi[jetIdx1], jetSmearEn[jetIdx1])")
                .Define("jet2",                     "hp::P4En(jetSmearPt[jetIdx2], jetEta[jetIdx2], jetPhi[jetIdx2], jetSmearEn[jetIdx2])")
                .Define("dijet",                    "jet1 + jet2")
                .Define("zepen",                    "H.Eta() - ((jet2.Eta() + jet1.Eta()) * 0.5)");
    return nf;
}


template <typename T>
auto AddWeights(T &df, string era, int year, bool isMC){
    if (!isMC) return df;

    // set up the puwei calculator
    PUWeightCalculator* puCalc[3] = {NULL, NULL, NULL};
    for (int i = 0; i < 3; i++){
        puCalc[i] = new PUWeightCalculator();
    }
    puCalc[0]->Init(PUfile(year, "nominal").c_str());
    puCalc[1]->Init(PUfile(year, "up").c_str());
    puCalc[2]->Init(PUfile(year, "down").c_str());

    auto get_pu = [&, puCalc](int run, ROOT::RVec<float>& puTrue){
        float puwei = puCalc[0]->GetWeight(run, puTrue[1]);
        return puwei;
    };
    auto get_pu_up = [&, puCalc](int run, ROOT::RVec<float>& puTrue){
        float puwei_up = puCalc[1]->GetWeight(run, puTrue[1]);
        return puwei_up;
    };
    auto get_pu_do = [&, puCalc](int run, ROOT::RVec<float>& puTrue){
        float puwei_down = puCalc[2]->GetWeight(run, puTrue[1]);
        return puwei_down;
    };

    // load the histogram of Fall17 V2 electron MVA ID SFs
    map<string, string> fFallMap;
    fFallMap["2016_preVFP"] = "./external/SFfiles/Fall17EleID/egammaEffi.txt_Ele_wp90iso_preVFP_EGM2D.root";
    fFallMap["2016_postVFP"] = "./external/SFfiles/Fall17EleID/egammaEffi.txt_Ele_wp90iso_postVFP_EGM2D.root";
    fFallMap["2017"] = "./external/SFfiles/Fall17EleID/egammaEffi.txt_EGM2D_MVA90iso_UL17.root";
    fFallMap["2018"] = "./external/SFfiles/Fall17EleID/egammaEffi.txt_Ele_wp90iso_EGM2D.root";
    TFile* fFall = new TFile(fFallMap[era.c_str()].c_str(), "READ");
    TH2F* hFall = fFall->Get<TH2F>("EGamma_SF2D");

    auto get_FallSFs = [&, hFall](float eleSCEta, float eleCalibPt){

        float SF = hFall->GetBinContent(hFall->FindBin(eleSCEta, eleCalibPt));
        return SF;
    };
    auto get_FallSFsErr = [&, hFall](float eleSCEta, float eleCalibPt){
        float SF_err = hFall->GetBinError(hFall->FindBin(eleSCEta, eleCalibPt));
        return SF_err;
    };

    // load the histogram of Merged electron MVA ID SFs
    TFile* fMerg = new TFile("./external/SFfiles/MergedID/MergedIDSFs_combined.root", "READ");
    TH2D* hMerg = fMerg->Get<TH2D>("hSF_2D");

    auto get_MergSFs = [&, hMerg](float eleSCEta, float eleCalibPt){
        int binx = hMerg->GetXaxis()->FindBin(abs(eleSCEta));
        int biny = hMerg->GetYaxis()->FindBin(eleCalibPt);
        int Nbiny = hMerg->GetYaxis()->GetNbins();
        if (biny > Nbiny)
            biny = Nbiny; // if pt > the last pt bin of hSF, use the sf of the last pt bin

        float SF = hMerg->GetBinContent(binx, biny);
        return SF;
    };
    auto get_MergSFsErr = [&, hMerg](float eleSCEta, float eleCalibPt){
        int binx = hMerg->GetXaxis()->FindBin(abs(eleSCEta));
        int biny = hMerg->GetYaxis()->FindBin(eleCalibPt);
        int Nbiny = hMerg->GetYaxis()->GetNbins();
        if (biny > Nbiny)
            biny = Nbiny; // if pt > the last pt bin of hSF, use the sf of the last pt bin

        float SF_err = hMerg->GetBinError(hMerg->FindBin(binx, biny));
        return SF_err;
    };

    // load the histogram of conversion safe electron veto SFs
    map<string, string> fVetoMap;
    fVetoMap["2016_preVFP"] = "./external/SFfiles/eleVeto/CSEV_SummaryPlot_UL16_preVFP.root";
    fVetoMap["2016_postVFP"] = "./external/SFfiles/eleVeto/CSEV_SummaryPlot_UL16_postVFP.root";
    fVetoMap["2017"] = "./external/SFfiles/eleVeto/CSEV_SummaryPlot_UL17.root";
    fVetoMap["2018"] = "./external/SFfiles/eleVeto/CSEV_SummaryPlot_UL18.root";
    TFile* fVeto = new TFile(fVetoMap[era.c_str()].c_str(), "READ");
    TH1F* hVeto = fVeto->Get<TH1F>("MVAID/SF_CSEV_MVAID");

    auto get_VetoSFs = [&, hVeto](float phoSCEta, float phoR9){
        float SF = 1.;
        if (fabs(phoSCEta) < 1.4446 && phoR9 > 0.96)
            SF = hVeto->GetBinContent(2);
        else if (fabs(phoSCEta) < 1.4446 && phoR9 < 0.96)
            SF = hVeto->GetBinContent(3);
        else if (fabs(phoSCEta) > 1.566 && phoR9 > 0.96)
            SF = hVeto->GetBinContent(5);
        else
            SF = hVeto->GetBinContent(6); // fabs(phoSCEta) > 1.566 && phoR9 < 0.96

        return SF;
    };
    auto get_VetoSFsErr = [&, hVeto](float phoSCEta, float phoR9){
        float SF_err = 1.;
        if (fabs(phoSCEta) < 1.4446 && phoR9 > 0.96)
            SF_err = hVeto->GetBinError(2);
        else if (fabs(phoSCEta) < 1.4446 && phoR9 < 0.96)
            SF_err = hVeto->GetBinError(3);
        else if (fabs(phoSCEta) > 1.566 && phoR9 > 0.96)
            SF_err = hVeto->GetBinError(5);
        else
            SF_err = hVeto->GetBinError(6); // fabs(phoSCEta) > 1.566 && phoR9 < 0.96

        return SF_err;
    };

    // load the histogram of photon ID SFs
    map<string, string> fPhoMap;
    fPhoMap["2016_preVFP"] = "./external/SFfiles/phoID/egammaEffi.txt_EGM2D_Pho_wp90_UL16.root";
    fPhoMap["2016_postVFP"] = "./external/SFfiles/phoID/egammaEffi.txt_EGM2D_Pho_MVA90_UL16_postVFP.root";
    fPhoMap["2017"] = "./external/SFfiles/phoID/egammaEffi.txt_EGM2D_PHO_MVA90_UL17.root";
    fPhoMap["2018"] = "./external/SFfiles/phoID/egammaEffi.txt_EGM2D_Pho_wp90.root_UL18.root.root";
    TFile* fPho = new TFile(fPhoMap[era.c_str()].c_str(), "READ");
    TH2F* hPho = fPho->Get<TH2F>("EGamma_SF2D");

    auto get_phoSFs = [&, hPho](float phoSCEta, float phoCalibEt){
        float SF = hPho->GetBinContent(hPho->FindBin(phoSCEta, phoCalibEt));
        return SF;
    };
    auto get_phoSFsErr = [&, hPho](float phoSCEta, float phoCalibEt){
        float SF_err = hPho->GetBinError(hPho->FindBin(phoSCEta, phoCalibEt));
        return SF_err;
    };

    auto nf = df.Define("puwei",            get_pu,             {"run", "puTrue"})
                .Define("puwei_up",         get_pu_up,          {"run", "puTrue"})
                .Define("puwei_down",       get_pu_do,          {"run", "puTrue"})
                .Define("Fall17IDSF",       get_FallSFs,        {"eleSCEta_lep1", "eleCalibPt_lep1"})
                .Define("Fall17IDSFErr",    get_FallSFsErr,     {"eleSCEta_lep1", "eleCalibPt_lep1"})
                .Define("Fall17IDSFUp",     "Fall17IDSF + Fall17IDSFErr")
                .Define("Fall17IDSFDo",     "Fall17IDSF - Fall17IDSFErr")
                .Define("MergedIDSF",       get_MergSFs,        {"eleSCEta_lep1", "eleCalibPt_lep1"})
                .Define("MergedIDSFErr",    get_MergSFsErr,     {"eleSCEta_lep1", "eleCalibPt_lep1"})
                .Define("MergedIDSFUp",     "MergedIDSF + MergedIDSFErr")
                .Define("MergedIDSFDo",     "MergedIDSF - MergedIDSFErr")
                .Define("eleVetoSF",        get_VetoSFs,        {"phoSCEta_lep1", "phoR9_lep1"})
                .Define("eleVetoSFErr",     get_VetoSFsErr,     {"phoSCEta_lep1", "phoR9_lep1"})
                .Define("eleVetoSFUp",      "eleVetoSF + eleVetoSFErr")
                .Define("eleVetoSFDo",      "eleVetoSF - eleVetoSFErr")
                .Define("phoIDSF",          get_phoSFs,         {"phoSCEta_lep1", "phoCalibEt_lep1"})
                .Define("phoIDSFErr",       get_phoSFsErr,      {"phoSCEta_lep1", "phoCalibEt_lep1"})
                .Define("phoIDSFUp",        "phoIDSF + phoIDSFErr")
                .Define("phoIDSFDo",        "phoIDSF - phoIDSFErr")
                .Define("genwei",           "if (genWeight > 0) return 1.; else return -1.;")
                .Define("L1PF",             "L1ECALPrefire")
                .Define("L1PFUp",           "L1ECALPrefireUp")
                .Define("L1PFDo",           "L1ECALPrefireDown")
                .Define("elewei",           "if (isResolved) return Fall17IDSF; else return MergedIDSF;")
                .Define("weight",           "mcwei * puwei * genwei * elewei * eleVetoSF * phoIDSF * phoConvWei * L1PF");
    return nf;
}


//====================================//
// Main function to do the analysis!! //
//====================================//
void rdfxAna(string infile, string outfile, int year, string era, int thread, bool isMC, string sysOption){
    TStopwatch time;
    time.Start();

    cout << "[INFO] Read_File(): " << infile.c_str() << endl;
    cout << "[INFO] Save_File(): " << outfile.c_str() << endl;

    // TMVA::Reader is not thread safety, so multi-thread can only be used for data :(
    // TMVA::Experimental::RReader allows multi-threading, but the performance is even worse than single thread :(
    // Not sure why
    if (!isMC){
        if (thread == -1){
            ROOT::EnableImplicitMT();
        }
        else{
            ROOT::EnableImplicitMT(thread);
        }
    }

    //=======================================================//
    // Load TTree into RDataFrame and then do the selections //
    //=======================================================//
    ROOT::RDataFrame df("ggNtuplizer/EventTree", infile.c_str());

    cout << "[INFO] Total number of events: " << df.Count().GetValue() << endl;
    cout << "[INFO] Pool size: " << df.GetNSlots() << endl;

    auto df1 = df.Define("isMC", [&, isMC]{return isMC;});
    auto df2 = FindGoodPho(df1, year);
    auto df3 = FindGoodEle(df2);
    auto df4 = FilterGoodEvents(df3, year);
    auto df5 = CatGoodEvents(df4, year, isMC, sysOption.c_str());
    auto df6 = CheckGenerator(df5, isMC);
    auto df7 = DefineFinalVars(df6);
    auto df8 = AddWeights(df7, era.c_str(), year, isMC);
    auto dfFin = df8;

    //========================================================//
    // Visualize the selection results (cut flow, event count)//
    //========================================================//
    cout << "[INFO] Cut flow:" << std::endl;
    auto report = dfFin.Report();
    report->Print();

    // cout the usefull information
    int M2 = dfFin.Filter("isM2").Count().GetValue();
    int M1 = dfFin.Filter("isM1").Count().GetValue();
    int Re = dfFin.Filter("isResolved").Count().GetValue();
    cout << "[INFO] Number of events:" << std::endl;
    cout << fixed << showpoint;
    cout << ">>> Merged-2Gsf: " << M2 << "(" << setprecision(2) << M2*100./(M2+M1+Re) << "%)" << endl;
    cout << ">>> Merged-1Gsf: " << M1 << "(" << setprecision(2) << M1*100./(M2+M1+Re) << "%)" << endl;
    cout << ">>> Resolved:    " << Re << "(" << setprecision(2) << Re*100./(M2+M1+Re) << "%)" << endl;

    if (isMC){
        int gM2 = dfFin.Filter("isTrueM2").Count().GetValue();
        int gM1 = dfFin.Filter("isTrueM1").Count().GetValue();
        int gRe = dfFin.Filter("isTrueResolved").Count().GetValue();
        cout << "[INFO] Purity of the selections:" << std::endl;
        cout << ">>> Merged-2Gsf: " << setprecision(2) << gM2*100./M2 << "%" << endl;
        cout << ">>> Merged-1Gsf: " << setprecision(2) << gM1*100./M1 << "%" << endl;
        cout << ">>> Resolved:    " << setprecision(2) << gRe*100./Re << "%" << endl;
    }

    //====================================================//
    // Define the final variables to save to the miniTree //
    //====================================================//
    // final vars to save
    vector<string> defColNames = dfFin.GetDefinedColumnNames();
    vector<string> Vars = {
        "category", "isM2", "isM1", "isHVbf", "isLVbf", "isBst", "isEBHR9", "isEBLR9", "isEE", "isResolved",
        "run", "lumis", "event",
        "dR_ele1g", "dR_ele2g", "CMS_higgs_mass", "zepen"
    };
    for (int i = 0; i < defColNames.size(); i++){
        size_t foundLep1 = defColNames[i].find("_lep1");
        size_t foundLep2 = defColNames[i].find("_lep2");

        if (foundLep1 != string::npos || foundLep2 != string::npos)
            Vars.push_back(defColNames[i]);
        if (dfFin.GetColumnType(defColNames[i]) == "TLorentzVector")
            Vars.push_back(defColNames[i]);
    }

    if (isMC){
        // concatenate Vars and weiVars
        vector<string> weiVars = {
            "puwei", "puwei_up", "puwei_down",
            "Fall17IDSF", "Fall17IDSFUp", "Fall17IDSFDo",
            "MergedIDSF", "MergedIDSFUp", "MergedIDSFDo",
            "eleVetoSF", "eleVetoSFUp", "eleVetoSFDo",
            "phoIDSF", "phoIDSFUp", "phoIDSFDo",
            "genwei", "mcwei", "L1PF", "L1PFUp", "L1PFDo", "weight",
            "nReco1MatchedGen", "nReco2MatchedGen", "isTrueM2", "isTrueM1", "isTrueResolved"
        };
        Vars.insert(Vars.end(), weiVars.begin(), weiVars.end());
    }
    vector<string> sysVars = {"category", "CMS_higgs_mass", "weight"};
    vector<string> finalVars = (sysOption == "Nominal") ? Vars : sysVars;

    // save the minitree
    dfFin.Snapshot("outTree", outfile.c_str(), finalVars);

    cout << "[INFO] Time taken: " << endl;
    time.Stop();
    time.Print();

    cout << endl;
}