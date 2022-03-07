#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "./pluginsV2/puweicalc.h"
#include "./pluginsV2/EleSelections.h"
#include "./pluginsV2/JetSelections.h"
#include "./pluginsV2/PhoR9SSCorrection.h"
#include "./pluginsV2/GenMatching.h"


// TODO: 
//! [1] add uncertaities 
//*     1) Hgg preselection
//! [2] some SFs need to be added
////    1) Fall17 V2 MVA ID
//*     2) Merged ID
//*     3) electron conversion veto and hits
//*     4) di-photn and dielectron trigger SFs 
//*     5) Hgg preselection


using namespace std;
using namespace ROOT::VecOps;

namespace Helper{
    template <typename T>
    bool HasMC(T &df){
        vector<string> colNames = df.GetColumnNames();
        if (count(colNames.begin(), colNames.end(), "nMC")) return true;
        else return false;
    }

    // Function to get the index vector sorted by pT
    // * Reference: https://root.cern/doc/master/vo006__IndexManipulation_8C.html
    ROOT::RVec<int> getIdx(ROOT::RVec<int>& isgood, ROOT::RVec<float>& pt){
        ROOT::RVec<int> idx_select = Nonzero(isgood);
        ROOT::RVec<int> idx_sort = Reverse(Argsort(pt));
        ROOT::RVec<int> idx = Intersect(idx_sort, idx_select);

        return idx;
    }

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
}


template <typename T>
auto FindGoodPho(T &df, int year) {
    const float mva_EB = (year == 2016) ? 0.2 : -0.02;
    const float mva_EE = (year == 2016) ? 0.2 : -0.26;

    auto nf = df.Define("isEBPho",      Form("abs(phoSCEta) < 1.4442 && phoIDMVA > %f", mva_EB))
                .Define("isEEPho",      Form("abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5 && phoIDMVA > %f", mva_EE))
                .Define("phoCuts",      "phoEleVeto != 0 && phoCalibEt > 15")
                .Define("isGoodPho",    "(isEBPho || isEEPho) && phoCuts")
                .Define("phoIdx1",      "if (Sum(isGoodPho) > 0) return Helper::getIdx(isGoodPho, phoCalibEt)[0]; else return -1;");
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

                .Define("m2Idx1",       "if (Sum(isM2GoodEle) > 0) return Helper::getIdx(isM2GoodEle, eleCalibPt)[0]; else return -1;")
                .Define("m1Idx1",       "if (Sum(isM1GoodEle) > 0) return Helper::getIdx(isM1GoodEle, eleCalibPt)[0]; else return -1;")
                .Define("reIdx1",       "if (Sum(isReGoodEle) > 1) return Helper::getIdx(isReGoodEle, eleCalibPt)[0]; else return -1;")
                .Define("reIdx2",       "if (Sum(isReGoodEle) > 1) return Helper::getIdx(isReGoodEle, eleCalibPt)[1]; else return -1;");
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

    const float pT1cut = (year == 2016) ? 33. : 33.;
    const float pT2cut = (year == 2016) ? 21. : 25.;

    auto nf = df.Filter("Sum(isM2GoodEle) > 0 || Sum(isM1GoodEle) ||  Sum(isReGoodEle) > 1", "good ele")
                .Filter("Sum(isGoodPho) > 0", "good pho")

                .Define("isDiPhoHLT",   "((HLTPho >> 14) & 1) == 1")
                .Define("isDiEleHLT",   isDiEleTrigger)
                .Define("isM2",         Form("Sum(isM2GoodEle) > 0 && isDiPhoHLT && TMath::Max(eleCalibPt[m2Idx1], phoCalibEt[phoIdx1]) > %f && TMath::Min(eleCalibPt[m2Idx1], phoCalibEt[phoIdx1]) > %f", pT1cut, pT2cut))
                .Define("isM1",         Form("isM2 != 1 && Sum(isM1GoodEle) > 0 && isDiPhoHLT && TMath::Max(eleCalibPt[m1Idx1], phoCalibEt[phoIdx1]) > %f && TMath::Min(eleCalibPt[m1Idx1], phoCalibEt[phoIdx1]) > %f", pT1cut, pT2cut))
                .Define("isResolved",   "isM2 != 1 && isM1 != 1 && Sum(isReGoodEle) > 1 && isDiEleHLT && TMath::Max(eleCalibPt[reIdx1], eleCalibPt[reIdx2]) > 25. && TMath::Min(eleCalibPt[reIdx1], eleCalibPt[reIdx2]) > 15.")
                
                .Filter("isM2 || isM1 || isResolved", "HLT/PtCut")

                // kinematic selections: construct physics object
                .Define("eleidx1",      "if (isM2) return m2Idx1; else if (isM1) return m1Idx1; else return reIdx1;")
                .Define("eleidx2",      "if (isM2) return m2Idx1; else if (isM1) return m1Idx1; else return reIdx2;")
                .Define("ele1",         "TLorentzVector v; v.SetPtEtaPhiE(eleCalibEn[eleidx1]/cosh(eleEta[eleidx1]), eleEta[eleidx1], elePhi[eleidx1], eleCalibEn[eleidx1]); return v;")
                .Define("ele2",         "TLorentzVector v; v.SetPtEtaPhiE(eleCalibEn[eleidx2]/cosh(eleEta[eleidx2]), eleEta[eleidx2], elePhi[eleidx2], eleCalibEn[eleidx2]); return v;")
                .Define("pho",          "TLorentzVector v; v.SetPtEtaPhiM(phoCalibEt[phoIdx1], phoEta[phoIdx1], phoPhi[phoIdx1], 0.); return v;")
                
                // only for merged electrons
                .Define("gsf1",         "TLorentzVector v; v.SetPtEtaPhiM(eleTrkPt[eleidx1], eleTrkEta[eleidx1], eleTrkPhi[eleidx1], 0.000511); return v;")
                .Define("gsf2",         "TLorentzVector v; v.SetPtEtaPhiM(eleSubTrkPt[eleidx1], eleSubTrkEta[eleidx1], eleSubTrkPhi[eleidx1], 0.000511); return v;")
                
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

                .Define("category",     "Helper::make_cat(isM2, isM1, isResolved, isHVbf, isLVbf, isBst, isEE, isEBHR9, isEBLR9)");
    return nf;
} 


template <typename T>
auto CheckGenerator(T &df, bool isMC){
    if (!isMC) return df;
    
    auto nf = df.Define("nReco1MatchedGen", "GenMatch(ele1, nMC, mcPt, mcEta, mcPhi, mcMass, mcPID, mcStatusFlag)")
                .Define("nReco2MatchedGen", "GenMatch(ele2, nMC, mcPt, mcEta, mcPhi, mcMass, mcPID, mcStatusFlag)")
                .Define("isTrueM2",         "isM2 && nReco1MatchedGen > 1")
                .Define("isTrueM1",         "isM1 && nReco1MatchedGen > 1")
                .Define("isTrueResolved",   "isResolved && nReco1MatchedGen == 1 && nReco2MatchedGen == 1");
    return nf;
}


template <typename T>
auto DefineFinalVars(T &df){
    auto nf = df.Define("CMS_higgs_mass",       "H.M()")

                .Define("eleCalibPt1",          "eleCalibPt[eleidx1]")
                .Define("eleEta1",              "eleEta[eleidx1]")
                .Define("elePhi1",              "elePhi[eleidx1]")
                .Define("eleSCEn1",             "eleSCEn[eleidx1]")
                .Define("eleSCEta1",            "eleSCEta[eleidx1]")
                .Define("eleSCPhi1",            "eleSCPhi[eleidx1]")
                .Define("eleEn1",               "eleEn[eleidx1]")
                .Define("eleCalibEn1",          "eleCalibEn[eleidx1]")
                .Define("eleEcalEn1",           "eleEcalEn[eleidx1]")
                .Define("eleGsfDeltaR1",        "gsfDeltaR[eleidx1]")
                .Define("eleGsfPtRatio1",       "gsfPtRatio[eleidx1]")
                .Define("eleGsfRelPtRatio1",    "gsfRelPtRatio[eleidx1]")

                .Define("eleCalibPt2",          "eleCalibPt[eleidx2]")
                .Define("eleEta2",              "eleEta[eleidx2]")
                .Define("elePhi2",              "elePhi[eleidx2]")
                .Define("eleSCEn2",             "eleSCEn[eleidx2]")
                .Define("eleSCEta2",            "eleSCEta[eleidx2]")
                .Define("eleSCPhi2",            "eleSCPhi[eleidx2]")
                .Define("eleEn2",               "eleEn[eleidx2]")
                .Define("eleCalibEn2",          "eleCalibEn[eleidx2]")
                .Define("eleEcalEn2",           "eleEcalEn[eleidx2]")
                .Define("eleGsfDeltaR2",        "gsfDeltaR[eleidx2]")
                .Define("eleGsfPtRatio2",       "gsfPtRatio[eleidx2]")
                .Define("eleGsfRelPtRatio2",    "gsfRelPtRatio[eleidx2]")

                .Define("phoCalibEt1",          "phoCalibEt[phoIdx1]")
                .Define("phoR91",               "phoR9[phoIdx1]")
                .Define("phoCorrR95x51",        "phoCorrR95x5")
                .Define("phoEta1",              "phoEta[phoIdx1]")
                .Define("phoPhi1",              "phoPhi[phoIdx1]")
                .Define("phoSCEta1",            "phoSCEta[phoIdx1]")
                .Define("phoSCPhi1",            "phoSCPhi[phoIdx1]")

                .Define("dR_ele1g",             "ele1.DeltaR(pho)")
                .Define("dR_ele2g",             "ele2.DeltaR(pho)")

                .Define("jet1",                 "TLorentzVector v; v.SetPtEtaPhiE(jetSmearPt[jetIdx1], jetEta[jetIdx1], jetPhi[jetIdx1], jetSmearEn[jetIdx1]); return v;")
                .Define("jet2",                 "TLorentzVector v; v.SetPtEtaPhiE(jetSmearPt[jetIdx2], jetEta[jetIdx2], jetPhi[jetIdx2], jetSmearEn[jetIdx2]); return v;")
                .Define("dijet",                "jet1 + jet2")
                .Define("zepen",                "H.Eta() - ((jet2.Eta() + jet1.Eta()) * 0.5)")

                .Define("eleScale_stat_up1",    "eleScale_stat_up[eleidx1]")
                .Define("eleScale_stat_dn1",    "eleScale_stat_dn[eleidx1]")
                .Define("eleScale_syst_up1",    "eleScale_syst_up[eleidx1]")
                .Define("eleScale_syst_dn1",    "eleScale_syst_dn[eleidx1]")
                .Define("eleScale_gain_up1",    "eleScale_gain_up[eleidx1]")
                .Define("eleScale_gain_dn1",    "eleScale_gain_dn[eleidx1]")
                .Define("eleResol_rho_up1",     "eleResol_rho_up[eleidx1]")
                .Define("eleResol_rho_dn1",     "eleResol_rho_dn[eleidx1]")
                .Define("eleResol_phi_up1",     "eleResol_phi_up[eleidx1]")
                .Define("eleResol_phi_dn1",     "eleResol_phi_dn[eleidx1]")

                .Define("eleScale_stat_up2",    "eleScale_stat_up[eleidx2]")
                .Define("eleScale_stat_dn2",    "eleScale_stat_dn[eleidx2]")
                .Define("eleScale_syst_up2",    "eleScale_syst_up[eleidx2]")
                .Define("eleScale_syst_dn2",    "eleScale_syst_dn[eleidx2]")
                .Define("eleScale_gain_up2",    "eleScale_gain_up[eleidx2]")
                .Define("eleScale_gain_dn2",    "eleScale_gain_dn[eleidx2]")
                .Define("eleResol_rho_up2",     "eleResol_rho_up[eleidx2]")
                .Define("eleResol_rho_dn2",     "eleResol_rho_dn[eleidx2]")
                .Define("eleResol_phi_up2",     "eleResol_phi_up[eleidx2]")
                .Define("eleResol_phi_dn2",     "eleResol_phi_dn[eleidx2]")

                .Define("phoScale_stat_up1",    "phoScale_stat_up[phoIdx1]")
                .Define("phoScale_stat_dn1",    "phoScale_stat_dn[phoIdx1]")
                .Define("phoScale_syst_up1",    "phoScale_syst_up[phoIdx1]")
                .Define("phoScale_syst_dn1",    "phoScale_syst_dn[phoIdx1]")
                .Define("phoScale_gain_up1",    "phoScale_gain_up[phoIdx1]")
                .Define("phoScale_gain_dn1",    "phoScale_gain_dn[phoIdx1]")
                .Define("phoResol_rho_up1",     "phoResol_rho_up[phoIdx1]")
                .Define("phoResol_rho_dn1",     "phoResol_rho_dn[phoIdx1]")
                .Define("phoResol_phi_up1",     "phoResol_phi_up[phoIdx1]")
                .Define("phoResol_phi_dn1",     "phoResol_phi_dn[phoIdx1]");

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
                .Define("Fall17IDSF",       get_FallSFs,        {"eleSCEta1", "eleCalibPt1"})
                .Define("Fall17IDSFErr",    get_FallSFsErr,     {"eleSCEta1", "eleCalibPt1"})
                .Define("Fall17IDSFUp",     "Fall17IDSF + Fall17IDSFErr")
                .Define("Fall17IDSFDo",     "Fall17IDSF - Fall17IDSFErr")

                .Define("eleVetoSF",        get_VetoSFs,        {"phoSCEta1", "phoR91"})
                .Define("eleVetoSFErr",     get_VetoSFsErr,     {"phoSCEta1", "phoR91"})
                .Define("eleVetoSFUp",      "eleVetoSF + eleVetoSFErr")
                .Define("eleVetoSFDo",      "eleVetoSF - eleVetoSFErr")

                .Define("phoIDSF",          get_phoSFs,         {"phoSCEta1", "phoCalibEt1"})
                .Define("phoIDSFErr",       get_phoSFsErr,      {"phoSCEta1", "phoCalibEt1"})
                .Define("phoIDSFUp",        "phoIDSF + phoIDSFErr")
                .Define("phoIDSFDo",        "phoIDSF - phoIDSFErr")

                .Define("genwei",           "if (genWeight > 0) return 1.; else return -1.;")
                .Define("L1PF",             "L1ECALPrefire")
                .Define("L1PFUp",           "L1ECALPrefireUp")
                .Define("L1PFDo",           "L1ECALPrefireDown")
                .Define("weight",           "mcwei * puwei * genwei * Fall17IDSF * eleVetoSF * phoIDSF * phoConvWei * L1PF");
    
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
    vector<string> Vars = {
        "category", "isM2", "isM1", "isHVbf", "isLVbf", "isBst", "isEBHR9", "isEBLR9", "isEE", "isResolved",
        "run", "lumis", "event", 
        "eleCalibPt1", "eleEta1", "elePhi1", "eleSCEn1", "eleSCEta1", "eleSCPhi1", "eleEn1", "eleCalibEn1", "eleEcalEn1", "eleGsfDeltaR1", "eleGsfPtRatio1", "eleGsfRelPtRatio1", 
        "eleCalibPt2", "eleEta2", "elePhi2", "eleSCEn2", "eleSCEta2", "eleSCPhi2", "eleEn2", "eleCalibEn2", "eleEcalEn2", "eleGsfDeltaR2", "eleGsfPtRatio2", "eleGsfRelPtRatio2", 
        "phoCalibEt1", "phoR91", "phoCorrR95x51", "phoEta1", "phoPhi1", "phoSCEta1", "phoSCPhi1",
        "dR_ele1g", "dR_ele2g",
        "CMS_higgs_mass",
        "eleScale_stat_up1", "eleScale_stat_dn1", "eleScale_syst_up1", "eleScale_syst_dn1", "eleScale_gain_up1", "eleScale_gain_dn1", "eleResol_rho_up1", "eleResol_rho_dn2", "eleResol_phi_up1", "eleResol_phi_dn1",
        "eleScale_stat_up2", "eleScale_stat_dn2", "eleScale_syst_up2", "eleScale_syst_dn2", "eleScale_gain_up2", "eleScale_gain_dn2", "eleResol_rho_up2", "eleResol_rho_dn2", "eleResol_phi_up2", "eleResol_phi_dn2",
        "phoScale_stat_up1", "phoScale_stat_dn1", "phoScale_syst_up1", "phoScale_syst_dn1", "phoScale_gain_up1", "phoScale_gain_dn1", "phoResol_rho_up1", "phoResol_rho_dn1", "phoResol_phi_up1", "phoResol_phi_dn1",
        "zepen",
        "ele1", "ele2", "gsf1", "gsf2", "diEle", "pho", "H", "Hgsf", "jet1", "jet2", "dijet"
    };
    if (isMC){
        // concatenate Vars and weiVars
        vector<string> weiVars = {
            "puwei", "puwei_up", "puwei_down", 
            "Fall17IDSF", "Fall17IDSFUp", "Fall17IDSFDo", 
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