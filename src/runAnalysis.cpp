#include <iostream>
#include <filesystem> // c++17
#include <algorithm>
#include <vector>
#include <map>

#include "yaml-cpp/yaml.h"
#include "boost/program_options.hpp"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TMath.h"

#include "WeightHandler.h"
#include "EnCalibrater.h"
#include "Utilities.h"
#include "PhoR9SSCorrection.h"
#include "PhotonSel.h"
#include "GsfTracks.h"
#include "ElectronSel.h"
#include "MergedIDPred.h"
#include "Categorizer.h"
#include "Generator.h"

#define CGREEN "\033[0;32m"
#define CEND "\033[0m"
#define cprintf(X) printf("%s%s%s", CGREEN, X, CEND)

namespace fs = std::filesystem;
namespace po = boost::program_options;
template <typename T>
using Vec = const ROOT::RVec<T>&;


ROOT::RDF::RNode FindGoodPho(ROOT::RDF::RNode df, bool isMC, std::string vary) {

    // systematic variation of photon energy
    auto ff = df;
    if (!isMC)                            ff = ff.Define("phoCalibEt_NEW", "phoCalibEt"); // data
    if (isMC && vary == "Nominal")        ff = ff.Define("phoCalibEt_NEW", "phoCalibEt"); // nominal mc
    if (isMC && vary == "PhoScaleStatUp") ff = ff.Define("phoCalibEt_NEW", "phoScale_stat_up/cosh(phoEta)");
    if (isMC && vary == "PhoScaleSystUp") ff = ff.Define("phoCalibEt_NEW", "phoScale_syst_up/cosh(phoEta)");
    if (isMC && vary == "PhoScaleGainUp") ff = ff.Define("phoCalibEt_NEW", "phoScale_gain_up/cosh(phoEta)");
    if (isMC && vary == "PhoScaleStatDo") ff = ff.Define("phoCalibEt_NEW", "phoScale_stat_dn/cosh(phoEta)");
    if (isMC && vary == "PhoScaleSystDo") ff = ff.Define("phoCalibEt_NEW", "phoScale_syst_dn/cosh(phoEta)");
    if (isMC && vary == "PhoScaleGainDo") ff = ff.Define("phoCalibEt_NEW", "phoScale_gain_dn/cosh(phoEta)");
    if (isMC && vary == "PhoSigmaPhiUp")  ff = ff.Define("phoCalibEt_NEW", "phoResol_phi_up/cosh(phoEta)");
    if (isMC && vary == "PhoSigmaRhoUp")  ff = ff.Define("phoCalibEt_NEW", "phoResol_rho_up/cosh(phoEta)");
    if (isMC && vary == "PhoSigmaRhoDo")  ff = ff.Define("phoCalibEt_NEW", "phoResol_rho_dn/cosh(phoEta)");

    // Fall17MVA ID 90% WP
    auto nf = ff.Define("isEBPho",              "abs(phoSCEta) < 1.4442 && phoIDMVA > -0.02")
                .Define("isEEPho",              "abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5 && phoIDMVA > -0.26")
                .Define("phoCuts",              "phoEleVeto != 0 && phoCalibEt_NEW > 15")
                .Define("isHggPho",             phoSel::HggPresel,          {"nPho", "rhoAll", "phoSCEta", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"})

                // good mpho for merged category (an additional Hgg preselection is required)
                // good rpho for resolved category
                .Define("isGoodRPho",           "(isEBPho || isEEPho) && phoCuts")
                .Define("isGoodMPho",           "isGoodRPho && isHggPho")
                .Define("RphoIdx1",             [](Vec<int> good, Vec<float> pt){if (ROOT::VecOps::Sum(good) > 0) return utils::getIdx(good, pt)[0]; else return -1;}, {"isGoodRPho", "phoCalibEt_NEW"})
                .Define("MphoIdx1",             [](Vec<int> good, Vec<float> pt){if (ROOT::VecOps::Sum(good) > 0) return utils::getIdx(good, pt)[0]; else return -1;}, {"isGoodMPho", "phoCalibEt_NEW"});
    return nf;
}


ROOT::RDF::RNode FindGoodEle(ROOT::RDF::RNode df, const YAML::Node cfg, bool isMC, std::string vary){
    // merged ID prediction
    auto MergedID_ph = cfg["external_files"]["mergedID_model"].as<std::map<std::string, std::string>>();
    auto MergedID_wp = cfg["external_files"]["mergedID_wp"].as<std::map<std::string, float>>();
    auto ff = MergedIDPred(df, "eleClass", MergedID_ph);
    
    ff = ff.Define("eleCalibEt",           "eleCalibEn/cosh(eleEta)");

    // systematic variation of electron energy
    if (!isMC)                            ff = ff.Define("eleHDALRegPt_NEW", "eleCalibPt").Define("eleCalibPt_NEW", "eleCalibPt"); // data
    if (isMC && vary == "Nominal")        ff = ff.Define("eleHDALRegPt_NEW", "eleCalibPt").Define("eleCalibPt_NEW", "eleCalibPt"); // nominal mc
    // for resolved electron and merged 1gsf electron
    if (isMC && vary == "EleScaleStatUp") ff = ff.Define("eleCalibPt_NEW", "eleScale_stat_up/cosh(eleEta)").Define("eleCalibEt_NEW", "eleCalibEt");
    if (isMC && vary == "EleScaleSystUp") ff = ff.Define("eleCalibPt_NEW", "eleScale_syst_up/cosh(eleEta)").Define("eleCalibEt_NEW", "eleCalibEt");
    if (isMC && vary == "EleScaleGainUp") ff = ff.Define("eleCalibPt_NEW", "eleScale_gain_up/cosh(eleEta)").Define("eleCalibEt_NEW", "eleCalibEt");
    if (isMC && vary == "EleScaleStatDo") ff = ff.Define("eleCalibPt_NEW", "eleScale_stat_dn/cosh(eleEta)").Define("eleCalibEt_NEW", "eleCalibEt");
    if (isMC && vary == "EleScaleSystDo") ff = ff.Define("eleCalibPt_NEW", "eleScale_syst_dn/cosh(eleEta)").Define("eleCalibEt_NEW", "eleCalibEt");
    if (isMC && vary == "EleScaleGainDo") ff = ff.Define("eleCalibPt_NEW", "eleScale_gain_dn/cosh(eleEta)").Define("eleCalibEt_NEW", "eleCalibEt");
    if (isMC && vary == "EleSigmaPhiUp")  ff = ff.Define("eleCalibPt_NEW", "eleResol_phi_up/cosh(eleEta)") .Define("eleCalibEt_NEW", "eleCalibEt");
    if (isMC && vary == "EleSigmaRhoUp")  ff = ff.Define("eleCalibPt_NEW", "eleResol_rho_up/cosh(eleEta)") .Define("eleCalibEt_NEW", "eleCalibEt");
    if (isMC && vary == "EleSigmaRhoDo")  ff = ff.Define("eleCalibPt_NEW", "eleResol_rho_dn/cosh(eleEta)") .Define("eleCalibEt_NEW", "eleCalibEt");
    // for merged 2gsf electron
    // if (isMC && vary == "EleHDALScaleUp") ff = ff.Define("eleCalibEt_NEW", "eleHDALCalibEt * (1 + eleHDALScale)").Define("eleCalibPt_NEW", "eleCalibPt");
    // if (isMC && vary == "EleHDALScaleDo") ff = ff.Define("eleCalibEt_NEW", "eleHDALCalibEt * (1 - eleHDALScale)").Define("eleCalibPt_NEW", "eleCalibPt");
    // if (isMC && vary == "EleHDALSmearUp") ff = ff.Define("eleCalibEt_NEW", "eleHDALCalibEt_2d[1]").Define("eleCalibPt_NEW", "eleCalibPt");
    // if (isMC && vary == "EleHDALSmearDo") ff = ff.Define("eleCalibEt_NEW", "eleHDALCalibEt_2d[2]").Define("eleCalibPt_NEW", "eleCalibPt");

    // good electron selections
    std::string isM2GoodEle = utils::joinAND({
        Form("eleClass == 0 && ((isEBEle && eleXGBID > %f) || (isEEEle && eleXGBID > %f))", MergedID_wp["M2EB"], MergedID_wp["M2EE"]),
        "nGsfMatchToReco >= 2",
        "isHggEle",
        "eleConvVeto == 1 && eleTrkMissHits < 1 && eleSubTrkMissHits < 1"
    });

    std::string isM1GoodEle = utils::joinAND({
        Form("eleClass == 0 && ((isEBEle && eleXGBID > %f) || (isEEEle && eleXGBID > %f))", MergedID_wp["M1EB"], MergedID_wp["M1EE"]),
        "nGsfMatchToReco == 1",
        "isHggEle",
        "eleConvVeto == 1 && eleTrkMissHits < 1"
    });
    auto nf = ff.Define("isEBEle",              "abs(eleSCEta) < 1.4442")
                .Define("isEEEle",              "abs(eleSCEta) > 1.566 && abs(eleSCEta) < 2.5")
                .Define("isHggEle",             eleSel::HggPresel,          {"nEle", "eleSCEta", "eleSCPhi", "nPho", "rhoAll", "phoSCEta", "phoSCPhi", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"})
                .Define("isMVAEle",             eleSel::Fall17V2ID,         {"nEle", "elePt", "eleSCEta", "eleIDMVAIso"})
                .Define("isM2GoodEle",          isM2GoodEle)
                .Define("isM1GoodEle",          isM1GoodEle)
                .Define("isReGoodEle",          "eleClass >= 1 && isMVAEle && (isEBEle || isEEEle)")
                .Define("m2Idx1",               [](Vec<int> good, Vec<float> pt){if (ROOT::VecOps::Sum(good) > 0) return utils::getIdx(good, pt)[0]; else return -1;}, {"isM2GoodEle", "eleHDALRegPt_NEW"})
                .Define("m1Idx1",               [](Vec<int> good, Vec<float> pt){if (ROOT::VecOps::Sum(good) > 0) return utils::getIdx(good, pt)[0]; else return -1;}, {"isM1GoodEle", "eleCalibPt_NEW"})
                .Define("reIdx1",               [](Vec<int> good, Vec<float> pt){if (ROOT::VecOps::Sum(good) > 1) return utils::getIdx(good, pt)[0]; else return -1;}, {"isReGoodEle", "eleCalibPt_NEW"})
                .Define("reIdx2",               [](Vec<int> good, Vec<float> pt){if (ROOT::VecOps::Sum(good) > 1) return utils::getIdx(good, pt)[1]; else return -1;}, {"isReGoodEle", "eleCalibPt_NEW"});
    return nf;
}


ROOT::RDF::RNode FilterGoodEvents(ROOT::RDF::RNode df){

    // event selections for merged-2gsf
    std::string evtSelM2 = utils::joinAND({
        "ROOT::VecOps::Sum(isM2GoodEle) > 0", // at least one good merged-2gsf electron
        "ROOT::VecOps::Sum(isGoodMPho) > 0", // at least one good photon
        "isDiPhoHLT",
        "TMath::Max(eleHDALRegPt_NEW[m2Idx1], phoCalibEt_NEW[MphoIdx1]) > 35.",
        "TMath::Min(eleHDALRegPt_NEW[m2Idx1], phoCalibEt_NEW[MphoIdx1]) > 25."
    });

    // event selections for merged-1gsf
    std::string evtSelM1 = utils::joinAND({
        "isM2 != 1",
        "ROOT::VecOps::Sum(isM1GoodEle) > 0", // at least one good merged-1gsf electron
        "ROOT::VecOps::Sum(isGoodMPho) > 0", // at least one good photon
        "isDiPhoHLT",
        "TMath::Max(eleCalibPt_NEW[m1Idx1], phoCalibEt_NEW[MphoIdx1]) > 35.",
        "TMath::Min(eleCalibPt_NEW[m1Idx1], phoCalibEt_NEW[MphoIdx1]) > 25."
    });

    // event selections for resolved
    std::string evtSelRe = utils::joinAND({
        "isM2 != 1",
        "isM1 != 1",
        "ROOT::VecOps::Sum(isReGoodEle) > 1", // at least two good resolved electrons
        "ROOT::VecOps::Sum(isGoodRPho) > 0", // at least one good photon
        "isDiEleHLT",
        "eleCalibPt_NEW[reIdx1] > 25",
        "eleCalibPt_NEW[reIdx2] > 15"
    });

    auto nf = df.Define("isM2",                 evtSelM2)
                .Define("isM1",                 evtSelM1)
                .Define("isRe",                 evtSelRe)

                .Filter("isM2 || isM1 || isRe", "event with good leps")

                // kinematic selections: construct physics object
                .Define("eleIdx1",              "if (isM2) return m2Idx1; else if (isM1) return m1Idx1; else return reIdx1;")
                .Define("eleIdx2",              "if (isM2) return m2Idx1; else if (isM1) return m1Idx1; else return reIdx2;")
                .Define("phoIdx1",              "if (isRe) return RphoIdx1; else return MphoIdx1")

                .Define("eleP4Pt",              "if (isM2) return eleHDALRegPt_NEW; else return eleCalibPt_NEW;")

                .Define("ele1",                 "ROOT::Math::PtEtaPhiMVector v(eleP4Pt[eleIdx1], eleEta[eleIdx1], elePhi[eleIdx1], M_ELE); return v;")
                .Define("ele2",                 "ROOT::Math::PtEtaPhiMVector v(eleP4Pt[eleIdx2], eleEta[eleIdx2], elePhi[eleIdx2], M_ELE); return v;")
                .Define("pho",                  "ROOT::Math::PtEtaPhiMVector v(phoCalibEt_NEW[phoIdx1], phoEta[phoIdx1], phoPhi[phoIdx1], 0.); return v;")

                // only for merged electrons
                .Define("gsf1",                 "ROOT::Math::PtEtaPhiMVector v(eleTrkPt[eleIdx1], eleTrkEta[eleIdx1], eleTrkPhi[eleIdx1], M_ELE); return v;")
                .Define("gsf2",                 "ROOT::Math::PtEtaPhiMVector v(eleSubTrkPt[eleIdx1], eleSubTrkEta[eleIdx1], eleSubTrkPhi[eleIdx1], M_ELE); return v;")
                .Define("diEle",                "if (isRe) return (ele1 + ele2); else return (gsf1 + gsf2);")
                .Define("H",                    "if (isRe) return (ele1 + ele2 + pho); else return (ele1 + pho);")
                .Define("Hgsf",                 "if (isRe) return H; else return (gsf1 + gsf2 + pho);")

                .Filter("if (isM2) return (gsf1.Pt()+gsf2.Pt()) > 44; else return true;", "gsfPt > 44 GeV")
                .Filter("if (isM1) return true; else return diEle.M() < 50;", "Mee < 50 GeV")
                .Filter("if (isRe) return (diEle.Pt()/H.M()) > 0.3 && (pho.Pt()/H.M()) > 0.3; else return (ele1.Pt()/H.M()) > 0.3 && (pho.Pt()/H.M()) > 0.3", "pt mass ratio cut")
                .Filter("ROOT::VecOps::DeltaR(ele1.Eta(), pho.Eta(), ele1.Phi(), pho.Phi()) > 1 && ROOT::VecOps::DeltaR(ele2.Eta(), pho.Eta(), ele2.Phi(), pho.Phi()) > 1", "dR cut")
                .Filter("H.M() > 110. && H.M() < 170.", "three body mass cut");
    return nf;
}


// DiJet tagged (VBF tagged) -> events with two jets passing:
//*     High VBF: Mass(jet1, jet2) > 500, Low VBF: 360 < Mass(jet1, jet2) <= 500
//*     1) Loose Jet ID, jetPt > 30, |jetEta| < 4.7 (per jet sel)
//*     2) dR(jet, leps) > 0.4, dR(jet, pho) > 0.4 (per jet sel)
//*     3) |dEta(jet1, jet2)| > 3.5
//*     4) Eta(uug) - ((jetEta1 + jetEta2) * 0.5) < 2.5 (zepen)
//*     5) |dPhi(uug, jj)| > 2.4
// Boosted tagged: pT(uug) > 60 GeV
// Untagged: EBHR9, EBLR9, EE
ROOT::RDF::RNode CatGoodEvents(ROOT::RDF::RNode df, std::string era, bool isMC, const YAML::Node cfg, std::string vary){
    auto ff = df;

    // perform the correction on the photon labeled by "phoIdx1"
    if (isMC){
        bool doCorr = cfg["external_files"]["showershape_corr"]["doCorr"].as<bool>();
        if (doCorr){
            auto ss_path = cfg["external_files"]["showershape_corr"]["phoR9"].as<std::map<std::string, std::string>>();
            ff = doSSCorrections(df, "phoR9CorrFull5x5_Lead", era, ss_path);
        }
        else
            ff = ff.Define("phoR9CorrFull5x5_Lead", "phoR9Full5x5[phoIdx1]");
    }
    else

    // create a branch for data
    if (!isMC)
        ff = ff.Define("phoCorrR9Full5x5_Lead", "phoCorrR9Full5x5[phoIdx1]");

    // systematic variation of jet energy
    if (!isMC)                      ff = ff.Define("jetSmearPt", "jetPt")                         .Define("jetSmearEn", "jetEn"); // data
    if (isMC && vary == "Nominal")  ff = ff.Define("jetSmearPt", "jetP4Smear * jetPt")            .Define("jetSmearEn", "jetP4Smear * jetEn"); // nominal mc
    if (isMC && vary == "JERUp")    ff = ff.Define("jetSmearPt", "jetP4SmearUp * jetPt")          .Define("jetSmearEn", "jetP4SmearUp * jetEn");
    if (isMC && vary == "JERDo")    ff = ff.Define("jetSmearPt", "jetP4SmearDo * jetPt")          .Define("jetSmearEn", "jetP4SmearDo * jetEn");
    if (isMC && vary == "JECUp")    ff = ff.Define("jetSmearPt", "jetP4Smear * jetPt + jetJECUnc").Define("jetSmearEn", "jetP4Smear * jetEn + jetJECUnc");
    if (isMC && vary == "JECDo")    ff = ff.Define("jetSmearPt", "jetP4Smear * jetPt - jetJECUnc").Define("jetSmearEn", "jetP4Smear * jetEn - jetJECUnc");

    std::string diJetSel = utils::joinAND({
        "jetSmearPt > 30.",
        "abs(jetEta) < 4.7",
        "passJetID",
        "dR_ele1_jet > 0.4",
        "dR_ele2_jet > 0.4",
        "dR_pho_jet > 0.4"
    });
    auto nf = ff.Define("dR_ele1_jet",          cat::dRVector,             {"ele1", "jetEta", "jetPhi"})
                .Define("dR_ele2_jet",          cat::dRVector,             {"ele2", "jetEta", "jetPhi"})
                .Define("dR_pho_jet",           cat::dRVector,             {"pho", "jetEta", "jetPhi"})
                .Define("passJetID",            cat::CutBasedJet,          {"nJet", "jetEta", "jetNHF", "jetNEF", "jetNNP", "jetNCH", "jetCHF", "jetCEF", "jetMUF", "era"})
                .Define("isGoodJet",            diJetSel)
                .Define("vbftag",               cat::VBFtag,               {"isGoodJet", "jetSmearPt", "jetEta", "jetPhi", "jetSmearEn", "H"})
                .Define("isHVbf",               "vbftag[0] == 1")
                .Define("isLVbf",               "vbftag[0] == 2")
                .Define("jetIdx1",              "vbftag[1]")
                .Define("jetIdx2",              "vbftag[2]")
                .Define("isBst",                "H.Pt() > 60.")
                .Define("isEBHR9",              "abs(phoSCEta[phoIdx1]) < 1.4442 && phoR9CorrFull5x5_Lead > 0.96")
                .Define("isEBLR9",              "abs(phoSCEta[phoIdx1]) < 1.4442 && phoR9CorrFull5x5_Lead <= 0.96")
                .Define("isEE",                 "abs(phoSCEta[phoIdx1]) > 1.566 && abs(phoSCEta[phoIdx1]) < 2.5")
                .Define("category",             cat::makeCat,               {"isM2", "isM1", "isRe", "isHVbf", "isLVbf", "isBst", "isEE", "isEBHR9", "isEBLR9"});
    return nf;
}


ROOT::RDF::RNode CheckGenerator(ROOT::RDF::RNode df, bool isMC){
    if (!isMC)
        return df;

    // match gen particles from prompt final state, hard process and from higgs
    auto nf = df.Define("genIdx_reco1",         gen::GenMatch,              {"ele1", "nMC", "mcEta", "mcPhi", "mcPID", "mcMomPID", "mcGMomPID", "mcStatusFlag"})
                .Define("genIdx_reco2",         gen::GenMatch,              {"ele2", "nMC", "mcEta", "mcPhi", "mcPID", "mcMomPID", "mcGMomPID", "mcStatusFlag"})
                .Define("nReco1MatchedGen",     "genIdx_reco1.size()") // number of generator particles can be matched to selected reco electron
                .Define("nReco2MatchedGen",     "genIdx_reco2.size()")

                .Define("isTrueM2",             "eleIdx1 == eleIdx2 && nGsfMatchToReco[eleIdx1] >= 2 && nReco1MatchedGen >= 2")
                .Define("isTrueM1",             "eleIdx1 == eleIdx2 && nGsfMatchToReco[eleIdx1] == 1 && nReco1MatchedGen >= 2")
                .Define("isTrueRe",             "eleIdx1 != eleIdx2 && nReco1MatchedGen == 1 && nReco2MatchedGen == 1 && genIdx_reco1[0] != genIdx_reco2[0]");
    return nf;
}


ROOT::RDF::RNode DefineFinalVars(ROOT::RDF::RNode df){
    auto nf = df.Define("CMS_higgs_mass",           "(float) H.M()")

                .Define("eleP4Pt_Lead",             "eleP4Pt[eleIdx1]")
                .Define("eleXGBID_Lead",            "eleXGBID[eleIdx1]")
                .Define("eleCalibPt_Lead",          "eleCalibPt_NEW[eleIdx1]")
                .Define("elePt_Lead",               "elePt[eleIdx1]")
                .Define("eleEta_Lead",              "eleEta[eleIdx1]")
                .Define("elePhi_Lead",              "elePhi[eleIdx1]")
                .Define("eleSCEn_Lead",             "eleSCEn[eleIdx1]")
                .Define("eleSCEta_Lead",            "eleSCEta[eleIdx1]")
                .Define("eleSCPhi_Lead",            "eleSCPhi[eleIdx1]")
                .Define("eleEn_Lead",               "eleEn[eleIdx1]")
                .Define("eleCalibEn_Lead",          "eleCalibEn[eleIdx1]")
                .Define("eleEcalEn_Lead",           "eleEcalEn[eleIdx1]")
                .Define("eleR9Full5x5_Lead",        "eleR9Full5x5[eleIdx1]")
                .Define("eleGsfDeltaR_Lead",        "gsfDeltaR[eleIdx1]")
                .Define("eleGsfPtRatio_Lead",       "gsfPtRatio[eleIdx1]")
                .Define("eleGsfRelPtRatio_Lead",    "gsfRelPtRatio[eleIdx1]")

                .Define("eleP4Pt_subLead",          "eleP4Pt[eleIdx2]")
                .Define("eleXGBID_subLead",         "eleXGBID[eleIdx2]")
                .Define("eleCalibPt_subLead",       "eleCalibPt_NEW[eleIdx2]")
                .Define("elePt_subLead",            "elePt[eleIdx2]")
                .Define("eleEta_subLead",           "eleEta[eleIdx2]")
                .Define("elePhi_subLead",           "elePhi[eleIdx2]")
                .Define("eleSCEn_subLead",          "eleSCEn[eleIdx2]")
                .Define("eleSCEta_subLead",         "eleSCEta[eleIdx2]")
                .Define("eleSCPhi_subLead",         "eleSCPhi[eleIdx2]")
                .Define("eleEn_subLead",            "eleEn[eleIdx2]")
                .Define("eleCalibEn_subLead",       "eleCalibEn[eleIdx2]")
                .Define("eleEcalEn_subLead",        "eleEcalEn[eleIdx2]")
                .Define("eleR9Full5x5_subLead",     "eleR9Full5x5[eleIdx2]")
                .Define("eleGsfDeltaR_subLead",     "gsfDeltaR[eleIdx2]")
                .Define("eleGsfPtRatio_subLead",    "gsfPtRatio[eleIdx2]")
                .Define("eleGsfRelPtRatio_subLead", "gsfRelPtRatio[eleIdx2]")

                .Define("phoCalibEt_Lead",          "phoCalibEt_NEW[phoIdx1]")
                .Define("phoEt_Lead",               "phoEt[phoIdx1]")
                .Define("phoR9Full5x5_Lead",        "phoR9Full5x5[phoIdx1]")
                .Define("phoR9_Lead",               "phoR9[phoIdx1]")
                .Define("phoEta_Lead",              "phoEta[phoIdx1]")
                .Define("phoPhi_Lead",              "phoPhi[phoIdx1]")
                .Define("phoSCEta_Lead",            "phoSCEta[phoIdx1]")
                .Define("phoSCPhi_Lead",            "phoSCPhi[phoIdx1]");
    return nf;
}


void xAna(const std::string inpath, const std::string outpath, std::string era, bool isMC, const YAML::Node cfg, int range){
    TStopwatch time_j;
    time_j.Start();

    printf(" [+] Read ggNtuple from:\n");
    printf("     - %s\n", inpath.c_str());
    ROOT::RDF::RNode df = ROOT::RDataFrame("skimTree", inpath);
    if (range != -1)
        df = df.Range(0, range);

    if (isMC)
        df = wei::ApplyWeight(df, cfg);

    auto df1 = df.Define("era", [&, era](){return era;});
    // auto df2 = FilterMinimum(df1, era, cfg, isMC);
    auto df2 = FindGoodPho(df1, isMC, "Nominal");
    // auto df4 = FindGSFTracks(df3);
    auto df3 = FindGoodEle(df2, cfg, isMC, "Nominal");
    auto df4 = FilterGoodEvents(df3);
    auto df5 = CatGoodEvents(df4, era, isMC, cfg, "Nominal");
    auto df6 = CheckGenerator(df5, isMC);
    auto df7 = DefineFinalVars(df6);
    auto dfFin = df7;
    if (isMC){
        dfFin = wei::GetScaleFactors(dfFin, cfg);
        dfFin = wei::GetFinalWeights(dfFin);
    }
    else
        dfFin = df7.Define("weight", "(float) 1."); // assign as 1 for data

    // book lazy cut flow report and useful information (not yet execute)
    // also book the yield wrt different normalization weights (for systematic yield unc)
    ROOT::RDF::RResultPtr<ULong64_t> tmp_count;
    auto report = dfFin.Report();
    auto M2_reco = dfFin.Filter("isM2").Count();
    auto M1_reco = dfFin.Filter("isM1").Count();
    auto Re_reco = dfFin.Filter("isRe").Count();
    auto M2_gen = (isMC) ? dfFin.Filter("isTrueM2").Count() : tmp_count;
    auto M1_gen = (isMC) ? dfFin.Filter("isTrueM1").Count() : tmp_count;
    auto Re_gen = (isMC) ? dfFin.Filter("isTrueRe").Count() : tmp_count;

    // save the minitree (trigger the event loop)
    // determine the coulumns to save
    std::vector<std::string> defColNames = dfFin.GetDefinedColumnNames();
    std::vector<std::string> Vars = {
        "category", "run", "event", "rho", "weight",
        "CMS_higgs_mass", "isM2", "isM1", "isRe"
    };
    std::string p4Type = "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >";
    for (size_t i = 0; i < defColNames.size(); i++){
        size_t foundLead = defColNames[i].find("_Lead");
        size_t foundSubLead = defColNames[i].find("_subLead");
        bool foundP4 = dfFin.GetColumnType(defColNames[i]) == p4Type;

        if (foundLead != std::string::npos || foundSubLead != std::string::npos || foundP4)
            Vars.push_back(defColNames[i]);
    }
    if (isMC){
        // concatenate Vars and weiVars(defined in wei::ApplyWeight)
        std::vector<std::string> weiVars = {
            "puwei", "puwei_up", "puwei_down",
            "L1ECALPrefire", "L1ECALPrefireUp", "L1ECALPrefireDown",
            "HLTSF", "HLTSFUp", "HLTSFDo",
            "genwei", "mcwei",

            "weight_puweiUp", "weight_puweiDo", "weight_l1pfUp", "weight_l1pfDo", "weight_hltUp", "weight_hltDo",
            "weight_recoEUp", "weight_recoEDo", "weight_eleIDUp", "weight_eleIDDo",
            "weight_phoIDUp", "weight_phoIDDo", "weight_csevUp", "weight_csevDo",

            "isTrueM2", "isTrueM1", "isTrueRe"
        };
        Vars.insert(Vars.end(), weiVars.begin(), weiVars.end());
    }

    printf(" [+] Save miniTree in:\n");
    printf("     - %s\n", outpath.c_str());
    dfFin.Snapshot("miniTree", outpath.c_str(), Vars);

    // print cut flow and useful information
    TString cutflow = utils::printReport(report);
    int m2 = *M2_reco;
    int m1 = *M1_reco;
    int re = *Re_reco;
    TString event_elements(" [+] Number of events:\n");
    event_elements += TString::Format("     - Merged-2Gsf: %d (%.2f%%)\n", m2, (m2*100./(m2+m1+re)));
    event_elements += TString::Format("     - Merged-1Gsf: %d (%.2f%%)\n", m1, (m1*100./(m2+m1+re)));
    event_elements += TString::Format("     - Resolved   : %d (%.2f%%)\n", re, (re*100./(m2+m1+re)));
    printf(event_elements.Data());
    TString purity(" [+] Purity:\n");
    if (isMC){
        purity += TString::Format("     - Merged-2Gsf: %.2f%%\n", *M2_gen*100./m2);
        purity += TString::Format("     - Merged-1Gsf: %.2f%%\n", *M1_gen*100./m1);
        purity += TString::Format("     - Resolved   : %.2f%%\n", *Re_gen*100./re);
        printf(purity.Data());
    }

    // variation
    if (cfg["variation"] && isMC){
        std::vector<std::string> uncVars = {"category", "CMS_higgs_mass", "weight"};
        ROOT::RDF::RSnapshotOptions opts;
        opts.fMode = "update";

        printf(" [+] Process miniTree under different variations\n");
        auto variation = cfg["variation"].as<std::vector<std::string>>();
        for (size_t v = 0; v < variation.size(); v++){
            if (variation[v].find("PhoScale") != std::string::npos || variation[v].find("PhoSigma") != std::string::npos){
                // check if "PhoScale" or "PhoSigma" in variation
                printf("     - calculate the effect of %s\n", variation[v].c_str());
                auto df3_phoShape = FindGoodPho(df2, isMC, variation[v]);
                // auto df4_phoShape = FindGSFTracks(df3_phoShape);
                auto df5_phoShape = FindGoodEle(df3_phoShape, cfg, isMC, "Nominal");
                auto df6_phoShape = FilterGoodEvents(df5_phoShape);
                auto df7_phoShape = CatGoodEvents(df6_phoShape, era, isMC, cfg, "Nominal");
                auto df8_phoShape = CheckGenerator(df7_phoShape, isMC);
                auto df9_phoShape = DefineFinalVars(df8_phoShape);
                auto dfFin_phoShape = df9_phoShape;
                dfFin_phoShape = wei::GetScaleFactors(dfFin_phoShape, cfg);
                dfFin_phoShape = wei::GetFinalWeights(dfFin_phoShape);
                dfFin_phoShape.Snapshot(Form("miniTree_%s", variation[v].c_str()), outpath.c_str(), uncVars, opts);
            }
            if (variation[v].find("Ele") != std::string::npos){
                printf("     - calculate the effect of %s\n", variation[v].c_str());
                auto df3_eleShape = FindGoodPho(df2, isMC, "Nominal");
                // auto df4_eleShape = FindGSFTracks(df3_eleShape);
                auto df5_eleShape = FindGoodEle(df3_eleShape, cfg, isMC, variation[v]);
                auto df6_eleShape = FilterGoodEvents(df5_eleShape);
                auto df7_eleShape = CatGoodEvents(df6_eleShape, era, isMC, cfg, "Nominal");
                auto df8_eleShape = CheckGenerator(df7_eleShape, isMC);
                auto df9_eleShape = DefineFinalVars(df8_eleShape);
                auto dfFin_eleShape = df9_eleShape;
                dfFin_eleShape = wei::GetScaleFactors(dfFin_eleShape, cfg);
                dfFin_eleShape = wei::GetFinalWeights(dfFin_eleShape);
                dfFin_eleShape.Snapshot(Form("miniTree_%s", variation[v].c_str()), outpath.c_str(), uncVars, opts);
            }
            if (variation[v].find("JER") != std::string::npos || variation[v].find("JEC") != std::string::npos){
                printf("     - calculate the effect of %s\n", variation[v].c_str());
                auto df7_Jet = CatGoodEvents(df6, era, isMC, cfg, variation[v]);
                auto df8_Jet = CheckGenerator(df7_Jet, isMC);
                auto df9_Jet = DefineFinalVars(df8_Jet);
                auto dfFin_Jet = df9_Jet;
                dfFin_Jet = wei::GetScaleFactors(dfFin_Jet, cfg);
                dfFin_Jet = wei::GetFinalWeights(dfFin_Jet);
                dfFin_Jet.Snapshot(Form("miniTree_%s", variation[v].c_str()), outpath.c_str(), uncVars, opts);
            }
        }
    }

    // save the cutflow and purity to the root file
    std::unique_ptr<TFile> fout(TFile::Open(outpath.c_str(), "UPDATE"));
    fout->WriteObject(&cutflow, "cutflow");
    fout->WriteObject(&event_elements, "event_elements");
    if (isMC)
        fout->WriteObject(&purity, "purity");

    time_j.Stop();
    auto htime_j = utils::GetHumanTime(time_j.RealTime());
    printf(" [+] Job's done:\n");
    printf("     - Total time: %d hours, %d mins, %d secs\n", htime_j[0], htime_j[1], htime_j[2]);
}


int main(int argc, char** argv){
    TStopwatch time;
    time.Start();

    std::string config_path;
    int range = -1;
    po::options_description desc{"Options"};
    desc.add_options()
        ("help,h",                                  "Higgs Dalitz decay analysis script")
        ("config,c",    po::value(&config_path),    "configuration file (required)")
        ("range,r",     po::value(&range),          "maximal number of events to process; -1 means all; turn off the MT");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        return 0;
    }
    if (vm.count("config") < 1)
        throw std::invalid_argument("the option '--config' is required but missing");

    //=======================================================//
    //    Load config file and setup the common parameters   //
    //=======================================================//
    const YAML::Node cfg = YAML::LoadFile(config_path);
    auto readpath = cfg["skimTree_path"].as<std::vector<std::string>>();
    auto savepath = cfg["miniTree_path"].as<std::vector<std::string>>();
    if (readpath.size() != savepath.size())
        throw std::runtime_error("Number of skimTree_path != number of miniTree_path");

    // check the directory to save the miniTree exist or not. if not, create it
    std::string direc = fs::path(savepath[0]).parent_path();
    if (!fs::exists(direc)){
        system(Form("mkdir -p %s", direc.c_str()));
    }
    auto nthreads = cfg["threads"].as<int>();
    auto isMC = cfg["isMC"].as<bool>();
    auto era = cfg["era"].as<std::string>();
    auto cfg_extfiles = cfg["external_files"];
    std::string range_text = (range == -1) ? "(-1 means all)" : "";
    std::string thread_txt = (range == -1) ? Form("%d", nthreads) : "1 (turn off MT if range is specified)";
    printf("\n******************************************************\n");
    printf("                  Common parameters                   \n");
    printf("******************************************************\n");
    printf("config files: %s\n", config_path.c_str());
    printf("era: %s\n", era.c_str());
    printf("isMC: %s\n", isMC ? "true" : "false");
    printf("nthreads: %s\n", thread_txt.c_str());
    printf("max events: %d %s\n", range, range_text.c_str());
    printf("******************************************************\n");
    printf("\n");

    //=======================================================//
    //   Run xAna over all of the ntuples list in datasets   //
    //=======================================================//
    cprintf("Process the following ntuples sequentially\n");
    for (size_t i = 0; i < readpath.size(); i++){
        cprintf(Form(" %ld. %s\n", i+1, readpath[i].c_str()));
    }
    printf("\n");

    if (range == -1 && nthreads != 1)
        ROOT::EnableImplicitMT(nthreads);
    for (size_t i = 0; i < readpath.size(); i++){
        printf("*************************  %ld  ************************\n", i+1);
        xAna(readpath[i], savepath[i], era, isMC, cfg, range);
    }

    time.Stop();
    auto htime = utils::GetHumanTime(time.RealTime());
    printf("******************************************************\n");
    printf("All done: %d hours, %d mins, %d secs\n", htime[0], htime[1], htime[2]);
    printf("\n");

    return 0;
}