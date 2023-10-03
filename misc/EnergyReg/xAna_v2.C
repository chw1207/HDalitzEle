#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <Math/Vector4D.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <boost/algorithm/string/join.hpp> // join
#include <boost/algorithm/string.hpp> // contains
#include "interface/roccor/RoccoR.cc"
#include "interface/PhoSelections.h"
#include "interface/MyRobustScaler.h"

// shared library libHDalitzEle.so
#include "PUWeightCalculator.h"
#include "Utilities.h"  // getIdx, P4Vector
#include "XGBReader.h"
#include "GsfTracks.h"


// helper function to make the RVec of dR between a set of lorentz vectors and a lorentz vector;
template <typename Lorentz_t>
ROOT::RVec<float> RVecDeltaR(const ROOT::RVec<Lorentz_t>& v1, const Lorentz_t& v2){
    auto vec = ROOT::VecOps::Map(v1, [&](Lorentz_t v){
        return (float) ROOT::VecOps::DeltaR(v.Eta(), v2.Eta(), v.Phi(), v2.Phi());
    });
    return vec;
}

void xAna_v2(const std::string infile, const std::string outfile, const std::string ext, const std::string era, const float xs, const float lumi){
    TStopwatch time;
    time.Start();

    //=======================================================//
    // Load TTree into RDataFrame and then do the selections //
    //=======================================================//
    std::cout << "[INFO] Read_File(): " << infile.c_str() << std::endl;
    ROOT::EnableImplicitMT(20);
    ROOT::RDF::RNode df = ROOT::RDataFrame("ggNtuplizer/EventTree", infile.c_str());
    const bool isMC = df.HasColumn("nMC");

    //=======================================================//
    //              Calculate the mc weight                  //
    //=======================================================//
    // set up the puwei calculator
    PUWeightCalculator puCalc[3];
    float mcwei = 1;
    if (isMC){
        // calculate the mc weight
        auto pos = df.Filter("genWeight >  0",  "positive event").Count();
        auto neg = df.Filter("genWeight <= 0",  "negative event").Count();
        const int totalev = pos.GetValue() - neg.GetValue();
        mcwei = ((xs * lumi) / totalev);
        std::cout << "[INFO] Number of events with genwei = " << totalev << std::endl;
        std::cout << "[INFO] MC weight = " << mcwei << std::endl;

        if (era == "UL2016preVFP"){
            puCalc[0].Init("/data4/cmkuo/tools/pileup/UL2016preVFP/PUreweight/PUreweight_13TeV_2016preVFP_GoldenJSON_69200nb.root");
            puCalc[1].Init("/data4/cmkuo/tools/pileup/UL2016preVFP/PUreweight/PUreweight_13TeV_2016preVFP_GoldenJSON_72383nb.root");
            puCalc[2].Init("/data4/cmkuo/tools/pileup/UL2016preVFP/PUreweight/PUreweight_13TeV_2016preVFP_GoldenJSON_66016nb.root");
        }
        else if (era == "UL2016postVFP"){
            puCalc[0].Init("/data4/cmkuo/tools/pileup/UL2016postVFP/PUreweight/PUreweight_13TeV_2016postVFP_GoldenJSON_69200nb.root");
            puCalc[1].Init("/data4/cmkuo/tools/pileup/UL2016postVFP/PUreweight/PUreweight_13TeV_2016postVFP_GoldenJSON_72383nb.root");
            puCalc[2].Init("/data4/cmkuo/tools/pileup/UL2016postVFP/PUreweight/PUreweight_13TeV_2016postVFP_GoldenJSON_66016nb.root");
        }
        else if (era == "UL2017"){
            puCalc[0].Init("/data4/cmkuo/tools/pileup/UL2017/PUreweight/PUreweight_13TeV_2017_GoldenJSON_69200nb.root");
            puCalc[1].Init("/data4/cmkuo/tools/pileup/UL2017/PUreweight/PUreweight_13TeV_2017_GoldenJSON_72383nb.root");
            puCalc[2].Init("/data4/cmkuo/tools/pileup/UL2017/PUreweight/PUreweight_13TeV_2017_GoldenJSON_66016nb.root");
        }
        else{
            puCalc[0].Init("/data4/cmkuo/tools/pileup/UL2018/PUreweight/PUreweight_13TeV_2018_GoldenJSON_69200nb.root");
            puCalc[1].Init("/data4/cmkuo/tools/pileup/UL2018/PUreweight/PUreweight_13TeV_2018_GoldenJSON_72383nb.root");
            puCalc[2].Init("/data4/cmkuo/tools/pileup/UL2018/PUreweight/PUreweight_13TeV_2018_GoldenJSON_66016nb.root");
        }

        // define the columns related to weights
        df = df.Define("mcwei",     [&](){return mcwei;})
               .Define("puwei",     [&](const int run, 
                                        const ROOT::RVec<float>& puTrue){
                                            const float puwei = puCalc[0].GetWeight(run, puTrue[1]);
                                            return puwei;
                                        }, {"run", "puTrue"})
               .Define("puwei_up",  [&](const int run, 
                                        const ROOT::RVec<float>& puTrue){
                                            const float puwei_up = puCalc[1].GetWeight(run, puTrue[1]);
                                            return puwei_up;
                                        }, {"run", "puTrue"})
               .Define("puwei_do",  [&](const int run, 
                                        const ROOT::RVec<float>& puTrue){
                                            const float puwei_do = puCalc[2].GetWeight(run, puTrue[1]);
                                            return puwei_do;
                                        }, {"run", "puTrue"})
               .Define("genwei",    "if (genWeight > 0) return (float) 1.; else return (float) -1.;")
               .Define("poswei",    "(float) puwei * mcwei")
               .Define("plotwei",   "(float) puwei * genwei")
               .Define("wei",       "(float) puwei * mcwei * genwei");
    }

    //=======================================================//
    //        Minimum event selections  (event df: ef)       //
    //=======================================================//
    const std::string trigHLT = (ext == "DoubleMuTrig") ? "((HLTEleMuX >> 14) & 1) == 1 || ((HLTEleMuX >> 15) & 1) == 1 || ((HLTEleMuX >> 41) & 1) == 1 || ((HLTEleMuX >> 42) & 1) == 1"
                                                        : "((HLTEleMuX >> 19) & 1) == 1 || ((HLTEleMuX >> 20) & 1) == 1";
    auto ef = df.Filter(trigHLT, "pass HLT")
                .Define("isDoubleMuTrig",   [&](){return (ext == "DoubleMuTrig") ? (int) 1 : (int) 0;})
                .Filter("isPVGood == 1",    "good Vtx")
                .Filter("(nPho > 0) && (nEle > 0) && (nConv > 0) && (nMu > 1) ", "non zero p");

    //=======================================================//
    //            Muon selections  (muon df: mf)             //
    //=======================================================//
    // rochester muon correction
    std::map<std::string, std::string> rcName;
    rcName["UL2016preVFP"]  = "/data4/chenghan/external/RoccorFiles/RoccoR2016aUL.txt";
    rcName["UL2016postVFP"] = "/data4/chenghan/external/RoccorFiles/RoccoR2016bUL.txt";
    rcName["UL2017"]        = "/data4/chenghan/external/RoccorFiles/RoccoR2017UL.txt";
    rcName["UL2018"]        = "/data4/chenghan/external/RoccorFiles/RoccoR2018UL.txt";
    RoccoR rc(rcName[era]);
    TRandom3 rgen; // random generator

    // trigger thresholds
    const float leadPtCut  = (ext == "DoubleMuTrig") ? 20 : 30;
    const float trailPtCut = 10;

    auto mf = ef.Define("muSF",             [&](const int nMu,
                                                const ROOT::RVec<int>& muCharge,
                                                const ROOT::RVec<int>& muTrkLayers,
                                                const ROOT::RVec<float>& muPt,
                                                const ROOT::RVec<float>& muEta,
                                                const ROOT::RVec<float>& muPhi){
                                                    ROOT::RVec<float> scale_factor(nMu);
                                                    ROOT::RVec<float> scale_factor_up(nMu);
                                                    ROOT::RVec<float> scale_factor_do(nMu);
                                                    
                                                    for (int i = 0; i < nMu; i++){
                                                        float SF = 1., SFerr = 0;
                                                        float rndm = rgen.Rndm();
                                                        
                                                        // compute the correction factor
                                                        if (isMC){
                                                            SF = rc.kSmearMC(muCharge[i], muPt[i], muEta[i], muPhi[i], muTrkLayers[i], rndm);
                                                            SFerr = rc.kSmearMCerror(muCharge[i], muPt[i], muEta[i], muPhi[i], muTrkLayers[i], rndm);
                                                        }    
                                                        else{
                                                            SF = rc.kScaleDT(muCharge[i], muPt[i], muEta[i], muPhi[i]);
                                                            SFerr = rc.kScaleDTerror(muCharge[i], muPt[i], muEta[i], muPhi[i]);
                                                        } 

                                                        scale_factor[i] = SF;
                                                        scale_factor_up[i] = SF + SFerr;
                                                        scale_factor_do[i] = SF - SFerr;
                                                    }

                                                    std::vector<ROOT::RVec<float>> scale_factor_vec = {scale_factor, scale_factor_up, scale_factor_do};
                                                    return scale_factor_vec;
                                                }, {"nMu", "muCharge", "muTrkLayers", "muPt", "muEta", "muPhi"})
                .Define("muCorrPt",         "muPt * muSF[0]")
                .Define("muCorrPt_up",      "muPt * muSF[1]")
                .Define("muCorrPt_do",      "muPt * muSF[2]")
                .Define("isMediumPrompt",   "auto v = ROOT::VecOps::Map(muIDbit, [](int bit){return ((bit >> 2) & 1);}); return v;")
                .Define("isGoodMuon",       Form("abs(muEta) < 2.4 && isMediumPrompt && muCorrPt > %f", trailPtCut))
                .Filter("ROOT::VecOps::Sum(isGoodMuon) > 1", "good mus") // at least two good muons
                .Define("mu_good_Idx",      utils::getIdx,              {"isGoodMuon", "muCorrPt"})
                .Define("mu1Idx",           "mu_good_Idx[0]")
                .Define("mu2Idx",           "mu_good_Idx[1]")
                .Define("mu1",              "ROOT::Math::PtEtaPhiMVector v(muCorrPt[mu1Idx], muEta[mu1Idx], muPhi[mu1Idx], 0.105658); return v;")
                .Define("mu2",              "ROOT::Math::PtEtaPhiMVector v(muCorrPt[mu2Idx], muEta[mu2Idx], muPhi[mu2Idx], 0.105658); return v;")
                .Define("mu1_old",          "ROOT::Math::PtEtaPhiMVector v(muPt[mu1Idx], muEta[mu1Idx], muPhi[mu1Idx], 0.105658); return v;")
                .Define("mu2_old",          "ROOT::Math::PtEtaPhiMVector v(muPt[mu2Idx], muEta[mu2Idx], muPhi[mu2Idx], 0.105658); return v;")
                .Filter(Form("mu1.Pt() > %f", leadPtCut), "mu1 pt cut")
                .Filter("(muCharge[mu1Idx] * muCharge[mu2Idx]) < 0", "+/- charge")
                .Filter("(mu1 + mu2).M() > 35", "dimu mass");

    //=======================================================//
    //            Photon selections  (photon df: pf)         //
    //=======================================================//
    /*  Photon slections
        Hgg preselection:
            https://arxiv.org/pdf/2012.06888.pdf
            https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/Taggers/python/flashggPreselectedDiPhotons_cfi.py#L72
        
        conversion vertex matching:
            https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/MicroAOD/plugins/LegacyVertexSelector.cc#L499-L579
            https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/CommonTools/Egamma/src/ConversionTools.cc#L103-L127
            https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/CommonTools/Egamma/interface/ConversionTools.h#L48-L52
    */
    // set up the scaler and reader for merged electron energy regression
    // MyRobustScaler scaler[2] = {
    //     std::string("/data4/chenghan/external/RegressionV4/XGBRegression_RobustScaling_EGMRegTarget_addTrk_EB/RobustScaler.root"),
    //     std::string("/data4/chenghan/external/RegressionV4/XGBRegression_RobustScaling_EGMRegTarget_addTrk_EE/RobustScaler.root")
    // };
    XGBReader readerEB("/data4/chenghan/external/RegressionFinal/XGBRegression_NoRobustScaling_EGMRegTarget_EB/XGB_Regression.txt");
    XGBReader readerEE("/data4/chenghan/external/RegressionFinal/XGBRegression_NoRobustScaling_EGMRegTarget_EE/XGB_Regression.txt");

    std::vector<std::string> feature_EB_vec = {
        "rho",
        "(float) nVtx",
        "eleSCEta[eleIdx]",
        "eleSCPhi[eleIdx]",
        "eleSCRawEn[eleIdx]",
        "eleCalibPt[eleIdx]",

        "eledEtaAtVtx[eleIdx]",
        "eledPhiAtVtx[eleIdx]",
        "elePtError[eleIdx]",
        "eleHoverE[eleIdx]",
        "eleEoverP[eleIdx]",
        "eleEoverPout[eleIdx]",
        "eleEoverPInv[eleIdx]",

        "eleSCEtaWidth[eleIdx]",
        "eleSCPhiWidth[eleIdx]",
        "eleSigmaIEtaIEtaFull5x5[eleIdx]",
        "eleSigmaIPhiIPhiFull5x5[eleIdx]",
        "eleR9Full5x5[eleIdx]",
        "eleBrem[eleIdx]",

        "gsfPtSum[eleIdx]",
        "gsfPtRatio[eleIdx]",
        "diTrkPt[eleIdx]",
        "gsfDeltaR[eleIdx]"
    };
    std::string feature_EB = boost::algorithm::join(feature_EB_vec, ", ");
    std::string feature_EE = feature_EB + ", eleESEnToRawE[eleIdx]";

    auto pf = mf.Define("isEBPho",          "abs(phoSCEta) < 1.4442")
                .Define("isEEPho",          "abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5")
                .Define("M_PHO",            "(float) 0.")
                .Define("phoP4",            utils::P4Vector,            {"phoCalibEt", "phoEta", "phoPhi", "M_PHO"})
                .Define("isMVAPho",         "(isEBPho && phoIDMVA > -0.02) || (isEEPho && phoIDMVA > -0.26)")
                .Define("isFSR",            "FSRSelec(mu1, mu2, phoP4)")
                .Define("isHggPho",         "HggPreSelection(mu1_old, mu2_old, phoP4, rhoAll, nPho, phoSCEta, phoPFChIso, phoPFPhoIso, phoTrkIsoHollowConeDR03, phoR9Full5x5, phoEt, phoSigmaIEtaIEtaFull5x5, phoHoverE)")    
                .Define("isGoodPho",        "(isEBPho || isEEPho) && isFSR && isHggPho")
                .Filter("ROOT::VecOps::Sum(isGoodPho) > 0", "good pho") // at least one good photon
                .Define("phoIdx",           "GetZPho(mu1, mu2, phoP4, isGoodPho)") // get the photon idx which makes z mass cloest to 91.18
                .Define("pho",              "phoP4[phoIdx]")

                //conversion vetex matching
                .Define("sc",               "ROOT::Math::RhoEtaPhiVectorF v(phoSCE[phoIdx]/cosh(fabs(phoSCEta[phoIdx])), phoSCEta[phoIdx], phoSCPhi[phoIdx]); return v;")
                .Define("vtxsP4",           "ROOT::VecOps::Construct<ROOT::Math::XYZVectorF>(convVtxX, convVtxY, convVtxZ)")
                .Define("vtxPsP4",          "ROOT::VecOps::Construct<ROOT::Math::XYZVectorF>(convFitPairPX, convFitPairPY, convFitPairPZ)")
                .Define("isGoodVtx",        "convNTrks == 2 && convFitProb < 1e-6 && hypot(convFitPairPX, convFitPairPY) > 10.")
                .Define("convIdx",          "GetConv(sc, vtxsP4, vtxPsP4, isGoodVtx)")
                .Filter("convIdx != -1",    "match conv") // require to match a conversion vertex
                .Define("convVtx",          "vtxsP4[convIdx]")
                .Define("convVtxP",         "vtxPsP4[convIdx]")

                // electron matching and track association
                .Define("M_ELE",            "(float) 0.000511")
                .Define("hasMatchedEle",    "(eleSCEta == phoSCEta[phoIdx]) && (eleSCPhi == phoSCPhi[phoIdx])")
                .Filter("ROOT::VecOps::Sum(hasMatchedEle) == 1", "match ele") // require to match a electron
                .Define("eleIdx",           "ROOT::VecOps::Nonzero(hasMatchedEle)[0]")
                .Define("ele",              "ROOT::Math::PtEtaPhiMVector v(elePt[eleIdx], eleEta[eleIdx], elePhi[eleIdx], M_ELE); return v;")
                .Define("ele_calib",        "ROOT::Math::PtEtaPhiMVector v(eleCalibPt[eleIdx], eleEta[eleIdx], elePhi[eleIdx], M_ELE); return v;")
                .Define("eleESEnToRawE",    "(eleESEnP1+eleESEnP2)/eleSCRawEn")
                .Define("isMainGSF",        gsf::IsMainGSF,             {"event", "nGSFTrk", "gsfD0", "gsfDz", "nEle", "eleD0", "eleDz"})
                .Define("ambGSF",           gsf::TrkEleAssociation,     {"nGSFTrk", "gsfD0", "gsfDz", "nEle", "eleD0", "eleDz", "isMainGSF"})
                .Define("nGsfMatchToReco",  gsf::CalcNGsfMatchToReco,   {"nEle", "ambGSF"})
                .Define("eleTrkIdx",        gsf::FindMainGSF,           {"nEle", "ambGSF"})
                .Define("eleSubTrkIdx",     gsf::FindSubGSF_dRMin,      {"nEle", "ambGSF", "gsfEta", "gsfPhi", "gsfCharge"})
                .Define("eleTrkPt",         gsf::MatchIndexF,           {"nEle", "eleTrkIdx", "gsfPt"})
                .Define("eleTrkEta",        gsf::MatchIndexF,           {"nEle", "eleTrkIdx", "gsfEta"})
                .Define("eleTrkPhi",        gsf::MatchIndexF,           {"nEle", "eleTrkIdx", "gsfPhi"})
                .Define("eleTrkCharge",     gsf::MatchIndexI,           {"nEle", "eleTrkIdx", "gsfCharge"})
                .Define("eleSubTrkPt",      gsf::MatchIndexF,           {"nEle", "eleSubTrkIdx", "gsfPt"})
                .Define("eleSubTrkEta",     gsf::MatchIndexF,           {"nEle", "eleSubTrkIdx", "gsfEta"})
                .Define("eleSubTrkPhi",     gsf::MatchIndexF,           {"nEle", "eleSubTrkIdx", "gsfPhi"})
                .Define("eleSubTrkCharge",  gsf::MatchIndexI,           {"nEle", "eleSubTrkIdx", "gsfCharge"})
                .Define("eleTrk1",          utils::P4Vector,            {"eleTrkPt", "eleTrkEta", "eleTrkPhi", "M_ELE"})
                .Define("eleTrk2",          utils::P4Vector,            {"eleSubTrkPt", "eleSubTrkEta", "eleSubTrkPhi", "M_ELE"})
                .Define("diTrkPt",          "auto v = ROOT::VecOps::Map(eleTrk1, eleTrk2, [](ROOT::Math::PtEtaPhiMVector trk1, ROOT::Math::PtEtaPhiMVector trk2){return (float) (trk1+trk2).Pt();}); return v;")
                .Define("gsfPtRatio",       gsf::GetTrkPtRatio,         {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfDeltaR",        gsf::GetTrkdR,              {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfPtSum",         gsf::GetTrkPtSum,           {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfRelPtRatio",    gsf::GetTrkRelPtRatio,      {"nEle", "eleSCRawEn", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})

                // electron energy regression
                .Define("eleRegVals",       Form("std::vector<float> v; if (fabs(eleSCEta[eleIdx]) < 1.479) v = {%s}; else v = {%s}; return v;", feature_EB.c_str(), feature_EE.c_str()))
                .Define("elePt_Lead",       "elePt[eleIdx]")
                .Define("eleSCEta_Lead",    "eleSCEta[eleIdx]")
                .Define("eleHDALRegPt_Lead",[&](const float elePt_Lead,
                                                const float eleSCEta_Lead,
                                                const std::vector<float>& eleRegVals){
                                                    // int iBE = (fabs(eleSCEta_Lead) < 1.479) ? 0 : 1;
                                                    // std::vector<float> scaled_features = scaler[iBE].Transform(eleRegVals);
                                                    // float reg = reader[iBE].Compute(scaled_features)[0];
                                                    float reg = 1.;
                                                    if (fabs(eleSCEta_Lead) < 1.479)
                                                        reg = readerEB.Compute(eleRegVals)[0];
                                                    else
                                                        reg = readerEE.Compute(eleRegVals)[0];
                                                    return elePt_Lead * reg;
                                                }, {"elePt_Lead", "eleSCEta_Lead", "eleRegVals"})
                .Define("ele_reg",          "ROOT::Math::PtEtaPhiMVector v(eleHDALRegPt_Lead, eleEta[eleIdx], elePhi[eleIdx], M_ELE); return v;")
                .Define("pho_reg",          "ROOT::Math::PtEtaPhiMVector v(eleHDALRegPt_Lead, pho.Eta(), pho.Phi(), pho.M()); return v;")
                .Define("Z",                "mu1 + mu2 + pho")
                .Define("Z_reg",            "mu1 + mu2 + pho_reg")
                .Define("dimu",             "mu1 + mu2")
                .Define("Zmass",            "(float) Z.M()")
                .Define("Zmass_reg",        "(float) Z_reg.M()");

    //=======================================================//
    //                   Generater matching                  //
    //=======================================================//
    if (isMC){
        pf = pf.Define("gensP4",           "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(mcPt, mcEta, mcPhi, mcMass)")
               .Define("mcPhoDeltaR",      "RVecDeltaR(gensP4, pho)")
               .Define("genIdx",           "ROOT::VecOps::ArgMin(mcPhoDeltaR)")
               .Define("gen",              "gensP4[genIdx]");
    }

    //=======================================================//
    //         Final columns definition (final df: ff)       //
    //=======================================================//
    auto ff = pf.Define("muCorrPt_up_Lead",     "muCorrPt_up[mu1Idx]")
                .Define("muCorrPt_do_Lead",     "muCorrPt_do[mu1Idx]")
                .Define("muCorrPt_up_subLead",  "muCorrPt_up[mu2Idx]")
                .Define("muCorrPt_do_subLead",  "muCorrPt_do[mu2Idx]")

                .Define("phoEleVeto_Lead",      "phoEleVeto[phoIdx]")
                .Define("phoSCEta_Lead",        "phoSCEta[phoIdx]")
                .Define("phoSCPhi_Lead",        "phoSCPhi[phoIdx]")
                .Define("phoSCE_Lead",          "phoSCE[phoIdx]")
                .Define("isMVAPho_Lead",        "isMVAPho[phoIdx]")
                .Define("isHggPho_Lead",        "isHggPho[phoIdx]")
                .Define("phoCalibE_Lead",       "phoCalibE[phoIdx]")
                .Define("phoCalibEt_Lead",      "phoCalibEt[phoIdx]")
                .Define("isEBPho_Lead",         "isEBPho[phoIdx]")
                .Define("isEEPho_Lead",         "isEEPho[phoIdx]")
                .Define("phoSCEtaWidth_Lead",   "phoSCEtaWidth[phoIdx]")
                .Define("phoSCPhiWidth_Lead",   "phoSCPhiWidth[phoIdx]")
                .Define("phoSCBrem_Lead",       "phoSCBrem[phoIdx]")
                .Define("phoHoverE_Lead",       "phoHoverE[phoIdx]")
                .Define("phoR9Full5x5_Lead",    "phoR9Full5x5[phoIdx]")

                .Define("convVtxRadius_Lead",   "convVtxRadius[convIdx]")
                .Define("convFitPairP_Lead",    "hypot(convFitPairPX[convIdx], convFitPairPY[convIdx], convFitPairPZ[convIdx])")
                .Define("convFitPairPt_Lead",   "hypot(convFitPairPX[convIdx], convFitPairPY[convIdx])")
                .Define("convD0_Lead",          "convD0[convIdx]")
                .Define("convDz_Lead",          "convDz[convIdx]")
                .Define("convL0_Lead",          "convL0[convIdx]")
                .Define("convLz_Lead",          "convLz[convIdx]")

                .Define("eleCalibPt_Lead",      "eleCalibPt[eleIdx]")
                .Define("nGsfMatchToReco_Lead", "nGsfMatchToReco[eleIdx]")
                .Define("eleConvVeto_Lead",     "eleConvVeto[eleIdx]")
                .Define("eleSCPhi_Lead",        "eleSCPhi[eleIdx]")
                .Define("eleSCEn_Lead",         "eleSCEn[eleIdx]")
                .Define("eleCalibEn_Lead",      "eleCalibEn[eleIdx]")
                .Define("eleSCRawEn_Lead",      "eleSCRawEn[eleIdx]")
                .Define("eledEtaAtVtx_Lead",    "eledEtaAtVtx[eleIdx]")
                .Define("eledPhiAtVtx_Lead",    "eledPhiAtVtx[eleIdx]")
                .Define("elePtError_Lead",      "elePtError[eleIdx]")
                .Define("eleHoverE_Lead",       "eleHoverE[eleIdx]")
                .Define("eleEoverP_Lead",       "eleEoverP[eleIdx]")
                .Define("eleEoverPout_Lead",    "eleEoverPout[eleIdx]")
                .Define("eleEoverPInv_Lead",    "eleEoverPInv[eleIdx]")
                .Define("eleSCEtaWidth_Lead",   "eleSCEtaWidth[eleIdx]")
                .Define("eleSCPhiWidth_Lead",   "eleSCPhiWidth[eleIdx]")
                .Define("eleSigmaIEtaIEtaFull5x5_Lead",      "eleSigmaIEtaIEtaFull5x5[eleIdx]")
                .Define("eleSigmaIPhiIPhiFull5x5_Lead",      "eleSigmaIPhiIPhiFull5x5[eleIdx]")
                .Define("eleR9Full5x5_Lead",    "eleR9Full5x5[eleIdx]")
                .Define("eleBrem_Lead",         "eleBrem[eleIdx]")
                .Define("elePFChIso_Lead",      "elePFChIso[eleIdx]")
                .Define("elePFPhoIso_Lead",     "elePFPhoIso[eleIdx]")
                .Define("elePFNeuIso_Lead",     "elePFNeuIso[eleIdx]")
                .Define("gsfPtRatio_Lead",      "gsfPtRatio[eleIdx]")
                .Define("gsfDeltaR_Lead",       "gsfDeltaR[eleIdx]")
                .Define("gsfRelPtRatio_Lead",   "gsfRelPtRatio[eleIdx]")
                .Define("gsfPtSum_Lead",        "gsfPtSum[eleIdx]")
                .Define("eleTrk1_Lead",         "eleTrk1[eleIdx]")
                .Define("eleTrk2_Lead",         "eleTrk2[eleIdx]");
    if (isMC){
        ff = ff.Define("mcFromHardProcess_Lead",    "(mcStatusFlag[genIdx] >> 0) & 1")
               .Define("mcPromptFinalState_Lead",   "(mcStatusFlag[genIdx] >> 1) & 1")
               .Define("mcPID_Lead",                "mcPID[genIdx]")
               .Define("mcMomPID_Lead",             "mcMomPID[genIdx]")
               .Define("mcGMomPID_Lead",            "mcGMomPID[genIdx]")
               .Define("mcPhoDeltaR_Lead",          "ROOT::VecOps::DeltaR(gen.Eta(), pho.Eta(), gen.Phi(), pho.Phi())")
               .Define("mcPhoPtDiff_Lead",          "(gen.Pt() - pho.Pt())/gen.Pt()");
    }

    //=======================================================//
    //                     save the results                  //
    //=======================================================//
    auto report = ff.Report();
    std::vector<std::string> colNames = ff.GetDefinedColumnNames();
    std::vector<std::string> outcolNames = {"event", "Zmass", "Zmass_reg", "isDoubleMuTrig"};;
    for (size_t i = 0; i < colNames.size(); i++){
        bool weiColNames  = boost::algorithm::contains(colNames[i], "wei");
        bool defColNames  = boost::algorithm::contains(colNames[i], "Lead");
        bool lvecColTypes = ff.GetColumnType(colNames[i]) == "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >";
        bool dvecColTypes = ff.GetColumnType(colNames[i]) == "ROOT::Math::DisplacementVector3D<ROOT::Math::CylindricalEta3D<float>,ROOT::Math::DefaultCoordinateSystemTag>"
                            || ff.GetColumnType(colNames[i]) == "ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>";
        if (defColNames || lvecColTypes || dvecColTypes || weiColNames)
            outcolNames.push_back(colNames[i]);
    }

    std::cout << "[INFO] Save_File(): " << outfile.c_str() << std::endl;
    ff.Snapshot("miniTree", outfile, outcolNames);

    std::cout << "[INFO] Cut flow:" << std::endl;
    report->Print();

    //=======================================================//
    //           All done, print the execution time          //
    //=======================================================//
    std::cout << "[INFO] Time taken: " << std::endl;
    time.Stop();
    time.Print();
    std::cout << std::endl;
}