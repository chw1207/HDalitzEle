#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <Math/Vector4D.h>
#include <TStopwatch.h>
#include <TRandom3.h>

#include "PUWeightCalculator.h"
#include "GsfTracks.h"
#include "Utilities.h"
#include "ElectronSel.h" // hgg preselection

#include "./interface/help.h"

template <typename T>
using Vec = const ROOT::RVec<T>&;
using namespace ROOT::VecOps;


std::map<std::string, float> xs = {
    {"UL2016preVFP",    6404. * 1000},
    {"UL2016postVFP",   6404. * 1000},
    {"UL2017",          6435. * 1000},
    {"UL2018",          6529. * 1000},
};


std::map<std::string, float> lumis = {
    {"UL2016preVFP",     19.52},
    {"UL2016postVFP",    16.81},
    {"UL2017",           41.48},
    {"UL2018",           59.82}
};


ROOT::RDF::RNode AddWeights(ROOT::RDF::RNode df, const std::string era){
    auto all = df.Count();
    auto pos = df.Filter("genWeight > 0",   "positive event").Count();
    auto neg = df.Filter("genWeight <= 0",  "negative event").Count();
    const int totalev = pos.GetValue() - neg.GetValue();
    const float procXS = xs.at(era);
    const float instwei = 1.;
    const float luminosity = lumis.at(era);
    const float mcwei = ((procXS * luminosity) / totalev);
    fmt::print("[INFO] Number of events with genwei = {}, without genwei = {}\n", totalev, all.GetValue());
    fmt::print("[INFO] MC weight = {}\n", mcwei);

    // set up the puwei calculator
    std::vector<std::shared_ptr<PUWeightCalculator>> puCalc = {
        std::make_shared<PUWeightCalculator>(),
        std::make_shared<PUWeightCalculator>(),
        std::make_shared<PUWeightCalculator>()
    };
    if (era == "UL2016preVFP"){
        puCalc[0]->Init("/data4/cmkuo/tools/pileup/UL2016preVFP/PUreweight/PUreweight_13TeV_2016preVFP_GoldenJSON_69200nb.root");
        puCalc[1]->Init("/data4/cmkuo/tools/pileup/UL2016preVFP/PUreweight/PUreweight_13TeV_2016preVFP_GoldenJSON_72383nb.root");
        puCalc[2]->Init("/data4/cmkuo/tools/pileup/UL2016preVFP/PUreweight/PUreweight_13TeV_2016preVFP_GoldenJSON_66016nb.root");
    }
    else if (era == "UL2016postVFP"){
        puCalc[0]->Init("/data4/cmkuo/tools/pileup/UL2016postVFP/PUreweight/PUreweight_13TeV_2016postVFP_GoldenJSON_69200nb.root");
        puCalc[1]->Init("/data4/cmkuo/tools/pileup/UL2016postVFP/PUreweight/PUreweight_13TeV_2016postVFP_GoldenJSON_72383nb.root");
        puCalc[2]->Init("/data4/cmkuo/tools/pileup/UL2016postVFP/PUreweight/PUreweight_13TeV_2016postVFP_GoldenJSON_66016nb.root");
    }
    else if (era == "UL2017"){
        puCalc[0]->Init("/data4/cmkuo/tools/pileup/UL2017/PUreweight/PUreweight_13TeV_2017_GoldenJSON_69200nb.root");
        puCalc[1]->Init("/data4/cmkuo/tools/pileup/UL2017/PUreweight/PUreweight_13TeV_2017_GoldenJSON_72383nb.root");
        puCalc[2]->Init("/data4/cmkuo/tools/pileup/UL2017/PUreweight/PUreweight_13TeV_2017_GoldenJSON_66016nb.root");
    }
    else{
        puCalc[0]->Init("/data4/cmkuo/tools/pileup/UL2018/PUreweight/PUreweight_13TeV_2018_GoldenJSON_69200nb.root");
        puCalc[1]->Init("/data4/cmkuo/tools/pileup/UL2018/PUreweight/PUreweight_13TeV_2018_GoldenJSON_72383nb.root");
        puCalc[2]->Init("/data4/cmkuo/tools/pileup/UL2018/PUreweight/PUreweight_13TeV_2018_GoldenJSON_66016nb.root");
    }

    auto get_pu = [&, puCalc](const int run, const ROOT::RVec<float>& puTrue){
        const float puwei = puCalc[0]->GetWeight(run, puTrue[1]);
        return puwei;
    };
    auto get_pu_up = [&, puCalc](const int run, const ROOT::RVec<float>& puTrue){
        const float puwei_up = puCalc[1]->GetWeight(run, puTrue[1]);
        return puwei_up;
    };
    auto get_pu_do = [&, puCalc](const int run, const ROOT::RVec<float>& puTrue){
        const float puwei_down = puCalc[2]->GetWeight(run, puTrue[1]);
        return puwei_down;
    };

    auto nf = df.Define("mcwei",            [&, mcwei]{return mcwei;})
                .Define("procXS",           [&, procXS]{return procXS;})
                .Define("instwei",          [&, instwei]{return instwei;})
                .Define("puwei",            get_pu,             {"run", "puTrue"})
                .Define("puwei_up",         get_pu_up,          {"run", "puTrue"})
                .Define("puwei_down",       get_pu_do,          {"run", "puTrue"})
                .Define("genwei",           "if (genWeight > 0) return 1.; else return -1.;")
                .Define("poswei",           "puwei * mcwei")
                .Define("wei",              "puwei * mcwei * genwei");
    return nf;
}


ROOT::RDF::RNode FindGen(ROOT::RDF::RNode df){
    auto nf = df.Define("hardProc",         "Helper::GenType(mcStatusFlag, 0)")
                .Define("isPrompt",         "Helper::GenType(mcStatusFlag, 1)")

                // single electron from Z
                .Define("isZEle",           "abs(mcPID) == 11 && mcMomPID == 23")
                .Define("isGenEle",         "isZEle && hardProc && isPrompt && mcEta < 2.5")

                // kick out the events without 2 gen electrons
                .Filter("Sum(isGenEle) == 2", "2 gen ele")

                // gen particle index
                .Define("genEleIdx1",       [](Vec<int> good, Vec<float> pt){return utils::getIdx(good, pt)[0];}, {"isGenEle", "mcPt"})
                .Define("genEleIdx2",       [](Vec<int> good, Vec<float> pt){return utils::getIdx(good, pt)[1];}, {"isGenEle", "mcPt"})

                // P4 of gen particle
                .Define("M_ELE",            "(float) 0.000511")
                .Define("GenEle_Lead",      "ROOT::Math::PtEtaPhiMVector v(mcPt[genEleIdx1], mcEta[genEleIdx1], mcPhi[genEleIdx1], M_ELE); return v;")
                .Define("GenEle_subLead",   "ROOT::Math::PtEtaPhiMVector v(mcPt[genEleIdx2], mcEta[genEleIdx2], mcPhi[genEleIdx2], M_ELE); return v;")
                .Define("diGenEle",         "GenEle_Lead + GenEle_subLead")

                .Define("GenElePtRatio",    "GenEle_subLead.Pt() / GenEle_Lead.Pt()")
                .Define("GenEleDeltaR",     "DeltaR(GenEle_Lead.Eta(), GenEle_subLead.Eta(), GenEle_Lead.Phi(), GenEle_subLead.Phi())")
                .Define("GenElePtSum",      "GenEle_Lead.Pt() + GenEle_subLead.Pt()");

    return nf;
}

ROOT::RVec<float> MatchIndex(
    const int nPar,
    const ROOT::RVec<int>& match_idx,
    const ROOT::RVec<float>& var2Par
){
    ROOT::RVec<float> v(nPar);
    for (int i = 0; i < nPar; i++){
        int idx = match_idx[i];
        v[i] = (idx != -1 && var2Par.size() > 0) ? var2Par[idx] : (float) 0;
    }
    return v;
}

ROOT::RDF::RNode GenToReco(ROOT::RDF::RNode df){
    auto match_pho = [](
        const ROOT::RVec<float>& eleSCEta,
        const ROOT::RVec<float>& phoSCEta
    ){
        ROOT::RVec<int> v(eleSCEta.size());
        for (size_t i = 0; i < eleSCEta.size(); i++){
            int phoIdx = -1;
            for (size_t j = 0; j < phoSCEta.size(); j++){
                if (phoSCEta[j] == eleSCEta[i]){
                    phoIdx = j;
                    break;
                }
            }
            v[i] = phoIdx;
        }
        return v;
    };

    auto nf = df.Define("eleCalibEt",       "eleCalibEn/cosh(eleEta)")
                .Define("recoEleIdx1",      "Helper::RecoIdx(GenEle_Lead, eleCalibEt, eleEta, elePhi, M_ELE)")
                .Define("recoEleIdx2",      "Helper::RecoIdx(GenEle_subLead, eleCalibEt, eleEta, elePhi, M_ELE)")

                // matching condition
                .Define("Gen2Reco2",        "recoEleIdx1 != recoEleIdx2 && recoEleIdx1 != -1 && recoEleIdx2 != -1") // resolved
                .Define("Gen2Reco1",        "recoEleIdx1 == recoEleIdx2 && recoEleIdx1 != -1 && recoEleIdx2 != -1") // merged

                // only event with proper reco particles left
                .Filter("(Gen2Reco2 || Gen2Reco1)", "2 reco ele")

                // P4 of reco particle
                .Define("RecoEle_Lead",     "ROOT::Math::PtEtaPhiMVector v(eleCalibEt[recoEleIdx1], eleEta[recoEleIdx1], elePhi[recoEleIdx1], M_ELE); return v;")
                .Define("RecoEle_subLead",  "ROOT::Math::PtEtaPhiMVector v(eleCalibEt[recoEleIdx2], eleEta[recoEleIdx2], elePhi[recoEleIdx2], M_ELE); return v;")
                .Define("elePresel",        eleSel::HggPresel,          {"nEle", "eleSCEta", "eleSCPhi", "nPho", "rhoAll", "phoSCEta", "phoSCPhi", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"})

                .Define("matchPhoIdx",          match_pho,                  {"eleSCEta", "phoSCEta"})
                .Define("eleE1x3Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE1x3Full5x5)")
                .Define("eleE2x2Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE2x2Full5x5)")
                .Define("eleE2x5Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE2x5Full5x5)")
                .Define("eleE3x3Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE3x3Full5x5)")
                .Define("eleE5x5Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE5x5Full5x5)");
    return nf;
}


ROOT::RDF::RNode RecoToGSF(ROOT::RDF::RNode df){
    auto nf = df.Define("isMainGSF",            gsf::IsMainGSF,             {"event", "nGSFTrk", "gsfD0", "gsfDz", "nEle", "eleD0", "eleDz"})
                .Define("ambGSF",               gsf::TrkEleAssociation,     {"nGSFTrk", "gsfD0", "gsfDz", "nEle", "eleD0", "eleDz", "isMainGSF"})
                .Define("nGsfMatchToReco",      gsf::CalcNGsfMatchToReco,   {"nEle", "ambGSF"})
                .Define("eleTrkIdx",            gsf::FindMainGSF,           {"nEle", "ambGSF"})
                .Define("eleSubTrkIdx",         gsf::FindSubGSF_dRMin,      {"nEle", "ambGSF", "gsfEta", "gsfPhi", "gsfCharge"})
                .Define("eleSubPtTrkIdx",       gsf::FindSubGSF_PtMax,      {"nEle", "ambGSF", "gsfPt", "gsfCharge"})

                .Define("eleTrkPt",             gsf::MatchIndexF,           {"nEle", "eleTrkIdx", "gsfPt"})
                .Define("eleTrkEta",            gsf::MatchIndexF,           {"nEle", "eleTrkIdx", "gsfEta"})
                .Define("eleTrkPhi",            gsf::MatchIndexF,           {"nEle", "eleTrkIdx", "gsfPhi"})
                .Define("eleTrkCharge",         gsf::MatchIndexI,           {"nEle", "eleTrkIdx", "gsfCharge"})
                .Define("eleTrkLayers",         gsf::MatchIndexI,           {"nEle", "eleTrkIdx", "gsfLayers"})
                .Define("eleTrkMissHits",       gsf::MatchIndexI,           {"nEle", "eleTrkIdx", "gsfMissHits"})
                .Define("eleTrkD0",             gsf::MatchIndexF,           {"nEle", "eleTrkIdx", "gsfD0"})
                .Define("eleTrkDz",             gsf::MatchIndexF,           {"nEle", "eleTrkIdx", "gsfDz"})

                .Define("eleSubTrkPt",          gsf::MatchIndexF,           {"nEle", "eleSubTrkIdx", "gsfPt"})
                .Define("eleSubTrkEta",         gsf::MatchIndexF,           {"nEle", "eleSubTrkIdx", "gsfEta"})
                .Define("eleSubTrkPhi",         gsf::MatchIndexF,           {"nEle", "eleSubTrkIdx", "gsfPhi"})
                .Define("eleSubTrkCharge",      gsf::MatchIndexI,           {"nEle", "eleSubTrkIdx", "gsfCharge"})
                .Define("eleSubTrkLayers",      gsf::MatchIndexI,           {"nEle", "eleSubTrkIdx", "gsfLayers"})
                .Define("eleSubTrkMissHits",    gsf::MatchIndexI,           {"nEle", "eleSubTrkIdx", "gsfMissHits"})
                .Define("eleSubTrkD0",          gsf::MatchIndexF,           {"nEle", "eleSubTrkIdx", "gsfD0"})
                .Define("eleSubTrkDz",          gsf::MatchIndexF,           {"nEle", "eleSubTrkIdx", "gsfDz"})

                .Define("eleSubPtTrkPt",        gsf::MatchIndexF,           {"nEle", "eleSubPtTrkIdx", "gsfPt"})
                .Define("eleSubPtTrkEta",       gsf::MatchIndexF,           {"nEle", "eleSubPtTrkIdx", "gsfEta"})
                .Define("eleSubPtTrkPhi",       gsf::MatchIndexF,           {"nEle", "eleSubPtTrkIdx", "gsfPhi"})
                .Define("eleSubPtTrkCharge",    gsf::MatchIndexI,           {"nEle", "eleSubPtTrkIdx", "gsfCharge"})
                .Define("eleSubPtTrkLayers",    gsf::MatchIndexI,           {"nEle", "eleSubPtTrkIdx", "gsfLayers"})
                .Define("eleSubPtTrkMissHits",  gsf::MatchIndexI,           {"nEle", "eleSubPtTrkIdx", "gsfMissHits"})
                .Define("eleSubPtTrkD0",        gsf::MatchIndexF,           {"nEle", "eleSubPtTrkIdx", "gsfD0"})
                .Define("eleSubPtTrkDz",        gsf::MatchIndexF,           {"nEle", "eleSubPtTrkIdx", "gsfDz"})

                .Define("eleTrk1",              utils::P4Vector,            {"eleTrkPt", "eleTrkEta", "eleTrkPhi", "M_ELE"})
                .Define("eleTrk2",              utils::P4Vector,            {"eleSubTrkPt", "eleSubTrkEta", "eleSubTrkPhi", "M_ELE"})
                .Define("gsfPtRatio",           gsf::GetTrkPtRatio,         {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfDeltaR",            gsf::GetTrkdR,              {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfPtSum",             gsf::GetTrkPtSum,           {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfRelPtRatio",        gsf::GetTrkRelPtRatio,      {"nEle", "eleCalibPt", "nGsfMatchToReco", "eleTrk1", "eleTrk2"});
    return nf;
}


ROOT::RDF::RNode DefineFinalVars(ROOT::RDF::RNode df){
    auto nf = df.Define("mcPt_Lead",                        "mcPt[genEleIdx1]")
                .Define("mcEta_Lead",                       "mcEta[genEleIdx1]")
                .Define("mcPhi_Lead",                       "mcPhi[genEleIdx1]")
                .Define("mcVtx_Lead",                       "mcVtx[genEleIdx1]")
                .Define("mcVty_Lead",                       "mcVty[genEleIdx1]")
                .Define("mcVtz_Lead",                       "mcVtz[genEleIdx1]")

                // leading reco variables
                .Define("eleCalibEt_Lead",                  "eleCalibEt[recoEleIdx1]")
                .Define("elePresel_Lead",                   "elePresel[recoEleIdx1]")
                .Define("eleCharge_Lead",                   "eleCharge[recoEleIdx1]")
                .Define("eleChargeConsistent_Lead",         "eleChargeConsistent[recoEleIdx1]")
                .Define("eleConvVeto_Lead",                 "eleConvVeto[recoEleIdx1]")
                .Define("eleEcalDrivenSeed_Lead",           "eleEcalDrivenSeed[recoEleIdx1]")
                .Define("eleEn_Lead",                       "eleEn[recoEleIdx1]")
                .Define("eleSCEn_Lead",                     "eleSCEn[recoEleIdx1]")
                .Define("eleEcalEn_Lead",                   "eleEcalEn[recoEleIdx1]")
                .Define("eleESEnP1_Lead",                   "eleESEnP1[recoEleIdx1]")
                .Define("eleESEnP2_Lead",                   "eleESEnP2[recoEleIdx1]")
                .Define("eleD0_Lead",                       "eleD0[recoEleIdx1]")
                .Define("eleDz_Lead",                       "eleDz[recoEleIdx1]")
                .Define("eleSIP_Lead",                      "eleSIP[recoEleIdx1]")
                .Define("elePt_Lead",                       "elePt[recoEleIdx1]")
                .Define("elePtError_Lead",                  "elePtError[recoEleIdx1]")
                .Define("eleEta_Lead",                      "eleEta[recoEleIdx1]")
                .Define("elePhi_Lead",                      "elePhi[recoEleIdx1]")
                .Define("eleR9_Lead",                       "eleR9[recoEleIdx1]")
                .Define("eleR9Full5x5_Lead",                "eleR9Full5x5[recoEleIdx1]")
                .Define("eleCalibPt_Lead",                  "eleCalibPt[recoEleIdx1]")
                .Define("eleDiffCalibOriPt_Lead",           "(eleCalibPt[recoEleIdx1] - elePt[recoEleIdx1])/elePt[recoEleIdx1]")
                .Define("eleCalibEn_Lead",                  "eleCalibEn[recoEleIdx1]")
                .Define("eleSCEta_Lead",                    "eleSCEta[recoEleIdx1]")
                .Define("eleSCPhi_Lead",                    "eleSCPhi[recoEleIdx1]")
                .Define("eleSCRawEn_Lead",                  "eleSCRawEn[recoEleIdx1]")
                .Define("eleSCEtaWidth_Lead",               "eleSCEtaWidth[recoEleIdx1]")
                .Define("eleSCPhiWidth_Lead",               "eleSCPhiWidth[recoEleIdx1]")
                .Define("eleHoverE_Lead",                   "eleHoverE[recoEleIdx1]")
                .Define("eleEoverP_Lead",                   "eleEoverP[recoEleIdx1]")
                .Define("eleEoverPout_Lead",                "eleEoverPout[recoEleIdx1]")
                .Define("eleEoverPInv_Lead",                "eleEoverPInv[recoEleIdx1]")
                .Define("eleBrem_Lead",                     "eleBrem[recoEleIdx1]")
                .Define("eledEtaAtVtx_Lead",                "eledEtaAtVtx[recoEleIdx1]")
                .Define("eledPhiAtVtx_Lead",                "eledPhiAtVtx[recoEleIdx1]")
                .Define("eleSigmaIEtaIEtaFull5x5_Lead",     "eleSigmaIEtaIEtaFull5x5[recoEleIdx1]")
                .Define("eleSigmaIPhiIPhiFull5x5_Lead",     "eleSigmaIPhiIPhiFull5x5[recoEleIdx1]")
                .Define("eleESEffSigmaRR_Lead",             "eleESEffSigmaRR[recoEleIdx1]")
                .Define("elePFChIso_Lead",                  "elePFChIso[recoEleIdx1]")
                .Define("elePFPhoIso_Lead",                 "elePFPhoIso[recoEleIdx1]")
                .Define("elePFNeuIso_Lead",                 "elePFNeuIso[recoEleIdx1]")
                .Define("elePFPUIso_Lead",                  "elePFPUIso[recoEleIdx1]")
                .Define("elePFClusEcalIso_Lead",            "elePFClusEcalIso[recoEleIdx1]")
                .Define("elePFClusHcalIso_Lead",            "elePFClusHcalIso[recoEleIdx1]")
                .Define("eleIDMVAIso_Lead",                 "eleIDMVAIso[recoEleIdx1]")
                .Define("eleIDMVANoIso_Lead",               "eleIDMVANoIso[recoEleIdx1]")
                .Define("eleTrkdxy_Lead",                   "eleTrkdxy[recoEleIdx1]")
                .Define("eleKFHits_Lead",                   "eleKFHits[recoEleIdx1]")
                .Define("eleKFChi2_Lead",                   "eleKFChi2[recoEleIdx1]")
                .Define("eleGSFChi2_Lead",                  "eleGSFChi2[recoEleIdx1]")
                .Define("eleESEnToRawE_Lead",               "(eleESEnP1[recoEleIdx1] + eleESEnP2[recoEleIdx1])/eleSCRawEn[recoEleIdx1]")
                .Define("nGsfMatchToReco_Lead",             "nGsfMatchToReco[recoEleIdx1]")
                .Define("eleTrk_Lead",                      "ROOT::Math::PtEtaPhiMVector v(eleTrkPt[recoEleIdx1], eleTrkEta[recoEleIdx1], eleTrkPhi[recoEleIdx1], M_ELE); return v;")
                .Define("eleTrkPt_Lead",                    "eleTrkPt[recoEleIdx1]")
                .Define("eleTrkEta_Lead",                   "eleTrkEta[recoEleIdx1]")
                .Define("eleTrkPhi_Lead",                   "eleTrkPhi[recoEleIdx1]")
                .Define("eleTrkCharge_Lead",                "eleTrkCharge[recoEleIdx1]")
                .Define("eleTrkLayers_Lead",                "eleTrkLayers[recoEleIdx1]")
                .Define("eleTrkMissHits_Lead",              "eleTrkMissHits[recoEleIdx1]")
                .Define("eleTrkD0_Lead",                    "eleTrkD0[recoEleIdx1]")
                .Define("eleTrkDz_Lead",                    "eleTrkDz[recoEleIdx1]")
                .Define("eleSubTrk_Lead",                   "ROOT::Math::PtEtaPhiMVector v(eleSubTrkPt[recoEleIdx1], eleSubTrkEta[recoEleIdx1], eleSubTrkPhi[recoEleIdx1], M_ELE); return v;")
                .Define("eleSubTrkPt_Lead",                 "eleSubTrkPt[recoEleIdx1]")
                .Define("eleSubTrkEta_Lead",                "eleSubTrkEta[recoEleIdx1]")
                .Define("eleSubTrkPhi_Lead",                "eleSubTrkPhi[recoEleIdx1]")
                .Define("eleSubTrkCharge_Lead",             "eleSubTrkCharge[recoEleIdx1]")
                .Define("eleSubTrkLayers_Lead",             "eleSubTrkLayers[recoEleIdx1]")
                .Define("eleSubTrkMissHits_Lead",           "eleSubTrkMissHits[recoEleIdx1]")
                .Define("eleSubTrkD0_Lead",                 "eleSubTrkD0[recoEleIdx1]")
                .Define("eleSubTrkDz_Lead",                 "eleSubTrkDz[recoEleIdx1]")
                .Define("eleSubPtTrk_Lead",                 "ROOT::Math::PtEtaPhiMVector v(eleSubPtTrkPt[recoEleIdx1], eleSubPtTrkEta[recoEleIdx1], eleSubPtTrkPhi[recoEleIdx1], M_ELE); return v;")
                .Define("eleSubPtTrkPt_Lead",               "eleSubPtTrkPt[recoEleIdx1]")
                .Define("eleSubPtTrkEta_Lead",              "eleSubPtTrkEta[recoEleIdx1]")
                .Define("eleSubPtTrkPhi_Lead",              "eleSubPtTrkPhi[recoEleIdx1]")
                .Define("eleSubPtTrkCharge_Lead",           "eleSubPtTrkCharge[recoEleIdx1]")
                .Define("eleSubPtTrkLayers_Lead",           "eleSubPtTrkLayers[recoEleIdx1]")
                .Define("eleSubPtTrkMissHits_Lead",         "eleSubPtTrkMissHits[recoEleIdx1]")
                .Define("eleSubPtTrkD0_Lead",               "eleSubPtTrkD0[recoEleIdx1]")
                .Define("eleSubPtTrkDz_Lead",               "eleSubPtTrkDz[recoEleIdx1]")
                .Define("gsfPtSum_Lead",                    "eleTrk_Lead.Pt() + eleSubTrk_Lead.Pt()")
                .Define("gsfPtRatio_Lead",                  "eleSubTrk_Lead.Pt() / eleTrk_Lead.Pt()")
                .Define("gsfDeltaR_Lead",                   "DeltaR(eleTrk_Lead.Eta(), eleSubTrk_Lead.Eta(), eleTrk_Lead.Phi(), eleSubTrk_Lead.Phi())")
                .Define("gsfMissHitsSum_Lead",              "eleTrkMissHits[recoEleIdx1] + eleSubTrkMissHits[recoEleIdx1]")
                .Define("gsfRelPtRatio_Lead",               "if (nGsfMatchToReco_Lead > 1) return (eleTrk_Lead + eleSubTrk_Lead).Pt() / eleSCRawEn_Lead; else return eleTrk_Lead.Pt() / eleSCRawEn_Lead;")

                .Define("eleE1x3Full5x5_Lead",            "eleE1x3Full5x5[recoEleIdx1]")
                .Define("eleE2x2Full5x5_Lead",            "eleE2x2Full5x5[recoEleIdx1]")
                .Define("eleE2x5Full5x5_Lead",            "eleE2x5Full5x5[recoEleIdx1]")
                .Define("eleE3x3Full5x5_Lead",            "eleE3x3Full5x5[recoEleIdx1]")
                .Define("eleE5x5Full5x5_Lead",            "eleE5x5Full5x5[recoEleIdx1]")
                // .Define("eleEmax_Lead",                   "eleEmax[recoEleIdx1]")
                // .Define("eleE2nd_Lead",                   "eleE2nd[recoEleIdx1]")
                // .Define("eleEtop_Lead",                   "eleEtop[recoEleIdx1]")
                // .Define("eleEleft_Lead",                  "eleEleft[recoEleIdx1]")
                // .Define("eleEright_Lead",                 "eleEright[recoEleIdx1]")
                // .Define("eleEbottom_Lead",                "eleEbottom[recoEleIdx1]")
                // .Define("eleIDbit_Lead",                  "eleIDbit[recoEleIdx1]")

                // trailing gen variables
                .Define("mcPt_subLead",                     "mcPt[genEleIdx2]")
                .Define("mcEta_subLead",                    "mcEta[genEleIdx2]")
                .Define("mcPhi_subLead",                    "mcPhi[genEleIdx2]")
                .Define("mcVtx_subLead",                    "mcVtx[genEleIdx2]")
                .Define("mcVty_subLead",                    "mcVty[genEleIdx2]")
                .Define("mcVtz_subLead",                    "mcVtz[genEleIdx2]")

                // trailing reco variables
                .Define("eleCalibEt_subLead",               "eleCalibEt[recoEleIdx2]")
                .Define("elePresel_subLead",                "elePresel[recoEleIdx2]")
                .Define("eleCharge_subLead",                "eleCharge[recoEleIdx2]")
                .Define("eleChargeConsistent_subLead",      "eleChargeConsistent[recoEleIdx2]")
                .Define("eleConvVeto_subLead",              "eleConvVeto[recoEleIdx2]")
                .Define("eleEcalDrivenSeed_subLead",        "eleEcalDrivenSeed[recoEleIdx2]")
                .Define("eleEn_subLead",                    "eleEn[recoEleIdx2]")
                .Define("eleSCEn_subLead",                  "eleSCEn[recoEleIdx2]")
                .Define("eleEcalEn_subLead",                "eleEcalEn[recoEleIdx2]")
                .Define("eleESEnP1_subLead",                "eleESEnP1[recoEleIdx2]")
                .Define("eleESEnP2_subLead",                "eleESEnP2[recoEleIdx2]")
                .Define("eleD0_subLead",                    "eleD0[recoEleIdx2]")
                .Define("eleDz_subLead",                    "eleDz[recoEleIdx2]")
                .Define("eleSIP_subLead",                   "eleSIP[recoEleIdx2]")
                .Define("elePt_subLead",                    "elePt[recoEleIdx2]")
                .Define("elePtError_subLead",               "elePtError[recoEleIdx2]")
                .Define("eleEta_subLead",                   "eleEta[recoEleIdx2]")
                .Define("elePhi_subLead",                   "elePhi[recoEleIdx2]")
                .Define("eleR9_subLead",                    "eleR9[recoEleIdx2]")
                .Define("eleR9Full5x5_subLead",             "eleR9Full5x5[recoEleIdx2]")
                .Define("eleCalibPt_subLead",               "eleCalibPt[recoEleIdx2]")
                .Define("eleDiffCalibOriPt_subLead",        "(eleCalibPt[recoEleIdx2] - elePt[recoEleIdx2])/elePt[recoEleIdx2]")
                .Define("eleCalibEn_subLead",               "eleCalibEn[recoEleIdx2]")
                .Define("eleSCEta_subLead",                 "eleSCEta[recoEleIdx2]")
                .Define("eleSCPhi_subLead",                 "eleSCPhi[recoEleIdx2]")
                .Define("eleSCRawEn_subLead",               "eleSCRawEn[recoEleIdx2]")
                .Define("eleSCEtaWidth_subLead",            "eleSCEtaWidth[recoEleIdx2]")
                .Define("eleSCPhiWidth_subLead",            "eleSCPhiWidth[recoEleIdx2]")
                .Define("eleHoverE_subLead",                "eleHoverE[recoEleIdx2]")
                .Define("eleEoverP_subLead",                "eleEoverP[recoEleIdx2]")
                .Define("eleEoverPout_subLead",             "eleEoverPout[recoEleIdx2]")
                .Define("eleEoverPInv_subLead",             "eleEoverPInv[recoEleIdx2]")
                .Define("eleBrem_subLead",                  "eleBrem[recoEleIdx2]")
                .Define("eledEtaAtVtx_subLead",             "eledEtaAtVtx[recoEleIdx2]")
                .Define("eledPhiAtVtx_subLead",             "eledPhiAtVtx[recoEleIdx2]")
                .Define("eleSigmaIEtaIEtaFull5x5_subLead",  "eleSigmaIEtaIEtaFull5x5[recoEleIdx2]")
                .Define("eleSigmaIPhiIPhiFull5x5_subLead",  "eleSigmaIPhiIPhiFull5x5[recoEleIdx2]")
                .Define("eleESEffSigmaRR_subLead",          "eleESEffSigmaRR[recoEleIdx2]")
                .Define("elePFChIso_subLead",               "elePFChIso[recoEleIdx2]")
                .Define("elePFPhoIso_subLead",              "elePFPhoIso[recoEleIdx2]")
                .Define("elePFNeuIso_subLead",              "elePFNeuIso[recoEleIdx2]")
                .Define("elePFPUIso_subLead",               "elePFPUIso[recoEleIdx2]")
                .Define("elePFClusEcalIso_subLead",         "elePFClusEcalIso[recoEleIdx2]")
                .Define("elePFClusHcalIso_subLead",         "elePFClusHcalIso[recoEleIdx2]")
                .Define("eleIDMVAIso_subLead",              "eleIDMVAIso[recoEleIdx2]")
                .Define("eleIDMVANoIso_subLead",            "eleIDMVANoIso[recoEleIdx2]")
                .Define("eleTrkdxy_subLead",                "eleTrkdxy[recoEleIdx2]")
                .Define("eleKFHits_subLead",                "eleKFHits[recoEleIdx2]")
                .Define("eleKFChi2_subLead",                "eleKFChi2[recoEleIdx2]")
                .Define("eleGSFChi2_subLead",               "eleGSFChi2[recoEleIdx2]")
                .Define("eleESEnToRawE_subLead",            "(eleESEnP1[recoEleIdx2] + eleESEnP2[recoEleIdx2])/eleSCRawEn[recoEleIdx2]")
                .Define("nGsfMatchToReco_subLead",          "nGsfMatchToReco[recoEleIdx2]")
                .Define("eleTrk_subLead",                   "ROOT::Math::PtEtaPhiMVector v(eleTrkPt[recoEleIdx2], eleTrkEta[recoEleIdx2], eleTrkPhi[recoEleIdx2], M_ELE); return v;")
                .Define("eleTrkPt_subLead",                 "eleTrkPt[recoEleIdx2]")
                .Define("eleTrkEta_subLead",                "eleTrkEta[recoEleIdx2]")
                .Define("eleTrkPhi_subLead",                "eleTrkPhi[recoEleIdx2]")
                .Define("eleTrkCharge_subLead",             "eleTrkCharge[recoEleIdx2]")
                .Define("eleTrkLayers_subLead",             "eleTrkLayers[recoEleIdx2]")
                .Define("eleTrkMissHits_subLead",           "eleTrkMissHits[recoEleIdx2]")
                .Define("eleTrkD0_subLead",                 "eleTrkD0[recoEleIdx2]")
                .Define("eleTrkDz_subLead",                 "eleTrkDz[recoEleIdx2]")
                .Define("eleSubTrk_subLead",                "ROOT::Math::PtEtaPhiMVector v(eleSubTrkPt[recoEleIdx2], eleSubTrkEta[recoEleIdx2], eleSubTrkPhi[recoEleIdx2], M_ELE); return v;")
                .Define("eleSubTrkPt_subLead",              "eleSubTrkPt[recoEleIdx2]")
                .Define("eleSubTrkEta_subLead",             "eleSubTrkEta[recoEleIdx2]")
                .Define("eleSubTrkPhi_subLead",             "eleSubTrkPhi[recoEleIdx2]")
                .Define("eleSubTrkCharge_subLead",          "eleSubTrkCharge[recoEleIdx2]")
                .Define("eleSubTrkLayers_subLead",          "eleSubTrkLayers[recoEleIdx2]")
                .Define("eleSubTrkMissHits_subLead",        "eleSubTrkMissHits[recoEleIdx2]")
                .Define("eleSubTrkD0_subLead",              "eleSubTrkD0[recoEleIdx2]")
                .Define("eleSubTrkDz_subLead",              "eleSubTrkDz[recoEleIdx2]")
                .Define("eleSubPtTrk_subLead",              "ROOT::Math::PtEtaPhiMVector v(eleSubPtTrkPt[recoEleIdx2], eleSubPtTrkEta[recoEleIdx2], eleSubPtTrkPhi[recoEleIdx2], M_ELE); return v;")
                .Define("eleSubPtTrkPt_subLead",            "eleSubPtTrkPt[recoEleIdx2]")
                .Define("eleSubPtTrkEta_subLead",           "eleSubPtTrkEta[recoEleIdx2]")
                .Define("eleSubPtTrkPhi_subLead",           "eleSubPtTrkPhi[recoEleIdx2]")
                .Define("eleSubPtTrkCharge_subLead",        "eleSubPtTrkCharge[recoEleIdx2]")
                .Define("eleSubPtTrkLayers_subLead",        "eleSubPtTrkLayers[recoEleIdx2]")
                .Define("eleSubPtTrkMissHits_subLead",      "eleSubPtTrkMissHits[recoEleIdx2]")
                .Define("eleSubPtTrkD0_subLead",            "eleSubPtTrkD0[recoEleIdx2]")
                .Define("eleSubPtTrkDz_subLead",            "eleSubPtTrkDz[recoEleIdx2]")
                .Define("gsfPtSum_subLead",                 "eleTrk_subLead.Pt() + eleSubTrk_subLead.Pt()")
                .Define("gsfPtRatio_subLead",               "eleSubTrk_subLead.Pt() / eleTrk_subLead.Pt()")
                .Define("gsfDeltaR_subLead",                "DeltaR(eleTrk_subLead.Eta(), eleSubTrk_subLead.Eta(), eleTrk_subLead.Phi(), eleSubTrk_subLead.Phi())")
                .Define("gsfMissHitsSum_subLead",           "eleTrkMissHits[recoEleIdx1] + eleSubTrkMissHits[recoEleIdx1]")
                .Define("gsfRelPtRatio_subLead",            "if (nGsfMatchToReco_subLead > 1) return (eleTrk_subLead + eleSubTrk_subLead).Pt() / eleSCRawEn_subLead; else return eleTrk_subLead.Pt() / eleSCRawEn_subLead;")

                .Define("eleE1x3Full5x5_subLead",            "eleE1x3Full5x5[recoEleIdx2]")
                .Define("eleE2x2Full5x5_subLead",            "eleE2x2Full5x5[recoEleIdx2]")
                .Define("eleE2x5Full5x5_subLead",            "eleE2x5Full5x5[recoEleIdx2]")
                .Define("eleE3x3Full5x5_subLead",            "eleE3x3Full5x5[recoEleIdx2]")
                .Define("eleE5x5Full5x5_subLead",            "eleE5x5Full5x5[recoEleIdx2]")
                // .Define("eleEmax_subLead",                   "eleEmax[recoEleIdx2]")
                // .Define("eleE2nd_subLead",                   "eleE2nd[recoEleIdx2]")
                // .Define("eleEtop_subLead",                   "eleEtop[recoEleIdx2]")
                // .Define("eleEleft_subLead",                  "eleEleft[recoEleIdx2]")
                // .Define("eleEright_subLead",                 "eleEright[recoEleIdx2]")
                // .Define("eleEbottom_subLead",                "eleEbottom[recoEleIdx2]")
                // .Define("eleIDbit_subLead",                  "eleIDbit[recoEleIdx2]")

                // event categorization
                // cat == 1 -> Resolved, cat == 2 -> M2, cat == 3 -> M1
                .Define("category",                         "Helper::make_cat(Gen2Reco2, Gen2Reco1, nGsfMatchToReco_Lead)")
                .Define("zMass",                            "(RecoEle_Lead + RecoEle_subLead).M()");
    return nf;
}


void xAnaGenDYJets(std::string infile, std::string outfile, std::string era){
    TStopwatch time;
    time.Start();

    fmt::print("[INFO] Read_File(): {}\n", infile);
    fmt::print("[INFO] Save_File(): {}\n", outfile);

    //=======================================================//
    // Load TTree into RDataFrame                            //
    //=======================================================//
    ROOT::EnableImplicitMT(25);
    ROOT::RDF::RNode df = ROOT::RDataFrame("ggNtuplizer/EventTree", infile.c_str());

    auto df1 = AddWeights(df, era);
    auto df2 = FindGen(df1);
    auto df3 = GenToReco(df2);
    auto df4 = RecoToGSF(df3);
    auto df5 = DefineFinalVars(df4);
    auto dfFin = df5;

    // book the calculation firt
    auto report = dfFin.Report();
    auto M2 = dfFin.Filter("category == 2").Count();
    auto M1 = dfFin.Filter("category == 3").Count();
    auto Re = dfFin.Filter("category == 1").Count();

    //====================================================//
    // Define the final variables to save to the miniTree //
    //====================================================//
    vector<string> defColNames = dfFin.GetDefinedColumnNames();
    vector<string> Vars;
    for (int i = 0; i < defColNames.size(); i++){
        size_t foundLep1 = defColNames[i].find("_Lead");
        size_t foundLep2 = defColNames[i].find("_subLead");

        if (foundLep1 != string::npos || foundLep2 != string::npos)
            Vars.push_back(defColNames[i]);
        if (dfFin.GetColumnType(defColNames[i]) == "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >")
            Vars.push_back(defColNames[i]);
    }
    vector<string> extralVars = {
        "event", "mcwei", "puwei", "genwei", "procXS", "instwei", "wei", "rho", "rhoAll", "nVtx", "nGoodVtx", "isPVGood",
        "category", "zMass"
    };
    Vars.insert(Vars.end(), extralVars.begin(), extralVars.end());
    dfFin.Snapshot("miniTree", outfile.c_str(), Vars);

    //========================================================//
    // Visualize the selection results (cut flow, event count)//
    //========================================================//
    fmt::print("[INFO] Cut flow: \n");
    report->Print();

    int m2 = *M2;
    int m1 = *M1;
    int re = *Re;
    fmt::print("[INFO] Number of events:\n");
    std::cout << fixed << showpoint;
    std::cout << "Merged-2Gsf: " << setprecision(2) << m2*100./(m2+m1+re) << "% (" << m2 << "/" << m2+m1+re << ")" << std::endl;
    std::cout << "Merged-1Gsf: " << setprecision(2) << m1*100./(m2+m1+re) << "% (" << m1 << "/" << m2+m1+re << ")" << std::endl;
    std::cout << "Resolved:    " << setprecision(2) << re*100./(m2+m1+re) << "% (" << re << "/" << m2+m1+re << ")" << std::endl;

    time.Stop();
    time.Print();
}