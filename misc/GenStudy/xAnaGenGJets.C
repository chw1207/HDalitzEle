#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/ranges.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "Math/Vector4D.h"
#include "TStopwatch.h"
#include "TRandom3.h"

#include "PUWeightCalculator.h"
#include "GsfTracks.h"
#include "Utilities.h"
#include "ElectronSel.h" // hgg preselection
#include "./interface/help.h"

template <typename T>
using Vec = const ROOT::RVec<T>&;
using namespace ROOT::VecOps;


// https://twiki.cern.ch/twiki/bin/viewauth/CMS/SummaryTable1G25ns#Gamma_jets
std::map<std::string, float> xs = { 
    {"HT40to100",        20790 * 1000.},
    {"HT100to200",       9238 * 1000.},
    {"HT200to400",       2305 * 1000.},
    {"HT400to600",       274.4 * 1000.},
    {"HT600toInf",       93.46 * 1000.}
};


std::map<std::string, double> lumis = {
    {"UL2016preVFP",     19.52},
    {"UL2016postVFP",    16.81},
    {"UL2017",           41.48},
    {"UL2018",           59.82}
};


ROOT::RDF::RNode AddWeights(ROOT::RDF::RNode df, const std::string era, const std::string proc){
    // calculate the mc weight
    auto all = df.Count();
    auto pos = df.Filter("genWeight > 0",   "positive event").Count();
    auto neg = df.Filter("genWeight <= 0",  "negative event").Count();
    const int totalev = pos.GetValue() - neg.GetValue();
    const double procXS = xs.at(proc);
    const double instwei = procXS/xs.at("HT40to100");
    const double luminosity = lumis.at(era);
    const double mcwei = ((procXS * luminosity) / totalev);
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


ROOT::RDF::RNode RecoToGen(ROOT::RDF::RNode df){
    // Find reco electron matched to gen photon 
    auto nf = df.Filter("nEle > 0 && nConv > 0", "1pho 1conv")
                .Define("M_ELE",            "(float) 0.000511")
                .Define("eleCalibEt",       "eleCalibEn/cosh(eleEta)")
                .Define("recoEle",          utils::P4Vector,     {"eleCalibEt", "eleEta", "elePhi", "M_ELE"})

                .Define("MatchPhoIdxVec",   "Helper::getGenPhoIdxVec(recoEle, mcPt, mcEta, mcPhi, mcMass, mcPID)")
                .Define("hasMathedGenPho",  "MatchPhoIdxVec != -1")
                .Define("MatchConvIdxVec",  "Helper::FindConvVtxWrtEle(nEle, eleSCEta, eleSCPhi, eleSCEn, nConv, convNTrks, convVtxX, convVtxY, convVtxZ, convFitPairPX, convFitPairPY, convFitPairPZ, convFitProb)")
                .Define("hasMathedConvVtx", "MatchConvIdxVec != -1")

                .Define("isConvEle",        "hasMathedGenPho && hasMathedConvVtx")
                .Filter("Sum(isConvEle) > 0", "1 convEle")

                .Define("recoEleIdx",       [](Vec<int> good, Vec<float> pt){return utils::getIdx(good, pt)[0];}, {"isConvEle", "eleCalibEt"})
                .Define("genPhoIdx",        "MatchPhoIdxVec[recoEleIdx]")
                .Define("convIdx",          "MatchConvIdxVec[recoEleIdx]")
                .Define("RecoEle_Lead",     "recoEle[recoEleIdx]")
                .Define("GenPho_Lead",      "ROOT::Math::PtEtaPhiMVector v(mcPt[genPhoIdx], mcEta[genPhoIdx], mcPhi[genPhoIdx], mcMass[genPhoIdx]); return v;")
                
                .Define("elePresel",        eleSel::HggPresel,          {"nEle", "eleSCEta", "eleSCPhi", "nPho", "rhoAll", "phoSCEta", "phoSCPhi", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"});
                // .Define("target",           "GenPho_Lead.Pt()/RecoEle_Lead.Pt()");

    return nf;
}

// ROOT::RDF::RNode GenToReco(ROOT::RDF::RNode df){
//     auto nf = df.Filter("nEle > 0", "1 ele")
//                 .Define("isGenPho", "mcPID == 22")
//                 .Define("isConvEle",[ ](const ROOT::RVec<int>& isGenPho,
//                                         const ROOT::RVec<float>& mcEta,
//                                         const ROOT::RVec<float>& mcPhi,
//                                         const ROOT::RVec<int>& eleEta,
//                                         const ROOT::RVec<int>& elePhi){
//                                             ROOT::RVec<int> v(eleEta.size());
//                                             for (size_t i = 0; i < eleEta.size(); i++){
//                                                 int matched = 0;
//                                                 for (size_t j = 0; j < mcEta.size(); j++){
//                                                     if (!isGenPh[j])
//                                                         continue;
//                                                     float dr = ROOT::VecOps::DeltaR(eleEta[i], mcEta[j], elePhi[i], mcPhi[j]);
//                                                     if (dr < 0.1)
//                                                         matched = 1;
//                                                 }
//                                                 v[i] = matched;
//                                             }
//                                             return v;
//                                         }, {"isGenPho", "mcEta", "mcPhi", "eleEta", "elePhi"})
//                 .Filter("Sum(isConvEle) > 0", "1 convEle")
//                 .Define("M_ELE",            "(float) 0.000511")
//                 .Define("eleCalibEt",       "eleCalibEn/cosh(eleEta)")
//                 .Define("recoEle",          utils::P4Vector,     {"eleCalibEt", "eleEta", "elePhi", "M_ELE"})
//                 .Define("MatchPhoIdxVec",   "Helper::getGenPhoIdxVec(recoEle, mcPt, mcEta, mcPhi, mcMass, mcPID)")
//                 .Define("hasMathedGenPho",  "MatchPhoIdxVec != -1")
//                 .Define("MatchConvIdxVec",  "Helper::FindConvVtxWrtEle(nEle, eleSCEta, eleSCPhi, eleSCEn, nConv, convNTrks, convVtxX, convVtxY, convVtxZ, convFitPairPX, convFitPairPY, convFitPairPZ, convFitProb)")
//                 .Define("hasMathedConvVtx", "MatchConvIdxVec != -1")
//                 .Define("recoEleIdx",       [](Vec<int> good, Vec<float> pt){return utils::getIdx(good, pt)[0];}, {"isConvEle", "eleCalibEt"})
// }

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


ROOT::RDF::RNode RecoToGSF(ROOT::RDF::RNode df){
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
                .Define("eleTrkConvVeto",       gsf::MatchIndexI,           {"nEle", "eleTrkIdx", "gsfConvVeto"})

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
                .Define("eleSubTrkConvVeto",    gsf::MatchIndexI,           {"nEle", "eleSubTrkIdx", "gsfConvVeto"})


                .Define("eleTrk1",              utils::P4Vector,            {"eleTrkPt", "eleTrkEta", "eleTrkPhi", "M_ELE"})
                .Define("eleTrk2",              utils::P4Vector,            {"eleSubTrkPt", "eleSubTrkEta", "eleSubTrkPhi", "M_ELE"})
                .Define("gsfPtRatio",           gsf::GetTrkPtRatio,         {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfDeltaR",            gsf::GetTrkdR,              {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfPtSum",             gsf::GetTrkPtSum,           {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfRelPtRatio",        gsf::GetTrkRelPtRatio,      {"nEle", "eleCalibPt", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})

                .Define("matchPhoIdx",          match_pho,                  {"eleSCEta", "phoSCEta"})
                .Define("eleE1x3Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE1x3Full5x5)")
                .Define("eleE2x2Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE2x2Full5x5)")
                .Define("eleE2x5Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE2x5Full5x5)")
                .Define("eleE3x3Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE3x3Full5x5)")
                .Define("eleE5x5Full5x5",       "MatchIndex(nEle, matchPhoIdx, phoE5x5Full5x5)");
                // .Define("eleEmax",              gsf::MatchIndexF,           {"nEle", "matchPhoIdx", "phoEmax"})
                // .Define("eleE2nd",              gsf::MatchIndexF,           {"nEle", "matchPhoIdx", "phoE2nd"})
                // .Define("eleEtop",              gsf::MatchIndexF,           {"nEle", "matchPhoIdx", "phoEtop"})
                // .Define("eleEleft",             gsf::MatchIndexF,           {"nEle", "matchPhoIdx", "phoEleft"})
                // .Define("eleEright",            gsf::MatchIndexF,           {"nEle", "matchPhoIdx", "phoEright"})
                // .Define("eleEbottom",           gsf::MatchIndexF,           {"nEle", "matchPhoIdx", "phoEbottom"});
    return nf;
}


ROOT::RDF::RNode DefineFinalVars(ROOT::RDF::RNode df){
    auto nf = df.Define("mcPt_Lead",                        "mcPt[genPhoIdx]")
                .Define("mcEta_Lead",                       "mcEta[genPhoIdx]")
                .Define("mcPhi_Lead",                       "mcPhi[genPhoIdx]")
                .Define("mcVtx_Lead",                       "mcVtx[genPhoIdx]")
                .Define("mcVty_Lead",                       "mcVty[genPhoIdx]")
                .Define("mcVtz_Lead",                       "mcVtz[genPhoIdx]")  
                
                // leading reco variables
                .Define("eleCalibEt_Lead",                  "eleCalibEt[recoEleIdx]")
                .Define("elePresel_Lead",                   "elePresel[recoEleIdx]")
                .Define("eleCharge_Lead",                   "eleCharge[recoEleIdx]")
                .Define("eleChargeConsistent_Lead",         "eleChargeConsistent[recoEleIdx]")
                .Define("eleConvVeto_Lead",                 "eleConvVeto[recoEleIdx]")
                .Define("eleEcalDrivenSeed_Lead",           "eleEcalDrivenSeed[recoEleIdx]")
                .Define("eleEn_Lead",                       "eleEn[recoEleIdx]")
                .Define("eleSCEn_Lead",                     "eleSCEn[recoEleIdx]")
                .Define("eleEcalEn_Lead",                   "eleEcalEn[recoEleIdx]")
                .Define("eleESEnP1_Lead",                   "eleESEnP1[recoEleIdx]")
                .Define("eleESEnP2_Lead",                   "eleESEnP2[recoEleIdx]")
                .Define("eleD0_Lead",                       "eleD0[recoEleIdx]")
                .Define("eleDz_Lead",                       "eleDz[recoEleIdx]")
                .Define("eleSIP_Lead",                      "eleSIP[recoEleIdx]")
                .Define("elePt_Lead",                       "elePt[recoEleIdx]")
                .Define("elePtError_Lead",                  "elePtError[recoEleIdx]")
                .Define("eleEta_Lead",                      "eleEta[recoEleIdx]")
                .Define("elePhi_Lead",                      "elePhi[recoEleIdx]")
                .Define("eleR9_Lead",                       "eleR9[recoEleIdx]")
                .Define("eleR9Full5x5_Lead",                "eleR9Full5x5[recoEleIdx]")
                .Define("eleCalibPt_Lead",                  "eleCalibPt[recoEleIdx]")
                .Define("eleDiffCalibOriPt_Lead",           "(eleCalibPt[recoEleIdx] - elePt[recoEleIdx])/elePt[recoEleIdx]")
                .Define("eleCalibEn_Lead",                  "eleCalibEn[recoEleIdx]")
                .Define("eleSCEta_Lead",                    "eleSCEta[recoEleIdx]")
                .Define("eleSCPhi_Lead",                    "eleSCPhi[recoEleIdx]")
                .Define("eleSCRawEn_Lead",                  "eleSCRawEn[recoEleIdx]")
                .Define("eleSCEtaWidth_Lead",               "eleSCEtaWidth[recoEleIdx]")
                .Define("eleSCPhiWidth_Lead",               "eleSCPhiWidth[recoEleIdx]")
                .Define("eleHoverE_Lead",                   "eleHoverE[recoEleIdx]")
                .Define("eleEoverP_Lead",                   "eleEoverP[recoEleIdx]")
                .Define("eleEoverPout_Lead",                "eleEoverPout[recoEleIdx]")
                .Define("eleEoverPInv_Lead",                "eleEoverPInv[recoEleIdx]")
                .Define("eleBrem_Lead",                     "eleBrem[recoEleIdx]")
                .Define("eledEtaAtVtx_Lead",                "eledEtaAtVtx[recoEleIdx]")
                .Define("eledPhiAtVtx_Lead",                "eledPhiAtVtx[recoEleIdx]")
                .Define("eleSigmaIEtaIEtaFull5x5_Lead",     "eleSigmaIEtaIEtaFull5x5[recoEleIdx]")
                .Define("eleSigmaIPhiIPhiFull5x5_Lead",     "eleSigmaIPhiIPhiFull5x5[recoEleIdx]")
                .Define("eleESEffSigmaRR_Lead",             "eleESEffSigmaRR[recoEleIdx]")
                .Define("elePFChIso_Lead",                  "elePFChIso[recoEleIdx]")
                .Define("elePFPhoIso_Lead",                 "elePFPhoIso[recoEleIdx]")
                .Define("elePFNeuIso_Lead",                 "elePFNeuIso[recoEleIdx]")
                .Define("elePFPUIso_Lead",                  "elePFPUIso[recoEleIdx]")
                .Define("elePFClusEcalIso_Lead",            "elePFClusEcalIso[recoEleIdx]")
                .Define("elePFClusHcalIso_Lead",            "elePFClusHcalIso[recoEleIdx]")
                .Define("eleIDMVAIso_Lead",                 "eleIDMVAIso[recoEleIdx]")
                .Define("eleIDMVANoIso_Lead",               "eleIDMVANoIso[recoEleIdx]")
                .Define("eleTrkdxy_Lead",                   "eleTrkdxy[recoEleIdx]")
                .Define("eleKFHits_Lead",                   "eleKFHits[recoEleIdx]")
                .Define("eleKFChi2_Lead",                   "eleKFChi2[recoEleIdx]")
                .Define("eleGSFChi2_Lead",                  "eleGSFChi2[recoEleIdx]")
                .Define("eleESEnToRawE_Lead",               "(eleESEnP1[recoEleIdx] + eleESEnP2[recoEleIdx])/eleSCRawEn[recoEleIdx]")
                .Define("nGsfMatchToReco_Lead",             "nGsfMatchToReco[recoEleIdx]")
                .Define("eleTrk_Lead",                      "ROOT::Math::PtEtaPhiMVector v(eleTrkPt[recoEleIdx], eleTrkEta[recoEleIdx], eleTrkPhi[recoEleIdx], M_ELE); return v;")
                .Define("eleTrkPt_Lead",                    "eleTrkPt[recoEleIdx]")
                .Define("eleTrkEta_Lead",                   "eleTrkEta[recoEleIdx]")
                .Define("eleTrkPhi_Lead",                   "eleTrkPhi[recoEleIdx]")
                .Define("eleTrkCharge_Lead",                "eleTrkCharge[recoEleIdx]")
                .Define("eleTrkLayers_Lead",                "eleTrkLayers[recoEleIdx]")
                .Define("eleTrkMissHits_Lead",              "eleTrkMissHits[recoEleIdx]")
                .Define("eleTrkD0_Lead",                    "eleTrkD0[recoEleIdx]")
                .Define("eleTrkDz_Lead",                    "eleTrkDz[recoEleIdx]")
                .Define("eleSubTrk_Lead",                   "ROOT::Math::PtEtaPhiMVector v(eleSubTrkPt[recoEleIdx], eleSubTrkEta[recoEleIdx], eleSubTrkPhi[recoEleIdx], M_ELE); return v;")
                .Define("eleSubTrkPt_Lead",                 "eleSubTrkPt[recoEleIdx]")
                .Define("eleSubTrkEta_Lead",                "eleSubTrkEta[recoEleIdx]")
                .Define("eleSubTrkPhi_Lead",                "eleSubTrkPhi[recoEleIdx]")
                .Define("eleSubTrkCharge_Lead",             "eleSubTrkCharge[recoEleIdx]")
                .Define("eleSubTrkLayers_Lead",             "eleSubTrkLayers[recoEleIdx]")
                .Define("eleSubTrkMissHits_Lead",           "eleSubTrkMissHits[recoEleIdx]")
                .Define("eleSubTrkD0_Lead",                 "eleSubTrkD0[recoEleIdx]")
                .Define("eleSubTrkDz_Lead",                 "eleSubTrkDz[recoEleIdx]")
                .Define("eleSubPtTrk_Lead",                 "ROOT::Math::PtEtaPhiMVector v(eleSubPtTrkPt[recoEleIdx], eleSubPtTrkEta[recoEleIdx], eleSubPtTrkPhi[recoEleIdx], M_ELE); return v;")
                .Define("eleSubPtTrkPt_Lead",               "eleSubPtTrkPt[recoEleIdx]")
                .Define("eleSubPtTrkEta_Lead",              "eleSubPtTrkEta[recoEleIdx]")
                .Define("eleSubPtTrkPhi_Lead",              "eleSubPtTrkPhi[recoEleIdx]")
                .Define("eleSubPtTrkCharge_Lead",           "eleSubPtTrkCharge[recoEleIdx]")
                .Define("eleSubPtTrkLayers_Lead",           "eleSubPtTrkLayers[recoEleIdx]")
                .Define("eleSubPtTrkMissHits_Lead",         "eleSubPtTrkMissHits[recoEleIdx]")
                .Define("eleSubPtTrkD0_Lead",               "eleSubPtTrkD0[recoEleIdx]")
                .Define("eleSubPtTrkDz_Lead",               "eleSubPtTrkDz[recoEleIdx]")
                .Define("gsfPtSum_Lead",                    "eleTrk_Lead.Pt() + eleSubTrk_Lead.Pt()")
                .Define("gsfPtRatio_Lead",                  "eleSubTrk_Lead.Pt() / eleTrk_Lead.Pt()")
                .Define("gsfDeltaR_Lead",                   "DeltaR(eleTrk_Lead.Eta(), eleSubTrk_Lead.Eta(), eleTrk_Lead.Phi(), eleSubTrk_Lead.Phi())")
                .Define("gsfMissHitsSum_Lead",              "eleTrkMissHits[recoEleIdx] + eleSubTrkMissHits[recoEleIdx]")
                .Define("gsfRelPtRatio_Lead",               "if (nGsfMatchToReco_Lead > 1) return (eleTrk_Lead + eleSubTrk_Lead).Pt() / eleSCRawEn_Lead; else return eleTrk_Lead.Pt() / eleSCRawEn_Lead;")

                .Define("eleIDbit_Lead",                    "eleIDbit[recoEleIdx]")

                .Define("convVtxRadius_Lead",               "convVtxRadius[convIdx]")
                .Define("convFitPairP_Lead",                "sqrt(pow(convFitPairPX[convIdx], 2) + pow(convFitPairPY[convIdx], 2) + pow(convFitPairPZ[convIdx], 2))")
                .Define("convD0_Lead",                      "convD0[convIdx]")
                .Define("convDz_Lead",                      "convDz[convIdx]")
                .Define("convL0_Lead",                      "convL0[convIdx]")
                .Define("convLz_Lead",                      "convLz[convIdx]")
                
                .Define("eleE1x3Full5x5_Lead",            "eleE1x3Full5x5[recoEleIdx]")
                .Define("eleE2x2Full5x5_Lead",            "eleE2x2Full5x5[recoEleIdx]")
                .Define("eleE2x5Full5x5_Lead",            "eleE2x5Full5x5[recoEleIdx]")
                .Define("eleE3x3Full5x5_Lead",            "eleE3x3Full5x5[recoEleIdx]")
                .Define("eleE5x5Full5x5_Lead",            "eleE5x5Full5x5[recoEleIdx]")
                
                .Define("eleSubTrkConvVeto_Lead",                "eleSubTrkConvVeto[recoEleIdx]")
                .Define("eleTrkConvVeto_Lead",                   "eleTrkConvVeto[recoEleIdx]");
                // .Define("eleEmax_Lead",                   "eleEmax[recoEleIdx]")
                // .Define("eleE2nd_Lead",                   "eleE2nd[recoEleIdx]")
                // .Define("eleEtop_Lead",                   "eleEtop[recoEleIdx]")
                // .Define("eleEleft_Lead",                  "eleEleft[recoEleIdx]")
                // .Define("eleEright_Lead",                 "eleEright[recoEleIdx]")
                // .Define("eleEbottom_Lead",                "eleEbottom[recoEleIdx]")
                // .Define("eleIDbit_Lead",                  "eleIDbit[recoEleIdx]");
    return nf;
}


void xAnaGenGJets(std::string infile, std::string outfile, std::string era, std::string proc){
    TStopwatch time;
    time.Start();

    fmt::print("[INFO] Read_File(): {}\n", infile);
    fmt::print("[INFO] Save_File(): {}\n", outfile);

    //=======================================================//
    // Load TTree into RDataFrame                            //
    //=======================================================//
    ROOT::EnableImplicitMT(25);
    ROOT::RDF::RNode df = ROOT::RDataFrame("ggNtuplizer/EventTree", infile.c_str());
    fmt::print("[INFO] Process: {} \n", proc);

    auto df1 = AddWeights(df, era, proc);
    auto df2 = RecoToGen(df1);
    auto df3 = RecoToGSF(df2);
    auto df4 = DefineFinalVars(df3);
    auto dfFin = df4;

    // book the calculation firt
    auto report = dfFin.Report();

    //====================================================//
    // Define the final variables to save to the miniTree //
    //====================================================//
    vector<string> defColNames = dfFin.GetDefinedColumnNames();
    vector<string> Vars;
    for (int i = 0; i < defColNames.size(); i++){
        size_t foundLep1 = defColNames[i].find("_Lead");
        // size_t foundLep2 = defColNames[i].find("_subLead");

        if (foundLep1 != string::npos)
            Vars.push_back(defColNames[i]);
        if (dfFin.GetColumnType(defColNames[i]) == "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >")
            Vars.push_back(defColNames[i]);
    }
    vector<string> extralVars = {
        "event", "mcwei", "puwei", "genwei", "procXS", "instwei", "wei", "rho", "rhoAll", "nVtx", "nGoodVtx", "isPVGood"
    };
    Vars.insert(Vars.end(), extralVars.begin(), extralVars.end());
    dfFin.Snapshot("miniTree", outfile.c_str(), Vars);
    // dfFin.Snapshot("miniTree", "test_qcd.root", Vars);

    fmt::print("[INFO] Cut flow: \n");
    report->Print();

    time.Stop();
    time.Print();
}