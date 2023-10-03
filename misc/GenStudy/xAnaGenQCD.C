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


std::map<std::string, double> xs = {
    {"HT50to100",       (187300000 * 1000.)},
    {"HT100to200",      (23590000  * 1000.)},
    {"HT200to300",      (1555000   * 1000.)},
    {"HT300to500",      (324500    * 1000.)},
    {"HT500to700",      (30310     * 1000.)},
    {"HT700to1000",     (6444      * 1000.)},
    {"HT1000to1500",    (1127      * 1000.)},
    {"HT1500to2000",    (109.8     * 1000.)},
    {"HT2000toInf",     (21.98     * 1000.)}
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
    const double instwei = procXS/xs.at("HT50to100");
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


ROOT::RDF::RNode RecotoGen(ROOT::RDF::RNode df){
    auto nf = df.Filter("nEle > 0",         "nEle > 0")

                .Define("M_ELE",            "(float) 0.000511")
                .Define("eleCalibEt",       "eleCalibEn/cosh(eleEta)")
                .Define("RecoEle",          utils::P4Vector,     {"eleCalibEt", "eleEta", "elePhi", "M_ELE"})
                .Define("GenIdx",           "Helper::getGenIdxVec(RecoEle, mcPt, mcEta, mcPhi, mcMass, mcPID, mcMomPID)")
                .Define("isGoodEle",        "GenIdx != -1") // can match to a gen electron

                .Filter("Sum(isGoodEle) > 0", "1 good ele")

                .Define("recoEleIdx1",      [](Vec<int> good, Vec<float> pt){return utils::getIdx(good, pt)[0];}, {"isGoodEle", "eleCalibEt"})
                .Define("RecoEle_Lead",     "RecoEle[recoEleIdx1]")
                .Define("genEleIdx1",       "GenIdx[recoEleIdx1]")
                .Define("GenEle_Lead",      "ROOT::Math::PtEtaPhiMVector v(mcPt[genEleIdx1], mcEta[genEleIdx1], mcPhi[genEleIdx1], mcMass[genEleIdx1]); return v;")
                .Define("elePresel",        eleSel::HggPresel,          {"nEle", "eleSCEta", "eleSCPhi", "nPho", "rhoAll", "phoSCEta", "phoSCPhi", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"});
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
                .Define("eleIDbit_Lead",                  "eleIDbit[recoEleIdx1]");
    return nf;
}


void xAnaGenQCD(std::string infile, std::string outfile, std::string era, std::string proc){
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
    auto df2 = RecotoGen(df1);
    auto df3 = RecoToGSF(df2);
    auto df4 = DefineFinalVars(df3);
    auto dfFin = df4;

    // book the calculation firt
    auto report = dfFin.Report();
    // auto text = dfFin.Display({"nMC", "nEle", "GenIdx", "recoEleIdx1", "eleCalibEt"}, 20)->AsString();
    // fmt::print("{}\n", text);

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