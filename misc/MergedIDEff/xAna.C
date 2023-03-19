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
#include "interface/roccor/RoccoR.cc"
#include "interface/MuSelections.h"
#include "interface/PhoSelections.h"
#include "MergedIDPred.h"
#include "PUWeightCalculator.h"
#include "GsfTracks.h"
#include "Utilities.h"

template <typename T>
using Vec = const ROOT::RVec<T>&;
using namespace ROOT::VecOps;

ROOT::RDF::RNode AddWeights(ROOT::RDF::RNode df, const bool isMC, const std::string era, const float xs, const float lumi){
    if (!isMC)
        return df;

    // calculate the mc weight
    auto pos = df.Filter("genWeight > 0",   "positive event").Count();
    auto neg = df.Filter("genWeight <= 0",  "negative event").Count();
    const int totalev = pos.GetValue() - neg.GetValue();
    const float mcwei = ((xs * lumi) / totalev);
    std::cout << "[INFO] Number of events with genwei = " << totalev << std::endl;
    std::cout << "[INFO] MC weight = " << mcwei << std::endl;

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
                .Define("puwei",            get_pu,             {"run", "puTrue"})
                .Define("puwei_up",         get_pu_up,          {"run", "puTrue"})
                .Define("puwei_down",       get_pu_do,          {"run", "puTrue"})
                .Define("genwei",           "if (genWeight > 0) return 1.; else return -1.;")
                .Define("poswei",           "puwei * mcwei")
                .Define("wei",              "puwei * mcwei * genwei");
    return nf;
}


ROOT::RDF::RNode FilterMinimum(ROOT::RDF::RNode df, const std::string ext){
    const std::string trig2mu = "((HLTEleMuX >> 14) & 1) == 1 || ((HLTEleMuX >> 15) & 1) == 1 || ((HLTEleMuX >> 41) & 1) == 1 || ((HLTEleMuX >> 42) & 1) == 1";
    const std::string trig1mu = "((HLTEleMuX >> 19) & 1) == 1 || ((HLTEleMuX >> 20) & 1) == 1";
    std::string trig = "";

    auto ff = df;
    if (ext == "DoubleMuTrig"){
        std::cout << "[INFO] HLT bit(double muon trigger): 14, 15, 41, 42" << std::endl;
        ff = ff.Filter(trig2mu, "Pass HLT").Define("isDoubleMuTrig", "(int) 1");
    }
    else {
        std::cout << "[INFO] HLT bit(single muon trigger): 19, 20" << std::endl;
        ff = ff.Filter(trig1mu, "Pass HLT").Define("isDoubleMuTrig", "(int) 0");
    }

    auto nf = ff.Filter("isPVGood == 1", "Good Vtx")
                .Filter("(nPho > 0) && (nEle > 0) && (nConv > 0) && (nMu > 1) ", "Nonzero");
    return nf;
}


ROOT::RDF::RNode FindGSFTracks(ROOT::RDF::RNode df) {
    auto nf = df.Define("M_ELE",                "(float) 0.000511")
                .Define("isMainGSF",            gsf::IsMainGSF,             {"event", "nGSFTrk", "gsfD0", "gsfDz", "nEle", "eleD0", "eleDz"})
                .Define("ambGSF",               gsf::TrkEleAssociation,     {"nGSFTrk", "gsfD0", "gsfDz", "nEle", "eleD0", "eleDz", "isMainGSF"})
                .Define("nGsfMatchToReco",      gsf::CalcNGsfMatchToReco,   {"nEle", "ambGSF"})
                .Define("eleTrkIdx",            gsf::FindMainGSF,           {"nEle", "ambGSF"})
                .Define("eleSubTrkIdx",         gsf::FindSubGSF_dRMin,      {"nEle", "ambGSF", "gsfEta", "gsfPhi", "gsfCharge"})

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

                .Define("eleTrk1",              utils::P4Vector,            {"eleTrkPt", "eleTrkEta", "eleTrkPhi", "M_ELE"})
                .Define("eleTrk2",              utils::P4Vector,            {"eleSubTrkPt", "eleSubTrkEta", "eleSubTrkPhi", "M_ELE"})
                .Define("gsfPtRatio",           gsf::GetTrkPtRatio,         {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfDeltaR",            gsf::GetTrkdR,              {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfPtSum",             gsf::GetTrkPtSum,           {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfRelPtRatio",        gsf::GetTrkRelPtRatio,      {"nEle", "eleSCRawEn", "nGsfMatchToReco", "eleTrk1", "eleTrk2"});
    return nf;
}


ROOT::RDF::RNode FindGoodMus(ROOT::RDF::RNode df, const std::string era) {
    std::map<std::string, std::string> rcName;
    rcName["UL2016preVFP"] = "interface/roccor/RoccoR2016aUL.txt";
    rcName["UL2016postVFP"] = "interface/roccor/RoccoR2016bUL.txt";
    rcName["UL2017"] = "interface/roccor/RoccoR2017UL.txt";
    rcName["UL2018"] = "interface/roccor/RoccoR2018UL.txt";

    auto rc = std::make_shared<RoccoR>(rcName[era.c_str()].c_str());
    std::random_device rd; // high quality random seed generator
    gRandom->SetSeed(rd());

    auto getRoccoR = [&, rc](
        const bool isMC,
        const Long64_t event,
        const int nMu,
        Vec<int>& muCharge,
        Vec<float>& muPt,
        Vec<float>& muEta,
        Vec<float>& muPhi,
        Vec<int>& muTrkLayers
    ){
        ROOT::RVec<float> sf;
        sf.clear();

        for (int i = 0; i < nMu; i++){
            if (isMC)
                sf.push_back(rc->kSmearMC(muCharge[i], muPt[i], muEta[i], muPhi[i], muTrkLayers[i], gRandom->Rndm(), 0, 0));
            else
                sf.push_back(rc->kScaleDT(muCharge[i], muPt[i], muEta[i], muPhi[i], 0, 0));
        }

        return sf;
    };

    // events with at least two tag muons passing the medium prompt muon ID are selected.
    // https://iopscience.iop.org/article/10.1088/1748-0221/13/06/P06015/pdf
    auto nf = df.Define("muSF",                 getRoccoR,      {"isMC", "event", "nMu", "muCharge", "muPt", "muEta", "muPhi", "muTrkLayers"})
                .Define("muCorrPt",             "muPt * muSF")
                .Define("isMediumPrompt",       "PassMuonID(muIDbit, 2)")
                .Define("isGoodMuon",           "abs(muEta) < 2.4 && isMediumPrompt")

                .Filter("Sum(isGoodMuon) > 1", "good mu")

                .Define("mu1Idx",               [](Vec<int> good, Vec<float> pt){if (Sum(good) > 1) return utils::getIdx(good, pt)[0]; else return -1;}, {"isGoodMuon", "muCorrPt"})
                .Define("mu2Idx",               [](Vec<int> good, Vec<float> pt){if (Sum(good) > 1) return utils::getIdx(good, pt)[1]; else return -1;}, {"isGoodMuon", "muCorrPt"})
                .Define("mu1",                  "ROOT::Math::PtEtaPhiMVector v(muCorrPt[mu1Idx], muEta[mu1Idx], muPhi[mu1Idx], 105.658*0.001); return v;")
                .Define("mu2",                  "ROOT::Math::PtEtaPhiMVector v(muCorrPt[mu2Idx], muEta[mu2Idx], muPhi[mu2Idx], 105.658*0.001); return v;")
                .Define("mu1_old",              "ROOT::Math::PtEtaPhiMVector v(muPt[mu1Idx], muEta[mu1Idx], muPhi[mu1Idx], 105.658*0.001); return v;")
                .Define("mu2_old",              "ROOT::Math::PtEtaPhiMVector v(muPt[mu2Idx], muEta[mu2Idx], muPhi[mu2Idx], 105.658*0.001); return v;")

                .Filter("(muCharge[mu1Idx] * muCharge[mu2Idx]) < 0", "+/- charge")
                .Filter("(mu1 + mu2).M() > 35", "dimu mass");
    return nf;
}


ROOT::RDF::RNode FindGoodPho(ROOT::RDF::RNode df, std::string era) {
    // Loose MVA ID
    float mva_EB = -0.02, mva_EE = -0.26; 

    auto nf = df.Define("M_PHO",                "(float) 0.")
                .Define("phoP4",                utils::P4Vector,                {"phoCalibEt", "phoEta", "phoPhi", "M_PHO"})

                // .Define("isEBPho",              Form("abs(phoSCEta) < 1.4442 && phoIDMVA > %f", mva_EB))
                // .Define("isEEPho",              Form("abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5 && phoIDMVA > %f", mva_EE))

                .Define("isEBPho",              "abs(phoSCEta) < 1.4442")
                .Define("isEEPho",              "abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5")

                .Define("isFSR",                "FSRSelec(mu1, mu2, phoP4)")
                .Define("isHggPho",             "HggPreSelection(mu1_old, mu2_old, phoP4, rhoAll, nPho, phoSCEta, phoPFChIso, phoPFPhoIso, phoTrkIsoHollowConeDR03, phoR9Full5x5, phoEt, phoSigmaIEtaIEtaFull5x5, phoHoverE)")
                .Define("isGoodPho",            "(isEBPho || isEEPho) && isFSR && phoCalibEt > 25. && isHggPho")

                .Filter("Sum(isGoodPho) > 0", "good pho")

                // pick the photon that can make 3 body mass closest to Z mass
                // HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v:
                // https://github.com/cmkuo/ggAnalysis/blob/106X/ggNtuplizer/plugins/ggNtuplizer_trigger.cc#L261-L263
                .Define("phoIdx",               "GetZPho(mu1, mu2, phoP4, isGoodPho)")
                .Define("pho",                  "phoP4[phoIdx]")
                // .Define("passDiPhoTrg",         "((phoFiredDoubleTrgs[phoIdx] >> 1) & 1) == 1 && ((phoFiredDoubleTrgs[phoIdx] >> 2) & 1) == 1")
                // .Define("passDiPhoTrg",         "((HLTPho >> 14) & 1) == 1")

                // check if fsr photon can match to a conversion vertex
                .Define("convIdx",              "ConvMatch(phoSCEta[phoIdx], phoSCPhi[phoIdx], phoSCE[phoIdx], nConv, convNTrks, convVtxX, convVtxY, convVtxZ, convFitPairPX, convFitPairPY, convFitPairPZ, convFitProb)")
                .Filter("convIdx != -1", "match conv");
    return nf;
}


ROOT::RDF::RNode FindFSRMu(ROOT::RDF::RNode df, std::string era, std::string ext) {
    float leadPt = 20, subleadPt = 10;
    if (ext == "SingleMuTrig"){
        if (era.find("2017") != std::string::npos || era.find("2018") != std::string::npos)
            leadPt = 30;
        else
            leadPt = 27;
    }
    auto nf = df.Define("fsrmu",                "if (DeltaR(mu1.Eta(), pho.Eta(), mu1.Phi(), pho.Phi()) > DeltaR(mu2.Eta(), pho.Eta(), mu2.Phi(), pho.Phi())) return mu2; else return mu1")
                .Define("nfsrmu",               "if (DeltaR(mu1.Eta(), pho.Eta(), mu1.Phi(), pho.Phi()) > DeltaR(mu2.Eta(), pho.Eta(), mu2.Phi(), pho.Phi())) return mu1; else return mu2")

                .Filter(Form("nfsrmu.Pt() > %f && fsrmu.Pt() > %f", leadPt, subleadPt), "HLT pT cut");
    return nf;
}


ROOT::RDF::RNode FindMatchEle(ROOT::RDF::RNode df) {
    std::map<std::string, std::string> MergedID_path;
    MergedID_path["M2EB"] = "/data4/chenghan/external/MergedID/Output_Merged2GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EB/XGB/XGB_modelXGB.txt";
    MergedID_path["M2EE"] = "/data4/chenghan/external/MergedID/Output_Merged2GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EE/XGB/XGB_modelXGB.txt";
    MergedID_path["M1EB"] = "/data4/chenghan/external/MergedID/Output_Merged1GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EB/XGB/XGB_modelXGB.txt";
    MergedID_path["M1EE"] = "/data4/chenghan/external/MergedID/Output_Merged1GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EE/XGB/XGB_modelXGB.txt";

    auto ff = MergedIDPred(df, "eleClass", MergedID_path);

    auto nf = ff.Define("eleIdx",               "GetMatchEle(phoSCEta[phoIdx], eleSCEta)")

                .Filter("eleIdx != -1", "match ele")

                .Define("ele",                  "ROOT::Math::PtEtaPhiMVector v(eleCalibPt[eleIdx], eleEta[eleIdx], elePhi[eleIdx], M_ELE); return v;");
    return nf;
}


ROOT::RDF::RNode DefineFinalVars(ROOT::RDF::RNode df){
    auto nf = df.Define("Z",                    "mu1 + mu2 + pho")
                .Define("dimu",                 "mu1 + mu2")
                .Define("zMass",                "Z.M()")
                .Define("phoEleVeto_Lead",      "phoEleVeto[phoIdx]")
                .Define("phoSCEta_Lead",        "phoSCEta[phoIdx]")
                .Define("phoSCPhi_Lead",        "phoSCPhi[phoIdx]")
                .Define("phoSCE_Lead",          "phoSCE[phoIdx]")
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
                .Define("convFitPairP_Lead",    "sqrt(pow(convFitPairPX[convIdx], 2) + pow(convFitPairPY[convIdx], 2) + pow(convFitPairPZ[convIdx], 2))")
                .Define("convFitPairPt_Lead",   "sqrt(pow(convFitPairPX[convIdx], 2) + pow(convFitPairPY[convIdx], 2))")
                .Define("convD0_Lead",          "convD0[convIdx]")
                .Define("convDz_Lead",          "convDz[convIdx]")
                .Define("convL0_Lead",          "convL0[convIdx]")
                .Define("convLz_Lead",          "convLz[convIdx]")

                .Define("nGsfMatchToReco_Lead", "nGsfMatchToReco[eleIdx]")
                .Define("eleConvVeto_Lead",     "eleConvVeto[eleIdx]")
                .Define("eleClass_Lead",        "eleClass[eleIdx]")
                .Define("eleXGBID_Lead",        "eleXGBID[eleIdx]")
                .Define("eleSCEta_Lead",        "eleSCEta[eleIdx]")
                .Define("eleSCPhi_Lead",        "eleSCPhi[eleIdx]")
                .Define("eleSCEn_Lead",         "eleSCEn[eleIdx]")
                .Define("eleCalibEn_Lead",      "eleCalibEn[eleIdx]")
                .Define("eleTrkMissHits_Lead",  "eleTrkMissHits[eleIdx]")
                .Define("eleSubTrkMissHits_Lead", "eleSubTrkMissHits[eleIdx]")

                .Define("eleCalibPt_Lead",      "eleCalibPt[eleIdx]")
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
                .Define("eleTrk_Lead",          "eleTrk1[eleIdx]")
                .Define("eleTrk_subLead",       "eleTrk2[eleIdx]");
    return nf;
}


void xAna(const std::string infile, const std::string outfile, const std::string ext, const std::string era, const bool isMC, const float xs, const float lumi){
    TStopwatch time;
    time.Start();

    std::cout << "[INFO] Read_File(): " << infile.c_str() << std::endl;
    ROOT::EnableImplicitMT(20);

    //=======================================================//
    // Load TTree into RDataFrame and then do the selections //
    //=======================================================//
    ROOT::RDF::RNode df = ROOT::RDataFrame("ggNtuplizer/EventTree", infile.c_str());

    auto df1 = df.Define("isMC", [&, isMC]{return isMC;});
    auto df2 = AddWeights(df1, isMC, era, xs ,lumi);
    auto df3 = FilterMinimum(df2, ext);
    auto df4 = FindGSFTracks(df3);
    auto df5 = FindGoodMus(df4, era);
    auto df6 = FindGoodPho(df5, era);
    auto df7 = FindFSRMu(df6, era, ext);
    auto df8 = FindMatchEle(df7);
    auto df9 = DefineFinalVars(df8);
    auto dfFin = df9;

    std::vector<std::string> defColNames = dfFin.GetDefinedColumnNames();
    std::vector<std::string> Vars = {"event", "zMass", "isDoubleMuTrig", "passDiPhoTrg"};
    for (int i = 0; i < defColNames.size(); i++){
        size_t foundLead = defColNames[i].find("_Lead");
        bool foundP4 = dfFin.GetColumnType(defColNames[i]) == "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >";

        if (foundLead != string::npos || foundP4)
            Vars.push_back(defColNames[i]);
    }
    if (isMC){
        // concatenate Vars and weiVars
        vector<string> weiVars = {
            "puwei", "puwei_up", "puwei_down",
            "genwei", "mcwei", "poswei", "wei",
        };
        Vars.insert(Vars.end(), weiVars.begin(), weiVars.end());
    }
    std::cout << "[INFO] Save_File(): " << outfile.c_str() << std::endl;
    auto d = dfFin.Report();
    dfFin.Snapshot("miniTree", outfile, Vars);

    std::cout << "[INFO] Cut flow:" << std::endl;
    d->Print();

    std::cout << "[INFO] Time taken: " << std::endl;
    time.Stop();
    time.Print();

    std::cout << std::endl;
}