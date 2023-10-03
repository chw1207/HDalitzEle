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
#include "Utilities.h"
#include "PUWeightCalculator.h"
#include "./interface/help.h"

template <typename T>
using Vec = const ROOT::RVec<T>&;
using namespace ROOT::VecOps;


// SM Higgs XS: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV
std::map<std::string, float> xsbr = { // xs times branching ratio (fb)
    {"ggF_120GeV",       (52.22  * 1000) * 3.80E-5},
    {"VBF_120GeV",       (3.935  * 1000) * 3.80E-5},
    {"WH_120GeV",        (1.565  * 1000) * 3.80E-5},
    {"ZH_120GeV",        (0.9939 * 1000) * 3.80E-5},
    {"ttH_120GeV",       (0.5697 * 1000) * 3.80E-5},
    {"bbH_120GeV",       (0.5534 * 1000) * 3.80E-5},

    {"ggF_125GeV",       (48.58  * 1000) * 3.92E-5},
    {"VBF_125GeV",       (3.782  * 1000) * 3.92E-5},
    {"WH_125GeV",        (1.373  * 1000) * 3.92E-5},
    {"ZH_125GeV",        (0.8839 * 1000) * 3.92E-5},
    {"ttH_125GeV",       (0.5071 * 1000) * 3.92E-5},
    {"bbH_125GeV",       (0.4880 * 1000) * 3.92E-5},

    {"ggF_130GeV",       (45.31  * 1000) * 3.90E-5},
    {"VBF_130GeV",       (3.637  * 1000) * 3.90E-5},
    {"WH_130GeV",        (1.209  * 1000) * 3.90E-5},
    {"ZH_130GeV",        (0.7899 * 1000) * 3.90E-5},
    {"ttH_130GeV",       (0.4539 * 1000) * 3.90E-5},
    {"bbH_130GeV",       (0.4304 * 1000) * 3.90E-5}
};


std::map<std::string, float> lumis = {
    {"UL2016preVFP",     19.52},
    {"UL2016postVFP",    16.81},
    {"UL2017",           41.48},
    {"UL2018",           59.82}
};


ROOT::RDF::RNode AddWeights(ROOT::RDF::RNode df, const std::string era, const std::string mode, int HiggsMass){
    // calculate the mc weight
    auto all = df.Count();
    auto pos = df.Filter("genWeight > 0",   "positive event").Count();
    auto neg = df.Filter("genWeight <= 0",  "negative event").Count();
    const int totalev = pos.GetValue() - neg.GetValue();
    const float procXS = xsbr.at(fmt::format("{}_{}GeV", mode, HiggsMass));
    const float instwei = procXS/xsbr.at(fmt::format("ggF_{}GeV", HiggsMass));
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


ROOT::RDF::RNode FindGenParticle(ROOT::RDF::RNode df){
    auto nf = df.Define("hardProc",         "Helper::GenType(mcStatusFlag, 0)")
                .Define("isPrompt",         "Helper::GenType(mcStatusFlag, 1)")

                .Define("isDALMu",          "abs(mcPID) == 13 && mcMomPID == 25") // higgs dalitz muon
                .Define("isHZGMu",          "abs(mcPID) == 13 && mcMomPID == 23 && mcGMomPID == 25") // higgs to zg muon
                .Define("isGenMu",          "(isDALMu || isHZGMu) && hardProc && isPrompt && mcEta < 2.5")

                .Define("isGenPho",         "abs(mcPID) == 22 && mcMomPID == 25 && hardProc && isPrompt && mcEta < 2.5")

                // kick out the events without 2 gen muons and 1 gen photon
                .Filter("Sum(isGenMu) == 2", "2 gen mu")
                .Filter("Sum(isGenPho) == 1", "1 gen pho") // may loose some events due to the photon internal conversion

                // gen particle index
                .Define("genMuIdx1",        [](Vec<int> good, Vec<float> pt){return utils::getIdx(good, pt)[0];}, {"isGenMu", "mcPt"})
                .Define("genMuIdx2",        [](Vec<int> good, Vec<float> pt){return utils::getIdx(good, pt)[1];}, {"isGenMu", "mcPt"})
                .Define("genPhoIdx1",       [](Vec<int> good, Vec<float> pt){return utils::getIdx(good, pt)[0];}, {"isGenPho", "mcPt"})

                .Define("M_Mu",             "(float) 0.1057")
                .Define("GenMu_Lead",       "ROOT::Math::PtEtaPhiMVector v(mcPt[genMuIdx1], mcEta[genMuIdx1], mcPhi[genMuIdx1], mcMass[genMuIdx1]); return v;")
                .Define("GenMu_subLead",    "ROOT::Math::PtEtaPhiMVector v(mcPt[genMuIdx2], mcEta[genMuIdx2], mcPhi[genMuIdx2], mcMass[genMuIdx2]); return v;")
                .Define("GenPho_Lead",      "ROOT::Math::PtEtaPhiMVector v(mcPt[genPhoIdx1], mcEta[genPhoIdx1], mcPhi[genPhoIdx1], 0); return v;")
                .Define("diGenMu",          "GenMu_Lead + GenMu_subLead")

                .Define("GenMuPtRatio",     "GenMu_subLead.Pt() / GenMu_Lead.Pt()")
                .Define("GenMuDeltaR",      "DeltaR(GenMu_Lead.Eta(), GenMu_subLead.Eta(), GenMu_Lead.Phi(), GenMu_subLead.Phi())")
                .Define("GenMuPtSum",       "GenMu_Lead.Pt() + GenMu_subLead.Pt()");
    return nf;
}


ROOT::RDF::RNode GenToReco(ROOT::RDF::RNode df){
    auto nf = df.Define("recoMuIdx1",                      "Helper::RecoIdx(GenMu_Lead, muPt, muEta, muPhi, M_Mu)")
                .Define("recoMuIdx2",                      "Helper::RecoIdx(GenMu_subLead, muPt, muEta, muPhi, M_Mu)")
                .Define("recoPhoIdx1",                     "Helper::RecoIdx(GenPho_Lead, phoCalibEt, phoEta, phoPhi, 0.)")

                // only event with proper reco particles left
                .Filter("recoMuIdx1 != -1 && recoMuIdx2 != -1", "2 reco mu")
                .Filter("recoPhoIdx1 != -1", "1 reco pho")

                // P4 of reco particle
                .Define("RecoMu_Lead",                     "ROOT::Math::PtEtaPhiMVector v(muPt[recoMuIdx1], muEta[recoMuIdx1], muPhi[recoMuIdx1], M_Mu); return v;")
                .Define("RecoMu_subLead",                  "ROOT::Math::PtEtaPhiMVector v(muPt[recoMuIdx2], muEta[recoMuIdx2], muPhi[recoMuIdx2], M_Mu); return v;")
                .Define("RecoPho_Lead",                    "ROOT::Math::PtEtaPhiMVector v(phoCalibEt[recoPhoIdx1], phoEta[recoPhoIdx1], phoPhi[recoPhoIdx1], 0.); return v;");
    return nf;
}


ROOT::RDF::RNode DefineFinalVars(ROOT::RDF::RNode df){
                // leading gen variables
    auto nf = df.Define("mcPt_Lead",                        "mcPt[genMuIdx1]")
                .Define("mcEta_Lead",                       "mcEta[genMuIdx1]")
                .Define("mcPhi_Lead",                       "mcPhi[genMuIdx1]")
                .Define("mcVtx_Lead",                       "mcVtx[genMuIdx1]")
                .Define("mcVty_Lead",                       "mcVty[genMuIdx1]")
                .Define("mcVtz_Lead",                       "mcVtz[genMuIdx1]")

                // trailing gen variables
                .Define("mcPt_subLead",                     "mcPt[genMuIdx2]")
                .Define("mcEta_subLead",                    "mcEta[genMuIdx2]")
                .Define("mcPhi_subLead",                    "mcPhi[genMuIdx2]")
                .Define("mcVtx_subLead",                    "mcVtx[genMuIdx2]")
                .Define("mcVty_subLead",                    "mcVty[genMuIdx2]")
                .Define("mcVtz_subLead",                    "mcVtz[genMuIdx2]")

                .Define("higgsMass",                        "(RecoMu_Lead + RecoMu_subLead + RecoPho_Lead).M()");

    return nf;
}



void xAnaGenMu(std::string infile, std::string outfile, std::string era, std::string mode, int HiggsMass){
    TStopwatch time;
    time.Start();

    fmt::print("[INFO] Read_File(): {}\n", infile);
    fmt::print("[INFO] Save_File(): {}\n", outfile);

    //=======================================================//
    // Load TTree into RDataFrame                            //
    //=======================================================//
    ROOT::EnableImplicitMT(20);
    ROOT::RDF::RNode df = ROOT::RDataFrame("skimTree", infile.c_str());
    fmt::print("[INFO] Process: {}_{} GeV\n", mode, HiggsMass);

    auto df1 = AddWeights(df, era, mode, HiggsMass);
    auto df2 = FindGenParticle(df1);
    auto df3 = GenToReco(df2);
    auto df4 = DefineFinalVars(df3);
    auto dfFin = df4;

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
    auto report = dfFin.Report();
    dfFin.Snapshot("miniTree", outfile.c_str(), Vars);

    //========================================================//
    // Visualize the selection results (cut flow, event count)//
    //========================================================//
    fmt::print("[INFO] Cut flow:\n");
    report->Print();

    time.Stop();
    time.Print();
}