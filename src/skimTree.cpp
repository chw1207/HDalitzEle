#include <iostream>
#include <fstream>
#include <cmath>
#include <mutex>
#include <string>
#include <vector>
#include <limits> // numeric_limits
#include <algorithm>
#include <filesystem>
#include "yaml-cpp/yaml.h"
#include "boost/program_options.hpp"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "Math/Vector4D.h"
#include "TStopwatch.h"
#include "TString.h"
#include "MergedEnCorrector.h"
#include "Utilities.h"
#include "GsfTracks.h"
#define CGREEN "\033[0;32m"
#define CEND "\033[0m"
#define cprintf(X) printf("%s%s%s", CGREEN, X, CEND)

namespace fs = std::filesystem;
namespace po = boost::program_options;
template <typename T>
using Vec = const ROOT::RVec<T>&;

//* Script to skim the ntuples:  ../bin/skimTree --config ../config/UL2017_SignalMC.yaml
//  For signal:  
//  1. remove events with diElectron mass > 60 GeV at LHE level.
//  2. split the events into two parts. one for ananlysis (skimTree) the other one for ID and regression training (trainTree)
//  3. for skimTree, several preselections are applied
//      a. HLTs (double photons or double electrons)
//      b. good primary vertex
//      c. at least one electron and one photn in the event
//  4. for skimTree, energy regression branch "eleHDALRegPt" and gsf track related branches are added (check src/MergedEnCorrector.cpp)
//  5. To-Do: add R9 correction and Residual correction 
//  For data:  
//  1. apply the same preselections as signal
//  2. branch "eleHDALRegPt" and gsf track related branches are also added


float CalcMLL(Vec<int> lhePID, Vec<float> lhePx, Vec<float> lhePy, Vec<float> lhePz){
    // MLL cut at 60 GeV
    auto LHE_Eles = abs(lhePID) == 11; // gen electron
    float mass = std::numeric_limits<float>::max();
    if (ROOT::VecOps::Nonzero(LHE_Eles).size() == 2){
        auto LHEElePx = lhePx[LHE_Eles];
        auto LHEElePy = lhePy[LHE_Eles];
        auto LHEElePz = lhePz[LHE_Eles];
        ROOT::Math::PxPyPzMVector ele1(LHEElePx[0], LHEElePy[0], LHEElePz[0], 0.000511);
        ROOT::Math::PxPyPzMVector ele2(LHEElePx[1], LHEElePy[1], LHEElePz[1], 0.000511);
        mass = (ele1 + ele2).M();
    }
    return mass;
}


bool IsPhoIntConv(const int nMC, Vec<int> mcPID, Vec<int> mcMomPID, Vec<unsigned short> mcStatusFlag){
    // remove internal conversion photon
    int realpho = 0;
    for (int i = 0; i < nMC; i++){
        bool hardProc = (mcStatusFlag[i] >> 0 & 1) == 1;
        bool isPrompt = (mcStatusFlag[i] >> 1 & 1) == 1;
        bool isRealPho = mcPID[i] == 22 && mcMomPID[i] == 25 && hardProc && isPrompt;
        if (isRealPho)
            realpho += 1;
    }
    bool isPhoIntConv = (realpho == 0) ? true : false; // no real photon
    return isPhoIntConv;
}


ROOT::RDF::RNode FindGSFTracks(ROOT::RDF::RNode df) {
    auto make_ditrkPt = [](Vec<ROOT::Math::PtEtaPhiMVector> eleTrk1, Vec<ROOT::Math::PtEtaPhiMVector> eleTrk2){
        ROOT::RVec<float> ditrkPt(eleTrk1.size());
        for (size_t i = 0; i < eleTrk1.size(); i++){
            ditrkPt[i] = (float)(eleTrk1[i]+eleTrk2[i]).Pt();
        }
        return ditrkPt;
    };

    // match the gsf tracks to the associated electron
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
                .Define("diTrkPt",              make_ditrkPt,               {"eleTrk1", "eleTrk2"})
                .Define("gsfPtRatio",           gsf::GetTrkPtRatio,         {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfDeltaR",            gsf::GetTrkdR,              {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfPtSum",             gsf::GetTrkPtSum,           {"nEle", "nGsfMatchToReco", "eleTrk1", "eleTrk2"})
                .Define("gsfRelPtRatio",        gsf::GetTrkRelPtRatio,      {"nEle", "eleSCRawEn", "nGsfMatchToReco", "eleTrk1", "eleTrk2"});
    return nf;
}


int main(int argc, char** argv){
    TStopwatch time;
	time.Start();

    std::string config_path;
    int range = -1;
    po::options_description desc{"Options"};
    desc.add_options()
        ("help,h",                                  "Higgs Dalitz decay skim script")
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
    auto ntuple_path = cfg["ntuple_path"].as<std::vector<std::string>>();
    auto skimTree_path = cfg["skimTree_path"].as<std::vector<std::string>>();
    if (ntuple_path.size() != skimTree_path.size())
        throw std::runtime_error("Number of ntuple_path != number of skimTree_path");

    // check the directory to save the miniTree exist or not. if not, create it
    std::string direc = std::filesystem::path(skimTree_path[0]).parent_path();
    if (!std::filesystem::exists(direc)){
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
    //          Run skimming over all of the ntuples         //
    //=======================================================//
   if (range == -1 && nthreads != 1)
        ROOT::EnableImplicitMT(nthreads); // enable multithreading
    cprintf("Start to skim the ggNtuples!\n");
    cprintf("Skim the following ntuples sequentially\n");
    for (size_t i = 0; i < ntuple_path.size(); i++){
        cprintf(Form(" %ld. %s\n", i+1, ntuple_path[i].c_str()));
    }
    for (size_t i = 0; i < ntuple_path.size(); i++){
        printf("********************************   %ld    ********************************\n", i+1);
        printf(" [+] Read ggNtuple from:\n");
        printf("     - %s\n", ntuple_path[i].c_str());
        printf(" [+] Save skimTree/trainTree in:\n");
        printf("     - %s\n", skimTree_path[i].c_str());
        
        // load the ntuples to dataframe
        ROOT::RDF::RNode df = ROOT::RDataFrame("ggNtuplizer/EventTree", ntuple_path[i]);
        if (range != -1)
            df = df.Range(0, range);
        
        // HLT trigger
        //* Di-photon trigger (for merged)
        //*     bit position: https://github.com/cmkuo/ggAnalysis/blob/106X/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc#L364
        //*     1) HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v (2016)
        //*     2) HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v (2017 and 2018)
        //* Di-electron trigger (for resolved)
        //*     bit position: https://github.com/cmkuo/ggAnalysis/blob/106X/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc#L297
        //*     bit position: https://github.com/cmkuo/ggAnalysis/blob/106X/ggNtuplizer/plugins/ggNtuplizer_globalEvent.cc#L332
        //*     1) HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v (2016)
        //*     2) HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v (2017 and 2018)   
        std::string diEleTrigger = "((HLTEleMuX >> 40 & 1) == 1) || ((HLTEleMuX >> 5 & 1) == 1)";
        if (era.find("2016") != std::string::npos)
            diEleTrigger = "(HLTEleMuX >> 40 & 1) == 1"; 

        // energy regression models
        auto model_path = cfg["external_files"]["energy_correction"]["energy_reg"].as<std::vector<std::string>>();

        if (isMC){
            auto df1 = df.Define("diLepMCMass",     CalcMLL,      {"lhePID", "lhePx", "lhePy", "lhePz"})
                         .Define("HasIntPho",       IsPhoIntConv, {"nMC", "mcPID", "mcMomPID", "mcStatusFlag"})
                         .Filter("diLepMCMass < 60 && HasIntPho == 0", "mc preslections");

            // calculate the weight wrt cross section
            auto df2 = df1.Filter("event % 2 == 1", "analysis"); // == 1 for analysis, == 0 for ID training and regression
            printf(" [+] Weight application for analysis tree:\n");
            auto pos = df2.Filter("genWeight > 0",   "positive event").Count();
            auto neg = df2.Filter("genWeight <= 0",  "negative event").Count();
            const int totalev = pos.GetValue() - neg.GetValue();
            const float xs = cfg["cross_section"].as<std::vector<float>>().at(i);
            const float lumi = cfg["luminosity"].as<float>();
            const float mcwei = ((xs * lumi) / totalev);
            printf("     - XS = %f fb, Luminosity = %f fb^-1, Nevents = %d \n", xs, lumi, totalev);
            printf("     - MC weight (XS * Luminosity / Nevents): %f \n", mcwei);
            
            // save the analysis tree 
            auto df_skim  = df1.Define("mcwei",   [&, mcwei](){return mcwei;})
                               .Define("isDiPhoHLT",           "((HLTPho >> 14) & 1) == 1")
                               .Define("isDiEleHLT",           diEleTrigger)
                               .Define("eleESEnToRawE",        "(eleESEnP1+eleESEnP2)/eleSCRawEn")
                               .Filter("isDiPhoHLT || isDiEleHLT", "diPho || diEle HLT")
                               .Filter("isPVGood == 1", "Good Vtx")
                               .Filter("nEle > 0 && nPho > 0 && nGSFTrk > 0", "none zero particles");           
            auto df_skim_1 = FindGSFTracks(df_skim);   
            auto df_skim_2 = doHDALRegressionXGB(df_skim_1, model_path);               

            // trigger the event loop and save the skimtree
            std::vector<std::string> out_variables;
            std::vector<std::string> variables = df_skim_2.GetColumnNames();
            for (size_t j = 0; j < variables.size(); j++){
                if (variables[j] == "ambGSF" || variables[j] == "isMainGSF" || variables[j] == "eleTrk1" || variables[j] == "eleTrk2")
                    continue;
                out_variables.push_back(variables[j]);
            }
            auto report = df_skim_2.Report(); 
            df_skim_2.Snapshot("skimTree", skimTree_path[i], out_variables);
            TString cutflow = utils::printReport(report);

            // save the training tree 
            ROOT::RDF::RSnapshotOptions opts;
            opts.fMode = "update";
            auto colNames = df1.GetColumnNames();
            auto df_train = df1.Filter("event % 2 == 0", "training"); // == 1 for analysis, == 0 for ID training and regression
            df_train.Snapshot("trainTree", skimTree_path[i], colNames, opts);
            printf("\n");

        }
        else{ // data
            auto df1  = df.Define("isDiPhoHLT",           "((HLTPho >> 14) & 1) == 1")
                          .Define("isDiEleHLT",           diEleTrigger)
                          .Define("eleESEnToRawE",        "(eleESEnP1+eleESEnP2)/eleSCRawEn")
                          .Filter("isDiPhoHLT || isDiEleHLT", "diPho || diEle HLT")
                          .Filter("isPVGood == 1", "Good Vtx")
                          .Filter("nEle > 0 && nPho > 0 && (nGSFTrk > 0)", "none zero particles");
            auto df2 = FindGSFTracks(df1); 
            auto df3 = doHDALRegressionXGB(df2, model_path); 
            auto report = df3.Report(); 

            // trigger the event loop and save the skimtree
            std::vector<std::string> out_variables;
            std::vector<std::string> variables = df3.GetColumnNames();
            for (size_t j = 0; j < variables.size(); j++){
                if (variables[j] == "ambGSF" || variables[j] == "isMainGSF" || variables[j] == "eleTrk1" || variables[j] == "eleTrk2")
                    continue;
                out_variables.push_back(variables[j]);
            }
            df3.Snapshot("skimTree", skimTree_path[i], out_variables);
            TString cutflow = utils::printReport(report);
        }
    }
    
    time.Stop();
    time.Print();

    return 0;
}