#include <iostream>
#include <filesystem> // c++17
#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include "yaml-cpp/yaml.h"
#include "boost/algorithm/string.hpp"
#include "boost/algorithm/cxx11/any_of.hpp"
#include "boost/program_options.hpp"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RSnapshotOptions.hxx"
#include "Math/Vector4D.h"
#include "TChain.h"
#include "TInterpreter.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TCollection.h"
#include "TMath.h"
#include "TObjString.h"

#include "WeightHandler.h"
#include "PUWeightCalculator.h"
#include "Utilities.h"
#include "GsfTracks.h"
#include "PhotonSel.h"
#include "ElectronSel.h"
#include "XGBReader.h"
#include "Categorizer.h"
#include "Generator.h"
#include "EnCalibrater.h"

#define CGREEN  "\033[0;32m"
#define CRED    "\033[1;31m"
#define CYELLOW "\033[1;33m"
#define CEND    "\033[0m"
#define Info(X)  printf("%s%s%s", CGREEN,  X, CEND)
#define Warn(X)  printf("%s%s%s", CYELLOW, X, CEND)
#define Error(X) printf("%s%s%s", CEND,    X, CEND)

namespace fs = std::filesystem;
namespace po = boost::program_options;
namespace ba = boost::algorithm;

// global variables
int range = -1;
bool saveTrainTree = false;

std::string config_path;
YAML::Node cfg;
std::vector<std::string> eras = {"UL2016preVFP", "UL2016postVFP", "UL2017", "UL2018"};
std::string era;
int nthreads = 1; 
bool isMC;
std::vector<std::string> readpaths;
std::vector<std::string> savepaths;
std::vector<std::string> scorpaths;
bool doVariation = false;
std::vector<std::string> variations; // booked by config (should be one of the key of allowed_variations)
std::map<std::string, std::string> allowed_variations = {
    // shower shape correction
    {"PhoNoR9Corr",     "phoR9Full5x5"},

    // photon energy 
    {"PhoScaleStatUp",  "phoScale_stat_up"},
    {"PhoScaleSystUp",  "phoScale_syst_up"},
    {"PhoScaleGainUp",  "phoScale_gain_up"},
    {"PhoScaleStatDo",  "phoScale_stat_dn"},
    {"PhoScaleSystDo",  "phoScale_syst_dn"},
    {"PhoScaleGainDo",  "phoScale_gain_dn"},
    {"PhoSigmaPhiUp",   "phoResol_phi_up"},
    {"PhoSigmaRhoUp",   "phoResol_rho_up"},
    {"PhoSigmaRhoDo",   "phoResol_rho_dn"},

    // merged electron electron calibration
    {"EleHDALScaleUp",  "+ eleHDALScaleErr"},
    {"EleHDALScaleDo",  "- eleHDALScaleErr"},
    {"EleHDALSmearUp",  "eleHDALCalibPtUp"},
    {"EleHDALSmearDo",  "eleHDALCalibPtDo"},

    // official electron energy calibration
    {"EleScaleStatUp",  "eleScale_stat_up"},
    {"EleScaleSystUp",  "eleScale_syst_up"},
    {"EleScaleGainUp",  "eleScale_syst_up"},
    {"EleScaleStatDo",  "eleScale_stat_dn"},
    {"EleScaleSystDo",  "eleScale_syst_dn"},
    {"EleScaleGainDo",  "eleScale_gain_dn"},
    {"EleSigmaPhiUp",   "eleResol_phi_up"},
    {"EleSigmaRhoUp",   "eleResol_rho_up"},
    {"EleSigmaRhoDo",   "eleResol_rho_dn"},
    
    // jet energy
    {"JERUp",           "jetP4SmearUp"},
    {"JERDo",           "jetP4SmearDo"},
    {"JECUp",           " + jetJECUnc"},
    {"JECDo",           " - jetJECUnc"},
};
std::vector<std::string> categories = {
    "Merged-2Gsf-Loose",
    "Merged-2Gsf-VBF", 
    "Merged-2Gsf-BST", 
    "Merged-2Gsf-EBHR9", 
    "Merged-2Gsf-EBLR9", 
    "Merged-2Gsf-EE"
};
std::vector<std::string> scale_factors = { // required scale factors, if values are not provided by config, they will be assigned as 1.  
    "RecoEleGt20",
    "MissHitsTrk",
    "MissHitsSubTrk",
    "ConvVetoTrk",
    "MergedEleID",
    "HggPreselForEle",
    "HggPreselForPho",
    "HggPhoCSEV",
    "HggPhoID",
    "DiPhoHLTSeedForEle",
    "DiPhoHLTUnSeedForEle",
    "DiPhoHLTSeedForPho",
    "DiPhoHLTUnSeedForPho"
};


void ArgumentParser(int argc, char** argv){
    po::options_description desc{"Options"};
    desc.add_options()
        ("help,h",                                      "Higgs Dalitz decay analysis script to select a merged elctron and a photon")
        ("config,c",        po::value(&config_path),    "configuration file (required)")
        ("range,r",         po::value(&range),          "maximal number of events to process; -1 means all; if specified MT will be turned off")
        ("trainTreeOnly",   po::value(&saveTrainTree),  "produce the tree for ID training and regression (only for MC)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) {
        std::cout << desc << std::endl;
        gSystem->Exit(0);
    }
    if (vm.count("config") < 1)
        throw std::invalid_argument("the option '--config' is required but missing");

    // Load config file and setup the common parameters 
    cfg = YAML::LoadFile(config_path);
    if (cfg["era"]){
        era = cfg["era"].as<std::string>();
        if (!ba::any_of_equal(eras, era))
            throw std::invalid_argument(Form("era: %s is not supported", era.c_str()));
    }
    else{
        throw std::runtime_error("the node 'era' cannot be found in config file");
    }

    if (cfg["threads"])
        nthreads = cfg["threads"].as<int>();

    if (cfg["isMC"])
        isMC = cfg["isMC"].as<bool>();
    else
        throw std::runtime_error("the node 'isMC' cannot be found in config file");

    if (cfg["ntuple_path"] && cfg["miniTree_path"]){
        readpaths = cfg["ntuple_path"].as<std::vector<std::string>>();
        savepaths = cfg["miniTree_path"].as<std::vector<std::string>>();
        if (readpaths.size() != savepaths.size())
            throw std::runtime_error("number of ntuple_path != number of miniTree_path");
    }
    else{
        throw std::runtime_error("the node 'ntuple_path' or 'miniTree_path' cannot be found in config file");
    }

    doVariation = (cfg["variation"]) ? true : false;
    if (doVariation) 
        variations = cfg["variation"].as<std::vector<std::string>>();
}


void DeclareMVAFunc(){
    // helper function to construct ID features vector (load to interpreter)
    gInterpreter->Declare(std::string(
        "template<std::size_t N, typename T, typename... Args_t>"
        "ROOT::RVec<vector<T>> MakeMVAVals(const ROOT::RVec<Args_t>&... args){"
        "    const std::size_t sizes[] = {args.size()...};"
        "    std::vector<std::vector<T>> outtrans(sizes[0], std::vector<T>(N, 0));"
        "    for(int i = 0; i < sizes[0]; i++) { "  
        "        outtrans[i] = {args[i]...};"
        "    }"
        "    return outtrans;"
        "}"
    ).c_str());
}


// assuming that the only difference of the miniTree path for merged_HLT and merged category is "HLT".
// eg. 
//      merged:     /data4/chenghan/electron/miniTree_merged/UL2017/miniTree_HDalitz_ggF_eeg_125_UL2017.root
//      merged_HLT: /data4/chenghan/electron/miniTree_merged_HLT/UL2017/miniTree_HDalitz_ggF_eeg_125_UL2017.root
std::vector<Long64_t> GetMergedEvents(std::string repath){
    TString MergedPath(repath);
    MergedPath.ReplaceAll("merged_HLT", "merged");
    
    std::vector<Long64_t> v; 
    TString MergedDir = gSystem->GetDirName(MergedPath.Data());
    if (!fs::exists(MergedDir.Data())){
        Warn(Form("Fail to find the root files for merged category: %s\n", MergedDir.Data()));
        return v;
    }

    auto df_merged = ROOT::RDataFrame("miniTree", MergedPath.Data(), {"event"});
    auto evCol = df_merged.Take<Long64_t>("event");

    printf(" [+] Find out the merged events in :\n");
    printf("     - %s\n", MergedPath.Data());
    printf("     - %ld events are found\n", (*evCol).size());
    return *evCol;
}


// main analysis for Heeg, where ee merged.
bool xAna(std::string inpath, std::string outpath, int iset){
    TStopwatch time_iset;
    time_iset.Start();

    //=======================================================//
    //              Load the files into RDataframe           //
    //=======================================================//
    printf(" [+] Read ggNtuple from:\n");
    printf("     - %s\n", inpath.c_str());
    TChain chain1("ggNtuplizer/EventTree");
    TChain chain2("CorrTree");
    bool isSignalMC = ba::contains(inpath, "Dalitz");

    TString InPathStr(inpath);
    TString InDir = gSystem->GetDirName(InPathStr.Data());
    TString SSInDir = gSystem->GetDirName(InPathStr.Data());
    bool sspath_found = false;
    if (isSignalMC){ // find the path containing shower shape corrections
        // eg: ntuple : "/data5/ggNtuples/V10_06_30_02/job_UL17_Dalitz_eeg_VBF_m130"
        //     sscorr : "/data5/ggNtuples/V10_06_30_02_sscorr/job_UL17_Dalitz_eeg_VBF_m130"
        TString tok;
        int position = 0;
        while (InDir.Tokenize(tok, position, "/")){
            if (tok.Contains("V10")){
                SSInDir.ReplaceAll(tok, tok+"_sscorr");
                break;
            }
        }
        if (fs::exists(SSInDir.Data())){
            printf(" [+] Shower shape corrected path found:\n");
            printf("     - %s\n", SSInDir.Data());
            sspath_found = true;
        }
        else 
            throw std::runtime_error("The path which contains shower shape corrected files cannot be found");
    }

    if (InPathStr.EndsWith("*.root")){ // all the files
        TSystemDirectory SD("SD", InDir.Data());
        TList* FileList = SD.GetListOfFiles();
        if (FileList){
            TIter Next(FileList);
            while (auto File = (TSystemFile*) Next()){
                TString fName = File->GetName();
                if (File->IsDirectory() || !fName.EndsWith(".root"))
                    continue;
                chain1.Add(Form("%s/%s", InDir.Data(), fName.Data()));
                if (sspath_found){
                    fName.ReplaceAll(".root", "_corr.root");
                    const char* fullSSPath = Form("%s/%s", SSInDir.Data(), fName.Data());
                    if (fs::exists(fullSSPath))
                        chain2.Add(fullSSPath);
                    else
                        throw std::runtime_error(Form("File doesn't exist: %s", fullSSPath));
                }
            }
            delete FileList;
        }
        if (sspath_found)
            chain1.AddFriend(&chain2); 
    }
    else{ // one file
        chain1.Add(InPathStr.Data());
        if (sspath_found){
            TString SSPathStr(InPathStr);
            SSPathStr.ReplaceAll(InDir.Data(), SSInDir.Data());
            SSPathStr.ReplaceAll(".root", "_corr.root");
            if (fs::exists(SSPathStr.Data()))
                chain2.Add(SSPathStr.Data());
            else
                throw std::runtime_error(Form("File doesn't exist: %s", SSPathStr.Data()));
            chain1.AddFriend(&chain2); 
        }
    }

    auto df = ROOT::RDataFrame(chain1);
    ROOT::RDF::RNode df_analysis(df);
    if (range != -1)
        df_analysis = df_analysis.Range(0, range);

    //=======================================================//
    //                 weight calculation for MC             //
    //=======================================================//
    PUWeightCalculator puCalc[3];
    float mcwei = 1.;
    if (isMC){
        auto pufiles = cfg["external_files"]["pileup_rewei"].as<std::map<std::string, std::string>>();

        // remove events with dielectron mass > 60 and photon is internally converted 
        auto nf =  (isSignalMC) ?  ROOT::RDF::AsRNode(df.Define("diLepMCMass",     gen::CalcLHEMee,   {"lhePID", "lhePx", "lhePy", "lhePz"})
                                                        .Define("HasIntPho",       gen::IsPhoIntConv, {"mcPID", "mcMomPID", "mcStatusFlag"})
                                                        .Filter("diLepMCMass < 60 && HasIntPho == 0", "mc preslections"))
                                :  ROOT::RDF::AsRNode(df);
                                       

        if (saveTrainTree){
            // save the tree for training
            auto skimTree_path = cfg["skimTree_path"].as<std::vector<std::string>>();
            std::string direc = fs::path(skimTree_path[iset]).parent_path();
            if (!fs::exists(direc)){
                int status = gSystem->mkdir(direc.c_str(), true);
                if (status == -1)
                    throw::std::runtime_error(Form("Fail to create directory: %s", direc.c_str()));
            }
            printf(" [+] Save trainTree in:\n");
            printf("     - %s\n", skimTree_path[iset].c_str());
            auto df_train = nf.Filter("event % 2 == 0", "training"); 
            df_train.Snapshot("trainTree", skimTree_path[iset]);
            return 0;
        }
        
        // calculate the mc weight wrt luminosity and puwei
        puCalc[0].Init(pufiles["nominal"].c_str());
        puCalc[1].Init(pufiles["up"].c_str());
        puCalc[2].Init(pufiles["down"].c_str());
        auto get_pu = [&](const int run, const ROOT::RVec<float>& puTrue){
            float puwei = puCalc[0].GetWeight(run, puTrue[0]);
            return puwei;
        };
        auto get_pu_up = [&](const int run, const ROOT::RVec<float>& puTrue){
            float puwei_up = puCalc[1].GetWeight(run, puTrue[0]);
            return puwei_up;
        };
        auto get_pu_do = [&](const int run, const ROOT::RVec<float>& puTrue){
            float puwei_down = puCalc[2].GetWeight(run, puTrue[0]);
            return puwei_down;
        };
        auto af = (isSignalMC) ? nf.Filter("event % 2 == 1",  "analysis") : nf; // == 1 for analysis, == 0 for ID training and regression
        const float xs   = cfg["cross_section"].as<std::vector<float>>()[iset]; 
        const float lumi = cfg["luminosity"].as<float>();
        auto pos = af.Filter("genWeight > 0",   "positive event").Count();
        auto neg = af.Filter("genWeight <= 0",  "negative event").Count();

        // define the weights in the dataframe
        printf(" [+] Weight calculation for analysis tree:\n");
        const int totalev = pos.GetValue() - neg.GetValue();
        mcwei = ((xs * lumi) / totalev);
        printf("     - XS = %f fb, Luminosity = %f fb^-1, Nevents = %d \n", xs, lumi, totalev);
        printf("     - MC weight (XS * Luminosity / Nevents): %f \n", mcwei);
        df_analysis = af.Define("mcwei",            [&](){return mcwei;})
                        .Define("genwei",           "if (genWeight > 0) return 1.; else return -1.;")
                        .Define("puwei",            get_pu,             {"run", "puTrue"})
                        .Define("puwei_up",         get_pu_up,          {"run", "puTrue"})
                        .Define("puwei_down",       get_pu_do,          {"run", "puTrue"});
    }
    else
        df_analysis = df_analysis.Define("weight", "(float) 1.");

    // =======================================================//
    //  minimum event filter (trigger, non zero gamma and e)  //
    // =======================================================//
    //! Note: do not apply DoublePhoHLT on MC (because L1 seed is not well simulated)
    auto merged_ev = GetMergedEvents(outpath);
    std::sort(merged_ev.begin(), merged_ev.end());
    df_analysis = df_analysis.Filter([&](Long64_t ev){
        if (std::binary_search(merged_ev.begin(), merged_ev.end(), ev))
            return false;
        else 
            return true;
    }, {"event"}, "remove merged events");

    // std::string DiEle27HLT = (isMC) ? "true" : "((HLTPho >> 13) & 1) == 1";
    auto mf = df_analysis.Define("era",                 [&](){return era;})
                         .Filter("((HLTPho >> 13) & 1) == 1",  "trigger cut") // HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v
                         .Filter("isPVGood == 1", "Good Vtx")
                         .Filter("nEle > 0 && nPho > 0", "1 ele, 1 pho")
                         .Define("phoS4Full5x5",        "phoE2x2Full5x5/phoE5x5Full5x5")
                         .Define("phoESEnToRawE",       "(phoESEnP1+phoESEnP2)/phoSCRawE"); 

    //=======================================================//
    //                  photon selections                    //
    //=======================================================//
    XGBReader hggReader[2]; // Hgg ID readers (only need for data) 
    std::vector<std::string> hggvals_EB;
    std::vector<std::string> hggvals_EE;
    if (!isSignalMC){
        auto hggModelFiles = cfg["external_files"]["HggPhoID_model"].as<std::map<std::string, std::string>>();
        hggReader[0].Init(hggModelFiles["EB"]);
        hggvals_EB = {
            "phoSCRawE",
            "phoR9Full5x5",
            "phoSigmaIEtaIEtaFull5x5",
            "phoSCEtaWidth",
            "phoSCPhiWidth",
            "phoSigmaIEtaIPhiFull5x5",
            "phoS4Full5x5",
            "phoPFPhoIso",
            "phoPFChIso",
            "phoPFChWorstIso",
            "phoSCEta",
            "phoRho"
        };
        std::string hggvals_EB_str = ba::join(hggvals_EB, ", ");

        hggReader[1].Init(hggModelFiles["EE"]);
        hggvals_EE = {
            "phoSCRawE",
            "phoR9Full5x5",
            "phoSigmaIEtaIEtaFull5x5",
            "phoSCEtaWidth",
            "phoSCPhiWidth",
            "phoSigmaIEtaIPhiFull5x5",
            "phoS4Full5x5",
            "phoPFPhoIso",
            "phoPFChIso",
            "phoPFChWorstIso",
            "phoSCEta",
            "phoRho",
            "phoESEffSigmaRR",
            "phoESEnToRawE"
        };
        std::string hggvals_EE_str = ba::join(hggvals_EE, ", ");
        mf = mf.Define("phoRho",            "ROOT::RVec<float> v(nPho, rho); return v;")
               .Define("phoCorrR9Full5x5",  "phoR9Full5x5")     
               .Define("phoMVAEBVals",      Form("MakeMVAVals<%d, float>(%s)", (int)hggvals_EB.size(), hggvals_EB_str.c_str()))
               .Define("phoMVAEEVals",      Form("MakeMVAVals<%d, float>(%s)", (int)hggvals_EE.size(), hggvals_EE_str.c_str()))
               .Define("phoCorrHggIDMVA",   [&](const ROOT::RVec<float>& phoSCEta,
                                                const ROOT::RVec<std::vector<float>>& phoMVAEBVals,
                                                const ROOT::RVec<std::vector<float>>& phoMVAEEVals){
                                                    ROOT::RVec<float> mva(phoSCEta.size());
                                                    for (size_t i = 0; i < phoSCEta.size(); i++){
                                                        float raw_mva;
                                                        if (fabs(phoSCEta[i]) < 1.479)
                                                            raw_mva = hggReader[0].Compute(phoMVAEBVals[i])[0];
                                                        else
                                                            raw_mva = hggReader[1].Compute(phoMVAEEVals[i])[0];

                                                        // convert to the MVA score between -1 and 1
                                                        mva[i] = 2.0 / (1.0 + TMath::Exp(-2.0 * raw_mva)) - 1;
                                                    }
                                                    return mva;
                                                }, {"phoSCEta", "phoMVAEBVals", "phoMVAEEVals"}); 
    }
    auto pf = mf.Define("isEBPho",                  "abs(phoSCEta) < 1.4442")
                .Define("isEEPho",                  "abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5")  
                .Define("isHggPho",                 phoSel::HggPresel,  {"nPho", "rhoAll", "phoSCEta", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"})
                // .Define("isGoodPho",                "(isEBPho || isEEPho) && phoEleVeto == 1 && isHggPho && phoCorrHggIDMVA > -0.9")
                .Define("isGoodPho",                "(isEBPho || isEEPho) && phoEleVeto == 1 && phoCorrHggIDMVA > -0.9")
                .Filter("ROOT::VecOps::Sum(isGoodPho) > 0", "event with good pho")
                .Define("phoIdx1",                  [ ](const ROOT::RVec<int>& isGoodPho,
                                                        const ROOT::RVec<float>& phoCalibEt){
                                                            // get the photon with highest et in good photons 
                                                            return utils::getIdx(isGoodPho, phoCalibEt)[0];
                                                        }, {"isGoodPho", "phoCalibEt"})
                .Define("phoCalibEt_Lead",          "phoCalibEt[phoIdx1]")
                .Define("phoSCEta_Lead",            "phoSCEta[phoIdx1]")
                .Define("phoCorrR9Full5x5_Lead",    "phoCorrR9Full5x5[phoIdx1]")
                .Define("phoCorrHggIDMVA_Lead",     "phoCorrHggIDMVA[phoIdx1]")
                .Define("phoEt_Lead",               "phoEt[phoIdx1]");

    // systematic variations for photon et
    std::vector<std::string> pho_sys;
    std::vector<std::string> pho_sys_bra;
    std::string r9sys = "";
    std::string r9sys_bra = "phoR9Full5x5[phoIdx1]";
    if (doVariation){
        for (size_t i = 0; i < variations.size(); i++){
            if (ba::contains(variations[i], "PhoScale") || ba::contains(variations[i], "PhoSigma")){
                pho_sys.emplace_back(variations[i]);
                pho_sys_bra.emplace_back(allowed_variations.at(variations[i])+"[phoIdx1]/cosh(phoEta[phoIdx1])");
            }
            if (ba::contains(variations[i], "PhoNoR9")){
                r9sys = variations[i];
            }
        }
        if (pho_sys.size() > 0 || r9sys != ""){
            printf(" [+] Book systematic variations for photon:\n");
            for (size_t i = 0; i < pho_sys.size(); i++){
                printf("     - %-21s: %s\n", pho_sys[i].c_str(), pho_sys_bra[i].c_str());
            }
            if (r9sys != "")
                printf("     - %-21s: %s\n", r9sys.c_str(), r9sys_bra.c_str());
            
            std::string pho_all_var = ba::join(pho_sys_bra, ", ");
            std::string varied_str  = Form("ROOT::RVec<float> v = {%s}; return v;", pho_all_var.c_str());
            pf = pf.Vary("phoCalibEt_Lead", varied_str, pho_sys, "pho");
            if (r9sys != "")
                pf = pf.Vary("phoCorrR9Full5x5_Lead", Form("ROOT::RVec<float> v = {%s}; return v;", r9sys_bra.c_str()), {r9sys}, "pho");
        }
    }
    pf = pf.Define("pho", "ROOT::Math::PtEtaPhiMVector v(phoCalibEt_Lead, phoEta[phoIdx1], phoPhi[phoIdx1], 0.); return v;");  

    //=======================================================//
    //                electron selections                    //
    //=======================================================//
    auto gf = gsf::DefineGSFColumns(pf); // create gsf tracks related columns

    // initialize the reader for merged electron ID
    XGBReader m2EBReader(cfg["external_files"]["mergedID_model"]["M2EB"].as<std::string>());
    XGBReader m2EEReader(cfg["external_files"]["mergedID_model"]["M2EE"].as<std::string>());
    float m2EBWPLoose = cfg["external_files"]["mergedID_wp"]["M2EBLoose"].as<float>();
    float m2EEWPLoose = cfg["external_files"]["mergedID_wp"]["M2EELoose"].as<float>();   
    std::vector<std::string> m2vals = {
        "eleRho",
        "eleSCEta",
        "eleSCRawEn",
        "eledEtaAtVtx",
        "eledPhiAtVtx",
        "elePtError",
        "eleHoverE",
        "eleEoverP",
        "eleEoverPout",
        "eleEoverPInv",
        "eleSCEtaWidth",
        "eleSCPhiWidth",
        "eleSigmaIEtaIEtaFull5x5",
        "eleSigmaIPhiIPhiFull5x5",
        "eleR9Full5x5",
        "eleBrem",
        "elePFChIso",
        "elePFPhoIso",
        "elePFNeuIso",
        "gsfPtRatio",
        "gsfDeltaR",
        "gsfRelPtRatio"
    };
    std::string m2vals_str = ba::join(m2vals, ", ");
    
    // merged electron energy regression
    auto energy_correction = cfg["external_files"]["energy_correction"];
    auto fregression    = energy_correction["energy_reg"].as<std::vector<std::string>>();
    auto fscale         = energy_correction["resi_scale"].as<std::vector<std::string>>();
    auto fsmear         = energy_correction["resi_smear"].as<std::vector<std::string>>();

    XGBReader regEBReader(fregression[0]);
    XGBReader regEEReader(fregression[1]);
    std::vector<std::string> regVals_EB = {
        "eleRho",
        "eleNVtx",
        "eleSCEta",
        "eleSCPhi",
        "eleSCRawEn",
        "eleCalibPt",
        "eledEtaAtVtx",
        "eledPhiAtVtx",
        "elePtError",
        "eleHoverE",
        "eleEoverP",
        "eleEoverPout",
        "eleEoverPInv",
        "eleSCEtaWidth",
        "eleSCPhiWidth",
        "eleSigmaIEtaIEtaFull5x5",
        "eleSigmaIPhiIPhiFull5x5",
        "eleR9Full5x5",
        "eleBrem",
        "gsfPtSum",
        "gsfPtRatio",
        "diTrkPt",
        "gsfDeltaR"
    };
    std::string regVals_EB_str = boost::algorithm::join(regVals_EB, ", ");
    std::string regVals_EE_str = regVals_EB_str + ", eleESEnToRawE";

    EnCalibrater calib;
    calib.SetScaleFiles(fscale);
    calib.SetSmearFiles(fsmear);

    auto ef = gf.Define("isEBEle",              "abs(eleSCEta) < 1.4442")
                .Define("isEEEle",              "abs(eleSCEta) > 1.566 && abs(eleSCEta) < 2.5")
                .Define("eleRho",               "ROOT::RVec<float> v(nEle, rho); return v;")
                .Define("eleNVtx",              "ROOT::RVec<float> v(nEle, nVtx); return v;")
                .Define("eleESEnToRawE",        "(eleESEnP1 + eleESEnP2) / eleSCRawEn")
                .Define("eleM2MVAVals",         Form("MakeMVAVals<%d, float>(%s)", (int)m2vals.size(),      m2vals_str.c_str()))
                .Define("eleRegEBMVAVals",      Form("MakeMVAVals<%d, float>(%s)", (int)regVals_EB.size(),  regVals_EB_str.c_str()))
                .Define("eleRegEEMVAVals",      Form("MakeMVAVals<%d, float>(%s)", (int)regVals_EB.size()+1,regVals_EE_str.c_str()))
                .Define("eleM2IDMVAs",          [&](const ROOT::RVec<float>& eleSCEta,
                                                    const ROOT::RVec<int>& nGsfMatchToReco,
                                                    const ROOT::RVec<std::vector<float>>& eleM2MVAVals){
                                                        ROOT::RVec<std::vector<float>> mvas;
                                                        for (size_t i = 0; i < eleSCEta.size(); i++){
                                                            std::vector<float> scores(3, -999);
                                                            if (nGsfMatchToReco[i] > 1){
                                                                if (fabs(eleSCEta[i]) < 1.479)
                                                                    scores = m2EBReader.Compute(eleM2MVAVals[i]);
                                                                else
                                                                    scores = m2EEReader.Compute(eleM2MVAVals[i]);
                                                            }
                                                            mvas.emplace_back(scores);
                                                        }
                                                        return mvas;
                                                    }, {"eleSCEta", "nGsfMatchToReco", "eleM2MVAVals"})
                .Define("isM2IDEle",            [&](const ROOT::RVec<float>& eleSCEta,
                                                    const ROOT::RVec<int>& nGsfMatchToReco, 
                                                    const ROOT::RVec<std::vector<float>>& eleM2IDMVAs){
                                                        ROOT::RVec<int> v(eleM2IDMVAs.size());
                                                        for (size_t i = 0; i < eleM2IDMVAs.size(); i++){
                                                            int class_with_max_score = TMath::LocMax(eleM2IDMVAs[i].size(), eleM2IDMVAs[i].data());
                                                            float wp = (fabs(eleSCEta[i]) < 1.479) ? m2EBWPLoose : m2EEWPLoose;
                                                            v[i] = (nGsfMatchToReco[i] > 1 && class_with_max_score == 0 && eleM2IDMVAs[i][0] > wp) ? 1 : 0;
                                                        }
                                                        return v;
                                                    }, {"eleSCEta", "nGsfMatchToReco", "eleM2IDMVAs"})
                .Define("eleHDALRegPt",         [&](const ROOT::RVec<float>& eleSCEta,
                                                    const ROOT::RVec<float>& elePt,
                                                    const ROOT::RVec<int>& nGsfMatchToReco,
                                                    const ROOT::RVec<std::vector<float>>& eleRegEBMVAVals,
                                                    const ROOT::RVec<std::vector<float>>& eleRegEEMVAVals){
                                                        ROOT::RVec<float> pt(elePt.size());
                                                        for (size_t i = 0; i < elePt.size(); i++){
                                                            float regSF = 1.;
                                                            if (nGsfMatchToReco[i] > 1){
                                                                if (fabs(eleSCEta[i]) < 1.479)
                                                                    regSF = regEBReader.Compute(eleRegEBMVAVals[i])[0];
                                                                else
                                                                    regSF = regEEReader.Compute(eleRegEEMVAVals[i])[0];
                                                            }
                                                            pt[i] = regSF * elePt[i];
                                                        }
                                                        return pt;
                                                    }, {"eleSCEta", "elePt", "nGsfMatchToReco", "eleRegEBMVAVals", "eleRegEEMVAVals"})
                // .Define("isHEEPEle",            [ ](const ROOT::RVec<UShort_t>& eleIDbit){
                //                                         auto v = ROOT::VecOps::Map(eleIDbit, [&](UShort_t idx){return (((eleIDbit[i] >> 4) & 1) == 1) ? (int) 1 : (int) 0; });
                //                                         return v;
                //                                     }, {"eleIDbit"})
                .Define("isHggEle",             eleSel::HggPresel,          {"nEle", "eleSCEta", "eleSCPhi", "nPho", "rhoAll", "phoSCEta", "phoSCPhi", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"})
                .Define("isPVTrk",              "(isEBEle && abs(eleTrkD0) < 0.02 && abs(eleSubTrkD0) < 0.02 && abs(eleTrkDz) < 0.1 && abs(eleSubTrkDz) < 0.1) || (isEEEle && abs(eleTrkD0) < 0.05 && abs(eleSubTrkD0) < 0.05 && abs(eleTrkDz) < 0.2 && abs(eleSubTrkDz) < 0.2)")
                .Define("isNonCovTrk",          "eleTrkMissHits < 1 && eleSubTrkMissHits < 1")
                // .Define("isGoodM2Ele",          "(isEBEle || isEEEle) && isM2IDEle && isHggEle && eleConvVeto == 1 && isPVTrk && isNonCovTrk")
                .Define("isGoodM2Ele",          "(isEBEle || isEEEle) && isM2IDEle && eleConvVeto == 1 && isPVTrk && isNonCovTrk")
                .Filter("ROOT::VecOps::Sum(isGoodM2Ele) > 0", "event with good ele");

    auto ef_scale = calib.GetScaleRDF(ef, isMC, "eleHDALRegPt", "eleSCEta", "eleHDALScalePt");
    auto ef_smear = calib.GetSmearRDF(ef_scale, isMC, "eleHDALScalePt", "eleSCEta", "eleHDALCalibPt");
    auto ef_calib = ef_smear.Define("eleIdx1",                  [ ](const ROOT::RVec<int>& isGoodM2Ele,
                                                                    const ROOT::RVec<float>& eleHDALCalibPt){
                                                                        return utils::getIdx(isGoodM2Ele, eleHDALCalibPt)[0];
                                                                    }, {"isGoodM2Ele", "eleHDALCalibPt"})
                            .Define("eleMathcedPhoIdx1",        [ ](const int eleIdx1,
                                                                    const ROOT::RVec<float> eleSCEta,
                                                                    const ROOT::RVec<float> eleSCPhi,
                                                                    const ROOT::RVec<float> phoSCEta,
                                                                    const ROOT::RVec<float> phoSCPhi){
                                                                        int idx = -1;
                                                                        for (size_t i = 0; i < phoSCEta.size(); i++){
                                                                            if ((eleSCEta[eleIdx1] == phoSCEta[i]) && (eleSCPhi[eleIdx1] == phoSCPhi[i])){
                                                                                idx = i;
                                                                                break;
                                                                            }
                                                                        }
                                                                        return idx;
                                                                    }, {"eleIdx1", "eleSCEta", "eleSCPhi", "phoSCEta", "phoSCPhi"})
                            .Define("eleMatchPhoCorrR9_Lead",   "phoCorrR9Full5x5[eleMathcedPhoIdx1]") // in order to apply the hgg preselected sf
                            .Define("eleMatchPhoSCEta_Lead",    "phoSCEta[eleMathcedPhoIdx1]")         // in order to apply the hgg preselected sf
                            .Define("eleHDALCalibPt_Lead",      "eleHDALCalibPt[eleIdx1]")
                            .Define("elePt_Lead",               "elePt[eleIdx1]")
                            .Define("eleCalibPt_Lead",          "eleCalibPt[eleIdx1]")
                            .Define("eleTrkPt_Lead",            "eleTrkPt[eleIdx1]")
                            .Define("eleSubTrkPt_Lead",         "eleSubTrkPt[eleIdx1]")
                            .Define("eleSCEta_Lead",            "eleSCEta[eleIdx1]")
                            .Define("eleAbsSCEta_Lead",         "fabs(eleSCEta[eleIdx1])")
                            .Define("eleTrkEta_Lead",           "eleTrkEta[eleIdx1]")
                            .Define("eleSubTrkEta_Lead",        "eleSubTrkEta[eleIdx1]")
                            .Define("eleM2IDMVA_Lead",          "eleM2IDMVAs[eleIdx1][0]")
                            .Define("eleHDALRegPt_Lead",        "eleHDALRegPt[eleIdx1]");

    // systematic variations for electron Pt
    std::vector<std::string> ele_sys;
    std::vector<std::string> ele_sys_bra;
    if (doVariation){
        for (size_t i = 0; i < variations.size(); i++){
            if (ba::contains(variations[i], "EleHDALScale")){
                ele_sys.emplace_back(variations[i]);
                ele_sys_bra.emplace_back("eleHDALCalibPt_Lead "+allowed_variations.at(variations[i])+"[eleIdx1] * eleHDALRegPt[eleIdx1]");
            }
            if (ba::contains(variations[i], "EleHDALSmear")){
                ele_sys.emplace_back(variations[i]);
                ele_sys_bra.emplace_back(allowed_variations.at(variations[i])+"[eleIdx1]");
            }
            if (ba::contains(variations[i], "EleScale") || ba::contains(variations[i], "EleSigma")){
                ele_sys.emplace_back(variations[i]);
                ele_sys_bra.emplace_back(allowed_variations.at(variations[i])+"[eleIdx1]/cosh(eleEta[eleIdx1])");
            }
        }
        if (ele_sys.size() > 0){
            printf(" [+] Book systematic variations for electron:\n");
            for (size_t i = 0; i < ele_sys.size(); i++){
                printf("     - %-21s: %s\n", ele_sys[i].c_str(), ele_sys_bra[i].c_str());
            }
            std::string ele_all_var = ba::join(ele_sys_bra, ", ");
            std::string varied_str  = Form("ROOT::RVec<float> v = {%s}; return v;", ele_all_var.c_str());
            ef_calib = ef_calib.Vary("eleHDALCalibPt_Lead", varied_str, ele_sys, "ele");
        }
    }
    // use the mass reconstructed by two gsf track to be merged electron mass
    ef_calib =  ef_calib.Define("gsf1",                     "eleTrk1[eleIdx1]")
                        .Define("gsf2",                     "eleTrk2[eleIdx1]")
                        .Define("ele1",                     "ROOT::Math::PtEtaPhiMVector v(eleHDALCalibPt_Lead, eleEta[eleIdx1], elePhi[eleIdx1], (float) (gsf1+gsf2).M()); return v;");

    //=======================================================//
    //            kinematic event selections                 //
    //=======================================================//
    auto kf = ef_calib.Define("H",                "ele1 + pho")
                    //   .Filter("TMath::Max(phoCalibEt_Lead, eleHDALCalibPt_Lead) > 35 && TMath::Min(phoCalibEt_Lead, eleHDALCalibPt_Lead) > 25", "trigger threshold")
                      .Filter("TMath::Max(phoCalibEt_Lead, eleHDALCalibPt_Lead) > 30 && TMath::Min(phoCalibEt_Lead, eleHDALCalibPt_Lead) > 30", "trigger threshold")
                      .Filter("(gsf1+gsf2).M() < 50", "Mgg < 50 GeV")
                      .Filter("(gsf1+gsf2).M() > 11 || (gsf1+gsf2).M() < 8", "reject Upsilon")
                      .Filter("(gsf1+gsf2).M() > 3.5 || (gsf1+gsf2).M() < 2.5", "reject Jpsi")
                      .Filter("(gsf1.Pt() + gsf2.Pt()) > 44", "gsfPtSum > 44 GeV")
                      .Filter("(ele1.Pt()/H.M()) > 0.3 && (pho.Pt()/H.M()) > 0.3", "pt mass ratio cut")
                      .Filter("ROOT::VecOps::DeltaR(ele1.Eta(), pho.Eta(), ele1.Phi(), pho.Phi()) > 1", "dR cut")
                      .Filter("H.M() > 110. && H.M() < 170.", "three body mass cut");
    //! End of the selections

    //=======================================================//
    //                      categorization                   //
    //=======================================================//
    // DiJet tagged (VBF tagged) -> events with two jets passing:
    //*     VBF: Mass(jet1, jet2) > 500
    //*     1) Loose Jet ID, jetPt > 30, |jetEta| < 4.7 (per jet sel)
    //*     2) dR(jet, leps) > 0.4, dR(jet, pho) > 0.4 (per jet sel)
    //*     3) |dEta(jet1, jet2)| > 3.5
    //*     4) Eta(eeg) - ((jetEta1 + jetEta2) * 0.5) < 2.5 (zeppen)
    //*     5) |dPhi(eeg, jj)| > 2.4
    // Boosted tagged: pT(eeg) > 60 GeV
    // Untagged: EBHR9, EBLR9, EE
    
    if (isMC){
        kf = kf.Define("jetSmearPt", "jetP4Smear * jetPt")
               .Define("jetSmearEn", "jetP4Smear * jetEn");
    }
    else{
        kf = kf.Define("jetSmearPt", "jetPt")
               .Define("jetSmearEn", "jetEn");
    }

    // systematic variations for jet pt & en
    std::vector<std::string> jet_sys;
    std::vector<std::string> jet_sys_bra1; // pt
    std::vector<std::string> jet_sys_bra2; // en
    if (doVariation){
        for (size_t i = 0; i < variations.size(); i++){
            if (ba::contains(variations[i], "JER")){
                jet_sys.emplace_back(variations[i]);
                jet_sys_bra1.emplace_back(allowed_variations.at(variations[i])+" * jetPt"); // eg. jetP4SmearUp * jetP
                jet_sys_bra2.emplace_back(allowed_variations.at(variations[i])+" * jetEn");
            }
            if (ba::contains(variations[i], "JEC")){
                jet_sys.emplace_back(variations[i]);
                jet_sys_bra1.emplace_back("jetP4Smear * jetPt"+allowed_variations.at(variations[i])); // eg. jetP4Smear * jetPt + jetJECUnc
                jet_sys_bra2.emplace_back("jetP4Smear * jetEn"+allowed_variations.at(variations[i]));
            }
        }
        if (jet_sys.size() > 0){
            printf(" [+] Book systematic variations for jet:\n");
            for (size_t i = 0; i < jet_sys.size(); i++){
                printf("     - %-21s: %s && %s\n", jet_sys[i].c_str(), jet_sys_bra1[i].c_str(), jet_sys_bra2[i].c_str());
            }
            std::string jet_all_var1 = ba::join(jet_sys_bra1, ", ");
            std::string jet_all_var2 = ba::join(jet_sys_bra2, ", ");
            std::string varied_str  = Form("ROOT::RVec<ROOT::RVec<ROOT::RVec<float>>> v = {{%s}, {%s}}; return v;", jet_all_var1.c_str(), jet_all_var2.c_str());
            kf = kf.Vary({"jetSmearPt", "jetSmearEn"}, varied_str, jet_sys, "jet");
        }
        // .Vary({"jetPt", "jetEn"},   [ ](const ROOT::RVec<float>& jetPt,
        //                                 const ROOT::RVec<float>& jetEn){
        //                                     ROOT::RVec<ROOT::RVec<ROOT::RVec<float>>> v = {{jetPt*1.1}, {jetEn*1.1}};
        //                                     return v;
        //                                 }, {"jetPt", "jetEn"}, {"test"}, "jet")
    }
    
    float m2EBWPTight = cfg["external_files"]["mergedID_wp"]["M2EBTight"].as<float>();
    float m2EEWPTight = cfg["external_files"]["mergedID_wp"]["M2EETight"].as<float>();  
    auto cf = kf.Define("dR_ele1_jet",          cat::dRVector,             {"ele1", "jetEta", "jetPhi"})
                .Define("dR_pho_jet",           cat::dRVector,             {"pho", "jetEta", "jetPhi"})
                .Define("passJetID",            cat::CutBasedJet,          {"nJet", "jetEta", "jetNHF", "jetNEF", "jetNNP", "jetNCH", "jetCHF", "jetCEF", "jetMUF", "era"})
                .Define("isGoodJet",            "jetSmearPt > 30. && abs(jetEta) < 4.7 && passJetID && dR_ele1_jet > 0.4 && dR_pho_jet > 0.4")
                .Define("vbftag",               cat::VBFtag,               {"isGoodJet", "jetSmearPt", "jetEta", "jetPhi", "jetSmearEn", "H"})
                .Define("isVbf",                "vbftag[0] == 1")
                .Define("jetIdx1",              "vbftag[1]")
                .Define("jetIdx2",              "vbftag[2]")
                .Define("isBst",                "H.Pt() > 60.")
                .Define("isEBHR9",              "abs(phoSCEta[phoIdx1]) < 1.4442 && phoCorrR9Full5x5_Lead > 0.96")
                .Define("isEBLR9",              "abs(phoSCEta[phoIdx1]) < 1.4442 && phoCorrR9Full5x5_Lead <= 0.96")
                .Define("isEE",                 "abs(phoSCEta[phoIdx1]) > 1.566  && abs(phoSCEta[phoIdx1]) < 2.5")
                .Define("isWPTight",            Form("(isEBEle[eleIdx1] && eleM2IDMVAs[eleIdx1][0] > %f) || (isEEEle[eleIdx1] && eleM2IDMVAs[eleIdx1][0] > %f)", m2EBWPTight, m2EEWPTight))
                .Define("category",             [&](const bool isWPTight,
                                                    const bool isVbf,
                                                    const bool isBst,
                                                    const bool isEBHR9,
                                                    const bool isEBLR9,
                                                    const bool isEE){
                                                        int cat = 0;
                                                        if (!isWPTight)     return cat = 1; // tag the events with relatively low quality of M2
                                                        else if (isVbf)     return cat = 2;
                                                        else if (isBst)     return cat = 3;
                                                        else if (isEBHR9)   return cat = 4;
                                                        else if (isEBLR9)   return cat = 5;
                                                        else if (isEE)      return cat = 6;
                                                        else throw std::runtime_error("No proper category for M2");
                                                    }, {"isWPTight", "isVbf", "isBst", "isEBHR9", "isEBLR9", "isEE"})
                .Define("category_str",         [&](int category){return categories.at(category-1);},      {"category"});

    //=======================================================//
    //                  scale factors for MC                 //
    //=======================================================//
    ROOT::RDF::RNode df_sfs = cf;
    std::map<std::string, std::shared_ptr<WeightHandler>> Weis;
    bool apply_sfs = (isMC && cfg["external_files"]["scale_factors"]) ? true : false;
    if (apply_sfs){
        printf(" [+] Apply the scale factors:\n");

        auto sfcfgs = cfg["external_files"]["scale_factors"];
        for (size_t i = 0; i < scale_factors.size(); i++){
            auto SFName = scale_factors[i];
            if (sfcfgs[SFName]){
                Weis[SFName] = std::make_shared<WeightHandler>();
                Weis[SFName]->Init(SFName, sfcfgs[SFName]);
            }
            else{
                Warn(Form("     - WARNING: Scale factors of %s are required but missing in config file, which are assigned as 1.\n", SFName.c_str()));
                Weis[SFName] = nullptr;
            }
        }

        std::vector<std::string> SFs;
        for (auto it = Weis.begin(); it != Weis.end(); it++){
            if (ba::contains(it->first, "DiPhoHLT")) // a little bit tricky
                continue;
            
            if (it->second){
                if (ba::contains(it->first, "Trk")){ // in bins of trk pt and trk eta
                    std::string xbin = (ba::contains(it->first, "SubTrk")) ? "eleSubTrkPt_Lead"  : "eleTrkPt_Lead";
                    std::string ybin = (ba::contains(it->first, "SubTrk")) ? "eleSubTrkEta_Lead" : "eleTrkEta_Lead";
                    df_sfs = it->second->GetRDF(df_sfs, it->first+"SF_Lead", it->first+"SFErr_Lead", {xbin, ybin});
                    
                    printf("     - %-21s: in bins of %-22s && %-22s\n", (it->first).c_str(), xbin.c_str(), ybin.c_str());
                }
                if (ba::contains(it->first, "Ele")){ 
                    std::string xbin = "elePt_Lead";    // in bins of sceta and pt
                    std::string ybin = "eleSCEta_Lead";
                    if (ba::contains(it->first, "Merged")) { 
                        xbin = "eleCalibPt_Lead";
                        ybin = "eleAbsSCEta_Lead";
                    }
                    if (ba::contains(it->first, "HggPresel")) {
                        xbin = "eleMatchPhoSCEta_Lead";
                        ybin = "eleMatchPhoCorrR9_Lead";
                    }
                    df_sfs = it->second->GetRDF(df_sfs, it->first+"SF_Lead", it->first+"SFErr_Lead", {xbin, ybin});

                    printf("     - %-21s: in bins of %-22s && %-22s\n", (it->first).c_str(), xbin.c_str(), ybin.c_str());
                }
                if (ba::contains(it->first, "Pho")){ // in bins of sceta and pt
                    std::string xbin = "phoSCEta_Lead";
                    std::string ybin = "phoCorrR9Full5x5_Lead";
                    df_sfs = it->second->GetRDF(df_sfs, it->first+"SF_Lead", it->first+"SFErr_Lead", {xbin, ybin});

                    printf("     - %-21s: in bins of %-22s && %-22s\n", (it->first).c_str(), xbin.c_str(), ybin.c_str());
                }
            }
            else{
                df_sfs = df_sfs.Define(it->first+"SF_Lead",    "(float) 1.")
                               .Define(it->first+"SFErr_Lead", "(float) 0.");
            }

            SFs.emplace_back(it->first+"SF_Lead");
        }
        if (Weis["DiPhoHLTSeedForPho"] && Weis["DiPhoHLTSeedForEle"] && Weis["DiPhoHLTUnSeedForPho"] && Weis["DiPhoHLTUnSeedForEle"]){
            df_sfs = df_sfs.Define("DiPhoHLTSeedSF_Lead",       [&](const float phoEt_Lead,
                                                                    const float phoSCEta_Lead,
                                                                    const float phoCorrR9Full5x5_Lead,
                                                                    const float elePt_Lead,
                                                                    const float eleSCEta_Lead){
                                                                        float SF = 1.;
                                                                        if (phoEt_Lead > elePt_Lead)
                                                                            SF = Weis["DiPhoHLTSeedForPho"]->GetSFFromHgg(phoSCEta_Lead, phoCorrR9Full5x5_Lead, phoEt_Lead);
                                                                        else
                                                                            SF = Weis["DiPhoHLTSeedForEle"]->GetSFFromEGM(eleSCEta_Lead, elePt_Lead);
                                                                        return SF;
                                                                    }, {"phoEt_Lead", "phoSCEta_Lead", "phoCorrR9Full5x5_Lead", "elePt_Lead", "eleSCEta_Lead"})
                        .Define("DiPhoHLTSeedSFErr_Lead",       [&](const float phoEt_Lead,
                                                                    const float phoSCEta_Lead,
                                                                    const float phoCorrR9Full5x5_Lead,
                                                                    const float elePt_Lead,
                                                                    const float eleSCEta_Lead){
                                                                        float SFerr = 0;
                                                                        if (phoEt_Lead > elePt_Lead)
                                                                            SFerr = Weis["DiPhoHLTSeedForPho"]->GetSFErrFromHgg(phoSCEta_Lead, phoCorrR9Full5x5_Lead, phoEt_Lead);
                                                                        else
                                                                            SFerr = Weis["DiPhoHLTSeedForEle"]->GetSFErrFromEGM(eleSCEta_Lead, elePt_Lead);
                                                                        return SFerr;
                                                                    }, {"phoEt_Lead", "phoSCEta_Lead", "phoCorrR9Full5x5_Lead", "elePt_Lead", "eleSCEta_Lead"})
                        .Define("DiPhoHLTUnSeedSF_Lead",        [&](const float phoEt_Lead,
                                                                    const float phoSCEta_Lead,
                                                                    const float phoCorrR9Full5x5_Lead,
                                                                    const float elePt_Lead,
                                                                    const float eleSCEta_Lead){
                                                                        float SF = 1.;
                                                                        if (phoEt_Lead < elePt_Lead)
                                                                            SF = Weis["DiPhoHLTUnSeedForPho"]->GetSFFromHgg(phoSCEta_Lead, phoCorrR9Full5x5_Lead, phoEt_Lead);
                                                                        else
                                                                            SF = Weis["DiPhoHLTUnSeedForEle"]->GetSFFromEGM(eleSCEta_Lead, elePt_Lead);
                                                                        return SF;
                                                                    }, {"phoEt_Lead", "phoSCEta_Lead", "phoCorrR9Full5x5_Lead", "elePt_Lead", "eleSCEta_Lead"})
                        .Define("DiPhoHLTUnSeedSFErr_Lead",     [&](const float phoEt_Lead,
                                                                    const float phoSCEta_Lead,
                                                                    const float phoCorrR9Full5x5_Lead,
                                                                    const float elePt_Lead,
                                                                    const float eleSCEta_Lead){
                                                                        float SFerr = 0;
                                                                        if (phoEt_Lead < elePt_Lead)
                                                                            SFerr = Weis["DiPhoHLTUnSeedForPho"]->GetSFErrFromHgg(phoSCEta_Lead, phoCorrR9Full5x5_Lead, phoEt_Lead);
                                                                        else
                                                                            SFerr = Weis["DiPhoHLTUnSeedForEle"]->GetSFErrFromEGM(eleSCEta_Lead, elePt_Lead);
                                                                        return SFerr;
                                                                    }, {"phoEt_Lead", "phoSCEta_Lead", "phoCorrR9Full5x5_Lead", "elePt_Lead", "eleSCEta_Lead"});
            printf("     - %-21s: in bins of photon and electron pt eta r9\n", "DiPhoHLTSeed");
            printf("     - %-21s: in bins of photon and electron pt eta r9\n", "DiPhoHLTUnSeed");
            SFs.emplace_back("DiPhoHLTSeedSF_Lead");
            SFs.emplace_back("DiPhoHLTUnSeedSF_Lead");
        }
        else{
            df_sfs = df_sfs.Define("DiPhoHLTSeedSF_Lead",       "(float) 1.")
                           .Define("DiPhoHLTSeedSFErr_Lead",    "(float) 0.")
                           .Define("DiPhoHLTUnSeedSF_Lead",     "(float) 1.")
                           .Define("DiPhoHLTUnSeedSFErr_Lead",  "(float) 0.");
            SFs.emplace_back("DiPhoHLTSeedSF_Lead");
            SFs.emplace_back("DiPhoHLTUnSeedSF_Lead");
        }

        // apply the final weights
        std::vector<std::string> SFsUp;
        std::vector<std::string> SFsDo;
        for (size_t i = 0; i < SFs.size(); i++){
            TString SFErrStr(SFs[i]);
            SFErrStr.ReplaceAll("SF", "SFErr");
            std::string SFUp_str = Form("(%s + %s)", SFs[i].c_str(), SFErrStr.Data());
            std::string SFDo_str = Form("(%s - %s)", SFs[i].c_str(), SFErrStr.Data());
            SFsUp.emplace_back(SFUp_str);
            SFsDo.emplace_back(SFDo_str);
        }

        std::map<std::string, std::string> Weights_str;
        Weights_str["nominal"] = "mcwei * genwei * puwei         * L1ECALPrefire";
        Weights_str["puweiUp"] = "mcwei * genwei * puwei_up      * L1ECALPrefire";
        Weights_str["puweiDo"] = "mcwei * genwei * puwei_down    * L1ECALPrefire";
        Weights_str["L1PreUp"] = "mcwei * genwei * puwei         * L1ECALPrefireUp";
        Weights_str["L1PreDo"] = "mcwei * genwei * puwei         * L1ECALPrefireDown";
        Weights_str["HLTUp"]   = "mcwei * genwei * puwei         * L1ECALPrefire";
        Weights_str["HLTDo"]   = "mcwei * genwei * puwei         * L1ECALPrefire";
        Weights_str["EleIDUp"] = "mcwei * genwei * puwei         * L1ECALPrefire";
        Weights_str["EleIDDo"] = "mcwei * genwei * puwei         * L1ECALPrefire";
        Weights_str["PhoIDUp"] = "mcwei * genwei * puwei         * L1ECALPrefire";
        Weights_str["PhoIDDo"] = "mcwei * genwei * puwei         * L1ECALPrefire";
        for (auto it = Weights_str.begin(); it != Weights_str.end(); it++){
            for (size_t i = 0; i < SFs.size(); i++){
                std::string wei_str = Form(" * %s", SFs[i].c_str());
                if (it->first == "HLTUp"   && (ba::contains(SFs[i], "DiPhoHLTSeedSF") || ba::contains(SFs[i], "DiPhoHLTUnSeedSF")))
                    wei_str = Form(" * %s", SFsUp[i].c_str());
                
                if (it->first == "HLTDo"   && (ba::contains(SFs[i], "DiPhoHLTSeedSF") || ba::contains(SFs[i], "DiPhoHLTUnSeedSF")))
                    wei_str = Form(" * %s", SFsDo[i].c_str());
                
                if (it->first == "EleIDUp" && (ba::contains(SFs[i], "Ele")  || ba::contains(SFs[i], "Trk")))
                    wei_str = Form(" * %s", SFsUp[i].c_str());
                
                if (it->first == "EleIDDo" && (ba::contains(SFs[i], "Ele")  || ba::contains(SFs[i], "Trk")))
                    wei_str = Form(" * %s", SFsDo[i].c_str());
                
                if (it->first == "PhoIDUp" && (ba::contains(SFs[i], "Pho")  && !ba::contains(SFs[i], "HLT")))
                    wei_str = Form(" * %s", SFsUp[i].c_str());
                
                if (it->first == "PhoIDDo" && (ba::contains(SFs[i], "Pho")  && !ba::contains(SFs[i], "HLT")))
                    wei_str = Form(" * %s", SFsDo[i].c_str());
                
                it->second += wei_str;
            }
        }
        
        for (auto it = Weights_str.begin(); it != Weights_str.end(); it++){
            std::string wei_col = (it->first == "nominal") ? "weight" : "weight_"+it->first;
            df_sfs = df_sfs.Define(wei_col, it->second);
        }
    }

    //=======================================================//
    //             check generater infomation                //
    //=======================================================//
    if (isSignalMC){
        df_sfs = df_sfs.Define("genIdx_reco1",              gen::MatchedGenEle,  {"ele1", "nMC", "mcEta", "mcPhi", "mcPID", "mcMomPID", "mcGMomPID", "mcStatusFlag"})
                       .Define("nReco1MatchedGen",          "genIdx_reco1.size()")  // number of generator electrons can be matched to selected reco electron
                       .Define("isTrueM2",                  "nGsfMatchToReco[eleIdx1] > 1 && nReco1MatchedGen > 1")
                       .Define("mcPIDForEle_Lead",          "if (nReco1MatchedGen > 0) return mcPID[genIdx_reco1[0]]; else return (int) 0;")
                       .Define("mcPtForEle_Lead",           "if (nReco1MatchedGen > 0) return mcPt[genIdx_reco1[0]]; else return (float) 0;")
                       .Define("mcMassForEle_Lead",         "if (nReco1MatchedGen > 0) return mcMass[genIdx_reco1[0]]; else return (float) 0;")
                       .Define("mcEtaForEle_Lead",          "if (nReco1MatchedGen > 0) return mcEta[genIdx_reco1[0]]; else return (float) 0;")
                       .Define("mcPhiForEle_Lead",          "if (nReco1MatchedGen > 0) return mcPhi[genIdx_reco1[0]]; else return (float) 0;")
                       .Define("mcMomPIDForEle_Lead",       "if (nReco1MatchedGen > 0) return mcMomPID[genIdx_reco1[0]]; else return (int) 0;")
                       .Define("mcMomPtForEle_Lead",        "if (nReco1MatchedGen > 0) return mcMomPt[genIdx_reco1[0]]; else return (float) 0;")
                       .Define("mcMomMassForEle_Lead",      "if (nReco1MatchedGen > 0) return mcMomMass[genIdx_reco1[0]]; else return (float) 0;")
                       .Define("mcMomEtaForEle_Lead",       "if (nReco1MatchedGen > 0) return mcMomEta[genIdx_reco1[0]]; else return (float) 0;")
                       .Define("mcMomPhiForEle_Lead",       "if (nReco1MatchedGen > 0) return mcMomPhi[genIdx_reco1[0]]; else return (float) 0;")
                       .Define("mcGMomPIDForEle_Lead",      "if (nReco1MatchedGen > 0) return mcGMomPID[genIdx_reco1[0]]; else return (int) 0;")
                       .Define("mcPIDForEle_subLead",       "if (nReco1MatchedGen > 1) return mcPID[genIdx_reco1[1]]; else return (int) 0;")
                       .Define("mcPtForEle_subLead",        "if (nReco1MatchedGen > 1) return mcPt[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMassForEle_subLead",      "if (nReco1MatchedGen > 1) return mcMass[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcEtaForEle_subLead",       "if (nReco1MatchedGen > 1) return mcEta[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcPhiForEle_subLead",       "if (nReco1MatchedGen > 1) return mcPhi[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMomPIDForEle_subLead",    "if (nReco1MatchedGen > 1) return mcMomPID[genIdx_reco1[1]]; else return (int) 0;")
                       .Define("mcMomPtForEle_subLead",     "if (nReco1MatchedGen > 1) return mcMomPt[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMomMassForEle_subLead",   "if (nReco1MatchedGen > 1) return mcMomMass[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMomEtaForEle_subLead",    "if (nReco1MatchedGen > 1) return mcMomEta[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMomPhiForEle_subLead",    "if (nReco1MatchedGen > 1) return mcMomPhi[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcGMomPIDForEle_subLead",   "if (nReco1MatchedGen > 1) return mcGMomPID[genIdx_reco1[1]]; else return (int) 0;")
                       // gen photon 
                       .Define("genPhoIdx",                 gen::MatchedGenPho,  {"pho", "nMC", "mcEta", "mcPhi", "mcPID", "mcMomPID", "mcStatusFlag"})
                       .Define("mcPIDForPho_Lead",          "if (genPhoIdx != -1) return mcPID[genPhoIdx]; else return (int) 0;")
                       .Define("mcPtForPho_Lead",           "if (genPhoIdx != -1) return mcPt[genPhoIdx]; else return (float) 0;")
                       .Define("mcMassForPho_Lead",         "if (genPhoIdx != -1) return mcMass[genPhoIdx]; else return (float) 0;")
                       .Define("mcEtaForPho_Lead",          "if (genPhoIdx != -1) return mcEta[genPhoIdx]; else return (float) 0;")
                       .Define("mcPhiForPho_Lead",          "if (genPhoIdx != -1) return mcPhi[genPhoIdx]; else return (float) 0;")
                       .Define("mcMomPIDForPho_Lead",       "if (genPhoIdx != -1) return mcMomPID[genPhoIdx]; else return (int) 0;")
                       .Define("mcMomPtForPho_Lead",        "if (genPhoIdx != -1) return mcMomPt[genPhoIdx]; else return (float) 0;")
                       .Define("mcMomMassForPho_Lead",      "if (genPhoIdx != -1) return mcMomMass[genPhoIdx]; else return (float) 0;")
                       .Define("mcMomEtaForPho_Lead",       "if (genPhoIdx != -1) return mcMomEta[genPhoIdx]; else return (float) 0;")
                       .Define("mcMomPhiForPho_Lead",       "if (genPhoIdx != -1) return mcMomPhi[genPhoIdx]; else return (float) 0;")
                       .Define("mcGMomPIDForPho_Lead",      "if (genPhoIdx != -1) return mcGMomPID[genPhoIdx]; else return (int) 0;");
    }
    else if (isMC){
        df_sfs = df_sfs.Define("genIdx_reco1",              gen::FindGenParticle,  {"ele1", "nMC", "mcEta", "mcPhi", "mcPt"})
                       .Define("mcPIDForEle_Lead",          "if (genIdx_reco1 != -1) return mcPID[genIdx_reco1]; else return (int) 0;")
                       .Define("mcPtForEle_Lead",           "if (genIdx_reco1 != -1) return mcPt[genIdx_reco1]; else return (float) 0;")
                       .Define("mcMassForEle_Lead",         "if (genIdx_reco1 != -1) return mcMass[genIdx_reco1]; else return (float) 0;")
                       .Define("mcEtaForEle_Lead",          "if (genIdx_reco1 != -1) return mcEta[genIdx_reco1]; else return (float) 0;")
                       .Define("mcPhiForEle_Lead",          "if (genIdx_reco1 != -1) return mcPhi[genIdx_reco1]; else return (float) 0;")
                       .Define("mcMomPIDForEle_Lead",       "if (genIdx_reco1 != -1) return mcMomPID[genIdx_reco1]; else return (int) 0;")
                       .Define("mcMomPtForEle_Lead",        "if (genIdx_reco1 != -1) return mcMomPt[genIdx_reco1]; else return (float) 0;")
                       .Define("mcMomMassForEle_Lead",      "if (genIdx_reco1 != -1) return mcMomMass[genIdx_reco1]; else return (float) 0;")
                       .Define("mcMomEtaForEle_Lead",       "if (genIdx_reco1 != -1) return mcMomEta[genIdx_reco1]; else return (float) 0;")
                       .Define("mcMomPhiForEle_Lead",       "if (genIdx_reco1 != -1) return mcMomPhi[genIdx_reco1]; else return (float) 0;")
                       .Define("mcGMomPIDForEle_Lead",      "if (genIdx_reco1 != -1) return mcGMomPID[genIdx_reco1]; else return (int) 0;")

                       // gen photon 
                       .Define("genPhoIdx",                 gen::FindGenParticle,  {"pho", "nMC", "mcEta", "mcPhi", "mcPt"})
                       .Define("mcPIDForPho_Lead",          "if (genPhoIdx != -1) return mcPID[genPhoIdx]; else return (int) 0;")
                       .Define("mcPtForPho_Lead",           "if (genPhoIdx != -1) return mcPt[genPhoIdx]; else return (float) 0;")
                       .Define("mcMassForPho_Lead",         "if (genPhoIdx != -1) return mcMass[genPhoIdx]; else return (float) 0;")
                       .Define("mcEtaForPho_Lead",          "if (genPhoIdx != -1) return mcEta[genPhoIdx]; else return (float) 0;")
                       .Define("mcPhiForPho_Lead",          "if (genPhoIdx != -1) return mcPhi[genPhoIdx]; else return (float) 0;")
                       .Define("mcMomPIDForPho_Lead",       "if (genPhoIdx != -1) return mcMomPID[genPhoIdx]; else return (int) 0;")
                       .Define("mcMomPtForPho_Lead",        "if (genPhoIdx != -1) return mcMomPt[genPhoIdx]; else return (float) 0;")
                       .Define("mcMomMassForPho_Lead",      "if (genPhoIdx != -1) return mcMomMass[genPhoIdx]; else return (float) 0;")
                       .Define("mcMomEtaForPho_Lead",       "if (genPhoIdx != -1) return mcMomEta[genPhoIdx]; else return (float) 0;")
                       .Define("mcMomPhiForPho_Lead",       "if (genPhoIdx != -1) return mcMomPhi[genPhoIdx]; else return (float) 0;")
                       .Define("mcGMomPIDForPho_Lead",      "if (genPhoIdx != -1) return mcGMomPID[genPhoIdx]; else return (int) 0;");
    }

    //=======================================================//
    //             define the final branches                 //
    //=======================================================//
    auto df_Final  = df_sfs.Define("CMS_higgs_mass",                "(float) H.M()")
                           // electron variables 
                           .Define("eleCharge_Lead",                "eleCharge[eleIdx1]")
                           .Define("eleEn_Lead",                    "eleEn[eleIdx1]")
                           .Define("eleSCEn_Lead",                  "eleSCEn[eleIdx1]")
                           .Define("eleEcalEn_Lead",                "eleEcalEn[eleIdx1]")
                           .Define("eleESEnToRawE_Lead",            "eleESEnToRawE[eleIdx1]")
                           .Define("eleD0_Lead",                    "eleD0[eleIdx1]")
                           .Define("eleDz_Lead",                    "eleDz[eleIdx1]")
                           .Define("eleSIP_Lead",                   "eleSIP[eleIdx1]")
                           .Define("elePtError_Lead",               "elePtError[eleIdx1]")
                           .Define("eleEta_Lead",                   "eleEta[eleIdx1]")
                           .Define("elePhi_Lead",                   "elePhi[eleIdx1]")
                           .Define("eleSCPhi_Lead",                 "eleSCPhi[eleIdx1]")
                           .Define("eleSCRawEn_Lead",               "eleSCRawEn[eleIdx1]")
                           .Define("eleSCEtaWidth_Lead",            "eleSCEtaWidth[eleIdx1]")
                           .Define("eleSCPhiWidth_Lead",            "eleSCPhiWidth[eleIdx1]")
                           .Define("eleHoverE_Lead",                "eleHoverE[eleIdx1]")
                           .Define("eleEoverP_Lead",                "eleEoverP[eleIdx1]")
                           .Define("eleEoverPout_Lead",             "eleEoverPout[eleIdx1]")
                           .Define("eleEoverPInv_Lead",             "eleEoverPInv[eleIdx1]")
                           .Define("eleBrem_Lead",                  "eleBrem[eleIdx1]")
                           .Define("eledEtaAtVtx_Lead",             "eledEtaAtVtx[eleIdx1]")
                           .Define("eledPhiAtVtx_Lead",             "eledPhiAtVtx[eleIdx1]")
                           .Define("eleSigmaIEtaIEtaFull5x5_Lead",  "eleSigmaIEtaIEtaFull5x5[eleIdx1]")
                           .Define("eleSigmaIPhiIPhiFull5x5_Lead",  "eleSigmaIPhiIPhiFull5x5[eleIdx1]")
                           .Define("eleConvVeto_Lead",              "eleConvVeto[eleIdx1]")
                           .Define("eleMissHits_Lead",              "eleMissHits[eleIdx1]")
                           .Define("eleESEffSigmaRR_Lead",          "eleESEffSigmaRR[eleIdx1]")
                           .Define("elePFChIso_Lead",               "elePFChIso[eleIdx1]")
                           .Define("elePFPhoIso_Lead",              "elePFPhoIso[eleIdx1]")
                           .Define("elePFNeuIso_Lead",              "elePFNeuIso[eleIdx1]")
                           .Define("elePFPUIso_Lead",               "elePFPUIso[eleIdx1]")
                           .Define("eleIDMVAIso_Lead",              "eleIDMVAIso[eleIdx1]")
                           .Define("eleIDMVANoIso_Lead",            "eleIDMVANoIso[eleIdx1]")
                           .Define("eleR9Full5x5_Lead",             "eleR9Full5x5[eleIdx1]")
                           .Define("eleGSFChi2_Lead",               "eleGSFChi2[eleIdx1]")
                           .Define("eleFiredSingleTrgs_Lead",       "eleFiredSingleTrgs[eleIdx1]")
                           .Define("eleFiredDoubleTrgs_Lead",       "eleFiredDoubleTrgs[eleIdx1]")
                           .Define("eleFiredL1Trgs_Lead",           "eleFiredL1Trgs[eleIdx1]")
                           .Define("eleTrkPhi_Lead",                "eleTrkPhi[eleIdx1]")
                           .Define("eleTrkCharge_Lead",             "eleTrkCharge[eleIdx1]")
                           .Define("eleTrkLayers_Lead",             "eleTrkLayers[eleIdx1]")
                           .Define("eleTrkMissHits_Lead",           "eleTrkMissHits[eleIdx1]")
                           .Define("eleTrkD0_Lead",                 "eleTrkD0[eleIdx1]")
                           .Define("eleTrkDz_Lead",                 "eleTrkDz[eleIdx1]")
                           .Define("eleSubTrkPhi_Lead",             "eleSubTrkPhi[eleIdx1]")
                           .Define("eleSubTrkCharge_Lead",          "eleSubTrkCharge[eleIdx1]")
                           .Define("eleSubTrkLayers_Lead",          "eleSubTrkLayers[eleIdx1]")
                           .Define("eleSubTrkMissHits_Lead",        "eleSubTrkMissHits[eleIdx1]")
                           .Define("eleSubTrkD0_Lead",              "eleSubTrkD0[eleIdx1]")
                           .Define("eleSubTrkDz_Lead",              "eleSubTrkDz[eleIdx1]")
                           .Define("eleGsfDeltaR_Lead",             "gsfDeltaR[eleIdx1]")
                           .Define("eleGsfPtRatio_Lead",            "gsfPtRatio[eleIdx1]")
                           .Define("eleGsfRelPtRatio_Lead",         "gsfRelPtRatio[eleIdx1]")
                           // photon variables 
                           .Define("phoE_Lead",                     "phoE[phoIdx1]")
                           .Define("phoSigmaE_Lead",                "phoSigmaE[phoIdx1]")
                           .Define("phoSCE_Lead",                   "phoSCE[phoIdx1]")
                           .Define("phoSCRawE_Lead",                "phoSCRawE[phoIdx1]")
                           .Define("phoESEnToRawE_Lead",            "phoESEnToRawE[phoIdx1]")
                           .Define("phoS4Full5x5_Lead",             "phoS4Full5x5[phoIdx1]")
                           .Define("phoSCPhi_Lead",                 "phoSCPhi[phoIdx1]")
                           .Define("phoSCEtaWidth_Lead",            "phoSCEtaWidth[phoIdx1]")
                           .Define("phoSCPhiWidth_Lead",            "phoSCPhiWidth[phoIdx1]")
                           .Define("phoSCBrem_Lead",                "phoSCBrem[phoIdx1]")
                           .Define("phohasPixelSeed_Lead",          "phohasPixelSeed[phoIdx1]")
                           .Define("phoEleVeto_Lead",               "phoEleVeto[phoIdx1]")
                           .Define("phoESEffSigmaRR_Lead",          "phoESEffSigmaRR[phoIdx1]")
                           .Define("phoSigmaIEtaIEtaFull5x5_Lead",  "phoSigmaIEtaIEtaFull5x5[phoIdx1]")
                           .Define("phoSigmaIEtaIPhiFull5x5_Lead",  "phoSigmaIEtaIPhiFull5x5[phoIdx1]")
                           .Define("phoSigmaIPhiIPhiFull5x5_Lead",  "phoSigmaIPhiIPhiFull5x5[phoIdx1]")
                           .Define("phoR9Full5x5_Lead",             "phoR9Full5x5[phoIdx1]")
                           .Define("phoPFChIso_Lead",               "phoPFChIso[phoIdx1]")
                           .Define("phoPFPhoIso_Lead",              "phoPFPhoIso[phoIdx1]")
                           .Define("phoPFNeuIso_Lead",              "phoPFNeuIso[phoIdx1]")
                           .Define("phoPFChWorstIso_Lead",          "phoPFChWorstIso[phoIdx1]")
                           .Define("phoTrkIsoHollowConeDR03_Lead",  "phoTrkIsoHollowConeDR03[phoIdx1]")
                           .Define("phoIDMVA_Lead",                 "phoIDMVA[phoIdx1]")
                           .Define("phoFiredSingleTrgs_Lead",       "phoFiredSingleTrgs[phoIdx1]")
                           .Define("phoFiredDoubleTrgs_Lead",       "phoFiredDoubleTrgs[phoIdx1]")
                           .Define("phoFiredL1Trgs_Lead",           "phoFiredL1Trgs[phoIdx1]")
                           // jet variables
                           .Define("jetSmearPt_Lead",               "if (isVbf) return jetSmearPt[jetIdx1]; else return (float) 0.;")
                           .Define("jetSmearEn_Lead",               "if (isVbf) return jetSmearEn[jetIdx1]; else return (float) 0.;")
                           .Define("jetEta_Lead",                   "if (isVbf) return jetEta[jetIdx1]; else return (float) 0.;")
                           .Define("jetPhi_Lead",                   "if (isVbf) return jetPhi[jetIdx1]; else return (float) 0.;")
                           .Define("jetMt_Lead",                    "if (isVbf) return jetMt[jetIdx1]; else return (float) 0.;")
                           .Define("jetArea_Lead",                  "if (isVbf) return jetArea[jetIdx1]; else return (float) 0.;")
                           .Define("jetLeadTrackPt_Lead",           "if (isVbf) return jetLeadTrackPt[jetIdx1]; else return (float) 0.;")
                           .Define("jetLeadTrackEta_Lead",          "if (isVbf) return jetLeadTrackEta[jetIdx1]; else return (float) 0.;")
                           .Define("jetLeadTrackPhi_Lead",          "if (isVbf) return jetLeadTrackPhi[jetIdx1]; else return (float) 0.;")
                           .Define("jetLepTrackPID_Lead",           "if (isVbf) return jetLepTrackPID[jetIdx1]; else return (int) 0;")
                           .Define("jetLepTrackPt_Lead",            "if (isVbf) return jetLepTrackPt[jetIdx1]; else return (float) 0.;")
                           .Define("jetLepTrackEta_Lead",           "if (isVbf) return jetLepTrackEta[jetIdx1]; else return (float) 0.;")
                           .Define("jetLepTrackPhi_Lead",           "if (isVbf) return jetLepTrackPhi[jetIdx1]; else return (float) 0.;")
                           .Define("jetCHF_Lead",                   "if (isVbf) return jetCHF[jetIdx1]; else return (float) 0.;") // chargedHadronEnergyFraction
                           .Define("jetNHF_Lead",                   "if (isVbf) return jetNHF[jetIdx1]; else return (float) 0.;") // neutralHadronEnergyFraction
                           .Define("jetCEF_Lead",                   "if (isVbf) return jetCEF[jetIdx1]; else return (float) 0.;") // chargedEmEnergyFraction
                           .Define("jetNEF_Lead",                   "if (isVbf) return jetNEF[jetIdx1]; else return (float) 0.;") // neutralEmEnergyFraction
                           .Define("jetNCH_Lead",                   "if (isVbf) return jetNCH[jetIdx1]; else return (int) 0;")    // chargedMultiplicity
                           .Define("jetNNP_Lead",                   "if (isVbf) return jetNNP[jetIdx1]; else return (int) 0;")    // neutralMultiplicity
                           .Define("jetSmearPt_subLead",            "if (isVbf) return jetSmearPt[jetIdx2]; else return (float) 0.;")
                           .Define("jetSmearEn_subLead",            "if (isVbf) return jetSmearEn[jetIdx2]; else return (float) 0.;")
                           .Define("jetEta_subLead",                "if (isVbf) return jetEta[jetIdx2]; else return (float) 0.;")
                           .Define("jetPhi_subLead",                "if (isVbf) return jetPhi[jetIdx2]; else return (float) 0.;")
                           .Define("jetMt_subLead",                 "if (isVbf) return jetMt[jetIdx2]; else return (float) 0.;")
                           .Define("jetArea_subLead",               "if (isVbf) return jetArea[jetIdx2]; else return (float) 0.;")
                           .Define("jetLeadTrackPt_subLead",        "if (isVbf) return jetLeadTrackPt[jetIdx2]; else return (float) 0.;")
                           .Define("jetLeadTrackEta_subLead",       "if (isVbf) return jetLeadTrackEta[jetIdx2]; else return (float) 0.;")
                           .Define("jetLeadTrackPhi_subLead",       "if (isVbf) return jetLeadTrackPhi[jetIdx2]; else return (float) 0.;")
                           .Define("jetLepTrackPID_subLead",        "if (isVbf) return jetLepTrackPID[jetIdx2]; else return (int) 0;")
                           .Define("jetLepTrackPt_subLead",         "if (isVbf) return jetLepTrackPt[jetIdx2]; else return (float) 0.;")
                           .Define("jetLepTrackEta_subLead",        "if (isVbf) return jetLepTrackEta[jetIdx2]; else return (float) 0.;")
                           .Define("jetLepTrackPhi_subLead",        "if (isVbf) return jetLepTrackPhi[jetIdx2]; else return (float) 0.;")
                           .Define("jetCHF_subLead",                "if (isVbf) return jetCHF[jetIdx2]; else return (float) 0.;") // chargedHadronEnergyFraction
                           .Define("jetNHF_subLead",                "if (isVbf) return jetNHF[jetIdx2]; else return (float) 0.;") // neutralHadronEnergyFraction
                           .Define("jetCEF_subLead",                "if (isVbf) return jetCEF[jetIdx2]; else return (float) 0.;") // chargedEmEnergyFraction
                           .Define("jetNEF_subLead",                "if (isVbf) return jetNEF[jetIdx2]; else return (float) 0.;") // neutralEmEnergyFraction
                           .Define("jetNCH_subLead",                "if (isVbf) return jetNCH[jetIdx2]; else return (int) 0;")    // chargedMultiplicity
                           .Define("jetNNP_subLead",                "if (isVbf) return jetNNP[jetIdx2]; else return (int) 0;");   // neutralMultiplicity

    //=======================================================//
    //                  extract the results                  //
    //=======================================================//
    auto report = df_Final.Report();
    std::vector<ROOT::RDF::RResultPtr<TH1D>> nominal_hists;
    for (size_t i = 0; i < categories.size(); i++){
        nominal_hists.emplace_back(df_Final.Filter(Form("category == %ld", i+1)).Histo1D({"Hmass", "", 60, 110, 170}, "CMS_higgs_mass", "weight"));
    }
    std::vector<ROOT::RDF::Experimental::RResultMap<TH1D>> hists;
    if (doVariation){
        for (size_t i = 0; i < nominal_hists.size(); i++){
            hists.emplace_back(ROOT::RDF::Experimental::VariationsFor(nominal_hists[i]));
        }
    }

    std::vector<std::string> defColNames = df_Final.GetDefinedColumnNames();
    std::vector<std::string> Cols = {"category", "category_str", "run", "event", "rho", "nVtx", "CMS_higgs_mass"};
    if (isSignalMC)
        Cols.emplace_back("isTrueM2");
    std::string p4Type = "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> >";
    for (size_t i = 0; i < defColNames.size(); i++){
        bool foundLead = ba::contains(defColNames[i], "_Lead");
        bool foundsubLead = ba::contains(defColNames[i], "_subLead");
        bool foundWeight = ba::contains(defColNames[i], "wei");
        bool foundP4 = df_Final.GetColumnType(defColNames[i]) == p4Type;
        if (foundLead || foundsubLead || foundWeight || foundP4)
            Cols.emplace_back(defColNames[i]);
    }
    std::sort(Cols.begin(), Cols.end());

    printf(" [+] Save miniTree in:\n");
    printf("     - %s\n", outpath.c_str());
    TString outDir = gSystem->GetDirName(outpath.c_str());
    if (!fs::exists(outDir.Data())){
        int status = gSystem->mkdir(outDir.Data(), true);
        if (status == -1)
            throw::std::runtime_error(Form("Fail to create directory: %s", outDir.Data()));
    }
    df_Final.Snapshot("miniTree", outpath.c_str(), Cols);    

    TString cutflow = utils::printReport(report);
    auto cutResult = std::make_shared<TObjString>(cutflow.Data());

    // save the report and histograms
    std::unique_ptr<TFile> fout(TFile::Open(outpath.c_str(), "UPDATE"));
    fout->cd();
    cutResult->Write("cutflow");
    if (doVariation){
        for(size_t i = 0; i < hists.size(); i++){
            for (auto it = hists[i].begin(); it != hists[i].end(); it++){
                TString var(it->first);
                if (var.Contains("pho:")) var.ReplaceAll("pho:", "");
                if (var.Contains("ele:")) var.ReplaceAll("ele:", "");
                if (var.Contains("jet:")) var.ReplaceAll("jet:", "");
                
                it->second->SetName(Form(("cat%ld_"+var).Data(), i+1));
                it->second->Write();
            }
        }
    }

    time_iset.Stop();
    auto htime_iset = utils::GetHumanTime(time_iset.RealTime());
    printf(" [+] Job's done:\n");
    printf("     - Total time: %d hours, %d mins, %d secs\n", htime_iset[0], htime_iset[1], htime_iset[2]);

    return 0;
}


int main(int argc, char** argv){
    TStopwatch time;
    time.Start();

    ArgumentParser(argc, argv);

    std::string range_text = (range == -1) ? "(-1 means all)" : "";
    std::string thread_txt = (range == -1) ? Form("%d", nthreads) : "1 (turn off MT if range is specified)";
    utils::printParameters(config_path, era, thread_txt, range);
    
    //=======================================================//
    //   Run xAna over all of the ntuples list in datasets   //
    //=======================================================//
    Info("Process the following ntuples sequentially\n");
    for (size_t i = 0; i < readpaths.size(); i++){
        Info(Form(" %ld. %s\n", i+1, readpaths[i].c_str()));
    }
    printf("\n");

    if (range == -1 && nthreads != 1)
        ROOT::EnableImplicitMT(nthreads);

    DeclareMVAFunc();
    for (size_t i = 0; i < readpaths.size(); i++){
        printf("*************************  %ld  ************************\n", i+1);
        xAna(readpaths[i], savepaths[i], i);
    }

    time.Stop();
    auto htime = utils::GetHumanTime(time.RealTime());
    printf("******************************************************\n");
    printf("All done: %d hours, %d mins, %d secs\n", htime[0], htime[1], htime[2]);
    printf("\n");

    return 0;
}