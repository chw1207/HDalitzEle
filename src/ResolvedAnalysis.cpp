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
#include "TMVASafeReader.h"

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
    {"PhoScaleGainUp",  "phoScale_syst_up"},
    {"PhoScaleStatDo",  "phoScale_stat_dn"},
    {"PhoScaleSystDo",  "phoScale_syst_dn"},
    {"PhoScaleGainDo",  "phoScale_gain_dn"},
    {"PhoSigmaPhiUp",   "phoResol_phi_up"},
    {"PhoSigmaRhoUp",   "phoResol_rho_up"},
    {"PhoSigmaRhoDo",   "phoResol_rho_dn"},

    // electron energy
    {"EleScaleStatUp",  "eleScale_stat_up"},
    {"EleScaleSystUp",  "eleScale_syst_up"},
    {"EleScaleGainUp",  "eleScale_syst_up"},
    {"EleScaleStatDo",  "eleScale_stat_dn"},
    {"EleScaleSystDo",  "eleScale_syst_dn"},
    {"EleScaleGainDo",  "eleScale_gain_dn"},
    {"EleSigmaPhiUp",   "eleResol_phi_up"},
    {"EleSigmaRhoUp",   "eleResol_rho_up"},
    {"EleSigmaRhoDo",   "eleResol_rho_dn"}
};
std::vector<std::string> categories = {
    // "Resolved-EBHR9", 
    // "Resolved-EBLR9", 
    // "Resolved-EE"
    "Resolved"
};
std::vector<std::string> scale_factors = { // required scale factors, if values are not provided by config, they will be assigned as 1.  
    "RecoEleGt20",
    "RecoEleLt20",
    "Fall17EleID",
    "SingleEleHLT",
    // "HggPreselForPho",
    "HggPhoCSEV",
    "Fall17PhoID"
};

int N_merged_cat = 5; //! number of merged categories 


void ArgumentParser(int argc, char** argv){
    po::options_description desc{"Options"};
    desc.add_options()
        ("help,h",                                      "Higgs Dalitz decay analysis script to select two elctrons and a photon")
        ("config,c",        po::value(&config_path),    "configuration file (required)")
        ("range,r",         po::value(&range),          "maximal number of events to process; -1 means all; turn off the MT")
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


// assuming that the only difference of the miniTree path for resolved and merged category is "resolved" -> "merged".
// eg. 
//      merged: /data4/chenghan/electron/miniTree_merged/UL2017/miniTree_HDalitz_ggF_eeg_125_UL2017.root
//      resolved: /data4/chenghan/electron/miniTree_resolved/UL2017/miniTree_HDalitz_ggF_eeg_125_UL2017.root
std::vector<Long64_t> GetMergedEvents(std::string repath){
    TString MergedPath(repath);
    MergedPath.ReplaceAll("resolved", "merged");
    if (MergedPath.Contains("SingleEle"))
        MergedPath.ReplaceAll("SingleEle", "DoubleEG");

    std::vector<Long64_t> v; 
    TString MergedDir = gSystem->GetDirName(MergedPath.Data());
    if (!fs::exists(MergedDir.Data())){
        Warn(Form(" [+] Fail to find the root files for merged category: %s\n", MergedDir.Data()));
        return v;
    }
    if (!fs::exists(MergedPath.Data()) && !MergedPath.EndsWith("*.root")){
        Warn(Form(" [+] Fail to find the root files for merged category: %s\n", MergedPath.Data()));
        return v;
    }
        
    auto df_merged = ROOT::RDataFrame("miniTree", MergedPath.Data(), {"event"});
    auto evCol = df_merged.Take<Long64_t>("event");

    printf(" [+] Find out the merged events in :\n");
    printf("     - %s\n", MergedPath.Data());
    printf("     - %ld events are found\n", (*evCol).size());
    return *evCol;
}


// analysis
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
        auto nf = df.Define("diLepMCMass",     gen::CalcLHEMee,   {"lhePID", "lhePx", "lhePy", "lhePz"})
                    .Define("HasIntPho",       gen::IsPhoIntConv, {"mcPID", "mcMomPID", "mcStatusFlag"})
                    .Filter("diLepMCMass < 60 && HasIntPho == 0", "mc preslections");

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
        auto af = nf.Filter("event % 2 == 1",  "analysis"); // == 1 for analysis, == 0 for ID training and regression
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
    // get the event numbers selected as merged category
    auto merged_ev = GetMergedEvents(outpath);
    std::sort(merged_ev.begin(), merged_ev.end());
    df_analysis = df_analysis.Filter([&](Long64_t ev){
        if (std::binary_search(merged_ev.begin(), merged_ev.end(), ev))
            return false;
        else 
            return true;
    }, {"event"}, "remove merged events");

    std::string SingleEleHLT;
    std::string DiEleHLT = "(HLTEleMuX >> 40  & 1) == 1";
    if (ba::contains(era, "2016")) SingleEleHLT = "(HLTEleMuX >> 4  & 1) == 1"; // HLT_Ele27_WPTight_Gsf_v
    if (ba::contains(era, "2017")) SingleEleHLT = "(HLTEleMuX >> 54 & 1) == 1"; // HLT_Ele32_WPTight_Gsf_L1DoubleEG_v
    if (ba::contains(era, "2018")) SingleEleHLT = "(HLTEleMuX >> 55 & 1) == 1"; // HLT_Ele32_WPTight_Gsf_v
    auto mf = df_analysis.Define("era",                 [&](){return era;})
                         .Filter(SingleEleHLT, "trigger cut")
                         .Filter("isPVGood == 1", "Good Vtx")
                         .Filter("nEle > 1 && nPho > 0", "2 eles, 1 pho")
                         .Define("phoS4Full5x5",        "phoE2x2Full5x5/phoE5x5Full5x5")
                         .Define("phoESEnToRawE",       "(phoESEnP1+phoESEnP2)/phoSCRawE"); 

    //=======================================================//
    //                  photon selections                    //
    //=======================================================//
    // official photon ID WP
    //https://github.com/cms-sw/cmssw/blob/809c232871e781eff61634dfb284c2611032fb43/RecoEgamma/PhotonIdentification/python/Identification/mvaPhotonID_Fall17_94X_V2_cff.py#L22
    // std::vector<std::unique_ptr<TMVASafeReader>> hggReader(2); // Hgg ID readers (only need for data or bkg mc) 
    // std::vector<std::string> hggvals_EB;
    // std::vector<std::string> hggvals_EE;
    if (!isSignalMC){
        // auto hggModelFiles = cfg["external_files"]["HggPhoID_model"].as<std::map<std::string, std::string>>();
        // hggReader[0] = std::make_unique<TMVASafeReader>(hggModelFiles["EB"], nthreads);
        // hggvals_EB = {
        //     "phoSCRawE",
        //     "phoR9Full5x5",
        //     "phoSigmaIEtaIEtaFull5x5",
        //     "phoSCEtaWidth",
        //     "phoSCPhiWidth",
        //     "phoSigmaIEtaIPhiFull5x5",
        //     "phoS4Full5x5",
        //     "phoPFPhoIso",
        //     "phoPFChIso",
        //     "phoPFChWorstIso",
        //     "phoSCEta",
        //     "phoRho"
        // };
        // std::string hggvals_EB_str = ba::join(hggvals_EB, ", ");

        // hggReader[1] = std::make_unique<TMVASafeReader>(hggModelFiles["EE"], nthreads);
        // hggvals_EE = {
        //     "phoSCRawE",
        //     "phoR9Full5x5",
        //     "phoSigmaIEtaIEtaFull5x5",
        //     "phoSCEtaWidth",
        //     "phoSCPhiWidth",
        //     "phoSigmaIEtaIPhiFull5x5",
        //     "phoS4Full5x5",
        //     "phoPFPhoIso",
        //     "phoPFChIso",
        //     "phoPFChWorstIso",
        //     "phoSCEta",
        //     "phoRho",
        //     "phoESEffSigmaRR",
        //     "phoESEnToRawE"
        // };
        // std::string hggvals_EE_str = ba::join(hggvals_EE, ", ");
        mf = mf.Define("phoRho",                "ROOT::RVec<float>(nPho, rho)")
               .Define("phoCorrIDMVA",          "phoIDMVA")   
               .Define("phoCorrR9Full5x5",      "phoR9Full5x5");  
            //    .Define("phoMVAEEVals",          Form("MakeMVAVals<%d, float>(%s)", (int)hggvals_EE.size(), hggvals_EE_str.c_str()))  
            //    .Define("phoMVAEBVals",          Form("MakeMVAVals<%d, float>(%s)", (int)hggvals_EB.size(), hggvals_EB_str.c_str()))
            //    .DefineSlot("phoCorrHggIDMVA",   [&](unsigned int slot,
            //                                         const ROOT::RVec<float>& phoSCEta,
            //                                         const ROOT::RVec<std::vector<float>>& phoMVAEBVals,
            //                                         const ROOT::RVec<std::vector<float>>& phoMVAEEVals){
            //                                             ROOT::RVec<float> mva(phoSCEta.size());
            //                                             for (size_t i = 0; i < phoSCEta.size(); i++){
            //                                                 mva[i] = (fabs(phoSCEta[i]) < 1.479) ? hggReader[0]->Compute(phoMVAEBVals[i], slot)[0]
            //                                                                                      : hggReader[1]->Compute(phoMVAEEVals[i], slot)[0];
            //                                             }
            //                                             return mva;
            //                                         }, {"phoSCEta", "phoMVAEBVals", "phoMVAEEVals"});
    }
    auto pf = mf.Define("isEBPho",                  "abs(phoSCEta) < 1.4442")
                .Define("isEEPho",                  "abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5")  
                .Define("isHggPho",                 phoSel::HggPresel,  {"nPho", "rhoAll", "phoSCEta", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"})
                // .Define("isGoodPho",                "(isEBPho || isEEPho) && phoEleVeto == 1 && isHggPho && phoCorrHggIDMVA > -0.9")
                .Define("isGoodPho",                "phoEleVeto == 1 && ((phoIDMVA > -0.02 && isEBPho) || (phoIDMVA > -0.26 && isEEPho))")
                .Filter("ROOT::VecOps::Sum(isGoodPho) > 0", "event with good pho")
                .Define("phoIdx1",                  [ ](const ROOT::RVec<int>& isGoodPho,
                                                        const ROOT::RVec<float>& phoCalibEt){
                                                            // get the photon with highest et in good photons 
                                                            return utils::getIdx(isGoodPho, phoCalibEt)[0];
                                                        }, {"isGoodPho", "phoCalibEt"})
                .Define("phoCalibEt_Lead",          "phoCalibEt[phoIdx1]")
                .Define("phoSCEta_Lead",            "phoSCEta[phoIdx1]")
                .Define("phoCorrR9Full5x5_Lead",    "phoCorrR9Full5x5[phoIdx1]")
                // .Define("phoCorrHggIDMVA_Lead",     "phoCorrHggIDMVA[phoIdx1]")
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
    pf = pf.Define("pho", "ROOT::Math::PtEtaPhiMVector(phoCalibEt_Lead, phoEta[phoIdx1], phoPhi[phoIdx1], 0.)");  

    //=======================================================//
    //                electron selections                    //
    //=======================================================//
    auto ef = pf.Define("eleESEnToRawE",        "(eleESEnP1 + eleESEnP2) / eleSCRawEn")
                .Define("isEBEle",              "abs(eleSCEta) < 1.4442")
                .Define("isEEEle",              "abs(eleSCEta) > 1.566 && abs(eleSCEta) < 2.5")
                .Define("isMVAEle",             eleSel::Fall17V2ID,         {"nEle", "elePt", "eleSCEta", "eleIDMVAIso"}) // WP90
                .Define("isGoodEle",            "((isEBEle && abs(eleD0) < 0.02 && abs(eleDz) < 0.1) && isMVAEle) ||"
                                                "((isEEEle && abs(eleD0) < 0.05 && abs(eleDz) < 0.2) && isMVAEle)")
                .Filter("ROOT::VecOps::Sum(isGoodEle) > 1", "event with good eles")
                .Define("eleIdx",               utils::getIdx,              {"isGoodEle", "eleCalibPt"})
                .Define("eleIdx1",              "eleIdx[0]")
                .Define("eleIdx2",              "eleIdx[1]")
                .Define("elePt_Lead",           "elePt[eleIdx1]")
                .Define("eleSCEta_Lead",        "eleSCEta[eleIdx1]")
                .Define("elePt_subLead",        "eleCalibPt[eleIdx2]")
                .Define("eleSCEta_subLead",     "eleSCEta[eleIdx2]")
                .Define("eleCalibPt_Lead",      "eleCalibPt[eleIdx1]")
                .Define("eleCalibPt_subLead",   "eleCalibPt[eleIdx2]");
    
    // systematic variations for electron pt
    std::vector<std::string> ele_sys;
    std::vector<std::string> ele_sys_bra1;
    std::vector<std::string> ele_sys_bra2;
    if (doVariation){
        for (size_t i = 0; i < variations.size(); i++){
            if (ba::contains(variations[i], "Ele")){
                ele_sys.emplace_back(variations[i]);
                ele_sys_bra1.emplace_back(allowed_variations.at(variations[i])+"[eleIdx1]/cosh(eleEta[eleIdx1])");
                ele_sys_bra2.emplace_back(allowed_variations.at(variations[i])+"[eleIdx2]/cosh(eleEta[eleIdx2])");
            }
        }
        if (ele_sys.size() > 0){
            printf(" [+] Book systematic variations for electron:\n");
            for (size_t i = 0; i < ele_sys.size(); i++){
                printf("     - %-21s: %s, %s\n", ele_sys[i].c_str(), ele_sys_bra1[i].c_str(), ele_sys_bra2[i].c_str());
            }

            std::string ele_all_var1 = ba::join(ele_sys_bra1, ", ");
            std::string ele_all_var2 = ba::join(ele_sys_bra2, ", ");
            std::string varied_str  = Form("ROOT::RVec<ROOT::RVec<float>> v = {{%s}, {%s}}; return v;", ele_all_var1.c_str(), ele_all_var2.c_str());
            ef = ef.Vary({"eleCalibPt_Lead", "eleCalibPt_subLead"}, varied_str, ele_sys, "ele");
        }
    }
    ef = ef.Define("ele1",  "ROOT::Math::PtEtaPhiMVector(eleCalibPt_Lead, eleEta[eleIdx1], elePhi[eleIdx1], 0.000511)")
           .Define("ele2",  "ROOT::Math::PtEtaPhiMVector(eleCalibPt_subLead, eleEta[eleIdx2], elePhi[eleIdx2], 0.000511)");

    //=======================================================//
    //             kinematic event selections                //
    //=======================================================//
    auto kf = ef.Define("diEle",            "ele1 + ele2")
                .Define("H",                "ele1 + ele2 + pho")
                .Filter("pho.Pt() > 15", "pho pt cut")
                .Filter("ele1.Pt() > 35 && ele2.Pt() > 7", "trigger threshold")
                .Filter("eleCharge[eleIdx1] * eleCharge[eleIdx2] < 0", "opposite charge")
                .Filter("diEle.M() < 50",  "Mee < 50 GeV")
                .Filter("diEle.M() > 11 || diEle.M() < 8", "reject Upsilon")
                .Filter("diEle.M() > 3.5 || diEle.M() < 2.5", "reject Jpsi")
                .Filter("(diEle.Pt()/H.M()) > 0.3 && (pho.Pt()/H.M()) > 0.3", "pt mass ratio cut")
                .Filter("ROOT::VecOps::DeltaR(ele1.Eta(), pho.Eta(), ele1.Phi(), pho.Phi()) > 1", "dR(e1,pho) > 1")
                .Filter("ROOT::VecOps::DeltaR(ele2.Eta(), pho.Eta(), ele2.Phi(), pho.Phi()) > 1", "dR(e2,pho) > 1")
                .Filter("H.M() > 110. && H.M() < 170.", "three body mass cut");
    //! End of the selections

    //=======================================================//
    //                      categorization                   //
    //=======================================================//
    auto cf = kf.Define("isEBHR9",              "abs(phoSCEta[phoIdx1]) < 1.4442 && phoCorrR9Full5x5_Lead > 0.96")
                .Define("isEBLR9",              "abs(phoSCEta[phoIdx1]) < 1.4442 && phoCorrR9Full5x5_Lead <= 0.96")
                .Define("isEE",                 "abs(phoSCEta[phoIdx1]) > 1.566  && abs(phoSCEta[phoIdx1]) < 2.5")
                .Define("category",             [&](const bool isEBHR9,
                                                    const bool isEBLR9,
                                                    const bool isEE){
                                                        // int cat = (N_merged_cat+1);
                                                        // if (isEBHR9)        return cat;
                                                        // else if (isEBLR9)   return cat + 1;
                                                        // else if (isEE)      return cat + 2;
                                                        // else throw std::runtime_error("No proper category for Re");
                                                        int cat = (N_merged_cat+1);
                                                        return cat;
                                                    }, {"isEBHR9", "isEBLR9", "isEE"})
                .Define("category_str",         [&](int category){return categories.at(category-(N_merged_cat+1));},      {"category"});

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
            if (it->second){
                if (ba::contains(it->first, "RecoEle")) // a little bit tricky
                    continue;

                if (ba::contains(it->first, "Ele")){ 
                    df_sfs = it->second->GetRDF(df_sfs, it->first+"SF_Lead", it->first+"SFErr_Lead", {"eleSCEta_Lead", "elePt_Lead"});
                    df_sfs = it->second->GetRDF(df_sfs, it->first+"SF_subLead", it->first+"SFErr_subLead", {"eleSCEta_subLead", "elePt_subLead"});
                    printf("     - %-21s: in bins of %-22s && %-22s\n", (it->first).c_str(), "eleSCEta_(sub)Lead", "elePt_(sub)Lead");

                    SFs.emplace_back(it->first+"SF_Lead");
                    SFs.emplace_back(it->first+"SF_subLead");
                }
                if (ba::contains(it->first, "Pho")){ // in bins of sceta and pt
                    if (ba::contains(it->first, "Hgg")){
                        df_sfs = it->second->GetRDF(df_sfs, it->first+"SF_Lead", it->first+"SFErr_Lead", {"phoSCEta_Lead", "phoCorrR9Full5x5_Lead"});
                        printf("     - %-21s: in bins of %-22s && %-22s\n", (it->first).c_str(), "phoSCEta_Lead", "phoCorrR9Full5x5_Lead");
                        SFs.emplace_back(it->first+"SF_Lead");
                    }
                    else{
                        df_sfs = it->second->GetRDF(df_sfs, it->first+"SF_Lead", it->first+"SFErr_Lead", {"phoSCEta_Lead", "phoCalibEt_Lead"});
                        printf("     - %-21s: in bins of %-22s && %-22s\n", (it->first).c_str(), "phoSCEta_Lead", "phoCalibEt_Lead");
                        SFs.emplace_back(it->first+"SF_Lead");
                    }
                }
            }
            else{
                df_sfs = df_sfs.Define(it->first+"SF_Lead",    "(float) 1.")
                               .Define(it->first+"SFErr_Lead", "(float) 0.");
                SFs.emplace_back(it->first+"SF_Lead");
            }
        }
        if (Weis["RecoEleGt20"] && Weis["RecoEleLt20"]){
            df_sfs = df_sfs.Define("RecoEleSF_Lead",            [&](const float elePt_Lead,
                                                                    const float eleSCEta_Lead){
                                                                    float SF = (elePt_Lead > 20) ? Weis["RecoEleGt20"]->GetSFFromEGM(eleSCEta_Lead, elePt_Lead)
                                                                                                 : Weis["RecoEleLt20"]->GetSFFromEGM(eleSCEta_Lead, elePt_Lead);
                                                                    return SF;
                                                                }, {"eleSCEta_Lead", "elePt_Lead"})
                           .Define("RecoEleSFErr_Lead",         [&](const float elePt_Lead,
                                                                    const float eleSCEta_Lead){
                                                                    float SFerr = (elePt_Lead > 20) ? Weis["RecoEleGt20"]->GetSFErrFromEGM(eleSCEta_Lead, elePt_Lead)
                                                                                                    : Weis["RecoEleLt20"]->GetSFErrFromEGM(eleSCEta_Lead, elePt_Lead);
                                                                    return SFerr;
                                                                }, {"eleSCEta_Lead", "elePt_Lead"})
                           .Define("RecoEleSF_subLead",         [&](const float elePt_subLead,
                                                                    const float eleSCEta_subLead){
                                                                    float SF = (elePt_subLead > 20) ? Weis["RecoEleGt20"]->GetSFFromEGM(eleSCEta_subLead, elePt_subLead)
                                                                                                    : Weis["RecoEleLt20"]->GetSFFromEGM(eleSCEta_subLead, elePt_subLead);
                                                                    return SF;
                                                                }, {"eleSCEta_subLead", "elePt_subLead"})
                           .Define("RecoEleSFErr_subLead",      [&](const float elePt_subLead,
                                                                    const float eleSCEta_subLead){
                                                                    float SFerr = (elePt_subLead > 20) ? Weis["RecoEleGt20"]->GetSFErrFromEGM(eleSCEta_subLead, elePt_subLead)
                                                                                                       : Weis["RecoEleLt20"]->GetSFErrFromEGM(eleSCEta_subLead, elePt_subLead);
                                                                    return SFerr;
                                                                }, {"eleSCEta_subLead", "elePt_subLead"});
            printf("     - %-21s: in bins of electron pt and sceta\n", "RecoEle");
            SFs.emplace_back("RecoEleSF_Lead");
            SFs.emplace_back("RecoEleSF_subLead");
        }
        else{
            df_sfs = df_sfs.Define("RecoEleSF_Lead",       "(float) 1.")
                           .Define("RecoEleSFErr_Lead",    "(float) 0.")
                           .Define("RecoEleSF_subLead",     "(float) 1.")
                           .Define("RecoEleSFErr_subLead",  "(float) 0.");
            SFs.emplace_back("RecoEleSF_Lead");
            SFs.emplace_back("RecoEleSF_subLead");
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
                if (it->first == "HLTUp" && ba::contains(SFs[i], "HLT"))
                    wei_str = Form(" * %s", SFsUp[i].c_str());
                
                if (it->first == "HLTDo" && ba::contains(SFs[i], "HLT"))
                    wei_str = Form(" * %s", SFsDo[i].c_str());
                
                if (it->first == "EleIDUp" && ba::contains(SFs[i], "Ele") && !ba::contains(SFs[i], "HLT"))
                    wei_str = Form(" * %s", SFsUp[i].c_str());
                
                if (it->first == "EleIDDo" && ba::contains(SFs[i], "Ele") && !ba::contains(SFs[i], "HLT"))
                    wei_str = Form(" * %s", SFsDo[i].c_str());
                
                if (it->first == "PhoIDUp" && (ba::contains(SFs[i], "Pho")))
                    wei_str = Form(" * %s", SFsUp[i].c_str());
                
                if (it->first == "PhoIDDo" && (ba::contains(SFs[i], "Pho")))
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
                       .Define("genIdx_reco2",              gen::MatchedGenEle,  {"ele2", "nMC", "mcEta", "mcPhi", "mcPID", "mcMomPID", "mcGMomPID", "mcStatusFlag"})
                       .Define("nReco1MatchedGen",          "genIdx_reco1.size()")  // number of generator electrons can be matched to selected reco electron
                       .Define("nReco2MatchedGen",          "genIdx_reco2.size()")  // number of generator electrons can be matched to selected reco electron
                       .Define("isTrueRe",                  "nReco1MatchedGen == 1 && nReco2MatchedGen == 1 && genIdx_reco1[0] != genIdx_reco2[0]")
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
                       .Define("mcPIDForEle_subLead",       "if (nReco2MatchedGen > 0) return mcPID[genIdx_reco1[1]]; else return (int) 0;")
                       .Define("mcPtForEle_subLead",        "if (nReco2MatchedGen > 0) return mcPt[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMassForEle_subLead",      "if (nReco2MatchedGen > 0) return mcMass[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcEtaForEle_subLead",       "if (nReco2MatchedGen > 0) return mcEta[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcPhiForEle_subLead",       "if (nReco2MatchedGen > 0) return mcPhi[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMomPIDForEle_subLead",    "if (nReco2MatchedGen > 0) return mcMomPID[genIdx_reco1[1]]; else return (int) 0;")
                       .Define("mcMomPtForEle_subLead",     "if (nReco2MatchedGen > 0) return mcMomPt[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMomMassForEle_subLead",   "if (nReco2MatchedGen > 0) return mcMomMass[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMomEtaForEle_subLead",    "if (nReco2MatchedGen > 0) return mcMomEta[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcMomPhiForEle_subLead",    "if (nReco2MatchedGen > 0) return mcMomPhi[genIdx_reco1[1]]; else return (float) 0;")
                       .Define("mcGMomPIDForEle_subLead",   "if (nReco2MatchedGen > 0) return mcGMomPID[genIdx_reco1[1]]; else return (int) 0;")
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
    //=======================================================//
    //             define the final branches                 //
    //=======================================================//
    auto df_Final  = df_sfs.Define("CMS_higgs_mass",                    "(float) H.M()")
                           // electron variables 
                           .Define("eleCharge_Lead",                    "eleCharge[eleIdx1]")
                           .Define("eleEn_Lead",                        "eleEn[eleIdx1]")
                           .Define("eleSCEn_Lead",                      "eleSCEn[eleIdx1]")
                           .Define("eleEcalEn_Lead",                    "eleEcalEn[eleIdx1]")
                           .Define("eleESEnToRawE_Lead",                "eleESEnToRawE[eleIdx1]")
                           .Define("eleD0_Lead",                        "eleD0[eleIdx1]")
                           .Define("eleDz_Lead",                        "eleDz[eleIdx1]")
                           .Define("eleSIP_Lead",                       "eleSIP[eleIdx1]")
                           .Define("elePtError_Lead",                   "elePtError[eleIdx1]")
                           .Define("eleEta_Lead",                       "eleEta[eleIdx1]")
                           .Define("elePhi_Lead",                       "elePhi[eleIdx1]")
                           .Define("eleSCPhi_Lead",                     "eleSCPhi[eleIdx1]")
                           .Define("eleSCRawEn_Lead",                   "eleSCRawEn[eleIdx1]")
                           .Define("eleSCEtaWidth_Lead",                "eleSCEtaWidth[eleIdx1]")
                           .Define("eleSCPhiWidth_Lead",                "eleSCPhiWidth[eleIdx1]")
                           .Define("eleHoverE_Lead",                    "eleHoverE[eleIdx1]")
                           .Define("eleEoverP_Lead",                    "eleEoverP[eleIdx1]")
                           .Define("eleEoverPout_Lead",                 "eleEoverPout[eleIdx1]")
                           .Define("eleEoverPInv_Lead",                 "eleEoverPInv[eleIdx1]")
                           .Define("eleBrem_Lead",                      "eleBrem[eleIdx1]")
                           .Define("eledEtaAtVtx_Lead",                 "eledEtaAtVtx[eleIdx1]")
                           .Define("eledPhiAtVtx_Lead",                 "eledPhiAtVtx[eleIdx1]")
                           .Define("eleSigmaIEtaIEtaFull5x5_Lead",      "eleSigmaIEtaIEtaFull5x5[eleIdx1]")
                           .Define("eleSigmaIPhiIPhiFull5x5_Lead",      "eleSigmaIPhiIPhiFull5x5[eleIdx1]")
                           .Define("eleConvVeto_Lead",                  "eleConvVeto[eleIdx1]")
                           .Define("eleMissHits_Lead",                  "eleMissHits[eleIdx1]")
                           .Define("eleESEffSigmaRR_Lead",              "eleESEffSigmaRR[eleIdx1]")
                           .Define("elePFChIso_Lead",                   "elePFChIso[eleIdx1]")
                           .Define("elePFPhoIso_Lead",                  "elePFPhoIso[eleIdx1]")
                           .Define("elePFNeuIso_Lead",                  "elePFNeuIso[eleIdx1]")
                           .Define("elePFPUIso_Lead",                   "elePFPUIso[eleIdx1]")
                           .Define("eleIDMVAIso_Lead",                  "eleIDMVAIso[eleIdx1]")
                           .Define("eleIDMVANoIso_Lead",                "eleIDMVANoIso[eleIdx1]")
                           .Define("eleR9Full5x5_Lead",                 "eleR9Full5x5[eleIdx1]")
                           .Define("eleGSFChi2_Lead",                   "eleGSFChi2[eleIdx1]")
                           .Define("eleFiredSingleTrgs_Lead",           "eleFiredSingleTrgs[eleIdx1]")
                           .Define("eleFiredDoubleTrgs_Lead",           "eleFiredDoubleTrgs[eleIdx1]")
                           .Define("eleFiredL1Trgs_Lead",               "eleFiredL1Trgs[eleIdx1]")
                           .Define("eleCharge_subLead",                 "eleCharge[eleIdx2]")
                           .Define("eleEn_subLead",                     "eleEn[eleIdx2]")
                           .Define("eleSCEn_subLead",                   "eleSCEn[eleIdx2]")
                           .Define("eleEcalEn_subLead",                 "eleEcalEn[eleIdx2]")
                           .Define("eleESEnToRawE_subLead",             "eleESEnToRawE[eleIdx2]")
                           .Define("eleD0_subLead",                     "eleD0[eleIdx2]")
                           .Define("eleDz_subLead",                     "eleDz[eleIdx2]")
                           .Define("eleSIP_subLead",                    "eleSIP[eleIdx2]")
                           .Define("elePtError_subLead",                "elePtError[eleIdx2]")
                           .Define("eleEta_subLead",                    "eleEta[eleIdx2]")
                           .Define("elePhi_subLead",                    "elePhi[eleIdx2]")
                           .Define("eleSCPhi_subLead",                  "eleSCPhi[eleIdx2]")
                           .Define("eleSCRawEn_subLead",                "eleSCRawEn[eleIdx2]")
                           .Define("eleSCEtaWidth_subLead",             "eleSCEtaWidth[eleIdx2]")
                           .Define("eleSCPhiWidth_subLead",             "eleSCPhiWidth[eleIdx2]")
                           .Define("eleHoverE_subLead",                 "eleHoverE[eleIdx2]")
                           .Define("eleEoverP_subLead",                 "eleEoverP[eleIdx2]")
                           .Define("eleEoverPout_subLead",              "eleEoverPout[eleIdx2]")
                           .Define("eleEoverPInv_subLead",              "eleEoverPInv[eleIdx2]")
                           .Define("eleBrem_subLead",                   "eleBrem[eleIdx2]")
                           .Define("eledEtaAtVtx_subLead",              "eledEtaAtVtx[eleIdx2]")
                           .Define("eledPhiAtVtx_subLead",              "eledPhiAtVtx[eleIdx2]")
                           .Define("eleSigmaIEtaIEtaFull5x5_subLead",   "eleSigmaIEtaIEtaFull5x5[eleIdx2]")
                           .Define("eleSigmaIPhiIPhiFull5x5_subLead",   "eleSigmaIPhiIPhiFull5x5[eleIdx2]")
                           .Define("eleConvVeto_subLead",               "eleConvVeto[eleIdx2]")
                           .Define("eleMissHits_subLead",               "eleMissHits[eleIdx2]")
                           .Define("eleESEffSigmaRR_subLead",           "eleESEffSigmaRR[eleIdx2]")
                           .Define("elePFChIso_subLead",                "elePFChIso[eleIdx2]")
                           .Define("elePFPhoIso_subLead",               "elePFPhoIso[eleIdx2]")
                           .Define("elePFNeuIso_subLead",               "elePFNeuIso[eleIdx2]")
                           .Define("elePFPUIso_subLead",                "elePFPUIso[eleIdx2]")
                           .Define("eleIDMVAIso_subLead",               "eleIDMVAIso[eleIdx2]")
                           .Define("eleIDMVANoIso_subLead",             "eleIDMVANoIso[eleIdx2]")
                           .Define("eleR9Full5x5_subLead",              "eleR9Full5x5[eleIdx2]")
                           .Define("eleGSFChi2_subLead",                "eleGSFChi2[eleIdx2]")
                           .Define("eleFiredSingleTrgs_subLead",        "eleFiredSingleTrgs[eleIdx2]")
                           .Define("eleFiredDoubleTrgs_subLead",        "eleFiredDoubleTrgs[eleIdx2]")
                           .Define("eleFiredL1Trgs_subLead",            "eleFiredL1Trgs[eleIdx2]")
                           // photon variables 
                           .Define("phoE_Lead",                         "phoE[phoIdx1]")
                           .Define("phoSigmaE_Lead",                    "phoSigmaE[phoIdx1]")
                           .Define("phoSCE_Lead",                       "phoSCE[phoIdx1]")
                           .Define("phoSCRawE_Lead",                    "phoSCRawE[phoIdx1]")
                           .Define("phoESEnToRawE_Lead",                "phoESEnToRawE[phoIdx1]")
                           .Define("phoS4Full5x5_Lead",                 "phoS4Full5x5[phoIdx1]")
                           .Define("phoSCPhi_Lead",                     "phoSCPhi[phoIdx1]")
                           .Define("phoSCEtaWidth_Lead",                "phoSCEtaWidth[phoIdx1]")
                           .Define("phoSCPhiWidth_Lead",                "phoSCPhiWidth[phoIdx1]")
                           .Define("phoSCBrem_Lead",                    "phoSCBrem[phoIdx1]")
                           .Define("phohasPixelSeed_Lead",              "phohasPixelSeed[phoIdx1]")
                           .Define("phoEleVeto_Lead",                   "phoEleVeto[phoIdx1]")
                           .Define("phoESEffSigmaRR_Lead",              "phoESEffSigmaRR[phoIdx1]")
                           .Define("phoSigmaIEtaIEtaFull5x5_Lead",      "phoSigmaIEtaIEtaFull5x5[phoIdx1]")
                           .Define("phoSigmaIEtaIPhiFull5x5_Lead",      "phoSigmaIEtaIPhiFull5x5[phoIdx1]")
                           .Define("phoSigmaIPhiIPhiFull5x5_Lead",      "phoSigmaIPhiIPhiFull5x5[phoIdx1]")
                           .Define("phoR9Full5x5_Lead",                 "phoR9Full5x5[phoIdx1]")
                           .Define("phoPFChIso_Lead",                   "phoPFChIso[phoIdx1]")
                           .Define("phoPFPhoIso_Lead",                  "phoPFPhoIso[phoIdx1]")
                           .Define("phoPFNeuIso_Lead",                  "phoPFNeuIso[phoIdx1]")
                           .Define("phoPFChWorstIso_Lead",              "phoPFChWorstIso[phoIdx1]")
                           .Define("phoTrkIsoHollowConeDR03_Lead",      "phoTrkIsoHollowConeDR03[phoIdx1]")
                           .Define("phoIDMVA_Lead",                     "phoIDMVA[phoIdx1]")
                           .Define("phoFiredSingleTrgs_Lead",           "phoFiredSingleTrgs[phoIdx1]")
                           .Define("phoFiredDoubleTrgs_Lead",           "phoFiredDoubleTrgs[phoIdx1]")
                           .Define("phoFiredL1Trgs_Lead",               "phoFiredL1Trgs[phoIdx1]");
                           
    //=======================================================//
    //                  extract the results                  //
    //=======================================================//
    auto report = df_Final.Report();
    std::vector<ROOT::RDF::RResultPtr<TH1D>> nominal_hists;
    for (size_t i = 0; i < categories.size(); i++){
        nominal_hists.emplace_back(df_Final.Filter(Form("category == %ld", i+(N_merged_cat+1))).Histo1D({"Hmass", "", 60, 110, 170}, "CMS_higgs_mass", "weight"));
    }
    std::vector<ROOT::RDF::Experimental::RResultMap<TH1D>> hists;
    if (doVariation){
        for (size_t i = 0; i < nominal_hists.size(); i++){
            hists.emplace_back(ROOT::RDF::Experimental::VariationsFor(nominal_hists[i]));
        }
    }

    std::vector<std::string> defColNames = df_Final.GetDefinedColumnNames();
    std::vector<std::string> Cols = {"category", "category_str", "run", "event", "rho", "nVtx", "lumis", "CMS_higgs_mass"};
    if (isSignalMC)
        Cols.emplace_back("isTrueRe");
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
                
                it->second->SetName(Form(("cat%ld_"+var).Data(), i+(N_merged_cat+1)));
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