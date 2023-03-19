R__LOAD_LIBRARY($HDalitzEle_LOC/lib/libHDalitzEle.so)
R__ADD_INCLUDE_PATH($HDalitzEle_LOC/include)
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"

#include "ElectronSel.h"
#include "PhotonSel.h"

using namespace std;
void test_diphoHLT(){
    // ROOT::EnableImplicitMT(10);


    // auto df = ROOT::RDataFrame("ggNtuplizer/EventTree", "/data5/ggNtuples/V10_06_30_02/job_UL16_Dalitz_eeg_m125_preVFP/*.root");
    // auto nf = df.Define("isGooEle", )
    
    
    
    // auto nf = df.Define("isHggPho",             phoSel::HggPresel,          {"nPho", "rhoAll", "phoSCEta", "phoPFChIso", "phoPFPhoIso",    "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"})
    //             .Define("isHggEle",             eleSel::HggPresel,          {"nEle", "eleSCEta", "eleSCPhi", "nPho", "rhoAll", "phoSCEta", "phoSCPhi", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"})
    //             .Define("isDiPhoHLT", "((HLTPho >> 14) & 1) == 1")
    //             .Filter("Sum(isHggPho) > 1 && Sum(isHggEle) != 0", "presel")
    //             .Filter("isDiPhoHLT", "hlt");
    
    
    
    auto df = ROOT::RDataFrame("miniTree", "~/HDalitzEle/misc/GenStudy/miniTree/UL2016preVFP/miniTree_HDalitz_ggF_eeg_125_UL2016preVFP.root");
    auto nf = df
                // .Filter("category == 2 && diGenEle.Pt() > (125 * 0.3) && elePresel_Lead == 1 && phoPresel_Lead == 1", "presel")
                .Filter("category == 2 && diGenEle.Pt() > (125 * 0.3) && phoPresel_Lead == 1 && RecoPho_Lead.Pt() > 35", "presel")
                .Filter("isDiPhoHLTSeed_g == 1", "hlt");
                // .Filter("isDiPhoHLTUnseed_g && isDiPhoHLTSeed_gs", "hlt");

    // auto nf = df.Filter("((HLTPho >> 14) & 1) == 1", "dipho")
    //             .Define("isHggPho",             phoSel::HggPresel,          {"nPho", "rhoAll", "phoSCEta", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"})
    //             .Define("isHggEle",             eleSel::HggPresel,          {"nEle", "eleSCEta", "eleSCPhi", "nPho", "rhoAll", "phoSCEta", "phoSCPhi", "phoPFChIso", "phoPFPhoIso", "phoTrkIsoHollowConeDR03", "phoR9Full5x5", "phoEt", "phoSigmaIEtaIEtaFull5x5", "phoHoverE"});
                // .Filter("Sum(isHggPho) == 0 || Sum(isHggEle) == 0", "test");
    
    nf.Report()->Print();

    // std::string display = nf.Display({"event", "isHggEle", "eleSCEta", "eleSCPhi", "isHggPho", "phoSCEta", "phoSCPhi", "isDiPhoHLT"}, 10)->AsString();
    // cout << display << endl;
}