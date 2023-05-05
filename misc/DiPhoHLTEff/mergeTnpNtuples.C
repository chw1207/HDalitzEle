#include <iostream>
#include <vector>
#include <string>
#include "TROOT.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TStopwatch.h"
#include "TSystem.h"
#include "XGBReader.h"
#include "boost/algorithm/string/join.hpp"
using namespace std;

// This script is mainly used to merge the tnp ntuples.
// It also adds merged electron ID prediction.

void mergeTnpNtuples(string infiles="/data4/chenghan/tnpTuples-unseedDiPhoHLT/UL2017_DY_LO/*.root", string outfile="/data4/chenghan/tnpTuples-unseedDiPhoHLT/UL2017_merged/DY_LO.root"){
    cout << "------- start -------" << endl;
    TStopwatch time;
    time.Start();

    XGBReader rM2EB("/data4/chenghan/external/MergedID/Output_Merged2GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EB/XGB/XGB_modelXGB.txt");
    XGBReader rM2EE("/data4/chenghan/external/MergedID/Output_Merged2GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EE/XGB/XGB_modelXGB.txt");

    vector<string> vals = {
        "event_rho",
        "el_sc_eta",
        "el_sc_rawE",
        "el_dEtaIn",
        "el_dPhiIn",
        "el_pterr",
        "el_gghoe",
        "el_ep",
        "el_eelepout",
        "el_EoverPInv",
        "el_etaW",
        "el_phiW",
        "el_5x5_sieie",
        "el_5x5_sipip",
        "el_5x5_r9",
        "el_fbrem",
        "el_chIso",
        "el_phoIso",
        "el_neuIso",
        "el_tksPtRatio",
        "el_tksdr",
        "el_tksRelPtRatio"
    };
    string vals_str = boost::algorithm::join(vals, ", ");
    
    cout << "[INFO] Read_Files(): " << infiles << endl;
    ROOT::EnableImplicitMT(10);
    auto df = ROOT::RDataFrame("tnpEleTrig/fitter_tree", infiles);

    auto nf = df.Define("el_mergedVal",         Form("vector<float> v = {%s}; return v;", vals_str.c_str()))
                .Define("el_mergedMVA",         [&](const float el_sc_eta,
                                                    const vector<float>& el_mergedVal){
                                                        if (fabs(el_sc_eta) < 1.479)
                                                            return rM2EB.Compute(el_mergedVal);
                                                        else
                                                            return rM2EE.Compute(el_mergedVal);
                                                    }, {"el_sc_eta", "el_mergedVal"})
                .Define("MERGEDMVA_WP",         "(fabs(el_sc_eta) < 1.479) ? (float) 0.4263 : (float) 0.4205")
                .Define("el_mergedMVAC0",       "el_mergedMVA[0]")
                .Define("el_mergedMVAC1",       "el_mergedMVA[1]")
                .Define("el_mergedMVAC2",       "el_mergedMVA[2]")
                .Define("el_maxMergedCl",       "ROOT::RVec<float> v(el_mergedMVA.data(), el_mergedMVA.size()); return (int) ROOT::VecOps::ArgMax(v);")
                .Define("el_passMergedMVA",     "el_ntks > 1 && el_maxMergedCl == 0 && (el_mergedMVAC0 > MERGEDMVA_WP)");
    
    auto colNames = nf.GetColumnNames();
    vector<string> outcolNames;
    for (int i = 0; i < colNames.size(); i++){
        if ((colNames[i] == "el_mergedMVA") || (colNames[i] == "el_mergedVal"))
            continue;
        outcolNames.push_back(colNames[i]);
    }
    
    cout << "[INFO] Save_File(): " << outfile << endl;
    nf.Snapshot("tnpEleTrig/fitter_tree", outfile, outcolNames);
    ROOT::DisableImplicitMT();

    cout << "[INFO] Time taken: " << endl;
    time.Stop();
    time.Print();

    cout << "------- done  -------" << endl;
    cout << endl;
}
