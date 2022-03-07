#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "./plugins/MuSelections.h"
#include "./plugins/PhoSelections.h"
#include "./plugins/puweicalc.h"


using namespace std;
using namespace ROOT::VecOps;


namespace Helper{
    // Function to get the index vector sorted by pT
    // * Reference: https://root.cern/doc/master/vo006__IndexManipulation_8C.html
    ROOT::RVec<int> getIdx(ROOT::RVec<int>& isgood, ROOT::RVec<float>& pt){
        ROOT::RVec<int> idx_select = Nonzero(isgood);
        ROOT::RVec<int> idx_sort = Reverse(Argsort(pt));
        ROOT::RVec<int> idx = Intersect(idx_sort, idx_select);

        return idx;
    }


    ROOT::RVec<TLorentzVector> P4Vector(ROOT::RVec<float>& pt, ROOT::RVec<float>& eta, ROOT::RVec<float>& phi, float m){
        ROOT::RVec<TLorentzVector> vec;
        vec.clear();

        for (int i = 0; i < pt.size(); ++i) {
            TLorentzVector p4;
            p4.SetPtEtaPhiM(pt.at(i), eta.at(i), phi.at(i), m);
            vec.push_back(p4);
        }

        return vec;
    }
}


template <typename T>
auto FindGoodMus(T &df) {
    auto nf = df.Define("isMediumID",   "PassMuonID(muIDbit, 2)")
                .Define("isGoodMuon",   "abs(muEta) < 2.4 && muPt > 10. && isMediumID")

                .Filter("Sum(isGoodMuon) > 1", "good muon")

                .Define("mu1Idx",       "Helper::getIdx(isGoodMuon, muPt)[0]")
                .Define("mu2Idx",       "Helper::getIdx(isGoodMuon, muPt)[1]")
                .Define("mu1",          "TLorentzVector v; v.SetPtEtaPhiM(muPt[mu1Idx], muEta[mu1Idx], muPhi[mu1Idx], 105.658*0.001); return v;")
                .Define("mu2",          "TLorentzVector v; v.SetPtEtaPhiM(muPt[mu2Idx], muEta[mu2Idx], muPhi[mu2Idx], 105.658*0.001); return v;")

                .Filter("muPt[mu1Idx] > 20 && muPt[mu2Idx] > 10", "HLT pTCut")
                .Filter("(muCharge[mu1Idx] * muCharge[mu2Idx]) < 0", "+/- charge");
    return nf;
}


template <typename T>
auto FindGoodPho(T &df) {
    auto nf = df.Define("phoP4",        "Helper::P4Vector(phoCalibEt, phoEta, phoPhi, 0.)")
                .Define("isEBPho",      "abs(phoSCEta) < 1.4442")
                .Define("isEEPho",      "abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5")
                .Define("isFSR",        "FSRSelec(mu1, mu2, phoP4)")
                .Define("isGoodPho",    "(isEBPho || isEEPho) && isFSR && phoCalibEt > 15.")
                .Define("isHggPho",     "HggPreSelection(rhoAll, nPho, phoSCEta, phoPFChIso, phoPFPhoIso, phoTrkIsoHollowConeDR03, phoR9Full5x5, phoCalibEt, phoSigmaIEtaIEtaFull5x5, phoHoverE)")
            
                .Filter("Sum(isGoodPho) > 0", "good phon")

                .Define("phoIdx",       "GetZPho(mu1, mu2, phoP4, isGoodPho)")
                .Define("pho",          "phoP4[phoIdx]")
                .Define("convIdx",      "ConvMatch(phoSCEta[phoIdx], phoSCPhi[phoIdx], phoSCE[phoIdx], nConv, convNTrks, convVtxX, convVtxY, convVtxZ, convFitPairPX, convFitPairPY, convFitPairPZ, convFitProb)")
                
                .Filter("convIdx != -1 && convVtxRadius[convIdx] < 16", "rconv < 16");
    return nf;
}


template <typename T>
auto FindFSRMu(T &df) {
    auto nf = df.Define("fsrmu",        "if (mu1.DeltaR(pho) > mu2.DeltaR(pho)) return mu2; else return mu1")
                .Define("nfsrmu",       "if (mu1.DeltaR(pho) > mu2.DeltaR(pho)) return mu1; else return mu2")
            
                .Filter("nfsrmu.Pt() > 20", "non-fsr mu");
    return nf;
}


template <typename T>
auto FindMatchEle(T &df) {
    auto nf = df.Define("eleIdx",       "GetMatchEle(phoSCEta[phoIdx], eleSCEta)")
                
                .Filter("eleIdx != -1", "match ele")
                .Filter("nGsfMatchToReco[eleIdx] > 0", "match gsf")
                
                .Define("ele",          "TLorentzVector v; v.SetPtEtaPhiM(eleCalibPt[eleIdx], eleEta[eleIdx], elePhi[eleIdx], 0.511*0.001); return v;");
    return nf;
}


template <typename T>
auto DefineFinalVars(T &df){
    auto nf = df.Define("Z",                "mu1 + mu2 + pho")
                .Define("dimu",             "mu1 + mu2")
                .Define("phoEleVeto1",      "phoEleVeto[phoIdx]")
                .Define("phoSCEta1",        "phoSCEta[phoIdx]")
                .Define("phoSCPhi1",        "phoSCPhi[phoIdx]")
                .Define("phoSCE1",          "phoSCE[phoIdx]")
                .Define("isHggPho1",        "isHggPho[phoIdx]")
                .Define("phoCalibE1",       "phoCalibE[phoIdx]")
                .Define("isEBPho1",         "isEBPho[phoIdx]")
                .Define("isEEPho1",         "isEEPho[phoIdx]")

                .Define("eleConvVeto1",     "eleConvVeto[eleIdx]")
                .Define("eleClass1",        "eleClass[eleIdx]")
                .Define("eleSCEta1",        "eleSCEta[eleIdx]")
                .Define("eleSCPhi1",        "eleSCPhi[eleIdx]")
                .Define("eleSCEn1",         "eleSCEn[eleIdx]")
                .Define("eleCalibEn1",      "eleCalibEn[eleIdx]")

                .Define("convVtxRadius1",   "convVtxRadius[convIdx]")
                .Define("convFitPairP1",    "sqrt(pow(convFitPairPX[convIdx], 2) + pow(convFitPairPY[convIdx], 2) + pow(convFitPairPZ[convIdx], 2))")
                .Define("convD01",          "convD0[convIdx]")
                .Define("convDz1",          "convDz[convIdx]")
                .Define("convL01",          "convL0[convIdx]")
                .Define("convLz1",          "convLz[convIdx]");        
    return nf;
}


template <typename T>
auto AddWeights(T &df, string era, int year, bool isMC){
    if (!isMC) return df;

    // set up the puwei calculator 
    PUWeightCalculator* puCalc[3] = {NULL, NULL, NULL};
    for (int i = 0; i < 3; i++){
        puCalc[i] = new PUWeightCalculator();
    }
    puCalc[0]->Init(PUfile(year, "nominal").c_str()); 
    puCalc[1]->Init(PUfile(year, "up").c_str());
    puCalc[2]->Init(PUfile(year, "down").c_str());

    auto get_pu = [&, puCalc](int run, ROOT::RVec<float>& puTrue){
        float puwei = puCalc[0]->GetWeight(run, puTrue[1]);
        return puwei;
    };
    auto get_pu_up = [&, puCalc](int run, ROOT::RVec<float>& puTrue){
        float puwei_up = puCalc[1]->GetWeight(run, puTrue[1]);
        return puwei_up;
    };
    auto get_pu_do = [&, puCalc](int run, ROOT::RVec<float>& puTrue){
        float puwei_down = puCalc[2]->GetWeight(run, puTrue[1]);
        return puwei_down;
    };

    auto nf = df.Define("puwei",            get_pu,             {"run", "puTrue"})
                .Define("puwei_up",         get_pu_up,          {"run", "puTrue"})
                .Define("puwei_down",       get_pu_do,          {"run", "puTrue"})
                .Define("genwei",           "if (genWeight > 0) return 1.; else return -1.;");

    return nf;
}


void rdfZmumug(string infile = "testSkimZg.root", string outfile = "test.root", int year = 2017, string era = "2017", bool isMC = true){
    TStopwatch time;
    time.Start();

    cout << "[INFO] Read_File(): " << infile.c_str() << endl;
    cout << "[INFO] Save_File(): " << outfile.c_str() << endl;

    ROOT::EnableImplicitMT();

    //=======================================================//
    // Load TTree into RDataFrame and then do the selections //
    //=======================================================//
    ROOT::RDataFrame df("ggNtuplizer/EventTree", infile.c_str());

    cout << "[INFO] Total number of events: " << df.Count().GetValue() << endl;
    cout << "[INFO] Pool size: " << df.GetNSlots() << endl;

    auto df1 = df.Define("isMC", [&, isMC]{return isMC;});
    auto df2 = FindGoodMus(df1);
    auto df3 = FindGoodPho(df2);
    auto df4 = FindFSRMu(df3);
    auto df5 = FindMatchEle(df4);
    auto df6 = DefineFinalVars(df5);
    auto df7 = AddWeights(df6, era.c_str(), year, isMC);
    auto dfFin = df7;

    //========================================================//
    // Visualize the selection results (cut flow, event count)//
    //========================================================//
    cout << "[INFO] Cut flow:" << std::endl;
    auto report = dfFin.Report();
    report->Print();

    //====================================================//
    // Define the final variables to save to the miniTree //
    //====================================================//
    vector<string> Vars = {
        "mu1", "mu2", "fsrmu", "nfsrmu", "Z", "dimu", "pho", "ele",
        "phoEleVeto1", "phoSCEta1", "phoSCPhi1", "phoSCE1", "isHggPho1", "phoCalibE1", "isEBPho1",
        "eleConvVeto1", "eleClass1", "eleSCEta1", "eleSCPhi1", "eleSCEn1", "eleCalibEn1",
        "convVtxRadius1", "convFitPairP1", "convD01", "convDz1", "convL01", "convLz1"
    };
    if (isMC){
        // concatenate Vars and weiVars
        vector<string> weiVars = {
            "puwei", "puwei_up", "puwei_down",  
            "genwei","mcwei"
        };
        Vars.insert(Vars.end(), weiVars.begin(), weiVars.end());
    }

    // save the minitree
    dfFin.Snapshot("outTree", outfile.c_str(), Vars);

    cout << "[INFO] Time taken: " << endl;
    time.Stop();
    time.Print();
}