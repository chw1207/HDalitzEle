#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "./plugins/roccor/RoccoR.cc"
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
auto FindGoodMus(T &df, string era) {
    map <string, string> rcName;
    rcName["2016_preVFP"] = "./plugins/roccor/RoccoR2016aUL.txt";
    rcName["2016_postVFP"] = "./plugins/roccor/RoccoR2016bUL.txt";
    rcName["2017"] = "./plugins/roccor/RoccoR2017UL.txt";
    rcName["2018"] = "./plugins/roccor/RoccoR2018UL.txt";

    RoccoR* rc = new RoccoR(rcName[era.c_str()].c_str());
    auto getRoccoR = [&, rc](bool isMC, long long event, int nMu, ROOT::RVec<int>& muCharge, ROOT::RVec<float>& muPt, ROOT::RVec<float>& muEta, ROOT::RVec<float>& muPhi, ROOT::RVec<int>& muTrkLayers){
        ROOT::RVec<float> sf;
        sf.clear();

        TRandom* rd = new TRandom3(event); // event as seed (just for reproducibility)
        float u[nMu];
        rd->RndmArray(nMu, u);
        for (int i = 0; i < nMu; i++){
            if (isMC)
                sf.push_back(rc->kSmearMC(muCharge[i], muPt[i], muEta[i], muPhi[i], muTrkLayers[i], u[i], 0, 0));
            else
                sf.push_back(rc->kScaleDT(muCharge[i], muPt[i], muEta[i], muPhi[i], 0, 0));
        }

        return sf;
    };

    auto nf = df.Define("muSF",                 getRoccoR,      {"isMC", "event", "nMu", "muCharge", "muPt", "muEta", "muPhi", "muTrkLayers"})
                .Define("muCorrPt",             "muPt * muSF")
                .Define("isMediumID",           "PassMuonID(muIDbit, 2)")
                .Define("isGoodMuon",           "abs(muEta) < 2.4 && muCorrPt > 10. && isMediumID")

                .Filter("Sum(isGoodMuon) > 1", "good muon")

                .Define("mu1Idx",               "Helper::getIdx(isGoodMuon, muCorrPt)[0]")
                .Define("mu2Idx",               "Helper::getIdx(isGoodMuon, muCorrPt)[1]")
                .Define("mu1",                  "TLorentzVector v; v.SetPtEtaPhiM(muCorrPt[mu1Idx], muEta[mu1Idx], muPhi[mu1Idx], 105.658*0.001); return v;")
                .Define("mu2",                  "TLorentzVector v; v.SetPtEtaPhiM(muCorrPt[mu2Idx], muEta[mu2Idx], muPhi[mu2Idx], 105.658*0.001); return v;")

                .Filter("muCorrPt[mu1Idx] > 20 && muCorrPt[mu2Idx] > 10", "HLT pTCut")
                .Filter("(muCharge[mu1Idx] * muCharge[mu2Idx]) < 0", "+/- charge");
    return nf;
}


template <typename T>
auto FindGoodPho(T &df) {
    auto nf = df.Define("phoP4",                "Helper::P4Vector(phoCalibEt, phoEta, phoPhi, 0.)")
                .Define("isEBPho",              "abs(phoSCEta) < 1.4442")
                .Define("isEEPho",              "abs(phoSCEta) > 1.566 && abs(phoSCEta) < 2.5")
                .Define("isFSR",                "FSRSelec(mu1, mu2, phoP4)")
                .Define("isGoodPho",            "(isEBPho || isEEPho) && isFSR && phoCalibEt > 15.")
                .Define("isHggPho",             "HggPreSelection(rhoAll, nPho, phoSCEta, phoPFChIso, phoPFPhoIso, phoTrkIsoHollowConeDR03, phoR9Full5x5, phoCalibEt, phoSigmaIEtaIEtaFull5x5, phoHoverE)")

                .Filter("Sum(isGoodPho) > 0", "good phon")

                .Define("phoIdx",               "GetZPho(mu1, mu2, phoP4, isGoodPho)")
                .Define("pho",                  "phoP4[phoIdx]")
                .Define("convIdx",              "ConvMatch(phoSCEta[phoIdx], phoSCPhi[phoIdx], phoSCE[phoIdx], nConv, convNTrks, convVtxX, convVtxY, convVtxZ, convFitPairPX, convFitPairPY, convFitPairPZ, convFitProb)")

                .Filter("convIdx != -1", "match conv");
    return nf;
}


template <typename T>
auto FindFSRMu(T &df) {
    auto nf = df.Define("fsrmu",                "if (mu1.DeltaR(pho) > mu2.DeltaR(pho)) return mu2; else return mu1")
                .Define("nfsrmu",               "if (mu1.DeltaR(pho) > mu2.DeltaR(pho)) return mu1; else return mu2")

                .Filter("nfsrmu.Pt() > 20", "non-fsr mu");
    return nf;
}


template <typename T>
auto FindMatchEle(T &df) {
    auto nf = df.Define("eleIdx",               "GetMatchEle(phoSCEta[phoIdx], eleSCEta)")

                .Filter("eleIdx != -1", "match ele")

                .Define("ele",                  "TLorentzVector v; v.SetPtEtaPhiM(eleCalibPt[eleIdx], eleEta[eleIdx], elePhi[eleIdx], 0.511*0.001); return v;");
    return nf;
}


template <typename T>
auto DefineFinalVars(T &df){
    auto nf = df.Define("Z",                    "mu1 + mu2 + pho")
                .Define("dimu",                 "mu1 + mu2")
                .Define("zMass",                "Z.M()")
                .Define("phoEleVeto_lep1",      "phoEleVeto[phoIdx]")
                .Define("phoSCEta_lep1",        "phoSCEta[phoIdx]")
                .Define("phoSCPhi_lep1",        "phoSCPhi[phoIdx]")
                .Define("phoSCE_lep1",          "phoSCE[phoIdx]")
                .Define("isHggPho_lep1",        "isHggPho[phoIdx]")
                .Define("phoCalibE_lep1",       "phoCalibE[phoIdx]")
                .Define("phoCalibEt_lep1",      "phoCalibEt[phoIdx]")
                .Define("isEBPho_lep1",         "isEBPho[phoIdx]")
                .Define("isEEPho_lep1",         "isEEPho[phoIdx]")

                .Define("phoSCEtaWidth_lep1",   "phoSCEtaWidth[phoIdx]")
                .Define("phoSCPhiWidth_lep1",   "phoSCPhiWidth[phoIdx]")
                .Define("phoSCBrem_lep1",       "phoSCBrem[phoIdx]")
                .Define("phoHoverE_lep1",       "phoHoverE[phoIdx]")
                .Define("phoR9Full5x5_lep1",    "phoR9Full5x5[phoIdx]")

                .Define("nGsfMatchToReco_lep1", "nGsfMatchToReco[eleIdx]")
                .Define("eleConvVeto_lep1",     "eleConvVeto[eleIdx]")
                .Define("eleClass_lep1",        "eleClass[eleIdx]")
                .Define("eleXGBID_lep1",        "eleXGBID[eleIdx]")
                .Define("eleSCEta_lep1",        "eleSCEta[eleIdx]")
                .Define("eleSCPhi_lep1",        "eleSCPhi[eleIdx]")
                .Define("eleSCEn_lep1",         "eleSCEn[eleIdx]")
                .Define("eleCalibEn_lep1",      "eleCalibEn[eleIdx]")

                .Define("convVtxRadius_lep1",   "convVtxRadius[convIdx]")
                .Define("convFitPairP_lep1",    "sqrt(pow(convFitPairPX[convIdx], 2) + pow(convFitPairPY[convIdx], 2) + pow(convFitPairPZ[convIdx], 2))")
                .Define("convD0_lep1",          "convD0[convIdx]")
                .Define("convDz_lep1",          "convDz[convIdx]")
                .Define("convL0_lep1",          "convL0[convIdx]")
                .Define("convLz_lep1",          "convLz[convIdx]")

                .Define("eleCalibPt_lep1",      "eleCalibPt[eleIdx]")
                .Define("eleSCRawEn_lep1",      "eleSCRawEn[eleIdx]")

                .Define("eledEtaAtVtx_lep1",    "eledEtaAtVtx[eleIdx]")
                .Define("eledPhiAtVtx_lep1",    "eledPhiAtVtx[eleIdx]")
                .Define("elePtError_lep1",      "elePtError[eleIdx]")
                .Define("eleHoverE_lep1",       "eleHoverE[eleIdx]")
                .Define("eleEoverP_lep1",       "eleEoverP[eleIdx]")
                .Define("eleEoverPout_lep1",    "eleEoverPout[eleIdx]")
                .Define("eleEoverPInv_lep1",    "eleEoverPInv[eleIdx]")

                .Define("eleSCEtaWidth_lep1",   "eleSCEtaWidth[eleIdx]")
                .Define("eleSCPhiWidth_lep1",   "eleSCPhiWidth[eleIdx]")
                .Define("eleSigmaIEtaIEtaFull5x5_lep1",      "eleSigmaIEtaIEtaFull5x5[eleIdx]")
                .Define("eleSigmaIPhiIPhiFull5x5_lep1",      "eleSigmaIPhiIPhiFull5x5[eleIdx]")
                .Define("eleR9Full5x5_lep1",    "eleR9Full5x5[eleIdx]")
                .Define("eleBrem_lep1",         "eleBrem[eleIdx]")

                .Define("elePFChIso_lep1",      "elePFChIso[eleIdx]")
                .Define("elePFPhoIso_lep1",     "elePFPhoIso[eleIdx]")
                .Define("elePFNeuIso_lep1",     "elePFNeuIso[eleIdx]")

                .Define("gsfPtRatio_lep1",      "gsfPtRatio[eleIdx]")
                .Define("gsfDeltaR_lep1",       "gsfDeltaR[eleIdx]")
                .Define("gsfRelPtRatio_lep1",   "gsfRelPtRatio[eleIdx]")
                .Define("gsfPtSum_lep1",        "gsfPtSum[eleIdx]");
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
                .Define("genwei",           "if (genWeight > 0) return 1.; else return -1.;")
                .Define("wei1",             "puwei * mcwei")
                .Define("wei2",             "puwei * mcwei * genwei");

    return nf;
}


void rdfZmumug(string infile, string outfile, int year, string era, bool isMC){
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
    auto df2 = FindGoodMus(df1, era.c_str());
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
    vector<string> defColNames = dfFin.GetDefinedColumnNames();
    vector<string> Vars = {"event", "zMass"};
    for (int i = 0; i < defColNames.size(); i++){
        size_t foundLep1 = defColNames[i].find("_lep1");
        size_t foundLep2 = defColNames[i].find("_lep2");

        if (foundLep1 != string::npos || foundLep2 != string::npos)
            Vars.push_back(defColNames[i]);
        if (dfFin.GetColumnType(defColNames[i]) == "TLorentzVector")
            Vars.push_back(defColNames[i]);
    }

    if (isMC){
        // concatenate Vars and weiVars
        vector<string> weiVars = {
            "puwei", "puwei_up", "puwei_down",
            "genwei","mcwei", "wei1", "wei2"
        };
        Vars.insert(Vars.end(), weiVars.begin(), weiVars.end());
    }

    // save the minitree
    dfFin.Snapshot("outTree", outfile.c_str(), Vars);

    cout << "[INFO] Time taken: " << endl;
    time.Stop();
    time.Print();

    cout << endl;
}