#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TStopwatch.h"
#include "TLorentzVector.h"
#include "./plugins/help.h"
#include "./plugins/puweicalc.h"
#include "../../pluginsV2/skim_utilities.h" // add gsf information to electron
#include "../../pluginsV2/EleSelections.h" // HggPreSelection


template <typename T>
auto RecotoGSF(T &df){
    auto nf = df.Filter("nEle > 0",                         "nEle > 0")
                .Define("RecoEle",                          "P4Vector(eleCalibPt, eleEta, elePhi, 0.000511)")
                .Define("isMainGSF",                        "IsMainGSF(eleD0, gsfD0)")
                .Define("ambiguousGSF",                     "TrackElectronMatching(eleD0, gsfD0, isMainGSF)")
                .Define("nGsfMatchToReco",                  "CalcNGsfMatchToReco(nEle, ambiguousGSF)")
                .Define("eleTrkIdx",                        "FindMainGSF(nEle, ambiguousGSF)")
                .Define("eleSubTrkIdx",                     "FindSubGSF_dRMin(nEle, ambiguousGSF, gsfEta, gsfPhi)")
                .Define("eleSubPtTrkIdx",                   "FindSubGSF_PtMax(nEle, ambiguousGSF, gsfPt)")

                // main gsf track variables
                .Define("eleTrkPt",                         "MatchIdex(eleTrkIdx, gsfPt)")
                .Define("eleTrkEta",                        "MatchIdex(eleTrkIdx, gsfEta)")
                .Define("eleTrkPhi",                        "MatchIdex(eleTrkIdx, gsfPhi)")
                .Define("eleTrkCharge",                     "MatchIdex(eleTrkIdx, gsfCharge)")
                .Define("eleTrkLayers",                     "MatchIdex(eleTrkIdx, gsfLayers)")
                .Define("eleTrkMissHits",                   "MatchIdex(eleTrkIdx, gsfMissHits)")
                .Define("eleTrkD0",                         "MatchIdex(eleTrkIdx, gsfD0)")
                .Define("eleTrkDz",                         "MatchIdex(eleTrkIdx, gsfDz)")

                // second gsf track variables(min dR track to main gsf track among the associated gsf tracks except for main gsf track)
                .Define("eleSubTrkPt",                      "MatchIdex(eleSubTrkIdx, gsfPt)")
                .Define("eleSubTrkEta",                     "MatchIdex(eleSubTrkIdx, gsfEta)")
                .Define("eleSubTrkPhi",                     "MatchIdex(eleSubTrkIdx, gsfPhi)")
                .Define("eleSubTrkCharge",                  "MatchIdex(eleSubTrkIdx, gsfCharge)")
                .Define("eleSubTrkLayers",                  "MatchIdex(eleSubTrkIdx, gsfLayers)")
                .Define("eleSubTrkMissHits",                "MatchIdex(eleSubTrkIdx, gsfMissHits)")
                .Define("eleSubTrkD0",                      "MatchIdex(eleSubTrkIdx, gsfD0)")
                .Define("eleSubTrkDz",                      "MatchIdex(eleSubTrkIdx, gsfDz)")

                // second gsf track variables(max pT track among the associated gsf tracks except for main gsf track)
                .Define("eleSubPtTrkPt",                    "MatchIdex(eleSubPtTrkIdx, gsfPt)")
                .Define("eleSubPtTrkEta",                   "MatchIdex(eleSubPtTrkIdx, gsfEta)")
                .Define("eleSubPtTrkPhi",                   "MatchIdex(eleSubPtTrkIdx, gsfPhi)")
                .Define("eleSubPtTrkCharge",                "MatchIdex(eleSubPtTrkIdx, gsfCharge)")
                .Define("eleSubPtTrkLayers",                "MatchIdex(eleSubPtTrkIdx, gsfLayers)")
                .Define("eleSubPtTrkMissHits",              "MatchIdex(eleSubPtTrkIdx, gsfMissHits)")
                .Define("eleSubPtTrkD0",                    "MatchIdex(eleSubPtTrkIdx, gsfD0)")
                .Define("eleSubPtTrkDz",                    "MatchIdex(eleSubPtTrkIdx, gsfDz)");

    return nf;
}


template <typename T>
auto RecotoGen(T &df){
    auto nf = df.Define("GenIdx",                           "Helper::GenIdxVec(RecoEle, mcPt, mcEta, mcPhi, mcMass, mcPID, mcMomPID)")
                .Define("isGoodEle",                        "GenIdx != -1")

                .Filter("Sum(isGoodEle) > 0", "1 good ele")

                .Define("recoEleIdx1",                      "Helper::getIdx(isGoodEle, eleCalibPt)[0]")
                .Define("RecoEle_lep1",                     "TLorentzVector v; v.SetPtEtaPhiM(eleCalibPt[recoEleIdx1], eleEta[recoEleIdx1], elePhi[recoEleIdx1], 0.000511); return v;")
                .Define("genEleIdx1",                       "GenIdx[recoEleIdx1]")
                .Define("GenEle_lep1",                      "TLorentzVector v; v.SetPtEtaPhiM(mcPt[genEleIdx1], mcEta[genEleIdx1], mcPhi[genEleIdx1], 0.000511); return v;")

                // conversion vtx matching
                .Define("convVtxIdx1",                      "Helper::RecoConvMatch(eleSCEta[recoEleIdx1], eleSCPhi[recoEleIdx1], eleSCEn[recoEleIdx1], nConv, convNTrks, convVtxX, convVtxY, convVtxZ, convFitPairPX, convFitPairPY, convFitPairPZ, convFitProb)")

                // Hgg preselection
                .Define("elePresel",                        "HggPreSelection(rhoAll, nEle, eleSCEta, nPho, phoSCEta, phoPFChIso, phoPFPhoIso, phoTrkIsoHollowConeDR03, phoR9Full5x5, phoCalibEt, phoSigmaIEtaIEtaFull5x5, phoHoverE)");

    return nf;
}


template <typename T>
auto DefineFinalVars(T &df){
                // leading gen variables
    auto nf = df.Define("mcPt_lep1",                        "mcPt[genEleIdx1]")
                .Define("mcEta_lep1",                       "mcEta[genEleIdx1]")
                .Define("mcPhi_lep1",                       "mcPhi[genEleIdx1]")
                .Define("mcVtx_lep1",                       "mcVtx[genEleIdx1]")
                .Define("mcVty_lep1",                       "mcVty[genEleIdx1]")
                .Define("mcVtz_lep1",                       "mcVtz[genEleIdx1]")

                // leading reco variables
                .Define("elePresel_lep1",                   "elePresel[recoEleIdx1]")
                .Define("eleCharge_lep1",                   "eleCharge[recoEleIdx1]")
                .Define("eleChargeConsistent_lep1",         "eleChargeConsistent[recoEleIdx1]")
                .Define("eleConvVeto_lep1",                 "eleConvVeto[recoEleIdx1]")
                .Define("eleEcalDrivenSeed_lep1",           "eleEcalDrivenSeed[recoEleIdx1]")
                .Define("eleEn_lep1",                       "eleEn[recoEleIdx1]")
                .Define("eleSCEn_lep1",                     "eleSCEn[recoEleIdx1]")
                .Define("eleEcalEn_lep1",                   "eleEcalEn[recoEleIdx1]")
                .Define("eleESEnP1_lep1",                   "eleESEnP1[recoEleIdx1]")
                .Define("eleESEnP2_lep1",                   "eleESEnP2[recoEleIdx1]")
                .Define("eleD0_lep1",                       "eleD0[recoEleIdx1]")
                .Define("eleDz_lep1",                       "eleDz[recoEleIdx1]")
                .Define("eleSIP_lep1",                      "eleSIP[recoEleIdx1]")
                .Define("elePt_lep1",                       "elePt[recoEleIdx1]")
                .Define("elePtError_lep1",                  "elePtError[recoEleIdx1]")
                .Define("eleEta_lep1",                      "eleEta[recoEleIdx1]")
                .Define("elePhi_lep1",                      "elePhi[recoEleIdx1]")
                .Define("eleR9_lep1",                       "eleR9[recoEleIdx1]")
                .Define("eleR9Full5x5_lep1",                "eleR9Full5x5[recoEleIdx1]")
                .Define("eleCalibPt_lep1",                  "eleCalibPt[recoEleIdx1]")
                .Define("eleDiffCalibOriPt_lep1",           "(eleCalibPt[recoEleIdx1] - elePt[recoEleIdx1])/elePt[recoEleIdx1]")
                .Define("eleCalibEn_lep1",                  "eleCalibEn[recoEleIdx1]")
                .Define("eleSCEta_lep1",                    "eleSCEta[recoEleIdx1]")
                .Define("eleSCPhi_lep1",                    "eleSCPhi[recoEleIdx1]")
                .Define("eleSCRawEn_lep1",                  "eleSCRawEn[recoEleIdx1]")
                .Define("eleSCEtaWidth_lep1",               "eleSCEtaWidth[recoEleIdx1]")
                .Define("eleSCPhiWidth_lep1",               "eleSCPhiWidth[recoEleIdx1]")
                .Define("eleHoverE_lep1",                   "eleHoverE[recoEleIdx1]")
                .Define("eleEoverP_lep1",                   "eleEoverP[recoEleIdx1]")
                .Define("eleEoverPout_lep1",                "eleEoverPout[recoEleIdx1]")
                .Define("eleEoverPInv_lep1",                "eleEoverPInv[recoEleIdx1]")
                .Define("eleBrem_lep1",                     "eleBrem[recoEleIdx1]")
                .Define("eledEtaAtVtx_lep1",                "eledEtaAtVtx[recoEleIdx1]")
                .Define("eledPhiAtVtx_lep1",                "eledPhiAtVtx[recoEleIdx1]")
                .Define("eleSigmaIEtaIEtaFull5x5_lep1",     "eleSigmaIEtaIEtaFull5x5[recoEleIdx1]")
                .Define("eleSigmaIPhiIPhiFull5x5_lep1",     "eleSigmaIPhiIPhiFull5x5[recoEleIdx1]")
                .Define("eleESEffSigmaRR_lep1",             "eleESEffSigmaRR[recoEleIdx1]")
                .Define("elePFChIso_lep1",                  "elePFChIso[recoEleIdx1]")
                .Define("elePFPhoIso_lep1",                 "elePFPhoIso[recoEleIdx1]")
                .Define("elePFNeuIso_lep1",                 "elePFNeuIso[recoEleIdx1]")
                .Define("elePFPUIso_lep1",                  "elePFPUIso[recoEleIdx1]")
                .Define("elePFClusEcalIso_lep1",            "elePFClusEcalIso[recoEleIdx1]")
                .Define("elePFClusHcalIso_lep1",            "elePFClusHcalIso[recoEleIdx1]")
                .Define("eleIDMVAIso_lep1",                 "eleIDMVAIso[recoEleIdx1]")
                .Define("eleIDMVANoIso_lep1",               "eleIDMVANoIso[recoEleIdx1]")
                .Define("eleTrkdxy_lep1",                   "eleTrkdxy[recoEleIdx1]")
                .Define("eleKFHits_lep1",                   "eleKFHits[recoEleIdx1]")
                .Define("eleKFChi2_lep1",                   "eleKFChi2[recoEleIdx1]")
                .Define("eleGSFChi2_lep1",                  "eleGSFChi2[recoEleIdx1]")
                .Define("eleESEnToRawE_lep1",               "(eleESEnP1[recoEleIdx1] + eleESEnP2[recoEleIdx1])/eleSCRawEn[recoEleIdx1]")
                .Define("nGsfMatchToReco_lep1",             "nGsfMatchToReco[recoEleIdx1]")
                .Define("eleTrk_lep1",                      "TLorentzVector v; v.SetPtEtaPhiM(eleTrkPt[recoEleIdx1], eleTrkEta[recoEleIdx1], eleTrkPhi[recoEleIdx1], 0.000511); return v;")
                .Define("eleTrkPt_lep1",                    "eleTrkPt[recoEleIdx1]")
                .Define("eleTrkEta_lep1",                   "eleTrkEta[recoEleIdx1]")
                .Define("eleTrkPhi_lep1",                   "eleTrkPhi[recoEleIdx1]")
                .Define("eleTrkCharge_lep1",                "eleTrkCharge[recoEleIdx1]")
                .Define("eleTrkLayers_lep1",                "eleTrkLayers[recoEleIdx1]")
                .Define("eleTrkMissHits_lep1",              "eleTrkMissHits[recoEleIdx1]")
                .Define("eleTrkD0_lep1",                    "eleTrkD0[recoEleIdx1]")
                .Define("eleTrkDz_lep1",                    "eleTrkDz[recoEleIdx1]")
                .Define("eleSubTrk_lep1",                   "TLorentzVector v; v.SetPtEtaPhiM(eleSubTrkPt[recoEleIdx1], eleSubTrkEta[recoEleIdx1], eleSubTrkPhi[recoEleIdx1], 0.000511); return v;")
                .Define("eleSubTrkPt_lep1",                 "eleSubTrkPt[recoEleIdx1]")
                .Define("eleSubTrkEta_lep1",                "eleSubTrkEta[recoEleIdx1]")
                .Define("eleSubTrkPhi_lep1",                "eleSubTrkPhi[recoEleIdx1]")
                .Define("eleSubTrkCharge_lep1",             "eleSubTrkCharge[recoEleIdx1]")
                .Define("eleSubTrkLayers_lep1",             "eleSubTrkLayers[recoEleIdx1]")
                .Define("eleSubTrkMissHits_lep1",           "eleSubTrkMissHits[recoEleIdx1]")
                .Define("eleSubTrkD0_lep1",                 "eleSubTrkD0[recoEleIdx1]")
                .Define("eleSubTrkDz_lep1",                 "eleSubTrkDz[recoEleIdx1]")
                .Define("eleSubPtTrk_lep1",                 "TLorentzVector v; v.SetPtEtaPhiM(eleSubPtTrkPt[recoEleIdx1], eleSubPtTrkEta[recoEleIdx1], eleSubPtTrkPhi[recoEleIdx1], 0.000511); return v;")
                .Define("eleSubPtTrkPt_lep1",               "eleSubPtTrkPt[recoEleIdx1]")
                .Define("eleSubPtTrkEta_lep1",              "eleSubPtTrkEta[recoEleIdx1]")
                .Define("eleSubPtTrkPhi_lep1",              "eleSubPtTrkPhi[recoEleIdx1]")
                .Define("eleSubPtTrkCharge_lep1",           "eleSubPtTrkCharge[recoEleIdx1]")
                .Define("eleSubPtTrkLayers_lep1",           "eleSubPtTrkLayers[recoEleIdx1]")
                .Define("eleSubPtTrkMissHits_lep1",         "eleSubPtTrkMissHits[recoEleIdx1]")
                .Define("eleSubPtTrkD0_lep1",               "eleSubPtTrkD0[recoEleIdx1]")
                .Define("eleSubPtTrkDz_lep1",               "eleSubPtTrkDz[recoEleIdx1]")
                .Define("gsfPtSum_lep1",                    "eleTrk_lep1.Pt() + eleSubTrk_lep1.Pt()")
                .Define("gsfPtRatio_lep1",                  "eleSubTrk_lep1.Pt() / eleTrk_lep1.Pt()")
                .Define("gsfDeltaR_lep1",                   "eleTrk_lep1.DeltaR(eleSubTrk_lep1)")
                .Define("gsfMissHitsSum_lep1",              "eleTrkMissHits[recoEleIdx1] + eleSubTrkMissHits[recoEleIdx1]")
                .Define("gsfRelPtRatio_lep1",               "if (nGsfMatchToReco_lep1 > 1) return (eleTrk_lep1 + eleSubTrk_lep1).Pt() / eleCalibPt_lep1; else return eleTrk_lep1.Pt() / eleCalibPt_lep1;")
                .Define("convMatched_lep1",                 "if (nConv > 0 && convVtxIdx1 != -1) return 1; else return 0;")
                .Define("convVtxRadius_lep1",               "if (nConv > 0 && convVtxIdx1 != -1) return convVtxRadius[convVtxIdx1]; else return (float) -999;")
                .Define("convD0_lep1",                      "if (nConv > 0 && convVtxIdx1 != -1) return convD0[convVtxIdx1]; else return (float) -999;")
                .Define("convDz_lep1",                      "if (nConv > 0 && convVtxIdx1 != -1) return convDz[convVtxIdx1]; else return (float) -999;")
                .Define("convL0_lep1",                      "if (nConv > 0 && convVtxIdx1 != -1) return convL0[convVtxIdx1]; else return (float) -999;")
                .Define("convLz_lep1",                      "if (nConv > 0 && convVtxIdx1 != -1) return convLz[convVtxIdx1]; else return (float) -999;")

                .Define("diTrk",                            "eleTrk_lep1 + eleSubTrk_lep1")
                .Define("diTrkPtMax",                       "eleTrk_lep1 + eleSubPtTrk_lep1");

    return nf;
}


void rdfGENQCD(vector<string> infile, string outfile, int year, string era, string proc){
    TStopwatch time;
    time.Start();

    cout << "[INFO] " << infile.size() << " files are found!" << endl;
    cout << "[INFO] Save_File(): " << outfile.c_str() << endl;

    //=======================================================//
    // Load TTree into RDataFrame                            //
    //=======================================================//
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("ggNtuplizer/EventTree", infile);
    cout << "[INFO] Process: " << proc << " in " << era << endl;

    //========================================================//
    // Calculate the seral weight (add into the RDataFrame)   //
    //========================================================//
    // calculate mc weight
    auto pos = df.Filter("genWeight > 0",   "positive event").Count();
    auto neg = df.Filter("genWeight <= 0",  "negative event").Count();
    const int totalev = pos.GetValue() - neg.GetValue();

    // XS
    map<string, float> xsMap; // fb
    xsMap["HT50to100"] = 187300000 * 1000.;
    xsMap["HT100to200"] = 23590000 * 1000.;
    xsMap["HT200to300"] = 1555000 * 1000.;
    xsMap["HT300to500"] = 324500 * 1000.;
    xsMap["HT500to700"] = 30310 * 1000.;
    xsMap["HT700to1000"] = 6444 * 1000.;
    xsMap["HT1000to1500"] = 1127 * 1000.;
    xsMap["HT1500to2000"] = 109.8 * 1000.;
    xsMap["HT2000toInf"] = 21.98 * 1000.;

    // luminosity
    map<int, float> lumiMap;
    lumiMap[2016] = 35.917;
    lumiMap[2017] = 41.525;
    lumiMap[2018] = 59.725;

    const float procXS = xsMap[proc.c_str()];
    const float mcwei = ((procXS * lumiMap[year]) / totalev);
    cout << scientific;
    cout << "[INFO] Number of events with genwei = " << totalev << endl;
    cout << "[INFO] MC weight = " << mcwei << endl;

    // set up the puwei calculator
    PUWeightCalculator* puCalc = new PUWeightCalculator();;
    puCalc->Init(PUfile(year, "nominal").c_str());

    auto get_pu = [&, puCalc](int run, ROOT::RVec<float>& puTrue){
        const float puwei = puCalc->GetWeight(run, puTrue[1]);
        return puwei;
    };

    float instwei = procXS/xsMap["HT50to100"];

    auto wf = df.Define("puwei",    get_pu,     {"run", "puTrue"})
                .Define("mcwei",    [&, mcwei]{return mcwei;})
                .Define("procXS",   [&, procXS]{return procXS;})
                .Define("instwei",  [&, instwei]{return instwei;})
                .Define("genwei",   "if (genWeight > 0) return (float) 1.; else return (float) -1.;")
                .Define("wei",      "mcwei * genwei * puwei");


    //=======================================================//
    // Do the selections on the RDataFrame                   //
    //=======================================================//
    auto df1 = RecotoGSF(wf);
    auto df2 = RecotoGen(df1);
    auto df3 = DefineFinalVars(df2);
    auto dfFin = df3.Cache();

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
    vector<string> Vars;
    for (int i = 0; i < defColNames.size(); i++){
        size_t foundLep1 = defColNames[i].find("_lep1");

        if (foundLep1 != string::npos)
            Vars.push_back(defColNames[i]);
    }
    vector<string> extralVars = {
        "mcwei", "puwei", "genwei", "procXS", "instwei", "wei", "rho", "rhoAll", "nVtx", "nGoodVtx", "isPVGood",
        "diTrk", "diTrkPtMax"
    };
    Vars.insert(Vars.end(), extralVars.begin(), extralVars.end());

    dfFin.Snapshot("outTree", outfile.c_str(), Vars);

    cout << "[INFO] Time taken: " << endl;
    time.Stop();
    time.Print();

    cout << endl;
}