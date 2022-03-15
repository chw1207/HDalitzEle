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

using namespace std;
using namespace ROOT::VecOps;

// fb
map<string, float> XS_HDalitz(){
    map<string, float> XSmap;
    XSmap["ggF_125GeV"] = 48.58 * 1000 * 8.10E-5;
    XSmap["VBF_125GeV"] = 3.782 * 1000 * 8.10E-5;
    XSmap["WH_125GeV"] = 1.373 * 1000 * 8.10E-5;
    XSmap["ZH_125GeV"] = 0.8839 * 1000 * 8.10E-5;

    XSmap["ggF_120GeV"] = 52.22 * 1000 * 7.88E-5;
    XSmap["VBF_120GeV"] = 3.935 * 1000 * 7.88E-5;
    XSmap["WH_120GeV"] = 1.565 * 1000 * 7.88E-5;
    XSmap["ZH_120GeV"] = 0.9939 * 1000 * 7.88E-5;

    XSmap["ggF_130GeV"] = 45.31 * 1000 * 8.02E-5;
    XSmap["VBF_130GeV"] = 3.637 * 1000 * 8.02E-5;
    XSmap["WH_130GeV"] = 1.209 * 1000 * 8.02E-5;
    XSmap["ZH_130GeV"] = 0.4539 * 1000 * 8.02E-5;
    
    return XSmap;
}


float luminosity(int year){
    float lumi = 1.;

    if (year == 2016)
        lumi = 35.917;
    else if (year == 2017)
        lumi = 41.525;
    else if (year == 2018)
        lumi = 59.725;
    else{
        cout << "[ERROR] No supported luminosity value of year: " << year << endl;
        exit(-1);
    }

    return lumi;
}


template <typename T>
auto FindGen(T &df){
    auto nf = df.Define("hardProc",                         "Helper::GenType(mcStatusFlag, 0)")
                .Define("isPrompt",                         "Helper::GenType(mcStatusFlag, 1)")
                
                .Define("isDALEle",                         "abs(mcPID) == 11 && mcMomPID == 25") // higgs dalitz electron
                .Define("isHZGEle",                         "abs(mcPID) == 11 && mcMomPID == 23 && mcGMomPID == 25") // higgs to zg electron
                .Define("isGenEle",                         "(isDALEle || isHZGEle) && hardProc && isPrompt && mcEta < 2.5")

                .Define("isGenPho",                         "abs(mcPID) == 22 && mcMomPID == 25 && hardProc && isPrompt && mcEta < 2.5")
                
                // kick out the events without 2 gen electrons and 1 gen photon
                .Filter("Sum(isGenEle) == 2", "2 gen ele")
                .Filter("Sum(isGenPho) == 1", "1 gen pho")
                
                // gen particle index
                .Define("genEleIdx1",                       "Helper::getIdx(isGenEle, mcPt)[0]")
                .Define("genEleIdx2",                       "Helper::getIdx(isGenEle, mcPt)[1]")
                .Define("genPhoIdx1",                       "Helper::getIdx(isGenPho, mcPt)[0]")
                
                // P4 of gen particle
                .Define("GenEle_lep1",                      "TLorentzVector v; v.SetPtEtaPhiM(mcPt[genEleIdx1], mcEta[genEleIdx1], mcPhi[genEleIdx1], 0.000511); return v;")
                .Define("GenEle_lep2",                      "TLorentzVector v; v.SetPtEtaPhiM(mcPt[genEleIdx2], mcEta[genEleIdx2], mcPhi[genEleIdx2], 0.000511); return v;")
                .Define("GenPho_lep1",                      "TLorentzVector v; v.SetPtEtaPhiM(mcPt[genPhoIdx1], mcEta[genPhoIdx1], mcPhi[genPhoIdx1], 0.000511); return v;")
                .Define("diGenEle",                         "GenEle_lep1 + GenEle_lep2")
                
                .Define("GenElePtRatio",                    "GenEle_lep2.Pt() / GenEle_lep1.Pt()")
                .Define("GenEleDeltaR",                     "GenEle_lep1.DeltaR(GenEle_lep2)")
                .Define("GenElePtSum",                      "GenEle_lep1.Pt() + GenEle_lep2.Pt()");
    
    return nf;
}


template <typename T>
auto GentoReco(T &df){
    auto nf = df.Define("recoEleIdx1",                      "Helper::RecoIdx(GenEle_lep1, eleCalibPt, eleEta, elePhi, 0.000511)")
                .Define("recoEleIdx2",                      "Helper::RecoIdx(GenEle_lep2, eleCalibPt, eleEta, elePhi, 0.000511)")
                .Define("recoPhoIdx1",                      "Helper::RecoIdx(GenPho_lep1, phoCalibEt, phoEta, phoPhi, 0.0)")

                // matching condition
                .Define("Gen2Reco2",                        "recoEleIdx1 != recoEleIdx2 && recoEleIdx1 != -1 && recoEleIdx2 != -1") // resolved 
                .Define("Gen2Reco1",                        "recoEleIdx1 == recoEleIdx2 && recoEleIdx1 != -1 && recoEleIdx2 != -1") // merged 

                // only event with proper reco particles left
                .Filter("(Gen2Reco2 || Gen2Reco1) && recoPhoIdx1 != -1", "good reco") 
                
                // P4 of reco particle 
                .Define("RecoEle_lep1",                     "TLorentzVector v; v.SetPtEtaPhiM(eleCalibPt[recoEleIdx1], eleEta[recoEleIdx1], elePhi[recoEleIdx1], 0.000511); return v;")
                .Define("RecoEle_lep2",                     "TLorentzVector v; v.SetPtEtaPhiM(eleCalibPt[recoEleIdx2], eleEta[recoEleIdx2], elePhi[recoEleIdx2], 0.000511); return v;")
                .Define("RecoPho_lep1",                     "TLorentzVector v; v.SetPtEtaPhiM(phoCalibEt[recoPhoIdx1], phoEta[recoPhoIdx1], phoPhi[recoPhoIdx1], 0.0); return v;");

    return nf;
}


template <typename T>
auto RecotoGSF(T &df){
    auto nf = df.Define("isMainGSF",                        "IsMainGSF(eleD0, gsfD0)")
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
                .Define("eleSubPtTrkDz",                    "MatchIdex(eleSubPtTrkIdx, gsfDz)")

                // conversion vtx matching
                .Define("convVtxIdx1",                      "Helper::RecoConvMatch(eleSCEta[recoEleIdx1], eleSCPhi[recoEleIdx1], eleSCEn[recoEleIdx1], nConv, convNTrks, convVtxX, convVtxY, convVtxZ, convFitPairPX, convFitPairPY, convFitPairPZ, convFitProb)")
                .Define("convVtxIdx2",                      "Helper::RecoConvMatch(eleSCEta[recoEleIdx2], eleSCPhi[recoEleIdx2], eleSCEn[recoEleIdx2], nConv, convNTrks, convVtxX, convVtxY, convVtxZ, convFitPairPX, convFitPairPY, convFitPairPZ, convFitProb)")

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
                
                // trailing gen variables
                .Define("mcPt_lep2",                        "mcPt[genEleIdx2]")
                .Define("mcEta_lep2",                       "mcEta[genEleIdx2]")
                .Define("mcPhi_lep2",                       "mcPhi[genEleIdx2]")
                .Define("mcVtx_lep2",                       "mcVtx[genEleIdx2]")
                .Define("mcVty_lep2",                       "mcVty[genEleIdx2]")
                .Define("mcVtz_lep2",                       "mcVtz[genEleIdx2]")

                // trailing reco variables
                .Define("elePresel_lep2",                   "elePresel[recoEleIdx2]")
                .Define("eleCharge_lep2",                   "eleCharge[recoEleIdx2]")
                .Define("eleChargeConsistent_lep2",         "eleChargeConsistent[recoEleIdx2]")
                .Define("eleConvVeto_lep2",                 "eleConvVeto[recoEleIdx2]")
                .Define("eleEcalDrivenSeed_lep2",           "eleEcalDrivenSeed[recoEleIdx2]")
                .Define("eleEn_lep2",                       "eleEn[recoEleIdx2]")
                .Define("eleSCEn_lep2",                     "eleSCEn[recoEleIdx2]")
                .Define("eleEcalEn_lep2",                   "eleEcalEn[recoEleIdx2]")
                .Define("eleESEnP1_lep2",                   "eleESEnP1[recoEleIdx2]")
                .Define("eleESEnP2_lep2",                   "eleESEnP2[recoEleIdx2]")
                .Define("eleD0_lep2",                       "eleD0[recoEleIdx2]")
                .Define("eleDz_lep2",                       "eleDz[recoEleIdx2]")
                .Define("eleSIP_lep2",                      "eleSIP[recoEleIdx2]")
                .Define("elePt_lep2",                       "elePt[recoEleIdx2]")
                .Define("elePtError_lep2",                  "elePtError[recoEleIdx2]")
                .Define("eleEta_lep2",                      "eleEta[recoEleIdx2]")
                .Define("elePhi_lep2",                      "elePhi[recoEleIdx2]")
                .Define("eleR9_lep2",                       "eleR9[recoEleIdx2]")
                .Define("eleR9Full5x5_lep2",                "eleR9Full5x5[recoEleIdx2]")
                .Define("eleCalibPt_lep2",                  "eleCalibPt[recoEleIdx2]")
                .Define("eleDiffCalibOriPt_lep2",           "(eleCalibPt[recoEleIdx2] - elePt[recoEleIdx2])/elePt[recoEleIdx2]")
                .Define("eleCalibEn_lep2",                  "eleCalibEn[recoEleIdx2]")
                .Define("eleSCEta_lep2",                    "eleSCEta[recoEleIdx2]")
                .Define("eleSCPhi_lep2",                    "eleSCPhi[recoEleIdx2]")
                .Define("eleSCRawEn_lep2",                  "eleSCRawEn[recoEleIdx2]")
                .Define("eleSCEtaWidth_lep2",               "eleSCEtaWidth[recoEleIdx2]")
                .Define("eleSCPhiWidth_lep2",               "eleSCPhiWidth[recoEleIdx2]")
                .Define("eleHoverE_lep2",                   "eleHoverE[recoEleIdx2]")
                .Define("eleEoverP_lep2",                   "eleEoverP[recoEleIdx2]")
                .Define("eleEoverPout_lep2",                "eleEoverPout[recoEleIdx2]")
                .Define("eleEoverPInv_lep2",                "eleEoverPInv[recoEleIdx2]")
                .Define("eleBrem_lep2",                     "eleBrem[recoEleIdx2]")
                .Define("eledEtaAtVtx_lep2",                "eledEtaAtVtx[recoEleIdx2]")
                .Define("eledPhiAtVtx_lep2",                "eledPhiAtVtx[recoEleIdx2]")
                .Define("eleSigmaIEtaIEtaFull5x5_lep2",     "eleSigmaIEtaIEtaFull5x5[recoEleIdx2]")
                .Define("eleSigmaIPhiIPhiFull5x5_lep2",     "eleSigmaIPhiIPhiFull5x5[recoEleIdx2]")
                .Define("eleESEffSigmaRR_lep2",             "eleESEffSigmaRR[recoEleIdx2]")
                .Define("elePFChIso_lep2",                  "elePFChIso[recoEleIdx2]")
                .Define("elePFPhoIso_lep2",                 "elePFPhoIso[recoEleIdx2]")
                .Define("elePFNeuIso_lep2",                 "elePFNeuIso[recoEleIdx2]")
                .Define("elePFPUIso_lep2",                  "elePFPUIso[recoEleIdx2]")
                .Define("elePFClusEcalIso_lep2",            "elePFClusEcalIso[recoEleIdx2]")
                .Define("elePFClusHcalIso_lep2",            "elePFClusHcalIso[recoEleIdx2]")
                .Define("eleIDMVAIso_lep2",                 "eleIDMVAIso[recoEleIdx2]")
                .Define("eleIDMVANoIso_lep2",               "eleIDMVANoIso[recoEleIdx2]")
                .Define("eleTrkdxy_lep2",                   "eleTrkdxy[recoEleIdx2]")
                .Define("eleKFHits_lep2",                   "eleKFHits[recoEleIdx2]")
                .Define("eleKFChi2_lep2",                   "eleKFChi2[recoEleIdx2]")
                .Define("eleGSFChi2_lep2",                  "eleGSFChi2[recoEleIdx2]")
                .Define("eleESEnToRawE_lep2",               "(eleESEnP1[recoEleIdx2] + eleESEnP2[recoEleIdx2])/eleSCRawEn[recoEleIdx2]")
                .Define("nGsfMatchToReco_lep2",             "nGsfMatchToReco[recoEleIdx2]")
                .Define("eleTrk_lep2",                      "TLorentzVector v; v.SetPtEtaPhiM(eleTrkPt[recoEleIdx2], eleTrkEta[recoEleIdx2], eleTrkPhi[recoEleIdx2], 0.000511); return v;")
                .Define("eleTrkPt_lep2",                    "eleTrkPt[recoEleIdx2]")
                .Define("eleTrkEta_lep2",                   "eleTrkEta[recoEleIdx2]")
                .Define("eleTrkPhi_lep2",                   "eleTrkPhi[recoEleIdx2]")
                .Define("eleTrkCharge_lep2",                "eleTrkCharge[recoEleIdx2]")
                .Define("eleTrkLayers_lep2",                "eleTrkLayers[recoEleIdx2]")
                .Define("eleTrkMissHits_lep2",              "eleTrkMissHits[recoEleIdx2]")
                .Define("eleTrkD0_lep2",                    "eleTrkD0[recoEleIdx2]")
                .Define("eleTrkDz_lep2",                    "eleTrkDz[recoEleIdx2]")
                .Define("eleSubTrk_lep2",                   "TLorentzVector v; v.SetPtEtaPhiM(eleSubTrkPt[recoEleIdx2], eleSubTrkEta[recoEleIdx2], eleSubTrkPhi[recoEleIdx2], 0.000511); return v;")
                .Define("eleSubTrkPt_lep2",                 "eleSubTrkPt[recoEleIdx2]")
                .Define("eleSubTrkEta_lep2",                "eleSubTrkEta[recoEleIdx2]")
                .Define("eleSubTrkPhi_lep2",                "eleSubTrkPhi[recoEleIdx2]")
                .Define("eleSubTrkCharge_lep2",             "eleSubTrkCharge[recoEleIdx2]")
                .Define("eleSubTrkLayers_lep2",             "eleSubTrkLayers[recoEleIdx2]")
                .Define("eleSubTrkMissHits_lep2",           "eleSubTrkMissHits[recoEleIdx2]")
                .Define("eleSubTrkD0_lep2",                 "eleSubTrkD0[recoEleIdx2]")
                .Define("eleSubTrkDz_lep2",                 "eleSubTrkDz[recoEleIdx2]")
                .Define("eleSubPtTrk_lep2",                 "TLorentzVector v; v.SetPtEtaPhiM(eleSubPtTrkPt[recoEleIdx2], eleSubPtTrkEta[recoEleIdx2], eleSubPtTrkPhi[recoEleIdx2], 0.000511); return v;")
                .Define("eleSubPtTrkPt_lep2",               "eleSubPtTrkPt[recoEleIdx2]")
                .Define("eleSubPtTrkEta_lep2",              "eleSubPtTrkEta[recoEleIdx2]")
                .Define("eleSubPtTrkPhi_lep2",              "eleSubPtTrkPhi[recoEleIdx2]")
                .Define("eleSubPtTrkCharge_lep2",           "eleSubPtTrkCharge[recoEleIdx2]")
                .Define("eleSubPtTrkLayers_lep2",           "eleSubPtTrkLayers[recoEleIdx2]")
                .Define("eleSubPtTrkMissHits_lep2",         "eleSubPtTrkMissHits[recoEleIdx2]")
                .Define("eleSubPtTrkD0_lep2",               "eleSubPtTrkD0[recoEleIdx2]")
                .Define("eleSubPtTrkDz_lep2",               "eleSubPtTrkDz[recoEleIdx2]")
                .Define("gsfPtSum_lep2",                    "eleTrk_lep2.Pt() + eleSubTrk_lep2.Pt()")
                .Define("gsfPtRatio_lep2",                  "eleSubTrk_lep2.Pt() / eleTrk_lep2.Pt()")
                .Define("gsfDeltaR_lep2",                   "eleTrk_lep2.DeltaR(eleSubTrk_lep2)")
                .Define("gsfMissHitsSum_lep2",              "eleTrkMissHits[recoEleIdx1] + eleSubTrkMissHits[recoEleIdx1]")
                .Define("gsfRelPtRatio_lep2",               "if (nGsfMatchToReco_lep2 > 1) return (eleTrk_lep2 + eleSubTrk_lep2).Pt() / eleCalibPt_lep2; else return eleTrk_lep2.Pt() / eleCalibPt_lep2;")
                .Define("convMatched_lep2",                 "if (nConv > 0 && convVtxIdx2 != -1) return 1; else return 0;")
                .Define("convVtxRadius_lep2",               "if (nConv > 0 && convVtxIdx2 != -1) return convVtxRadius[convVtxIdx2]; else return (float) -999;")
                .Define("convD0_lep2",                      "if (nConv > 0 && convVtxIdx2 != -1) return convD0[convVtxIdx2]; else return (float) -999;")
                .Define("convDz_lep2",                      "if (nConv > 0 && convVtxIdx2 != -1) return convDz[convVtxIdx2]; else return (float) -999;")
                .Define("convL0_lep2",                      "if (nConv > 0 && convVtxIdx2 != -1) return convL0[convVtxIdx2]; else return (float) -999;")
                .Define("convLz_lep2",                      "if (nConv > 0 && convVtxIdx2 != -1) return convLz[convVtxIdx2]; else return (float) -999;")

                // comparison between different methods of picking the second gsf track (only for merged case)
                .Define("diTrk",                            "eleTrk_lep1 + eleSubTrk_lep1")
                .Define("diTrkPtMax",                       "eleTrk_lep1 + eleSubPtTrk_lep1")
                .Define("meeRatio",                         "diTrk.M() / diGenEle.M()")
                .Define("meeRatioPtMax",                    "diTrkPtMax.M() / diGenEle.M()")
                
                // event categorization
                // cat == 1 -> Resolved, cat == 2 -> M2, cat == 3 -> M1
                .Define("category",                         "Helper::make_cat(Gen2Reco2, Gen2Reco1, nGsfMatchToReco_lep1)");

    return nf;
}


void rdfGEN(string infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_m125/*.root", string outfile = "testGen.root", int year = 2017, string era = "2017", string proc = "HDalitz", string prod = "ggF", int HiggsMass = 125){
    TStopwatch time;
    time.Start();

    cout << "[INFO] Read_File(): " << infile.c_str() << endl;
    cout << "[INFO] Save_File(): " << outfile.c_str() << endl;

    //=======================================================//
    // Load TTree into RDataFrame                            //
    //=======================================================//
    ROOT::EnableImplicitMT();
    ROOT::RDataFrame df("ggNtuplizer/EventTree", infile.c_str());
    cout << "[INFO] Process: " << proc << " " << prod << "_" << HiggsMass << "GeV" << endl;

    //========================================================//
    // Calculate the seral weight (add into the RDataFrame)   //
    //========================================================//
    // calculate mc weight
    auto pos = df.Filter("genWeight > 0",   "positive event").Count();
    auto neg = df.Filter("genWeight <= 0",  "negative event").Count();
    const int totalev = pos.GetValue() - neg.GetValue();

    map<string, float> XSmap = XS_HDalitz();
    const float procXS = XSmap[Form("%s_%dGeV", prod.c_str(), HiggsMass)];
    const float mcwei = ((procXS * luminosity(year)) / totalev);
    cout << "[INFO] Number of events with genwei = " << totalev << endl;
    cout << "[INFO] MC weight = " << mcwei << endl;

    // set up the puwei calculator 
    PUWeightCalculator* puCalc = new PUWeightCalculator();;
    puCalc->Init(PUfile(year, "nominal").c_str()); 

    auto get_pu = [&, puCalc](int run, ROOT::RVec<float>& puTrue){
        const float puwei = puCalc->GetWeight(run, puTrue[1]);
        return puwei;
    };

    auto wf = df.Define("puwei",    get_pu,     {"run", "puTrue"})
                .Define("mcwei",    [&]{return mcwei;})
                .Define("procXS",   [&]{return procXS;})
                .Define("instwei",  [&]{return XSmap[Form("ggF_%dGeV", HiggsMass)];})
                .Define("genwei",   "if (genWeight > 0) return (float) 1.; else return (float) -1.;");

    //=======================================================//
    // Do the selections on the RDataFrame                   //
    //=======================================================//
    auto df1 = FindGen(wf);
    auto df2 = GentoReco(df1);
    auto df3 = RecotoGSF(df2);
    auto df4 = DefineFinalVars(df3);
    auto dfFin = df4;

    //========================================================//
    // Visualize the selection results (cut flow, event count)//
    //========================================================//
    cout << "[INFO] Cut flow:" << std::endl;
    auto report = dfFin.Report();
    report->Print();

    int M2 = dfFin.Filter("category == 2").Count().GetValue();
    int M1 = dfFin.Filter("category == 3").Count().GetValue();
    int Re = dfFin.Filter("category == 1").Count().GetValue();
    cout << "[INFO] Number of events:" << std::endl;
    cout << fixed << showpoint;
    cout << ">>> Merged-2Gsf: " << setprecision(2) << M2*100./(M2+M1+Re) << "%(" << M2 << "/" << M2+M1+Re << ")" << endl;
    cout << ">>> Merged-1Gsf: " << setprecision(2) << M1*100./(M2+M1+Re) << "%(" << M1 << "/" << M2+M1+Re << ")" << endl;
    cout << ">>> Resolved:    " << setprecision(2) << Re*100./(M2+M1+Re) << "%(" << Re << "/" << M2+M1+Re << ")" << endl;

    //====================================================//
    // Define the final variables to save to the miniTree //
    //====================================================//
    vector<string> defColNames = dfFin.GetDefinedColumnNames();
    vector<string> Vars;
    for (int i = 0; i < defColNames.size(); i++){
        size_t foundLep1 = defColNames[i].find("_lep1");
        size_t foundLep2 = defColNames[i].find("_lep2");

        if (foundLep1 != string::npos || foundLep2 != string::npos) 
            Vars.push_back(defColNames[i]);
    }
    vector<string> extralVars = {
        "mcwei", "puwei", "genwei", "procXS", "instwei",
        "diTrk", "diTrkPtMax", "meeRatio", "meeRatioPtMax", "category"
    };
    Vars.insert(Vars.end(), extralVars.begin(), extralVars.end());

    dfFin.Snapshot("outTree", outfile.c_str(), Vars);

    cout << "[INFO] Time taken: " << endl;
    time.Stop();
    time.Print();
}