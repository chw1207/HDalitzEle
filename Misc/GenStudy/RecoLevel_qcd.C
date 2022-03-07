#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <Riostream.h>
#include <stdlib.h>
#include <map>

#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

#include "/home/chenghan/Analysis/Dalitz/muon/plugins/untuplizerv8.h"
#include "/home/chenghan/Analysis/Dalitz/muon/plugins/Utilities.h"

#include "./interface/ElectronSelection.h"
#include "./interface/MC.h"
#include "./interface/Event_qcd.h" // This should be declared before Minitree.h
#include "./interface/Minitree_qcd.h"
#include "./interface/processbar.h"
using namespace std;

bool debug = false;


vector<Int_t> RecoGsfMatch(TreeReader &data, TLorentzVector recoele)
{
    Int_t nGSFTrk = data.GetInt("nGSFTrk");
    float *gsfPt = data.GetPtrFloat("gsfPt");
    float *gsfEta = data.GetPtrFloat("gsfEta");
    float *gsfPhi = data.GetPtrFloat("gsfPhi");

    vector<Int_t> v_matchedidx;
    v_matchedidx.clear();
    for (int i = 0; i < nGSFTrk; i++)
    {
        if (debug)
            printf("[RecoGsfMatch] nGsf=%d; iGsf=%d; Gsf (pT,eta,phi)=(%f,%f,%f); |dEta|=%f; |dPhi|=%f\n", nGSFTrk, i, gsfPt[i], gsfEta[i], gsfPhi[i], fabs(deltaEta(gsfEta[i], recoele.Eta())), fabs(deltaPhi(gsfPhi[i], recoele.Phi())));

        if (gsfPt[i]/recoele.Pt() > 5) continue; // remove the gsf track with unreaonable high Pt

        // Reference: https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/RecoEgamma/EgammaElectronProducers/plugins/GsfElectronProducer.cc#L188-L191
        if (fabs(deltaEta(gsfEta[i], recoele.Eta())) > 0.02 ||
            fabs(deltaPhi(gsfPhi[i], recoele.Phi())) > 0.15)
            continue;

        v_matchedidx.push_back(i);
    }

    return v_matchedidx;
}

vector<Int_t> RecoBCMatch(TreeReader &data, TLorentzVector recoele)
{
    Int_t nBC       = data.GetInt("nBC");
    float *bcEta    = data.GetPtrFloat("bcEta");
    float *bcPhi    = data.GetPtrFloat("bcPhi");

    vector<Int_t> v_matchedidx;
    v_matchedidx.clear();
    for (int i = 0; i < nBC; i++)
    {
        // Reference: https://github.com/cms-sw/cmssw/blob/fa58f0d5d3345359eac47f4ae964b7d7379d5fbb/RecoEcal/EgammaClusterProducers/python/particleFlowSuperClusterECALBox_cfi.py#L48-L52
        if (fabs(deltaEta(bcEta[i], recoele.Eta())) > 0.04 ||
            fabs(deltaPhi(bcPhi[i], recoele.Phi())) > 0.28)
            continue;

        v_matchedidx.push_back(i);
    }

    return v_matchedidx;
}

void RecoLevel_qcd(const char *Inpath, const char *outpath, const char *runyear, const char *proc, const char *prod)
{
    TreeReader data(Inpath);

    Minitree_qcd outTree(outpath);
    outTree.SetMinitree("outTree");

    Long64_t totalEvents = 0;
    double procXS = 0;
    Float_t mcwei_ = 0;
    Float_t instwei_ = 1.;
    Float_t Lumi = luminosity(runyear);

    std::map<string, double> xs_QCD = XS_QCD();
    procXS = xs_QCD[prod];
    mcwei_ = mcwei(data, xs_QCD[prod], Lumi, 1, totalEvents);

    if (strcmp(prod, "HT100to200") != 0)
        instwei_ = procXS / xs_QCD["HT100to200"];
    else
        instwei_ = 1.;

    printf("[INFO] The physics process of the sample is %s %s, with number of events (w. genwei) = %lld, with mcwei = %f\n", proc, prod, totalEvents, mcwei_);
    
    Int_t showev = int(TMath::Power(10, int(TMath::Log10(data.GetEntriesFast() / 100.))));

    // Long64_t evt2print = 310060;
    // for (Long64_t ev = 0; ev < 1000; ev++)
    for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++)
    {
        if (!debug && ev % showev == 0)
            printprocess(ev, showev, data.GetEntriesFast());

        Event_qcd event = {};

        data.GetEntry(ev);

        Float_t GenWeight = data.GetFloat("genWeight");
        Float_t genwei_ = (GenWeight > 0.) ? 1. : -1.;

        ULong64_t HLTEleMuX = (ULong64_t)data.GetLong64("HLTEleMuX");
        ULong64_t HLTEleMuXIsPrescaled = (ULong64_t)data.GetLong64("HLTEleMuXIsPrescaled");
        ULong64_t HLTPho = (ULong64_t)data.GetLong64("HLTPho");
        ULong64_t HLTPhoIsPrescaled = (ULong64_t)data.GetLong64("HLTPhoIsPrescaled");
        Float_t rho = data.GetFloat("rho");
        Float_t rhoCentral = data.GetFloat("rhoCentral");
        Int_t nVtx = data.GetInt("nVtx");
        Int_t nGoodVtx = data.GetInt("nGoodVtx");
        Bool_t isPVGood = data.GetBool("isPVGood");

        // GEN information
        int nMC = data.GetInt("nMC"); // MC
        int *mcPID = data.GetPtrInt("mcPID");
        int *mcMomPID = data.GetPtrInt("mcMomPID");
        float *mcMomMass = data.GetPtrFloat("mcMomMass");
        float *mcMomPt = data.GetPtrFloat("mcMomPt");
        float *mcMomEta = data.GetPtrFloat("mcMomEta");
        float *mcMomPhi = data.GetPtrFloat("mcMomPhi");
        int *mcGMomPID = data.GetPtrInt("mcGMomPID");
        float *mcPt = data.GetPtrFloat("mcPt");
        float *mcEta = data.GetPtrFloat("mcEta");
        float *mcPhi = data.GetPtrFloat("mcPhi");
        float *mcMass = data.GetPtrFloat("mcMass");

        // RECO information
        int nEle = data.GetInt("nEle");
        float *eleCalibPt = data.GetPtrFloat("eleCalibPt");
        float *eleEta = data.GetPtrFloat("eleEta");
        float *elePhi = data.GetPtrFloat("elePhi");
        float Mele = 0.000510999;

        vector<int> ele_ind;
        ele_ind.clear();
        map<int, int> ele2gen_ind;
        ele2gen_ind.clear();
        for (int i = 0; i < nEle; i++){
            // printf("---------------------------------------------\n");
            TLorentzVector tmpele;
            tmpele.SetPtEtaPhiM(eleCalibPt[i], eleEta[i], elePhi[i], Mele);
            // printf("[RECO] elePt: %f, eleEta: %f, elePhi: %f\n", tmpele.Pt(), tmpele.Eta(), tmpele.Phi());
            
            int temp_index = -1;
            float temp_ptratio = 999.;
            for (int j = 0; j < nMC; j++){
                
                TLorentzVector genobj;
                genobj.SetPtEtaPhiM(mcPt[j], mcEta[j], mcPhi[j], mcMass[j]);
                // printf("[GEN] Pt: %f, Eta: %f, Phi: %f, dR: %f\n", genobj.Pt(), genobj.Eta(), genobj.Phi(), tmpele.DeltaR(genobj));

                if (tmpele.DeltaR(genobj) > 0.1)
                    continue;

                if (fabs(mcPID[j]) != 11)
                    continue;
                // printf("PID: %d, MomPID: %d, dR: %f\n", mcPID[j], mcMomPID[j], tmpele.DeltaR(genobj));
                
                if (fabs((tmpele.Pt() / genobj.Pt()) - 1.) < temp_ptratio){
                    temp_ptratio = fabs((tmpele.Pt() / genobj.Pt()) - 1.);
                    temp_index = j;
                    continue;
                }
                else
                    continue;
                
            }
            if (temp_index == -1) // no matched gen electron
                continue;
            if (mcMomPID[temp_index] == 22) // reject converted photon 
                continue;
            // printf("ind, %d, PID: %d, MomPID: %d\n", temp_index, mcPID[temp_index], mcMomPID[temp_index]);

            ele_ind.push_back(i);
            ele2gen_ind[i] = temp_index;
        }

        if (ele_ind.size() < 1)
            continue;

        for (size_t i = 0; i < ele_ind.size(); i++){
            int RecoEleIdx = i;
            TLorentzVector recoele;
            recoele.SetPtEtaPhiM(eleCalibPt[RecoEleIdx], eleEta[RecoEleIdx], elePhi[RecoEleIdx], Mele);
            event.procXS = procXS;
            event.mcwei = mcwei_;
            event.genwei = genwei_;
            event.instwei = instwei_;
            event.totalEvents = totalEvents;
            event.HLTEleMuX = HLTEleMuX;
            event.HLTEleMuXIsPrescaled = HLTEleMuXIsPrescaled;
            event.HLTPho = HLTPho;
            event.HLTPhoIsPrescaled = HLTPhoIsPrescaled;
            event.rho = rho;
            event.rhoCentral = rhoCentral;
            event.nVtx = nVtx;
            event.nGoodVtx = nGoodVtx;
            event.isPVGood = isPVGood;
            event.iGen = ele2gen_ind[RecoEleIdx];
            event.iReco = RecoEleIdx;
            event.GsfIdx = RecoGsfMatch(data, recoele);
            event.BCIdx = RecoBCMatch(data, recoele);

            if (event.GsfIdx.size() < 1)
                continue;
            if (event.BCIdx.size() < 1)
                continue;

            outTree.FillMinitree(data, event);
        }
    }
    
}