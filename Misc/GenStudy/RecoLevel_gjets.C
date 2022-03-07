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
#include "./interface/Event_gjets.h" // This should be declared before Minitree.h
#include "./interface/Minitree_gjets.h"
#include "./interface/processbar.h"

using namespace std;

bool debug = false;

Int_t Gen2Reco(TreeReader &data, TLorentzVector genobj, const char *recotype)
{
    // RECO-level information
    Int_t nReco = 0;
    float *recoPt = 0;
    float *recoEta = 0;
    float *recoPhi = 0;
    Float_t Mparticle = 0;
    if (strcmp(recotype, "ele") == 0)
    {
        nReco = data.GetInt("nEle");
        recoPt = data.GetPtrFloat("eleCalibPt");
        recoEta = data.GetPtrFloat("eleEta");
        recoPhi = data.GetPtrFloat("elePhi");
        Mparticle = 0.000510999;
    }
    else if (strcmp(recotype, "pho") == 0)
    {
        nReco = data.GetInt("nPho");
        recoPt = data.GetPtrFloat("phoCalibEt");
        recoEta = data.GetPtrFloat("phoEta");
        recoPhi = data.GetPtrFloat("phoPhi");
        Mparticle = 0.0;
    }
    else if (strcmp(recotype, "mu") == 0)
    {
        nReco = data.GetInt("nMu");
        recoPt = data.GetPtrFloat("muPt");
        recoEta = data.GetPtrFloat("muEta");
        recoPhi = data.GetPtrFloat("muPhi");
        Mparticle = 0.105658;
    }
    else
    {
        printf("[INFO] recotype [%s] is not ele/pho/mu. Please check! Use ele as default.\n", recotype);
        nReco = data.GetInt("nEle");
        recoPt = data.GetPtrFloat("eleCalibPt");
        recoEta = data.GetPtrFloat("eleEta");
        recoPhi = data.GetPtrFloat("elePhi");
        Mparticle = 0.000510999;
    }

    if (debug)
        printf("[Gen2Reco] nReco=%d\n", nReco);
    int temp_index = -1;
    float temp_ptratio = 999.;
    for (Int_t i = 0; i < nReco; i++)
    {
        TLorentzVector tmpele;
        tmpele.SetPtEtaPhiM(recoPt[i], recoEta[i], recoPhi[i], Mparticle);

        if (debug)
        {
            printf("[Gen2Reco] genobj (pT,eta,phi)=(%f,%f,%f); iReco=%d (pT,eta,phi)=(%f,%f,%f); dR=%f, (pTratio-1)=%f\n", genobj.Pt(), genobj.Eta(), genobj.Phi(), i, tmpele.Pt(), tmpele.Eta(), tmpele.Phi(), genobj.DeltaR(tmpele), fabs((tmpele.Pt() / genobj.Pt()) - 1.));
        }

        if (tmpele.DeltaR(genobj) > 0.1)
            continue;

        if (fabs((tmpele.Pt() / genobj.Pt()) - 1.) < temp_ptratio)
        {
            temp_ptratio = fabs((tmpele.Pt() / genobj.Pt()) - 1.);
            temp_index = i;
            continue;
        }
        else
            continue;
    }

    return temp_index;
}

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

void RecoLevel_gjets(const char *Inpath, const char *outpath, const char *runyear, const char *proc, const char *prod, Int_t lepId)
{
    TreeReader data(Inpath);

    Minitree_gjets outTree(outpath);
    outTree.SetMinitree("outTree");

    Long64_t totalEvents = 0;
    Float_t procXS = 0;
    Float_t mcwei_ = 0;
    Float_t instwei_ = 1.;
    Float_t Lumi = luminosity(runyear);

    std::map<string, Float_t> xs_GJets = XS_GJets();
    procXS = xs_GJets[prod];
    mcwei_ = mcwei(data, xs_GJets[prod], Lumi, 1, totalEvents);
    if (strcmp(prod, "pt20_MGG_40to80") == 0 || strcmp(prod, "pt20to40_MGG_80toInf") == 0 || strcmp(prod, "pt40_MGG_80toInf") == 0)
        instwei_ = procXS / (154500 * 0.0239 * 1000);
    else
        instwei_ = 1.;

    printf("[INFO] The physics process of the sample is %s %s, with number of events (w. genwei) = %lld, with mcwei = %f\n", proc, prod, totalEvents, mcwei_);

    Int_t showev = int(TMath::Power(10, int(TMath::Log10(data.GetEntriesFast() / 100.))));

    // Long64_t evt2print = 310060;
    // for (Long64_t ev = evt2print; ev < evt2print + 1; ev++)
    for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++)
    {
        if (!debug && ev % showev == 0)
            printprocess(ev, showev, data.GetEntriesFast());

        Event_gjets event = {};

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

        for (Int_t i = 0; i < nMC; i++)
        {
            if (mcPID[i] == 22)
            {
                if (debug)
                    printf("[INFO] ev = %lld; mcPID = %d; mcMomPID = %d; mcGMomPID = %d; mcMomPt = %f \n", ev, mcPID[i], mcMomPID[i], mcGMomPID[i], mcMomPt[i]);

                TLorentzVector tmpobj;
                tmpobj.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0.0);
                Int_t RecoEleIdx = Gen2Reco(data, tmpobj, "ele"); // converted photon -> in reco electron collection

                if (RecoEleIdx == -1)
                    continue;

                float *eleCalibPt = data.GetPtrFloat("eleCalibPt");
                float *eleEta = data.GetPtrFloat("eleEta");
                float *elePhi = data.GetPtrFloat("elePhi");
                TLorentzVector recoele;
                recoele.SetPtEtaPhiM(eleCalibPt[RecoEleIdx], eleEta[RecoEleIdx], elePhi[RecoEleIdx], genlepmass(lepId));

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
                event.iGen = i;
                event.iReco = RecoEleIdx;
                event.GsfIdx = RecoGsfMatch(data, recoele);
                event.BCIdx = RecoBCMatch(data, recoele);

                if (event.GsfIdx.size() < 1)
                    continue;
                if (event.BCIdx.size() < 1)
                    continue;

                if (debug)
                    printf("[INFO] ev=%lld; iGen = %d; iReco=%d; number of Gsf tracks that match to RECO ele=%zu\n", ev, i, RecoEleIdx, event.GsfIdx.size());

                outTree.FillMinitree(data, event);
            }
            else
                continue;
        }
        if (debug)
            printf("------------------------------------------------\n");
    }
}