#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <Riostream.h>
#include <stdlib.h>
#include <cmath>
#include <map>

#include <TLorentzVector.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

#include "/home/chenghan/Analysis/Dalitz/muon/plugins/untuplizerv8.h"
#include "/home/chenghan/Analysis/Dalitz/muon/plugins/Utilities.h"

#include "./interface/ElectronSelection.h"
#include "./interface/MC.h"
#include "./interface/GenSelection.h"
#include "./interface/RecoMatching.h"
#include "./interface/Event.h" // This should be declared before Minitree.h
#include "./interface/Minitree.h"
#include "./interface/processbar.h"

using namespace std;

bool debug = false;

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


//  Reference: https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/MicroAOD/plugins/LegacyVertexSelector.cc#L499-L579
Int_t RecoConvMatch(TreeReader &data, int recoele_idx)
{
    Int_t  nConv            = data.GetInt("nConv");
    Int_t* convNTrks        = data.GetPtrInt("convNTrks");
    float* convVtxX         = data.GetPtrFloat("convVtxX");
    float* convVtxY         = data.GetPtrFloat("convVtxY");
    float* convVtxZ         = data.GetPtrFloat("convVtxZ");
    float* convTrksPin0X    = data.GetPtrFloat("convTrksPin0X");
    float* convTrksPin0Y    = data.GetPtrFloat("convTrksPin0Y");
    float* convTrksPin0Z    = data.GetPtrFloat("convTrksPin0Z");
    float* convFitPairPX    = data.GetPtrFloat("convFitPairPX");
    float* convFitPairPY    = data.GetPtrFloat("convFitPairPY");
    float* convFitPairPZ    = data.GetPtrFloat("convFitPairPZ");
    float* convFitProb      = data.GetPtrFloat("convFitProb");

    float* eleSCEta         = data.GetPtrFloat("eleSCEta");
    float* eleSCPhi         = data.GetPtrFloat("eleSCPhi");
    float* eleSCEn          = data.GetPtrFloat("eleSCEn");

    TVector3 SC;
    SC.SetPtEtaPhi(eleSCEn[recoele_idx]/cosh(eleSCEta[recoele_idx]), eleSCEta[recoele_idx], eleSCPhi[recoele_idx]);

    float mindR_2Trks = 999, mindR_1Trks = 999;
    int selected_conversion_index_2Trks = -1, selected_conversion_index_1Trks = -1;
    for (int i = 0; i < nConv; i++)
    {
        TVector3 VtxtoSC;
        VtxtoSC.SetXYZ(SC.x() - convVtxX[i], SC.y() - convVtxY[i], SC.z() - convVtxZ[i]);

        if (convNTrks[i] == 2)
        {
            float pairPt = sqrt(convFitPairPX[i]*convFitPairPX[i] + convFitPairPY[i]*convFitPairPY[i]);
            if (pairPt < 10.)
                continue;
            if (convFitProb[i] < 1e-6)
                continue;

            TVector3 RefPairMo;
            RefPairMo.SetXYZ(convFitPairPX[i], convFitPairPY[i], convFitPairPZ[i]);

            float dR = VtxtoSC.DeltaR(RefPairMo);
            if(dR < mindR_2Trks) 
            {
                mindR_2Trks = dR;
                selected_conversion_index_2Trks = i;
            }
        }
        
        if (convNTrks[i] == 1)
        {
            TVector3 RefPairMo;
            RefPairMo.SetXYZ(convTrksPin0X[i], convTrksPin0Y[i], convTrksPin0Z[i]);

            float dR = VtxtoSC.DeltaR(RefPairMo);
            if(dR < mindR_1Trks) 
            {
                mindR_1Trks = dR;
                selected_conversion_index_1Trks = i;
            }
        }
    }

    if (mindR_2Trks < 0.1)
        return selected_conversion_index_2Trks;

    if (mindR_1Trks < 0.1)
        return selected_conversion_index_1Trks;

    return -1;
}


void RecoLevel(const char *Inpath, const char *outpath, const char *runyear, const char *proc, const char *prod, Int_t HiggsMass, Int_t lepId)
{
    TreeReader data(Inpath);
    TTree *tmptree = data.GetTree();
    printf("[TEST] total number of events = %lli, number of events with nVtx>1 = %lli\n", tmptree->GetEntries(), tmptree->GetEntries("nVtx>1"));

    Minitree outTree(outpath);
    outTree.SetHist();
    outTree.SetMinitree("outTree");

    Long64_t totalEvents = 0;
    Float_t instwei_ = 1.;
    Float_t procXS = 0.;
    Float_t mcwei_ = 0.;
    Float_t Lumi = luminosity(runyear);
    if (strcmp(proc, "HDalitz") == 0)
    {
        std::map<string, Float_t> map_ProdXS = XS_HDalitz();
        mcwei_ = mcwei(data, map_ProdXS[Form("%s_%dGeV", prod, HiggsMass)], Lumi, 1, totalEvents);
        procXS = map_ProdXS[Form("%s_%dGeV", prod, HiggsMass)];
        instwei_ = procXS / map_ProdXS[Form("ggF_%dGeV", HiggsMass)]; // Used in XGBoost training: relative cross-section (ggF as 1, for it being the production with the largest XS)
    }
    else if (strcmp(proc, "DYJetsToLL") == 0)
    {
        Float_t prodXS = 6077.22 * 1000.; // Reference: https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat13TeV
        mcwei_ = mcwei(data, prodXS, Lumi, 1, totalEvents);
        procXS = prodXS;
        instwei_ = 1.;
    }
    else
    {
        cout << "Process [" << proc << "] is not available. Please check!" << endl;
        gSystem->Exit(0);
    }

    printf("[INFO] The physics process of the sample is %s %s, with number of events (w. genwei) = %lld, with mcwei = %f\n", proc, prod, totalEvents, mcwei_);

    std::map<string, Int_t> count_case;
    count_case["Total"] = 0;
    count_case["2Gen2Reco->2Reco2Gsf->2Gsf2Gen"] = 0; // Resolved
    count_case["2Gen1Reco->1Reco2Gsf->2Gsf2Gen"] = 0; // Merged-2Gsf
    // count_case["2Gen1Reco->1Reco2Gsf->1Gsf2Gen"] = 0; // Merged-1MissinGsf
    count_case["2Gen1Reco->1Reco1Gsf->1Gsf2Gen"] = 0; // Merged-1MissinGsf
    count_case["Others"] = 0;                         // Not properly reconstructed
    std::map<string, Int_t> count_category;
    count_category["Total"] = 0;
    count_category["Resolved"] = 0;
    count_category["Merged2Gsf"] = 0;
    count_category["Merged1MissingGsf"] = 0;
    count_category["NPR"] = 0; // Not properly reconstructed

    Int_t showev = int(TMath::Power(10, int(TMath::Log10(data.GetEntriesFast() / 100.))));

    std::map<string, Float_t> map_expyield;
    map_expyield["Total"] = 0.;
    map_expyield["Resolved"] = 0.;
    map_expyield["Merged2Gsf"] = 0.;
    map_expyield["Merged1MissingGsf"] = 0.;
    map_expyield["NPR"] = 0.;

    for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++)
    {
        if (ev % showev == 0)
            printprocess(ev, showev, data.GetEntriesFast());

        Event event = {};

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
        float *mcPt = data.GetPtrFloat("mcPt");
        float *mcEta = data.GetPtrFloat("mcEta");
        float *mcPhi = data.GetPtrFloat("mcPhi");
        float *mcMomMass = data.GetPtrFloat("mcMomMass");
        float *mcMomPt = data.GetPtrFloat("mcMomPt");
        float *mcMomEta = data.GetPtrFloat("mcMomEta");
        float *mcMomPhi = data.GetPtrFloat("mcMomPhi");
        int nEle = data.GetInt("nEle"); // MC

        vector<Int_t> vec_LepIdx, vec_PhoIdx;
        vector<TLorentzVector> vec_Lep, vec_Pho;
        vec_LepIdx.clear();
        vec_PhoIdx.clear();
        vec_Lep.clear();
        vec_Pho.clear();

        GenSelection(data, vec_LepIdx, vec_PhoIdx, lepId, proc);
        for (size_t i = 0; i < vec_LepIdx.size(); i++)
        {
            TLorentzVector tmpobj;
            tmpobj.SetPtEtaPhiM(mcPt[vec_LepIdx[i]], mcEta[vec_LepIdx[i]], mcPhi[vec_LepIdx[i]], genlepmass(lepId));
            vec_Lep.push_back(tmpobj);
        }
        for (size_t i = 0; i < vec_PhoIdx.size(); i++)
        {
            TLorentzVector tmpobj;
            tmpobj.SetPtEtaPhiM(mcPt[vec_PhoIdx[i]], mcEta[vec_PhoIdx[i]], mcPhi[vec_PhoIdx[i]], 0.0);
            vec_Pho.push_back(tmpobj);
        }

        bool cut = false;
        if (strcmp(proc, "HDalitz") == 0)
        {
            cut = (nEle < 1) ||
                  (vec_Lep.size() != 2 || vec_Pho.size() != 1);
        }
        else if (strcmp(proc, "DYJetsToLL") == 0)
        {
            cut = (nEle < 1 || vec_Lep.size() != 2);
        }
        else if (strcmp(proc, "gjets") == 0)
        {
            cut = (nEle < 1 || vec_Lep.size() != 2);
        }
        else
        {
            cout << "Process [" << proc << "] is not available. Please check!" << endl;
            gSystem->Exit(0);
        }

        if (cut)
            continue;
        
        count_case["Total"] += 1;
        count_category["Total"] += 1;

        // check the order of vec_lep -> make sure the leading lepton has index 0 and trailing lepton has 1
        if (vec_Lep[0].Pt() < vec_Lep[1].Pt())
        {
            std::reverse(vec_Lep.begin(), vec_Lep.end());
            std::reverse(vec_LepIdx.begin(), vec_LepIdx.end());
        }

        Int_t GenReco_case = -1, RecoGsf_case = -1, GsfGen_case = -1;
        std::map<string, Int_t> GenMatchReco_idx = genreco_match(data, vec_Lep, GenReco_case);
        std::map<string, vector<Int_t>> RecoMatchGsf_idx = recogsf_match(data, GenMatchReco_idx, RecoGsf_case);
        std::map<string, Int_t> GenMatchGSF_idx = gsfgen_match(data, vec_Lep, RecoMatchGsf_idx, GsfGen_case);

        if (debug)
            debuginfo(data, ev, vec_Lep, GenMatchReco_idx, RecoMatchGsf_idx, GenMatchGSF_idx, GenReco_case, RecoGsf_case, GsfGen_case);

        //BC matching 
        float *elePt = data.GetPtrFloat("eleCalibPt");
        float *eleEta = data.GetPtrFloat("eleEta");
        float *elePhi = data.GetPtrFloat("elePhi");

        TLorentzVector tmpReco_lep1, tmpReco_lep2;
        tmpReco_lep1.SetPtEtaPhiM(elePt[GenMatchReco_idx["GenLep1"]], eleEta[GenMatchReco_idx["GenLep1"]], elePhi[GenMatchReco_idx["GenLep1"]], 0.000510999);
        tmpReco_lep2.SetPtEtaPhiM(elePt[GenMatchReco_idx["GenLep2"]], eleEta[GenMatchReco_idx["GenLep2"]], elePhi[GenMatchReco_idx["GenLep2"]], 0.000510999);
        event.BCIdxLep1 = RecoBCMatch(data, tmpReco_lep1);
        event.BCIdxLep2 = RecoBCMatch(data, tmpReco_lep2);

        if (event.BCIdxLep1.size() < 1)
            continue;

        if (event.BCIdxLep2.size() < 1)
            continue;


        // Conversion matching
        event.ConvIdxLep1 = RecoConvMatch(data, GenMatchReco_idx["GenLep1"]);
        event.ConvIdxLep2 = RecoConvMatch(data, GenMatchReco_idx["GenLep2"]);


        const char *Category_str = "";
        event.category = Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str);
        event.GenReco_case = GenReco_case;
        event.RecoGsf_case = RecoGsf_case;
        event.GsfGen_case = GsfGen_case;
        Category_count(GenReco_case, RecoGsf_case, GsfGen_case, count_case, count_category);

        event.HLTEleMuX = HLTEleMuX;
        event.HLTEleMuXIsPrescaled = HLTEleMuXIsPrescaled;
        event.HLTPho = HLTPho;
        event.HLTPhoIsPrescaled = HLTPhoIsPrescaled;
        event.rho = rho;
        event.rhoCentral = rhoCentral;
        event.nVtx = nVtx;
        event.nGoodVtx = nGoodVtx;
        event.isPVGood = isPVGood;
        event.totalEvents = totalEvents;
        event.instwei = instwei_;

        event.iGenLep[0] = vec_LepIdx[0];
        event.iGenLep[1] = vec_LepIdx[1];
        event.iGenPho = vec_PhoIdx[0];
        event.iele[0] = GenMatchReco_idx["GenLep1"];
        event.iele[1] = GenMatchReco_idx["GenLep2"];
        event.GsfIdxLep1 = RecoMatchGsf_idx["GenLep1"];
        event.GsfIdxLep2 = RecoMatchGsf_idx["GenLep2"];
        event.igsf[0] = GenMatchGSF_idx["GenLep1"];
        event.igsf[1] = GenMatchGSF_idx["GenLep2"];

        map_expyield["Total"] += genwei_;
        if (Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str) == 1)
            map_expyield["Resolved"] += genwei_;
        if (Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str) == 2)
            map_expyield["Merged2Gsf"] += genwei_;
        if (Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str) == 3)
            map_expyield["Merged1MissingGsf"] += genwei_;
        if (Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str) == 4)
            map_expyield["NPR"] += genwei_;

        if (Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str) != 1 &&
            Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str) != 2 &&
            Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str) != 3 &&
            Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str) != 4)
            cout << "[WARNING] category " << Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str) << " isn't designed category. Please check!" << endl;

        event.mcwei = mcwei_;
        event.procXS = procXS;
        event.genwei = genwei_;
        outTree.FillMinitree(data, event);
    }

    for (std::map<string, Int_t>::iterator it = count_case.begin(); it != count_case.end(); ++it)
    {
        printf("%s: %d\n", it->first.c_str(), count_case[it->first]);
    }

    for (std::map<string, Int_t>::iterator it = count_category.begin(); it != count_category.end(); ++it)
    {
        printf("%s: %d\n", it->first.c_str(), count_category[it->first]);
    }
    cout << proc << " " << prod << " Total: " << map_expyield["Total"] * mcwei_ << "; Resolved: " << map_expyield["Resolved"] * mcwei_ << "; Merged2Gsf: " << map_expyield["Merged2Gsf"] * mcwei_ << "; Merged1Gsf: " << map_expyield["Merged1MissingGsf"] * mcwei_ << "; NPR: " << map_expyield["NPR"] * mcwei_ << "; sum = " << map_expyield["Resolved"] * mcwei_ + map_expyield["Merged2Gsf"] * mcwei_ + map_expyield["Merged1MissingGsf"] * mcwei_ + map_expyield["NPR"] * mcwei_ << endl;
    if (strcmp(proc, "HDalitz") == 0)

    cout << "||================================================================================================================||" << endl;
}
