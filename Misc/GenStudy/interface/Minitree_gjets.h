#include <TLorentzVector.h>

class Minitree_gjets
{
    TFile *fo;
    TTree *_tree;

    ULong64_t HLTEleMuX, HLTEleMuXIsPrescaled, HLTPho, HLTPhoIsPrescaled;
    Float_t rho, rhoCentral;
    Int_t nVtx, nGoodVtx;
    Bool_t isPVGood;
    // Generator-level information
    Float_t mcPt_, mcEta_, mcPhi_, mcVtx_, mcVty_, mcVtz_;
    Int_t mcMomPID_, mcGMomPID_;
    UShort_t mcStatusFlag_;
    Float_t procXS, mcwei, genwei;
    Long64_t totalEvents;
    Float_t instwei;
    // Reco-level information
    Bool_t elePresel;
    Int_t eleCharge_, eleChargeConsistent_, eleConvVeto_, eleMissHits_, eleEcalDrivenSeed_;
    Float_t eleEn_, eleSCEn_, eleEcalEn_, eleESEnP1_, eleESEnP2_, eleD0_, eleDz_, eleSIP_, elePt_, elePtError_, eleEta_, elePhi_, eleR9_, eleCalibPt_, eleDiffCalibOriPt_, eleCalibEn_, eleSCEta_, eleSCPhi_, eleSCRawEn_, eleSCEtaWidth_, eleSCPhiWidth_, eleHoverE_, eleEoverP_, eleEoverPout_, eleEoverPInv_, eleBrem_, eledEtaAtVtx_, eledPhiAtVtx_, eleSigmaIEtaIEtaFull5x5_, eleSigmaIPhiIPhiFull5x5_, eleESEffSigmaRR_, elePFChIso_, elePFPhoIso_, elePFNeuIso_, elePFPUIso_, elePFClusEcalIso_, elePFClusHcalIso_, eleIDMVAIso_, eleIDMVANoIso_, eleR9Full5x5_, eleTrkdxy_, eleKFHits_, eleKFChi2_, eleGSFChi2_;
    ULong64_t eleFiredSingleTrgs_, eleFiredDoubleTrgs_, eleFiredL1Trgs_;
    UShort_t eleIDbit_;
    // Gsf track information
    vector<Int_t> GsfIdx, gsfMissHits_reco;
    vector<Float_t> gsfPt_reco, gsfEta_reco, gsfPhi_reco;
    // Float_t gsfPt_, gsfEta_, gsfPhi_;
    Int_t nGsfMatchToReco;
    Float_t gsfPtSum, gsfPtRatio, gsfDeltaR; // Stored only when there are at least 2 gsf tracks being associated to a reco electron
    Int_t gsfMissHitsSum;

    float gsfPtSum2, gsfRelPtRatio, gsfPerpEleSum, gsfDiffPtRatio;
    
    // BC information
    float circularity;

    // preshower variables
    float eleESEnToRawE_; 
    // float eleESEffSigmaRR

public:
    Minitree_gjets(){};
    Minitree_gjets(const char *outpath)
    {
        fo = TFile::Open(outpath, "RECREATE");
        cout << "outpath [" << outpath << "] is opened." << endl;
        fo->cd();
    };

    ~Minitree_gjets()
    {
        _tree->Write("", TObject::kOverwrite);

        delete _tree;
        delete fo;
    };

    void Refresh();
    void SetMinitree(TString treename);
    void FillMinitree(TreeReader &data, Event_gjets &event);
};

void Minitree_gjets::Refresh()
{
    GsfIdx.clear();
    gsfPt_reco.clear();
    gsfEta_reco.clear();
    gsfPhi_reco.clear();
    gsfMissHits_reco.clear();
}

void Minitree_gjets::SetMinitree(TString treename)
{
    _tree = new TTree(treename, "minitree");

    // Category
    _tree->Branch("HLTEleMuX", &HLTEleMuX);
    _tree->Branch("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled);
    _tree->Branch("HLTPho", &HLTPho);
    _tree->Branch("HLTPhoIsPrescaled", &HLTPhoIsPrescaled);
    _tree->Branch("rho", &rho);
    _tree->Branch("rhoCentral", &rhoCentral);
    _tree->Branch("nVtx", &nVtx);
    _tree->Branch("nGoodVtx", &nGoodVtx);
    _tree->Branch("isPVGood", &isPVGood);
    _tree->Branch("totalEvents", &totalEvents);
    _tree->Branch("instwei", &instwei);
    // Generator information
    _tree->Branch("procXS", &procXS);
    _tree->Branch("mcwei", &mcwei);
    _tree->Branch("genwei", &genwei);
    _tree->Branch("mcPt", &mcPt_);
    _tree->Branch("mcEta", &mcEta_);
    _tree->Branch("mcPhi", &mcPhi_);
    _tree->Branch("mcVtx", &mcVtx_);
    _tree->Branch("mcVty", &mcVty_);
    _tree->Branch("mctz", &mcVtz_);
    _tree->Branch("mcMomPID", &mcMomPID_);
    _tree->Branch("mcGMomPID", &mcGMomPID_);
    _tree->Branch("mcStatusFlag", &mcStatusFlag_);
    // Reco information
    _tree->Branch("elePresel", &elePresel);
    _tree->Branch("eleCharge", &eleCharge_);
    _tree->Branch("eleChargeConsistent", &eleChargeConsistent_);
    _tree->Branch("eleConvVeto", &eleConvVeto_);
    _tree->Branch("eleMissHits", &eleMissHits_);
    _tree->Branch("eleEcalDrivenSeed", &eleEcalDrivenSeed_);
    _tree->Branch("eleEn", &eleEn_);
    _tree->Branch("eleSCEn", &eleSCEn_);
    _tree->Branch("eleEcalEn", &eleEcalEn_);
    _tree->Branch("eleESEnP1", &eleESEnP1_);
    _tree->Branch("eleESEnP2", &eleESEnP2_);
    _tree->Branch("eleD0", &eleD0_);
    _tree->Branch("eleDz", &eleDz_);
    _tree->Branch("eleSIP", &eleSIP_);
    _tree->Branch("elePt", &elePt_);
    _tree->Branch("elePtError", &elePtError_);
    _tree->Branch("eleEta", &eleEta_);
    _tree->Branch("elePhi", &elePhi_);
    _tree->Branch("eleR9", &eleR9_);
    _tree->Branch("eleCalibPt", &eleCalibPt_);
    _tree->Branch("eleDiffCalibOriPt", &eleDiffCalibOriPt_);
    _tree->Branch("eleCalibEn", &eleCalibEn_);
    _tree->Branch("eleSCEta", &eleSCEta_);
    _tree->Branch("eleSCPhi", &eleSCPhi_);
    _tree->Branch("eleSCRawEn", &eleSCRawEn_);
    _tree->Branch("eleSCEtaWidth", &eleSCEtaWidth_);
    _tree->Branch("eleSCPhiWidth", &eleSCPhiWidth_);
    _tree->Branch("eleHoverE", &eleHoverE_);
    _tree->Branch("eleEoverP", &eleEoverP_);
    _tree->Branch("eleEoverPout", &eleEoverPout_);
    _tree->Branch("eleEoverPInv", &eleEoverPInv_);
    _tree->Branch("eleBrem", &eleBrem_);
    _tree->Branch("eledEtaAtVtx", &eledEtaAtVtx_);
    _tree->Branch("eledPhiAtVtx", &eledPhiAtVtx_);
    _tree->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
    _tree->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);
    _tree->Branch("eleESEffSigmaRR", &eleESEffSigmaRR_);
    _tree->Branch("elePFChIso", &elePFChIso_);
    _tree->Branch("elePFPhoIso", &elePFPhoIso_);
    _tree->Branch("elePFNeuIso", &elePFNeuIso_);
    _tree->Branch("elePFPUIso", &elePFPUIso_);
    _tree->Branch("elePFClusEcalIso", &elePFClusEcalIso_);
    _tree->Branch("elePFClusHcalIso", &elePFClusHcalIso_);
    _tree->Branch("eleIDMVAIso", &eleIDMVAIso_);
    _tree->Branch("eleIDMVANoIso", &eleIDMVANoIso_);
    _tree->Branch("eleR9Full5x5", &eleR9Full5x5_);
    _tree->Branch("eleTrkdxy", &eleTrkdxy_);
    _tree->Branch("eleKFHits", &eleKFHits_);
    _tree->Branch("eleKFChi2", &eleKFChi2_);
    _tree->Branch("eleGSFChi2", &eleGSFChi2_);
    _tree->Branch("eleFiredSingleTrgs", &eleFiredSingleTrgs_);
    _tree->Branch("eleFiredDoubleTrgs", &eleFiredDoubleTrgs_);
    _tree->Branch("eleFiredL1Trgs", &eleFiredL1Trgs_);
    _tree->Branch("eleIDbit", &eleIDbit_);
    // Gsf track information
    _tree->Branch("GsfIdx", &GsfIdx);
    _tree->Branch("gsfPt_reco", &gsfPt_reco);
    _tree->Branch("gsfEta_reco", &gsfEta_reco);
    _tree->Branch("gsfPhi_reco", &gsfPhi_reco);
    _tree->Branch("gsfMissHits_reco", &gsfMissHits_reco);
    // _tree->Branch("gsfPt", &gsfPt_);
    // _tree->Branch("gsfEta", &gsfEta_);
    // _tree->Branch("gsfPhi", &gsfPhi_);
    _tree->Branch("nGsfMatchToReco", &nGsfMatchToReco);
    _tree->Branch("gsfPtSum", &gsfPtSum);
    _tree->Branch("gsfPtRatio", &gsfPtRatio);
    _tree->Branch("gsfDeltaR", &gsfDeltaR);
    _tree->Branch("gsfMissHitsSum", &gsfMissHitsSum);
    _tree->Branch("gsfPtSum2", &gsfPtSum2);
    _tree->Branch("gsfRelPtRatio", &gsfRelPtRatio);
    _tree->Branch("gsfPerpEleSum", &gsfPerpEleSum);
    _tree->Branch("gsfDiffPtRatio", &gsfDiffPtRatio);
    _tree->Branch("circularity", &circularity);
    _tree->Branch("eleESEnToRawE", &eleESEnToRawE_);
}

void Minitree_gjets::FillMinitree(TreeReader &data, Event_gjets &event)
{
    Refresh();
    // GEN information
    int nMC = data.GetInt("nMC"); // MC
    int *mcPID = data.GetPtrInt("mcPID");
    int *mcMomPID = data.GetPtrInt("mcMomPID");
    float *mcMomMass = data.GetPtrFloat("mcMomMass");
    int *mcGMomPID = data.GetPtrInt("mcGMomPID");
    float *mcPt = data.GetPtrFloat("mcPt");
    float *mcEta = data.GetPtrFloat("mcEta");
    float *mcPhi = data.GetPtrFloat("mcPhi");
    float *mcVtx = data.GetPtrFloat("mcVtx");
    float *mcVty = data.GetPtrFloat("mcVty");
    float *mcVtz = data.GetPtrFloat("mcVtz");
    UShort_t *mcStatusFlag = (UShort_t *)data.GetPtrShort("mcStatusFlag");
    // RECO information
    Int_t nEle = data.GetInt("nEle");
    int *eleCharge = data.GetPtrInt("eleCharge");
    int *eleChargeConsistent = data.GetPtrInt("eleChargeConsistent");
    float *eleEn = data.GetPtrFloat("eleEn");
    float *eleSCEn = data.GetPtrFloat("eleSCEn");
    float *eleEcalEn = data.GetPtrFloat("eleEcalEn");
    float *eleESEnP1 = data.GetPtrFloat("eleESEnP1");
    float *eleESEnP2 = data.GetPtrFloat("eleESEnP2");
    float *eleD0 = data.GetPtrFloat("eleD0");
    float *eleDz = data.GetPtrFloat("eleDz");
    float *eleSIP = data.GetPtrFloat("eleSIP");
    float *elePt = data.GetPtrFloat("elePt");
    float *elePtError = data.GetPtrFloat("elePtError");
    float *eleEta = data.GetPtrFloat("eleEta");
    float *elePhi = data.GetPtrFloat("elePhi");
    float *eleR9 = data.GetPtrFloat("eleR9");
    float *eleCalibPt = data.GetPtrFloat("eleCalibPt");
    float *eleCalibEn = data.GetPtrFloat("eleCalibEn");
    float *eleSCEta = data.GetPtrFloat("eleSCEta");
    float *eleSCPhi = data.GetPtrFloat("eleSCPhi");
    float *eleSCRawEn = data.GetPtrFloat("eleSCRawEn");
    float *eleSCEtaWidth = data.GetPtrFloat("eleSCEtaWidth");
    float *eleSCPhiWidth = data.GetPtrFloat("eleSCPhiWidth");
    float *eleHoverE = data.GetPtrFloat("eleHoverE");
    float *eleEoverP = data.GetPtrFloat("eleEoverP");
    float *eleEoverPout = data.GetPtrFloat("eleEoverPout");
    float *eleEoverPInv = data.GetPtrFloat("eleEoverPInv");
    float *eleBrem = data.GetPtrFloat("eleBrem");
    float *eledEtaAtVtx = data.GetPtrFloat("eledEtaAtVtx");
    float *eledPhiAtVtx = data.GetPtrFloat("eledPhiAtVtx");
    float *eleSigmaIEtaIEtaFull5x5 = data.GetPtrFloat("eleSigmaIEtaIEtaFull5x5");
    float *eleSigmaIPhiIPhiFull5x5 = data.GetPtrFloat("eleSigmaIPhiIPhiFull5x5");
    int *eleConvVeto = data.GetPtrInt("eleConvVeto");
    int *eleMissHits = data.GetPtrInt("eleMissHits");
    float *eleESEffSigmaRR = data.GetPtrFloat("eleESEffSigmaRR");
    float *elePFChIso = data.GetPtrFloat("elePFChIso");
    float *elePFPhoIso = data.GetPtrFloat("elePFPhoIso");
    float *elePFNeuIso = data.GetPtrFloat("elePFNeuIso");
    float *elePFPUIso = data.GetPtrFloat("elePFPUIso");
    float *elePFClusEcalIso = data.GetPtrFloat("elePFClusEcalIso");
    float *elePFClusHcalIso = data.GetPtrFloat("elePFClusHcalIso");
    float *eleIDMVAIso = data.GetPtrFloat("eleIDMVAIso");
    float *eleIDMVANoIso = data.GetPtrFloat("eleIDMVANoIso");
    float *eleR9Full5x5 = data.GetPtrFloat("eleR9Full5x5");
    int *eleEcalDrivenSeed = data.GetPtrInt("eleEcalDrivenSeed");
    float *eleTrkdxy = data.GetPtrFloat("eleTrkdxy");
    float *eleKFHits = data.GetPtrFloat("eleKFHits");
    float *eleKFChi2 = data.GetPtrFloat("eleKFChi2");
    float *eleGSFChi2 = data.GetPtrFloat("eleGSFChi2");
    ULong64_t *eleFiredSingleTrgs = (ULong64_t *)data.GetPtrLong64("eleFiredSingleTrgs");
    ULong64_t *eleFiredDoubleTrgs = (ULong64_t *)data.GetPtrLong64("eleFiredDoubleTrgs");
    ULong64_t *eleFiredL1Trgs = (ULong64_t *)data.GetPtrLong64("eleFiredL1Trgs");
    UShort_t *eleIDbit = (UShort_t *)data.GetPtrShort("eleIDbit");
    // Gsf track information
    Int_t nGSFTrk = data.GetInt("nGSFTrk");
    float *gsfPt = data.GetPtrFloat("gsfPt");
    float *gsfEta = data.GetPtrFloat("gsfEta");
    float *gsfPhi = data.GetPtrFloat("gsfPhi");
    int *gsfMissHits = data.GetPtrInt("gsfMissHits");

    // BC information 
    float *bcR15 = data.GetPtrFloat("bcR15");

    HLTEleMuX = event.HLTEleMuX;
    HLTEleMuXIsPrescaled = event.HLTEleMuXIsPrescaled;
    HLTPho = event.HLTPho;
    HLTPhoIsPrescaled = event.HLTPhoIsPrescaled;
    rho = event.rho;
    rhoCentral = event.rhoCentral;
    nVtx = event.nVtx;
    nGoodVtx = event.nGoodVtx;
    isPVGood = event.isPVGood;
    totalEvents = event.totalEvents;
    instwei = event.instwei;

    // Generator information
    procXS = event.procXS;
    mcwei = event.mcwei;
    genwei = event.genwei;
    mcPt_ = mcPt[event.iGen];
    mcEta_ = mcEta[event.iGen];
    mcPhi_ = mcPhi[event.iGen];
    mcVtx_ = mcVtx[event.iGen];
    mcVty_ = mcVty[event.iGen];
    mcVtz_ = mcVtz[event.iGen];
    mcMomPID_ = mcMomPID[event.iGen];
    mcGMomPID_ = mcGMomPID[event.iGen];
    mcStatusFlag_ = mcStatusFlag[event.iGen];

    elePresel = preselection(data, event.iReco);
    eleCharge_ = eleCharge[event.iReco];
    eleChargeConsistent_ = eleChargeConsistent[event.iReco];
    eleConvVeto_ = eleConvVeto[event.iReco];
    eleMissHits_ = eleMissHits[event.iReco];
    eleEcalDrivenSeed_ = eleEcalDrivenSeed[event.iReco];
    eleEn_ = eleEn[event.iReco];
    eleSCEn_ = eleSCEn[event.iReco];
    eleEcalEn_ = eleEcalEn[event.iReco];
    eleESEnP1_ = eleESEnP1[event.iReco];
    eleESEnP2_ = eleESEnP2[event.iReco];
    eleD0_ = eleD0[event.iReco];
    eleDz_ = eleDz[event.iReco];
    eleSIP_ = eleSIP[event.iReco];
    elePt_ = elePt[event.iReco];
    elePtError_ = elePtError[event.iReco];
    eleEta_ = eleEta[event.iReco];
    elePhi_ = elePhi[event.iReco];
    eleR9_ = eleR9[event.iReco];
    eleCalibPt_ = eleCalibPt[event.iReco];
    eleDiffCalibOriPt_ = (eleCalibPt[event.iReco] - elePt[event.iReco]) / elePt[event.iReco];
    eleCalibEn_ = eleCalibEn[event.iReco];
    eleSCEta_ = eleSCEta[event.iReco];
    eleSCPhi_ = eleSCPhi[event.iReco];
    eleSCRawEn_ = eleSCRawEn[event.iReco];
    eleSCEtaWidth_ = eleSCEtaWidth[event.iReco];
    eleSCPhiWidth_ = eleSCPhiWidth[event.iReco];
    eleHoverE_ = eleHoverE[event.iReco];
    eleEoverP_ = eleEoverP[event.iReco];
    eleEoverPout_ = eleEoverPout[event.iReco];
    eleEoverPInv_ = eleEoverPInv[event.iReco];
    eleBrem_ = eleBrem[event.iReco];
    eledEtaAtVtx_ = eledEtaAtVtx[event.iReco];
    eledPhiAtVtx_ = eledPhiAtVtx[event.iReco];
    eleSigmaIEtaIEtaFull5x5_ = eleSigmaIEtaIEtaFull5x5[event.iReco];
    eleSigmaIPhiIPhiFull5x5_ = eleSigmaIPhiIPhiFull5x5[event.iReco];
    eleESEffSigmaRR_ = eleESEffSigmaRR[event.iReco];
    elePFChIso_ = elePFChIso[event.iReco];
    elePFPhoIso_ = elePFPhoIso[event.iReco];
    elePFNeuIso_ = elePFNeuIso[event.iReco];
    elePFPUIso_ = elePFPUIso[event.iReco];
    elePFClusEcalIso_ = elePFClusEcalIso[event.iReco];
    elePFClusHcalIso_ = elePFClusHcalIso[event.iReco];
    eleIDMVAIso_ = eleIDMVAIso[event.iReco];
    eleIDMVANoIso_ = eleIDMVANoIso[event.iReco];
    eleR9Full5x5_ = eleR9Full5x5[event.iReco];
    eleTrkdxy_ = eleTrkdxy[event.iReco];
    eleKFHits_ = eleKFHits[event.iReco];
    eleKFChi2_ = eleKFChi2[event.iReco];
    eleGSFChi2_ = eleGSFChi2[event.iReco];
    eleFiredSingleTrgs_ = eleFiredSingleTrgs[event.iReco];
    eleFiredDoubleTrgs_ = eleFiredDoubleTrgs[event.iReco];
    eleFiredL1Trgs_ = eleFiredL1Trgs[event.iReco];
    eleIDbit_ = eleIDbit[event.iReco];

    eleESEnToRawE_ = (eleESEnP1[event.iReco] + eleESEnP2[event.iReco])/(eleSCRawEn[event.iReco]);

    GsfIdx = event.GsfIdx;
    if (event.GsfIdx.size() > 0)
    {
        for (size_t i = 0; i < event.GsfIdx.size(); i++)
        {
            gsfPt_reco.push_back(gsfPt[event.GsfIdx[i]]);
            gsfEta_reco.push_back(gsfEta[event.GsfIdx[i]]);
            gsfPhi_reco.push_back(gsfPhi[event.GsfIdx[i]]);
            gsfMissHits_reco.push_back(gsfMissHits[event.GsfIdx[i]]);
        }
    }
    else
    {
        gsfPt_reco.push_back(gsfPt[-1]);
        gsfEta_reco.push_back(gsfEta[-1]);
        gsfPhi_reco.push_back(gsfPhi[-1]);
        gsfMissHits_reco.push_back(gsfPhi[-1]);
    }

    nGsfMatchToReco = GsfIdx.size();

    // gsfPt_ = gsfPt[event.iGsf];
    // gsfEta_ = gsfEta[event.iGsf];
    // gsfPhi_ = gsfPhi[event.iGsf];

    if (event.GsfIdx.size() > 1)
    {
        // Float_t tmp_gsfPtSum = -1., tmp_gsfPtRatio = -1., tmp_dR = 999.;
        // Int_t tmp_gsfMissHitsSum = -1;
        // for (size_t i = 0; i < event.GsfIdx.size(); i++)
        // {
        //     TLorentzVector l1;
        //     l1.SetPtEtaPhiM(gsfPt[event.GsfIdx[i]], gsfEta[event.GsfIdx[i]], gsfPhi[event.GsfIdx[i]], 0.000510999);
        //     for (size_t j = i + 1; j < event.GsfIdx.size(); j++)
        //     {
        //         TLorentzVector l2;
        //         l2.SetPtEtaPhiM(gsfPt[event.GsfIdx[j]], gsfEta[event.GsfIdx[j]], gsfPhi[event.GsfIdx[j]], 0.000510999);
        //         if (l1.DeltaR(l2) < tmp_dR)
        //         {
        //             tmp_gsfPtSum = l1.Pt() + l2.Pt();
        //             tmp_gsfPtRatio = (l1.Pt() > l2.Pt()) ? (l2.Pt()/l1.Pt()) : (l1.Pt()/l2.Pt());
        //             tmp_dR = l1.DeltaR(l2);
        //             tmp_gsfMissHitsSum = gsfMissHits[event.GsfIdx[i]] + gsfMissHits[event.GsfIdx[j]];
        //         }
        //     }
        // }

        // gsfPtSum = tmp_gsfPtSum;
        // gsfPtRatio = tmp_gsfPtRatio;
        // gsfDeltaR = tmp_dR;
        // gsfMissHitsSum = tmp_gsfMissHitsSum;

        // two highest pT gsf trakcs (CH's version)
        float gsfPt_lead = 0.;
        int gsfPt_lead_ind = 0;
        for (size_t i = 0; i < event.GsfIdx.size(); i++){
            if(gsfPt_lead < gsfPt[event.GsfIdx[i]]) {
                gsfPt_lead = gsfPt[event.GsfIdx[i]];
                gsfPt_lead_ind = event.GsfIdx[i];
            }
        }

        float gsfPt_sublead = 0.;
        int gsfPt_sublead_ind = 0;
        for (size_t j = 0; j < event.GsfIdx.size(); j++){
            if (event.GsfIdx[j] == gsfPt_lead_ind) continue;
            if(gsfPt_sublead < gsfPt[event.GsfIdx[j]]) {
                gsfPt_sublead = gsfPt[event.GsfIdx[j]];
                gsfPt_sublead_ind = event.GsfIdx[j];
            }
        }
        
        TLorentzVector gsf1, gsf2, digsf, ele;
        gsf1.SetPtEtaPhiM(gsfPt[gsfPt_lead_ind], gsfEta[gsfPt_lead_ind], gsfPhi[gsfPt_lead_ind], 0.000511);
        gsf2.SetPtEtaPhiM(gsfPt[gsfPt_sublead_ind], gsfEta[gsfPt_sublead_ind], gsfPhi[gsfPt_sublead_ind], 0.000511);
        ele.SetPtEtaPhiE(eleCalibEn[event.iReco]/cosh(eleEta[event.iReco]), eleEta[event.iReco], elePhi[event.iReco], eleCalibEn[event.iReco]);
        digsf = gsf1 + gsf2;

        gsfPtSum = gsfPt_lead + gsfPt_sublead;
        gsfPtRatio = gsfPt_sublead / gsfPt_lead;
        gsfDeltaR = deltaR(gsfEta[gsfPt_lead_ind], gsfPhi[gsfPt_lead_ind], gsfEta[gsfPt_sublead_ind], gsfPhi[gsfPt_sublead_ind]);
        gsfMissHitsSum = gsfMissHits[gsfPt_lead_ind] + gsfMissHits[gsfPt_sublead_ind];

        gsfPtSum2 = (gsfPt_lead * gsfPt_lead) + (gsfPt_sublead * gsfPt_sublead);
        gsfRelPtRatio = digsf.Pt() / eleCalibPt[event.iReco];

        gsfPerpEleSum = (gsf1.Vect().Perp(ele.Vect()) + gsf2.Vect().Perp(ele.Vect()));

        gsfDiffPtRatio = (digsf.Pt() - eleCalibPt[event.iReco])/(digsf.Pt() + eleCalibPt[event.iReco]);
    }
    else if (event.GsfIdx.size() == 1){
        TLorentzVector gsf1, ele;
        gsf1.SetPtEtaPhiM(gsfPt[event.GsfIdx[0]], gsfEta[event.GsfIdx[0]], gsfPhi[event.GsfIdx[0]], 0.000511);
        ele.SetPtEtaPhiE(eleCalibEn[event.iReco]/cosh(eleEta[event.iReco]), eleEta[event.iReco], elePhi[event.iReco], eleCalibEn[event.iReco]);

        gsfPtSum = gsfPt[event.GsfIdx[0]];
        gsfPtRatio = -1.;
        gsfDeltaR = -999.;
        gsfMissHitsSum = gsfMissHits[event.GsfIdx[0]];

        gsfPtSum2 = (gsfPt[event.GsfIdx[0]] * gsfPt[event.GsfIdx[0]]);
        gsfRelPtRatio = gsfPt[event.GsfIdx[0]]/ eleCalibPt[event.iReco];
        gsfPerpEleSum = gsf1.Vect().Perp(ele.Vect());

        gsfDiffPtRatio = (gsfPt[event.GsfIdx[0]] - eleCalibPt[event.iReco])/(gsfPt[event.GsfIdx[0]] + eleCalibPt[event.iReco]);
    }
    else{ // no matched Gsf
        gsfPtSum = 0;
        gsfPtRatio = -1.;
        gsfDeltaR = -999.;
        gsfMissHitsSum = 0;

        gsfPtSum2 = 0;
        gsfRelPtRatio = -1;
        gsfPerpEleSum = -999;
        gsfDiffPtRatio = -999;
    }

    circularity = 1 - bcR15[event.BCIdx[0]];

    _tree->Fill();
}