#include <TLorentzVector.h>
#include <cmath>

class Minitree
{
    TFile *fo;
    TTree *_tree;

    // Category label: 1: Resolved, 2: Merged-2Gsf, 3: Merged-1 Missing track, 4: Not properly reconstructed
    Int_t category, GenReco_case, RecoGsf_case, GsfGen_case;
    ULong64_t HLTEleMuX, HLTEleMuXIsPrescaled, HLTPho, HLTPhoIsPrescaled;
    Float_t rho, rhoCentral;
    Int_t nVtx, nGoodVtx;
    Bool_t isPVGood;
    // Generator-level information
    Float_t mcPt_lep1, mcEta_lep1, mcPhi_lep1, mcVtx_lep1, mcVty_lep1, mcVtz_lep1, mcPt_lep2, mcEta_lep2, mcPhi_lep2, mcVtx_lep2, mcVty_lep2, mcVtz_lep2, mcPt_pho, mcEta_pho, mcPhi_pho;
    Float_t mcPt_ll, mcEta_ll, mcPhi_ll, mc_Mll, mcPt_llg, mcEta_llg, mcPhi_llg, mc_Mllg;
    Int_t mcMomPID_lep1, mcGMomPID_lep1, mcMomPID_lep2, mcGMomPID_lep2;
    Int_t mcMomPID_pho, mcGMomPID_pho;
    UShort_t mcStatusFlag_lep1, mcStatusFlag_lep2, mcStatusFlag_pho;
    Float_t mcwei, genwei;
    Float_t procXS;
    Long64_t totalEvents;
    Float_t instwei;
    // Reco-level information
    Bool_t elePresel_lep1, elePresel_lep2;
    Int_t eleCharge_lep1, eleCharge_lep2, eleChargeConsistent_lep1, eleChargeConsistent_lep2, eleConvVeto_lep1, eleConvVeto_lep2, eleMissHits_lep1, eleMissHits_lep2, eleEcalDrivenSeed_lep1, eleEcalDrivenSeed_lep2;
    Float_t eleEn_lep1, eleEn_lep2, eleSCEn_lep1, eleSCEn_lep2, eleEcalEn_lep1, eleEcalEn_lep2, eleESEnP1_lep1, eleESEnP1_lep2, eleESEnP2_lep1, eleESEnP2_lep2, eleD0_lep1, eleD0_lep2, eleDz_lep1, eleDz_lep2, eleSIP_lep1, eleSIP_lep2, elePt_lep1, elePt_lep2, elePtError_lep1, elePtError_lep2, eleEta_lep1, eleEta_lep2, elePhi_lep1, elePhi_lep2, eleR9_lep1, eleR9_lep2, eleCalibPt_lep1, eleCalibPt_lep2, eleDiffCalibOriPt_lep1, eleDiffCalibOriPt_lep2, eleCalibEn_lep1, eleCalibEn_lep2, eleSCEta_lep1, eleSCEta_lep2, eleSCPhi_lep1, eleSCPhi_lep2, eleSCRawEn_lep1, eleSCRawEn_lep2, eleSCEtaWidth_lep1, eleSCEtaWidth_lep2, eleSCPhiWidth_lep1, eleSCPhiWidth_lep2, eleHoverE_lep1, eleHoverE_lep2, eleEoverP_lep1, eleEoverP_lep2, eleEoverPout_lep1, eleEoverPout_lep2, eleEoverPInv_lep1, eleEoverPInv_lep2, eleBrem_lep1, eleBrem_lep2, eledEtaAtVtx_lep1, eledEtaAtVtx_lep2, eledPhiAtVtx_lep1, eledPhiAtVtx_lep2, eleSigmaIEtaIEtaFull5x5_lep1, eleSigmaIEtaIEtaFull5x5_lep2, eleSigmaIPhiIPhiFull5x5_lep1, eleSigmaIPhiIPhiFull5x5_lep2, eleESEffSigmaRR_lep1, eleESEffSigmaRR_lep2, elePFChIso_lep1, elePFChIso_lep2, elePFPhoIso_lep1, elePFPhoIso_lep2, elePFNeuIso_lep1, elePFNeuIso_lep2, elePFPUIso_lep1, elePFPUIso_lep2, elePFClusEcalIso_lep1, elePFClusEcalIso_lep2, elePFClusHcalIso_lep1, elePFClusHcalIso_lep2, eleIDMVAIso_lep1, eleIDMVAIso_lep2, eleIDMVANoIso_lep1, eleIDMVANoIso_lep2, eleR9Full5x5_lep1, eleR9Full5x5_lep2, eleTrkdxy_lep1, eleTrkdxy_lep2, eleKFHits_lep1, eleKFHits_lep2, eleKFChi2_lep1, eleKFChi2_lep2, eleGSFChi2_lep1, eleGSFChi2_lep2;
    ULong64_t eleFiredSingleTrgs_lep1, eleFiredSingleTrgs_lep2, eleFiredDoubleTrgs_lep1, eleFiredDoubleTrgs_lep2, eleFiredL1Trgs_lep1, eleFiredL1Trgs_lep2;
    UShort_t eleIDbit_lep1, eleIDbit_lep2;
    // Gsf track information
    vector<Int_t> GsfIdxLep1, GsfIdxLep2, gsfMissHits_reco1, gsfMissHits_reco2; //RECO-Gsf matching
    vector<Float_t> gsfPt_reco1, gsfPt_reco2, gsfEta_reco1, gsfEta_reco2, gsfPhi_reco1, gsfPhi_reco2; //RECO-Gsf matching
    Float_t gsfPt_Lep1, gsfPt_Lep2, gsfEta_Lep1, gsfEta_Lep2, gsfPhi_Lep1, gsfPhi_Lep2; // Gen-Gsf matching
    Int_t nGsfMatchToReco_lep1, nGsfMatchToReco_lep2; //RECO-Gsf matching
    Float_t gsfPtSum_lep1, gsfPtRatio_lep1, gsfDeltaR_lep1, gsfPtSum_lep2, gsfPtRatio_lep2, gsfDeltaR_lep2; //RECO-Gsf matching; Stored only when resolved and Merged-2Gsf
    Int_t gsfMissHitsSum_lep1, gsfMissHitsSum_lep2;
    
    // BC information
    Float_t circularity_lep1, circularity_lep2;

    // Conv information
    Int_t nConv_ = 0;
    Bool_t convMatched_lep1 = 0, convMatched_lep2 = 0;
    Float_t convVtxRadius_lep1 = 0, convVtxRadius_lep2 = 0;
    Float_t convP_lep1 = 0, convP_lep2 = 0;
    Int_t convNTrks_lep1 = 0, convNTrks_lep2 = 0;
    Float_t convD0_lep1 = 0, convD0_lep2 = 0, convDz_lep1 = 0, convDz_lep2 = 0, convL0_lep1 = 0, convL0_lep2 = 0, convLz_lep1 = 0, convLz_lep2 = 0;

    // preshower variables
    float eleESEnToRawE_lep1, eleESEnToRawE_lep2;

    float gsfPtSum2_lep1, gsfPtSum2_lep2;
    float gsfRelPtRatio_lep1, gsfRelPtRatio_lep2;
    float gsfPerpEleSum_lep1, gsfPerpEleSum_lep2;
    float gsfDiffPtRatio_lep1, gsfDiffPtRatio_lep2;

    vector<const char *> vhname_category{"Total", "Resolved", "Merged2Gsf", "Merged1MissingGsf", "NPR"};
    vector<TH1F *> v_hM_Mll;
    vector<TH1F *> v_hM_Mll_30to40;
    vector<TH1F *> v_hM_dR;
    vector<TH1F *> v_hM_Ptll;
    vector<TH1F *> v_hM_LeadLepPt;
    vector<TH1F *> v_hM_TrailLepPt;
    vector<TH1F *> v_hM_Ptllg;
    vector<TH1F *> v_hM_sumPtll;
    vector<const char *> vhname_category_showershape{"Resolved", "Resolved-leading", "Resolved-trailing", "Merged2Gsf", "Merged1MissingGsf"};
    vector<TH1F *> v_hM_eleEn;
    vector<TH1F *> v_hM_eleSCEn;
    vector<TH1F *> v_hM_eleEcalEn;
    vector<TH1F *> v_hM_eleESEnP1;
    vector<TH1F *> v_hM_eleESEnP2;
    vector<TH1F *> v_hM_eleD0;
    vector<TH1F *> v_hM_eleDz;
    vector<TH1F *> v_hM_eleSIP;
    vector<TH1F *> v_hM_elePt;
    vector<TH1F *> v_hM_elePtError;
    vector<TH1F *> v_hM_eleEta;
    vector<TH1F *> v_hM_elePhi;
    vector<TH1F *> v_hM_eleR9;
    vector<TH1F *> v_hM_eleCalibPt;
    vector<TH1F *> v_hM_diffCalibOriPt;
    vector<TH1F *> v_hM_eleCalibEn;
    vector<TH1F *> v_hM_eleSCEta;
    vector<TH1F *> v_hM_eleSCPhi;
    vector<TH1F *> v_hM_eleSCRawEn;
    vector<TH1F *> v_hM_eleSCEtaWidth;
    vector<TH1F *> v_hM_eleSCPhiWidth;
    vector<TH1F *> v_hM_eleHoverE;
    vector<TH1F *> v_hM_eleEoverP;
    vector<TH1F *> v_hM_eleEoverPout;
    vector<TH1F *> v_hM_eleEoverPInv;
    vector<TH1F *> v_hM_eleBrem;
    vector<TH1F *> v_hM_eledEtaAtVtx;
    vector<TH1F *> v_hM_eledPhiAtVtx;
    vector<TH1F *> v_hM_eleSigmaIEtaIEtaFull5x5;
    vector<TH1F *> v_hM_eleSigmaIPhiIPhiFull5x5;
    vector<TH1F *> v_hM_eleSigmaIEtaIEtaFull5x5_EB;
    vector<TH1F *> v_hM_eleSigmaIPhiIPhiFull5x5_EB;
    vector<TH1F *> v_hM_eleSigmaIEtaIEtaFull5x5_EE;
    vector<TH1F *> v_hM_eleSigmaIPhiIPhiFull5x5_EE;
    vector<TH1F *> v_hM_eleESEffSigmaRR;
    vector<TH1F *> v_hM_elePFChIso;
    vector<TH1F *> v_hM_elePFPhoIso;
    vector<TH1F *> v_hM_elePFNeuIso;
    vector<TH1F *> v_hM_elePFPUIso;
    vector<TH1F *> v_hM_elePFClusEcalIso;
    vector<TH1F *> v_hM_elePFClusHcalIso;
    vector<TH1F *> v_hM_eleIDMVAIso;
    vector<TH1F *> v_hM_eleIDMVANoIso;
    vector<TH1F *> v_hM_eleR9Full5x5;
    vector<TH1F *> v_hM_eleTrkdxy;
    vector<TH1F *> v_hM_eleKFHits;
    vector<TH1F *> v_hM_eleKFChi2;
    vector<TH1F *> v_hM_eleGSFChi2;

    TH1F *hM_dRSCs;
    TH1F *hM_dRSCs_SCEtaPhi;
    TH2F *hM_dRSCs_OrdSC;
    TH2F *hM_dEta_OrdSC;
    TH2F *hM_dPhi_OrdSC;
    TH2F *hM_corr_step12;
    TH2F *hM_corr_step23;
    TH1F *hM_eleGsfPtSum_resolved;
    TH1F *hM_eleGsfPtSum_merged2Gsf;
    TH2F *hM_ele1GsfPt_ele2GsfPt_resolved;
    TH2F *hM_ele1GsfPt_ele2GsfPt_merged2Gsf;

private:
    vector<TH1F *> Book1DHist(const char *basename, vector<const char *> vname, Int_t Nbins, Float_t minval, Float_t maxval)
    {
        vector<TH1F *> vhist;
        vhist.clear();
        for (size_t i = 0; i < vname.size(); i++)
        {
            TH1F *hM = new TH1F(Form("%s_%s", basename, vname[i]), Form("%s_%s", basename, vname[i]), Nbins, minval, maxval);
            vhist.push_back(hM);
        }

        return vhist;
    };
    void Fill1DHist(vector<TH1F *> &vhist, Int_t category, Float_t val, Float_t evtweight)
    {
        vhist[0]->Fill(val, evtweight);
        vhist[category]->Fill(val, evtweight);
    };
    void Fill1DHist_SS(vector<TH1F *> &vhist, Int_t category, Float_t val_ele1, Float_t val_ele2, Float_t evtweight) // For shower shape
    {
        // Resolved
        if (category == 1)
        {
            vhist[0]->Fill(val_ele1, evtweight); // Leading + Trailing
            vhist[0]->Fill(val_ele2, evtweight); // Leading + Trailing
            vhist[1]->Fill(val_ele1, evtweight); // Leading
            vhist[2]->Fill(val_ele2, evtweight); // Trailing
        }
        // Merged2Gsf & Merged1MissingGsf -> val_ele1 = val_ele2
        if (category == 2)
            vhist[3]->Fill(val_ele1, evtweight);
        if (category == 3)
            vhist[4]->Fill(val_ele1, evtweight);
    };
    void Fill1DHist_SS_EBEE(vector<TH1F *> &vhist_EB, vector<TH1F *> &vhist_EE, Int_t category, Float_t SCEta_ele1, Float_t SCEta_ele2, Float_t val_ele1, Float_t val_ele2, Float_t evtweight) // For shower shape, to saperate EB and EE
    {
        // Resolved
        if (category == 1)
        {
            if (fabs(SCEta_ele1) < 1.4442)
            {
                vhist_EB[0]->Fill(val_ele1, evtweight);
                vhist_EB[1]->Fill(val_ele1, evtweight);
            }
            if (fabs(SCEta_ele1) > 1.566 && fabs(SCEta_ele1) < 2.5)
            {
                vhist_EE[0]->Fill(val_ele1, evtweight);
                vhist_EE[1]->Fill(val_ele1, evtweight);
            }
            if (fabs(SCEta_ele2) < 1.4442)
            {
                vhist_EB[0]->Fill(val_ele2, evtweight);
                vhist_EB[2]->Fill(val_ele2, evtweight);
            }
            if (fabs(SCEta_ele2) > 1.566 && fabs(SCEta_ele2) < 2.5)
            {
                vhist_EE[0]->Fill(val_ele2, evtweight);
                vhist_EE[2]->Fill(val_ele2, evtweight);
            }
        }
        // Merged2Gsf & Merged1MissingGsf -> val_ele1 = val_ele2
        if (category == 2)
        {
            if (fabs(SCEta_ele1) < 1.4442)
                vhist_EB[3]->Fill(val_ele1, evtweight);
            if (fabs(SCEta_ele1) > 1.566 && fabs(SCEta_ele1) < 2.5)
                vhist_EE[3]->Fill(val_ele1, evtweight);
        }
        if (category == 3)
        {
            if (fabs(SCEta_ele1) < 1.4442)
                vhist_EB[4]->Fill(val_ele1, evtweight);
            if (fabs(SCEta_ele1) > 1.566 && fabs(SCEta_ele1) < 2.5)
                vhist_EE[4]->Fill(val_ele1, evtweight);
        }
    };
    Float_t deltaPhi(Float_t phi1, Float_t phi2)
    {

        Float_t dPhi = phi1 - phi2;
        if (dPhi > TMath::Pi())
            dPhi -= 2. * TMath::Pi();
        if (dPhi < -TMath::Pi())
            dPhi += 2. * TMath::Pi();

        return dPhi;
    }

public:
    Minitree(){};
    Minitree(const char *outpath)
    {
        fo = TFile::Open(outpath, "RECREATE");
        cout << "outpath [" << outpath << "] is opened." << endl;
        fo->cd();
    };

    ~Minitree()
    {
        _tree->Write("", TObject::kOverwrite);
        for (size_t i = 0; i < vhname_category.size(); i++)
        {
            v_hM_Mll[i]->Write();
            v_hM_Mll_30to40[i]->Write();
            v_hM_dR[i]->Write();
            v_hM_Ptll[i]->Write();
            v_hM_LeadLepPt[i]->Write();
            v_hM_TrailLepPt[i]->Write();
            v_hM_Ptllg[i]->Write();
            v_hM_sumPtll[i]->Write();
        }

        for (size_t i = 0; i < vhname_category_showershape.size(); i++)
        {
            v_hM_eleEn[i]->Write();
            v_hM_eleSCEn[i]->Write();
            v_hM_eleEcalEn[i]->Write();
            v_hM_eleESEnP1[i]->Write();
            v_hM_eleESEnP2[i]->Write();
            v_hM_eleD0[i]->Write();
            v_hM_eleDz[i]->Write();
            v_hM_eleSIP[i]->Write();
            v_hM_elePt[i]->Write();
            v_hM_elePtError[i]->Write();
            v_hM_eleEta[i]->Write();
            v_hM_elePhi[i]->Write();
            v_hM_eleR9[i]->Write();
            v_hM_eleCalibPt[i]->Write();
            v_hM_diffCalibOriPt[i]->Write();
            v_hM_eleCalibEn[i]->Write();
            v_hM_eleSCEta[i]->Write();
            v_hM_eleSCPhi[i]->Write();
            v_hM_eleSCRawEn[i]->Write();
            v_hM_eleSCEtaWidth[i]->Write();
            v_hM_eleSCPhiWidth[i]->Write();
            v_hM_eleHoverE[i]->Write();
            v_hM_eleEoverP[i]->Write();
            v_hM_eleEoverPout[i]->Write();
            v_hM_eleEoverPInv[i]->Write();
            v_hM_eleBrem[i]->Write();
            v_hM_eledEtaAtVtx[i]->Write();
            v_hM_eledPhiAtVtx[i]->Write();
            v_hM_eleSigmaIEtaIEtaFull5x5[i]->Write();
            v_hM_eleSigmaIPhiIPhiFull5x5[i]->Write();
            v_hM_eleSigmaIEtaIEtaFull5x5[i]->Write();
            v_hM_eleSigmaIPhiIPhiFull5x5[i]->Write();
            v_hM_eleSigmaIEtaIEtaFull5x5_EB[i]->Write();
            v_hM_eleSigmaIPhiIPhiFull5x5_EB[i]->Write();
            v_hM_eleSigmaIEtaIEtaFull5x5_EE[i]->Write();
            v_hM_eleSigmaIPhiIPhiFull5x5_EE[i]->Write();
            v_hM_eleESEffSigmaRR[i]->Write();
            v_hM_elePFChIso[i]->Write();
            v_hM_elePFPhoIso[i]->Write();
            v_hM_elePFNeuIso[i]->Write();
            v_hM_elePFPUIso[i]->Write();
            v_hM_elePFClusEcalIso[i]->Write();
            v_hM_elePFClusHcalIso[i]->Write();
            v_hM_eleIDMVAIso[i]->Write();
            v_hM_eleIDMVANoIso[i]->Write();
            v_hM_eleR9Full5x5[i]->Write();
            v_hM_eleTrkdxy[i]->Write();
            v_hM_eleKFHits[i]->Write();
            v_hM_eleKFChi2[i]->Write();
            v_hM_eleGSFChi2[i]->Write();
        }

        hM_dRSCs->Write();
        hM_dRSCs_SCEtaPhi->Write();
        hM_dRSCs_OrdSC->Write();
        hM_dEta_OrdSC->Write();
        hM_dPhi_OrdSC->Write();
        hM_corr_step12->Write();
        hM_corr_step23->Write();
        hM_eleGsfPtSum_resolved->Write();
        hM_eleGsfPtSum_merged2Gsf->Write();
        hM_ele1GsfPt_ele2GsfPt_resolved->Write();
        hM_ele1GsfPt_ele2GsfPt_merged2Gsf->Write();

        delete _tree;
        delete fo;
    };

    void SetHist();
    void Refresh();
    void SetMinitree(TString treename);
    void FillMinitree(TreeReader &data, Event &event);
};

void Minitree::SetHist()
{
    v_hM_Mll = Book1DHist("hM_Mll", vhname_category, 25, 0, 50);
    v_hM_Mll_30to40 = Book1DHist("hM_Mll_30to40", vhname_category, 20, 30, 40);
    v_hM_dR = Book1DHist("hM_dR", vhname_category, 30, 0, 3);
    v_hM_Ptll = Book1DHist("hM_Ptll", vhname_category, 50, 0, 200);
    v_hM_LeadLepPt = Book1DHist("hM_LeadLepPt", vhname_category, 50, 0, 200);
    v_hM_TrailLepPt = Book1DHist("hM_TrailLepPt", vhname_category, 25, 0, 100);
    v_hM_Ptllg = Book1DHist("hM_Ptllg", vhname_category, 50, 0, 200);
    v_hM_sumPtll = Book1DHist("v_hM_sumPtll", vhname_category, 50, 0, 200);
    // RECO object: Shower shape
    // See variable definitions in https://cmsdoxygen.web.cern.ch/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
    v_hM_eleEn = Book1DHist("hM_eleEn", vhname_category_showershape, 50, 0, 200);
    v_hM_eleSCEn = Book1DHist("hM_eleSCEn", vhname_category_showershape, 50, 0, 200);
    v_hM_eleEcalEn = Book1DHist("hM_eleEcalEn", vhname_category_showershape, 50, 0, 200);
    v_hM_eleESEnP1 = Book1DHist("hM_eleESEnP1", vhname_category_showershape, 60, 0, 60);
    v_hM_eleESEnP2 = Book1DHist("hM_eleESEnP2", vhname_category_showershape, 60, 0, 60);
    v_hM_eleD0 = Book1DHist("hM_eleD0", vhname_category_showershape, 60, 0, 3);
    v_hM_eleDz = Book1DHist("hM_eleDz", vhname_category_showershape, 50, 0, 5);
    v_hM_eleSIP = Book1DHist("hM_eleSIP", vhname_category_showershape, 50, 0, 10);
    v_hM_elePt = Book1DHist("hM_elePt", vhname_category_showershape, 50, 0, 200);
    v_hM_elePtError = Book1DHist("hM_elePtError", vhname_category_showershape, 50, 0, 20);
    v_hM_eleEta = Book1DHist("hM_eleEta", vhname_category_showershape, 60, -3, 3);
    v_hM_elePhi = Book1DHist("hM_elePhi", vhname_category_showershape, 70, -3.5, 3.5);
    v_hM_eleR9 = Book1DHist("hM_eleR9", vhname_category_showershape, 50, 0, 1);
    v_hM_eleCalibPt = Book1DHist("hM_eleCalibPt", vhname_category_showershape, 50, 0, 200);
    v_hM_diffCalibOriPt = Book1DHist("hM_diffCalibOriPt", vhname_category_showershape, 60, -0.6, 0.6);
    v_hM_eleCalibEn = Book1DHist("hM_eleCalibEn", vhname_category_showershape, 50, 0, 200);
    v_hM_eleSCEta = Book1DHist("hM_eleSCEta", vhname_category_showershape, 60, -3, 3);
    v_hM_eleSCPhi = Book1DHist("hM_eleSCPhi", vhname_category_showershape, 70, -3.5, 3.5);
    v_hM_eleSCRawEn = Book1DHist("hM_eleSCRawEn", vhname_category_showershape, 50, 0, 200);
    v_hM_eleSCEtaWidth = Book1DHist("hM_eleSCEtaWidth", vhname_category_showershape, 50, 0, 0.2); // https://github.com/cms-sw/cmssw/blob/master/RecoEcal/EgammaCoreTools/src/SuperClusterShapeAlgo.cc
    v_hM_eleSCPhiWidth = Book1DHist("hM_eleSCPhiWidth", vhname_category_showershape, 50, 0, 0.5); // https://github.com/cms-sw/cmssw/blob/master/RecoEcal/EgammaCoreTools/src/SuperClusterShapeAlgo.cc
    v_hM_eleHoverE = Book1DHist("hM_eleHoverE", vhname_category_showershape, 50, 0, 5);
    v_hM_eleEoverP = Book1DHist("hM_eleEoverP", vhname_category_showershape, 50, 0, 5);       // the supercluster energy / track momentum at the PCA to the beam spot
    v_hM_eleEoverPout = Book1DHist("hM_eleEoverPout", vhname_category_showershape, 50, 0, 2); // the electron cluster energy / track momentum at calo extrapolated from the outermost track state
    v_hM_eleEoverPInv = Book1DHist("hM_eleEoverPInv", vhname_category_showershape, 50, -1, 1);
    v_hM_eleBrem = Book1DHist("hM_eleBrem", vhname_category_showershape, 50, 0, 1);                 // the brem fraction from gsf fit: (track momentum in - track momentum out) / track momentum in
    v_hM_eledEtaAtVtx = Book1DHist("hM_eledEtaAtVtx", vhname_category_showershape, 120, -0.3, 0.3); // the supercluster eta - track eta position at calo extrapolated from the innermost track state
    v_hM_eledPhiAtVtx = Book1DHist("hM_eledPhiAtVtx", vhname_category_showershape, 120, -0.3, 0.3); // the supercluster phi - track phi position at calo extrapolated from the innermost track state
    v_hM_eleSigmaIEtaIEtaFull5x5 = Book1DHist("hM_eleSigmaIEtaIEtaFull5x5", vhname_category_showershape, 50, 0, 0.1);
    v_hM_eleSigmaIPhiIPhiFull5x5 = Book1DHist("hM_eleSigmaIPhiIPhiFull5x5", vhname_category_showershape, 50, 0, 0.1);
    v_hM_eleSigmaIEtaIEtaFull5x5_EB = Book1DHist("hM_eleSigmaIEtaIEtaFull5x5_EB", vhname_category_showershape, 50, 0, 0.1);
    v_hM_eleSigmaIPhiIPhiFull5x5_EB = Book1DHist("hM_eleSigmaIPhiIPhiFull5x5_EB", vhname_category_showershape, 50, 0, 0.1);
    v_hM_eleSigmaIEtaIEtaFull5x5_EE = Book1DHist("hM_eleSigmaIEtaIEtaFull5x5_EE", vhname_category_showershape, 50, 0, 0.1);
    v_hM_eleSigmaIPhiIPhiFull5x5_EE = Book1DHist("hM_eleSigmaIPhiIPhiFull5x5_EE", vhname_category_showershape, 50, 0, 0.1);
    v_hM_eleESEffSigmaRR = Book1DHist("hM_eleESEffSigmaRR", vhname_category_showershape, 50, 0, 10); // Preshower effective sigmaRR
    v_hM_elePFChIso = Book1DHist("hM_elePFChIso", vhname_category_showershape, 80, 0, 20);
    v_hM_elePFPhoIso = Book1DHist("hM_elePFPhoIso", vhname_category_showershape, 80, 0, 20);
    v_hM_elePFNeuIso = Book1DHist("hM_elePFNeuIso", vhname_category_showershape, 80, 0, 20);
    v_hM_elePFPUIso = Book1DHist("hM_elePFPUIso", vhname_category_showershape, 80, 0, 20);
    v_hM_elePFClusEcalIso = Book1DHist("hM_elePFClusEcalIso", vhname_category_showershape, 50, 0, 0.5);
    v_hM_elePFClusHcalIso = Book1DHist("hM_elePFClusHcalIso", vhname_category_showershape, 50, 0, 1.5);
    v_hM_eleIDMVAIso = Book1DHist("hM_eleIDMVAIso", vhname_category_showershape, 50, -1, 1);
    v_hM_eleIDMVANoIso = Book1DHist("hM_eleIDMVANoIso", vhname_category_showershape, 50, -1, 1);
    v_hM_eleR9Full5x5 = Book1DHist("hM_eleR9Full5x5", vhname_category_showershape, 50, 0, 1);
    v_hM_eleTrkdxy = Book1DHist("hM_eleTrkdxy", vhname_category_showershape, 50, 0, 1);
    v_hM_eleKFHits = Book1DHist("hM_eleKFHits", vhname_category_showershape, 50, 0, 1);
    v_hM_eleKFChi2 = Book1DHist("hM_eleKFChi2", vhname_category_showershape, 50, 0, 10);
    v_hM_eleGSFChi2 = Book1DHist("hM_eleGSFChi2", vhname_category_showershape, 50, 0, 10);

    hM_dRSCs = new TH1F("hM_dRSCs", "hM_dRSCs", 50, 0, 1);
    hM_dRSCs_SCEtaPhi = new TH1F("hM_dRSCs_SCEtaPhi", "hM_dRSCs_SCEtaPhi", 50, 0, 1);
    hM_dRSCs_OrdSC = new TH2F("hM_dRSCs_OrdSC", "hM_dRSCs_OrdSC", 200, 0, 4, 200, 0, 4);
    hM_dEta_OrdSC = new TH2F("hM_dEta_OrdSC", "hM_dEta_OrdSC", 200, 0, 4, 200, 0, 4);
    hM_dPhi_OrdSC = new TH2F("hM_dPhi_OrdSC", "hM_dPhi_OrdSC", 175, 0, 3.5, 175, 0, 3.5);
    hM_corr_step12 = new TH2F("hM_corr_step12", "hM_corr_step12", 4, 0.5, 4.5, 4, 0.5, 4.5); //GenToReco vs RecoToGsf
    hM_corr_step23 = new TH2F("hM_corr_step23", "hM_corr_step23", 4, 0.5, 4.5, 3, 0.5, 3.5); //RecoToGsf vs GsfToGen

    hM_eleGsfPtSum_resolved = new TH1F("hM_eleGsfPtSum_resolved", "hM_eleGsfPtSum_resolved", 50, 0, 200);
    hM_eleGsfPtSum_merged2Gsf = new TH1F("hM_eleGsfPtSum_merged2Gsf", "hM_eleGsfPtSum_merged2Gsf", 50, 0, 200);
    hM_ele1GsfPt_ele2GsfPt_resolved = new TH2F("hM_ele1GsfPt_ele2GsfPt_resolved", "hM_ele1GsfPt_ele2GsfPt_resolved", 50, 0, 100, 50, 0, 100);
    hM_ele1GsfPt_ele2GsfPt_merged2Gsf = new TH2F("hM_ele1GsfPt_ele2GsfPt_merged2Gsf", "hM_ele1GsfPt_ele2GsfPt_merged2Gsf", 50, 0, 100, 50, 0, 100);
}

void Minitree::Refresh()
{
    GsfIdxLep1.clear();
    GsfIdxLep2.clear();
    gsfPt_reco1.clear();
    gsfPt_reco2.clear();
    gsfEta_reco1.clear();
    gsfEta_reco2.clear();
    gsfPhi_reco1.clear();
    gsfPhi_reco2.clear();
    gsfMissHits_reco1.clear();
    gsfMissHits_reco2.clear();
}

void Minitree::SetMinitree(TString treename)
{
    _tree = new TTree(treename, "minitree");

    // Category
    _tree->Branch("category", &category);
    _tree->Branch("GenReco_case", &GenReco_case);
    _tree->Branch("RecoGsf_case", &RecoGsf_case);
    _tree->Branch("GsfGen_case", &GsfGen_case);
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
    _tree->Branch("mcPt_lep1", &mcPt_lep1);
    _tree->Branch("mcEta_lep1", &mcEta_lep1);
    _tree->Branch("mcPhi_lep1", &mcPhi_lep1);
    _tree->Branch("mcVtx_lep1", &mcVtx_lep1);
    _tree->Branch("mcVty_lep1", &mcVty_lep1);
    _tree->Branch("mctz_lep1", &mcVtz_lep1);
    _tree->Branch("mcPt_lep2", &mcPt_lep2);
    _tree->Branch("mcEta_lep2", &mcEta_lep2);
    _tree->Branch("mcPhi_lep2", &mcPhi_lep2);
    _tree->Branch("mcVtx_lep2", &mcVtx_lep2);
    _tree->Branch("mcVty_lep2", &mcVty_lep2);
    _tree->Branch("mcVtz_lep2", &mcVtz_lep2);
    _tree->Branch("mcPt_pho", &mcPt_pho);
    _tree->Branch("mcEta_pho", &mcEta_pho);
    _tree->Branch("mcPhi_pho", &mcPhi_pho);
    _tree->Branch("mcPt_ll", &mcPt_ll);
    _tree->Branch("mcEta_ll", &mcEta_ll);
    _tree->Branch("mcPhi_ll", &mcPhi_ll);
    _tree->Branch("mc_Mll", &mc_Mll);
    _tree->Branch("mcPt_llg", &mcPt_llg);
    _tree->Branch("mcEta_llg", &mcEta_llg);
    _tree->Branch("mcPhi_llg", &mcPhi_llg);
    _tree->Branch("mc_Mllg", &mc_Mllg);
    _tree->Branch("mcMomPID_lep1", &mcMomPID_lep1);
    _tree->Branch("mcGMomPID_lep1", &mcGMomPID_lep1);
    _tree->Branch("mcMomPID_lep2", &mcMomPID_lep2);
    _tree->Branch("mcGMomPID_lep2", &mcGMomPID_lep2);
    _tree->Branch("mcMomPID_pho", &mcMomPID_pho);
    _tree->Branch("mcGMomPID_pho", &mcGMomPID_pho);
    _tree->Branch("mcStatusFlag_lep1", &mcStatusFlag_lep1);
    _tree->Branch("mcStatusFlag_lep2", &mcStatusFlag_lep2);
    _tree->Branch("mcStatusFlag_pho", &mcStatusFlag_pho);
    // Reco information
    _tree->Branch("elePresel_lep1", &elePresel_lep1);
    _tree->Branch("elePresel_lep2", &elePresel_lep2);
    _tree->Branch("eleCharge_lep1", &eleCharge_lep1);
    _tree->Branch("eleCharge_lep2", &eleCharge_lep2);
    _tree->Branch("eleChargeConsistent_lep1", &eleChargeConsistent_lep1);
    _tree->Branch("eleChargeConsistent_lep2", &eleChargeConsistent_lep2);
    _tree->Branch("eleConvVeto_lep1", &eleConvVeto_lep1);
    _tree->Branch("eleConvVeto_lep2", &eleConvVeto_lep2);
    _tree->Branch("eleMissHits_lep1", &eleMissHits_lep1);
    _tree->Branch("eleMissHits_lep2", &eleMissHits_lep2);
    _tree->Branch("eleEcalDrivenSeed_lep1", &eleEcalDrivenSeed_lep1);
    _tree->Branch("eleEcalDrivenSeed_lep2", &eleEcalDrivenSeed_lep2);
    _tree->Branch("eleEn_lep1", &eleEn_lep1);
    _tree->Branch("eleEn_lep2", &eleEn_lep2);
    _tree->Branch("eleSCEn_lep1", &eleSCEn_lep1);
    _tree->Branch("eleSCEn_lep2", &eleSCEn_lep2);
    _tree->Branch("eleEcalEn_lep1", &eleEcalEn_lep1);
    _tree->Branch("eleEcalEn_lep2", &eleEcalEn_lep2);
    _tree->Branch("eleESEnP1_lep1", &eleESEnP1_lep1);
    _tree->Branch("eleESEnP1_lep2", &eleESEnP1_lep2);
    _tree->Branch("eleESEnP2_lep1", &eleESEnP2_lep1);
    _tree->Branch("eleESEnP2_lep2", &eleESEnP2_lep2);
    _tree->Branch("eleD0_lep1", &eleD0_lep1);
    _tree->Branch("eleD0_lep2", &eleD0_lep2);
    _tree->Branch("eleDz_lep1", &eleDz_lep1);
    _tree->Branch("eleDz_lep2", &eleDz_lep2);
    _tree->Branch("eleSIP_lep1", &eleSIP_lep1);
    _tree->Branch("eleSIP_lep2", &eleSIP_lep2);
    _tree->Branch("elePt_lep1", &elePt_lep1);
    _tree->Branch("elePt_lep2", &elePt_lep2);
    _tree->Branch("elePtError_lep1", &elePtError_lep1);
    _tree->Branch("elePtError_lep2", &elePtError_lep2);
    _tree->Branch("eleEta_lep1", &eleEta_lep1);
    _tree->Branch("eleEta_lep2", &eleEta_lep2);
    _tree->Branch("elePhi_lep1", &elePhi_lep1);
    _tree->Branch("elePhi_lep2", &elePhi_lep2);
    _tree->Branch("eleR9_lep1", &eleR9_lep1);
    _tree->Branch("eleR9_lep2", &eleR9_lep2);
    _tree->Branch("eleCalibPt_lep1", &eleCalibPt_lep1);
    _tree->Branch("eleCalibPt_lep2", &eleCalibPt_lep2);
    _tree->Branch("eleDiffCalibOriPt_lep1", &eleDiffCalibOriPt_lep1);
    _tree->Branch("eleDiffCalibOriPt_lep2", &eleDiffCalibOriPt_lep2);
    _tree->Branch("eleCalibEn_lep1", &eleCalibEn_lep1);
    _tree->Branch("eleCalibEn_lep2", &eleCalibEn_lep2);
    _tree->Branch("eleSCEta_lep1", &eleSCEta_lep1);
    _tree->Branch("eleSCEta_lep2", &eleSCEta_lep2);
    _tree->Branch("eleSCPhi_lep1", &eleSCPhi_lep1);
    _tree->Branch("eleSCPhi_lep2", &eleSCPhi_lep2);
    _tree->Branch("eleSCRawEn_lep1", &eleSCRawEn_lep1);
    _tree->Branch("eleSCRawEn_lep2", &eleSCRawEn_lep2);
    _tree->Branch("eleSCEtaWidth_lep1", &eleSCEtaWidth_lep1);
    _tree->Branch("eleSCEtaWidth_lep2", &eleSCEtaWidth_lep2);
    _tree->Branch("eleSCPhiWidth_lep1", &eleSCPhiWidth_lep1);
    _tree->Branch("eleSCPhiWidth_lep2", &eleSCPhiWidth_lep2);
    _tree->Branch("eleHoverE_lep1", &eleHoverE_lep1);
    _tree->Branch("eleHoverE_lep2", &eleHoverE_lep2);
    _tree->Branch("eleEoverP_lep1", &eleEoverP_lep1);
    _tree->Branch("eleEoverP_lep2", &eleEoverP_lep2);
    _tree->Branch("eleEoverPout_lep1", &eleEoverPout_lep1);
    _tree->Branch("eleEoverPout_lep2", &eleEoverPout_lep2);
    _tree->Branch("eleEoverPInv_lep1", &eleEoverPInv_lep1);
    _tree->Branch("eleEoverPInv_lep2", &eleEoverPInv_lep2);
    _tree->Branch("eleBrem_lep1", &eleBrem_lep1);
    _tree->Branch("eleBrem_lep2", &eleBrem_lep2);
    _tree->Branch("eledEtaAtVtx_lep1", &eledEtaAtVtx_lep1);
    _tree->Branch("eledEtaAtVtx_lep2", &eledEtaAtVtx_lep2);
    _tree->Branch("eledPhiAtVtx_lep1", &eledPhiAtVtx_lep1);
    _tree->Branch("eledPhiAtVtx_lep2", &eledPhiAtVtx_lep2);
    _tree->Branch("eleSigmaIEtaIEtaFull5x5_lep1", &eleSigmaIEtaIEtaFull5x5_lep1);
    _tree->Branch("eleSigmaIEtaIEtaFull5x5_lep2", &eleSigmaIEtaIEtaFull5x5_lep2);
    _tree->Branch("eleSigmaIPhiIPhiFull5x5_lep1", &eleSigmaIPhiIPhiFull5x5_lep1);
    _tree->Branch("eleSigmaIPhiIPhiFull5x5_lep2", &eleSigmaIPhiIPhiFull5x5_lep2);
    _tree->Branch("eleESEffSigmaRR_lep1", &eleESEffSigmaRR_lep1);
    _tree->Branch("eleESEffSigmaRR_lep2", &eleESEffSigmaRR_lep2);
    _tree->Branch("elePFChIso_lep1", &elePFChIso_lep1);
    _tree->Branch("elePFChIso_lep2", &elePFChIso_lep2);
    _tree->Branch("elePFPhoIso_lep1", &elePFPhoIso_lep1);
    _tree->Branch("elePFPhoIso_lep2", &elePFPhoIso_lep2);
    _tree->Branch("elePFNeuIso_lep1", &elePFNeuIso_lep1);
    _tree->Branch("elePFNeuIso_lep2", &elePFNeuIso_lep2);
    _tree->Branch("elePFPUIso_lep1", &elePFPUIso_lep1);
    _tree->Branch("elePFPUIso_lep2", &elePFPUIso_lep2);
    _tree->Branch("elePFClusEcalIso_lep1", &elePFClusEcalIso_lep1);
    _tree->Branch("elePFClusEcalIso_lep2", &elePFClusEcalIso_lep2);
    _tree->Branch("elePFClusHcalIso_lep1", &elePFClusHcalIso_lep1);
    _tree->Branch("elePFClusHcalIso_lep2", &elePFClusHcalIso_lep2);
    _tree->Branch("eleIDMVAIso_lep1", &eleIDMVAIso_lep1);
    _tree->Branch("eleIDMVAIso_lep2", &eleIDMVAIso_lep2);
    _tree->Branch("eleIDMVANoIso_lep1", &eleIDMVANoIso_lep1);
    _tree->Branch("eleIDMVANoIso_lep2", &eleIDMVANoIso_lep2);
    _tree->Branch("eleR9Full5x5_lep1", &eleR9Full5x5_lep1);
    _tree->Branch("eleR9Full5x5_lep2", &eleR9Full5x5_lep2);
    _tree->Branch("eleTrkdxy_lep1", &eleTrkdxy_lep1);
    _tree->Branch("eleTrkdxy_lep2", &eleTrkdxy_lep2);
    _tree->Branch("eleKFHits_lep1", &eleKFHits_lep1);
    _tree->Branch("eleKFHits_lep2", &eleKFHits_lep2);
    _tree->Branch("eleKFChi2_lep1", &eleKFChi2_lep1);
    _tree->Branch("eleKFChi2_lep2", &eleKFChi2_lep2);
    _tree->Branch("eleGSFChi2_lep1", &eleGSFChi2_lep1);
    _tree->Branch("eleGSFChi2_lep2", &eleGSFChi2_lep2);
    _tree->Branch("eleFiredSingleTrgs_lep1", &eleFiredSingleTrgs_lep1);
    _tree->Branch("eleFiredSingleTrgs_lep2", &eleFiredSingleTrgs_lep2);
    _tree->Branch("eleFiredDoubleTrgs_lep1", &eleFiredDoubleTrgs_lep1);
    _tree->Branch("eleFiredDoubleTrgs_lep2", &eleFiredDoubleTrgs_lep2);
    _tree->Branch("eleFiredL1Trgs_lep1", &eleFiredL1Trgs_lep1);
    _tree->Branch("eleFiredL1Trgs_lep2", &eleFiredL1Trgs_lep2);
    _tree->Branch("eleIDbit_lep1", &eleIDbit_lep1);
    _tree->Branch("eleIDbit_lep2", &eleIDbit_lep2);
    // Gsf track information
    _tree->Branch("GsfIdxLep1", &GsfIdxLep1);
    _tree->Branch("GsfIdxLep2", &GsfIdxLep2);
    _tree->Branch("gsfMissHits_reco1", &gsfMissHits_reco1);
    _tree->Branch("gsfMissHits_reco2", &gsfMissHits_reco2);
    _tree->Branch("gsfPt_reco1", &gsfPt_reco1);
    _tree->Branch("gsfPt_reco2", &gsfPt_reco2);
    _tree->Branch("gsfEta_reco1", &gsfEta_reco1);
    _tree->Branch("gsfEta_reco2", &gsfEta_reco2);
    _tree->Branch("gsfPhi_reco1", &gsfPhi_reco1);
    _tree->Branch("gsfPhi_reco2", &gsfPhi_reco2);
    _tree->Branch("gsfPt_Lep1", &gsfPt_Lep1);
    _tree->Branch("gsfPt_Lep2", &gsfPt_Lep2);
    _tree->Branch("gsfEta_Lep1", &gsfEta_Lep1);
    _tree->Branch("gsfEta_Lep2", &gsfEta_Lep2);
    _tree->Branch("gsfPhi_Lep1", &gsfPhi_Lep1);
    _tree->Branch("gsfPhi_Lep2", &gsfPhi_Lep2);
    _tree->Branch("nGsfMatchToReco_lep1", &nGsfMatchToReco_lep1);
    _tree->Branch("nGsfMatchToReco_lep2", &nGsfMatchToReco_lep2);
    _tree->Branch("gsfPtSum_lep1", &gsfPtSum_lep1);
    _tree->Branch("gsfPtRatio_lep1", &gsfPtRatio_lep1);
    _tree->Branch("gsfDeltaR_lep1", &gsfDeltaR_lep1);
    _tree->Branch("gsfMissHitsSum_lep1", &gsfMissHitsSum_lep1);
    _tree->Branch("gsfPtSum_lep2", &gsfPtSum_lep2);
    _tree->Branch("gsfPtRatio_lep2", &gsfPtRatio_lep2);
    _tree->Branch("gsfDeltaR_lep2", &gsfDeltaR_lep2);
    _tree->Branch("gsfMissHitsSum_lep2", &gsfMissHitsSum_lep2);
     _tree->Branch("gsfPtSum2_lep1", &gsfPtSum2_lep1);
    _tree->Branch("gsfPtSum2_lep2", &gsfPtSum2_lep2);
    _tree->Branch("gsfRelPtRatio_lep1", &gsfRelPtRatio_lep1);
    _tree->Branch("gsfRelPtRatio_lep2", &gsfRelPtRatio_lep2);
    _tree->Branch("gsfPerpEleSum_lep1", &gsfPerpEleSum_lep1);
    _tree->Branch("gsfPerpEleSum_lep2", &gsfPerpEleSum_lep2);
    _tree->Branch("gsfDiffPtRatio_lep1", &gsfDiffPtRatio_lep1);
    _tree->Branch("gsfDiffPtRatio_lep2", &gsfDiffPtRatio_lep2);

    //pre shower information
    _tree->Branch("eleESEnToRawE_lep1", &eleESEnToRawE_lep1);
    _tree->Branch("eleESEnToRawE_lep2", &eleESEnToRawE_lep2);

    // BC information
    _tree->Branch("circularity_lep1", &circularity_lep1);
    _tree->Branch("circularity_lep2", &circularity_lep2);

    // Conv information
    _tree->Branch("nConv", &nConv_);
    _tree->Branch("convMatched_lep1", &convMatched_lep1);
    _tree->Branch("convMatched_lep2", &convMatched_lep2);
    _tree->Branch("convVtxRadius_lep1", &convVtxRadius_lep1);
    _tree->Branch("convVtxRadius_lep2", &convVtxRadius_lep2);
    _tree->Branch("convP_lep1", &convP_lep1);
    _tree->Branch("convP_lep2", &convP_lep2);
    _tree->Branch("convNTrks_lep1", &convNTrks_lep1);
    _tree->Branch("convNTrks_lep2", &convNTrks_lep2);
    _tree->Branch("convD0_lep1", &convD0_lep1);
    _tree->Branch("convD0_lep2", &convD0_lep2);
    _tree->Branch("convDz_lep1", &convDz_lep1);
    _tree->Branch("convDz_lep2", &convDz_lep2);
    _tree->Branch("convL0_lep1", &convL0_lep1);
    _tree->Branch("convL0_lep2", &convL0_lep2);
    _tree->Branch("convLz_lep1", &convLz_lep1);
    _tree->Branch("convLz_lep2", &convLz_lep2);
   
}

void Minitree::FillMinitree(TreeReader &data, Event &event)
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

    // Conv information
    Int_t nConv             = data.GetInt("nConv");
    float *convVtxRadius    = data.GetPtrFloat("convVtxRadius");
    Int_t* convNTrks        = data.GetPtrInt("convNTrks");
    float* convTrksPin0X    = data.GetPtrFloat("convTrksPin0X");
    float* convTrksPin0Y    = data.GetPtrFloat("convTrksPin0Y");
    float* convTrksPin0Z    = data.GetPtrFloat("convTrksPin0Z");
    float* convFitPairPX    = data.GetPtrFloat("convFitPairPX");
    float* convFitPairPY    = data.GetPtrFloat("convFitPairPY");
    float* convFitPairPZ    = data.GetPtrFloat("convFitPairPZ");
    float* convD0           = data.GetPtrFloat("convD0");
    float* convDz           = data.GetPtrFloat("convDz");
    float* convL0           = data.GetPtrFloat("convL0");
    float* convLz           = data.GetPtrFloat("convLz");

    category = event.category;
    GenReco_case = event.GenReco_case;
    RecoGsf_case = event.RecoGsf_case;
    GsfGen_case = event.GsfGen_case;

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
    mcPt_lep1 = mcPt[event.iGenLep[0]];
    mcEta_lep1 = mcEta[event.iGenLep[0]];
    mcPhi_lep1 = mcPhi[event.iGenLep[0]];
    mcVtx_lep1 = mcVtx[event.iGenLep[0]];
    mcVty_lep1 = mcVty[event.iGenLep[0]];
    mcVtz_lep1 = mcVtz[event.iGenLep[0]];
    mcPt_lep2 = mcPt[event.iGenLep[1]];
    mcEta_lep2 = mcEta[event.iGenLep[1]];
    mcPhi_lep2 = mcPhi[event.iGenLep[1]];
    mcVtx_lep2 = mcVtx[event.iGenLep[1]];
    mcVty_lep2 = mcVty[event.iGenLep[1]];
    mcVtz_lep2 = mcVtz[event.iGenLep[1]];
    mcPt_pho = mcPt[event.iGenPho];
    mcEta_pho = mcEta[event.iGenPho];
    mcPhi_pho = mcPhi[event.iGenPho];
    mcMomPID_lep1 = mcMomPID[event.iGenLep[0]];
    mcGMomPID_lep1 = mcGMomPID[event.iGenLep[0]];
    mcMomPID_lep2 = mcMomPID[event.iGenLep[1]];
    mcGMomPID_lep2 = mcGMomPID[event.iGenLep[1]];
    mcMomPID_pho = mcMomPID[event.iGenPho];
    mcGMomPID_pho = mcGMomPID[event.iGenPho];
    mcStatusFlag_lep1 = mcStatusFlag[event.iGenLep[0]];
    mcStatusFlag_lep2 = mcStatusFlag[event.iGenLep[1]];
    mcStatusFlag_pho = mcStatusFlag[event.iGenPho];

    TLorentzVector GenLep1, GenLep2, GenPho, RecoLep1, RecoLep2, RecoSCLep1, RecoSCLep2, GsfLep1, GsfLep2;
    GenLep1.SetPtEtaPhiM(mcPt[event.iGenLep[0]], mcEta[event.iGenLep[0]], mcPhi[event.iGenLep[0]], 0.000510999);
    GenLep2.SetPtEtaPhiM(mcPt[event.iGenLep[1]], mcEta[event.iGenLep[1]], mcPhi[event.iGenLep[1]], 0.000510999);
    GenPho.SetPtEtaPhiM(mcPt[event.iGenPho], mcEta[event.iGenPho], mcPhi[event.iGenPho], 0.0);
    RecoLep1.SetPtEtaPhiM(elePt[event.iele[0]], eleEta[event.iele[0]], elePhi[event.iele[0]], 0.000510999);
    RecoLep2.SetPtEtaPhiM(elePt[event.iele[1]], eleEta[event.iele[1]], elePhi[event.iele[1]], 0.000510999);
    RecoSCLep1.SetPtEtaPhiM(elePt[event.iele[0]], eleSCEta[event.iele[0]], eleSCPhi[event.iele[0]], 0.000510999);
    RecoSCLep2.SetPtEtaPhiM(elePt[event.iele[1]], eleSCEta[event.iele[1]], eleSCPhi[event.iele[1]], 0.000510999);
    GsfLep1.SetPtEtaPhiM(gsfPt[event.igsf[0]], gsfEta[event.igsf[0]], gsfPhi[event.igsf[0]], 0.000510999);
    GsfLep2.SetPtEtaPhiM(gsfPt[event.igsf[1]], gsfEta[event.igsf[1]], gsfPhi[event.igsf[1]], 0.000510999);

    mcPt_ll = (GenLep1 + GenLep2).Pt();
    mcEta_ll = (GenLep1 + GenLep2).Eta();
    mcPhi_ll = (GenLep1 + GenLep2).Phi();
    mc_Mll = (GenLep1 + GenLep2).M();
    mcPt_llg = (GenLep1 + GenLep2 + GenPho).Pt();
    mcEta_llg = (GenLep1 + GenLep2 + GenPho).Eta();
    mcPhi_llg = (GenLep1 + GenLep2 + GenPho).Phi();
    mc_Mllg = (GenLep1 + GenLep2 + GenPho).M();

    elePresel_lep1 = preselection(data, event.iele[0]);
    elePresel_lep2 = preselection(data, event.iele[1]);
    eleCharge_lep1 = eleCharge[event.iele[0]];
    eleCharge_lep2 = eleCharge[event.iele[1]];
    eleChargeConsistent_lep1 = eleChargeConsistent[event.iele[0]];
    eleChargeConsistent_lep2 = eleChargeConsistent[event.iele[1]];
    eleConvVeto_lep1 = eleConvVeto[event.iele[0]];
    eleConvVeto_lep2 = eleConvVeto[event.iele[1]];
    eleMissHits_lep1 = eleMissHits[event.iele[0]];
    eleMissHits_lep2 = eleMissHits[event.iele[1]];
    eleEcalDrivenSeed_lep1 = eleEcalDrivenSeed[event.iele[0]];
    eleEcalDrivenSeed_lep2 = eleEcalDrivenSeed[event.iele[1]];
    eleEn_lep1 = eleEn[event.iele[0]];
    eleEn_lep2 = eleEn[event.iele[1]];
    eleSCEn_lep1 = eleSCEn[event.iele[0]];
    eleSCEn_lep2 = eleSCEn[event.iele[1]];
    eleEcalEn_lep1 = eleEcalEn[event.iele[0]];
    eleEcalEn_lep2 = eleEcalEn[event.iele[1]];
    eleESEnP1_lep1 = eleESEnP1[event.iele[0]];
    eleESEnP1_lep2 = eleESEnP1[event.iele[1]];
    eleESEnP2_lep1 = eleESEnP2[event.iele[0]];
    eleESEnP2_lep2 = eleESEnP2[event.iele[1]];
    eleD0_lep1 = eleD0[event.iele[0]];
    eleD0_lep2 = eleD0[event.iele[1]];
    eleDz_lep1 = eleDz[event.iele[0]];
    eleDz_lep2 = eleDz[event.iele[1]];
    eleSIP_lep1 = eleSIP[event.iele[0]];
    eleSIP_lep2 = eleSIP[event.iele[1]];
    elePt_lep1 = elePt[event.iele[0]];
    elePt_lep2 = elePt[event.iele[1]];
    elePtError_lep1 = elePtError[event.iele[0]];
    elePtError_lep2 = elePtError[event.iele[1]];
    eleEta_lep1 = eleEta[event.iele[0]];
    eleEta_lep2 = eleEta[event.iele[1]];
    elePhi_lep1 = elePhi[event.iele[0]];
    elePhi_lep2 = elePhi[event.iele[1]];
    eleR9_lep1 = eleR9[event.iele[0]];
    eleR9_lep2 = eleR9[event.iele[1]];
    eleCalibPt_lep1 = eleCalibPt[event.iele[0]];
    eleCalibPt_lep2 = eleCalibPt[event.iele[1]];
    eleDiffCalibOriPt_lep1 = (eleCalibPt[event.iele[0]] - elePt[event.iele[0]]) / elePt[event.iele[0]];
    eleDiffCalibOriPt_lep2 = (eleCalibPt[event.iele[1]] - elePt[event.iele[1]]) / elePt[event.iele[1]];
    eleCalibEn_lep1 = eleCalibEn[event.iele[0]];
    eleCalibEn_lep2 = eleCalibEn[event.iele[1]];
    eleSCEta_lep1 = eleSCEta[event.iele[0]];
    eleSCEta_lep2 = eleSCEta[event.iele[1]];
    eleSCPhi_lep1 = eleSCPhi[event.iele[0]];
    eleSCPhi_lep2 = eleSCPhi[event.iele[1]];
    eleSCRawEn_lep1 = eleSCRawEn[event.iele[0]];
    eleSCRawEn_lep2 = eleSCRawEn[event.iele[1]];
    eleSCEtaWidth_lep1 = eleSCEtaWidth[event.iele[0]];
    eleSCEtaWidth_lep2 = eleSCEtaWidth[event.iele[1]];
    eleSCPhiWidth_lep1 = eleSCPhiWidth[event.iele[0]];
    eleSCPhiWidth_lep2 = eleSCPhiWidth[event.iele[1]];
    eleHoverE_lep1 = eleHoverE[event.iele[0]];
    eleHoverE_lep2 = eleHoverE[event.iele[1]];
    eleEoverP_lep1 = eleEoverP[event.iele[0]];
    eleEoverP_lep2 = eleEoverP[event.iele[1]];
    eleEoverPout_lep1 = eleEoverPout[event.iele[0]];
    eleEoverPout_lep2 = eleEoverPout[event.iele[1]];
    eleEoverPInv_lep1 = eleEoverPInv[event.iele[0]];
    eleEoverPInv_lep2 = eleEoverPInv[event.iele[1]];
    eleBrem_lep1 = eleBrem[event.iele[0]];
    eleBrem_lep2 = eleBrem[event.iele[1]];
    eledEtaAtVtx_lep1 = eledEtaAtVtx[event.iele[0]];
    eledEtaAtVtx_lep2 = eledEtaAtVtx[event.iele[1]];
    eledPhiAtVtx_lep1 = eledPhiAtVtx[event.iele[0]];
    eledPhiAtVtx_lep2 = eledPhiAtVtx[event.iele[1]];
    eleSigmaIEtaIEtaFull5x5_lep1 = eleSigmaIEtaIEtaFull5x5[event.iele[0]];
    eleSigmaIEtaIEtaFull5x5_lep2 = eleSigmaIEtaIEtaFull5x5[event.iele[1]];
    eleSigmaIPhiIPhiFull5x5_lep1 = eleSigmaIPhiIPhiFull5x5[event.iele[0]];
    eleSigmaIPhiIPhiFull5x5_lep2 = eleSigmaIPhiIPhiFull5x5[event.iele[1]];
    eleESEffSigmaRR_lep1 = eleESEffSigmaRR[event.iele[0]];
    eleESEffSigmaRR_lep2 = eleESEffSigmaRR[event.iele[1]];
    elePFChIso_lep1 = elePFChIso[event.iele[0]];
    elePFChIso_lep2 = elePFChIso[event.iele[1]];
    elePFPhoIso_lep1 = elePFPhoIso[event.iele[0]];
    elePFPhoIso_lep2 = elePFPhoIso[event.iele[1]];
    elePFNeuIso_lep1 = elePFNeuIso[event.iele[0]];
    elePFNeuIso_lep2 = elePFNeuIso[event.iele[1]];
    elePFPUIso_lep1 = elePFPUIso[event.iele[0]];
    elePFPUIso_lep2 = elePFPUIso[event.iele[1]];
    elePFClusEcalIso_lep1 = elePFClusEcalIso[event.iele[0]];
    elePFClusEcalIso_lep2 = elePFClusEcalIso[event.iele[1]];
    elePFClusHcalIso_lep1 = elePFClusHcalIso[event.iele[0]];
    elePFClusHcalIso_lep2 = elePFClusHcalIso[event.iele[1]];
    eleIDMVAIso_lep1 = eleIDMVAIso[event.iele[0]];
    eleIDMVAIso_lep2 = eleIDMVAIso[event.iele[1]];
    eleIDMVANoIso_lep1 = eleIDMVANoIso[event.iele[0]];
    eleIDMVANoIso_lep2 = eleIDMVANoIso[event.iele[1]];
    eleR9Full5x5_lep1 = eleR9Full5x5[event.iele[0]];
    eleR9Full5x5_lep2 = eleR9Full5x5[event.iele[1]];
    eleTrkdxy_lep1 = eleTrkdxy[event.iele[0]];
    eleTrkdxy_lep2 = eleTrkdxy[event.iele[1]];
    eleKFHits_lep1 = eleKFHits[event.iele[0]];
    eleKFHits_lep2 = eleKFHits[event.iele[1]];
    eleKFChi2_lep1 = eleKFChi2[event.iele[0]];
    eleKFChi2_lep2 = eleKFChi2[event.iele[1]];
    eleGSFChi2_lep1 = eleGSFChi2[event.iele[0]];
    eleGSFChi2_lep2 = eleGSFChi2[event.iele[1]];
    eleFiredSingleTrgs_lep1 = eleFiredSingleTrgs[event.iele[0]];
    eleFiredSingleTrgs_lep2 = eleFiredSingleTrgs[event.iele[1]];
    eleFiredDoubleTrgs_lep1 = eleFiredDoubleTrgs[event.iele[0]];
    eleFiredDoubleTrgs_lep2 = eleFiredDoubleTrgs[event.iele[1]];
    eleFiredL1Trgs_lep1 = eleFiredL1Trgs[event.iele[0]];
    eleFiredL1Trgs_lep2 = eleFiredL1Trgs[event.iele[1]];
    eleIDbit_lep1 = eleIDbit[event.iele[0]];
    eleIDbit_lep2 = eleIDbit[event.iele[1]];

    eleESEnToRawE_lep1 = (eleESEnP1[event.iele[0]] + eleESEnP2[event.iele[0]])/(eleSCRawEn[event.iele[0]]);
    eleESEnToRawE_lep2 = (eleESEnP1[event.iele[1]] + eleESEnP2[event.iele[1]])/(eleSCRawEn[event.iele[1]]);

    GsfIdxLep1 = event.GsfIdxLep1;
    GsfIdxLep2 = event.GsfIdxLep2;

    if (event.GsfIdxLep1.size() > 0)
    {
        // cloest two gsf trakcs (HR's version)
        Float_t tmp_gsfPtSum = -1., tmp_gsfPtRatio = -1., tmp_dR = 999.;
        Int_t tmp_gsfMissHitsSum = -1;
        for (size_t i = 0; i < event.GsfIdxLep1.size(); i++)
        {
            gsfPt_reco1.push_back(gsfPt[event.GsfIdxLep1[i]]);
            gsfEta_reco1.push_back(gsfEta[event.GsfIdxLep1[i]]);
            gsfPhi_reco1.push_back(gsfPhi[event.GsfIdxLep1[i]]);
            gsfMissHits_reco1.push_back(gsfMissHits[event.GsfIdxLep1[i]]);

            TLorentzVector gsf1;
            gsf1.SetPtEtaPhiM(gsfPt[event.GsfIdxLep1[i]], gsfEta[event.GsfIdxLep1[i]], gsfPhi[event.GsfIdxLep1[i]], 0.000510999);
            for (size_t j = i + 1; j < event.GsfIdxLep1.size(); j++)
            {
                TLorentzVector gsf2;
                gsf2.SetPtEtaPhiM(gsfPt[event.GsfIdxLep1[j]], gsfEta[event.GsfIdxLep1[j]], gsfPhi[event.GsfIdxLep1[j]], 0.000510999);
                if (gsf1.DeltaR(gsf2) < tmp_dR)
                {
                    tmp_gsfPtSum = gsf1.Pt() + gsf2.Pt();
                    tmp_gsfPtRatio = (gsf1.Pt() > gsf2.Pt()) ? (gsf2.Pt()/gsf1.Pt()) : (gsf1.Pt()/gsf2.Pt()); // Make sure that the ratio is always smaller than 1
                    tmp_dR = gsf1.DeltaR(gsf2);
                    tmp_gsfMissHitsSum = gsfMissHits[event.GsfIdxLep1[i]] + gsfMissHits[event.GsfIdxLep1[j]];
                }
            }
        }
        // gsfPtSum_lep1 = tmp_gsfPtSum;
        // gsfPtRatio_lep1 = tmp_gsfPtRatio;
        // gsfDeltaR_lep1 = tmp_dR;
        // gsfMissHitsSum_lep1 = tmp_gsfMissHitsSum;


        // two highest pT gsf trakcs (CH's version)
        if (event.GsfIdxLep1.size() > 1){
            float gsfPt_lead = 0.;
            int gsfPt_lead_ind = 0;
            for (size_t i = 0; i < event.GsfIdxLep1.size(); i++){
                if(gsfPt_lead < gsfPt[event.GsfIdxLep1[i]]) {
                    gsfPt_lead = gsfPt[event.GsfIdxLep1[i]];
                    gsfPt_lead_ind = event.GsfIdxLep1[i];
                }
            }

            float gsfPt_sublead = 0.;
            int gsfPt_sublead_ind = 0;
            for (size_t j = 0; j < event.GsfIdxLep1.size(); j++){
                if (event.GsfIdxLep1[j] == gsfPt_lead_ind) continue;
                if(gsfPt_sublead < gsfPt[event.GsfIdxLep1[j]]) {
                    gsfPt_sublead = gsfPt[event.GsfIdxLep1[j]];
                    gsfPt_sublead_ind = event.GsfIdxLep1[j];
                }
            }
            
            gsfPtSum_lep1 = gsfPt_lead + gsfPt_sublead;
            gsfPtRatio_lep1 = gsfPt_sublead / gsfPt_lead;
            gsfDeltaR_lep1 = deltaR(gsfEta[gsfPt_lead_ind], gsfPhi[gsfPt_lead_ind], gsfEta[gsfPt_sublead_ind], gsfPhi[gsfPt_sublead_ind]);
            gsfMissHitsSum_lep1 = gsfMissHits[gsfPt_lead_ind] + gsfMissHits[gsfPt_sublead_ind];


            TLorentzVector gsf1, gsf2, digsf, ele;
            gsf1.SetPtEtaPhiM(gsfPt[gsfPt_lead_ind], gsfEta[gsfPt_lead_ind], gsfPhi[gsfPt_lead_ind], 0.000511);
            gsf2.SetPtEtaPhiM(gsfPt[gsfPt_sublead_ind], gsfEta[gsfPt_sublead_ind], gsfPhi[gsfPt_sublead_ind], 0.000511);
            ele.SetPtEtaPhiE(eleCalibEn[event.iele[0]]/cosh(eleEta[event.iele[0]]), eleEta[event.iele[0]], elePhi[event.iele[0]], eleCalibEn[event.iele[0]]);
            digsf = gsf1 + gsf2;

            gsfPtSum2_lep1 = (gsfPt_lead * gsfPt_lead) + (gsfPt_sublead * gsfPt_sublead);
            gsfRelPtRatio_lep1 = digsf.Pt() / eleCalibPt[event.iele[0]];
            gsfPerpEleSum_lep1 = (gsf1.Vect().Perp(ele.Vect()) + gsf2.Vect().Perp(ele.Vect()));

            gsfDiffPtRatio_lep1 = (digsf.Pt() - eleCalibPt[event.iele[0]]) / (digsf.Pt() + eleCalibPt[event.iele[0]]);
        }
        else{
            TLorentzVector gsf1, ele;
            gsf1.SetPtEtaPhiM(gsfPt[event.GsfIdxLep1[0]], gsfEta[event.GsfIdxLep1[0]], gsfPhi[event.GsfIdxLep1[0]], 0.000511);
            ele.SetPtEtaPhiE(eleCalibEn[event.iele[0]]/cosh(eleEta[event.iele[0]]), eleEta[event.iele[0]], elePhi[event.iele[0]], eleCalibEn[event.iele[0]]);

            gsfPtSum_lep1 = gsfPt[event.GsfIdxLep1[0]];
            gsfPtRatio_lep1 = -1.;
            gsfDeltaR_lep1 = -999.;
            gsfMissHitsSum_lep1 = gsfMissHits[event.GsfIdxLep1[0]];

            gsfPtSum2_lep1 = (gsfPt[event.GsfIdxLep1[0]] * gsfPt[event.GsfIdxLep1[0]]);
            gsfRelPtRatio_lep1 = gsfPt[event.GsfIdxLep1[0]] / eleCalibPt[event.iele[0]];
            gsfPerpEleSum_lep1 = gsf1.Vect().Perp(ele.Vect());

            gsfDiffPtRatio_lep1 = (gsfPt[event.GsfIdxLep1[0]] - eleCalibPt[event.iele[0]]) / (gsfPt[event.GsfIdxLep1[0]] + eleCalibPt[event.iele[0]]);
        }
        
    }
    else
    {
        gsfPt_reco1.push_back(gsfPt[-1]);
        gsfEta_reco1.push_back(gsfEta[-1]);
        gsfPhi_reco1.push_back(gsfPhi[-1]);
        gsfMissHits_reco1.push_back(gsfMissHits[-1]);
    }

    if (event.GsfIdxLep2.size() > 0)
    {
        // cloest two gsf trakcs (HR's version)
        Float_t tmp_gsfPtSum = -1., tmp_gsfPtRatio = -1., tmp_dR = 999.;
        Int_t tmp_gsfMissHitsSum = -1;
        for (size_t i = 0; i < event.GsfIdxLep2.size(); i++)
        {
            gsfPt_reco2.push_back(gsfPt[event.GsfIdxLep2[i]]);
            gsfEta_reco2.push_back(gsfEta[event.GsfIdxLep2[i]]);
            gsfPhi_reco2.push_back(gsfPhi[event.GsfIdxLep2[i]]);
            gsfMissHits_reco2.push_back(gsfMissHits[event.GsfIdxLep2[i]]);

            TLorentzVector gsf1;
            gsf1.SetPtEtaPhiM(gsfPt[event.GsfIdxLep2[i]], gsfEta[event.GsfIdxLep2[i]], gsfPhi[event.GsfIdxLep2[i]], 0.000510999);
            for (size_t j = i + 1; j < event.GsfIdxLep2.size(); j++)
            {
                TLorentzVector gsf2;
                gsf2.SetPtEtaPhiM(gsfPt[event.GsfIdxLep2[j]], gsfEta[event.GsfIdxLep2[j]], gsfPhi[event.GsfIdxLep2[j]], 0.000510999);
                if (gsf1.DeltaR(gsf2) < tmp_dR)
                {
                    tmp_gsfPtSum = gsf1.Pt() + gsf2.Pt();
                    tmp_gsfPtRatio = (gsf1.Pt() > gsf2.Pt()) ? (gsf2.Pt()/gsf1.Pt()) : (gsf1.Pt()/gsf2.Pt());
                    tmp_dR = gsf1.DeltaR(gsf2);
                    tmp_gsfMissHitsSum = gsfMissHits[event.GsfIdxLep2[i]] + gsfMissHits[event.GsfIdxLep2[j]];
                }
            }
        }
        // gsfPtSum_lep2= tmp_gsfPtSum;
        // gsfPtRatio_lep2 = tmp_gsfPtRatio;
        // gsfDeltaR_lep2 = tmp_dR;
        // gsfMissHitsSum_lep2 = tmp_gsfMissHitsSum;

        // two highest pT gsf trakcs (CH's version)
        if (event.GsfIdxLep2.size() > 1){
            float gsfPt_lead = 0.;
            int gsfPt_lead_ind = 0;
            for (size_t i = 0; i < event.GsfIdxLep2.size(); i++){
                if(gsfPt_lead < gsfPt[event.GsfIdxLep2[i]]) {
                    gsfPt_lead = gsfPt[event.GsfIdxLep2[i]];
                    gsfPt_lead_ind = event.GsfIdxLep2[i];
                }
            }

            float gsfPt_sublead = 0.;
            int gsfPt_sublead_ind = 0;
            for (size_t j = 0; j < event.GsfIdxLep2.size(); j++){
                if (event.GsfIdxLep2[j] == gsfPt_lead_ind) continue;
                if(gsfPt_sublead < gsfPt[event.GsfIdxLep2[j]]) {
                    gsfPt_sublead = gsfPt[event.GsfIdxLep2[j]];
                    gsfPt_sublead_ind = event.GsfIdxLep2[j];
                }
            }
            
            gsfPtSum_lep2 = gsfPt_lead + gsfPt_sublead;
            gsfPtRatio_lep2 = gsfPt_sublead / gsfPt_lead;
            gsfDeltaR_lep2 = deltaR(gsfEta[gsfPt_lead_ind], gsfPhi[gsfPt_lead_ind], gsfEta[gsfPt_sublead_ind], gsfPhi[gsfPt_sublead_ind]);
            gsfMissHitsSum_lep2 = gsfMissHits[gsfPt_lead_ind] + gsfMissHits[gsfPt_sublead_ind];

            TLorentzVector gsf1, gsf2, digsf, ele;
            gsf1.SetPtEtaPhiM(gsfPt[gsfPt_lead_ind], gsfEta[gsfPt_lead_ind], gsfPhi[gsfPt_lead_ind], 0.000511);
            gsf2.SetPtEtaPhiM(gsfPt[gsfPt_sublead_ind], gsfEta[gsfPt_sublead_ind], gsfPhi[gsfPt_sublead_ind], 0.000511);
            ele.SetPtEtaPhiE(eleCalibEn[event.iele[1]]/cosh(eleEta[event.iele[1]]), eleEta[event.iele[1]], elePhi[event.iele[1]], eleCalibEn[event.iele[1]]);
            digsf = gsf1 + gsf2;

            gsfPtSum2_lep2 = (gsfPt_lead * gsfPt_lead) + (gsfPt_sublead * gsfPt_sublead);
            gsfRelPtRatio_lep2 = digsf.Pt() / eleCalibPt[event.iele[1]];
            gsfPerpEleSum_lep2 = (gsf1.Vect().Perp(ele.Vect()) + gsf2.Vect().Perp(ele.Vect()));

            gsfDiffPtRatio_lep2 = (digsf.Pt() - eleCalibPt[event.iele[1]])  / (digsf.Pt() + eleCalibPt[event.iele[1]]);
        }
        else{
            TLorentzVector gsf1, ele;
            gsf1.SetPtEtaPhiM(gsfPt[event.GsfIdxLep2[0]], gsfEta[event.GsfIdxLep2[0]], gsfPhi[event.GsfIdxLep2[0]], 0.000511);
            ele.SetPtEtaPhiE(eleCalibEn[event.iele[1]]/cosh(eleEta[event.iele[1]]), eleEta[event.iele[1]], elePhi[event.iele[1]], eleCalibEn[event.iele[1]]);

            gsfPtSum_lep2 = gsfPt[event.GsfIdxLep2[0]];
            gsfPtRatio_lep2 = -1.;
            gsfDeltaR_lep2 = -999.;
            gsfMissHitsSum_lep2 = gsfMissHits[event.GsfIdxLep2[0]];

            gsfPtSum2_lep2 = (gsfPt[event.GsfIdxLep2[0]] * gsfPt[event.GsfIdxLep2[0]]);
            gsfRelPtRatio_lep2 = gsfPt[event.GsfIdxLep2[0]] / eleCalibPt[event.iele[1]];
            gsfPerpEleSum_lep2 = gsf1.Vect().Perp(ele.Vect());
            gsfDiffPtRatio_lep2 = (gsfPt[event.GsfIdxLep2[0]] - eleCalibPt[event.iele[1]]) / (gsfPt[event.GsfIdxLep2[0]] + eleCalibPt[event.iele[1]]);
        }
    }
    else
    {
        gsfPt_reco2.push_back(gsfPt[-1]);
        gsfEta_reco2.push_back(gsfEta[-1]);
        gsfPhi_reco2.push_back(gsfPhi[-1]);
        gsfMissHits_reco2.push_back(gsfMissHits[-1]);
    }

    gsfPt_Lep1 = gsfPt[event.igsf[0]];
    gsfPt_Lep2 = gsfPt[event.igsf[1]];
    gsfEta_Lep1 = gsfEta[event.igsf[0]];
    gsfEta_Lep2 = gsfEta[event.igsf[1]];
    gsfPhi_Lep1 = gsfPhi[event.igsf[0]];
    gsfPhi_Lep2 = gsfPhi[event.igsf[1]];

    nGsfMatchToReco_lep1 = GsfIdxLep1.size();
    nGsfMatchToReco_lep2 = GsfIdxLep2.size();

    circularity_lep1 = 1 - bcR15[event.BCIdxLep1[0]];
    circularity_lep2 = 1 - bcR15[event.BCIdxLep2[0]];

    // Conv information
    nConv_ = nConv;
    if (nConv > 0){
        convMatched_lep1 = (event.ConvIdxLep1 != -1);
        convMatched_lep2 = (event.ConvIdxLep2 != -1);
        convVtxRadius_lep1 = convVtxRadius[event.ConvIdxLep1];
        convVtxRadius_lep2 = convVtxRadius[event.ConvIdxLep2];
        convNTrks_lep1 = convNTrks[event.ConvIdxLep1];
        convNTrks_lep2 = convNTrks[event.ConvIdxLep2];
        convD0_lep1 = convD0[event.ConvIdxLep1];
        convD0_lep2 = convD0[event.ConvIdxLep2];
        convDz_lep1 = convDz[event.ConvIdxLep1];
        convDz_lep2 = convDz[event.ConvIdxLep2];
        convL0_lep1 = convL0[event.ConvIdxLep1];
        convL0_lep2 = convL0[event.ConvIdxLep2];
        convLz_lep1 = convLz[event.ConvIdxLep1];
        convLz_lep2 = convLz[event.ConvIdxLep2];

        if (convNTrks[event.ConvIdxLep1] == 2)
            convP_lep1 = sqrt(pow(convFitPairPX[event.ConvIdxLep1], 2) + pow(convFitPairPY[event.ConvIdxLep1], 2) + pow(convFitPairPZ[event.ConvIdxLep1], 2));
        else
            convP_lep1 = sqrt(pow(convTrksPin0X[event.ConvIdxLep1], 2) + pow(convTrksPin0Y[event.ConvIdxLep1], 2) + pow(convTrksPin0Z[event.ConvIdxLep1], 2));

        if (convNTrks[event.ConvIdxLep2] == 2)
            convP_lep2 = sqrt(pow(convFitPairPX[event.ConvIdxLep2], 2) + pow(convFitPairPY[event.ConvIdxLep2], 2) + pow(convFitPairPZ[event.ConvIdxLep2], 2));
        else
            convP_lep2 = sqrt(pow(convTrksPin0X[event.ConvIdxLep2], 2) + pow(convTrksPin0Y[event.ConvIdxLep2], 2) + pow(convTrksPin0Z[event.ConvIdxLep2], 2));
    }
    
    

    Fill1DHist(v_hM_Mll, event.category, (GenLep1 + GenLep2).M(), event.mcwei * event.genwei);
    Fill1DHist(v_hM_Mll_30to40, event.category, (GenLep1 + GenLep2).M(), event.mcwei * event.genwei);
    Fill1DHist(v_hM_dR, event.category, GenLep1.DeltaR(GenLep2), event.mcwei * event.genwei);
    Fill1DHist(v_hM_Ptll, event.category, (GenLep1 + GenLep2).Pt(), event.mcwei * event.genwei);
    Fill1DHist(v_hM_LeadLepPt, event.category, GenLep1.Pt(), event.mcwei * event.genwei);
    Fill1DHist(v_hM_TrailLepPt, event.category, GenLep2.Pt(), event.mcwei * event.genwei);
    Fill1DHist(v_hM_Ptllg, event.category, (GenLep1 + GenLep2 + GenPho).Pt(), event.mcwei * event.genwei);
    Fill1DHist(v_hM_sumPtll, event.category, GenLep1.Pt() + GenLep2.Pt(), event.mcwei * event.genwei);
    if (event.category == 1)
    {
        hM_dRSCs->Fill(RecoLep1.DeltaR(RecoLep2), event.mcwei * event.genwei);
        hM_dRSCs_SCEtaPhi->Fill(RecoSCLep1.DeltaR(RecoSCLep2), event.mcwei * event.genwei);
        hM_dRSCs_OrdSC->Fill(RecoLep1.DeltaR(RecoLep2), RecoSCLep1.DeltaR(RecoSCLep2), event.mcwei * event.genwei);
        hM_dEta_OrdSC->Fill(fabs(RecoLep1.Eta() - RecoLep2.Eta()), fabs(RecoSCLep1.Eta() - RecoSCLep2.Eta()), event.mcwei * event.genwei);
        hM_dPhi_OrdSC->Fill(fabs(deltaPhi(RecoLep1.Phi(), RecoLep2.Phi())), fabs(deltaPhi(RecoSCLep1.Phi(), RecoSCLep2.Phi())), event.mcwei * event.genwei);
        hM_eleGsfPtSum_resolved->Fill(gsfPt_Lep1 + gsfPt_Lep2, event.mcwei * event.genwei);
        hM_ele1GsfPt_ele2GsfPt_resolved->Fill(gsfPt[event.igsf[0]], gsfPt[event.igsf[1]], event.mcwei * event.genwei);
    }
    if (event.category == 2)
    {
        hM_ele1GsfPt_ele2GsfPt_merged2Gsf->Fill(gsfPt[event.igsf[0]], gsfPt[event.igsf[1]], event.mcwei * event.genwei);
        hM_eleGsfPtSum_merged2Gsf->Fill(gsfPtSum_lep1, event.mcwei * event.genwei); //gsfPtSum_lep1 is the same as gsfPtSum_lep2 for resolved category
    }

    hM_corr_step12->Fill(event.GenReco_case, event.RecoGsf_case, event.mcwei * event.genwei);
    hM_corr_step23->Fill(event.RecoGsf_case, event.GsfGen_case, event.mcwei * event.genwei);

    // Shower shape
    Fill1DHist_SS(v_hM_eleEn, event.category, eleEn[event.iele[0]], eleEn[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleSCEn, event.category, eleSCEn[event.iele[0]], eleSCEn[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleEcalEn, event.category, eleEcalEn[event.iele[0]], eleEcalEn[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleESEnP1, event.category, eleESEnP1[event.iele[0]], eleESEnP1[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleESEnP2, event.category, eleESEnP2[event.iele[0]], eleESEnP2[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleD0, event.category, eleD0[event.iele[0]], eleD0[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleDz, event.category, eleDz[event.iele[0]], eleDz[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleSIP, event.category, eleSIP[event.iele[0]], eleSIP[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_elePt, event.category, elePt[event.iele[0]], elePt[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_elePtError, event.category, elePtError[event.iele[0]], elePtError[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleEta, event.category, eleEta[event.iele[0]], eleEta[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_elePhi, event.category, elePhi[event.iele[0]], elePhi[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleR9, event.category, eleR9[event.iele[0]], eleR9[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleCalibPt, event.category, eleCalibPt[event.iele[0]], eleCalibPt[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_diffCalibOriPt, event.category, (eleCalibPt[event.iele[0]] - elePt[event.iele[0]]) / elePt[event.iele[0]], (eleCalibPt[event.iele[1]] - elePt[event.iele[1]]) / elePt[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleCalibEn, event.category, eleCalibEn[event.iele[0]], eleCalibEn[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleSCEta, event.category, eleSCEta[event.iele[0]], eleSCEta[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleSCPhi, event.category, eleSCPhi[event.iele[0]], eleSCPhi[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleSCRawEn, event.category, eleSCRawEn[event.iele[0]], eleSCRawEn[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleSCEtaWidth, event.category, eleSCEtaWidth[event.iele[0]], eleSCEtaWidth[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleSCPhiWidth, event.category, eleSCPhiWidth[event.iele[0]], eleSCPhiWidth[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleHoverE, event.category, eleHoverE[event.iele[0]], eleHoverE[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleEoverP, event.category, eleEoverP[event.iele[0]], eleEoverP[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleEoverPout, event.category, eleEoverPout[event.iele[0]], eleEoverPout[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleEoverPInv, event.category, eleEoverPInv[event.iele[0]], eleEoverPInv[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleBrem, event.category, eleBrem[event.iele[0]], eleBrem[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eledEtaAtVtx, event.category, eledEtaAtVtx[event.iele[0]], eledEtaAtVtx[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eledPhiAtVtx, event.category, eledPhiAtVtx[event.iele[0]], eledPhiAtVtx[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleSigmaIEtaIEtaFull5x5, event.category, eleSigmaIEtaIEtaFull5x5[event.iele[0]], eleSigmaIEtaIEtaFull5x5[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleSigmaIPhiIPhiFull5x5, event.category, eleSigmaIPhiIPhiFull5x5[event.iele[0]], eleSigmaIPhiIPhiFull5x5[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS_EBEE(v_hM_eleSigmaIEtaIEtaFull5x5_EB, v_hM_eleSigmaIEtaIEtaFull5x5_EE, event.category, eleSCEta[event.iele[0]], eleSCEta[event.iele[1]], eleSigmaIEtaIEtaFull5x5[event.iele[0]], eleSigmaIEtaIEtaFull5x5[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS_EBEE(v_hM_eleSigmaIPhiIPhiFull5x5_EB, v_hM_eleSigmaIPhiIPhiFull5x5_EE, event.category, eleSCEta[event.iele[0]], eleSCEta[event.iele[1]], eleSigmaIPhiIPhiFull5x5[event.iele[0]], eleSigmaIPhiIPhiFull5x5[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleESEffSigmaRR, event.category, eleESEffSigmaRR[event.iele[0]], eleESEffSigmaRR[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_elePFChIso, event.category, elePFChIso[event.iele[0]], elePFChIso[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_elePFPhoIso, event.category, elePFPhoIso[event.iele[0]], elePFPhoIso[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_elePFNeuIso, event.category, elePFNeuIso[event.iele[0]], elePFNeuIso[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_elePFPUIso, event.category, elePFPUIso[event.iele[0]], elePFPUIso[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_elePFClusEcalIso, event.category, elePFClusEcalIso[event.iele[0]], elePFClusEcalIso[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_elePFClusHcalIso, event.category, elePFClusHcalIso[event.iele[0]], elePFClusHcalIso[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleIDMVAIso, event.category, eleIDMVAIso[event.iele[0]], eleIDMVAIso[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleIDMVANoIso, event.category, eleIDMVANoIso[event.iele[0]], eleIDMVANoIso[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleR9Full5x5, event.category, eleR9Full5x5[event.iele[0]], eleR9Full5x5[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleTrkdxy, event.category, eleTrkdxy[event.iele[0]], eleTrkdxy[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleKFHits, event.category, eleKFHits[event.iele[0]], eleKFHits[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleKFChi2, event.category, eleKFChi2[event.iele[0]], eleKFChi2[event.iele[1]], event.mcwei * event.genwei);
    Fill1DHist_SS(v_hM_eleGSFChi2, event.category, eleGSFChi2[event.iele[0]], eleGSFChi2[event.iele[1]], event.mcwei * event.genwei);

    _tree->Fill();
}