# include "/home/chenghan/Analysis/Dalitz/electron/interface/CMS_lumi_mod.C"
using namespace std;
string LineName[2] = {"#202020","#ec0101"};

void Draw2Dhist(TH2F* h, string outName){
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);

    h->Scale(1. / h->Integral(-1, -1, -1, -1));
    
    float TickSize = 0.02, TitleSize = 0.04, LabelSize = TitleSize*1.1;
    TCanvas *c = new TCanvas("c", "c", 700, 700);
    c->cd();
    gPad->SetRightMargin(0.14);
    gPad->SetBottomMargin(0.12);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.12);

    h->GetXaxis()->SetTitle("P^{1st gsf}_{T} [GeV]");
    h->GetXaxis()->SetMoreLogLabels();
    h->GetXaxis()->SetTickSize(TickSize);
    h->GetXaxis()->SetTitleSize(TitleSize);
    h->GetXaxis()->SetLabelSize(LabelSize);
    h->GetXaxis()->SetTitleOffset(1.4);

    h->GetYaxis()->SetTitle("P^{ 2nd gsf}_{T} [GeV]");
    h->GetYaxis()->SetTickSize(TickSize);
    h->GetYaxis()->SetTitleSize(TitleSize);
    h->GetYaxis()->SetLabelSize(LabelSize);
    h->GetYaxis()->SetTitleOffset(1.3);

    h->GetZaxis()->SetLabelSize(0.03);
    h->GetZaxis()->SetRangeUser(0, 0.0015);
    h->GetZaxis()->SetNdivisions(510);
    h->Draw("COLZ");

    CMS_lumi(c, 4, 11, "41.5 fb^{-1}", 2017, true, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "");
    c->RedrawAxis();


    c->SaveAs(outName.c_str());
    c->Close();
}

void Draw1Dhist(vector<TH1F*> vh, string outName){
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);

    TH1F* h_sig = (TH1F *)vh[0]->Clone();
    for (int i = 1; i < 4; i ++){
        h_sig->Add(vh[i],1);
    }
    TH1F* h_bkg = (TH1F *)vh[4]->Clone();
    for (int i = 5; i < 8; i ++){
        h_bkg->Add(vh[i],1);
    }

    h_sig->Scale(1. / h_sig->Integral(-1, -1));
    h_bkg->Scale(1. / h_bkg->Integral(-1, -1));

    TCanvas* c = new TCanvas("c","c",800,800);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.07);
    c->SetLeftMargin(0.14);
    c->SetBottomMargin(0.15);
    // c->SetLogy();
    c->cd();
    
    float ymax = h_sig->GetBinContent(h_sig->GetMaximumBin()) * 1.2;
    h_sig->GetXaxis()->SetTitle("P^{1st gsf}_{T} + P^{2nd gsf}_{T} [GeV]");
    h_sig->GetXaxis()->SetTickSize(0.02);
    h_sig->GetXaxis()->SetTitleSize(0.04);
    h_sig->GetXaxis()->SetLabelSize(0.04);
    h_sig->GetXaxis()->SetLabelOffset(0.02);
    h_sig->GetXaxis()->SetTitleOffset(1.4);
    h_sig->GetYaxis()->SetTitle("a.u.");
    // h_sig->GetXaxis()->SetRangeUser(110, 140);
    h_sig->GetYaxis()->SetNdivisions(506);
    // h_sig->SetMaximum(ymax);
    h_sig->GetYaxis()->SetRangeUser(0, ymax);
    h_sig->GetYaxis()->SetTickSize(0.02);
    h_sig->GetYaxis()->SetTitleSize(0.04);
    h_sig->GetYaxis()->SetTitleOffset(1.6);
    h_sig->GetYaxis()->SetLabelSize(0.04);

    h_sig->SetLineColor(TColor::GetColor(LineName[0].c_str()));
    h_sig->SetLineWidth(2);
    h_sig->Draw("hist");

    h_bkg->SetLineColor(TColor::GetColor(LineName[1].c_str()));
    h_bkg->SetLineWidth(2);
    h_bkg->Draw("hist same");

    TLegend *legend = new TLegend(0.6, 0.8, 0.85, 0.9);//(x1,y1,x2,y2)
    legend->SetTextSize(0.04);
    legend->AddEntry(h_sig, "Signal", "l");
    legend->AddEntry(h_bkg, "Background", "l");
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->Draw();

    CMS_lumi(c, 4, 11, "41.5 fb^{-1}", 2017, true, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "");

    c->SaveAs(outName.c_str());
    c->Close();
}


void gsfPtSum(){

    vector<string> vfiles = {
        "./minitree/2017/Minitree_HDalitz_ggF_eeg_m125_2017_RECO.root",
        "./minitree/2017/Minitree_HDalitz_VBF_eeg_m125_2017_RECO.root",
        "./minitree/2017/Minitree_HDalitz_WH_eeg_m125_2017_RECO.root",
        "./minitree/2017/Minitree_HDalitz_ZH_eeg_m125_2017_RECO.root"
    }; 

    TChain* tree1 = new TChain("outTree");
    for(auto i : vfiles){
        tree1->Add(i.c_str());
    }

    TH1F* hsig = new TH1F("hsig", "", 100, 0, 200);
    TH2F* hgsf = new TH2F("hgsfPt", "", 100, 0, 100, 90, 0, 90);
    int category;
    float gsfPtSum_lep1, eleCalibPt_lep1;
    int eleEcalDrivenSeed_lep1;
    bool elePresel_lep1;
    float mcwei, genwei;
    float gsfPt_Lep1, gsfPt_Lep2;
    tree1->SetBranchStatus("*", 0);
    tree1->SetBranchStatus("category", 1);
    tree1->SetBranchStatus("gsfPtSum_lep1", 1);
    tree1->SetBranchStatus("elePresel_lep1", 1);
    tree1->SetBranchStatus("eleEcalDrivenSeed_lep1", 1);
    tree1->SetBranchStatus("eleCalibPt_lep1", 1);
    tree1->SetBranchStatus("mcwei", 1);
    tree1->SetBranchStatus("genwei", 1);
    tree1->SetBranchStatus("gsfPt_Lep1", 1);
    tree1->SetBranchStatus("gsfPt_Lep2", 1);

    tree1->SetBranchAddress("category", &category);
    tree1->SetBranchAddress("gsfPtSum_lep1", &gsfPtSum_lep1);
    tree1->SetBranchAddress("elePresel_lep1", &elePresel_lep1);
    tree1->SetBranchAddress("eleEcalDrivenSeed_lep1", &eleEcalDrivenSeed_lep1);
    tree1->SetBranchAddress("eleCalibPt_lep1", &eleCalibPt_lep1);
    tree1->SetBranchAddress("mcwei", &mcwei);
    tree1->SetBranchAddress("genwei", &genwei);
    tree1->SetBranchAddress("gsfPt_Lep1", &gsfPt_Lep1);
    tree1->SetBranchAddress("gsfPt_Lep2", &gsfPt_Lep2);

    int tev = tree1->GetEntries();
    for (int ev = 0; ev < tev; ev++){
        tree1->GetEntry(ev);

        if (category != 2) continue;
        if (elePresel_lep1 != 1) continue;
        if (eleCalibPt_lep1 < 25) continue;
        if (gsfPt_Lep1 < gsfPt_Lep2) continue;

        float wei = mcwei * genwei;
        hsig->Fill(gsfPtSum_lep1, wei);
        hgsf->Fill(gsfPt_Lep1, gsfPt_Lep2, wei);
    }

    Draw2Dhist(hgsf, "./gsfSumPt_2D_sinal.pdf");
    // Draw2Dhist(vgsfPt, "../plot/gsfSumPt_2D_bkg.pdf", true);

    // Draw1Dhist(vgsfPtSum, "./gsfSumPt_1D.pdf");
}