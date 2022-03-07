#include "/afs/cern.ch/work/h/hajheng/private/HDalitz/interface/tdrstyle.C"
#include "/afs/cern.ch/work/h/hajheng/private/HDalitz/interface/CMS_lumi_mod.C"

using namespace std;

Float_t Lumi_2017 = 41.5;
Float_t refXS = 48.58 + 3.782 + 1.373 + 0.8839;
Float_t BR_eeg = 8.07E-5;
Float_t BR_mmg = 3.83E-5;

void Pal_blue()
{
    static Int_t colors[1000];
    static Bool_t initialized = kFALSE;
    Double_t Red[5] = {250. / 255., 178. / 255., 141. / 255., 104. / 255., 67. / 255.};
    Double_t Green[5] = {250. / 255., 201. / 255., 174. / 255., 146. / 255., 120. / 255.};
    Double_t Blue[5] = {250. / 255., 230. / 255., 216. / 255., 202. / 255., 188. / 255.};
    Double_t Length[5] = {0.00, 0.25, 0.50, 0.75, 1.00};
    if (!initialized)
    {
        Int_t FI = TColor::CreateGradientColorTable(5, Length, Red, Green, Blue, 1000);
        for (int i = 0; i < 1000; i++)
            colors[i] = FI + i;
        initialized = kTRUE;
        return;
    }
    gStyle->SetPalette(1000, colors);
}

void Draw_single2Dhist(TH2F *hist, Int_t step, bool logz, const char *XaxisName, const char *YaxisName, const char *outname)
{
    gStyle->SetOptStat(0);
    // gStyle->SetPalette(kBlueGreenYellow);

    hist->Scale(100.); // Percentage
    vector<const char *> vbinlable_x_step12{"2Gen2Reco", "2Gen1Reco", "2Gen1Match1N", "2Gen2N"};
    vector<const char *> vbinlable_y_step12{"2Reco2Gsf", "1Reco2Gsf", "1Reco1Gsf", "0Reco0Gsf"};
    vector<const char *> vbinlable_x_step23{"2Reco2Gsf", "1Reco2Gsf", "1Reco1Gsf", "0Reco0Gsf"};
    vector<const char *> vbinlable_y_step23{"2Gsf2Gen", "1Gsf2Gen", "0Gsf2Gen"};

    Float_t TickSize = 0.03, AxisTitleSize = 0.05, AxisLabelSize = 0.04;
    Float_t LeftMargin = 0.27;

    TCanvas *c = new TCanvas("c", "", 700, 600);
    c->cd();
    if (logz)
        c->SetLogz();
    gPad->SetRightMargin(0.13);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(LeftMargin);
    gPad->SetBottomMargin(0.18);

    if (step == 12)
    {
        for (size_t i = 0; i < vbinlable_x_step12.size(); i++)
            hist->GetXaxis()->SetBinLabel(i + 1, vbinlable_x_step12[i]);
        for (size_t i = 0; i < vbinlable_y_step12.size(); i++)
            hist->GetYaxis()->SetBinLabel(i + 1, vbinlable_y_step12[i]);
    }
    else if (step == 23)
    {
        for (size_t i = 0; i < vbinlable_x_step23.size(); i++)
            hist->GetXaxis()->SetBinLabel(i + 1, vbinlable_x_step23[i]);
        for (size_t i = 0; i < vbinlable_y_step23.size(); i++)
            hist->GetYaxis()->SetBinLabel(i + 1, vbinlable_y_step23[i]);
    }
    else
    {
        cout << "[WARNING] step " << step << " is neither 12 nor 23. It's not correct! Please change!" << endl;
        for (size_t i = 0; i < vbinlable_x_step12.size(); i++)
            hist->GetXaxis()->SetBinLabel(i + 1, vbinlable_x_step12[i]);
        for (size_t i = 0; i < vbinlable_y_step12.size(); i++)
            hist->GetYaxis()->SetBinLabel(i + 1, vbinlable_y_step12[i]);
    }

    // hist->SetMinimum(1.);
    hist->GetXaxis()->SetTitle(XaxisName);
    hist->GetYaxis()->SetTitle(YaxisName);
    hist->GetXaxis()->SetTickSize(TickSize);
    hist->GetYaxis()->SetTickSize(TickSize);
    hist->GetXaxis()->SetTitleSize(AxisTitleSize);
    hist->GetYaxis()->SetTitleSize(AxisTitleSize);
    hist->GetXaxis()->SetLabelSize(AxisLabelSize * 1.3);
    hist->GetYaxis()->SetLabelSize(AxisLabelSize * 1.3);
    hist->GetXaxis()->SetTitleOffset(1.7);
    hist->GetYaxis()->SetTitleOffset(2.5);
    hist->GetZaxis()->SetLabelSize(AxisLabelSize);
    hist->SetContour(1000);
    hist->SetLineColor(1);
    hist->SetMarkerSize(1.8);
    gStyle->SetPaintTextFormat(".3g %%");
    hist->Draw("COLZ TEXT");
    static TExec *ex1 = new TExec("ex1", "Pal_blue();");
    ex1->Draw();
    hist->Draw("COLZ TEXT same");
    CMS_lumi(c, 5, 0, "", 2017, true, "Simulation", "", "");
    c->RedrawAxis();

    c->SaveAs(Form("%s.png", outname));
    c->SaveAs(Form("%s.pdf", outname));
}

std::vector<Float_t> getmcwei()
{
    vector<Float_t> vmcwei;
    vmcwei.clear();
    std::ifstream file("mcwei.txt");
    if (!file)
    {
        cerr << "Can't open file!\n";
        exit(1);
    }
    else
    {
        std::string line;
        while (std::getline(file, line))
        {
            std::stringstream sst(line);
            std::string a;
            std::vector<Float_t> vval;
            vval.clear();
            while (getline(sst, a, ' '))
            {
                vval.push_back(stof(a));
            }
            vmcwei.push_back(vval[0]);
        }
        return vmcwei;
    }
}

TH2F *gethist(const char *channel = "ele", bool normalize = true, const char *histname = "hM_corr")
{
    const char *chan = "";
    if (strcmp(channel, "ele") == 0)
    {
        chan = "eeg";
    }
    else if (strcmp(channel, "mu") == 0)
    {
        chan = "mmg";
    }
    else
    {
        cout << "[WARNING] " << channel << " is neither [ele] nor [mu]. Use default [eeg]. " << endl;
        chan = "eeg";
    }

    vector<const char *> v_prod{"ggF", "VBF", "WH", "ZH"};
    vector<Float_t> v_XSprod{48.58, 3.782, 1.373, 0.8839}; // unit: pb
    vector<Float_t> v_mcwei = getmcwei();
    TH2F *hist = new TH2F();

    for (size_t iprod = 0; iprod < v_prod.size(); iprod++)
    {
        TFile *fin = new TFile(Form("./minitree/Minitree_HDalitz_%s_%s_m125_2017_RECO.root", v_prod[iprod], chan), "READ");
        if (iprod == 0)
        {
            hist = (TH2F *)fin->Get(histname);
            // hist->Scale(v_mcwei[iprod]);
            cout << hist->Integral(-1,-1) << endl;
            hist->SetDirectory(0);
        }
        else
        {
            TH2F *hist_tmp = (TH2F *)fin->Get(histname);
            // hist_tmp->Scale(v_mcwei[iprod]);
            cout << hist_tmp->Integral(-1,-1) << endl;
            hist_tmp->SetDirectory(0);
            hist->Add(hist_tmp);
        }

        fin->Close();
    }

    if (normalize)
        hist->Scale(1. / hist->Integral(-1, -1));

    return hist;
}

void plot_corr()
{
    setTDRStyle();

    system("mkdir -p ./plots/Corr");

    vector<const char *> v_prod{"ggF", "VBF", "WH", "ZH"};
    for (size_t iprod = 0; iprod < v_prod.size(); iprod++)
    {
        TFile *f = new TFile(Form("./minitree/Minitree_HDalitz_%s_eeg_m125_2017_RECO.root", v_prod[iprod]), "READ");
        TH2F *hM_corr_step12 = (TH2F *)f->Get("hM_corr_step12");
        hM_corr_step12->SetDirectory(0);
        cout << hM_corr_step12->Integral(-1,-1) << endl;
        hM_corr_step12->Scale(1. / hM_corr_step12->Integral(-1, -1));
        // Draw_single2Dhist(TH2F *hist, bool logz, const char *XaxisName, const char *YaxisName, const char *outname)
        Draw_single2Dhist(hM_corr_step12, 12, false, "Gen to Reco", "Reco to Gsf", Form("./plots/Corr/MatchingCaseCorr_step12_%s", v_prod[iprod]));

        TH2F *hM_corr_step23 = (TH2F *)f->Get("hM_corr_step23");
        hM_corr_step23->SetDirectory(0);
        cout << hM_corr_step23->Integral(-1,-1) << endl;
        hM_corr_step23->Scale(1. / hM_corr_step23->Integral(-1, -1));
        Draw_single2Dhist(hM_corr_step23, 23, false, "Reco to Gsf", "Gsf to Gen", Form("./plots/Corr/MatchingCaseCorr_step23_%s", v_prod[iprod]));
        f->Close();
    }

    TH2F *hM_corr_step12_allprod = gethist("ele",true,"hM_corr_step12");
    Draw_single2Dhist(hM_corr_step12_allprod, 12, false, "Gen to Reco", "Reco to Gsf", "./plots/Corr/MatchingCaseCorr_step12_allprod");
    TH2F *hM_corr_step23_allprod = gethist("ele",true,"hM_corr_step23");
    Draw_single2Dhist(hM_corr_step23_allprod, 23, false, "Reco to Gsf", "Gsf to Gen", "./plots/Corr/MatchingCaseCorr_step23_allprod");
}