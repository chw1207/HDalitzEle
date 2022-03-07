// #include <fstream.h>
#include <stdlib.h>

#include "/home/chenghan/Analysis/Dalitz/electron/interface/tdrstyle.C"
#include "/home/chenghan/Analysis/Dalitz/electron/interface/CMS_lumi_mod.C"

Float_t Lumi_2017 = 41.5;
Float_t Lumi_2016 = 35.9;
Float_t Lumi_2018 = 59.7;
Float_t Lumi_FullRun2 = 35.9 + 41.5 + 59.7;     // invfb
Float_t refXS = 48.58 + 3.782 + 1.373 + 0.8839; // pb
Float_t BR_eeg = 8.10E-5;
Float_t BR_mmg = 3.90E-5;

std::vector<Float_t> getyield()
{
    Float_t yield_all = 0., yield_PR = 0.;
    std::ifstream file("yield_prod.txt");
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
            yield_all += vval[0];
            yield_PR += vval[1];
        }

        vector<Float_t> vyield;
        vyield.clear();
        vyield.push_back(yield_all);
        vyield.push_back(yield_PR);
        return vyield;
    }
}

void Draw_1DHistStack(vector<TH1F *> vhist, vector<const char *> vleg, const char *drawopt, bool addcum, bool cumhistPR, bool logy, const char *legsty, const char *XaxisName, const char *yaxisunit, const char *outname)
{
    Float_t TickSize = 0.03, AxisTitleSize = 0.05, AxisLabelSize = 0.05;

    setTDRStyle();
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBlueGreenYellow);

    THStack *hs = new THStack("hs", "hs");
    TH1F *hM_cumulative_tmp = new TH1F("hM_cumulative_tmp", "hM_cumulative_tmp", vhist[0]->GetNbinsX(), vhist[0]->GetBinLowEdge(1), vhist[0]->GetBinLowEdge(vhist[0]->GetNbinsX()) + vhist[0]->GetXaxis()->GetBinWidth(1));
    TH1F *hM_cumulative_tmp2 = new TH1F("hM_cumulative_tmp2", "hM_cumulative_tmp2", vhist[0]->GetNbinsX(), vhist[0]->GetBinLowEdge(1), vhist[0]->GetBinLowEdge(vhist[0]->GetNbinsX()) + vhist[0]->GetXaxis()->GetBinWidth(1));
    hM_cumulative_tmp->StatOverflows();
    for (int i = 0; i < vhist.size(); i++)
    {
        hM_cumulative_tmp->Add(vhist[i]);
        if (i > 0)
            hM_cumulative_tmp2->Add(vhist[i]);

        // vhist[i]->Scale(Lumi_FullRun2 / Lumi_2017);
        hs->Add(vhist[i]);
        vhist[i]->GetXaxis()->SetRange(1, vhist[i]->GetNbinsX() + 1);
    }

    // Normalize the histogram for better agreement
    std::vector<Float_t> vecyield = getyield();
    hM_cumulative_tmp->Scale((Lumi_FullRun2 / Lumi_2017) * (vecyield[0] / hM_cumulative_tmp->Integral(-1, -1)));
    hM_cumulative_tmp2->Scale((Lumi_FullRun2 / Lumi_2017) * (vecyield[1] / hM_cumulative_tmp2->Integral(-1, -1)));
    // hM_cumulative_tmp = (TH1F *)hM_cumulative_tmp->GetCumulative();
    TH1F *hM_cumulative = new TH1F("hM_cumulative", "hM_cumulative", vhist[0]->GetNbinsX() + 1, vhist[0]->GetBinLowEdge(1), vhist[0]->GetBinLowEdge(vhist[0]->GetNbinsX()) + 2 * vhist[0]->GetXaxis()->GetBinWidth(1));
    TH1F *hM_cumulative_PR = new TH1F("hM_cumulative_PR", "hM_cumulative_PR", vhist[0]->GetNbinsX() + 1, vhist[0]->GetBinLowEdge(1), vhist[0]->GetBinLowEdge(vhist[0]->GetNbinsX()) + 2 * vhist[0]->GetXaxis()->GetBinWidth(1));
    Float_t cum = 0., cum_PR = 0.;
    for (int i = 1; i <= hM_cumulative_tmp->GetNbinsX() + 1; i++)
    {
        cum += hM_cumulative_tmp->GetBinContent(i);
        hM_cumulative->SetBinContent(i, cum);
        cum_PR += hM_cumulative_tmp2->GetBinContent(i);
        hM_cumulative_PR->SetBinContent(i, cum_PR);
    }

    TCanvas *c = new TCanvas("c", "", 700, 600);
    c->cd();
    if (logy)
        c->SetLogy();
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);
    hs->Draw(drawopt);

    Float_t binwidth = vhist[0]->GetXaxis()->GetBinWidth(1);
    const char *yaxtitletext = "";
    if (strcmp(yaxisunit, "") == 0)
        yaxtitletext = Form("Expected events / %g", binwidth);
    else
        yaxtitletext = Form("Expected events / %g %s", binwidth, yaxisunit);

    hs->GetXaxis()->SetRangeUser(hs->GetHistogram()->GetBinLowEdge(1), hs->GetHistogram()->GetBinLowEdge(hs->GetHistogram()->GetNbinsX()) + 2 * binwidth);
    hs->GetXaxis()->SetTitle(XaxisName);
    hs->GetYaxis()->SetTitle(yaxtitletext);
    if (addcum)
        hs->SetMinimum(hs->GetMinimum() * 0.3);
    else
        hs->SetMinimum(hs->GetMinimum() * 0.5);
    if (addcum)
        hs->SetMaximum(hs->GetMaximum() * 750);
    else
        hs->SetMaximum(hs->GetMaximum() * 2);
    hs->GetXaxis()->SetTickSize(TickSize);
    hs->GetXaxis()->SetTitleSize(AxisTitleSize);
    hs->GetXaxis()->SetLabelSize(AxisLabelSize);
    hs->GetYaxis()->SetTickSize(TickSize);
    hs->GetYaxis()->SetTitleSize(AxisTitleSize);
    hs->GetYaxis()->SetLabelSize(AxisLabelSize);
    hs->GetXaxis()->SetTitleOffset(1.3);
    hs->GetYaxis()->SetTitleOffset(1.4);
    if (addcum)
    {
        gStyle->SetPaintTextFormat(".2f");
        if (cumhistPR)
        {
            hM_cumulative_PR->SetLineWidth(2);
            hM_cumulative_PR->SetLineColor(kBlack);
            hM_cumulative_PR->Draw("HISTtext same");
        }
        else
        {
            hM_cumulative->SetLineWidth(2);
            hM_cumulative->SetLineColor(kBlack);
            hM_cumulative->Draw("HISTtext same");
        }
    }
    // CMS_lumi(c, 4, 11, "41.5 fb^{-1}", 2017, true, "Simulation", "", ""); // CMS Simulation inside the frame
    CMS_lumi(c, 5, 0, Form("%.1f fb^{-1}", Lumi_FullRun2), 2017, true, "Simulation", "", ""); // CMS Simulation outside the frame
    c->Modified();

    Float_t legpos_x0 = 0.57, legpos_x1 = 0.85, legpos_y0 = 0.65, legpos_y1 = 0.89;
    if (addcum)
    {
        legpos_x0 = 0.17;
        legpos_x1 = 0.8;
        legpos_y0 = 0.75;
    }
    TLegend *l = new TLegend(legpos_x0, legpos_y0, legpos_x1, legpos_y1);
    l->SetTextSize(0.035);
    if (addcum)
    {
        l->SetNColumns(2);
        l->SetTextSize(0.03);
    }
    for (int i = 0; i < vhist.size(); i++)
    {
        l->AddEntry(vhist[i], vleg[i], legsty);
    }
    if (addcum)
    {
        if (cumhistPR)
            l->AddEntry(hM_cumulative_PR, "Cumulative (Resolved + Merged-2Gsf + Merged-1MissingGsf)", "l");
        else
            l->AddEntry(hM_cumulative, "Cumulative (All categories)", "l");
    }
    l->SetFillColor(0); //Set the background to be white
    l->SetLineColor(0);
    l->SetFillStyle(0);
    l->Draw("same");

    c->SaveAs(Form("%s.png", outname));
    c->SaveAs(Form("%s.pdf", outname));
}

void Draw_GrAsymErr(vector<TGraphAsymmErrors *> vgr, vector<const char *> vleg, bool logy, vector<Float_t> xlim, vector<Float_t> ylim, const char *legsty, const char *XaxisName, const char *YaxisName, const char *outname)
{
    Float_t TickSize = 0.03, AxisTitleSize = 0.05, AxisLabelSize = 0.05;

    setTDRStyle();
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c", "", 700, 600);
    c->cd();
    if (logy)
        c->SetLogy();
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.15);

    vgr[0]->GetXaxis()->SetLimits(xlim[0], xlim[1]);
    vgr[0]->SetMinimum(ylim[0]);
    vgr[0]->SetMaximum(ylim[1]);
    vgr[0]->GetXaxis()->SetTitle(XaxisName);
    vgr[0]->GetYaxis()->SetTitle(YaxisName);
    vgr[0]->GetXaxis()->SetTickSize(TickSize);
    vgr[0]->GetXaxis()->SetTitleSize(AxisTitleSize);
    vgr[0]->GetXaxis()->SetLabelSize(AxisLabelSize);
    vgr[0]->GetYaxis()->SetTickSize(TickSize);
    vgr[0]->GetYaxis()->SetTitleSize(AxisTitleSize);
    vgr[0]->GetYaxis()->SetLabelSize(AxisLabelSize);
    vgr[0]->GetXaxis()->SetTitleOffset(1.3);
    vgr[0]->GetYaxis()->SetTitleOffset(1.4);
    vgr[0]->Draw("AP");
    for (size_t i = 1; i < vgr.size(); i++)
    {
        vgr[i]->Draw("P same");
    }
    CMS_lumi(c, 5, 0, Form("%.1f fb^{-1}", Lumi_FullRun2), 2017, true, "Simulation", "", "");
    c->RedrawAxis();
    c->Modified();

    TLegend *l = new TLegend(0.35, 0.78, 0.93, 0.9);
    l->SetNColumns(2);
    l->SetTextSize(0.035);
    for (int i = 0; i < vgr.size(); i++)
    {
        l->AddEntry(vgr[i], vleg[i], legsty);
    }
    l->SetFillColor(0); //Set the background to be white
    l->SetLineColor(1);
    l->Draw("same");

    c->SaveAs(Form("%s.png", outname));
    c->SaveAs(Form("%s.pdf", outname));
}

TH1F *gethist(const char *channel, const char *histname)
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
    TH1F *hist = new TH1F();

    for (size_t iprod = 0; iprod < v_prod.size(); iprod++)
    {
        TFile *fin = new TFile(Form("./minitree/2017/Minitree_HDalitz_%s_%s_m125_2017_RECO.root", v_prod[iprod], chan), "READ");
        if (iprod == 0)
        {
            hist = (TH1F *)fin->Get(histname);
            // hist->Scale(v_XSprod[iprod] / refXS); // The original histogram is already mc-weighted
            hist->Scale(Lumi_FullRun2 / Lumi_2017); // The original histogram is mc-weighted to 2017 luminosity -> now project to Full-Run2
            hist->SetDirectory(0);
        }
        else
        {
            TH1F *hist_tmp = (TH1F *)fin->Get(histname);
            // hist_tmp->Scale(v_XSprod[iprod] / refXS); // The original histogram is already mc-weighted
            hist_tmp->Scale(Lumi_FullRun2 / Lumi_2017); // The original histogram is mc-weighted to 2017 luminosity -> now project to Full-Run2
            hist_tmp->SetDirectory(0);
            hist->Add(hist_tmp);
        }

        fin->Close();
    }

    return hist;
}

void plot_GenMatch()
{
    system("mkdir -p ./plots/GenMatchFrac");
    system("mkdir -p ./plots/GenMatchStack");

    vector<const char *> vcolor_frac{"#0b8457", "#cf3030", "#34699a", "#2b2726"};
    vector<const char *> vcolor_stack{"#91BFA5", "#C06463", "#8298CD", "#4D4C4C"};
    vector<const char *> vhname_category{"Total", "NPR", "Merged1MissingGsf", "Merged2Gsf", "Resolved"};
    vector<const char *> vleg{"NPR", "Merged-1MissingGsf", "Merged-2Gsf", "Resolved"};
    vector<const char *> vbasename{"hM_Mll", "hM_Mll_30to40", "hM_dR", "hM_Ptll", "hM_LeadLepPt", "hM_TrailLepPt", "hM_Ptllg"};
    vector<const char *> vXaistitle{"GEN M_{ee} (GeV)", "GEN M_{ee} (GeV)", "GEN #DeltaR(e,e)", "GEN p_{T}^{ee} (GeV)", "GEN p_{T}^{Leading e} (GeV)", "GEN p_{T}^{Trailing e} (GeV)", "GEN p_{T}^{ee#gamma} (GeV)"};
    vector<const char *> vYaxisunit{"GeV", "GeV", "", "GeV", "GeV", "GeV", "GeV"};
    vector<const char *> plotname{"Mll", "Mll_30to40", "DeltaRll", "Ptll", "LeadLepPt", "TrailLepPt", "Ptllg"};
    vector<vector<Float_t>> vxlim{{0, 50}, {30, 40}, {0, 3}, {0, 200}, {0, 200}, {0, 100}, {0, 200}};
    vector<vector<Float_t>> vylim{{0, 1.15}, {0, 1.15}, {0, 1.15}, {0, 1}, {0, 1}, {0, 1.05}, {0, 1}};

    for (size_t ihist = 0; ihist < vbasename.size(); ihist++)
    {
        vector<const char *> vhistname;
        vhistname.clear();
        for (size_t iname = 0; iname < vhname_category.size(); iname++)
        {
            vhistname.push_back(Form("%s_%s", vbasename[ihist], vhname_category[iname]));
        }

        vector<TGraphAsymmErrors *> v_gr;
        vector<TH1F *> v_hM;
        v_gr.clear();
        v_hM.clear();

        TH1F *hm_den = gethist("ele", vhistname[0]);
        Float_t Nevt_total = hm_den->Integral(-1, -1);

        for (size_t i = 1; i < vhistname.size(); i++)
        {
            TH1F *hm_tmp = gethist("ele", vhistname[i]);
            hm_tmp->SetDirectory(0);
            hm_tmp->StatOverflows();
            hm_tmp->GetXaxis()->SetRange(1, hm_tmp->GetNbinsX() + 1); // Overflow bin
            hm_tmp->SetLineWidth(0);
            hm_tmp->SetLineColor(kBlack);
            hm_tmp->SetFillColor(TColor::GetColor(vcolor_stack[i - 1]));
            // hm_tmp->Scale((refXS * 1000 * BR_eeg * Lumi_FullRun2) / Nevt_total);
            v_hM.push_back(hm_tmp);

            // if (i == 1)
            // hm_den->Scale((refXS * 1000 * BR_eeg * Lumi_FullRun2) / Nevt_total);
            TGraphAsymmErrors *frac = new TGraphAsymmErrors(hm_tmp, hm_den);
            frac->SetMarkerStyle(kFullCircle);
            frac->SetLineWidth(2);
            frac->SetLineColor(TColor::GetColor(vcolor_frac[i - 1]));
            frac->SetMarkerColor(TColor::GetColor(vcolor_frac[i - 1]));
            frac->SetMarkerSize(1.3);
            v_gr.push_back(frac);
        }
        Draw_1DHistStack(v_hM, vleg, "hist", false, false, true, "f", vXaistitle[ihist], vYaxisunit[ihist], Form("./plots/GenMatchStack/%s_GenMatch_Stack", plotname[ihist]));
        Draw_1DHistStack(v_hM, vleg, "hist", true, false, true, "f", vXaistitle[ihist], vYaxisunit[ihist], Form("./plots/GenMatchStack/%s_GenMatch_Stack_cumulative_all", plotname[ihist]));
        Draw_1DHistStack(v_hM, vleg, "hist", true, true, true, "f", vXaistitle[ihist], vYaxisunit[ihist], Form("./plots/GenMatchStack/%s_GenMatch_Stack_cumulative_PR", plotname[ihist]));
        Draw_GrAsymErr(v_gr, vleg, false, vxlim[ihist], vylim[ihist], "PELX0", vXaistitle[ihist], "Fraction", Form("./plots/GenMatchFrac/%s_GenMatch_frac", plotname[ihist]));
    }
}