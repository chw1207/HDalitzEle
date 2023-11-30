R__LOAD_LIBRARY($HDalitzEle_LOC/lib_ana/libHDalitzEle.so)
R__ADD_INCLUDE_PATH($HDalitzEle_LOC/include)
#include <iostream>
#include <string>
#include <vector>
#include "boost/filesystem.hpp"
#include "yaml-cpp/yaml.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "tdrstyle.h"
#include "CMS_lumi_modification.h"


using namespace std;
namespace bfs = boost::filesystem;

void draw_var(int iBE, ROOT::RDF::RNode df_data, ROOT::RDF::RNode df_GJets, ROOT::RDF::RNode df_QCD, ROOT::RDF::RNode df_TT, ROOT::RDF::RNode df_Zg, ROOT::RDF::RNode df_sig, string varName, string sel, vector<double> bins, string legName){
    // preprocess mc dataframe
    auto nf_GJets = df_GJets;
    if (!nf_GJets.HasColumn(varName))
        nf_GJets = nf_GJets.Define(varName, sel);
    
    auto nf_QCD = df_QCD;
    if (!nf_QCD.HasColumn(varName))
        nf_QCD = nf_QCD.Define(varName, sel);

    auto nf_TT = df_TT;
    if (!nf_TT.HasColumn(varName))
        nf_TT = nf_TT.Define(varName, sel);

    auto nf_Zg = df_Zg;
    if (!nf_Zg.HasColumn(varName))
        nf_Zg = nf_Zg.Define(varName, sel);

    auto nf_sig = df_sig;
    if (!nf_sig.HasColumn(varName))
        nf_sig = nf_sig.Define(varName, sel);

    // preprocess data dataframe
    auto nf_data = df_data;
    if (!nf_data.HasColumn(varName))
        nf_data = nf_data.Define(varName, sel).Filter("(CMS_higgs_mass > 110 && CMS_higgs_mass <= 120) || (CMS_higgs_mass >= 130 && CMS_higgs_mass < 170)");

    auto h_GJets = nf_GJets.Histo1D({Form("h_%s_GJets", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName, "weight");
    auto h_QCD = nf_QCD.Histo1D({Form("h_%s_QCD", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName, "weight");
    auto h_TT = nf_TT.Histo1D({Form("h_%s_TT", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName, "weight");
    auto h_Zg = nf_Zg.Histo1D({Form("h_%s_Zg", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName, "weight");
    auto h_sig = nf_sig.Histo1D({Form("h_%s_sig", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName, "weight");
    h_sig->Scale(50);
    auto h_data = nf_data.Histo1D({Form("h_%s_data", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName);
    // cout << h_data->Integral()/h_mc->Integral() << endl;
    // h_mc->Scale(h_data->Integral()/h_mc->Integral());
    

    setTDRStyle();
    auto canv = new TCanvas("c", "c", 800, 800);
    canv->cd();

    TPad* pad1 = new TPad("pad1", " ", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.05);
    pad1->SetTopMargin(0.08);
    pad1->SetRightMargin(0.05);
    pad1->SetLeftMargin(0.13);
    pad1->SetBottomMargin(0.03);
    pad1->Draw();    
    pad1->SetLogy();         
    pad1->cd();

    pad2->SetRightMargin(0.05);
    pad2->SetLeftMargin(0.13);
    pad2->SetTopMargin(0.06);
    pad2->SetBottomMargin(0.4);

    auto hs = new THStack("hs", " ");
    h_Zg->SetFillColor(TColor::GetColor(100, 192, 232));
    h_Zg->SetLineColor(TColor::GetColor(100, 192, 232));
    hs->Add(h_Zg.GetPtr());

    h_TT->SetFillColor(TColor::GetColor(222, 90, 106));
    h_TT->SetLineColor(TColor::GetColor(222, 90, 106));
    hs->Add(h_TT.GetPtr());

    h_QCD->SetFillColor(TColor::GetColor("#00A88F"));
    h_QCD->SetLineColor(TColor::GetColor("#00A88F"));
    hs->Add(h_QCD.GetPtr());

    h_GJets->SetFillColor(TColor::GetColor(248, 206, 104));
    h_GJets->SetLineColor(TColor::GetColor(248, 206, 104));
    hs->Add(h_GJets.GetPtr());

    h_data->GetXaxis()->SetLabelOffset(0.05);
    h_data->GetXaxis()->SetTickSize(0.03);
    h_data->GetYaxis()->SetTitle(Form("Events"));
    h_data->GetYaxis()->SetTitleSize(0.07);
    h_data->GetYaxis()->SetRangeUser(1, h_data->GetBinContent(h_data->GetMaximumBin()) * 100);
    h_data->GetYaxis()->SetTickSize(0.03);
    h_data->GetYaxis()->SetTitleSize(0.07);
    h_data->GetYaxis()->SetLabelSize(0.055);
    h_data->GetYaxis()->SetTitleOffset(0.8);

    h_data->SetMarkerColor(TColor::GetColor("#202020"));
    h_data->SetMarkerSize(1.1);
    h_data->SetMarkerStyle(20);
    h_data->SetLineColor(TColor::GetColor("#202020"));
    h_data->SetLineWidth(2);
    h_data->Draw("EP");
    // cout << h_data->Integral() << endl;;
    // cout << h_sig->Integral() << endl;;

    hs->Draw("hist same");

    h_sig->SetLineColor(TColor::GetColor("#1F3C88"));
    h_sig->SetLineWidth(3);
    h_sig->Draw("hist same");

    h_data->Draw("EP same");

    TString lumi_text("59.82 fb^{-1}");
    int year = 2018;
    TString extra_text("Work-in-progress");
    TString proc_text("H #rightarrow #gamma* #gamma #rightarrow ee#gamma");
    CMS_lumi(pad1 , 5, 10, lumi_text, year, true, extra_text, proc_text, "");
    pad1->RedrawAxis();
    pad1->Update();

    auto leg = new TLegend(0.69, 0.62, 0.9, 0.88);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->AddEntry(h_data.GetPtr(),  "Data", "LE1P");
    leg->AddEntry(h_sig.GetPtr(), "Signal #times 50", "f");
    leg->AddEntry(h_GJets.GetPtr(), "GJets", "f");
    leg->AddEntry(h_QCD.GetPtr(), "QCD", "f");
    leg->AddEntry(h_TT.GetPtr(), "TT", "f");
    leg->AddEntry(h_Zg.GetPtr(), "Zg", "f");
    leg->Draw("same");

    auto ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.05);
    ltx->DrawLatex(0.17, 0.7, (iBE == 0) ? "ECAL Barrel" : "ECAL Endcap");

    canv->cd();

    //---------- 2nd Pad's Setting(Draw the Scale Factors) ----------//
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    pad2->SetGridy();
    pad2->SetRightMargin(0.05);
    pad2->SetLeftMargin(0.13);
    pad2->SetTopMargin(0.06);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();
    
    auto h_ratio = (TH1F*) h_data.GetPtr()->Clone();
    auto hs_hist = (TH1F*) h_GJets.GetPtr()->Clone();
    hs_hist->Add(h_QCD.GetPtr(), 1);
    hs_hist->Add(h_TT.GetPtr(), 1);
    hs_hist->Add(h_Zg.GetPtr(), 1);
    h_ratio->Divide(hs_hist);

    h_ratio->GetXaxis()->SetTitle(legName.c_str());
    h_ratio->GetYaxis()->SetTitle("Data / MC");
    h_ratio->GetYaxis()->SetRangeUser(0, 6);

    h_ratio->SetMarkerColor(TColor::GetColor("#202020"));
    h_ratio->SetMarkerSize(1.1);
    h_ratio->SetMarkerStyle(20);
    h_ratio->SetLineColor(TColor::GetColor("#202020"));
    h_ratio->SetLineWidth(2);
    h_ratio->GetXaxis()->SetTitleSize(0.16);
    h_ratio->GetXaxis()->SetTitleOffset(1);
    h_ratio->GetXaxis()->SetLabelSize(0.12);
    h_ratio->GetXaxis()->SetLabelOffset(0.03);
    h_ratio->GetYaxis()->SetTitleSize(0.13);
    h_ratio->GetYaxis()->SetTitleOffset(0.4);
    h_ratio->GetYaxis()->SetLabelSize(0.12);
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->Draw("EP");

    string plot_path = (iBE == 0) ? "../plots/UL2018_eeg_bkg_EB" : "../plots/UL2018_eeg_bkg_EE";
    if (!bfs::exists(plot_path)){
        system(Form("mkdir -p %s", plot_path.c_str()));
    }
    canv->Print(Form("%s/%s.pdf", plot_path.c_str(), varName.c_str()));
    delete canv;
    delete leg;
}

void drawBkgComp(int iBE){
    // ROOT::EnableImplicitMT();

    string eta_cut = (iBE == 0) ? "abs(eleSCEta_Lead) < 1.4442" : "abs(eleSCEta_Lead) > 1.4442";

    vector<string> data_files = {
        "/data4/chenghan/electron/miniTree_HGGCheck_WPTight_merged/UL2018/miniTree_EGamma_Run2018A_UL2018.root",
        "/data4/chenghan/electron/miniTree_HGGCheck_WPTight_merged/UL2018/miniTree_EGamma_Run2018B_UL2018.root",
        "/data4/chenghan/electron/miniTree_HGGCheck_WPTight_merged/UL2018/miniTree_EGamma_Run2018C_UL2018.root",
        "/data4/chenghan/electron/miniTree_HGGCheck_WPTight_merged/UL2018/miniTree_EGamma_Run2018D_UL2018.root"
    };
    ROOT::RDF::RNode df_data    = ROOT::RDataFrame("miniTree", data_files).Filter("abs(eleSCEta_Lead) < 1.4442");
    ROOT::RDF::RNode df_GJets   = ROOT::RDataFrame("miniTree", "/data4/chenghan/electron/miniTree_HGGCheck_WPTight_merged/UL2018/miniTree_GJets_*.root").Filter("abs(eleSCEta_Lead) < 1.4442");
    ROOT::RDF::RNode df_QCD     = ROOT::RDataFrame("miniTree", "/data4/chenghan/electron/miniTree_HGGCheck_WPTight_merged/UL2018/miniTree_QCD_*.root").Filter("abs(eleSCEta_Lead) < 1.4442");
    ROOT::RDF::RNode df_TT      = ROOT::RDataFrame("miniTree", "/data4/chenghan/electron/miniTree_HGGCheck_WPTight_merged/UL2018/miniTree_TT_*.root").Filter("abs(eleSCEta_Lead) < 1.4442");
    ROOT::RDF::RNode df_Zg      = ROOT::RDataFrame("miniTree", "/data4/chenghan/electron/miniTree_HGGCheck_WPTight_merged/UL2018/miniTree_Zg_*.root").Filter("abs(eleSCEta_Lead) < 1.4442");
    ROOT::RDF::RNode df_sig     = ROOT::RDataFrame("miniTree", "/data4/chenghan/electron/miniTree_HGGCheck_WPTight_merged/UL2018/miniTree_HDalitz_*_125_UL2018.root").Filter("abs(eleSCEta_Lead) < 1.4442");

    map<string, string> var_sel = {
        {"rho",                 "rho"},
        {"eleSCRawEn",          "eleSCRawEn_Lead"},
        {"eledEtaAtVtx",        "eledEtaAtVtx_Lead"},
        {"eledPhiAtVtx",        "eledPhiAtVtx_Lead"},
        {"elePtError",          "elePtError_Lead"},
        {"eleHoverE",           "eleHoverE_Lead"},
        {"eleEoverP",           "eleEoverP_Lead"},
        {"eleEoverPout",        "eleEoverPout_Lead"},
        {"eleEoverPInv",        "eleEoverPInv_Lead"},
        {"eleSCEtaWidth",       "eleSCEtaWidth_Lead"},
        {"eleSCPhiWidth",       "eleSCPhiWidth_Lead"},
        {"eleSigmaIEtaIEtaFull5x5", "eleSigmaIEtaIEtaFull5x5_Lead"},
        {"eleSigmaIPhiIPhiFull5x5", "eleSigmaIPhiIPhiFull5x5_Lead"},
        {"eleR9Full5x5",        "eleR9Full5x5_Lead"},
        {"eleBrem",             "eleBrem_Lead"},
        {"elePFChIso",          "elePFChIso_Lead"},
        {"elePFPhoIso",         "elePFPhoIso_Lead"},
        {"elePFNeuIso",         "elePFNeuIso_Lead"},
        {"eleGsfPtRatio",       "eleGsfPtRatio_Lead"},
        {"eleGsfDeltaR",        "eleGsfDeltaR_Lead"},
        {"eleGsfRelPtRatio",    "eleGsfRelPtRatio_Lead"}
        // {"gsf2Pt",               "gsf2.Pt()"},
        // {"gsf1Pt",               "gsf1.Pt()"},
        // {"gsfPtSum",             "gsf2.Pt() + gsf1.Pt()"},
        // {"digsfM",               "(gsf2+gsf1).M()"},
        // {"EOverdiTrkP",          "eleSCRawEn_Lead/(gsf2+gsf1).Pt()"},
        // {"gsfPtRatio",           "eleGsfPtRatio_Lead"},
        // {"gsf1D0",               "eleTrkD0_Lead"},
        // {"gsf2D0",               "eleSubTrkD0_Lead"},
        // {"gsf1Dz",               "eleTrkDz_Lead"},
        // {"gsf2Dz",               "eleSubTrkDz_Lead"},
        // {"eleSIP",               "eleSIP_Lead"},
        // {"gsf1Layers", "eleTrkLayers_Lead"},
        // {"gsf2Layers", "eleSubTrkLayers_Lead"},
        // {"higgsMass",            "H.M()"},
        // {"elePFCHIso",           "elePFChIso_Lead"}
        // {"mu1Eta",              "mu1.Eta()"},
        // {"mu2Eta",              "mu2.Eta()"},
        // {"phoEt",               "pho.Pt()"},
        // {"dRmumu",              "dRmumu"},
        // {"dRmax",               "dRmax"},
        // {"dRmin",               "dRmin"},
        // {"muPtRatio",           "muPtRatio"},
        // {"phoEtRatio",          "phoEtRatio"},
        // {"phoCorrR9Full5x5",    "phoCorrR9Full5x5"},
    };

    map<string, vector<double>> var_bins = {
        {"rho",                     {35, 0, 70}},
        {"eleSCRawEn",              {40, 10, 250}},
        {"eledEtaAtVtx",            {30, -0.015, 0.015}},
        {"eledPhiAtVtx",            {30, -0.05, 0.05}},
        {"elePtError",              {25, 0, 5}},
        {"eleHoverE",               {25, 0, 0.1}},
        {"eleEoverP",               {40, 0, 4}},
        {"eleEoverPout",            {40, 0, 4}},
        {"eleEoverPInv",            {25, -0.03, 0.02}},
        {"eleSCEtaWidth",           {30, 0.002, 0.02}},
        {"eleSCPhiWidth",           {40, 0, 0.08}},
        {"eleSigmaIEtaIEtaFull5x5", {40, 0.005, 0.015}},
        {"eleSigmaIPhiIPhiFull5x5", {40, 0.005, 0.03}},
        {"eleR9Full5x5",            {25, 0.5, 1}},
        {"eleBrem",                 {25, 0, 1}},
        {"elePFChIso",              {40, 0, 10}},
        {"elePFPhoIso",             {40, 0, 10}},
        {"elePFNeuIso",             {40, 0, 10}},
        {"eleGsfPtRatio",           {40, 0, 1}},
        {"eleGsfDeltaR",            {25, 0, 0.1}},
        {"eleGsfRelPtRatio",        {40, 0.2, 1.8}}
        // {"gsf2Pt",               {25, 0,  100}},
        // {"gsf1Pt",               {25, 0,  140}},
        // {"gsfPtSum",             {40, 44, 244}},
        // {"digsfM",               {25, 0, 1}},
        // {"EOverdiTrkP",          {20, 0, 5}},
        // {"gsf1D0",               {20, 0, 0.05}},
        // {"gsf2D0",               {20, 0, 0.05}},
        // {"gsf1Dz",               {20, 0, 0.1}},
        // {"gsf2Dz",               {20, 0, 0.1}},
        // {"gsfPtRatio",           {20, 0, 1}},
        // {"gsf1Layers", {20, 0, 20}},
        // {"gsf2Layers", {20, 0, 20}},
        // {"higgsMass",               {30, 110, 170}},
        // {"eleSIP",               {20, 0, 5}},
        // {"elePFCHIso",           {20, 0, 20}},
        // {"mu1Eta",              {60, -3, 3}},
        // {"mu2Eta",              {60, -3, 3}},
        // {"phoEt",               {35, 0, 140}},
        // {"dRmumu",              {80, 0, 4}},
        // {"dRmax",               {80, 0, 4}},
        // {"dRmin",               {80, 0, 4}},
        // {"muPtRatio",           {40, 0, 20}},
        // {"phoEtRatio",          {50, 0, 5}},
        // {"phoCorrR9Full5x5",    {50, 0, 1}},
    };

    map<string, string> var_legend = {
        {"rho",                     "#rho"},
        {"eleSCRawEn",              "raw E_{SC} (GeV)"},
        {"eledEtaAtVtx",            "#Delta#eta at vtx."},
        {"eledPhiAtVtx",            "#Delta#phi at vtx."},
        {"elePtError",              "p_{T} error (GeV)"},
        {"eleHoverE",               "H/E"},
        {"eleEoverP",               "E/P"},
        {"eleEoverPout",            "E/P_{out}"},
        {"eleEoverPInv",            "1/E - 1/P"},
        {"eleSCEtaWidth",           "#eta_{SC} width"},
        {"eleSCPhiWidth",           "#phi_{SC} width"},
        {"eleSigmaIEtaIEtaFull5x5", "#sigma_{i#etai#eta}"},
        {"eleSigmaIPhiIPhiFull5x5", "#sigma_{i#phii#phi}"},
        {"eleR9Full5x5",            "R_{9}"},
        {"eleBrem",                 "f_{brem}"},
        {"elePFChIso",              "Iso_{charged} (GeV)"},
        {"elePFPhoIso",             "Iso_{photon} (GeV)"},
        {"elePFNeuIso",             "Iso_{neutral} (GeV)"},
        {"eleGsfPtRatio",           "p^{Trk2}_{T} / p^{Trk1}_{T}"},
        {"eleGsfDeltaR",            "#DeltaR(Trk1, Trk2)"},
        {"eleGsfRelPtRatio",        "p^{diTrk}_{T}/E^{raw}_{SC}"}
        // {"gsf2Pt",               "p^{Trk2}_{T} (GeV)"},
        // {"gsf1Pt",               "p^{Trk1}_{T} (GeV)"},
        // {"gsfPtSum",             "p^{Trk2}_{T} + p^{Trk1}_{T} (GeV)"},
        // {"digsfM",               "m_{diTrk} (GeV)"},
        // {"EOverdiTrkP",          "E^{SC}_{raw} / p^{diTrk}_{T}"},
        // {"gsf1D0",               "dxy^{Trk1} (cm)"},
        // {"gsf2D0",               "dxy^{Trk2} (cm)"},
        // {"gsf1Dz",               "dz^{Trk1} (cm)"},
        // {"gsf2Dz",               "dz^{Trk2} (cm)"},
        // {"gsfPtRatio",           "p^{Trk2}_{T} / p^{Trk1}_{T}"},
        // {"gsf1Layers",           "expected inner hits of Trk1"},
        // {"gsf2Layers",           "expected inner hits of Trk2"},
        // {"higgsMass",            "m_{ee#gamma} (GeV)"},
        // {"elePFCHIso",           "Iso_{CH} (GeV)"},
        // {"eleSIP",               "Significance of IP"}
        // {"mu1Eta",              "#eta^{#mu^{1}}"},
        // {"mu2Eta",              "#eta^{#mu^{2}}"},
        // {"phoEt",               "E^{#gamma}_{T} (GeV)"},
        // {"dRmumu",              "#DeltaR(#mu^{1},#mu^{2})"},
        // {"dRmax",               "#DeltaR_{max}(#mu,#gamma)"},
        // {"dRmin",               "#DeltaR_{min}(#mu,#gamma)"},
        // {"muPtRatio",           "p^{#mu^{1}}_{T} / p^{#mu^{2}}_{T}"},
        // {"phoEtRatio",          "p^{#mu#mu}_{T} / p^{#gamma}_{T}"},
        // {"phoCorrR9Full5x5",    "Full5x5 R_{9}"},
    };

    for (auto it = var_sel.begin(); it != var_sel.end(); it++){
        const string varName        = it->first;
        const string sel            = it->second;
        const vector<double> bins   = var_bins[varName];
        const string legName        = var_legend[varName];
        draw_var(iBE, df_data, df_GJets, df_QCD, df_TT, df_Zg, df_sig, varName, sel, bins, legName);
    }
}