#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include "boost/program_options.hpp"
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


void draw_var(ROOT::RDF::RNode df_mc, ROOT::RDF::RNode df_data, const std::string varName, const YAML::Node setting, const YAML::Node lumi_option, std::string plot_path){
    // preprocess mc dataframe
    auto nf_mc = df_mc.Filter(setting["selecs"].as<std::string>());
    if (!nf_mc.HasColumn(setting["column"].as<std::string>()))
        nf_mc = nf_mc.Define(varName, setting["column"].as<std::string>());
    else
        nf_mc = nf_mc.Alias(varName, setting["column"].as<std::string>());

    // preprocess data dataframe
    auto nf_data = df_data.Filter(setting["selecs"].as<std::string>());
    if (setting["addsel"])
        nf_data = nf_data.Filter(setting["addsel"].as<std::string>());
    if (!nf_data.HasColumn(setting["column"].as<std::string>()))
        nf_data = nf_data.Define(varName, setting["column"].as<std::string>());
    else
        nf_data = nf_data.Alias(varName, setting["column"].as<std::string>());

    auto bins = setting["numbin"].as<std::vector<float>>();
    auto h_mc = nf_mc.Histo1D({Form("h_%s_mc", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName, setting["mc_wei"].as<std::string>());
    if (setting["normalize"].as<bool>())
        h_mc->Scale(1./h_mc->Integral());
    else
        h_mc->Scale(setting["mc_scale"].as<int>());
    auto h_data = nf_data.Histo1D({Form("h_%s_data", varName.c_str()), "", (int)bins[0], bins[1], bins[2]}, varName);
    if (setting["normalize"].as<bool>())
        h_mc->Scale(1./h_data->Integral());

    setTDRStyle();
    auto canv = new TCanvas("c", "c", 800, 800);
    canv->cd();
    if (setting["isLogY"].as<bool>())
        canv->SetLogy();
    canv->SetTopMargin(0.06);
    canv->SetRightMargin(0.05);

    h_data->GetXaxis()->SetTitle(setting["x-axis"].as<std::string>().c_str());
    h_data->GetYaxis()->SetTitle(setting["y-axis"].as<std::string>().c_str());
    h_data->GetXaxis()->SetTitleSize(0.04);
    h_data->GetXaxis()->SetLabelSize(0.04);
    h_data->GetXaxis()->SetLabelOffset(0.015);
    h_data->GetXaxis()->SetTitleOffset(1.25);
    h_data->GetYaxis()->SetTitleSize(0.04);
    h_data->GetYaxis()->SetLabelSize(0.04);
    h_data->GetYaxis()->SetTitleOffset(2);
    if (setting["isLogY"].as<bool>())
        h_data->GetYaxis()->SetRangeUser(1, h_data->GetBinContent(h_data->GetMaximumBin()) * setting["y_scale"].as<float>());
    else
        h_data->GetYaxis()->SetRangeUser(0, h_data->GetBinContent(h_data->GetMaximumBin()) * setting["y_scale"].as<float>());

    h_data->SetMarkerColor(TColor::GetColor("#202020"));
    h_data->SetMarkerSize(1.2);
    h_data->SetMarkerStyle(20);
    h_data->SetLineColor(TColor::GetColor("#202020"));
    h_data->SetLineWidth(2);
    h_data->Draw("EP");

    h_mc->SetFillColor(TColor::GetColor(setting["mc_fcolor"].as<std::string>().c_str()));
    h_mc->SetLineColor(TColor::GetColor(setting["mc_lcolor"].as<std::string>().c_str()));
    h_mc->Draw("hist same");
    h_data->Draw("P same");
    canv->RedrawAxis();

    float text_posx = setting["latex_posx"].as<float>();
    float text_posy = setting["latex_posy"].as<float>();
    auto latex_text = setting["latex_text"].as<std::string>();
    auto ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.035);
    ltx->DrawLatex(text_posx, text_posy, latex_text.c_str());

    auto leg_posx = setting["leg_posx"].as<std::vector<float>>();
    auto leg_posy = setting["leg_posy"].as<std::vector<float>>();
    auto leg = new TLegend(leg_posx[0], leg_posy[0], leg_posx[1], leg_posy[1]);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->AddEntry(h_data.GetPtr(), "data", "LE1P");
    leg->AddEntry(h_mc.GetPtr(), Form("signal #times %d", setting["mc_scale"].as<int>()), "f");
    leg->Draw("same");

    TString lumi_text(lumi_option["lumi_text"].as<std::string>());
    int year = lumi_option["year"].as<int>();
    TString extra_text(lumi_option["extra_text"].as<std::string>());
    TString proc_text(lumi_option["proc_text"].as<std::string>());
    CMS_lumi((TCanvas*) canv, 5, 10, lumi_text, year, true, extra_text, proc_text, "");

    if (!std::filesystem::exists(plot_path)){
        system(Form("mkdir -p %s", plot_path.c_str()));
    }
    canv->Print(Form("%s/%s.pdf", plot_path.c_str(), varName.c_str()));

    delete canv;
    delete ltx;
    delete leg;
}


int main(int argc, char** argv){

    if (argc != 2)
        throw std::invalid_argument("Please specify the config file for plotting");
    std::string config_path(argv[1]);

    ROOT::EnableImplicitMT(10);

    const YAML::Node cfg = YAML::LoadFile(config_path);
    auto mc_files = cfg["mcTree_path"].as<std::vector<std::string>>();
    ROOT::RDF::RNode df_mc = ROOT::RDataFrame("miniTree", mc_files);

    auto data_files = cfg["dataTree_path"].as<std::vector<std::string>>();
    ROOT::RDF::RNode df_data = ROOT::RDataFrame("miniTree", data_files);

    auto plot_path = cfg["plot_path"].as<std::string>();
    const YAML::Node lumi_option = cfg["cms_lumi"];
    for (auto it = cfg["kinematics"].begin(); it != cfg["kinematics"].end(); it++){
        const auto varName = it->first.as<std::string>();
        const YAML::Node setting = it->second;
        draw_var(df_mc, df_data, varName, setting, lumi_option, plot_path);
    }

    return 0;
}