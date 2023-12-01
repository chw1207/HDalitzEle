#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <filesystem>
#include "yaml-cpp/yaml.h"
#include "fmt/core.h"
#include "fmt/format.h"
#include "fmt/ranges.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TString.h"
#include "TSystem.h"
#include "Utilities.h" // utils::sigmaEff
#include "TColor.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "tdrstyle.h"
#include "CMS_lumi_modification.h"

std::vector<std::string> categories;

// approximate median significance (AMS): https://www.pp.rhul.ac.uk/~cowan/stat/notes/medsigNote.pdf (so-called combined significance)
std::map<std::string, double> calc_sigmaEff(ROOT::RDF::RNode df_sig_cat, ROOT::RDF::RNode df_bkg_cat, std::string mass_column, std::string weight_column, std::vector<int> mass_range, float range_percentage){
    std::map<std::string, double> effsig;
    auto arr = df_sig_cat.Take<float>(mass_column);
    std::vector<float> mass((*arr).begin(), (*arr).end());
    std::vector<float> sigma1 = utils::sigmaEff(mass, 0.683);
    std::vector<float> sigma2 = utils::sigmaEff(mass, range_percentage);

    std::string sig_range = fmt::format("{} > {} && {} < {}", mass_column, sigma2[0], mass_column, sigma2[1]);
    auto sig_stats = df_sig_cat.Stats(mass_column, weight_column);
    auto sig_stats_range = df_sig_cat.Filter(sig_range).Stats(mass_column, weight_column);
    const double sig_yield = sig_stats->GetW();
    const double sig_yield_range = sig_stats_range->GetW();
    effsig["Resolution"] = sigma1[2];
    effsig["Signal"] = sig_yield;
    effsig["Signal_2sigma"] = sig_yield_range;

    // non-resonant background is estimated from data sideband region scaled to the same mass window as signal(2 sigma)
    // still need to think a way to do a beetter estimation....
    std::string sideband_range = fmt::format("({} > {} && {} < 120) || ({} < {} && {} > 130)", mass_column, mass_range[0], mass_column, mass_column, mass_range[1], mass_column);
    auto sideband_stats = df_bkg_cat.Filter(sideband_range).Stats(mass_column);
    auto total_data_stats = df_bkg_cat.Stats(mass_column);
    const double total_data_yield = total_data_stats->GetN(); // total selected data evtens
    const double scale = (sigma2[1] - sigma2[0]) /*2sigma window*/ / ((120.-mass_range[0]) + (mass_range[1]-130.)); /*sideband width*/
    const double sigma2_yield = sideband_stats->GetN() * scale;
    effsig["Events"] = total_data_yield;
    effsig["Sideband_2sigma"] = sigma2_yield;
    effsig["Ratio"] = sig_yield_range*100/(sigma2_yield+sig_yield_range); // %
    effsig["AMS"] = TMath::Sqrt(2 * ( (sig_yield_range + sigma2_yield) * TMath::Log( 1 + (sig_yield_range/sigma2_yield) ) - sig_yield_range));
    return effsig;
}


void makeSummaryPlot(std::vector<std::map<std::string, double>> cat_summary, std::string outName){
    setTDRStyle();
    auto canv = new TCanvas("c", "c", 1500, 900);
    canv->cd();

    // make lumi text
    TLatex* ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.065);
    ltx->DrawLatex(0.1, 0.95, "CMS preliminary");
   
    // make effective sigma plot
    canv->cd();

    TLatex* ltx3 = new TLatex();
    ltx3->SetNDC();

    TPad* pad1 = new TPad("pad1", "pad1", 0, 0, 0.45, 1);
    pad1->SetTopMargin(0.08);
    pad1->SetLeftMargin(0.4);
    pad1->SetBottomMargin(0.12);
    pad1->SetGridx();    
    pad1->Draw();
    pad1->cd();
    int Ncat = cat_summary.size();
    auto h1 = new TH1F("h1", "", Ncat, 0, Ncat);
    for(int i = 0; i < Ncat; i++){
        h1->SetBinContent(i+1, cat_summary[i].at("Resolution"));
        h1->GetXaxis()->SetBinLabel(i+1, categories[i].c_str());
    }
    h1->GetXaxis()->SetLabelSize(0.06);
    h1->GetYaxis()->SetTitle("#sigma_{eff} (GeV)");
    h1->GetYaxis()->SetRangeUser(1, 3);
    h1->GetYaxis()->SetNdivisions(505);
    h1->GetYaxis()->SetTickSize(0.03);
    h1->GetYaxis()->SetTitleSize(0.055);
    h1->GetYaxis()->SetLabelSize(0.04);
    h1->GetYaxis()->SetTitleOffset(0.95);

    h1->SetBarWidth(0.4);
    h1->SetBarOffset(0.3);
    h1->SetFillColor(TColor::GetColor("#109F65"));
    h1->Draw("hbar");

    ltx3->SetTextFont(61);
    ltx3->SetTextSize(0.07);
    ltx3->DrawLatex(0.22, 0.925, "CMS");
    ltx3->SetTextFont(52);
    ltx3->SetTextSize(0.05);
    ltx3->DrawLatex(0.4, 0.93, "Preliminary");
    ltx3->SetTextFont(42);
    ltx3->SetTextSize(0.05);
    pad1->RedrawAxis();

    // make purity plot
    canv->cd();
    TPad* pad2 = new TPad("pad2", "pad2", 0.45, 0, 0.71, 1);
    pad2->SetTopMargin(0.08);
    pad2->SetLeftMargin(0.02);
    pad2->SetBottomMargin(0.12);
    pad2->SetGridx();
    pad2->Draw();    
    pad2->cd();
    auto h2 = new TH1F("h2", "", Ncat, 0, Ncat);
    for(int i = 0; i < Ncat; i++){
        h2->SetBinContent(i+1, cat_summary[i].at("Ratio"));
        h2->GetXaxis()->SetBinLabel(i+1, "");
    }
    h2->GetYaxis()->SetTitle("S/(S+B) in #pm 2#sigma_{eff} (%)");
    // h2->GetYaxis()->SetRangeUser(0, 5.5);
    // h2->GetYaxis()->SetNdivisions(505);
    h2->GetYaxis()->SetTickSize(0.02);
    h2->GetYaxis()->SetTitleSize(0.075);
    h2->GetYaxis()->SetLabelSize(0.065);
    h2->GetYaxis()->SetLabelOffset(-0.012);
    h2->GetYaxis()->SetTitleOffset(0.73);

    h2->SetBarWidth(0.4);
    h2->SetBarOffset(0.3);
    h2->SetFillColor(TColor::GetColor("#0174BE"));
    h2->Draw("hbar");
    pad2->RedrawAxis();

    // make sensitivity plot
    canv->cd();
    TPad* pad3 = new TPad("pad3", "pad3", 0.71, 0, 1, 1);
    pad3->SetTopMargin(0.08);
    pad3->SetLeftMargin(0.02);
    pad3->SetRightMargin(0.04);
    pad3->SetBottomMargin(0.12);
    pad3->SetGridx();    
    pad3->Draw();    
    pad3->cd();
    auto h3 = new TH1F("h3", "", Ncat, 0, Ncat);
    for(int i = 0; i < Ncat; i++){
        h3->SetBinContent(i+1, cat_summary[i].at("AMS"));
        h3->GetXaxis()->SetBinLabel(i+1, "");
    }
    h3->GetYaxis()->SetTitle("#sqrt{2((S+B)ln(1+S/B)-S)} in #pm 2#sigma_{eff}");
    // h3->GetYaxis()->SetRangeUser(0, 0.5);
    h3->GetYaxis()->SetNdivisions(505);
    h3->GetYaxis()->SetTickSize(0.02);
    h3->GetYaxis()->SetTitleSize(0.072);
    h3->GetYaxis()->SetLabelSize(0.065);
    h3->GetYaxis()->SetLabelOffset(-0.013);
    h3->GetYaxis()->SetTitleOffset(0.72);

    h3->SetBarWidth(0.4);
    h3->SetBarOffset(0.3);
    h3->SetFillColor(TColor::GetColor("#FFC436"));
    h3->Draw("hbar");
    
    ltx3->SetTextFont(42);
    ltx3->SetTextSize(0.08);
    ltx3->DrawLatex(0.43, 0.93, "138 fb^{-1} (13TeV)");

    pad3->RedrawAxis();

    TString outDir = gSystem->GetDirName(outName.c_str());
    if (!std::filesystem::exists(outDir.Data()))
        gSystem->Exec(Form("mkdir -p %s", outDir.Data()));
    canv->Print(outName.c_str()); 
}


int main(int argc, char** argv){
    if (argc != 2)
        throw std::invalid_argument("Please specify the config file for plotting");
    std::string config_path(argv[1]);

    ROOT::EnableImplicitMT();
    // load the config and mini tree files
    const YAML::Node cfg = YAML::LoadFile(config_path);
    auto mc_files = cfg["mcTree_path"].as<std::vector<std::string>>();
    ROOT::RDF::RNode df_sig = ROOT::RDataFrame("miniTree", mc_files);

    auto data_files = cfg["dataTree_path"].as<std::vector<std::string>>();
    ROOT::RDF::RNode df_bkg = ROOT::RDataFrame("miniTree", data_files);

    auto mass_range = cfg["significance"]["mass_range"].as<std::vector<int>>();
    auto mass_column = cfg["significance"]["mass_column"].as<std::string>();
    auto weight_column = cfg["significance"]["weight_column"].as<std::string>();
    auto category = cfg["significance"]["category"];
    auto range = cfg["range"].as<std::vector<float>>(); 

    // calculate the combined significance and format to a table
    TString all_info;
    for (size_t i = 0; i < range.size(); i++){
        TString info = TString::Format("%-20s|%10s|%10s|%10s|%10s|%10s|%10s|%12s|\n", "category", "AMS", Form("f%.1f[%%]", range[i]), "sigma[GeV]", "Nsig",  Form("Nsig%.1f", range[i]), "Nbkg", Form("Nbkg%.1f", range[i]));

        std::vector<std::map<std::string, double>> cat_summary;
        info += TString::Format("====================================================================================================\n");
        for (auto it = category.begin(); it != category.end(); it++){
            auto cat_name = it->first.as<std::string>();
            auto cat_filter = it->second.as<std::string>();
            auto df_sig_cat = df_sig.Filter(cat_filter);
            auto df_bkg_cat = df_bkg.Filter(cat_filter);

            auto effsig = calc_sigmaEff(df_sig_cat, df_bkg_cat, mass_column, weight_column, mass_range, range[i]/100);
            info += TString::Format("%-20s|%10.2f|%10.2f|%10.2f|%10.2f|%10.2f|%10d|%12.2f|\n", cat_name.c_str(), effsig["AMS"], effsig["Ratio"], effsig["Resolution"], effsig["Signal"], effsig["Signal_2sigma"], (int)effsig["Events"], effsig["Sideband_2sigma"]);
            categories.emplace_back(cat_name);
            cat_summary.emplace_back(effsig);
        }
        info += TString::Format("====================================================================================================\n");
        info += TString::Format("                                                                                                    \n");
        printf(info.Data());
        if (range[i] == 95.5)
            makeSummaryPlot(cat_summary, Form("../plots/CatSummary_%.1f.pdf", range[i]));
        categories.clear();

        all_info += info;
    }

    // save the results to the file
    auto outName = cfg["AMS_file"].as<std::string>();
    printf("Save table in: %s\n", outName.c_str());
    std::filesystem::path outpath = outName;
    if (!std::filesystem::exists(outpath.parent_path()))
        system(Form("mkdir -p %s", outpath.parent_path().c_str()));
    std::ofstream outfile(outName);
    if (outfile.is_open()){
        outfile << all_info.Data();
        outfile.close();

    }     
    else
        throw std::runtime_error(Form("Failed to open the file: %s", outName.c_str()));
    
    return 0;
}