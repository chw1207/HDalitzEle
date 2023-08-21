#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
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
#include "Utilities.h" // utils::sigmaEff


// approximate median significance (AMS): https://www.pp.rhul.ac.uk/~cowan/stat/notes/medsigNote.pdf (so-called combined significance)
std::map<std::string, double> calc_sigmaEff(ROOT::RDF::RNode df_sig_cat, ROOT::RDF::RNode df_bkg_cat, std::string mass_column, std::string weight_column, std::vector<int> mass_range){
    std::map<std::string, double> effsig;
    auto arr = df_sig_cat.Take<float>(mass_column);
    std::vector<float> mass((*arr).begin(), (*arr).end());
    std::vector<float> sigma1 = utils::sigmaEff(mass, 0.683);
    std::vector<float> sigma2 = utils::sigmaEff(mass, 0.955);

    std::string sig_range = fmt::format("{} > {} && {} < {}", mass_column, sigma2[0], mass_column, sigma2[1]);
    auto sig_stats = df_sig_cat.Stats(mass_column, weight_column);
    auto sig_stats_range = df_sig_cat.Filter(sig_range).Stats(mass_column, weight_column);
    const double sig_yield = sig_stats->GetW();
    const double sig_yield_range = sig_stats_range->GetW();
    effsig["Resolution"] = sigma1[2];
    effsig["Signal"] = sig_yield;
    effsig["Signal_2sigma"] = sig_yield_range;

    // non-resonant background is estimated from data sideband region scaled to the same mass window as signal(2 sigma)
    std::string sideband_range = fmt::format("({} > {} && {} < {}) || ({} < {} && {} > {})", mass_column, mass_range[0], mass_column, sigma2[0], mass_column, mass_range[1], mass_column, sigma2[1]);
    auto sideband_stats = df_bkg_cat.Stats(mass_column);
    auto sideband_stats_range = df_bkg_cat.Filter(sideband_range).Stats(mass_column);
    const double sideband_yield = sideband_stats->GetN();
    const double sideband_yield_range = sideband_stats_range->GetN() * (sigma2[1] - sigma2[0])/(mass_range[1] - mass_range[0]);
    effsig["Sideband"] = sideband_yield;
    effsig["Sideband_2sigma"] = sideband_yield_range; // %
    effsig["Ratio"] = sig_yield_range*100/(sideband_yield_range+sig_yield_range);
    effsig["AMS"] = sqrt(2 * ( (sig_yield_range + sideband_yield_range) * log( 1 + (sig_yield_range/sideband_yield_range) ) - sig_yield_range));

    return effsig;
}


int main(int argc, char** argv){
    if (argc != 2)
        throw std::invalid_argument("Please specify the config file for plotting");
    std::string config_path(argv[1]);

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

    // calculate the combined significance and format to a table
    TString info = TString::Format("%-20s|%10s|%10s|%10s|%10s|%10s|%10s|%12s|\n", "category", "AMS", "f95.5[%]","sigma[GeV]", "Nsig", "Nsig95.5", "Nbkg", "Nbkg95.5");
    info += TString::Format("====================================================================================================\n");
    for (auto it = category.begin(); it != category.end(); it++){
        auto cat_name = it->first.as<std::string>();
        auto cat_filter = it->second.as<std::string>();
        auto df_sig_cat = df_sig.Filter(cat_filter);
        auto df_bkg_cat = df_bkg.Filter(cat_filter);

        auto effsig = calc_sigmaEff(df_sig_cat, df_bkg_cat, mass_column, weight_column, mass_range);
        info += TString::Format("%-20s|%10.2f|%10.2f|%10.2f|%10.2f|%10.2f|%10d|%12.2f|\n", cat_name.c_str(), effsig["AMS"], effsig["Ratio"], effsig["Resolution"], effsig["Signal"], effsig["Signal_2sigma"], (int)effsig["Sideband"], effsig["Sideband_2sigma"]);
    }
    printf(info.Data());

    // save the results to the file
    auto outName = cfg["AMS_file"].as<std::string>();
    printf("Save table in: %s\n", outName.c_str());
    std::filesystem::path outpath = outName;
    if (!std::filesystem::exists(outpath.parent_path()))
        system(Form("mkdir -p %s", outpath.parent_path().c_str()));
    std::ofstream outfile(outName);
    if (outfile.is_open()){
        outfile << info.Data();
        outfile.close();
    }
    else
        throw std::runtime_error(Form("Failed to open the file: %s", outName.c_str()));

    return 0;
}