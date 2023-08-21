#include <algorithm>
#include <iostream>
#include <filesystem> // C++17
#include "Utilities.h"

std::vector<std::string> utils::find_files(std::string dirname){
    std::vector<std::string> v;
    for (const auto &entry : std::filesystem::directory_iterator(dirname)){
        if (!std::filesystem::is_directory(entry) && entry.path().extension() == ".root")
            v.push_back(entry.path());
    }

    if (v.size() == 0)
        throw std::runtime_error(Form("No files found in %s", dirname.c_str()));

    std::string fcount = v.size() == 1 ? "file" : "files";
    return v;
}


ROOT::RVec<int> utils::getIdx(const ROOT::RVec<int>& isgood, const ROOT::RVec<float>& pt){

    ROOT::RVec<int> idx_select = ROOT::VecOps::Nonzero(isgood);
    ROOT::RVec<int> idx_sort = ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(pt));
    ROOT::RVec<int> idx = ROOT::VecOps::Intersect(idx_sort, idx_select);

    return idx;
}


ROOT::RVec<ROOT::Math::PtEtaPhiMVector> utils::P4Vector(
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& phi,
    float m
){
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> v(pt.size());

    for (size_t i = 0; i < pt.size(); ++i){
        v[i].SetCoordinates(pt[i], eta[i], phi[i], m);
    }

    return v;
}


std::string utils::joinAND(const std::vector<std::string>& vstr){
    std::string sel = "";
    for (size_t i = 0; i < vstr.size(); i++){
        sel.append(vstr[i]);
        if (i < vstr.size()-1)
            sel.append(" && ");
    }

    return sel;
}


TString utils::printReport(ROOT::RDF::RResultPtr<ROOT::RDF::RCutFlowReport> cutflow){
    TString cutReport(" [+] Cut flow report\n");
    for (auto &&cutInfo : cutflow){
        const auto &name = cutInfo.GetName();
        const auto pass = cutInfo.GetPass();
        const auto all = cutInfo.GetAll();
        const auto eff = cutInfo.GetEff();
        cutReport += TString::Format("     - %-21s: pass = %-10lldall = %-10lldeff = %3.2f%% \n", name.c_str(), pass, all, eff);
    }
    printf(cutReport.Data());

    return cutReport;
}


void utils::printParameters(std::string config, std::string era, std::string nthreads, int max_events){
    TString pars("******************************************************\n");
    pars += TString::Format("                  Common parameters                   \n");
    pars += TString::Format("******************************************************\n");
    pars += TString::Format("config files: %s\n", config.c_str());
    pars += TString::Format("era: %s\n", era.c_str());
    pars += TString::Format("nthreads: %s\n", nthreads.c_str());
    pars += TString::Format("max events: %d (-1 means all)\n", max_events);
    pars += TString::Format("******************************************************\n");
    pars += TString::Format("\n");

    printf(pars.Data());
}


std::vector<int> utils::GetHumanTime(double time){
    std::vector<int> t(3, 0);
    if (time < 0)
        return t;
    else{
        double real = time;

        int hours = (int) real / 3600.;
        real -= hours * 3600.;
        t[0] = hours;

        int mins = (int) real / 60.;
        real -= mins * 60.;
        t[1] = mins;

        int sec = (int) real;
        t[2] = sec;
    }

    return t;
}


int utils::FindBins(std::vector<float> bin, const float var){
    std::sort(bin.begin(), bin.end()); // ensure the values are sorted from low to high 
    auto it = std::lower_bound(bin.begin(), bin.end(), var);
    int pos = std::distance(bin.begin(), it) - 1;
    int max_pos = bin.size() - 2;
    if ((pos < 0) || (pos > max_pos)) // underflow or overflow
        return -1;
    return pos;
}


std::vector<float> utils::sigmaEff(std::vector<float> v, const float threshold){
    std::sort(v.begin(), v.end());

    int total = v.size();
    int max = (int)(threshold * total);

    std::vector<float>  start;
    std::vector<float>  stop;
    std::vector<float>  width;

    int i = 0;
    while (i != total-1){
        int count = 0;
        int j = i;
        while (j != total-1 && count < max){
            ++count;
            ++j;
        }
        if (j != total-1){
            start.push_back(v[i]);
            stop .push_back(v[j]);
            width.push_back(v[j] - v[i]);
        }
        ++i;
    }

    float minwidth = *min_element(width.begin(), width.end()) * 0.5;
    auto pos = min_element(width.begin(), width.end()) - width.begin();
    float xmin = start[pos];
    float xmax = stop[pos];

    std::vector<float> sigma = {xmin, xmax, minwidth};
    return sigma;
}