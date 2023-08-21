#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <string>
#include <vector>
#include "TString.h"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDF/RCutFlowReport.hxx"

namespace utils{
    std::vector<std::string> find_files(std::string dirname);


    // Function to get the index vector sorted by pT
    // Reference: https://root.cern/doc/master/vo006__IndexManipulation_8C.html
    ROOT::RVec<int> getIdx(const ROOT::RVec<int>& isgood, const ROOT::RVec<float>& pt);


    // Function to construct RVec of LorentzVector
    ROOT::RVec<ROOT::Math::PtEtaPhiMVector> P4Vector(
        const ROOT::RVec<float>& pt,
        const ROOT::RVec<float>& eta,
        const ROOT::RVec<float>& phi,
        float m
    );


    // concatenate strings(selections) with '&&'
    // example: vector<string> vstr = {"pt > 30", "mass < 50"} -> string str = "pt > 30 && mass < 50"
    std::string joinAND(const std::vector<std::string>& vstr);


    // print cut report gracefully and return the report as TString
    TString printReport(ROOT::RDF::RResultPtr<ROOT::RDF::RCutFlowReport> cutflow);

    // print parameters gracefully
    void printParameters(std::string config, std::string era, std::string nthreads, int max_events);

    // return vector -> {hours, mins, seconds}
    std::vector<int> GetHumanTime(double time);

    // find which bin a value fall into
    // useful when there are many else if statements
    // example:
    //      Origin:
    //         float aaa = 0.
    //         if (fabs(eta) >= 0. && fabs(eta) < 1.479)
    //              aaa = somthing_1
    //         else if (fabs(eta >= 1.479 && fabs(eta) < 2.5)
    //              aaa = somthing_2
    //      Using FindBins:
    //          vector<float> someting = {somthing_1, somthing_2};
    //          vector<float> etas = {0, 1.479, 2.5};
    //          float aaa = someting[FindBins(etas, fabs(eta))];
    // Note: FindBins return -1 when underflow and overflow
    int FindBins(std::vector<float> bin, const float var);


    // calculate the effective sigma of a vector under a given threshold
    // return vector -> {lower_bound, upper_bound, sigma}
    std::vector<float> sigmaEff(std::vector<float> v, const float threshold);
}
#endif