#ifndef WEIGHTHANDLER_H_
#define WEIGHTHANDLER_H_

#include <iostream>
#include <vector>
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TH2F.h"
#include "yaml-cpp/yaml.h"


// namespace wei{
//     // should be used after "ApplyWeight" and "GetScaleFactors"
//     ROOT::RDF::RNode GetFinalWeights(ROOT::RDF::RNode df, std::vector<std::string> column_SFs);
// }


class WeightHandler{
public:
    WeightHandler(){};
    ~WeightHandler();

    void Init(std::string name_, YAML::Node sfcfg_);
    float GetSFFromHgg(float sceta, float full5x5r9, float pt = 0);
    float GetSFErrFromHgg(float sceta, float full5x5r9, float pt = 0);
    float GetSFFromEGM(float x, float y);
    float GetSFErrFromEGM(float x, float y);

    // get the dataframe which has scale factor and it's error column
    ROOT::RDF::RNode GetRDF(ROOT::RDF::RNode df, std::string column, std::string column_err, std::vector<std::string> in_cloumns);

private:
    std::string name;
    YAML::Node sfcfg;
    
    // for egm sfs
    TFile* InFile = nullptr;
    TH2F* SFs = nullptr;

    // for hgg sfs 
    std::vector<float>              etabins;
    std::vector<std::vector<float>> r9bins;
    std::vector<std::vector<float>> values_2D; // 2D
    std::vector<std::vector<float>> uncert_2D; // 2D
    std::vector<float>              ptbins;                 // 3D
    std::vector<std::vector<std::vector<float>>> values_3D; // 3D
    std::vector<std::vector<std::vector<float>>> uncert_3D; // 3D
};



#endif