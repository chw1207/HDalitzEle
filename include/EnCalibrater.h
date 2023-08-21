#ifndef ENCALIBRATER_H_
#define ENCALIBRATER_H_

#include <vector>
#include <string>
#include <random> 
#include "TFile.h"
#include "TRandom3.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"


class EnCalibrater{
public:
    EnCalibrater(){};
    ~EnCalibrater();

    /*
    sepcify the correction files to read
    --> assuming the first file is for EB and the second file is for EE
    */
    void SetScaleFiles(std::vector<std::string> scale_files);
    void SetSmearFiles(std::vector<std::string> smear_files);

    /*
    get the dataframe with corrected column
    */
    ROOT::RDF::RNode GetScaleRDF(ROOT::RDF::RNode df, bool isMC, std::string column, std::string column_eta, std::string column_corr);
    ROOT::RDF::RNode GetSmearRDF(ROOT::RDF::RNode df, bool isMC, std::string column, std::string column_eta, std::string column_corr);

private:
    TH1F* getHist(std::string histFile, std::string histName);

    // for residual correction of energy scale 
    bool hasScale = false;
    TH1F* hscaleEB = nullptr;
    TH1F* hscaleEE = nullptr;
    
    // for residual correction of energy resolution 
    std::random_device rd;
    bool hasSmear = false;
    TH1F* hsmearEB = nullptr;
    TH1F* hsmearEE = nullptr;
    
    std::vector<std::shared_ptr<TRandom3>> rnds;
};

#endif