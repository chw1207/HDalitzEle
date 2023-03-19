#ifndef ENCALIBRATER_H_
#define ENCALIBRATER_H_

#include <vector>
#include <string>
#include "TFile.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"


class EnCalibrater{
public:
    EnCalibrater(){};

    /*
    sepcify the correction files to read
    --> assuming the first file is for EB and the second file is for EE
    */
    void SetCorrFiles(std::vector<std::string> corr_files);
    void SetScaleFiles(std::vector<std::string> scale_files);
    void SetSmearFiles(std::vector<std::string> smear_files);

    /*
    get the dataframe with corrected column
    */
    ROOT::RDF::RNode GetConvCorr(ROOT::RDF::RNode df, std::string column, std::string column_eta, std::string column_corr);
    ROOT::RDF::RNode GetScaleDif(ROOT::RDF::RNode df, std::string column, std::string column_eta, std::string column_corr, bool isMC);
    ROOT::RDF::RNode GetSmearDif(ROOT::RDF::RNode df, std::string column, std::string column_eta, std::string column_corr, bool isMC);

private:
    std::string _fcorrEB  = "";
    std::string _fcorrEE  = "";
    std::string _fscaleEB = "";
    std::string _fscaleEE = "";
    std::string _fsmearEB = "";
    std::string _fsmearEE = "";
};

#endif