#ifndef MERGEDENCORRECTOR_H_
#define MERGEDENCORRECTOR_H_
#include <vector>
#include <string>
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"


ROOT::RDF::RNode doHDALRegression(ROOT::RDF::RNode df, std::vector<std::string> path);

ROOT::RDF::RNode doHDALRegressionXGB(ROOT::RDF::RNode df, std::vector<std::string> path);
// ROOT::RDF::RNode doHDALRegressionV2(ROOT::RDF::RNode df, std::vector<std::string> path);

#endif