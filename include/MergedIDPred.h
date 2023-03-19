#ifndef MERGEDIDPRED_H_
#define MERGEDIDPRED_H_

#include <map>
#include <string>
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"


ROOT::RDF::RNode MergedIDPred(ROOT::RDF::RNode df, std::string predColumn, std::map<std::string, std::string> modelMap);

#endif