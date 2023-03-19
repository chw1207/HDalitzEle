#ifndef PHOR9SSCORRECTION_H_
#define PHOR9SSCORRECTION_H_

#include <map>
#include <string>
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"


ROOT::RDF::RNode doSSCorrections(ROOT::RDF::RNode df, std::string corrColumn, std::string era, std::map<std::string, std::string> weightMap);

#endif