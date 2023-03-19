#ifndef WEIGHTHANDLER_H_
#define WEIGHTHANDLER_H_

#include <iostream>
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "yaml-cpp/yaml.h"


namespace wei{
    ROOT::RDF::RNode ApplyWeight(ROOT::RDF::RNode df, const YAML::Node cfg);

    ROOT::RDF::RNode GetScaleFactors(ROOT::RDF::RNode df, const YAML::Node cfg);

    // should be used after "ApplyWeight" and "GetScaleFactors"
    ROOT::RDF::RNode GetFinalWeights(ROOT::RDF::RNode df);
}

#endif