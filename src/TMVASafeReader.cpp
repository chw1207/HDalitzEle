#include "TMVASafeReader.h"
#include "TString.h"


TMVASafeReader::TMVASafeReader(std::string path, int nslot/*=1*/){
    // Load config
    auto c = TMVA::Experimental::Internal::ParseXMLConfig(path);
    fVariables = c.variables;
    fExpressions = c.expressions;
    fAnalysisType = c.analysisType;
    fNumClasses = c.numClasses;
    fNslot = (nslot < 1) ? 1 : nslot;

    // Setup reader
    const int numVars = fVariables.size();
    fReaders.resize(fNslot);
    fInputs.resize(fNslot, std::vector<float>(numVars));
    for (int i = 0; i < fNslot; ++i) {
        fReaders[i] = std::make_unique<TMVA::Reader>("!Color:Silent");
        for (int j = 0; j < numVars; j++){
            fReaders[i]->AddVariable(TString(fExpressions[j]), &fInputs[i][j]);
        }
        fReaders[i]->BookMVA(name, path);
    }
}


std::vector<float> TMVASafeReader::Compute(const std::vector<float> &x, int slot/*=0*/){
    if (x.size() != fVariables.size())
        throw std::runtime_error("Size of input vector is not equal to number of variables.");

    if (slot > (fNslot-1))
        throw std::runtime_error("Slot number is larger than the total number of slots - 1.");

    if (slot < 0)
        throw std::runtime_error("Negative slot number is not reasonable.");

    // Copy over inputs to memory used by TMVA reader
    for (size_t i = 0; i < x.size(); i++) {
        fInputs[slot][i] = x[i];
    }

    // Evaluate TMVA model
    // Classification
    if (fAnalysisType == TMVA::Experimental::Internal::AnalysisType::Classification) {
        return std::vector<float>({static_cast<float>(fReaders[slot]->EvaluateMVA(name))});
    }
    // Regression
    else if (fAnalysisType == TMVA::Experimental::Internal::AnalysisType::Regression) {
        return fReaders[slot]->EvaluateRegression(name);
    }
    // Multiclass
    else if (fAnalysisType == TMVA::Experimental::Internal::AnalysisType::Multiclass) {
        return fReaders[slot]->EvaluateMulticlass(name);
    }
    // Throw error
    else {
        throw std::runtime_error("TMVASafeReader has undefined analysis type.");
        return std::vector<float>();
    }
}


std::vector<std::string> TMVASafeReader::GetVariableNames() { return fVariables; }