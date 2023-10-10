#include "TMVA/Reader.h"
#include "TMVA/RReader.hxx"
#include <memory> // std::unique_ptr
#include <vector>
#include <string>

// The classical class TMVA::Reader is not thread-safe, so it cannot be used directly when EnableImplicitMT is turned on.
// ROOT had developed a thread-safe class TMVA::Experimental::RReader by adding a mutex on TMVA::Reader.
// However, in this case, we loose the advantage of using multi-threading.
// Therefore, I write a simple class by creating multiple TMVA::Readers and using one instance per thread to ensure thread safety based on TMVA::Experimental::RReader.
/*
    int num_threads = 10;
    TMVASafeReader reader("/path/to/xml_file", num_threads);

    ROOT::EnableImplicitMT(num_threads);
    auto df = ROOT::RDataFrame("tree", "/path/to/root_file")
    auto df_mva = df.DefineSlot("mva",  [&](unsigned int slot,
                                            std::vector<float>& features){
                                                return reader.Compute(features, slot)[0];
                                            }, {....});
*/

class TMVASafeReader{
public:
    TMVASafeReader(std::string path, int nslot=1);
    
    std::vector<float> Compute(const std::vector<float> &x, int slot=0);
    std::vector<std::string> GetVariableNames();

private:
    std::vector<std::unique_ptr<TMVA::Reader>> fReaders;
    unsigned int fNumClasses;
    TMVA::Experimental::Internal::AnalysisType fAnalysisType;
    std::vector<std::string> fExpressions;
    std::vector<std::string> fVariables;
    std::vector<std::vector<float>> fInputs;
    int fNslot;
    const char *name = "RReader";
};

