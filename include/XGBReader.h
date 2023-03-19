#ifndef XGBREADER_H_
#define XGBREADER_H_

#include <iostream>
#include <vector>
// #include "TMVA/RTensor.hxx"
#include "xgboost/c_api.h"


class XGBReader{
public:
    XGBReader(std::string _fName);
    virtual ~XGBReader();

    std::vector<float> Compute(const std::vector<float>& features);
    // std::vector<std::vector<float>> Compute(std::vector<std::vector<float>>& features, bool verbose=false);
    // TMVA::Experimental::RTensor<float> Compute(TMVA::Experimental::RTensor<float>& features, bool verbose=false);

    std::vector<int> argMax(const std::vector<std::vector<float>>& scores);

    uint64_t GetNumFeature();

private:
    /* Row-major matrix */
    // struct _Matrix{
    //     float *data;
    //     size_t shape[2];
    // };

    // /* A custom data type for demo. */
    // typedef struct _Matrix *Matrix;

    // void Matrix_Create(Matrix *self, float const *data, size_t n_samples, size_t n_features);
    // void Matrix_Free(Matrix self);
    void safe_xgboost(uint64_t call);

    BoosterHandle booster;
    // DMatrixHandle dpred;

    // c_json_config
    const char* config = "{\"training\": false, \"type\": 0, \"iteration_begin\": 0, \"iteration_end\": 0, \"strict_shape\": true}";

    // prediction-results initialization
    // uint64_t out_dim = 0;
    // uint64_t const *out_shape = NULL;
    // float const *out_results = NULL;
};
#endif // XGBREADER_H_