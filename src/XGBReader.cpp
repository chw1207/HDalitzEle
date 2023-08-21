// #include <stddef.h>
// #include <stdlib.h>
#include <string>
#include <algorithm>
#include "TString.h"
#include "XGBReader.h"



XGBReader::XGBReader(std::string _fName){
    safe_xgboost(XGBoosterCreate(NULL, 0, &booster));
    safe_xgboost(XGBoosterSetParam(booster, "nthread", "1"));
    safe_xgboost(XGBoosterLoadModel(booster, _fName.c_str()));
}


XGBReader::XGBReader(){
    safe_xgboost(XGBoosterCreate(NULL, 0, &booster));
    safe_xgboost(XGBoosterSetParam(booster, "nthread", "1"));
}


XGBReader::~XGBReader(){
    // safe_xgboost(XGDMatrixFree(dpred));
    safe_xgboost(XGBoosterFree(booster));
}


void XGBReader::Init(std::string _fName){
    safe_xgboost(XGBoosterLoadModel(booster, _fName.c_str()));
}


std::vector<float> XGBReader::Compute(const std::vector<float>& features){
    DMatrixHandle dpred;
    uint64_t out_dim;
    uint64_t const *out_shape;
    float const *out_results;

    safe_xgboost(XGDMatrixCreateFromMat((float*)features.data(), 1, features.size(), -99999., &dpred));
    safe_xgboost(XGBoosterPredictFromDMatrix(booster, dpred, config, &out_shape, &out_dim, &out_results));
    // printf("dim: %lu, shape[0]: %lu, shape[1]: %lu\n", out_dim, out_shape[0], out_shape[1]);

    std::vector<float> score(out_results, out_results + out_shape[1]);
    safe_xgboost(XGDMatrixFree(dpred));

    return score;
}


// std::vector<std::vector<float>> XGBReader::Compute(const std::vector<std::vector<float>>& features, bool verbose/*=false*/){
//     const auto numEntries = features.size();
//     const auto numVars = features[0].size();

//     // fill the row-major matrix to perform the prediction
//     Matrix mfeatures;
//     Matrix_Create(&mfeatures, NULL, numEntries, numVars);
//     int k = 0;
//     for (size_t i = 0; i < numEntries; i++){
//         for (size_t j = 0; j < numVars; j++){
//             mfeatures->data[k] = features[i][j];
//             k++;
//         }
//     }

//     // perform the prediction
//     safe_xgboost(XGDMatrixCreateFromMat(mfeatures->data, numEntries, numVars, -99999., &dpred));
//     if (verbose == true){
//         bst_ulong row, col;
//         safe_xgboost(XGDMatrixNumRow(dpred, &row));
//         safe_xgboost(XGDMatrixNumCol(dpred, &col));
//         printf("[INFO] Vector is converted to %lux%lu DMatrix.\n", row, col);
//     }
//     safe_xgboost(XGBoosterPredictFromDMatrix(booster, dpred, config, &out_shape, &out_dim, &out_results));

//     // save the results to 2D vector
//     std::vector<std::vector<float>> score(numEntries, std::vector<float>(out_shape[1], 0.));
//     for (size_t i = 0; i < numEntries; i++){ // row
//         for (size_t j = 0; j < out_shape[1]; j++){ // column
//             score[i][j] = out_results[i * out_shape[1] + j];
//         }
//     }

//     // release the memory
//     Matrix_Free(mfeatures);

//     return score;
// }


// TMVA::Experimental::RTensor<float> XGBReader::Compute(const TMVA::Experimental::RTensor<float>& features, bool verbose/*=false*/){
//     const auto shape = features.GetShape();
//     const auto numEntries = shape[0];
//     const auto numVars = shape[1];

//     // fill the row-major matrix to perform the prediction
//     Matrix mfeatures;
//     Matrix_Create(&mfeatures, NULL, numEntries, numVars);
//     int k = 0;
//     for (size_t i = 0; i < numEntries; i++){
//         for (size_t j = 0; j < numVars; j++){
//             mfeatures->data[k] = features(i, j);
//             k++;
//         }
//     }

//     // perform the prediction
//     safe_xgboost(XGDMatrixCreateFromMat(mfeatures->data, numEntries, numVars, -99999., &dpred));
//     if (verbose == true){
//         bst_ulong row, col;
//         safe_xgboost(XGDMatrixNumRow(dpred, &row));
//         safe_xgboost(XGDMatrixNumCol(dpred, &col));
//         printf("[INFO] RTensor is converted to %lux%lu DMatrix.\n", row, col);
//     }
//     safe_xgboost(XGBoosterPredictFromDMatrix(booster, dpred, config, &out_shape, &out_dim, &out_results));

//     // save the results to RTensor
//     TMVA::Experimental::RTensor<float> score({numEntries, out_shape[1]});
//     for (size_t i = 0; i < numEntries; i++){ // row
//         for (size_t j = 0; j < out_shape[1]; j++){ // column
//             score(i, j) = out_results[i * out_shape[1] + j];
//         }
//     }

//     // release the memory
//     Matrix_Free(mfeatures);

//     return score;
// }


std::vector<int> XGBReader::argMax(const std::vector<std::vector<float>>& scores){
    std::vector<int> idx;
    for(size_t i = 0; i < scores.size(); i++){
        idx.push_back(std::max_element(scores[i].begin(), scores[i].end()) - scores[i].begin());
    }

    return idx;
}


uint64_t XGBReader::GetNumFeature(){
    uint64_t nf;
    safe_xgboost(XGBoosterGetNumFeature(booster, &nf));

    return nf;
}


// void XGBReader::Matrix_Create(Matrix *self, float const *data, size_t n_samples, size_t n_features){
//     if (self == NULL) {
//         fprintf(stderr, "Invalid pointer to %s\n", __func__);
//         exit(-1);
//     }

//     *self = (Matrix)malloc(sizeof(struct _Matrix));
//     (*self)->data = (float *)malloc(n_samples * n_features * sizeof(float));
//     (*self)->shape[0] = n_samples;
//     (*self)->shape[1] = n_features;

//     if (data != NULL)
//         memcpy((*self)->data, data, (*self)->shape[0] * (*self)->shape[1] * sizeof(float));
// }


// void XGBReader::Matrix_Free(Matrix self){
//     if (self != NULL){
//         if (self->data != NULL){
//             self->shape[0] = 0;
//             self->shape[1] = 0;
//             free(self->data);
//             self->data = NULL;
//         }
//         free(self);
//     }
// }


void XGBReader::safe_xgboost(uint64_t call){
    if (call != 0)
        throw std::runtime_error(XGBGetLastError());
}