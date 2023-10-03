#include <iostream>
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"

// a simple class to read the root file containing the center and scale of RobustScaler
// and then do the transformation
// do not have the function "fit"!

template <typename T>
void print_vector(const vector<T>& vec){
    std::cout << "{";
    for (size_t i = 0; i < vec.size(); i++){
        if (i == vec.size()-1)
            std::cout << vec[i];
        else
            std::cout << vec[i] << ", ";
    }
    std::cout << "}" << std::endl;
}


class MyRobustScaler{
public:
    MyRobustScaler(std::string filename);
    ~MyRobustScaler(){};

    std::vector<float> Transform(const std::vector<float>& input);
    void Print();

private:
    std::vector<double> center_vec; // The median value for each feature in the training set.
    std::vector<double> scale_vec;  // The (scaled) interquartile range for each feature in the training set.
};


MyRobustScaler::MyRobustScaler(std::string filename){
    auto file = new TFile(filename.c_str(), "READ");
    auto tree = (TTree*) file->Get("RobustScaler");

    double center_, scale_;
    tree->SetBranchAddress("center", &center_);
    tree->SetBranchAddress("scale", &scale_);

    for (int i = 0; i < tree->GetEntries(); i++){
        tree->GetEntry(i);
        center_vec.push_back(center_);
        scale_vec.push_back(scale_);
    }
    
    file->Close();
    delete file;
}


std::vector<float> MyRobustScaler::Transform(const std::vector<float>& input){
    std::vector<float> output(input.size());
    for (size_t i = 0; i < input.size(); i++) {
        output[i] = (input[i] - center_vec[i]) / scale_vec[i];
    }
    return output;
}


void MyRobustScaler::Print(){
    std::cout << "center_ is the median value for each feature in the training set." << std::endl;
    std::cout << "scale_ is the (scaled) interquartile range for each feature in the training set." << std::endl;
    std::cout << "center_: ";
    print_vector(center_vec);

    std::cout << "scale_: ";
    print_vector(scale_vec);
}