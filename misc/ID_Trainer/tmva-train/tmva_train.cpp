#include <iostream>
#include <string>
#include <vector>

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/DataSetInfo.h"
#include "TMVA/Config.h"
#include "TMVA/MethodDL.h"
#include "TMVA/Types.h"

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TStopwatch.h"
// DNN training using GPU
// usage: g++ tmva_train.cpp -o tmva_train -lTMVA -Wall $(root-config --glibs --cflags --libs)


void train(){
    auto input_file  = TFile::Open("/data4/chenghan/external/MergedID/Output_Merged2GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EE/XGB/train_tmva.root");
    auto output_file = TFile::Open("train_results_EE.root", "RECREATE");
    auto train_tree  = (TTree*) input_file->Get("train_tree");
    auto test_tree   = (TTree*) input_file->Get("test_tree");
    std::cout << "Classification          : Using input file: " << input_file->GetName() << std::endl;

    ROOT::EnableImplicitMT(20);

    // create factory
    auto factory    = new TMVA::Factory("MergedMVA_EE", output_file,
                                        "!V:!Silent:Color:DrawProgressBar:AnalysisType=multiclass");
    auto dataloader = new TMVA::DataLoader("dataset");
    std::vector<std::string> features = {
        "rho",
        "eleSCEta",
        "eleSCRawEn",

        "eledEtaAtVtx",
        "eledPhiAtVtx",
        "elePtError",
        "eleHoverE",
        "eleEoverP",
        "eleEoverPout",
        "eleEoverPInv",

        "eleSCEtaWidth",
        "eleSCPhiWidth",
        "eleSigmaIEtaIEtaFull5x5",
        "eleSigmaIPhiIPhiFull5x5",
        "eleR9Full5x5",
        "eleBrem",

        "elePFChIso",
        "elePFPhoIso",
        "elePFNeuIso",

        "gsfPtRatio",
        "gsfDeltaR",
        "gsfRelPtRatio",

        // "eleIDMVAIso"
    };
    for (size_t i = 0; i < features.size(); i++){
        dataloader->AddVariable(features[i].c_str());
    }

    // prepare training data
    // https://root-forum.cern.ch/t/how-to-specify-different-input-files-for-training-and-testing-events/33549/2
    dataloader->AddTree(train_tree, "sig0", 1.0, "EleType == 0", TMVA::Types::kTraining);
    dataloader->AddTree(train_tree, "bkg1", 1.0, "EleType == 1", TMVA::Types::kTraining);
    dataloader->AddTree(train_tree, "bkg2", 1.0, "EleType == 2", TMVA::Types::kTraining);
    dataloader->AddTree(test_tree,  "sig0", 1.0, "EleType == 0", TMVA::Types::kTesting);
    dataloader->AddTree(test_tree,  "bkg1", 1.0, "EleType == 1", TMVA::Types::kTesting);
    dataloader->AddTree(test_tree,  "bkg2", 1.0, "EleType == 2", TMVA::Types::kTesting);
    dataloader->SetWeightExpression("NewWt", "sig0");
    dataloader->SetWeightExpression("NewWt", "bkg1");
    dataloader->SetWeightExpression("NewWt", "bkg2");
    // TCut mycuts = "";
    // TCut mycutb = "";
    // std::string pre_opt = "SplitMode=Random:SplitSeed=42:NormMode=NumEvents:!V";
    // dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, pre_opt);
    
    // Deep Neural Network
    TString layout("Layout=RELU|64,RELU|128,RELU|128,RELU|64,LINEAR"); // DENSE
    // one can catenate several training strategies
    TString train_strategy( "TrainingStrategy=LearningRate=1e-5,Momentum=0.2,Repetitions=1,ConvergenceSteps=50,"
                            "BatchSize=50000,TestRepetitions=1,WeightDecay=1e-6,DropConfig=0.05,"
                            "MaxEpochs=2000,"
                            "Regularization=L2," // regularization is used to avoid overfitting
                            "Optimizer=Adam"
                          );
    TString dnn_opts("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=G,D,N:WeightInitialization=XAVIER:Architecture=GPU");
    dnn_opts.Append(":");
    dnn_opts.Append(layout);
    dnn_opts.Append(":");
    dnn_opts.Append(train_strategy);
    factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN", dnn_opts);

    // Gradient Boost Decision Tree
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG",
                        "!H:!V:NTrees=1500:BoostType=Grad:Shrinkage=0.01:nCuts=500:MaxDepth=7:"
                        // "UseBaggedGrad:BaggedSampleFraction=0.5:"
                        "PruneMethod=CostComplexity:PruneStrength=50"
                        );

    // trainning
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    output_file->Close();
    input_file->Close();

    delete factory;
    delete dataloader;
}


int main(int argc, char** argv){
    TStopwatch time;
    time.Start();

    train();

    time.Stop();
    time.Print();

    return 0;
}