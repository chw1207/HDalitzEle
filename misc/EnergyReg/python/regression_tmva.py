import os, sys
import ROOT
import time
from datetime import datetime
from colorPrint import *
from argparse import ArgumentParser
ROOT.TMVA.Tools.Instance()
ROOT.TMVA.PyMethodBase.PyInitialize()   


def get_parser():
    parser = ArgumentParser(description="Script for merged electron energy regression using TMVA")
    parser.add_argument("-r",  "--region", help="Which region to run. Options: [EB, EE]", default="", type=str)
    parser.add_argument("-e",  "--ext", help="External text pass to plot name and model name", default="", type=str)
    parser.add_argument("--fast", help="Train the model without evalution. only model will be produce", default=False, action="store_true")

    return parser


def train_TMVA():
    # clean the working directory
    outDir = "../reg_results/TMVARegression{}_{}".format(args.ext, args.region)
    if os.path.exists(outDir):
        os.system("rm -r {}".format(outDir))
    os.makedirs(outDir)
    
    # enable multi threadings
    ROOT.EnableImplicitMT(30)
    
    # setup the input and output file
    input_file  = ROOT.TFile("../reg_signal.root", "read")
    output_file = ROOT.TFile("{}/TMVAResults.root".format(outDir), "RECREATE")
    print("TMVARegression           : Using input file: {}".format(input_file.GetName()), flush=True)
    
    # create factory
    factory     = ROOT.TMVA.Factory("TMVARegression", output_file, 
                                    "V:!Silent:Color:DrawProgressBar:AnalysisType=Regression")
    dataloader  = ROOT.TMVA.DataLoader(outDir)
    
    features = [
        "rho",
        "nVtx",
        "eleSCEta_Lead",
        "eleSCPhi_Lead",
        "eleSCRawEn_Lead",
        # "elePt_Lead",
        "eleCalibPt_Lead",

        "eledEtaAtVtx_Lead",
        "eledPhiAtVtx_Lead",
        "elePtError_Lead",
        "eleHoverE_Lead",
        "eleEoverP_Lead",
        "eleEoverPout_Lead",
        "eleEoverPInv_Lead",

        "eleSCEtaWidth_Lead",
        "eleSCPhiWidth_Lead",
        "eleSigmaIEtaIEtaFull5x5_Lead",
        "eleSigmaIPhiIPhiFull5x5_Lead",
        "eleR9Full5x5_Lead",
        "eleBrem_Lead",

        "gsfPtSum_Lead",
        "gsfPtRatio_Lead",
        "diTrk.Pt()",
        "gsfDeltaR_Lead"
    ]
    if args.region == "EE":
        features.append("eleESEnToRawE_Lead")
        
    for f in features:
        dataloader.AddVariable(f)
    dataloader.AddTarget("diGenEle.Pt()/elePt_Lead", "target", "")
    
    # load the ttree
    regTree = input_file.Get("miniTree")
    # regTree  = inTree.CloneTree(10000)
    dataloader.AddRegressionTree(regTree, 1.)
    dataloader.SetWeightExpression("instwei * puwei", "Regression")
    
    # train test split
    # the efficiency of presel is roughly 30%
    # SplitMode=Random:SplitSeed=123456
    region_cut = "fabs(eleSCEta_Lead) < 1.479" if args.region == "EB" else "fabs(eleSCEta_Lead) >= 1.479"
    presel = ROOT.TCut("elePresel_Lead == 1 && eleCalibPt_Lead > 25 && category == 2 && {} && (diGenEle.Pt()/elePt_Lead) < 1.2 && (diGenEle.Pt()/elePt_Lead) > 0.8".format(region_cut))
    dataloader.PrepareTrainingAndTestTree(presel,
                                          "nTrain_Regression=0:nTest_Regression=0:SplitMode=Alternate:NormMode=None:!CalcCorrelations:V")
    
    # book the training methods
    # BDTG
    # factory.BookMethod( dataloader, ROOT.TMVA.Types.kBDT, "BDTG", 
    #                     ROOT.TString("!H:V:SeparationType=RegressionVariance:BoostType=Grad:"
    #                                  "NTrees=1500:nCuts=40:MaxDepth=7:Shrinkage=0.1:"
    #                                  "MinNodeSize=0.1%:"
    #                                  "UseBaggedBoost:BaggedSampleFraction=0.5:"
    #                                  "RegressionLossFunctionBDTG=AbsoluteDeviation"))
    
    # DNN
    # ROOT installed by conda do not support GPU
    Layout   = "Layout=RELU|128,RELU|512,RELU|512,RELU|512,LINEAR"
    Strategy = ROOT.TString("TrainingStrategy=LearningRate=1e-5,Momentum=0.2,ConvergenceSteps=100,"
                            "BatchSize=5000,TestRepetitions=1,DropConfig=0.0+0.1+0.1+0.1"
                            "Regularization=L2,WeightDecay=1e-4,"
                            "MaxEpochs=4000,"
                            "Optimizer=Adam,"
                            "Multithreading=True")
    DNN_opts = ROOT.TString("!H:V:ErrorStrategy=SUMOFSQUARES:WeightInitialization=XAVIER:VarTransform=G:Architecture=CPU:ValidationSize=20%")
    DNN_opts.Append(":")
    DNN_opts.Append(Layout)
    DNN_opts.Append(":")
    DNN_opts.Append(Strategy)
    DNN_base = factory.BookMethod(dataloader, ROOT.TMVA.Types.kDL, "DNN", DNN_opts)
    
    if args.fast:
        # Run TMVA minimally
        # TrainAllMethods runs the event loop to get the ResultsRegression
        # To get the model only 
        # https://root.cern/doc/master/MethodBase_8cxx_source.html#l00650
        # https://root.cern/doc/master/MethodBase_8cxx_source.html#l00744
        DNN_base.Train()
        dataloader.GetDataSetInfo().GetDataSet().SetCurrentType(ROOT.TMVA.Types.kTraining)
        DNN_base.WriteStateToFile()
    else:
        # Run TMVA
        factory.TrainAllMethods()
        factory.TestAllMethods()
        factory.EvaluateAllMethods()
    
    # close the files
    input_file.Close()
    output_file.Close()


if __name__ == "__main__": 
    start_time = time.time()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush=True)
    print(color.BLUE + "---Start to train the regression!---" + color.END, flush=True)
         
    parser = get_parser()
    args = parser.parse_args()
    
    train_TMVA()
    
    print(color.BLUE + "---All done!---" + color.END)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))), flush=True)