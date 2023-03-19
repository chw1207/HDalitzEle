import numpy as np
from glob import glob
#####################################################################
# All plots, models, config file will be stored here
OutputDirName  = "./results/Output_Merged2GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EE"
Parquet_signal = glob("./data/DataFrames-Merged-2Gsf-FullRun2UL-EE/*.parquet")
Parquet_bkg    = glob("./data/DataFrames-DYJets-FullRun2UL-*-EE/*.parquet") + glob("./data/DataFrames-QCD-FullRun2UL-EE/*.parquet")
Clfname        = "MergedID"
#####################################################################
# prepare sample
testsize = 0.4  # (0.2 means 20%)
CommonCut = " and ".join([
    "(eleCalibEt > 25.)",
    "(elePresel == 1)",       # H->gg preselection
    "(nGsfMatchToReco >= 2)"  # number of Gsf tracks associated to the electron
])
#####################################################################
Classes = ["Merged-2Gsf", "DYJets", "QCD"]
ClassColors = ["#377eb8", "#3A9679", "#E16262"]
MVA = "XGB"
saveROOT = True # save the training set to root file
doShapPlot = True # estimate the prediction contribution

featureplotparam_json = "FeaturePlotParam_MergedID_EE.json"
features = [
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
    
    # "elePFClusEcalIso",
    # "elePFClusHcalIso",
    # "eleIDMVAIso"
]

# in order to safe the memory, just to load the columns you want
# features, selections, rewei, ptvar, etavar
# columns = None to load all of the columns
columns = features+["Class", "eleCalibEt", "elePresel", "nVtx", "nGsfMatchToReco", "eleSCEta", "instwei"]

opt = True # do the hyper parameter searching or not
RandomState = 40
model_params = {
    "objective": "multi:softprob",
    "num_class": len(Classes),
    "eval_metric": "mlogloss",
    "tree_method": "gpu_hist",
    "random_state": RandomState
}

# if opt = False, this should be specified
hyper_params = {
    "colsample_bytree": 0.7225396062573325,
    "grow_policy": "lossguide",
    "learning_rate": 0.07862366155766044,
    "max_delta_step": 4,
    "max_depth": 9,
    "min_child_weight": 655.7809329936954,
    "min_split_loss": 1.7123973653711355,
    "num_class": 3,
    "reg_lambda": 1.4015566968649453,
    "subsample": 0.9187857742123362
}

# number of trees
num_boost_round = 2000
early_stopping_rounds = 50
#####################################################################
# Reweighting scheme
# ! Notice that "intwei" should be included in the input files
# * instwei: (relative xs) ggF as 1, for it being the production with the largest xs and VBF as xs_VBF/xs_ggF
# * if uniformwt = True, then instwei will be assigned as 1. (treat electrons from different processes equally)
# * Reweighing = Balanced(balanced reweight) | one of the class in "Classes"(pt-eta reweight)
uniformwt  = False
Reweighing = "Merged-2Gsf"
ptbins     = [25, 35, 50, 65, 85, 120, 160, 250, 20000]
etabins    = [-2.5, -2.2, -1.8, -1.479, 1.479, 1.8, 2.2, 2.5]
ptwtvar    = "eleCalibEt"
etawtvar   = "eleSCEta"