import numpy as np

#####################################################################
# All plots, models, config file will be stored here
OutputDirName = "./Results/Output_Merged1GsfID_hyperTune_Fall17_EE"
Pickle_signal = "./data/Merged-ID/Dataframe_MergedID_HDalitz_Fall17_EE.pkl"
Pickle_bkg = "./data/Merged-ID/Dataframe_MergedID_DYJets_QCD_Fall17_EE.pkl"
Clfname = "MergedID"
#####################################################################
# prepare sample
testsize = 0.2  #(0.2 means 20%)
CommonCut = " and ".join([
    "(eleCalibPt > 5.)",
    "(elePresel == 1)",
    "(nVtx > 1)",
    "(nGsfMatchToReco == 1)",
    "(Class != 'Merged-2Gsf')"
])
#####################################################################
Classes = ["Merged-1Gsf", "DYJets", "QCD"]
ClassColors = ["#377eb8", "#3A9679", "#E16262"]
MVA = "XGB"

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

    # "gsfPtRatio",
    # "gsfDeltaR",
    "gsfRelPtRatio"
]

opt = True # do the hyper parameter searching or not
RandomState = 42
model_params = {
    "objective": "multi:softprob",
    "num_class": len(Classes),
    "eval_metric": "mlogloss",
    "tree_method": "gpu_hist",
    "random_state": RandomState
}

# if opt = False, this should be specified
hyper_params = {
    "learning_rate": 0.02216257774390921,
    "max_delta_step": 10,
    "max_depth": 12,
    "min_child_weight": 134.62634063073708,
    "min_split_loss": 0.23112380376907937,
    "subsample": 0.7933362725563781
}

# number of trees
num_boost_round = 500
early_stopping_rounds = 10

#####################################################################
# Reweighting scheme
# ! Notice that "intwei" should be included in the input pkl files
# * instwei: (relative xs) ggF as 1, for it being the production with the largest xs and VBF as xs_VBF/xs_ggF
# * if uniformwt = True, then instwei will be assigned as 1. (treat electrons from different processes equally)
# * Reweighing = Balanced(balanced reweight) | one of the class in "Classes"(pt-eta reweight)
uniformwt = True
Reweighing = "Merged-1Gsf"
ptbins = [5, 15, 25, 35, 50, 65, 85, 120, 160, 200, 20000]
etabins = [-2.5, -2.2, -1.8, -1.479, 1.479, 1.8, 2.2, 2.5]
ptwtvar = "eleCalibPt"
etawtvar = "eleSCEta"