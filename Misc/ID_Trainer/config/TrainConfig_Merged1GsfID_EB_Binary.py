import numpy as np

#####################################################################
# All plots, models, config file will be stored here
OutputDirName = "./Results/Output_Merged1GsfID_dask_Fall17_EB_CheckALL_Binary"
Pickle_signal = "./data/Merged-ID/Dataframe_MergedID_HDalitz_Fall17_EB.pkl"
Pickle_bkg = "./data/Merged-ID/Dataframe_MergedID_DYJets_QCD_Fall17_EB.pkl"
Clfname = "MergedID"
#####################################################################
# prepare sample
testsize = 0.2  #(0.2 means 20%)
CommonCut = "(elePresel == 1) and (nVtx > 1) and (nGsfMatchToReco == 1) and (Class != 'Merged-2Gsf')"
#####################################################################
Classes = ["Signal", "Background"]
ClassColors = ["#377eb8", "#E16262"]
MVA = "XGB"

featureplotparam_json = "FeaturePlotParam_MergedID_EB.json"
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

RandomState = 42
param = {
    "objective": "binary:logistic",
    "eval_metric": "logloss",
    "tree_method": "gpu_hist",
    "random_state": RandomState,

    "learning_rate": 0.2,
    "max_delta_step": 7,
    "max_depth": 5,
    "min_child_weight": 100,
    "min_split_loss": 1,
    "subsample": 0.5
}

# number of trees
num_boost_round = 500
early_stopping_rounds = 10
#####################################################################
# Reweighting scheme
# ! Notice that "intwei" should be included in the input pkl files
# * instwei: (relative xs) ggF as 1, for it being the production with the largest xs and VBF as xs_VBF/xs_ggF
# * Reweighing = Balanced | one of the class in "Classes"
Reweighing = "Signal"
ptbins = [0., 15., 25., 35., 60., 20000]
etabins = np.arange(-1.6, 1.6, 0.2)
ptwtvar = "eleCalibPt"
etawtvar = "eleSCEta"
