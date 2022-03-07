import numpy as np

#####################################################################
# All plots, models, config file will be stored here
OutputDirName = "./Results/Output_Merged1GsfID_EB_2017_dask"
Pickle_signal = "./data/Merged-ID/Dataframe_MergedID_HDalitz_eeg_allMass_EB_2017.pkl"
Pickle_bkg = "./data/Merged-ID/Dataframe_MergedID_DYJets_GJets_QCD_EB_2017.pkl"
Clfname = "MergedID"
#####################################################################
# prepare sample
testsize = 0.2  #(0.2 means 20%)
CommonCut = "(elePresel == True) and (nVtx > 1) and (nGsfMatchToReco == 1) and (Class != 'Merged-2Gsf')"
#####################################################################
Classes = ["Merged-1Gsf", "DYJets", "GJets", "QCD"]
ClassColors = ["#377eb8", "#3A9679", "#E16262", "#787A91"]
MVA = "XGB"

featureplotparam_json = "FeaturePlotParam_MergedID_EB.json"
features = [
    "rho",
    "eleSCEta",
    "eleSCRawEn",

    "eleD0",
    "eleDz",
    "eleSIP",
    "elePtError",
    "eleSCPhiWidth",
    "eleEoverP",
    "eleEoverPout",
    "eleEoverPInv",
    "eledEtaAtVtx",
    "eledPhiAtVtx",
    "eleSigmaIEtaIEtaFull5x5",
    "elePFChIso",
    "elePFPhoIso",
    "eleR9Full5x5",
    "eleTrkdxy",
    "eleMissHits",
    "gsfRelPtRatio"
]

RandomState = 42
param = {
    "objective": "multi:softprob",
    "num_class": len(Classes),
    "eval_metric": "mlogloss",
    "tree_method": "gpu_hist",
    "random_state": RandomState,
    "eta": 0.2,
    "max_delta_step": 5,
    "max_depth": 5,
    "min_child_weight": 40,
    "gamma": 1,
    "subsample": 0.5
}

# number of trees
num_boost_round = 500
early_stopping_rounds = 20

#####################################################################
# Reweighting scheme
# ! Notice that "intwei" should be included in the input pkl files
# * instwei: (relative xs) ggF as 1, for it being the production with the largest xs and VBF as xs_VBF/xs_ggF 
Reweighing = "Nothing"
ptbins = [0, 10, 20, 30, 40, 50, 60, 100, 150] 
etabins = np.linspace(-1.5, 1.5, num = 11)
ptwtvar = "eleCalibPt"
etawtvar = "eleEta"