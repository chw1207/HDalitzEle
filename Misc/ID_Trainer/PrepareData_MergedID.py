import sys
from datetime import datetime
import os
import uproot as uproot4
import glob
import pandas as pd
import numpy as np
import ROOT
import pickle
import gc
from Tools.PlotTools import *
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")

# convert object dtype to float
# def convert_dtype(df):
#     cols = df.select_dtypes(include = [object]).columns
#     df[cols] = df[cols].astype(np.float32)
#     return df

# create merged categories df
def prepare_Merged_df(df, cat):
    if ((cat != "Merged-2Gsf") and (cat != "Merged-1Gsf")):
        print("No such category! [Merged-2Gsf | Merged-1Gsf]")
        sys.exit(1)
    
    sel = "category == 2" if cat == "Merged-2Gsf" else "category == 3"
    df_merged = df.query(sel)

    df_merged_lep1 = df_merged.filter(regex = "lep1") # variables of lep1
    df_merged_lep1.columns = df_merged_lep1.columns.str.replace("_lep1", "") # rename the columns
    df_merged_lep1[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]] = df_merged[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]].to_numpy()

    df_merged_new = df_merged_lep1.sample(frac = 1).reset_index(drop = True) # shuffle DataFrame rows
    df_merged_new = df_merged_new.dropna()
    df_merged_new["Class"] = cat

    return df_merged_new

def prepare_GJetsQCD_df(df, cat):
    if ((cat != "GJets") and (cat != "QCD")):
        print("No such category! [GJets | QCD]")
        sys.exit(1)

    df_new = df.sample(frac = 1, random_state = 42).reset_index(drop = True) # shuffle DataFrame rows
    df_new["Class"] = cat
    df_new = df_new.dropna()
    return df_new

def prepare_DYJets_df(df):
    df_resolved = df.query("category == 1")
    df_resolved_lep1 = df_resolved.filter(regex = "lep1") # variables of lep1
    df_resolved_lep2 = df_resolved.filter(regex = "lep2") # variables of lep1

    df_resolved_lep1.columns = df_resolved_lep1.columns.str.replace("_lep1", "") # rename the columns
    df_resolved_lep2.columns = df_resolved_lep2.columns.str.replace("_lep2", "") # rename the columns

    df_resolved_lep1[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]] = df_resolved[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]].to_numpy()
    df_resolved_lep2[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]] = df_resolved[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]].to_numpy()

    # combine two leptons 
    df_new = pd.concat([df_resolved_lep1, df_resolved_lep2], ignore_index = True, sort = False)
    df_new = df_new.sample(frac = 1, random_state = 42).reset_index(drop = True) # shuffle DataFrame rows
    df_new["Class"] = "DYJets"
    df_new = df_new.dropna()
    return df_new

def df2pickle(df, datasetname):
    datadir = "./data/Merged-ID"
    os.makedirs(datadir, exist_ok=True)

    region = ["EB","EE"]
    sel = ["(abs(eleSCEta) <= 1.479)", "(abs(eleSCEta) > 1.479 and abs(eleSCEta) <= 2.5)"]

    for i, reg in enumerate(region):
        df_tmp = df.query(sel[i])

        file = open("{}/Dataframe_MergedID_{}_{}_2017.pkl".format(datadir, datasetname, reg),"wb")
        pickle.dump(df_tmp, file)
        file.close()

    file = open("{}/Dataframe_MergedID_{}_2017.pkl".format(datadir, datasetname),"wb")
    pickle.dump(df, file)
    print("Save pkl file in {}/Dataframe_MergedID_{}_2017.pkl".format(datadir, datasetname))
    file.close()

def main():
    print("RDataFrame is loading HDalitz minitrees to pandas dataframe...")
    df_tmp_HDalitz = ROOT.RDataFrame("outTree", "{}/Minitree_HDalitz_*.root".format(treedir))
    df_tmp_HDalitz_array = df_tmp_HDalitz.AsNumpy(columns = feature_HDalitz_DY)
    df_HDalitz = pd.DataFrame(df_tmp_HDalitz_array)
    # df_HDalitz = convert_dtype(df_HDalitz)

    print("RDataFrame is loading GJets minitrees to pandas dataframe...")
    df_tmp_GJets = ROOT.RDataFrame("outTree", "{}/Minitree_gjet_*_MGG_*.root".format(treedir))
    df_tmp_GJets_array = df_tmp_GJets.AsNumpy(columns = feature_GJets)
    df_GJets = pd.DataFrame(df_tmp_GJets_array)
    # df_GJets = convert_dtype(df_GJets)

    print("RDataFrame is loading QCD minitrees to pandas dataframe...")
    df_tmp_QCD = ROOT.RDataFrame("outTree", "{}/Minitree_QCD_*.root".format(treedir))
    df_tmp_QCD_array = df_tmp_QCD.AsNumpy(columns = feature_QCD)
    df_QCD = pd.DataFrame(df_tmp_QCD_array)
    # df_QCD = convert_dtype(df_QCD)

    print("RDataFrame is loading DYJets minitrees to pandas dataframe...")
    df_tmp_DYJets = ROOT.RDataFrame("outTree", "{}/Minitree_DYJetsToLL_*.root".format(treedir))
    df_tmp_DYJets_array = df_tmp_DYJets.AsNumpy(columns = feature_HDalitz_DY)
    df_DYJets = pd.DataFrame(df_tmp_DYJets_array)
    # df_DYJets = convert_dtype(df_DYJets)

    print("Deal with HDalitz pandas dataframe ...")
    df_merged2gsf_HDalitz_new = prepare_Merged_df(df_HDalitz, "Merged-2Gsf")
    df_merged1gsf_HDalitz_new = prepare_Merged_df(df_HDalitz, "Merged-1Gsf")
    # combine merged2gsf and merged1gsf
    df_HDalitz_new = pd.concat([df_merged2gsf_HDalitz_new, df_merged1gsf_HDalitz_new], ignore_index = True, sort = False) 
    df_HDalitz_new = df_HDalitz_new.sample(frac = 1, random_state = 42).reset_index(drop = True) 
    df_HDalitz_new = df_HDalitz_new.dropna()
    print(color.RED + "[INFO] Does HDalitz dataframe have NAN value(s)? {}".format(df_HDalitz_new.isnull().values.any()) + color.END)

    print("Deal with GJets pandas dataframe ...")
    df_GJets_new = prepare_GJetsQCD_df(df_GJets, "GJets")
    print(color.RED + "[INFO] Does GJets dataframe have NAN value(s)? {}".format(df_GJets_new.isnull().values.any()) + color.END)

    print("Deal with QCD pandas dataframe ...")
    df_QCD_new = prepare_GJetsQCD_df(df_QCD, "QCD")
    print(color.RED + "[INFO] Does QCD dataframe have NAN value(s)? {}".format(df_QCD_new.isnull().values.any()) + color.END)

    print("Deal with DYJets pandas dataframe ...")
    df_DYJets_new = prepare_DYJets_df(df_DYJets)
    print (color.RED + "[INFO] Does DYJets dataframe have NAN value(s)? {}".format(df_DYJets_new.isnull().values.any()) + color.END)

    print("Deal with background dataframe(DY+GJets+QCD) dataframe ...")
    df_Background_new = pd.concat([df_GJets_new, df_QCD_new, df_DYJets_new], ignore_index = True, sort = False)
    df_Background_new = df_Background_new.sample(frac = 1, random_state = 42).reset_index(drop = True) # shuffle DataFrame rows
    df_Background_new = df_Background_new.dropna()
    print(color.RED + "[INFO] Does background dataframe(DY+GJets+QCD) dataframe have NAN value(s)? {}".format(df_Background_new.isnull().values.any()) + color.END)

    print("Save signal(HDalitz) pkl file ...")
    df2pickle(df_HDalitz_new, "HDalitz_eeg_allMass")

    print("Save background(DY+GJets+QCD) pkl file ...")
    df2pickle(df_Background_new, "DYJets_GJets_QCD")

if __name__ == "__main__" :
    print (color.BOLD + color.BLUE + "---Data preparation start!---" + color.END)

    treedir = "/home/chenghan/Analysis/Dalitz/electron/Misc/GenStudy/minitree/2017"
    ROOT.EnableImplicitMT(20) # Tell ROOT you want to go parallel

    # HDalitz and DYJets us the same set of feature list in the root file
    feature_HDalitz_DY = [
        "category",
        "rho",
        "nVtx",
        "nGoodVtx",
        "isPVGood",
        "instwei",
        "elePresel_lep1", "elePresel_lep2",
        "eleSCEta_lep1", "eleSCEta_lep2",
        "eleEn_lep1","eleEn_lep2",
        "eleSCEn_lep1","eleSCEn_lep2",
        "eleSCRawEn_lep1","eleSCRawEn_lep2",
        "eleEcalEn_lep1","eleEcalEn_lep2",
        "eleESEnP1_lep1","eleESEnP1_lep2",
        "eleESEnP2_lep1","eleESEnP2_lep2",
        "eleD0_lep1","eleD0_lep2",
        "eleDz_lep1","eleDz_lep2",
        "eleSIP_lep1","eleSIP_lep2",
        "elePtError_lep1","elePtError_lep2",
        "eleEta_lep1","eleEta_lep2",
        "eleR9_lep1","eleR9_lep2",
        "eleSCEtaWidth_lep1","eleSCEtaWidth_lep2",
        "eleSCPhiWidth_lep1","eleSCPhiWidth_lep2",
        "eleHoverE_lep1","eleHoverE_lep2",
        "eleEoverP_lep1","eleEoverP_lep2",
        "eleEoverPout_lep1","eleEoverPout_lep2",
        "eleEoverPInv_lep1","eleEoverPInv_lep2",
        "eleBrem_lep1","eleBrem_lep2",
        "eledEtaAtVtx_lep1","eledEtaAtVtx_lep2",
        "eledPhiAtVtx_lep1","eledPhiAtVtx_lep2",
        "eleSigmaIEtaIEtaFull5x5_lep1","eleSigmaIEtaIEtaFull5x5_lep2",
        "eleSigmaIPhiIPhiFull5x5_lep1","eleSigmaIPhiIPhiFull5x5_lep2",
        "elePFChIso_lep1","elePFChIso_lep2",
        "elePFPhoIso_lep1","elePFPhoIso_lep2",
        "elePFNeuIso_lep1","elePFNeuIso_lep2",
        "elePFPUIso_lep1","elePFPUIso_lep2",
        "eleIDMVAIso_lep1","eleIDMVAIso_lep2",
        "eleIDMVANoIso_lep1","eleIDMVANoIso_lep2",
        "eleR9Full5x5_lep1","eleR9Full5x5_lep2",
        "eleTrkdxy_lep1","eleTrkdxy_lep2",
        "eleMissHits_lep1","eleMissHits_lep2",
        "nGsfMatchToReco_lep1","nGsfMatchToReco_lep2", 
        "circularity_lep1", "circularity_lep2",
        "eleConvVeto_lep1", "eleConvVeto_lep2",
        "eleKFHits_lep1", "eleKFHits_lep2",
        "eleKFChi2_lep1", "eleKFChi2_lep2",
        "eleGSFChi2_lep1", "eleGSFChi2_lep1",
        "eleCalibPt_lep1", "eleCalibPt_lep2",
        "eleEcalDrivenSeed_lep1", "eleEcalDrivenSeed_lep2",
        
        "eleESEnToRawE_lep1", "eleESEnToRawE_lep2",
        "eleESEffSigmaRR_lep1", "eleESEffSigmaRR_lep2",
        
        "gsfPtSum_lep1", "gsfPtRatio_lep1", "gsfDeltaR_lep1", "gsfMissHitsSum_lep1", "gsfPtSum2_lep1", "gsfRelPtRatio_lep1", "gsfPerpEleSum_lep1",

        "gsfPtSum_lep2", "gsfPtRatio_lep2", "gsfDeltaR_lep2", "gsfMissHitsSum_lep2", "gsfPtSum2_lep2", "gsfRelPtRatio_lep2", "gsfPerpEleSum_lep2",

        "gsfDiffPtRatio_lep1", "gsfDiffPtRatio_lep2"
    ]

    feature_GJets = [
        "rho",
        "nVtx",
        "nGoodVtx",
        "isPVGood",
        "instwei",
        "elePresel",
        "eleSCEta",
        "eleEn",
        "eleSCEn",
        "eleSCRawEn",
        "eleEcalEn",
        "eleESEnP1",
        "eleESEnP2",
        "eleD0",
        "eleDz",
        "eleSIP",
        "elePtError",
        "eleEta",
        "eleR9",
        "eleSCEtaWidth",
        "eleSCPhiWidth",
        "eleHoverE",
        "eleEoverP",
        "eleEoverPout",
        "eleEoverPInv",
        "eleBrem",
        "eledEtaAtVtx",
        "eledPhiAtVtx",
        "eleSigmaIEtaIEtaFull5x5",
        "eleSigmaIPhiIPhiFull5x5",
        "elePFChIso",
        "elePFPhoIso",
        "elePFNeuIso",
        "elePFPUIso",
        "eleIDMVAIso",
        "eleIDMVANoIso",
        "eleR9Full5x5",
        "eleTrkdxy",
        "eleMissHits",
        "nGsfMatchToReco",
        "circularity", 
        "eleConvVeto",
        "eleKFHits",
        "eleKFChi2",
        "eleGSFChi2",
        "eleEcalDrivenSeed",
        "eleCalibPt",
        "eleESEnToRawE",
        "eleESEffSigmaRR",
        "gsfPtSum", "gsfPtRatio", "gsfDeltaR", "gsfMissHitsSum", "gsfPtSum2", "gsfRelPtRatio", "gsfPerpEleSum", "gsfDiffPtRatio"
    ] 

    feature_QCD = [
        "rho",
        "nVtx",
        "nGoodVtx",
        "isPVGood",
        "instwei",
        "elePresel",
        "eleSCEta",
        "eleEn",
        "eleSCEn",
        "eleSCRawEn",
        "eleEcalEn",
        "eleESEnP1",
        "eleESEnP2",
        "eleD0",
        "eleDz",
        "eleSIP",
        "elePtError",
        "eleEta",
        "eleR9",
        "eleSCEtaWidth",
        "eleSCPhiWidth",
        "eleHoverE",
        "eleEoverP",
        "eleEoverPout",
        "eleEoverPInv",
        "eleBrem",
        "eledEtaAtVtx",
        "eledPhiAtVtx",
        "eleSigmaIEtaIEtaFull5x5",
        "eleSigmaIPhiIPhiFull5x5",
        "elePFChIso",
        "elePFPhoIso",
        "elePFNeuIso",
        "elePFPUIso",
        "eleIDMVAIso",
        "eleIDMVANoIso",
        "eleR9Full5x5",
        "eleTrkdxy",
        "eleMissHits",
        "nGsfMatchToReco",
        "circularity", 
        "eleConvVeto",
        "eleKFHits",
        "eleKFChi2",
        "eleGSFChi2",
        "eleEcalDrivenSeed",
        "eleCalibPt",
        "eleESEnToRawE",
        "eleESEffSigmaRR",
        "gsfPtSum", "gsfPtRatio", "gsfDeltaR", "gsfMissHitsSum", "gsfPtSum2", "gsfRelPtRatio", "gsfPerpEleSum", "gsfDiffPtRatio"
    ] 

    state = main()

    print (color.BOLD + color.BLUE + "---Data preparation done---" + color.END)
    sys.exit(state)



# print (color.BOLD + color.BLUE + "---Data preparation start!---" + color.END)

# #### ----- feature list: Merged ID ----- ####
# # HDalitz and DYJets us the same set of feature list in the root file
# feature_HDalitz_DY = [
#     "category",
#     "rho",
#     "nVtx",
#     "nGoodVtx",
#     "isPVGood",
#     "instwei",
#     "elePresel_lep1", "elePresel_lep2",
#     "eleSCEta_lep1", "eleSCEta_lep2",
#     "eleEn_lep1","eleEn_lep2",
#     "eleSCEn_lep1","eleSCEn_lep2",
#     "eleSCRawEn_lep1","eleSCRawEn_lep2",
#     "eleEcalEn_lep1","eleEcalEn_lep2",
#     "eleESEnP1_lep1","eleESEnP1_lep2",
#     "eleESEnP2_lep1","eleESEnP2_lep2",
#     "eleD0_lep1","eleD0_lep2",
#     "eleDz_lep1","eleDz_lep2",
#     "eleSIP_lep1","eleSIP_lep2",
#     "elePtError_lep1","elePtError_lep2",
#     "eleEta_lep1","eleEta_lep2",
#     "eleR9_lep1","eleR9_lep2",
#     "eleSCEtaWidth_lep1","eleSCEtaWidth_lep2",
#     "eleSCPhiWidth_lep1","eleSCPhiWidth_lep2",
#     "eleHoverE_lep1","eleHoverE_lep2",
#     "eleEoverP_lep1","eleEoverP_lep2",
#     "eleEoverPout_lep1","eleEoverPout_lep2",
#     "eleEoverPInv_lep1","eleEoverPInv_lep2",
#     "eleBrem_lep1","eleBrem_lep2",
#     "eledEtaAtVtx_lep1","eledEtaAtVtx_lep2",
#     "eledPhiAtVtx_lep1","eledPhiAtVtx_lep2",
#     "eleSigmaIEtaIEtaFull5x5_lep1","eleSigmaIEtaIEtaFull5x5_lep2",
#     "eleSigmaIPhiIPhiFull5x5_lep1","eleSigmaIPhiIPhiFull5x5_lep2",
#     "elePFChIso_lep1","elePFChIso_lep2",
#     "elePFPhoIso_lep1","elePFPhoIso_lep2",
#     "elePFNeuIso_lep1","elePFNeuIso_lep2",
#     "elePFPUIso_lep1","elePFPUIso_lep2",
#     "eleIDMVAIso_lep1","eleIDMVAIso_lep2",
#     "eleIDMVANoIso_lep1","eleIDMVANoIso_lep2",
#     "eleR9Full5x5_lep1","eleR9Full5x5_lep2",
#     "eleTrkdxy_lep1","eleTrkdxy_lep2",
#     "eleMissHits_lep1","eleMissHits_lep2",
#     "nGsfMatchToReco_lep1","nGsfMatchToReco_lep2", 
#     "circularity_lep1", "circularity_lep2",
#     "eleConvVeto_lep1", "eleConvVeto_lep2",
#     "eleKFHits_lep1", "eleKFHits_lep2",
#     "eleKFChi2_lep1", "eleKFChi2_lep2",
#     "eleGSFChi2_lep1", "eleGSFChi2_lep1",
#     "eleCalibPt_lep1", "eleCalibPt_lep2",
#     "eleEcalDrivenSeed_lep1", "eleEcalDrivenSeed_lep2",
    
#     "eleESEnToRawE_lep1", "eleESEnToRawE_lep2",
#     "eleESEffSigmaRR_lep1", "eleESEffSigmaRR_lep2",
    
#     "gsfPtSum_lep1", "gsfPtRatio_lep1", "gsfDeltaR_lep1", "gsfMissHitsSum_lep1", "gsfPtSum2_lep1", "gsfRelPtRatio_lep1", "gsfPerpEleSum_lep1",

#     "gsfPtSum_lep2", "gsfPtRatio_lep2", "gsfDeltaR_lep2", "gsfMissHitsSum_lep2", "gsfPtSum2_lep2", "gsfRelPtRatio_lep2", "gsfPerpEleSum_lep2",

#     "gsfDiffPtRatio_lep1", "gsfDiffPtRatio_lep2"
# ]
 

# feature_GJets = [
#     "rho",
#     "nVtx",
#     "nGoodVtx",
#     "isPVGood",
#     "instwei",
#     "elePresel",
#     "eleSCEta",
#     "eleEn",
#     "eleSCEn",
#     "eleSCRawEn",
#     "eleEcalEn",
#     "eleESEnP1",
#     "eleESEnP2",
#     "eleD0",
#     "eleDz",
#     "eleSIP",
#     "elePtError",
#     "eleEta",
#     "eleR9",
#     "eleSCEtaWidth",
#     "eleSCPhiWidth",
#     "eleHoverE",
#     "eleEoverP",
#     "eleEoverPout",
#     "eleEoverPInv",
#     "eleBrem",
#     "eledEtaAtVtx",
#     "eledPhiAtVtx",
#     "eleSigmaIEtaIEtaFull5x5",
#     "eleSigmaIPhiIPhiFull5x5",
#     "elePFChIso",
#     "elePFPhoIso",
#     "elePFNeuIso",
#     "elePFPUIso",
#     "eleIDMVAIso",
#     "eleIDMVANoIso",
#     "eleR9Full5x5",
#     "eleTrkdxy",
#     "eleMissHits",
#     "nGsfMatchToReco",
#     "circularity", 
#     "eleConvVeto",
#     "eleKFHits",
#     "eleKFChi2",
#     "eleGSFChi2",
#     "eleEcalDrivenSeed",
#     "eleCalibPt",
#     "eleESEnToRawE",
#     "eleESEffSigmaRR",
#     "gsfPtSum", "gsfPtRatio", "gsfDeltaR", "gsfMissHitsSum", "gsfPtSum2", "gsfRelPtRatio", "gsfPerpEleSum", "gsfDiffPtRatio"
# ] 

# feature_QCD = [
#     "rho",
#     "nVtx",
#     "nGoodVtx",
#     "isPVGood",
#     "instwei",
#     "elePresel",
#     "eleSCEta",
#     "eleEn",
#     "eleSCEn",
#     "eleSCRawEn",
#     "eleEcalEn",
#     "eleESEnP1",
#     "eleESEnP2",
#     "eleD0",
#     "eleDz",
#     "eleSIP",
#     "elePtError",
#     "eleEta",
#     "eleR9",
#     "eleSCEtaWidth",
#     "eleSCPhiWidth",
#     "eleHoverE",
#     "eleEoverP",
#     "eleEoverPout",
#     "eleEoverPInv",
#     "eleBrem",
#     "eledEtaAtVtx",
#     "eledPhiAtVtx",
#     "eleSigmaIEtaIEtaFull5x5",
#     "eleSigmaIPhiIPhiFull5x5",
#     "elePFChIso",
#     "elePFPhoIso",
#     "elePFNeuIso",
#     "elePFPUIso",
#     "eleIDMVAIso",
#     "eleIDMVANoIso",
#     "eleR9Full5x5",
#     "eleTrkdxy",
#     "eleMissHits",
#     "nGsfMatchToReco",
#     "circularity", 
#     "eleConvVeto",
#     "eleKFHits",
#     "eleKFChi2",
#     "eleGSFChi2",
#     "eleEcalDrivenSeed",
#     "eleCalibPt",
#     "eleESEnToRawE",
#     "eleESEffSigmaRR",
#     "gsfPtSum", "gsfPtRatio", "gsfDeltaR", "gsfMissHitsSum", "gsfPtSum2", "gsfRelPtRatio", "gsfPerpEleSum", "gsfDiffPtRatio"
# ] 

# #### ----- Load miniTree ----- ####
# treedir = "/home/chenghan/Analysis/Dalitz/electron/Misc/GenStudy/minitree/2017"

# ROOT.EnableImplicitMT(10) # Tell ROOT you want to go parallel

# print ("RDataFrame is loading HDalitz minitrees to pandas dataframe...")
# df_tmp_HDalitz = ROOT.RDataFrame("outTree", "{}/Minitree_HDalitz_*.root".format(treedir))
# df_tmp_HDalitz_array = df_tmp_HDalitz.AsNumpy(columns = feature_HDalitz_DY)
# df_HDalitz = pd.DataFrame(df_tmp_HDalitz_array)
# df_HDalitz = convert_dtype(df_HDalitz)

# print ("RDataFrame is loading GJets minitrees to pandas dataframe...")
# df_tmp_GJets = ROOT.RDataFrame("outTree", "{}/Minitree_gjet_*_MGG_*.root".format(treedir))
# df_tmp_GJets_array = df_tmp_GJets.AsNumpy(columns = feature_GJets)
# df_GJets = pd.DataFrame(df_tmp_GJets_array)

# print ("RDataFrame is loading QCD minitrees to pandas dataframe...")
# df_tmp_QCD = ROOT.RDataFrame("outTree", "{}/Minitree_QCD_*.root".format(treedir))
# df_tmp_QCD_array = df_tmp_QCD.AsNumpy(columns = feature_QCD)
# df_QCD = pd.DataFrame(df_tmp_QCD_array)

# print ("RDataFrame is loading DYJets minitrees to pandas dataframe...")
# df_tmp_DYJets = ROOT.RDataFrame("outTree", "{}/Minitree_DYJetsToLL_*.root".format(treedir))
# df_tmp_DYJets_array = df_tmp_DYJets.AsNumpy(columns = feature_HDalitz_DY)
# df_DYJets = pd.DataFrame(df_tmp_DYJets_array)

# ### ----- HDalitz dataframe ----- ####
# import warnings
# if not sys.warnoptions:
#     warnings.simplefilter("ignore")

# print("Deal with HDalitz pandas dataframe ...")
# # merged2gsf
# df_merged2gsf_HDalitz = df_HDalitz.query("category == 2")
# df_merged2gsf_HDalitz_lep1 = df_merged2gsf_HDalitz.filter(regex = "lep1") # variables of lep1
# df_merged2gsf_HDalitz_lep1.columns = df_merged2gsf_HDalitz_lep1.columns.str.replace("_lep1", "") # rename the columns
# df_merged2gsf_HDalitz_lep1[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]] = df_merged2gsf_HDalitz[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]].to_numpy()
# df_merged2gsf_HDalitz_new = df_merged2gsf_HDalitz_lep1.sample(frac = 1).reset_index(drop = True) # shuffle DataFrame rows
# df_merged2gsf_HDalitz_new = df_merged2gsf_HDalitz_new.dropna()
# df_merged2gsf_HDalitz_new["Class"] = "Merged-2Gsf"

# # merged1gsf
# df_merged1gsf_HDalitz = df_HDalitz.query("category == 3")
# df_merged1gsf_HDalitz_lep1 = df_merged1gsf_HDalitz.filter(regex = "lep1") # variables of lep1
# df_merged1gsf_HDalitz_lep1.columns = df_merged1gsf_HDalitz_lep1.columns.str.replace("_lep1", "") # rename the columns
# df_merged1gsf_HDalitz_lep1[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]] = df_merged1gsf_HDalitz[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]].to_numpy()
# df_merged1gsf_HDalitz_new = df_merged1gsf_HDalitz_lep1.sample(frac = 1, random_state = 42).reset_index(drop = True) # shuffle DataFrame rows
# df_merged1gsf_HDalitz_new = df_merged1gsf_HDalitz_new.dropna()
# df_merged1gsf_HDalitz_new["Class"] = "Merged-1Gsf"

# # combine merged2gsf and merged1gsf
# df_HDalitz_new = pd.concat([df_merged2gsf_HDalitz_new, df_merged1gsf_HDalitz_new], ignore_index = True, sort = False) # concatenate the dataframe 
# df_HDalitz_new = df_HDalitz_new.sample(frac = 1, random_state = 42).reset_index(drop = True) # shuffle DataFrame rows
# df_HDalitz_new = df_HDalitz_new.dropna()
# print (color.RED + "[INFO] Does HDalitz dataframe have NAN value(s)? {}".format(df_HDalitz_new.isnull().values.any()) + color.END)

# #### ----- GJets dataframe ----- ####
# print("Deal with GJets pandas dataframe ...")
# df_GJets_new = df_GJets.sample(frac = 1, random_state = 42).reset_index(drop = True) # shuffle DataFrame rows
# df_GJets_new["Class"] = "GJets"
# df_GJets_new = df_GJets_new.dropna()
# print (color.RED + "[INFO] Does GJets dataframe have NAN value(s)? {}".format(df_GJets_new.isnull().values.any()) + color.END)

# #### ----- QCD dataframe ----- ####
# print("Deal with QCD pandas dataframe ...")
# df_QCD_new = df_QCD.sample(frac = 1, random_state = 42).reset_index(drop = True) # shuffle DataFrame rows
# df_QCD_new["Class"] = "QCD"
# df_QCD_new = df_QCD_new.dropna()
# print (color.RED + "[INFO] Does QCD dataframe have NAN value(s)? {}".format(df_QCD_new.isnull().values.any()) + color.END)

# # #### ----- DYJets dataframe ----- ####
# print("Deal with DYJets pandas dataframe ...")
# df_resolved_DYJets = df_DYJets.query("category == 1")
# df_resolved_DYJets_lep1 = df_resolved_DYJets.filter(regex = "lep1") # variables of lep1
# df_resolved_DYJets_lep2 = df_resolved_DYJets.filter(regex = "lep2") # variables of lep1

# df_resolved_DYJets_lep1.columns = df_resolved_DYJets_lep1.columns.str.replace("_lep1", "") # rename the columns
# df_resolved_DYJets_lep2.columns = df_resolved_DYJets_lep2.columns.str.replace("_lep2", "") # rename the columns

# df_resolved_DYJets_lep1[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]] = df_resolved_DYJets[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]].to_numpy()
# df_resolved_DYJets_lep2[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]] = df_resolved_DYJets[["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei"]].to_numpy()

# # combine two leptons 
# df_DYJets_new = pd.concat([df_resolved_DYJets_lep1, df_resolved_DYJets_lep2], ignore_index = True, sort = False)

# df_DYJets_new = df_DYJets_new.sample(frac = 0.5, random_state = 42).reset_index(drop = True) # shuffle DataFrame rows
# df_DYJets_new["Class"] = "DYJets"
# df_DYJets_new = df_DYJets_new.dropna()
# print (color.RED + "[INFO] Does DYJets dataframe have NAN value(s)? {}".format(df_DYJets_new.isnull().values.any()) + color.END)

# # #### ----- Background dataframe(DY+GJets+QCD) ----- ####
# print("Deal with background dataframe(DY+GJets+QCD) dataframe ...")
# df_Background_new = pd.concat([df_GJets_new, df_QCD_new, df_DYJets_new], ignore_index = True, sort = False)
# df_Background_new = df_Background_new.sample(frac = 1, random_state = 42).reset_index(drop = True) # shuffle DataFrame rows
# df_Background_new = df_Background_new.dropna()
# print (color.RED + "[INFO] Does background dataframe(DY+GJets+QCD) dataframe have NAN value(s)? {}".format(df_Background_new.isnull().values.any()) + color.END)

# #### ----- Save the dataframe to pkl files ----- ####
# def df2pickle(df, datasetname):
#     datadir = "./data/Merged-ID"
#     os.makedirs(datadir, exist_ok=True)

#     region = ["EB","EE"]
#     sel = ["(abs(eleSCEta) <= 1.479)", "(abs(eleSCEta) > 1.479 and abs(eleSCEta) <= 2.5)"]

#     for i, reg in enumerate(region):
#         df_tmp = df.query(sel[i])

#         file = open("{}/Dataframe_MergedID_{}_{}_2017.pkl".format(datadir, datasetname, reg),"wb")
#         pickle.dump(df_tmp, file)
#         file.close()

#     file = open("{}/Dataframe_MergedID_{}_2017.pkl".format(datadir, datasetname),"wb")
#     pickle.dump(df, file)
#     print("Save pkl file in {}/Dataframe_MergedID_{}_2017.pkl".format(datadir, datasetname))
#     file.close()
    
# print("Save signal(HDalitz) pkl file ...")
# df2pickle(df_HDalitz_new, "HDalitz_allprod_eeg_allMass")

# print("Save background(DY+GJets+QCD) pkl file ...")
# df2pickle(df_Background_new, "DYJets_GJets_QCD")

# print (color.BOLD + color.BLUE + "---Data preparation done---" + color.END)
