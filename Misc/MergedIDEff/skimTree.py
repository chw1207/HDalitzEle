import ROOT # version: 6.24/02
import os, sys
import time
import pickle
import uproot # version: 4.1.9
import numpy as np # version: 1.21.4
import pandas as pd # version: 1.3.2
from glob import glob
from datetime import datetime
from argparse import ArgumentParser
from plugins.colorPrint import *
import plugins.SampleConfig as cfg # get the sample information and Merged ID information
import xgboost as xgb # version: 1.4.2


def find_files(indir):
    Infile_list = sorted(glob(indir))
    f = ROOT.std.vector("string")()
    for i in Infile_list:
        f.push_back(i)
    return f


# For M1: classes [0, 1, 2] -> [1, 2, 3]
# For M2: classes [0, 1, 2] -> [0, 2, 3]
# Therefore the class of a electron should be in [0, 1, 2, 3] which corresponds to ["Merged-2Gsf", "Merged-1Gsf", "DYJets", "QCD"].
def convert_class(arr, Type):
    if Type not in ["Merged-1Gsf", "Merged-2Gsf"]:
        print("This type of model is not available!(Merged-1Gsf or Merged-2Gsf)")
        sys.exit(1)
    new_arr = [i+1 for i in arr] if Type == "Merged-1Gsf" else [0 if i == 0 else i+1 for i in arr]

    return np.asarray(new_arr)


# It contains the function to generate gsf tracks related branches and the function to calculate the phoConvWei
ROOT.gInterpreter.ProcessLine("""
# include "../../pluginsV2/skim_utilities.h"
""")


# Reference:
# [1]: RDataFrame class: https://root.cern/doc/master/classROOT_1_1RDataFrame.html
# [2]: Introduction: https://root-forum.cern.ch/t/rdataframe-a-modern-tool-to-manipulate-and-analyze-root-datasets/29384
#![NOTE]: When ImplicitMT is used to do the Snapshot, the entries are shuffled.
#!https://root-forum.cern.ch/t/snapshot-shuffles-data-invisibly-with-implicitmt/43473
class PreProcess():
    #____________________________________________________________________________________
    def __init__(self, inputlist, outname, xs, luminosity, isMC, ncpu):
        self.outname = outname
        self.inputlist = inputlist
        
        print(color.GREEN + "Add weights to the branches(takes time to calculate...): " + color.END, flush = True)
        MCwei = 1.
        
        if (ncpu == -1):
            ROOT.EnableImplicitMT()
        else:
            ROOT.EnableImplicitMT(ncpu)
            
        df1 = ROOT.RDataFrame("ggNtuplizer/EventTree", self.inputlist)
        
        if (isMC == True):
            pos = df1.Filter("genWeight > 0", "positive event to calculate mcwei").Count()
            neg = df1.Filter("genWeight <= 0", "negative event to calculate mcwei").Count()
            totalev = pos.GetValue() - neg.GetValue()
            MCwei = ((xs * luminosity)/totalev)
            print("# of total events with genweight = {}".format(totalev), flush = True)
            print("mcwei = {}, XS = {}, Lumi = {}".format(MCwei, xs, luminosity), flush = True)
            
        self.df = df1.Define("mcwei", str(MCwei))

    #____________________________________________________________________________________
    def runSkim(self):
        print(color.GREEN + "Start to skim the ggNtuple: " + color.END, flush = True)
        
        df2 = (
            self.df
            .Filter("((HLTEleMuX >> 14) & 1) == 1 || ((HLTEleMuX >> 15) & 1) == 1 || ((HLTEleMuX >> 41) & 1) == 1 || ((HLTEleMuX >> 42) & 1) == 1", "Pass HLT")
            .Filter("isPVGood == 1", "Good Vtx")
            .Filter("(nEle > 0) && (nPho > 0) && (nGSFTrk > 0) && (nMu > 1)", "Nonzero")
            
            .Define("TrkIdxVec",         "GetTrkIdx(nEle, eleCalibEn, eleEta, elePhi, nGSFTrk, gsfPt, gsfEta, gsfPhi)")
            .Define("eleTrkIdx",         "TrkIdxVec[0]")
            .Define("eleSubTrkIdx",      "TrkIdxVec[1]")
            .Define("nGsfMatchToReco",   "TrkIdxVec[2]")

            .Define("eleTrkPt",          "MatchIdex(eleTrkIdx, gsfPt)")
            .Define("eleTrkEta",         "MatchIdex(eleTrkIdx, gsfEta)")
            .Define("eleTrkPhi",         "MatchIdex(eleTrkIdx, gsfPhi)")
            .Define("eleTrkCharge",      "MatchIdex(eleTrkIdx, gsfCharge)")
            .Define("eleTrkLayers",      "MatchIdex(eleTrkIdx, gsfLayers)")
            .Define("eleTrkMissHits",    "MatchIdex(eleTrkIdx, gsfMissHits)")
            .Define("eleTrkD0",          "MatchIdex(eleTrkIdx, gsfD0)")
            .Define("eleTrkDz",          "MatchIdex(eleTrkIdx, gsfDz)")

            .Define("eleSubTrkPt",       "MatchIdex(eleSubTrkIdx, gsfPt)")
            .Define("eleSubTrkEta",      "MatchIdex(eleSubTrkIdx, gsfEta)")
            .Define("eleSubTrkPhi",      "MatchIdex(eleSubTrkIdx, gsfPhi)")
            .Define("eleSubTrkCharge",   "MatchIdex(eleSubTrkIdx, gsfCharge)")
            .Define("eleSubTrkLayers",   "MatchIdex(eleSubTrkIdx, gsfLayers)")
            .Define("eleSubTrkMissHits", "MatchIdex(eleSubTrkIdx, gsfMissHits)")
            .Define("eleSubTrkD0",       "MatchIdex(eleSubTrkIdx, gsfD0)")
            .Define("eleSubTrkDz",       "MatchIdex(eleSubTrkIdx, gsfDz)")

            .Define("eleTrk1",           "P4Vector(eleTrkPt, eleTrkEta, eleTrkPhi, 0.000511)")
            .Define("eleTrk2",           "P4Vector(eleSubTrkPt, eleSubTrkEta, eleSubTrkPhi, 0.000511)")
            .Define("gsfPtRatio",        "GetTrkPtRatio(nEle, nGsfMatchToReco, eleTrk1, eleTrk2)")
            .Define("gsfDeltaR",         "GetTrkdR(nEle, nGsfMatchToReco, eleTrk1, eleTrk2)")
            .Define("gsfRelPtRatio",     "GetTrkRelPtRatio(nEle, eleCalibEn, eleEta, nGsfMatchToReco, eleTrk1, eleTrk2)")
        )
        df2.Report().Print()

        branches = df2.GetColumnNames()
        branches_remain = ROOT.std.vector("string")()
        for i in branches:
            # useless branches
            if ((str(i)[:3] == "HLT") or (str(i)[:2] == "pf") or (str(i)[:2] == "bc")):
                continue
            # special branch type (not allowed by Snapshot)
            if (str(i) == "TrkIdxVec") or (str(i) == "eleTrk1") or (str(i) == "eleTrk2"):
                continue 
            branches_remain.push_back(i)
        
        print(color.GREEN + "Save skimmed tree in(takes time to execute the event loop...):  " + color.END, flush = True)
        print(self.outname, flush = True)
        df2.Snapshot("ggNtuplizer/EventTree", self.outname, branches_remain) # save the skimmed tree
        ROOT.DisableImplicitMT() #! MT should be closed here -> to do the prediction later
    
        #____________________________________________________________________________________
    # Add the prediction results to the existing skimmed tree
    # non-flat tree(event level) -> flat object level dataframe -> do the prediction -> convert back to non-flat tree -> add the prediction to the skimmed tree  
    def addPred(self, features, models): 
        # features is a dict containing {"M1": feature list for M1 ID, "M2": feature list for M2 ID} 
        # models is a dict containing {"M1EB": xgb model M1EB, "M2EB": xgb model M2EB, "M1EE": xgb model M1EE, "M2EB": xgb model M2EE} 
        print(color.GREEN + "Predict the classes of electrons by xgboost(takes time)..." + color.END, flush = True)
        print("Large dataframe(memory > 500 MB) will be split to small chunks to process", flush = True)

        # open the skimmed tree
        fout = ROOT.TFile(self.outname, "UPDATE") 
        tout = fout.Get("ggNtuplizer/EventTree")
        eleClass = ROOT.std.vector("float")()
        TB1 = tout.Branch("eleClass", eleClass)
        
        branches = list(set(features["M1"] + features["M2"] + ["nGsfMatchToReco"]))
        with uproot.open("{}:ggNtuplizer/EventTree".format(self.outname)) as tree:
            # split the dataframe based on the memory
            for df_flat, report in tree.iterate(branches, step_size = "500 MB", library = "pd", report = True): 
                print(report, flush = True)
                
                df_flat_0gsf = df_flat.query("nGsfMatchToReco == 0")
                df_flat_0gsf.insert(loc = 0, column = "eleClass", value = -1)

                # EB 1gsf prediction 
                df_flat_EB_1gsf = df_flat.query("(abs(eleSCEta) < 1.479) and (nGsfMatchToReco == 1)")
                x_EB_1gsf = xgb.DMatrix(df_flat_EB_1gsf.loc[:,features["M1"]].values)
                df_flat_EB_1gsf.insert(loc = 0, column = "eleClass", value = convert_class(models["M1EB"].predict(x_EB_1gsf).argmax(axis = 1), "Merged-1Gsf"))

                # EB 2gsf prediction 
                df_flat_EB_2gsf = df_flat.query("(abs(eleSCEta) < 1.479) and (nGsfMatchToReco >= 2)")
                x_EB_2gsf = xgb.DMatrix(df_flat_EB_2gsf.loc[:,features["M2"]].values)
                df_flat_EB_2gsf.insert(loc = 0, column = "eleClass", value = convert_class(models["M2EB"].predict(x_EB_2gsf).argmax(axis = 1), "Merged-2Gsf"))

                # EE 1gsf prediction 
                df_flat_EE_1gsf = df_flat.query("(abs(eleSCEta) >= 1.479) and (nGsfMatchToReco == 1)")
                x_EE_1gsf = xgb.DMatrix(df_flat_EE_1gsf.loc[:,features["M1"]].values)
                df_flat_EE_1gsf.insert(loc = 0, column = "eleClass", value = convert_class(models["M1EE"].predict(x_EE_1gsf).argmax(axis = 1), "Merged-1Gsf"))

                # EE 2gsf prediction 
                df_flat_EE_2gsf = df_flat.query("(abs(eleSCEta) >= 1.479) and (nGsfMatchToReco >= 2)")
                x_EE_2gsf = xgb.DMatrix(df_flat_EE_2gsf.loc[:,features["M2"]].values)
                df_flat_EE_2gsf.insert(loc = 0, column = "eleClass", value = convert_class(models["M2EE"].predict(x_EE_2gsf).argmax(axis = 1), "Merged-2Gsf"))

                df_new_EBEE = pd.concat([df_flat_0gsf, df_flat_EB_1gsf, df_flat_EB_2gsf, df_flat_EE_1gsf, df_flat_EE_2gsf], sort = False).sort_index()

                # fill the branch ("entry" is the index name which uproot creates)
                df_new = df_new_EBEE.groupby("entry").agg({"eleClass": lambda x: x.to_list()})
                arr_eleClass = df_new["eleClass"].to_numpy() # loop numpy array is much faster than loop the pandas series

                for i in range(len(arr_eleClass)): # loop entry
                    eleClass.clear()

                    for j in range(len(arr_eleClass[i])): # loop subentry
                        eleClass.push_back(arr_eleClass[i][j])

                    TB1.Fill()
                    
        fout.Write()
        fout.Close()



now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush = True)

ncpus = -1
features = cfg.features
models = {
    "M1EB": pickle.load(open(cfg.models["M1EB"], "rb")),
    "M2EB": pickle.load(open(cfg.models["M2EB"], "rb")),
    "M1EE": pickle.load(open(cfg.models["M1EE"], "rb")),
    "M2EE": pickle.load(open(cfg.models["M2EE"], "rb"))
}

pre = PreProcess("/data4/cmkuo/ggtree_mc_fall17_Zg.root", "testSkimZg.root", 106.17*1000., 41.525, True, ncpus)
pre.runSkim()
pre.addPred(features, models)

# isMC = True
# Era = "Zg_2017"
# path, outpath = cfg.MCSample[Era]["path"], cfg.MCSample[Era]["outpath"]
# lumi = cfg.MCSample[Era]["lumi"][0]
# xs = cfg.MCSample[Era]["xs"]
# production = cfg.MCSample[Era]["production"]

# for i in range(len(path)):
#     print(color.GREEN + "Processing MC production {}...".format(production[i]) + color.END, flush = True)
#     InFile_vector = find_files("{}/*.root".format(path[i]))
#     print("Find_files(): {} files are found in {}".format(InFile_vector.size(), path[i]), flush = True)

#     os.makedirs(outpath[i], exist_ok = True)
#     OutFile = "{}/skim.root".format(outpath[i])

#     pre = PreProcess(InFile_vector, OutFile, 2017, xs[i], lumi, isMC, ncpus)
#     pre.runSkim()
#     pre.addPred(features, models)
#     print("", flush = True)


# isMC = True
# Era = "Zg_2016"
# path, outpath = cfg.MCSample[Era]["path"], cfg.MCSample[Era]["outpath"]
# lumi = cfg.MCSample[Era]["lumi"][0]
# xs = cfg.MCSample[Era]["xs"]
# production = cfg.MCSample[Era]["production"]

# for i in range(len(path)):
#     print(color.GREEN + "Processing MC production {}...".format(production[i]) + color.END, flush = True)
#     InFile_vector = find_files("{}/*.root".format(path[i]))
#     print("Find_files(): {} files are found in {}".format(InFile_vector.size(), path[i]), flush = True)

#     os.makedirs(outpath[i], exist_ok = True)
#     OutFile = "{}/skim.root".format(outpath[i])

#     pre = PreProcess(InFile_vector, OutFile, 2016, xs[i], lumi, isMC, ncpus)
#     pre.runSkim()
#     pre.addPred(features, models)
#     print("", flush = True)

# isMC = True
# Era = "Zg_2018"
# path, outpath = cfg.MCSample[Era]["path"], cfg.MCSample[Era]["outpath"]
# lumi = cfg.MCSample[Era]["lumi"][0]
# xs = cfg.MCSample[Era]["xs"]
# production = cfg.MCSample[Era]["production"]

# for i in range(len(path)):
#     print(color.GREEN + "Processing MC production {}...".format(production[i]) + color.END, flush = True)
#     InFile_vector = find_files("{}/*.root".format(path[i]))
#     print("Find_files(): {} files are found in {}".format(InFile_vector.size(), path[i]), flush = True)

#     os.makedirs(outpath[i], exist_ok = True)
#     OutFile = "{}/skim.root".format(outpath[i])

#     pre = PreProcess(InFile_vector, OutFile, 2018, xs[i], lumi, isMC, ncpus)
#     pre.runSkim()
#     pre.addPred(features, models)
#     print("", flush = True)