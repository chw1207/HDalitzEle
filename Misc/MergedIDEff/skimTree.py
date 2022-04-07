import ROOT  # version: 6.24/02
import os
import sys
import time
import pickle
import uproot  # version: 4.1.9
import numpy as np  # version: 1.21.4
import pandas as pd  # version: 1.3.2
from glob import glob
from datetime import datetime
from argparse import ArgumentParser
from plugins.colorPrint import *
# get the sample information and Merged ID information
import plugins.SampleConfig as cfg
import xgboost as xgb  # version: 1.4.2


# This is the python script to skim the ggNtuple and add the XGBoost prediction results
def get_parser():
    parser = ArgumentParser(description="python script to skim the ggNtuple")
    parser.add_argument(
        "-r", "--run",
        help="sample to run [ Data | ZGToLLG | TTJets ]",
        type=str
    )
    parser.add_argument(
        "-e", "--era",
        help="era to run [ 2016_preVFP | 2016_postVFP | 2017 | 2018 ], (default = 2017)",
        default="2017",
        type=str
    )
    return parser


def find_files(indir):
    Infile_list = sorted(glob(indir))
    f = ROOT.std.vector("string")()
    for i in Infile_list:
        f.push_back(i)
    return f


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
    # ____________________________________________________________________________________
    def __init__(self, inputlist, outname, xs, luminosity, isMC, ncpu):
        self.outname = outname
        self.inputlist = inputlist

        fileName = self.outname.split("/")[-1]
        directory = self.outname.replace(fileName, "")
        if not os.path.exists(directory):
            print("Create directory: {}".format(directory), flush=True)
            os.makedirs(directory)

        MCwei = 1.
        if (ncpu == -1):
            ROOT.EnableImplicitMT()
        else:
            ROOT.EnableImplicitMT(ncpu)

        df1 = ROOT.RDataFrame("ggNtuplizer/EventTree", self.inputlist)

        if (isMC == True):
            print(color.GREEN + "Add weights to the branches(takes time to calculate...): " +
                  color.END, flush=True)
            pos = df1.Filter("genWeight > 0", "positive event").Count()
            neg = df1.Filter("genWeight <= 0", "negative event").Count()
            totalev = pos.GetValue() - neg.GetValue()
            MCwei = ((xs * luminosity)/totalev)
            print("# of total events with genweight = {}".format(totalev), flush=True)
            print("mcwei = {}, XS = {}, Lumi = {}".format(MCwei, xs, luminosity), flush=True)

        self.df = df1.Define("mcwei", str(MCwei))

    # ____________________________________________________________________________________
    def runSkim(self):
        print(color.GREEN + "Start to skim the ggNtuple: " + color.END, flush=True)

        df2 = (
            self.df
            # minimum selections on the events
            .Filter("((HLTEleMuX >> 14) & 1) == 1 || ((HLTEleMuX >> 15) & 1) == 1 || ((HLTEleMuX >> 41) & 1) == 1 || ((HLTEleMuX >> 42) & 1) == 1", "Pass HLT")
            .Filter("isPVGood == 1", "Good Vtx")
            .Filter("(nEle > 0) && (nPho > 0) && (nGSFTrk > 0) && (nMu > 1)", "Nonzero")

            # add GSF track information
            .Define("isMainGSF",         "IsMainGSF(eleD0, gsfD0)")
            .Define("ambiguousGSF",      "TrackElectronMatching(eleD0, gsfD0, isMainGSF)")
            .Define("nGsfMatchToReco",   "CalcNGsfMatchToReco(nEle, ambiguousGSF)")
            .Define("eleTrkIdx",         "FindMainGSF(nEle, ambiguousGSF)")
            .Define("eleSubTrkIdx",      "FindSubGSF_PtMax(nEle, ambiguousGSF, gsfPt)")

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
            .Define("gsfPtSum",          "GetTrkPtSum(nEle, nGsfMatchToReco, eleTrk1, eleTrk2)")
            .Define("gsfRelPtRatio",     "GetTrkRelPtRatio(nEle, eleCalibPt, nGsfMatchToReco, eleTrk1, eleTrk2)")
        )
        df2.Report().Print()

        branches = df2.GetColumnNames()
        branches_remain = ROOT.std.vector("string")()
        for i in branches:
            # useless branches
            if ((str(i)[:3] == "HLT") or (str(i)[:2] == "pf") or (str(i)[:2] == "bc")):
                continue
            # special branch type (not allowed by Snapshot)
            if (str(i) == "ambiguousGSF") or (str(i) == "eleTrk1") or (str(i) == "eleTrk2"):
                continue
            branches_remain.push_back(i)

        print(color.GREEN + "Save skimmed tree in(takes time to execute the event loop...):  " + color.END, flush=True)
        print(self.outname, flush=True)
        df2.Snapshot("ggNtuplizer/EventTree", self.outname, branches_remain)  # save the skimmed tree
        ROOT.DisableImplicitMT()  # ! MT should be closed here -> to do the prediction later

    # ____________________________________________________________________________________
    # For M1: classes [0, 1, 2] -> [1, 2, 3]
    # For M2: classes [0, 1, 2] -> [0, 2, 3]
    # Therefore the class of a electron should be in [0, 1, 2, 3] which corresponds to ["Merged-2Gsf", "Merged-1Gsf", "DYJets", "QCD"].
    def convert_class(self, arr, Type):
        if Type not in ["Merged-1Gsf", "Merged-2Gsf"]:
            print("This type of model is not available!(Merged-1Gsf or Merged-2Gsf)")
            sys.exit(1)

        new_arr = [i+1 for i in arr] if Type == "Merged-1Gsf" else [0 if i == 0 else i+1 for i in arr]

        return np.asarray(new_arr)

    # ____________________________________________________________________________________
    # Add the prediction results to the existing skimmed tree
    # non-flat tree(event level) -> flat object level dataframe -> do the prediction -> convert back to non-flat tree -> add the prediction to the skimmed tree
    def addPred(self, features, models):
        # features is a dict containing {"M1": feature list for M1 ID, "M2": feature list for M2 ID}
        # models is a dict containing {"M1EB": xgb model M1EB, "M2EB": xgb model M2EB, "M1EE": xgb model M1EE, "M2EB": xgb model M2EE}
        print(color.GREEN + "Predict the classes of electrons by xgboost(takes time)..." +
              color.END, flush=True)
        print("Large dataframe(memory > 500 MB) will be split to small chunks to process", flush=True)

        # open the skimmed tree
        fout = ROOT.TFile(self.outname, "UPDATE")
        tout = fout.Get("ggNtuplizer/EventTree")
        eleClass = ROOT.std.vector("float")()
        TB1 = tout.Branch("eleClass", eleClass)

        branches = list(
            set(features["M1"] + features["M2"] + ["nGsfMatchToReco"]))
        with uproot.open("{}:ggNtuplizer/EventTree".format(self.outname)) as tree:
            # split the dataframe based on the memory
            for df_flat, report in tree.iterate(branches, step_size="500 MB", library="pd", report=True):
                print(report, flush=True)

                df_flat_0gsf = df_flat.query("nGsfMatchToReco == 0")
                df_flat_0gsf.insert(loc=0, column="eleClass", value=-1)

                # EB 1gsf prediction
                df_flat_EB_1gsf = df_flat.query("(abs(eleSCEta) < 1.479) and (nGsfMatchToReco == 1)")
                x_EB_1gsf = xgb.DMatrix(df_flat_EB_1gsf.loc[:, features["M1"]].values)
                df_flat_EB_1gsf.insert(loc=0, column="eleClass", value=self.convert_class(models["M1EB"].predict(x_EB_1gsf).argmax(axis=1), "Merged-1Gsf"))

                # EB 2gsf prediction
                df_flat_EB_2gsf = df_flat.query("(abs(eleSCEta) < 1.479) and (nGsfMatchToReco >= 2)")
                x_EB_2gsf = xgb.DMatrix(df_flat_EB_2gsf.loc[:, features["M2"]].values)
                df_flat_EB_2gsf.insert(loc=0, column="eleClass", value=self.convert_class(models["M2EB"].predict(x_EB_2gsf).argmax(axis=1), "Merged-2Gsf"))

                # EE 1gsf prediction
                df_flat_EE_1gsf = df_flat.query("(abs(eleSCEta) >= 1.479) and (nGsfMatchToReco == 1)")
                x_EE_1gsf = xgb.DMatrix(df_flat_EE_1gsf.loc[:, features["M1"]].values)
                df_flat_EE_1gsf.insert(loc=0, column="eleClass", value=self.convert_class(models["M1EE"].predict(x_EE_1gsf).argmax(axis=1), "Merged-1Gsf"))

                # EE 2gsf prediction
                df_flat_EE_2gsf = df_flat.query("(abs(eleSCEta) >= 1.479) and (nGsfMatchToReco >= 2)")
                x_EE_2gsf = xgb.DMatrix(df_flat_EE_2gsf.loc[:, features["M2"]].values)
                df_flat_EE_2gsf.insert(loc=0, column="eleClass", value=self.convert_class(models["M2EE"].predict(x_EE_2gsf).argmax(axis=1), "Merged-2Gsf"))

                df_new_EBEE = pd.concat([df_flat_0gsf, df_flat_EB_1gsf, df_flat_EB_2gsf, df_flat_EE_1gsf, df_flat_EE_2gsf], sort=False).sort_index()

                # fill the branch ("entry" is the index name which uproot creates)
                df_new = df_new_EBEE.groupby("entry").agg({"eleClass": lambda x: x.to_list()})
                # loop numpy array is much faster than loop the pandas series
                arr_eleClass = df_new["eleClass"].to_numpy()

                for i in range(len(arr_eleClass)):  # loop entry
                    eleClass.clear()

                    for j in range(len(arr_eleClass[i])):  # loop subentry
                        eleClass.push_back(arr_eleClass[i][j])

                    TB1.Fill()

        fout.Write()
        fout.Close()


def main():
    if (sample == "Data"):
        isMC = False
        Era = "{}_{}".format(sample, era)
        path, outpath = cfg.DataSample[Era]["path"], cfg.DataSample[Era]["outpath"]
        lumi = cfg.DataSample[Era]["lumi"]
        xs = 1.
        run = cfg.DataSample[Era]["run"]

        for i in range(len(path)):
            print(color.GREEN + ">>> Processing data run {}...".format(run[i]) + color.END, flush=True)
            InFile_vector = find_files("{}/*.root".format(path[i]))
            print("Find_files(): {} files are found in {}".format(InFile_vector.size(), path[i]), flush=True)
            OutFile = "{}/skim.root".format(outpath[i])

            pre = PreProcess(InFile_vector, OutFile, xs, lumi[i], isMC, ncpus)
            pre.runSkim()
            pre.addPred(features, models)
            print("", flush=True)

    elif (sample == "ZGToLLG" or sample == "TTJets"):
        isMC = True
        Era = "{}_{}".format(sample, era)
        path, outpath = cfg.MCSample[Era]["path"], cfg.MCSample[Era]["outpath"]
        lumi = cfg.MCSample[Era]["lumi"][0]
        xs = cfg.MCSample[Era]["xs"]
        production = cfg.MCSample[Era]["production"]

        for i in range(len(path)):
            print(color.GREEN + ">>> Processing MC production {}...".format(production[i]) + color.END, flush=True)
            InFile_vector = find_files("{}/*.root".format(path[i]))
            print("Find_files(): {} files are found in {}".format(InFile_vector.size(), path[i]), flush=True)
            OutFile = "{}/skim.root".format(outpath[i])

            pre = PreProcess(InFile_vector, OutFile, xs[i], lumi, isMC, ncpus)
            pre.runSkim()
            pre.addPred(features, models)
            print("", flush=True)


#===============================================#
#            Set up the whole script            #
#  Sample information: pluginsV.SamplConfig.py  #
#===============================================#
if __name__ == "__main__":
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED+"Execution date and time = {}".format(dt_string) +
          color.END, flush=True)

    start_time = time.time()

    # get parser to determine which sample to run
    parser = get_parser()
    args = parser.parse_args()
    sample, era = args.run, args.era
    sample_list = ["Data", "ZGToLLG", "TTJets"]
    era_list = ["2016_preVFP", "2016_postVFP", "2017", "2018"]
    if (sample not in sample_list) or (era not in era_list):
        parser.print_help()
        sys.exit(1)

    print(color.BLUE + "---Start to preprocess the ggNtuple!---" +
          color.END, flush=True)
    print("", flush=True)

    ncpus = -1
    features = cfg.features
    models = {
        "M1EB": pickle.load(open(cfg.models["M1EB"], "rb")),
        "M2EB": pickle.load(open(cfg.models["M2EB"], "rb")),
        "M1EE": pickle.load(open(cfg.models["M1EE"], "rb")),
        "M2EE": pickle.load(open(cfg.models["M2EE"], "rb"))
    }

    main()

    print(color.BLUE + "---All done!---" + color.END, flush=True)
    seconds = time.time() - start_time
    print("Time Taken:", time.strftime(
        "%H:%M:%S", time.gmtime(seconds)), flush=True)
