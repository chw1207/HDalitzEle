import os, sys
import glob
import time
import pandas as pd
import numpy as np
import ROOT
import pickle
from Tools.PlotTools import *
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")


class Preparation():
    #____________________________________________________________________________________
    def __init__(self, inFile, sample):
        print(">>> Loading the minitrees: ", flush = True)
        print(">>> {}".format(inFile), flush = True)
        self.rdf = ROOT.RDataFrame("outTree", inFile)
        self.sample = sample

        branches_arr = np.array(self.rdf.GetColumnNames(), dtype = "object")
        self.branches_arr1 = [str(i) for i in branches_arr if ("_lep1" or "_lep2" in str(i)) and (self.rdf.GetColumnType(i) != "TLorentzVector") and ("mc" not in str(i))]
        self.branches_arr2 = ["rho", "nVtx", "nGoodVtx", "isPVGood", "instwei", "procXS", "category"]
        self.branches_arr3 = self.branches_arr2[:-1]

        self.branches = self.branches_arr1 + self.branches_arr2 if (self.sample != "QCD") else self.branches_arr1 + self.branches_arr3

    #____________________________________________________________________________________
    def rdf2df(self):
        print(">>> Processing {} pandas dataframe".format(self.sample), flush = True)
        if (self.sample == "HDalitz"):
            df = pd.DataFrame(self.rdf.AsNumpy(columns = self.branches))

            catDict = {"Merged-2Gsf": 2, "Merged-1Gsf": 3}
            dfsig_list = []
            for cat in catDict.keys():
                df = df.query("category == {}".format(catDict[cat]))
                df_lep1 = df.filter(regex = "lep1") # variables of lep1
                df_lep1.columns = df_lep1.columns.str.replace("_lep1", "") # rename the columns
                df_lep1[self.branches_arr3] = df[self.branches_arr3].to_numpy()

                df_lep1 = df_lep1.dropna()
                df_lep1["Class"] = cat # for multi-classification

                dfsig_list.append(df_lep1)

            df_new = pd.concat(dfsig_list, ignore_index = True, sort = False)

        elif (self.sample == "DYJets"):
            # small training sets
            # https://root-forum.cern.ch/t/use-rdataframe-to-produce-small-training-set/37920/3
            df = pd.DataFrame(self.rdf.AsNumpy(columns = self.branches))

            df = df.query("category == 1")
            lep1 = [i for i in self.branches_arr1 if "_lep1" in i] # braches for lep1
            lep2 = [i.replace("_lep1", "_lep2") for i in lep1] # branches for lep2

            # value is randomly picked from column lep1 or lep2
            # https://stackoverflow.com/a/59752385
            np.random.seed(1234)
            df_new = df[self.branches_arr3]
            idx = np.random.randint(0, 2, df_new.shape[0])
            for i in range(len(lep1)):
                v = df[[lep1[i], lep2[i]]].values
                df_new.insert(loc = 0, column = lep1[i].replace("_lep1", ""), value = np.take_along_axis(v, idx[:, None], 1))

            df_new = df_new.dropna()
            df_new["Class"] = "DYJets" # for multi-classification

        elif (self.sample == "QCD"):
            df = pd.DataFrame(self.rdf.AsNumpy(columns = self.branches))

            df_new = df.filter(regex = "lep1") # variables of lep1
            df_new.columns = df_new.columns.str.replace("_lep1", "") # rename the columns
            df_new[self.branches_arr3] = df[self.branches_arr3].to_numpy()

            df_new = df_new.dropna()
            df_new["Class"] = "QCD" # for multi-classification

        else:
            print("[ERROR] Sample name should be in [HDalitz, DYJets, QCD]", flush = True)
            sys.exit(-1)

        return df_new


def df2pickle(df, outName):
    datadir = "./data/Merged-ID"
    os.makedirs(datadir, exist_ok = True)

    region = ["EB", "EE"]
    sel = ["(abs(eleSCEta) <= 1.479)", "(abs(eleSCEta) > 1.479 and abs(eleSCEta) <= 2.5)"]

    for i, reg in enumerate(region):
        df_tmp = df.query(sel[i])
        print("Save pkl file in {}/{}_{}.pkl".format(datadir, outName, reg))
        with open("{}/{}_{}.pkl".format(datadir, outName, reg), "wb") as f:
            pickle.dump(df_tmp, f)


def main():
    SampleDict = {
        "HDalitz": "../GENstudy/minitree/2017/Minitree_HDalitz_*.root",
        "DYJets": "../GENstudy/minitree/2017/Minitree_DYJets_*.root",
        "QCD": "../GENstudy/minitree/2017/Minitree_QCD_*.root",
    }

    bkgData = []
    for key, value in SampleDict.items():
        if (key == "HDalitz"): # signal
            pre = Preparation(value, key)
            df_sig = pre.rdf2df()
            df2pickle(df_sig, "Dataframe_MergedID_{}_Fall17".format(key))

        else: # background
            pre = Preparation(value, key)
            df_bkgtmp = pre.rdf2df()
            bkgData.append(df_bkgtmp)
        print("", flush = True)

    df_bkg = pd.concat(bkgData, ignore_index = True, sort = False)
    df2pickle(df_bkg, "Dataframe_MergedID_DYJets_QCD_Fall17".format(key))


if __name__ == "__main__":
    start_time = time.time()

    print(color.BLUE + "---Start to prepare the dataframe for training!---" + color.END, flush = True)

    ROOT.EnableImplicitMT()
    main()

    print(color.BLUE + "---All done!---" + color.END, flush = True)

    seconds = time.time() - start_time
    print("Time Taken:", time.strftime("%H:%M:%S", time.gmtime(seconds)), flush = True)
