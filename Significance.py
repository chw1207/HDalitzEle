import ROOT
import time
import pandas as pd
import numpy as np
import os,sys
from tabulate import tabulate
from pluginsV2.sigmaEff import sigmaEff
from pluginsV2.colorPrint import *


def CalcEff(df_sig_cat, df_bkg_cat):
    arr = df_sig_cat.AsNumpy(columns = ["CMS_higgs_mass"])["CMS_higgs_mass"]

    effsig = {}
    xmin_1, xmax_1, effsig["sigma"] = sigmaEff(arr, threshold = 0.683)
    xmin_2, xmax_2, effsig["region"] = sigmaEff(arr, threshold = percentage)

    rangePerc = "CMS_higgs_mass > {} && CMS_higgs_mass < {}".format(xmin_2, xmax_2)
    effsig["sigPerc"] = df_sig_cat.Filter(rangePerc).Stats("CMS_higgs_mass", "weight").GetW()
    effsig["signal"] = df_sig_cat.Stats("CMS_higgs_mass", "weight").GetW()

    # non-resonant background is estimated from data sideband region scaled to the same mass window as signal(2 sigma)
    sidebandPerc = "(CMS_higgs_mass > 105 && CMS_higgs_mass < {}) || (CMS_higgs_mass < 170 && CMS_higgs_mass > {})".format(xmin_2, xmax_2)
    effsig["bkgPerc"] = df_bkg_cat.Filter(sidebandPerc).Stats("CMS_higgs_mass").GetN() * (xmax_2 - xmin_2)/(170 - 105)
    effsig["data"] = df_bkg_cat.Stats("CMS_higgs_mass").GetN()

    effsig["fPerc"] = effsig["sigPerc"]*100/effsig["bkgPerc"]
    effsig["ZPerc"] = 2*(((effsig["sigPerc"] + effsig["bkgPerc"])*np.log(1 + (effsig["sigPerc"]/effsig["bkgPerc"]))) - effsig["sigPerc"])

    return effsig


def main(df, percentage, topology):
    effsig_list = []
    if (topology == "Resolved"):
        colName = ["all"]
        topo = "category == 13"
        df_sig_cat = df["sig"].Filter(topo)
        df_bkg_cat = df["bkg"].Filter(topo)
        effsig_list.append(CalcEff(df_sig_cat, df_bkg_cat))
    else:
        colName = cats
        for i in range(len(cats)):
            if topology == "Merged-2Gsf":
                topo = "category == {}".format(i+1)
            elif topology == "Merged-1Gsf":
                topo = "category == {}".format(i+7)

            df_sig_cat = df["sig"].Filter(topo)
            df_bkg_cat = df["bkg"].Filter(topo)
            effsig_list.append(CalcEff(df_sig_cat, df_bkg_cat))

    df = pd.DataFrame(data = effsig_list, index = colName)
    df = df.round(2)

    print(color.GREEN + color.BOLD + "{} significance table".format(topology) + color.END)
    print(color.GREEN +  ">>> {}".format(percentage) + color.END)
    table = tabulate(df, headers = "keys", tablefmt = "orgtbl")
    print(table)

    print("[INFO] Save the significance table in {}/significance_{}_{}.txt".format(outDir, topology, percentage))
    with open("{}/significance_{}_{}.txt".format(outDir, topology, percentage), "w") as f:
        f.write(table)
    print("")


if __name__ == "__main__" :
    start_time = time.time()

    outDir = "./Significance"
    os.makedirs(outDir, exist_ok = True)

    ROOT.EnableImplicitMT()

    cats = ["HVBF", "LVBF", "BST", "EBHR9", "EBLR9", "EE"]

    fileDict = {
        "sig": "./miniTree/2017/miniTree_HDalitz_*_m125.root",
        "bkg": "./miniTree/2017/miniTree_Data_2017.root"
    }
    df = {}
    for key in fileDict.keys():
        df[key] = ROOT.RDataFrame("outTree", fileDict[key])

    for percentage in [0.955, 0.9]:
        for topology in ["Merged-2Gsf", "Merged-1Gsf", "Resolved"]:
            main(df, percentage, topology)

    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))

