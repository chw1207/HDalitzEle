import ROOT
import time
import pandas as pd
import numpy as np
import os,sys
from tabulate import tabulate
from pluginsV2.sigmaEff import sigmaEff
from pluginsV2.colorPrint import *


def main(df, percentage, topology):
    effsig_list = []
    for i in range(len(cats)):
        topo = ""
        if topology == "Merged-2Gsf":
            topo = "category == {}".format(i+1)
        elif topology == "Merged-1Gsf":
            topo = "category == {}".format(i+7)
        else:
            topo = "category == 13"

        if (topology == "Resolved") and (i != 0):
            continue


        effsig = {}
        df_tmp = df["sig"].Filter(topo)
        arr = df_tmp.AsNumpy(columns = ["CMS_higgs_mass"])["CMS_higgs_mass"]

        xmin_1, xmax_1, effsig["sigma"] = sigmaEff(arr, threshold = 0.683)
        xmin_2, xmax_2, effsig["region"] = sigmaEff(arr, threshold = percentage)
        
        rangePerc = "CMS_higgs_mass > {} && CMS_higgs_mass < {}".format(xmin_2, xmax_2)
        effsig["sigPerc"] = df_tmp.Filter(rangePerc).Stats("CMS_higgs_mass", "weight").GetW()
        effsig["signal"] = df_tmp.Stats("CMS_higgs_mass", "weight").GetW()

        # non-resonant background is estimated from data sideband region scaled to the same mass window as signal(2 sigma)
        sidebandPerc = "(CMS_higgs_mass > 105 && CMS_higgs_mass < {}) || (CMS_higgs_mass < 170 && CMS_higgs_mass > {})".format(xmin_2, xmax_2)
        effsig["bkgPerc"] = df["bkg"].Filter(topo).Filter(sidebandPerc).Stats("CMS_higgs_mass").GetN() * (xmax_2 - xmin_2)/(170 - 105)
        effsig["data"] = df["bkg"].Filter(topo).Stats("CMS_higgs_mass").GetN()
        
        effsig["fPerc"] = effsig["sigPerc"]*100/effsig["bkgPerc"]
        effsig["ZPerc"] = 2*(((effsig["sigPerc"] + effsig["bkgPerc"])*np.log(1 + (effsig["sigPerc"]/effsig["bkgPerc"]))) - effsig["sigPerc"])
        
        effsig_list.append(effsig)

    df = pd.DataFrame(data = effsig_list, index = cats)
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

    outDir = "./SignificanceNoGJets"
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

