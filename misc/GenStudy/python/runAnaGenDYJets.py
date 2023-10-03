import ROOT
import os, sys
import time
from datetime import datetime
from glob import glob
from colorPrint import *
from argparse import ArgumentParser


# ntuples for full run2
ntuples = {
    "UL2016preVFP" : "/data5/ggNtuples/V10_06_30_01/job_UL16_DYJetsToLL_m50_aMCatNLO_preVFP/*.root",
    "UL2016postVFP": "/data5/ggNtuples/V10_06_30_01/job_UL16_DYJetsToLL_m50_aMCatNLO_postVFP/*.root",
    "UL2017"       : "/data5/ggNtuples/V10_06_30_01/job_UL17_DYJetsToLL_m50_aMCatNLO/*.root",
    "UL2018"       : "/data5/ggNtuples/V10_06_30_01/job_UL18_DYJetsToLL_m50_aMCatNLO/*.root"
}
skimpath = {
    "UL2016preVFP" : "./miniTree/UL2016preVFP/miniTree_DYJets_UL2016preVFP.root",
    "UL2016postVFP": "./miniTree/UL2016postVFP/miniTree_DYJets_UL2016postVFP.root",
    "UL2017"       : "./miniTree/UL2017/miniTree_DYJets_UL2017.root",
    "UL2018"       : "./miniTree/UL2018/miniTree_DYJets_UL2018.root",
}


def main():
    start_time = time.time()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush=True)
    print(color.BLUE + "---Start to produce the generator mini trees!---" + color.END, flush=True)

    xAnaGenDYJets(ntuples[era], skimpath[era], era)

    print(color.BLUE + "---All done!---" + color.END, flush = True)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))


if __name__ == "__main__":
    era = sys.argv[1]
    eras = ["UL2016preVFP", "UL2016postVFP", "UL2017", "UL2018"]
    if era not in eras:
        print("[Error] Unknow era: {}, available era: {}".format(era, eras))
        sys.exit(-1)

    dirtomake = ["./miniTree/UL2016preVFP", "./miniTree/UL2016postVFP", "./miniTree/UL2017", "./miniTree/UL2018"]
    for d in dirtomake:
        os.makedirs(d, exist_ok=True)

    # setup the build directory
    ROOT.gROOT.SetBatch()
    ROOT.gSystem.SetBuildDir("build", ROOT.kTRUE)

    # link the compiled shared libarary
    loc = ROOT.gSystem.Getenv("HDalitzEle_LOC") # declare in env.sh
    if (loc == ""):
        print("[ERROR] Please source env.sh in advance!")
        sys.exit(-1)

    ROOT.gSystem.AddIncludePath("-I{}/include".format(loc))
    ROOT.gSystem.Load("{}/lib/libHDalitzEle.so".format(loc))

    # compile the macro and execute it
    ROOT.gROOT.ProcessLine(".L xAnaGenDYJets.C+O") # compile the macro with ACLiC
    from ROOT import xAnaGenDYJets

    main()