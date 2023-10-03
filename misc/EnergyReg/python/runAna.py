import ROOT
import time
import os, sys
import sampleCfg as cfg
from glob import glob
from colorPrint import *
from datetime import datetime

#! NOTE: Do not execute this script in python folder, otherwise you need to modify the
#!       path to save the miniTree. -> use runAll.sh to run

def mergeFiles(target, source):
    haddfileName = target
    Filelist = sorted(glob(source))

    # make sure the hadd target file is not in the glob file list
    if haddfileName in Filelist:
        Filelist.remove(haddfileName)

    filestr = ""
    for f in Filelist:
        filestr += " {}".format(f)
    os.system("hadd -f {}{}".format(haddfileName, filestr))
    print("", flush=True)


def main():

    # specify the path to store the minitrees
    outDir = "./miniTree/addDiMuCuts"
    os.makedirs("{}/{}".format(outDir, era), exist_ok=True)

    if (run == "ZGToLLG"):
        target = cfg.MCSample["ZGToLLG_{}".format(era)]

        # For double muon trigger
        for file in target["path"]:
            ROOT.xAna_v2(
                file,
                "{}/{}/miniTree_ZGToLLG_DoubleMuTrig_{}.root".format(outDir, era, era),
                "DoubleMuTrig", era,
                target["xs"][0], target["lumi"][0]
            )

        # For single muon trigger
        for file in target["path"]:
            ROOT.xAna_v2(
                file,
                "{}/{}/miniTree_ZGToLLG_SingleMuTrig_{}.root".format(outDir, era, era),
                "SingleMuTrig", era,
                target["xs"][0], target["lumi"][0]
            )


    if (run == "DYJets"):
        target = cfg.MCSample["DYJets_{}".format(era)]

        # For double muon trigger
        for file in target["path"]:
            ROOT.xAna_v2(
                file,
                "{}/{}/miniTree_DYJets_DoubleMuTrig_{}.root".format(outDir, era, era),
                "DoubleMuTrig", era,
                target["xs"][0], target["lumi"][0]
            )

        # For single muon trigger
        for file in target["path"]:
            ROOT.xAna_v2(
                file,
                "{}/{}/miniTree_DYJets_SingleMuTrig_{}.root".format(outDir, era, era),
                "SingleMuTrig", era,
                target["xs"][0], target["lumi"][0]
            )


    if (run == "TTJets"):
        target = cfg.MCSample["TTJets_{}".format(era)]

        # For double muon trigger
        for file in target["path"]:
            ROOT.xAna_v2(
                file,
                "{}/{}/miniTree_TTJets_DoubleMuTrig_{}.root".format(outDir, era, era),
                "DoubleMuTrig", era,
                target["xs"][0], target["lumi"][0]
            )

        # For single muon trigger
        for file in target["path"]:
            ROOT.xAna_v2(
                file,
                "{}/{}/miniTree_TTJets_SingleMuTrig_{}.root".format(outDir, era, era),
                "SingleMuTrig", era,
                target["xs"][0], target["lumi"][0]
            )


    if (run == "Data1Mu"):
        for i, file in enumerate(cfg.DataSample["Data1Mu_{}".format(era)]["path"]):
            dataset = cfg.DataSample["Data1Mu_{}".format(era)]["run"]
            ROOT.xAna_v2(
                file,
                "{}/{}/miniTree_Data_SingleMuTrig_{}.root".format(outDir, era, dataset[i]),
                "SingleMuTrig", era,
                1, 1
            )
        mergeFiles(
            "{}/{}/miniTree_Data_SingleMuTrig_{}.root".format(outDir, era, era),
            "{}/{}/miniTree_Data_SingleMuTrig_*.root".format(outDir, era)
        )


    if (run == "Data2Mu"):
        for i, file in enumerate(cfg.DataSample["Data2Mu_{}".format(era)]["path"]):
            dataset = cfg.DataSample["Data2Mu_{}".format(era)]["run"]
            ROOT.xAna_v2(
                file,
                "{}/{}/miniTree_Data_DoubleMuTrig_{}.root".format(outDir, era, dataset[i]),
                "DoubleMuTrig", era,
                1, 1
            )
        mergeFiles(
            "{}/{}/miniTree_Data_DoubleMuTrig_{}.root".format(outDir, era, era),
            "{}/{}/miniTree_Data_DoubleMuTrig_*.root".format(outDir, era)
        )


if __name__ == "__main__":
    # setup the build directory
    ROOT.gROOT.SetBatch()
    ROOT.gSystem.SetBuildDir("build", ROOT.kTRUE)

    # link the compiled shared libarary
    loc = ROOT.gSystem.Getenv("HDalitzEle_LOC") # declare in env.sh
    if (loc == ""):
        raise ValueError("Cannot fine HDalitzEle_LOC. Please source env.sh in advance!")

    ROOT.gSystem.AddIncludePath("-I{}/include".format(loc))
    ROOT.gSystem.Load("{}/lib/libHDalitzEle.so".format(loc))

    # compile the macro and execute it
    opt = ROOT.gSystem.GetFlagsOpt() # "-O3"
    opt += " -std=c++17 -Wall -Wextra"
    ROOT.gSystem.SetFlagsOpt(opt)
    print("Info in <TUnixSystem::ACLiC>: compilation option is {}".format(opt), flush=True)
    ROOT.gROOT.ProcessLine(".L xAna_v2.C+") # compile the macro with ACLiC

    start_time = time.time()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED + "Execution date and time = {}".format(dt_string) + color.END, flush=True)
    print(color.BLUE + "---Start to produce the mini trees!---" + color.END, flush=True)
    print("", flush=True)

    run = sys.argv[1]
    runs = ["ZGToLLG", "Data2Mu", "Data1Mu", "DYJets", "TTJets"]
    if run not in runs:
        print(runs)
        sys.exit(-1)

    era = sys.argv[2]
    eras = ["UL2016preVFP", "UL2016postVFP", "UL2017", "UL2018"]
    if era not in eras:
        print(eras)
        sys.exit(-1)

    main()

    print(color.BLUE + "---All done!---" + color.END, flush = True)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))