import ROOT
import time
import os, sys
from colorPrint import *
from datetime import datetime
from argparse import ArgumentParser


def get_parser():
    parser = ArgumentParser(description = "python script to run the xAnaGenQCD.C")
    parser.add_argument(
        "-e", "--era",
        help = "era to run [ UL2016preVFP | UL2016postVFP | UL2017 | UL2018 ]",
        type = str
    )
    parser.add_argument(
        "-p", "--prefix",
        help = "directory to find the ntuples",
        default = "/data5/ggNtuples/V10_06_30_01/",
        type = str
    )
    return parser


# find the way that Ming names the ggNtuples
# e.g. /data5/ggNtuples/V10_06_30_01/job_UL18_GJets_HT40to100/*.root
def get_name(mode):
    # remove the some chars in era (ex. UL2016preVFP -> UL16, UL2017 -> UL17)
    ext = args.era.replace("20", "").replace("preVFP", "").replace("postVFP", "")

    end = mode
    if args.era == "UL2016preVFP":
        end = end + "_preVFP"
    if args.era == "UL2016postVFP":
        end = end + "_postVFP"

    ntuple =  args.prefix + "job_{}_GJets_{}/*.root".format(ext, end)
    return ntuple


def main():
    # creating directory to save the mini trees
    # for i in ["./miniTree/UL2016preVFP", "./miniTree/UL2016postVFP", "./miniTree/UL2017", "./miniTree/UL2018"]:
    #     os.makedirs(i, exist_ok=True)

    # decay_modes = [
    #     "HT40to100",
    #     "HT100to200",
    #     "HT200to400",
    #     "HT400to600",
    #     "HT600toInf"
    # ]
    # for mode in decay_modes:
    #     infile = get_name(mode)
    #     skimfile = "./miniTree/{}/miniTree_GJets_{}_{}.root".format(args.era, mode, args.era)
    #     # start executing the xAna
    #     try:
    #         xAnaGenGJets(infile, skimfile, args.era, mode)
    #     except KeyboardInterrupt:
    #         print("KeyboardInterrupt exception is caught", flush=True)

    #     print("", flush=True)

    xAnaGenGJets("ggtree_mc.root", "ggtree_mc_skim.root", "UL2017", "HT600toInf")
    
    
if __name__ == "__main__" :

    parser = get_parser()
    args = parser.parse_args()

    eras = ["UL2016preVFP", "UL2016postVFP", "UL2017", "UL2018"]
    if args.era not in eras:
        print("[Error] Unknow era: {}, available era: {}".format(args.era, eras))
        sys.exit(-1)

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
    ROOT.gROOT.ProcessLine(".L xAnaGenGJets.C+O") # compile the macro with ACLiC
    from ROOT import xAnaGenGJets

    start_time = time.time()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush=True)
    print(color.BLUE + "---Start to produce the generator mini trees!---" + color.END, flush=True)

    main()

    print(color.BLUE + "---All done!---" + color.END, flush=True)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))