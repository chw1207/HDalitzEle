import ROOT
import time
import os, sys
from colorPrint import *
from datetime import datetime
from argparse import ArgumentParser


def get_parser():
    parser = ArgumentParser(description = "python script to run the xAnaGen.C")
    parser.add_argument(
        "-e", "--era",
        help = "era to run [ UL2016preVFP | UL2016postVFP | UL2017 | UL2018 ]",
        type = str
    )
    parser.add_argument(
        "-p", "--prefix",
        help = "directory to find the ntuples",
        default = "/data4/chenghan/muon/skimTree/",
        type = str
    )
    return parser


# find the way that Ming names the ggNtuples
# e.g. /data5/ggNtuples/V10_06_30_01/job_UL16_Dalitz_eeg_m125_preVFP/*.root
def get_name(mode, mpoint):
    # remove the some chars in era (ex. UL2016preVFP -> UL16, UL2017 -> UL17)
    ext = args.era.replace("20", "").replace("preVFP", "").replace("postVFP", "")

    end = "{}_m{}".format(mode, mpoint)
    if mode == "ggF":
        end = "m{}".format(mpoint)
    if args.era == "UL2016preVFP":
        end = end + "_preVFP"
    if args.era == "UL2016postVFP":
        end = end + "_postVFP"

    ntuple =  args.prefix + "job_{}_Dalitz_mmg_{}.root".format(ext, end)
    return ntuple


def main():
    # creating directory to save the mini trees
    for i in ["./miniTree/UL2016preVFP", "./miniTree/UL2016postVFP", "./miniTree/UL2017", "./miniTree/UL2018"]:
        os.makedirs(i, exist_ok=True)

    decay_modes = [
        "ggF",
        "VBF",
        "WH",
        "ZH",
        "ttH",
        "bbH"
    ]
    mass_points = [
        120,
        125,
        130
    ]
    for mode in decay_modes:
        for mpoint in mass_points:
            infile = get_name(mode, mpoint)
            skimfile = "./miniTree/{}/miniTree_HDalitz_{}_mmg_{}_{}.root".format(args.era, mode, mpoint, args.era)

            # start executing the xAna
            try:
                ROOT.xAnaGenMu(infile, skimfile, args.era, mode, 125)
            except KeyboardInterrupt:
                print("KeyboardInterrupt exception is caught", flush=True)
                sys.exit(-1)

            print("", flush=True)
            
            
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

    # # compile the macro and execute it
    ROOT.gROOT.ProcessLine(".L xAnaGenMu.C++") # compile the macro with ACLiC

    start_time = time.time()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush=True)
    print(color.BLUE + "---Start to produce the generator mini trees!---" + color.END, flush=True)

    main()

    print(color.BLUE + "---All done!---" + color.END, flush=True)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))