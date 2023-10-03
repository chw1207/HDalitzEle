import sys, os
import ROOT
import json
import numpy as np
from colorPrint import *
from argparse import ArgumentParser


# https://root-forum.cern.ch/t/adding-data-from-an-external-containe
def get_parser():
    parser = ArgumentParser(
        description="combine the events selected by single and double muon HLT."
    )
    parser.add_argument(
        "-r", "--run",
        help="sample to run [ Data | ZGToLLG | DYJets | TTJets ]",
        type=str
    )
    parser.add_argument(
        "-e", "--era",
        help="era to run [ UL2016preVFP | UL2016postVFP | UL2017 | UL2018], (default = UL2017)",
        default="UL2017",
        type=str
    )
    return parser


def main():
    readin = [
        "../miniTree/addDiMuCuts/{}/miniTree_{}_SingleMuTrig_{}.root".format(args.era, args.run, args.era),
        "../miniTree/addDiMuCuts/{}/miniTree_{}_DoubleMuTrig_{}.root".format(args.era, args.run, args.era)
    ]
    print(color.GREEN + "[INFO] Load files for {}, {}: ".format(args.run, args.era) + color.END)
    print(json.dumps(readin, indent=4))

    mrgFile = "../miniTree/addDiMuCuts/{}/miniTree_{}_{}.root".format(args.era, args.run, args.era)
    print("[INFO] Merged file: ")
    print("   >>> {}".format(mrgFile))

    df1 = ROOT.RDataFrame("miniTree", readin[0])
    arr1 = df1.AsNumpy(columns=["event"])["event"]
    
    #! NOTE: Name of tree should be the same between two files
    df2 = ROOT.RDataFrame("miniTree", readin[1])
    arr2 = df2.AsNumpy(columns=["event"])["event"]
    
    # find out the duplicated events
    arr, count = np.unique(np.concatenate([arr1, arr2]), return_counts=True)
    dupEv = arr[count > 1]
    before, after = len(arr2), len(arr2)+len(arr1)-len(dupEv),
    print("[INFO] Number of events for miniTree:")
    print("   >>> with double muon HLT: {}".format(before))
    print("   >>> after adding single muon HLT: {}".format(after))
    print("   >>> {}%% improvements".format(round((after-before)*100/before, 2)))

    @ROOT.Numba.Declare(["long"], "int")
    def FindDup(ev):
        isDupEV = 1 if ev in dupEv else 0
        return isDupEV
    
    un_file = "tmp_unique.root"
    print("[INFO] Create tmp file to save events from single muon file. (no duplicated events)")
    print("   >>> {}".format(un_file))
    df_uni = df1.Define("isDupEvent", "Numba::FindDup(event)").Filter("isDupEvent != 1")
    df_uni.Snapshot("miniTree", un_file)
    
    # merge files
    # ignore the cols: isDupEvent, compresion level f2 is used
    # https://root-forum.cern.ch/t/a-question-about-hadd/5471
    print("")
    os.system("hadd -f2 {} {} {}".format(mrgFile, readin[1], un_file))
    print("")

    print("[INFO] remove {}\n".format(un_file))
    os.system("rm {}".format(un_file))


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    main()