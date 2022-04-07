import ROOT
import os, sys
import numpy as np
from array import array
from argparse import ArgumentParser
from plugins.CMS_lumi import CMS_lumi


def get_parser():
    parser = ArgumentParser(description = "python script to make 1D efficiency plots")
    parser.add_argument(
        "-e", "--era", 
        help = "era to run [ 2016_preVFP | 2016_postVFP | 2017 | 2018 | combined ], (default = 2017)", 
        default = "2017",
        type = str
    )
    return parser

def Luminosity(era):
    lumi = 1.
    if (era == "2018"):
        lumi = 59.8
    elif (era == "2017"):
        lumi = 41.5
    elif (era == "2016_preVFP"):
        lumi = 19.3
    elif (era == "2016_postVFP"):
        lumi = 16.6
    else:
        lumi = 137.2
        
    return lumi


def main():
    # MC efficiency calculation
    eff_MC = ROOT.TEfficiency()
    for ev in tree_MC:
        if (ev.convVtxRadius_lep1 > 20 and ev.isHggPho_lep1 != 1):
            continue
        bpassed = ev.eleClass_lep1 == 0
        eff_MC.FillWeighted(bpassed, ev.wei2, ev.phoSCEta_lep1, ev.phoCalibEt_lep1)
    
    eff_MC.SetStatisticOption(ROOT.TEfficiency.kBUniform)
    eff_MC.SetConfidenceLevel(0.683)
    eff_MC.SetPosteriorMode(1)
        
    # Data efficiency calculation    
    eff_Da = ROOT.TEfficiency()
    for ev in tree_Da:
        if (ev.convVtxRadius_lep1 > 20 and ev.isHggPho_lep1 != 1):
            continue
        bpassed = ev.eleClass_lep1 == 0
        eff_Da.FillWeighted(bpassed, ev.wei2, ev.phoSCEta_lep1, ev.phoCalibEt_lep1)
    
    eff_Da.SetStatisticOption(ROOT.TEfficiency.kBUniform)
    eff_Da.SetConfidenceLevel(0.683)
    eff_Da.SetPosteriorMode(1)
        




if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    
    fDa = {
        "2016_preVFP": "./miniTree/2016_preVFP/miniTree_Data_2016_preVFP.root",
        "2016_postVFP": "./miniTree/2016_postVFP/miniTree_Data_2016_postVFP.root",
        "2017": "./miniTree/2017/miniTree_Data_2017.root",
        "2018": "./miniTree/2018/miniTree_Data_2018.root"
    }
    
    fMC = {
        "2016_preVFP": "./miniTree/2016_preVFP/miniTree_ZGToLLG_2016_preVFP.root",
        "2016_postVFP": "./miniTree/2016_postVFP/miniTree_ZGToLLG_2016_postVFP.root",
        "2017": "./miniTree/2017/miniTree_ZGToLLG_2017.root",
        "2018": "./miniTree/2018/miniTree_ZGToLLG_2018.root"
    }

    tree_MC, tree_Da = ROOT.TChain("outTree"), ROOT.TChain("outTree")
    if (args.era != "combined"):
        tree_MC.Add(fMC[args.era])
        tree_Da.Add(fDa[args.era])
    else:
        for f in list(fMC.values()):
            tree_MC.Add(f)
        for f in list(fDa.values()):
            tree_Da.Add(f)
    
    main()