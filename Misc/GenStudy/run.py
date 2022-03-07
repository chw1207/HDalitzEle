#!/usr/bin/env python
import argparse
from ROOT import *
import os
import sys

gROOT.SetBatch()

parser = argparse.ArgumentParser(
    description='python script to run xAna')
parser.add_argument('-r', '--runfile', dest='runfile', type=str, default='all',
                    help='which samples to run')
opt = parser.parse_args()

runHDalitz = False
runDYJetsToLL = False
runGJets = False
runAll = False
if opt.runfile == 'all':
    runAll = True
elif opt.runfile == 'dalitz':
    runHDalitz = True
elif opt.runfile == 'DYjets':
    runDYJetsToLL = True
elif opt.runfile == 'gjets':
    runGJets = True
else:
    print('[INFO] runfile argument [{}] is not currently supported. Please check!\n'.format(opt.runfile))


dirtomake = ['./minitree/', './minitree/2016/', './minitree/2017/', './minitree/2018/']
for d in dirtomake:
    os.makedirs(d, exist_ok=True)

gSystem.AddIncludePath('-Iexternal')
gSystem.SetBuildDir('tmpdir', kTRUE)
gROOT.ProcessLine('.L RecoLevel.C++')

if runHDalitz or runAll:
    RecoLevel("/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_m125/*.root", "test_EEG.root", "2017", "HDalitz", "ggF", 125, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_m125/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_ggF_eeg_m125_2017_RECO.root", "2017", "HDalitz", "ggF", 125, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_VBF_m125/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_VBF_eeg_m125_2017_RECO.root", "2017", "HDalitz", "VBF", 125, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_ZH_m125/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_ZH_eeg_m125_2017_RECO.root", "2017", "HDalitz", "ZH", 125, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_WH_m125/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_WH_eeg_m125_2017_RECO.root", "2017", "HDalitz", "WH", 125, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_m120/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_ggF_eeg_m120_2017_RECO.root", "2017", "HDalitz", "ggF", 120, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_VBF_m120/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_VBF_eeg_m120_2017_RECO.root", "2017", "HDalitz", "VBF", 120, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_ZH_m120/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_ZH_eeg_m120_2017_RECO.root", "2017", "HDalitz", "ZH", 120, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_WH_m120/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_WH_eeg_m120_2017_RECO.root", "2017", "HDalitz", "WH", 120, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_m120/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_ggF_eeg_m130_2017_RECO.root", "2017", "HDalitz", "ggF", 130, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_VBF_m120/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_VBF_eeg_m130_2017_RECO.root", "2017", "HDalitz", "VBF", 130, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_ZH_m120/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_ZH_eeg_m130_2017_RECO.root", "2017", "HDalitz", "ZH", 130, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_WH_m120/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_WH_eeg_m130_2017_RECO.root", "2017", "HDalitz", "WH", 130, 11)

if runDYJetsToLL or runAll:
    RecoLevel("/data6/ggNtuples/V10_02_10_07/job_fall17_DYJetsToLL_m50_aMCatNLO*/ggtree_mc_*.root", "./minitree/2017/Minitree_DYJetsToLL_2017_RECO.root", "2017", "DYJetsToLL", "", 125, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_summer16_DYJetsToLL_m50_aMCatNLO_ext2/ggtree_mc_*.root", "./minitree/2016/Minitree_DYJetsToLL_2016_RECO.root", "2016", "DYJetsToLL", "", 125, 11)
    # RecoLevel("/data6/ggNtuples/V10_02_10_07/job_autumn18_DYJetsToLL_m50_aMCatNLO*/ggtree_mc_*.root", "./minitree/2018/Minitree_DYJetsToLL_2018_RECO.root", "2018", "DYJetsToLL", "", 125, 11)
