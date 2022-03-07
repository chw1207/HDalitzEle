#!/usr/bin/env python
import argparse
from ROOT import *
import os
import sys

gROOT.SetBatch()

# parser = argparse.ArgumentParser(
#     description='python script to run xAna')
# parser.add_argument('-r', '--runfile', dest='runfile', type=str, default='all',
#                     help='which samples to run')
# opt = parser.parse_args()

# runHDalitz = False
# runDYJetsToLL = False
# runGJets = False
# runAll = False
# if opt.runfile == 'all':
#     runAll = True
# elif opt.runfile == 'GJets':
#     runGJets = True
# else:
#     print('[INFO] runfile argument [{}] is not currently supported. Please check!\n'.format(opt.runfile))


dirtomake = ["./minitree/", "./minitree/2016/", "./minitree/2017/", "./minitree/2018/", "./status"]
for d in dirtomake:
    os.makedirs(d, exist_ok=True)

gSystem.AddIncludePath("-Iexternal")
gSystem.SetBuildDir("tmpdir", kTRUE)
gROOT.ProcessLine(".L RecoLevel_qcd.C++")


# def runAna(proc):
#     sys.stdout = open("./status/status_qcd_{}.txt".format(proc), "w")
#     RecoLevel_qcd("/data6/ggNtuples/V10_02_10_07/job_fall17_QCD_{}*/ggtree_mc_*.root".format(proc), "./minitree/2017/Minitree_QCD_{}_2017_RECO.root".format(proc), "2017", "QCD", "{}".format(proc))

# procs = [
#     "HT100to200",
#     "HT200to300",
#     "HT300to500",
#     "HT500to700",
#     "HT700to1000",
#     "HT1000to1500",
#     "HT1500to2000",
#     "HT2000toInf"
# ]

# if __name__ == "__main__":

#     tmp = []
#     for p in procs:
#         process_tmp = Process(target = runAna, args = (p,))
#         process_tmp.start()
#         tmp.append(process_tmp)

#     for i in tmp:
#         i.join()


# pool = Pool()
# pool.map_async(runAna, procs)
# pool.close()
# pool.join()
# tmp = []
# for p in procs:
#     process_tmp = Process(target = runAna, args = (p,))
#     process_tmp.start()
#     tmp.append(process_tmp)




RecoLevel_qcd("/data6/ggNtuples/V10_02_10_07/job_fall17_QCD_HT100to200*/ggtree_mc_*.root", "./minitree/2017/Minitree_QCD_HT100to200_2017_RECO.root", "2017", "QCD", "HT100to200")
RecoLevel_qcd("/data6/ggNtuples/V10_02_10_07/job_fall17_QCD_HT200to300/ggtree_mc_*.root", "./minitree/2017/Minitree_QCD_HT200to300_2017_RECO.root", "2017", "QCD", "HT200to300")
RecoLevel_qcd("/data6/ggNtuples/V10_02_10_07/job_fall17_QCD_HT300to500/ggtree_mc_*.root", "./minitree/2017/Minitree_QCD_HT300to500_2017_RECO.root", "2017", "QCD", "HT300to500")
RecoLevel_qcd("/data6/ggNtuples/V10_02_10_07/job_fall17_QCD_HT500to700/ggtree_mc_*.root", "./minitree/2017/Minitree_QCD_HT500to700_2017_RECO.root", "2017", "QCD", "HT500to700")
RecoLevel_qcd("/data6/ggNtuples/V10_02_10_07/job_fall17_QCD_HT700to1000/ggtree_mc_*.root", "./minitree/2017/Minitree_QCD_HT700to1000_2017_RECO.root", "2017", "QCD", "HT700to1000")
RecoLevel_qcd("/data6/ggNtuples/V10_02_10_07/job_fall17_QCD_HT1000to1500/ggtree_mc_*.root", "./minitree/2017/Minitree_QCD_HT1000to1500_2017_RECO.root", "2017", "QCD", "HT1000to1500")
RecoLevel_qcd("/data6/ggNtuples/V10_02_10_07/job_fall17_QCD_HT1500to2000/ggtree_mc_*.root", "./minitree/2017/Minitree_QCD_HT1500to2000_2017_RECO.root", "2017", "QCD", "HT1500to2000")
RecoLevel_qcd("/data6/ggNtuples/V10_02_10_07/job_fall17_QCD_HT2000toInf/ggtree_mc_*.root", "./minitree/2017/Minitree_QCD_HT2000toInf_2017_RECO.root", "2017", "QCD", "HT2000toInf")
