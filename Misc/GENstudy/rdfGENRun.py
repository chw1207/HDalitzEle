import ROOT
import os, sys
import time
from datetime import datetime
from plugins.colorPrint import *
from argparse import ArgumentParser


dirtomake = ["./minitree/", "./minitree/2016/", "./minitree/2017/", "./minitree/2018/"]
for d in dirtomake:
    os.makedirs(d, exist_ok = True)

ROOT.gROOT.SetBatch()
ROOT.gSystem.AddIncludePath("-Iexternal")
ROOT.gSystem.SetBuildDir("tmpdir", ROOT.kTRUE)
ROOT.gROOT.ProcessLine(".L rdfGEN.C+")

start_time = time.time()  
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush = True)
print(color.BLUE + "---Start to produce the generator mini trees!---" + color.END, flush = True)

# signal sample for 2017 @ 125 GeV 
ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_m125/*.root",     outfile = "./minitree/2017/Minitree_HDalitz_ggF_eeg_m125_2017_RECO.root",   year = 2017,    era = "2017",   proc = "HDalitz",   prod = "ggF",   HiggsMass = 125)
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_VBF_m125/*.root", outfile = "./minitree/2017/Minitree_HDalitz_VBF_eeg_m125_2017_RECO.root",   year = 2017,    era = "2017",   proc = "HDalitz",   prod = "VBF",   HiggsMass = 125)
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_ZH_m125/*.root",  outfile = "./minitree/2017/Minitree_HDalitz_ZH_eeg_m125_2017_RECO.root",    year = 2017,    era = "2017",   proc = "HDalitz",   prod = "ZH",    HiggsMass = 125)
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_WH_m125/*.root",  outfile = "./minitree/2017/Minitree_HDalitz_WH_eeg_m125_2017_RECO.root",    year = 2017,    era = "2017",   proc = "HDalitz",   prod = "WH",    HiggsMass = 125)

# # signal sample for 2017 @ 120 GeV 
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_m120/*.root",     outfile = "./minitree/2017/Minitree_HDalitz_ggF_eeg_m120_2017_RECO.root",   year = 2017,    era = "2017",   proc = "HDalitz",   prod = "ggF",   HiggsMass = 120)
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_VBF_m120/*.root", outfile = "./minitree/2017/Minitree_HDalitz_VBF_eeg_m120_2017_RECO.root",   year = 2017,    era = "2017",   proc = "HDalitz",   prod = "VBF",   HiggsMass = 120)
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_ZH_m120/*.root",  outfile = "./minitree/2017/Minitree_HDalitz_ZH_eeg_m120_2017_RECO.root",    year = 2017,    era = "2017",   proc = "HDalitz",   prod = "ZH",    HiggsMass = 120)
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_WH_m120/*.root",  outfile = "./minitree/2017/Minitree_HDalitz_WH_eeg_m120_2017_RECO.root",    year = 2017,    era = "2017",   proc = "HDalitz",   prod = "WH",    HiggsMass = 120)

# # signal sample for 2017 @ 130 GeV 
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_m130/*.root",     outfile = "./minitree/2017/Minitree_HDalitz_ggF_eeg_m130_2017_RECO.root",   year = 2017,    era = "2017",   proc = "HDalitz",   prod = "ggF",   HiggsMass = 130)
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_VBF_m130/*.root", outfile = "./minitree/2017/Minitree_HDalitz_VBF_eeg_m130_2017_RECO.root",   year = 2017,    era = "2017",   proc = "HDalitz",   prod = "VBF",   HiggsMass = 130)
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_ZH_m130/*.root",  outfile = "./minitree/2017/Minitree_HDalitz_ZH_eeg_m130_2017_RECO.root",    year = 2017,    era = "2017",   proc = "HDalitz",   prod = "ZH",    HiggsMass = 130)
# ROOT.rdfGEN(infile = "/data6/ggNtuples/V10_02_10_08/job_fall17_Dalitz_eeg_WH_m130/*.root",  outfile = "./minitree/2017/Minitree_HDalitz_WH_eeg_m130_2017_RECO.root",    year = 2017,    era = "2017",   proc = "HDalitz",   prod = "WH",    HiggsMass = 130)

print(color.BLUE + "---All done!---" + color.END, flush = True)
seconds = time.time() - start_time
print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))