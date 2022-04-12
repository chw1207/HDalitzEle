import ROOT
import os, sys
import time
from datetime import datetime
from glob import glob
from plugins.colorPrint import *
from argparse import ArgumentParser


def find_files(indir):
    Infile_list = sorted(glob(indir))
    f = ROOT.std.vector("string")()
    for i in Infile_list:
        f.push_back(i)
    return f


dirtomake = ["./minitree/", "./minitree/2016/", "./minitree/2017/", "./minitree/2018/"]
for d in dirtomake:
    os.makedirs(d, exist_ok = True)

ROOT.gROOT.SetBatch()
ROOT.gSystem.AddIncludePath("-Iexternal")
ROOT.gSystem.SetBuildDir("tmpdir", ROOT.kTRUE)
ROOT.gROOT.ProcessLine(".L rdfGENDYJets.C++")

start_time = time.time()
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush = True)
print(color.BLUE + "---Start to produce the generator mini trees!---" + color.END, flush = True)

ROOT.rdfGENDYJets(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_newPU_DYJetsToLL_m50_aMCatNLO*/*.root"),     "./minitree/2017/Minitree_DYJets_2017_RECO.root",   2017,    "2017",   "DYJets")

print(color.BLUE + "---All done!---" + color.END, flush = True)
seconds = time.time() - start_time
print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))