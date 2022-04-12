import ROOT
import os, sys
import time
from glob import glob
from datetime import datetime
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
ROOT.gROOT.ProcessLine(".L rdfGENQCD.C+")

start_time = time.time()
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush = True)
print(color.BLUE + "---Start to produce the generator mini trees!---" + color.END, flush = True)

# fall 17
ROOT.rdfGENQCD(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_QCD_HT50to100/*.root"),     "./minitree/2017/Minitree_QCD_HT50to100_2017_RECO.root",    2017, "2017", "HT50to100")
ROOT.rdfGENQCD(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_QCD_HT100to200*/*.root"),   "./minitree/2017/Minitree_QCD_HT100to200_2017_RECO.root",   2017, "2017", "HT100to200")
ROOT.rdfGENQCD(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_QCD_HT200to300/*.root"),    "./minitree/2017/Minitree_QCD_HT200to300_2017_RECO.root",   2017, "2017", "HT200to300")
ROOT.rdfGENQCD(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_QCD_HT300to500/*.root"),    "./minitree/2017/Minitree_QCD_HT300to500_2017_RECO.root",   2017, "2017", "HT300to500")
ROOT.rdfGENQCD(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_QCD_HT500to700/*.root"),    "./minitree/2017/Minitree_QCD_HT500to700_2017_RECO.root",   2017, "2017", "HT500to700")
ROOT.rdfGENQCD(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_QCD_HT700to1000/*.root"),   "./minitree/2017/Minitree_QCD_HT700to1000_2017_RECO.root",  2017, "2017", "HT700to1000")
ROOT.rdfGENQCD(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_QCD_HT1000to1500/*.root"),  "./minitree/2017/Minitree_QCD_HT1000to1500_2017_RECO.root", 2017, "2017", "HT1000to1500")
ROOT.rdfGENQCD(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_QCD_HT1500to2000/*.root"),  "./minitree/2017/Minitree_QCD_HT1500to2000_2017_RECO.root", 2017, "2017", "HT1500to2000")
ROOT.rdfGENQCD(find_files("/data6/ggNtuples/V10_02_10_08/job_fall17_QCD_HT2000toInf/*.root"),   "./minitree/2017/Minitree_QCD_HT2000toInf_2017_RECO.root",  2017, "2017", "HT2000toInf")

print(color.BLUE + "---All done!---" + color.END, flush = True)
seconds = time.time() - start_time
print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))