import sys, os
import ROOT
from pprint import pprint
from glob import glob
from CMS_lumi import CMS_lumi

ROOT.gROOT.LoadMacro("../interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")

rdf_test = ROOT.RDataFrame("dataset_reg_EB/TestTree", "reg_results_EB_signal.root")\
               .Define("bias", "BDTG-diGenEle.Pt___D_eleCalibPt_Lead")
# rdf_test = ROOT.RDataFrame("dataset_reg_EB/TrainTree","reg_results_EB_signal.root")