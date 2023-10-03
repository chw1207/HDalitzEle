import ROOT
import os
import time
from CMS_lumi import CMS_lumi
from sigmaEff import sigmaEff
from glob import glob
import numpy as np

ROOT.EnableImplicitMT(10)
infile_sigs = glob("../miniTree/*/miniTree_HDalitz_*_eeg_*.root")
rdf_sigs = ROOT.RDataFrame("miniTree", infile_sigs).Define("weight", "mcwei * genwei")\
               .Filter("elePresel_Lead == 1 && category == 2 && nGsfMatchToReco_Lead >= 2")
    
sigs_before = rdf_sigs.Sum("weight").GetValue()
sigs_after  = rdf_sigs.Filter("eleConvVeto_Lead == 1 && eleTrkMissHits_Lead == 0 && eleSubTrkMissHits_Lead == 0").Sum("weight").GetValue()
print("{0}% of signal events are suppresed by conversion veto and missing hits".format((sigs_before-sigs_after)*100/sigs_before))
    
    
infile_gjet = glob("../miniTree/*/miniTree_GJets_*.root")
rdf_gjet = ROOT.RDataFrame("miniTree", infile_gjet).Define("weight", "mcwei * genwei")\
               .Filter("elePresel_Lead == 1 && nGsfMatchToReco_Lead >= 2")
               
gjet_before = rdf_gjet.Sum("weight").GetValue()
gjet_after  = rdf_gjet.Filter("eleConvVeto_Lead == 1 && eleTrkMissHits_Lead == 0 && eleSubTrkMissHits_Lead == 0").Sum("weight").GetValue()
print("{0}% of photon conversion events are suppresed by conversion veto and missing hits".format((gjet_before-gjet_after)*100/gjet_before))