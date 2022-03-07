# Datacard making script: uses output pkl file of makeYields_HDalitz.py script
# ! Systematic variations are yet included (To-do: assume to be 1 for now)
# to-do : split mass maslist argument

import os, sys
import re
import ROOT
import pandas as pd
import pickle
from argparse import ArgumentParser
from glob import glob
from collections import OrderedDict as od
# from systematics import theory_systematics, experimental_systematics, signal_shape_systematics
from tools.writeToDatacard_HDalitz import writePreamble, writeProcesses, writeSystematic, writeParamSystematic, writePdfIndex
from commonTools_HDalitz import *

from systematics_HDalitz import theory_systematics, experimental_systematics
from tools.calcSystematics_HDalitz import addConstantSyst


print(color.GREEN + "Preliminary work......" + color.END)
print(" --> Loading per category dataframes into single dataframe")
pkl_files = glob("./yields/*.pkl")
pkl_files.sort(key = str.lower)
df_data = pd.DataFrame()
for f_pkl_name in pkl_files:
    with open(f_pkl_name, "rb") as f_pkl: 
        df = pickle.load(f_pkl)
        df_data = pd.concat([df_data, df], ignore_index = True, axis = 0, sort = False)

# Theory:
print(" --> Adding theory systematics variations to dataFrame")
for s in theory_systematics:
    if (s["type"] == "constant"):
        df_data = addConstantSyst(df_data, s)

# Experimental:
print(" --> Adding experimental systematics variations to dataFrame")
# Add constant systematics to dataFrame
for s in experimental_systematics:
    if (s["type"] == "constant"): 
        df_data = addConstantSyst(df_data, s)

channel = "heeg"
massList = [120+i for i in range(11)]
cats = catNumMap.keys()

outdir = "./electron" # --> directory for data cards
if not os.path.exists(outdir):
    os.makedirs(outdir)

for cat in cats:
    for mass in massList:
        fdataName = "{outdir}/datacard_{channel}_runII_{category}_{mass_point}.txt".format(outdir = outdir, channel = channel, category = cat, mass_point = mass)
        print("[INFO] Creating the data card {}".format(fdataName))
        fdata = open(fdataName, "w")

        if not writePreamble(fdata): 
            print " --> [ERROR] in writing preamble. Leaving..."
            sys.exit(1)

        if not writeProcesses(fdata, df_data, cat, mass, outdir = outdir):
            print " --> [ERROR] in writing processes. Leaving..."
            sys.exit(1)

        for syst in theory_systematics:
            if not writeSystematic(fdata, df_data, syst, cat, mass, years = [2016, 2017, 2018]):
                print " --> [ERROR] in writing theory systematics. Leaving..."
                sys.exit(1)

        for syst in experimental_systematics:
            if not writeSystematic(fdata, df_data, syst, cat, mass, years = [2016, 2017, 2018]):
                print " --> [ERROR] in writing theory systematics. Leaving..."
                sys.exit(1)

        if not writeParamSystematic(fdata, df_data, cat, mass):
            print " --> [ERROR] in writing pram systematics. Leaving..."
            sys.exit(1)

        if not writePdfIndex(fdata, df_data, cat, ext = "DiPho_13TeV"):
            print " --> [ERROR] in writing pdf indices. Leaving..."
            sys.exit(1)

        fdata.close()


# fdataName = "./test2.txt"
# print("[INFO] Creating the data card {}".format(fdataName))
# fdata = open(fdataName, "w")

# if not writePreamble(fdata): 
#     print " --> [ERROR] in writing preamble. Leaving..."
#     sys.exit(1)

# if not writeProcesses(fdata, df_data, "Merged2Gsf_EE", 125):
#     print " --> [ERROR] in writing processes. Leaving..."
#     sys.exit(1)


# for syst in theory_systematics:
#     if not writeSystematic(fdata, df_data, syst, "Merged2Gsf_EE", 125, years = [2016, 2017, 2018]):
#         print " --> [ERROR] in writing theory systematics. Leaving..."
#         sys.exit(1)

# for syst in experimental_systematics:
#     if not writeSystematic(fdata, df_data, syst, "Merged2Gsf_EE", 125, years = [2016, 2017, 2018]):
#         print " --> [ERROR] in writing theory systematics. Leaving..."
#         sys.exit(1)

# if not writeParamSystematic(fdata, df_data, "Merged2Gsf_EE", 125):
#     print " --> [ERROR] in writing pram systematics. Leaving..."
#     sys.exit(1)

# if not writePdfIndex(fdata, df_data, "Merged2Gsf_EE"):
#     print " --> [ERROR] in writing pdf indices. Leaving..."
#     sys.exit(1)

# fdata.close()