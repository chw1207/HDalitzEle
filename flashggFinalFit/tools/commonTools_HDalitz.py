import os, sys
import re
import ROOT
import math
from glob import glob
from collections import OrderedDict as od

class color:
   PURPLE = "\033[95m"
   CYAN = "\033[96m"
   DARKCYAN = "\033[36m"
   BLUE = "\033[94m"
   GREEN = "\033[92m"
   YELLOW = "\033[93m"
   RED = "\033[91m"
   BOLD = "\033[1m"
   UNDERLINE = "\033[4m"
   END = "\033[0m"

# Function for iterating over ROOT argsets in workspace
# https://root-forum.cern.ch/t/iterating-over-rooargset/16331/2
def rooiter(x):
    iter = x.iterator()
    ret = iter.Next()
    while ret:
        yield ret
        ret = iter.Next()

def extractWSFileNames(_inputWSDir):
    state = False
    if not os.path.isdir(_inputWSDir):
        print("[ERROR] No such directory...")
        return False
    else:
        return glob("{}/*.root".format(_inputWSDir))

# asssume the structure of the file is workspace_HDalitz_sigMC_{mass}_{cat}_{proc}.root
def extractListOfProcs(_listOfWSFileNames):
    procs = []
    for pName in _listOfWSFileNames:
        f = pName.split("/")[-1] # extract the file name (without the directroy name)
        p = f.split("_")[-1].split(".root")[0] 
        if p not in procs:
            procs.append(p)
    procs.sort(key = str.lower) # Define list of procs(alphabetically ordered)
    return procs

def extractListOfMass(_listOfWSFileNames):
    mass = []
    for pName in _listOfWSFileNames:
        f = pName.split("/")[-1] # extract the file name (without the directroy name)
        m = f.split("_")[3]
        if m not in mass:
            mass.append(m) 
    mass.sort()
    return mass

# function to extract the directory form the full path
# "./test/test.root" -> "./test"
def extractDirectory(fileName):
    vdirs = fileName.split("/")
    del vdirs[-1]
    direc = "/".join(vdirs)
    return direc

# function to converte process to production mode in dataset name
procToDatacardNameMap = od()
procToDatacardNameMap["ggF"] = "ggH"
procToDatacardNameMap["VBF"] = "qqH"
procToDatacardNameMap["WH"] = "WH"
procToDatacardNameMap["ZH"] = "ZH"
def procToDatacardName(_proc):
    k = _proc.split("_")[0]
    if k in procToDatacardNameMap: 
        _proc = re.sub(k, procToDatacardNameMap[k], _proc)
    return _proc


# function to map the cat number (use to create datacard)
catNumMap = od({
    "Merged2Gsf_HVBF": 0,
    "Merged2Gsf_LVBF": 1,
    "Merged2Gsf_BST": 2,
    "Merged2Gsf_EBHR9": 3,
    "Merged2Gsf_EBLR9": 4,
    "Merged2Gsf_EE": 5,
    "Merged1Gsf_HVBF": 6, 
    "Merged1Gsf_LVBF": 7,
    "Merged1Gsf_BST": 8,
    "Merged1Gsf_EBHR9": 9, 
    "Merged1Gsf_EBLR9": 10,
    "Merged1Gsf_EE": 11,
    "Resolved": 12
})

M2Untag = ["Merged2Gsf_EBHR9", "Merged2Gsf_EBLR9", "Merged2Gsf_EE"]
M1Untag = ["Merged1Gsf_EBHR9", "Merged1Gsf_EBLR9", "Merged1Gsf_EE"]
M2tag = ["Merged2Gsf_HVBF", "Merged2Gsf_LVBF", "Merged2Gsf_BST"]
M1tag = ["Merged1Gsf_HVBF", "Merged1Gsf_LVBF", "Merged1Gsf_BST"]

# Base mass point list for signal samples
massBaseList = [120, 125, 130]
massList = [120+i for i in range(11)]
years = [2016, 2017, 2018]