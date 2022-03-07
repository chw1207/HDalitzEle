# Script to exrtract yields and some useful information for creating datacard
#  * Uses Pandas dataframe to store all proc x cat yields
#  ! Systematic variations are yet included (To-do)


import os, sys
import ROOT
import pandas as pd
import pickle
import math
from glob import glob
from argparse import ArgumentParser
from collections import OrderedDict as od
# from systematics import theory_systematics, experimental_systematics, signal_shape_systematics
# from commonObjects import *
from commonTools_HDalitz import *

def get_parser():
    parser = ArgumentParser(description = "Script to make dataframe for creating datacard")
    parser.add_argument("-sigWSDirMap", "--SigWSDirMap", 
        help = "Map. Format: year = signal model and dataset ws (separate years by comma)", 
        default = "2016:../Signal/workspace/Interpolation/2016,2017:../Signal/workspace/Interpolation/2017,2018:../Signal/workspace/Interpolation/2018", type = str
    )
    parser.add_argument("-bkgWsDir", "--bkgModelWSDir", help = "background model ws", default = "../Background/MultiPdf", type = str)
    parser.add_argument("-e", "--ext", help = "ext string", default = "DiPho_13TeV", type = str)
    parser.add_argument("-cat", "--category", help = "Analysis category", type = str)
    parser.add_argument("-dc", "--decayMode", help = "decay mode", default = "heeg", type = str)
    parser.add_argument("-my", "--mergeYears", help = "Merge category across years", default = True, type = bool)

    return parser


def main():
    # Extract years and signal ws directory
    SigWSDirMap = od()
    for i in SigWSDirMapstr.split(","):
        SigWSDirMap[i.split(":")[0]] = i.split(":")[1]

    years = SigWSDirMap.keys()

    procsMap, massMap = od(), od()
    for year, iWSDir in SigWSDirMap.iteritems():
        WSFileNamesList = extractWSFileNames(iWSDir)
        if not WSFileNamesList:
            sys.exit(1)
        procsMap[year] = extractListOfProcs(WSFileNamesList)
        massMap[year] = extractListOfMass(WSFileNamesList)

    # Initiate pandas dataframe
    column_title = ["year", "type", "procOriginal", "proc", "cat", "mass", "modelWSFile", "model", "rate"]
    df_data = pd.DataFrame(columns = column_title)
    
    # FILL DATAFRAME: signal
    print("[INFO] Adding signal to dataFrame")
    for year in years:
        for mass in massMap[year]:
            for proc in procsMap[year]:
                # _id = "{}_{}_{}".format(proc, year, cat) # Identifier

                _procOriginal = proc
                _proc = "{}_{}_{}".format(procToDatacardName(proc), year, decayMode)
                
                if mergeYears:
                    _cat = cat
                else:
                    _cat = "{}_{}".format(cat, year)
                
                _modelWSFile = glob("{}/*_{}_{}_{}.root".format(SigWSDirMap[year], mass, cat, proc))[0]
                _model = "{}:newSigPdf_{}_{}_{}".format(_inputWSName, _procOriginal, year, cat)
                # _model = "{}:SigPdf_{}".format(_inputWSName, cat)

                _rate = float(lumiMap[year]) * 1000

                df_data.loc[len(df_data)] = [year, "sig", _procOriginal, _proc, _cat, mass, _modelWSFile, _model, _rate]

    # FILL DATAFRAME: background
    print("[INFO] Adding background/data to dataFrame")
    _proc_bkg = "bkg_mass"
    _proc_data = "data_obs"
    if mergeYears:
        _cat = cat
        _modelWSFile = "{}/Multipdf.root".format(bkgModelWSDir)
        _model_bkg = "{}:CMS_higgs_{}_{}_bkgshape".format(_bkgWSName, _cat, _ext)
        _model_data = "{}:roohist_data_mass_{}".format(_bkgWSName, _cat)
        _mass = "-" # not needed for data/bkg
        df_data.loc[len(df_data)] = ["merged", "bkg", _proc_bkg, _proc_bkg, _cat, _mass, _modelWSFile, _model_bkg, 1] #!FIXEDME
        df_data.loc[len(df_data)] = ["merged", "data", _proc_data, _proc_data, _cat, _mass, _modelWSFile, _model_data, -1]
    else:
        for year in years:
            _cat = "{}_{}".format(cat, year)
            _catStripYear = cat
            _modelWSFile = "{}/{}/Multipdf_{}.root".format(bkgModelWSDir, year, year)
            _model_bkg = "{}:CMS_higgs_{}_{}_bkgshape".format(_bkgWSName, _catStripYear, _ext)
            _model_data = "{}:roohist_data_mass_{}".format(_bkgWSName, _catStripYear)
            _mass = "-" # not needed for data/bkg

            df_data.loc[len(df_data)] = [year, "bkg", _proc_bkg,  _proc_bkg, _cat, _mass, _modelWSFile, _model_bkg, 1] #!FIXEDME
            df_data.loc[len(df_data)] = [year, "data", _proc_data, _proc_data, _cat, _mass, _modelWSFile, _model_data, -1]

    # Yields: for each signal row in dataFrame extract the yield
    # Loop over signal rows in dataFrame: extract yields (nominal & systematic variations)
    # totalSignalRows = float(df_data[df_data["type"] == "sig"].shape[0])
    # print("[INFO] Adding signal yields to dataFrame")
    df_data["nominal_yield"] = "-"
    for ir, r in df_data[df_data["type"] == "sig"].iterrows():
        # open input WS file and extract workspace 
        fin = ROOT.TFile.Open(r.modelWSFile)
        if not fin:
            sys.exit(1)
        inputWS = ROOT.RooWorkspace()
        fin.GetObject(_inputWSName, inputWS)
        
        # Extract nominal yield
        # _yield = inputWS.var("ExpYield_{}".format(r.mass)).getVal() #!FIXEDME
        _yield = inputWS.var("ExpYield").getVal()
        df_data.at[ir, "nominal_yield"] = _yield

        # Remove the workspace and file from heap
        inputWS.Delete()
        fin.Close()

    # SAVE YIELDS DATAFRAME
    print ("[INFO] Saving yields dataframe: ./yields/{}.pkl".format(cat))
    if not os.path.isdir("./yields"):
        os.makedirs("./yields")
    with open("./yields/{}.pkl".format(cat), "wb") as fout:
        pickle.dump(df_data, fout)

    y = [round(df_data.query("type == 'sig' and mass == '{}'".format(m))["nominal_yield"].sum(), 3) for m in massList]
    print("----> yield list: {}".format(y))
    print("")

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    SigWSDirMapstr, bkgModelWSDir = args.SigWSDirMap, args.bkgModelWSDir
    cat, mergeYears, decayMode = args.category, args.mergeYears, args.decayMode
    _ext = args.ext
    
    lumiMap = {"2016": 36.33, "2017": 41.48, "2018": 59.35, "combined": 137.17, "merged": 137.17}
    _inputWSName = "w"
    _bkgWSName = "multipdf"

    print(color.GREEN + "Processing the category " + cat + "......" + color.END)
    state = main()
    sys.exit(state)
