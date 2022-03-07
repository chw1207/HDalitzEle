import ROOT
import sys
import os
import numpy as np
from argparse import ArgumentParser
from collections import OrderedDict as od
from root_numpy import array2tree, tree2array   
from commonTools_HDalitz import *

# Script to convert minitrees from ggNtuple to RooWorkspace (compatible for finalFits)
def get_parser():
    parser = ArgumentParser(description = "Script to convert data trees to RooWorkspace (compatible for finalFits)")
    parser.add_argument("-ic", "--inputConfig", help = "Input config: specify list of variables/analysis categories", type = str)
    parser.add_argument("-y",  "--year",        help = "Year",                      default = 2017,                   type = int)
    parser.add_argument("-m",  "--mass",        help = "mass point",                default = 125,                    type = int)
    parser.add_argument("-v",  "--verbosity",   help = "verbose information",       default = False,                  type = bool)

    return parser


# Function to add vars to workspace
def add_vars_to_workspace(_ws, _var_list):
    Mass_lower, Mass_upper = Conf["massBounds"][0], Conf["massBounds"][1]

    # Add intLumi var
    intLumi = ROOT.RooRealVar("intLumi", "intLumi", 1000., 0., 999999999.)
    intLumi.setConstant(True)
    getattr(_ws, "import")(intLumi)

    var_dict = od()
    for var in _var_list:
        if var == _var_list[0]:
            var_dict[var] = ROOT.RooRealVar(var, var, 125, Mass_lower, Mass_upper, "GeV") # initial, lower bound, upper bound
            var_dict[var].setBins(int((Mass_upper - Mass_lower) * 2)) # not sure is needed ir not, just keep it
        if var == _var_list[1]: 
            var_dict[var] = ROOT.RooRealVar(var, var, 0.)
        else:
            var_dict[var] = ROOT.RooRealVar(var, var, 1., -999999, 999999)
        
        getattr(_ws,"import")(var_dict[var])
    
    return var_dict.keys()


# Function to make RooArgSet
def make_argset(_ws, _varNames):
    _aset = ROOT.RooArgSet()
    for v in _varNames: 
        _aset.add(_ws.var(v))
    return _aset


def main():
    # read the miniTree
    for i, f in enumerate(Conf["inputTreeFiles"]):
        fin = ROOT.TFile.Open(f, "READ")
        if not fin:
            print(f)
            sys.exit(1)

        tree = ROOT.TTree()
        fin.GetObject(Conf["inputTreeName"], tree)

        # create root file to store ws
        outdir = extractDirectory(Conf["outputWSFile"][i])
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outFile = ROOT.TFile(Conf["outputWSFile"][i], "RECREATE")
        outFile.cd()
        print("[INFO] Save the root file in: {}".format(Conf["outputWSFile"][i]))

        # initiate workspace
        ws = ROOT.RooWorkspace("w", "work space for modelings")
        varNames = add_vars_to_workspace(ws, Conf["TreeVars"]) # Add variables to workspace 
        aset = make_argset(ws, varNames) # Make argset

        # Loop over cats
        for icat, cat in enumerate(Conf["cats"]):
            sel = "({} == {}) & ({} > {}) & ({} < {})".format(Conf["catBranch"], icat+1, Conf["TreeVars"][0], Conf["massBounds"][0], Conf["TreeVars"][0], Conf["massBounds"][1])
            array_cat = tree2array(tree, branches = Conf["TreeVars"], selection = sel)
            tree_cat = array2tree(array_cat)
            
            # Add dataset to worksapce
            dname = "set_%i_%s" %(Conf["MassPoint"], cat)
            dset = ROOT.RooDataSet(dname, dname, tree_cat, aset, "", Conf["TreeVars"][1])
            # rename the tree varables to the name used in the following modelings
            getattr(ws, "import")(dset, ROOT.RooFit.RenameVariable(Conf["TreeVars"][0], Conf["WSVars"][0]))

        ws.Write()
        if verbose:
            ws.Print()

        outFile.Close()


if __name__ == "__main__" :
    # extract information from config file:
    parser = get_parser()
    args = parser.parse_args()
    
    inputConfig = args.inputConfig.split(".")[0] # config_HDalitz.py -> config_HDalitz
    if inputConfig is None:
        print("Please specify the config file! eg. config_HDalitz.py")
        parser.print_help()
        sys.exit(1)
    exec("from "+inputConfig+" import Info")

    mass, year, verbose = args.mass, args.year, args.verbosity
    Conf = Info[mass][year] # extract the information from the config file

    print(color.GREEN + "Tree2WS -> year: {}, mass: {}GeV".format(year, mass) + color.END)
    print(color.GREEN + "Mass wundow [{}, {}] GeV".format(Conf["massBounds"][0], Conf["massBounds"][1]) + color.END)

    # turn on the silent mode of RooFit
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    ROOT.RooMsgService.instance().setSilentMode(True)
    
    # main script
    state = main()
    sys.exit(state)
