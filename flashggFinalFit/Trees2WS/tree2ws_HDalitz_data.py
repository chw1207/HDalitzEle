import ROOT
import sys
import os
from argparse import ArgumentParser
from collections import OrderedDict as od
from root_numpy import array2tree, tree2array
from commonTools_HDalitz import *


def get_parser():
    parser = ArgumentParser(description = "Script to convert data trees to RooWorkspace (compatible for finalFits)")
    parser.add_argument("-ic", "--inputConfig", help = "Input config: specify list of variables/analysis categories", type = str)
    parser.add_argument("-v",  "--verbosity",   help = "verbose information",       default = False,                  type = bool)

    return parser

# Function to add vars to workspace
def add_vars_to_workspace(_ws, _var_list):
    Mass_lower, Mass_upper = Conf["massBounds"][0], Conf["massBounds"][1]
    Nbins_2 = 0.5 * (Mass_upper - Mass_lower)  # binwidth: 2GeV
    Nbins_1 = 1. * (Mass_upper - Mass_lower)   # binwidth: 1GeV 
    Nbins_p5 = 2. * (Mass_upper - Mass_lower)  # binwidth: 0.5GeV
    var_dict = od()
  
    # Add intLumi var
    intLumi = ROOT.RooRealVar("intLumi", "intLumi", 1000., 0., 999999999.)
    intLumi.setConstant(True)
    getattr(_ws, "import")(intLumi)

    for var in _var_list:
        if var == "CMS_higgs_mass":
            var_dict[var] = ROOT.RooRealVar(var, var, 125, Mass_lower, Mass_upper, "GeV")# initial, lower bound, upper bound
            var_dict[var].setBins(int(Nbins_2))
        if var == "weight":
            var_dict[var] = ROOT.RooRealVar(var, var, 0.)
        
        getattr(_ws,"import")(var_dict[var])

    return var_dict.keys()


# Function to make RooArgSet
def make_argset(_ws, _varNames):
    _aset = ROOT.RooArgSet()
    for v in _varNames: 
        _aset.add(_ws.var(v))
    return _aset


def main():
    # ws separate per year
    for i, f in enumerate(Conf["inputTreeFiles"]):
        print(color.GREEN + "Proccessing the miniTree " + f + color.END)

        # open miniTree file and extract TTree 
        fin = ROOT.TFile.Open(f, "READ")
        if not fin:
            sys.exit(1)
        inTree = ROOT.TTree()
        fin.GetObject(Conf["inputTreeName"], inTree)

        # create root file to store ws
        outdir = extractDirectory(Conf["outputWSFiles"][i])
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        # create the output file
        fout = ROOT.TFile(Conf["outputWSFiles"][i], "RECREATE")
        fout.cd()
        print ("[INFO] Saving the ws file in: {}".format(Conf["outputWSFiles"][i]))

        ws = ROOT.RooWorkspace("w", "work space for background modeling")
        varNames = add_vars_to_workspace(ws, Conf["TreeVars"]) # Add variables to workspace
        aset = make_argset(ws, varNames) # Make argset

        # loop over all category
        for icat, cat in enumerate(Conf["cats"]):
            sel = "({} == {}) & ({} > {}) & ({} < {})".format(Conf["catBranch"], icat+1, Conf["TreeVars"][0], Conf["massBounds"][0], Conf["TreeVars"][0], Conf["massBounds"][1])
            array_cat = tree2array(inTree, branches = Conf["TreeVars"], selection = sel)
            tree_cat = array2tree(array_cat)

            print("[INFO] Total entries of {}: {}".format(cat, tree_cat.GetEntries()))
            # define dataset for cat
            dname = "data_obs_{}".format(cat)  
            dset = ROOT.RooDataSet(dname, dname, aset, Conf["TreeVars"][1])

            # event loop to set the values of variables
            for ev in range(tree_cat.GetEntries()):
                tree_cat.GetEntry(ev)

                mass = getattr(tree_cat, Conf["TreeVars"][0])

                ws.var(Conf["TreeVars"][0]).setVal(mass)
                dset.add(aset, 1.)

            getattr(ws, "import")(dset, ROOT.RooFit.RenameVariable(Conf["TreeVars"][0], Conf["WSVars"][0]))

        if verbose:
            ws.Print()
        ws.Write()
        ws.Delete()
        fout.Close()

    # merged ws
    if (Conf["MergeYears"] == True):
        print(color.GREEN + "Proccessing the miniTrees for full RunII." + color.END)
        chain = ROOT.TChain(Conf["inputTreeName"])
        for i in Conf["inputTreeFiles"]:
            chain.Add(i) 

        outFile = ROOT.TFile(Conf["outputWSFile"], "RECREATE")
        outFile.cd()
        print("Save the WS file in: {}".format(Conf["outputWSFile"]))

        ws = ROOT.RooWorkspace("w", "work space for background modeling")
        varNames = add_vars_to_workspace(ws, Conf["TreeVars"]) # Add variables to workspace
        aset = make_argset(ws, varNames) # Make argset

        for icat, cat in enumerate(Conf["cats"]):
            sel = "({} == {}) & ({} > {}) & ({} < {})".format(Conf["catBranch"], icat+1, Conf["TreeVars"][0], Conf["massBounds"][0], Conf["TreeVars"][0], Conf["massBounds"][1])

            array_cat = tree2array(chain, branches = Conf["TreeVars"], selection = sel)
            tree_cat = array2tree(array_cat)
            print("[INFO] Total entries of {}: {}".format(cat, tree_cat.GetEntries()))

            dname = "data_obs_{}".format(cat)  
            dset = ROOT.RooDataSet(dname, dname, aset, Conf["TreeVars"][1])

            for ev in range(tree_cat.GetEntries()):
                tree_cat.GetEntry(ev)

                ws.var(Conf["TreeVars"][0]).setVal(getattr(tree_cat, Conf["TreeVars"][0]))
                dset.add(aset, 1.)

            getattr(ws, "import")(dset, ROOT.RooFit.RenameVariable(Conf["TreeVars"][0], Conf["WSVars"][0]))

        ws.Write()
        if verbose:
            ws.Print()
        ws.Delete()
        outFile.Close()


if __name__ == "__main__" :
    # Extract information from config file:
    parser = get_parser()
    args = parser.parse_args()
    
    if args.inputConfig is None:
        print("Please specify the config file! eg. config_HDalitz_data")
        parser.print_help()
        sys.exit(1)

    inputConfig = args.inputConfig.split(".")[0] # config_HDalitz_data.py -> config_HDalitz_data
    exec("from "+inputConfig+" import Info as Conf")
    verbose = args.verbosity

    print(color.GREEN + "Mass wundow [{}, {}] GeV".format(Conf["massBounds"][0], Conf["massBounds"][1]) + color.END)

    # turn on the silent mode of RooFit
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    ROOT.RooMsgService.instance().setSilentMode(True)

    state = main()
    sys.exit(state)
    