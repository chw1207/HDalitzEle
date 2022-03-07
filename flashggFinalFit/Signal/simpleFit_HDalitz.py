import ROOT
import sys
import os
import numpy as np
import json
from argparse import ArgumentParser
from collections import OrderedDict as od
from tools.CMS_lumi import CMS_lumi
from commonTools_HDalitz import *

def get_parser():
    parser = ArgumentParser(description = "Script to fit the signal shape")
    parser.add_argument("-if", "--inputFile", help = "Input file contains work space",  type = str)
    parser.add_argument("-iw", "--inputWS",   help = "Input work space name",           default = "w", type = str)
    parser.add_argument("-id", "--inputdSet", help = "Input dataset in the work space", type = str)
    parser.add_argument("-proc", "--process", help = "Process of the signal sample(ggF, VBF, WH, ZH)", type = str)
    parser.add_argument("-y", "--year",       help = "year of the signal sample",       type = int)
    
    return parser

# extract the mass point and category name frome the dataset name
# compatible with tree2ws_HDalitz.py
def get_mass_catName(dname):
    vname = dname.split("_")
    mass = vname[1]

    vname = vname[2:]
    catName = "_".join(vname)
    return mass, catName

# plot fuction
def DrawSingleFit(CMS_higgs_mass, dataset, SigPdf, fitRes, outName, massLegName, massPoint, catName, process):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    Nbins = 170 - 105
    xframe = CMS_higgs_mass.frame(105, 170)
    
    dataset.plotOn(xframe, ROOT.RooFit.Name("set"), ROOT.RooFit.Binning(Nbins), ROOT.RooFit.MarkerStyle(ROOT.kFullCircle), ROOT.RooFit.MarkerSize(1.5), ROOT.RooFit.XErrorSize(10e-5), ROOT.RooFit.DataError(ROOT.RooDataSet.SumW2))
    SigPdf.plotOn(xframe, ROOT.RooFit.Name("sigFit"), ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#5893D4")), ROOT.RooFit.LineWidth(4))
    SigPdf.plotOn(xframe, ROOT.RooFit.Components("SigDCB"), ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#1C7747")), ROOT.RooFit.LineStyle(7), ROOT.RooFit.Name("DCB"), ROOT.RooFit.LineWidth(4))
    SigPdf.plotOn(xframe, ROOT.RooFit.Components("SigGauss"), ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#E23E57")), ROOT.RooFit.LineStyle(7), ROOT.RooFit.Name("Gauss"), ROOT.RooFit.LineWidth(4))
    dataset.plotOn(xframe, ROOT.RooFit.Name("set"), ROOT.RooFit.Binning(Nbins), ROOT.RooFit.MarkerStyle(ROOT.kFullCircle), ROOT.RooFit.MarkerSize(1.5), ROOT.RooFit.XErrorSize(10e-5), ROOT.RooFit.DataError(ROOT.RooDataSet.SumW2))
    
    xframe.SetTitle("")
    xframe.GetXaxis().SetTickSize(0.03)
    xframe.GetXaxis().SetTitleSize(0.04)
    xframe.GetXaxis().SetLabelSize(0.04)
    xframe.GetXaxis().SetLabelOffset(0.02)
    xframe.GetXaxis().SetTitleOffset(1.4)
    xframe.GetXaxis().SetTitle(massLegName)
    xframe.GetYaxis().SetTitle("Signal shape / ({} GeV)".format((170 - 105)/Nbins))
    xframe.GetYaxis().SetNdivisions(510)
    xframe.GetYaxis().SetTickSize(0.03)
    xframe.GetYaxis().SetTitleSize(0.04)
    xframe.GetYaxis().SetLabelSize(0.04)
    xframe.GetYaxis().SetTitleOffset(1.6)
    xframe.SetMaximum(xframe.GetMaximum() * 500.)
    xframe.SetMinimum(xframe.GetMaximum() * 0.00000001)

    c = ROOT.TCanvas("c", "", 900, 900)
    c.cd()
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.14)
    c.SetBottomMargin(0.12)
    c.SetLogy()
    xframe.Draw()

    CMS_lumi(c, 4, 11, "", year, True, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma ", "")
    c.Update()
    c.RedrawAxis()
    
    ltx = ROOT.TLatex()
    ltx.SetNDC()
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.037)
    ltx.DrawLatex(0.53, 0.86, "{}, {}".format(process, catName))

    leg1 = ROOT.TLegend(0.53, 0.7, 0.9, 0.83)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.037)
    leg1.SetFillColor(0)
    leg1.SetLineColor(0)
    leg1.AddEntry(xframe.findObject("set"), "Simulation", "ep")
    leg1.AddEntry(xframe.findObject("sigFit"), "Parametric model", "l")
    # leg1.AddEntry((TObject *)0, Form("-log(L) = %g", minNll), "")
    leg1.Draw()

    leg2 = ROOT.TLegend(0.61, 0.62, 0.9, 0.7)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.033)
    leg2.SetFillColor(0)
    leg2.SetLineColor(0)
    leg2.AddEntry(xframe.findObject("DCB"), "DCB", "l")
    leg2.AddEntry(xframe.findObject("Gauss"), "Gauss", "l")
    leg2.Draw("same")

    c.Print(outName)
    del c


#############################
##                         ##
##    Main fitting part    ##
##                         ##
#############################
ROOT.gROOT.SetBatch() # PyROOT does not display any graphics(root "-b" option)

parser = get_parser()
args = parser.parse_args()
infile, inWsName, indName = args.inputFile, args.inputWS, args.inputdSet
proc, year = args.process, args.year

# turn on the silent mode of RooFit
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
ROOT.RooMsgService.instance().setSilentMode(True)

# open ROOT file storing datasets
fin = ROOT.TFile.Open(infile, "READ")
if not fin:
    sys.exit(1)
inWS = ROOT.RooWorkspace()
fin.GetObject(inWsName, inWS)

dataName = [] # double check if the name of the dataset is in the ws.
for i in inWS.allData():
    dataName.append(i.GetName())
if indName not in dataName:
    print("RooDataSet should be the following!")
    print(json.dumps(dataName, indent = 4)) # print the list gracefully using the json structure
    sys.exit(1)
# data = inWS.data(indName)
oldData = inWS.data(indName)

mass, catName = get_mass_catName(indName) # mass = 120 or 125 or 130

# extract the fitting variable
CMS_higgs_mass = inWS.var("CMS_higgs_mass")
CMS_higgs_mass.setUnit("GeV")
CMS_higgs_mass.setMin(105) #! do not remove this line
CMS_higgs_mass.setMax(170) #! do not remove this line
CMS_higgs_mass.setRange("NormRange", 105, 170)

#! FIXEDME: No available signal samples for 2016 and 2018 for now, so they are mimicked by 2017 signal samples 
NewWei = ROOT.RooRealVar("NewWei", "NewWei", -10000, 10000)
rewei = [(35.9/41.525), 1., (59.7/41.525)] # 2016, 2017, 2018(luminosity is used to reweight)
data = oldData.emptyClone()
data.SetName("set")
for i in range(oldData.numEntries()):
    CMS_higgs_mass.setVal(oldData.get(i).getRealValue(CMS_higgs_mass.GetName()))
    NewWei.setVal(rewei[year - 2016] * oldData.weight())
    data.add(ROOT.RooArgSet(CMS_higgs_mass, NewWei), NewWei.getVal())

# Build DCB 
mean = ROOT.RooRealVar("mean", " ", float(mass), float(mass) - 3, float(mass) + 3, "GeV")
sigma_CB = ROOT.RooRealVar("sigma_CB", " ", 2.5 , 0.1, 5.)
alpha1 = ROOT.RooRealVar("alpha1", " ", 1, 0, 15) 
alpha2 = ROOT.RooRealVar("alpha2", " ", 2, 0, 15)
n1 = ROOT.RooRealVar("n1", " ", 5, 0, 10)
n2 = ROOT.RooRealVar("n2", " ", 3, 0, 10)
DCBPdf = ROOT.RooDoubleCB("SigDCB", "Double sided CB", CMS_higgs_mass, mean, sigma_CB, alpha1, n1, alpha2, n2)

# Build Gaussian
# Gaussian share the same mean as DCB
sigma_Gauss = ROOT.RooRealVar("sigma_Gauss", " ", 20., 18, 60)
GaussPdf = ROOT.RooGaussian("SigGauss", "Gaussian", CMS_higgs_mass, mean, sigma_Gauss)

# Build Final pdf
frac_DCB = ROOT.RooRealVar("frac_DCB", " ", 0.5, 0.3, 0.9999999) # fraction of DCB
SigPdf = ROOT.RooAddPdf("SigPdf", " ", DCBPdf, GaussPdf, frac_DCB)

# unbinned likelihood fit
outStat = "./status_check"
if not os.path.exists(outStat):
    os.makedirs(outStat)
Fitresult = SigPdf.fitTo(data, ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Range("NormRange"), ROOT.RooFit.Minimizer("Minuit","minimize"), ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.PrintLevel(-1))
if (Fitresult.status() != 0):
    print(color.YELLOW + "[WARNING] Fit status = {}".format(Fitresult.status()) + color.END)
    with open("{}/WARNSTATUS_{}_{}_{}.txt".format(outStat, mass, catName, proc), "w") as f:
        f.write("Refit the category\n")
        f.write("mass = {}, cat = {}, proc = {}, year = {}\n".format(mass, catName, proc, year))

# plot the fit results
outPlot = "./plots/{}".format(year)
if not os.path.exists(outPlot):
    os.makedirs(outPlot)
outName = "{}/singleFit_{}_{}_{}.pdf".format(outPlot, mass, catName, proc)
DrawSingleFit(CMS_higgs_mass, data, SigPdf, Fitresult, outName, "M_{ee#gamma} [GeV]", int(mass), catName, proc)

# save the fit results
outRes = "./FitResults/{}".format(year)
if not os.path.exists(outRes):
    os.makedirs(outRes)
print("[INFO]: Save RooFitResult in {}/fit_{}_{}_{}.root".format(outRes, mass, catName, proc))
fres = ROOT.TFile("{}/fit_{}_{}_{}.root".format(outRes, mass, catName, proc), "RECREATE")
fres.cd()
Fitresult.SetName("fitres")
Fitresult.Write()
fres.Close()
fres.Delete("*;*")

# build the work space
outws = "./workspace/{}".format(year)
if not os.path.exists(outws):
    os.makedirs(outws)
fws = ROOT.TFile("{}/workspace_HDalitz_sigMC_{}_{}_{}.root".format(outws, mass, catName, proc), "RECREATE")
fws.cd()
ws = ROOT.RooWorkspace("w")
getattr(ws, "import")(data)
getattr(ws, "import")(SigPdf)
params = SigPdf.getParameters(data)
ws.saveSnapshot("nominal_values", params)
print("[INFO]: Save workspace in {}/workspace_HDalitz_sigMC_{}_{}_{}.root".format(outws, mass, catName, proc))
ws.Write()
fws.Close()
fws.Delete("*;*")
ws.Delete() #! do not remove this line otherwise script will break when you import data(https://sft.its.cern.ch/jira/browse/ROOT-9890)
