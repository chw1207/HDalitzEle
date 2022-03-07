import ROOT
import sys
import os
import numpy as np
from argparse import ArgumentParser
from collections import OrderedDict as od
from tools.CMS_lumi import CMS_lumi
from commonTools_HDalitz import *

def get_parser():
    parser = ArgumentParser(description = "Script to plot the interpolated signal shape")
    parser.add_argument("-c", "--category", help = "category", type = str)
    parser.add_argument("-y", "--year", help = "year", type = int)
    parser.add_argument("-p", "--process", help = "process", type = str)
    parser.add_argument("-fr120", "--FitReslutsFileM120", help = "fit results of mass 120 GeV", type = str)
    parser.add_argument("-fr125", "--FitReslutsFileM125", help = "fit results of mass 125 GeV", type = str)
    parser.add_argument("-fr130", "--FitReslutsFileM130", help = "fit results of mass 130 GeV", type = str)
    parser.add_argument("-fw120", "--WSFileM120", help = "workspace of mass 120 GeV", type = str)
    parser.add_argument("-fw125", "--WSFileM125", help = "workspace of mass 125 GeV", type = str)
    parser.add_argument("-fw130", "--WSFileM130", help = "workspace of mass 130 GeV", type = str)
    
    return parser

def getNomFitres(fname, fwname, masspoint, cat, year, proc): # get the fit results of masses 120, 125 ,130
    Fitparams = od()

    f = ROOT.TFile.Open(fname, "READ")
    if not f:
        sys.exit(1)
    fitres = ROOT.RooFitResult()
    f.GetObject("fitres", fitres)

    # extract DCB parameters
    mean = fitres.floatParsFinal().find("mean")
    sigma_CB = fitres.floatParsFinal().find("sigma_CB")
    alpha1 = fitres.floatParsFinal().find("alpha1")
    alpha2 = fitres.floatParsFinal().find("alpha2")
    n1 = fitres.floatParsFinal().find("n1")
    n2 = fitres.floatParsFinal().find("n1")
    frac_DCB = fitres.floatParsFinal().find("frac_DCB")

    # extract Gauss parameters
    sigma_Gauss = fitres.floatParsFinal().find("sigma_Gauss")

    # build a dictionary to store the value
    Fitparams["mean_{}_{}_{}".format(masspoint, cat, proc)] = np.array([mean.getValV(), mean.getError()])
    Fitparams["sigma_CB_{}_{}_{}".format(masspoint, cat, proc)] = np.array([sigma_CB.getValV(), sigma_CB.getError()])
    Fitparams["alpha1_{}_{}_{}".format(masspoint, cat, proc)] = np.array([alpha1.getValV(), alpha1.getError()])
    Fitparams["alpha2_{}_{}_{}".format(masspoint, cat, proc)] = np.array([alpha2.getValV(), alpha2.getError()])
    Fitparams["n1_{}_{}_{}".format(masspoint, cat, proc)] = np.array([n1.getValV(), n1.getError()])
    Fitparams["n2_{}_{}_{}".format(masspoint, cat, proc)] = np.array([n2.getValV(), n2.getError()])
    Fitparams["frac_DCB_{}_{}_{}".format(masspoint, cat, proc)] = np.array([frac_DCB.getValV(), frac_DCB.getError()])
    Fitparams["sigma_Gauss_{}_{}_{}".format(masspoint, cat, proc)] = np.array([sigma_Gauss.getValV(), sigma_Gauss.getError()])

    f.Close()

    fw = ROOT.TFile.Open(fwname, "READ")
    if not fw:
        sys.exit(1)
    inWS = ROOT.RooWorkspace()
    fw.GetObject("w", inWS)
    data = inWS.data("set")
    Fitparams["ExpYield_{}_{}_{}".format(masspoint, cat, proc)] = np.array([data.sumEntries(), data.numEntries()])
    
    fw.Close()

    return Fitparams

# calculate the values using Interpolation
def getInterpolation(massHigh, massLow, mass, FitParamsHigh, FitParamsLow, cat, proc):
    params = od()
    
    a, b = 0, 0
    if ((mass == 120) or (mass == 125) or ((mass == 130))):
        a, b = 1, 0
    else:
        a = (massHigh - mass)/(massHigh - massLow)
        b = 1 - a

    mean_intp = a * FitParamsLow["mean_{}_{}_{}".format(int(massLow), cat, proc)][0] + b * FitParamsHigh["mean_{}_{}_{}".format(int(massHigh), cat, proc)][0]
    sigma_CB_intp = a * FitParamsLow["sigma_CB_{}_{}_{}".format(int(massLow), cat, proc)][0] + b * FitParamsHigh["sigma_CB_{}_{}_{}".format(int(massHigh), cat, proc)][0]
    alpha1_intp = a * FitParamsLow["alpha1_{}_{}_{}".format(int(massLow), cat, proc)][0] + b * FitParamsHigh["alpha1_{}_{}_{}".format(int(massHigh), cat, proc)][0]
    alpha2_intp = a * FitParamsLow["alpha2_{}_{}_{}".format(int(massLow), cat, proc)][0] + b * FitParamsHigh["alpha2_{}_{}_{}".format(int(massHigh), cat, proc)][0]
    n1_intp = a * FitParamsLow["n1_{}_{}_{}".format(int(massLow), cat, proc)][0] + b * FitParamsHigh["n1_{}_{}_{}".format(int(massHigh), cat, proc)][0]
    n2_intp = a * FitParamsLow["n2_{}_{}_{}".format(int(massLow), cat, proc)][0] + b * FitParamsHigh["n2_{}_{}_{}".format(int(massHigh), cat, proc)][0]
    frac_DCB_intp = a * FitParamsLow["frac_DCB_{}_{}_{}".format(int(massLow), cat, proc)][0] + b * FitParamsHigh["frac_DCB_{}_{}_{}".format(int(massHigh), cat, proc)][0]
    sigma_Gauss_intp = a * FitParamsLow["sigma_Gauss_{}_{}_{}".format(int(massLow), cat, proc)][0] + b * FitParamsHigh["sigma_Gauss_{}_{}_{}".format(int(massHigh), cat, proc)][0]

    norm_intp = a * FitParamsLow["ExpYield_{}_{}_{}".format(int(massLow), cat, proc)][0] + b * FitParamsHigh["ExpYield_{}_{}_{}".format(int(massHigh), cat, proc)][0] # used to normalize
    # print(norm_intp, FitParamsLow["ExpYield_{}_{}_{}".format(int(massHigh), cat, proc)][0], FitParamsHigh["ExpYield_{}_{}_{}".format(int(massHigh), cat, proc)][0], a, b)

    # calculate the errors using Interpolation
    mean_err_intp = np.sqrt(pow(a * FitParamsLow["mean_{}_{}_{}".format(int(massLow), cat, proc)][1], 2) + pow(b * FitParamsHigh["mean_{}_{}_{}".format(int(massHigh), cat, proc)][1], 2))
    sigma_CB_err_intp = np.sqrt(pow(a * FitParamsLow["mean_{}_{}_{}".format(int(massLow), cat, proc)][1], 2) + pow(b * FitParamsHigh["mean_{}_{}_{}".format(int(massHigh), cat, proc)][1], 2))
    alpha1_err_intp = np.sqrt(pow(a * FitParamsLow["mean_{}_{}_{}".format(int(massLow), cat, proc)][1], 2) + pow(b * FitParamsHigh["mean_{}_{}_{}".format(int(massHigh), cat, proc)][1], 2))
    alpha2_err_intp = np.sqrt(pow(a * FitParamsLow["mean_{}_{}_{}".format(int(massLow), cat, proc)][1], 2) + pow(b * FitParamsHigh["mean_{}_{}_{}".format(int(massHigh), cat, proc)][1], 2))
    n1_err_intp = np.sqrt(pow(a * FitParamsLow["mean_{}_{}_{}".format(int(massLow), cat, proc)][1], 2) + pow(b * FitParamsHigh["mean_{}_{}_{}".format(int(massHigh), cat, proc)][1], 2))
    n2_err_intp = np.sqrt(pow(a * FitParamsLow["mean_{}_{}_{}".format(int(massLow), cat, proc)][1], 2) + pow(b * FitParamsHigh["mean_{}_{}_{}".format(int(massHigh), cat, proc)][1], 2))
    frac_DCB_err_intp = np.sqrt(pow(a * FitParamsLow["mean_{}_{}_{}".format(int(massLow), cat, proc)][1], 2) + pow(b * FitParamsHigh["mean_{}_{}_{}".format(int(massHigh), cat, proc)][1], 2))
    sigma_Gauss_err_intp = np.sqrt(pow(a * FitParamsLow["mean_{}_{}_{}".format(int(massLow), cat, proc)][1], 2) + pow(b * FitParamsHigh["mean_{}_{}_{}".format(int(massHigh), cat, proc)][1], 2))

    params["mean"],         params["mean_err"]          = mean_intp,        mean_err_intp
    params["sigma_CB"],     params["sigma_CB_err"]      = sigma_CB_intp,    sigma_CB_err_intp
    params["alpha1"],       params["alpha1_err"]        = alpha1_intp,      alpha1_err_intp
    params["alpha2"],       params["alpha2_err"]        = alpha2_intp,      alpha2_err_intp
    params["n1"],           params["n1_err"]            = n1_intp,          n1_err_intp
    params["n2"],           params["n2_err"]            = n2_intp,          n2_err_intp
    params["frac_DCB"],     params["frac_DCB_err"]      = frac_DCB_intp,    frac_DCB_err_intp
    params["sigma_Gauss"],  params["sigma_Gauss_err"]   = sigma_Gauss_intp, sigma_Gauss_err_intp
    params["norm"]                                      = norm_intp

    return params

def get_electron_cat(catName): # Merged2Gsf_EBLR9 -> ["Merged2Gsf", "EBLR9"]
    vname = catName.split("_")
    return vname

#############################
##                         ##
##      Main script        ##
##                         ##
#############################
ROOT.gROOT.SetBatch() # PyROOT does not display any graphics(root "-b" option)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
ROOT.RooMsgService.instance().setSilentMode(True)

parser = get_parser()
args = parser.parse_args()
cat, year, proc = args.category, args.year, args.process
fr120, fr125, fr130 = args.FitReslutsFileM120, args.FitReslutsFileM125, args.FitReslutsFileM130
fw120, fw125, fw130 = args.WSFileM120, args.WSFileM125, args.WSFileM130

base_mass = [120, 125, 130]
mass_interval = 1
v_mass = np.array([120. + mass_interval * i for i in range(int(((base_mass[2] - base_mass[0])/mass_interval) + 1))])

# get the fit results of masses 120, 125 ,130
FitResM120 = getNomFitres(fr120, fw120, 120, cat, year, proc)
FitResM125 = getNomFitres(fr125, fw125, 125, cat, year, proc)
FitResM130 = getNomFitres(fr130, fw130, 130, cat, year, proc)

# ROOT plotting style
ROOT.gStyle.SetPadTickX(1)
ROOT.gStyle.SetPadTickY(1)
ROOT.gStyle.SetOptStat(0)

# create RooRealVar to fit
CMS_higgs_mass = ROOT.RooRealVar("CMS_higgs_mass", "CMS_higgs_mass", 105, 145, "GeV")
CMS_higgs_mass.setRange("NormRange", 105, 145)
CMS_higgs_mass.setMin(105) 
CMS_higgs_mass.setMax(145) 
xframe = CMS_higgs_mass.frame(110, 140)

for imass, mass in enumerate(v_mass):
    FitParamsHigh, FitParamsLow = od(), od ()

    if (mass == 120) or (mass == 125) or (mass == 130):
        if (mass == 120):
            FitParamsHigh,  FitParamsLow, massHigh, massLow = FitResM120, FitResM120, 120, 120
        if (mass == 125):
            FitParamsHigh,  FitParamsLow, massHigh, massLow = FitResM125, FitResM125, 125, 125
        if (mass == 130):
            FitParamsHigh,  FitParamsLow, massHigh, massLow = FitResM130, FitResM130, 130, 130
    elif mass < 125:
        FitParamsHigh, FitParamsLow, massHigh, massLow = FitResM125, FitResM120, 125, 120
    else:
        FitParamsHigh, FitParamsLow, massHigh, massLow = FitResM130, FitResM125, 130, 125

    intp_params = getInterpolation(massHigh, massLow, mass, FitParamsHigh, FitParamsLow, cat, proc)

    mean = ROOT.RooRealVar("mean_{}".format(cat), " ", intp_params["mean"], "GeV")
    sigma_CB = ROOT.RooRealVar("sigma_CB_{}".format(cat), " ", intp_params["sigma_CB"], "GeV")
    alpha1 = ROOT.RooRealVar("alpha1_{}".format(cat), " ", intp_params["alpha1"]) 
    alpha2 = ROOT.RooRealVar("alpha2_{}".format(cat), " ", intp_params["alpha2"])
    n1 = ROOT.RooRealVar("n1_{}".format(cat), " ", intp_params["n1"])
    n2 = ROOT.RooRealVar("n2_{}".format(cat), " ", intp_params["n2"])
    DCBPdf = ROOT.RooDoubleCBFast("SigDCB_{}".format(cat), "Double sided CB", CMS_higgs_mass, mean, sigma_CB, alpha1, n1, alpha2, n2)

    # Build Gaussian
    # Gaussian share the same mean as DCB
    sigma_Gauss = ROOT.RooRealVar("sigma_Gauss_{}".format(cat), " ", intp_params["sigma_Gauss"], "GeV")
    GaussPdf = ROOT.RooGaussian("SigGauss_{}".format(cat), "Gaussian", CMS_higgs_mass, mean, sigma_Gauss)

    # Build Final pdf
    frac_DCB = ROOT.RooRealVar("frac_DCB_{}".format(cat), " ", intp_params["frac_DCB"]) # fraction of DCB

    # set the errors of parameters for signal pdf
    mean.setError(intp_params["mean_err"])
    sigma_CB.setError(intp_params["sigma_CB_err"])
    alpha1.setError(intp_params["alpha1_err"])
    alpha2.setError(intp_params["alpha2_err"])
    n1.setError(intp_params["n1_err"])
    n2.setError(intp_params["n2_err"])
    sigma_Gauss.setError(intp_params["sigma_Gauss_err"])
    frac_DCB.setError(intp_params["frac_DCB_err"])

    sigpdf = ROOT.RooAddPdf("SigPdf_{}".format(cat), " ", DCBPdf, GaussPdf, frac_DCB)
    if (mass != 120) and (mass != 125) and (mass != 130):
        sigpdf.plotOn(
            xframe, ROOT.RooFit.Range("NormRange"), 
            ROOT.RooFit.Normalization(intp_params["norm"], ROOT.RooAbsReal.NumEvent), 
            ROOT.RooFit.LineColor(ROOT.TColor.GetColorPalette(imass * 20)), 
            ROOT.RooFit.LineStyle(7)
        )
    elif (mass == 120):
        sigpdf.plotOn(
            xframe, ROOT.RooFit.Range("NormRange"),
            ROOT.RooFit.Normalization(intp_params["norm"], ROOT.RooAbsReal.NumEvent), 
            ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#0F52BA")), 
            ROOT.RooFit.LineWidth(4), ROOT.RooFit.Name("120")
        )
    elif (mass == 125):
        sigpdf.plotOn(
            xframe, ROOT.RooFit.Range("NormRange"),
            ROOT.RooFit.Normalization(intp_params["norm"], ROOT.RooAbsReal.NumEvent), 
            ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#1C7747")), 
            ROOT.RooFit.LineWidth(4), ROOT.RooFit.Name("125")
        )
    elif (mass == 130):
        sigpdf.plotOn(
            xframe, ROOT.RooFit.Range("NormRange"),
            ROOT.RooFit.Normalization(intp_params["norm"], ROOT.RooAbsReal.NumEvent), 
            ROOT.RooFit.LineColor(ROOT.TColor.GetColor("#E23E57")), 
            ROOT.RooFit.LineWidth(4), ROOT.RooFit.Name("130")
        )

    # build the output work space and save the yield and pdf
    outws = "./workspace/Interpolation/{}".format(year)
    if not os.path.exists(outws):
        os.makedirs(outws)
    fws = ROOT.TFile("{}/workspace_HDalitz_sigMC_{}_{}_{}.root".format(outws, int(mass), cat, proc), "RECREATE")
    fws.cd()
    print("[INFO] Saving the ws file in: {}/workspace_HDalitz_sigMC_{}_{}_{}.root".format(outws, int(mass), cat, proc))
    ws = ROOT.RooWorkspace("w")
    ExpYield = ROOT.RooRealVar("ExpYield", "ExpYield", intp_params["norm"], "GeV")
    ExpYield.setConstant(True)
    getattr(ws, "import")(ExpYield)
    getattr(ws, "import")(sigpdf)

    # set the parameters of pdf to be constant
    ws.defineSet("SigPdfParams_{}".format(cat), ROOT.RooArgSet(
        ws.var("mean_{}".format(cat)),
        ws.var("sigma_CB_{}".format(cat)),
        ws.var("alpha1_{}".format(cat)),
        ws.var("alpha2_{}".format(cat)),
        ws.var("n1_{}".format(cat)),
        ws.var("n2_{}".format(cat)),
        ws.var("sigma_Gauss_{}".format(cat)),
        ws.var("frac_DCB_{}".format(cat)) 
    ))
    for _var in rooiter(ws.set("SigPdfParams_{}".format(cat))):
        _var.setConstant(True) 

    # for the systematic uncertainties
    ws.factory("CMS_heeg_scale_{}_{}_{}[1]".format(proc, year, cat)) # set the initial value to be 1
    ws.factory("CMS_heeg_res_{}_{}_{}[1]".format(proc, year, cat))
    ws.factory("prod::newMean_{}_{}_{}(mean_{}, CMS_heeg_scale_{}_{}_{})".format(proc, year, cat, cat, proc, year, cat))
    ws.factory("prod::newSigma_{}_{}_{}(sigma_CB_{}, CMS_heeg_scale_{}_{}_{})".format(proc, year, cat, cat, proc, year, cat))
    ws.factory("EDIT::newSigPdf_{}_{}_{}(SigPdf_{}, mean_{} = newMean_{}_{}_{}, sigma_CB_{} = newSigma_{}_{}_{})".format(proc, year, cat, cat, cat, proc, year, cat, cat, proc, year, cat))
    
    ws.Write()
    fws.Close()

# set up the canvas to draw
xframe.SetTitle("")

xframe.GetXaxis().SetTickSize(0.03)
xframe.GetXaxis().SetTitleSize(0.04)
xframe.GetXaxis().SetLabelSize(0.04)
xframe.GetXaxis().SetLabelOffset(0.02)
xframe.GetXaxis().SetTitleOffset(1.4)
xframe.GetXaxis().SetTitle("M_{ee#gamma} [GeV]")

xframe.GetYaxis().SetTitle("Signal shape")
xframe.GetYaxis().SetNdivisions(510)
xframe.GetYaxis().SetTickSize(0.03)
xframe.GetYaxis().SetTitleSize(0.04)
xframe.GetYaxis().SetLabelSize(0.04)
xframe.GetYaxis().SetTitleOffset(1.8)

xframe.SetMaximum(xframe.GetMaximum() * 1.4)
xframe.SetMinimum(0.00001)

c = ROOT.TCanvas("c", "", 900, 900)
c.cd()
c.SetRightMargin(0.05)
c.SetTopMargin(0.07)
c.SetLeftMargin(0.14)
c.SetBottomMargin(0.12)
xframe.Draw()

catProc = ROOT.TLatex()
catProc.SetTextFont(42)
catProc.SetNDC()
catProc.SetTextSize(0.037)
catProc.DrawLatex(0.53, 0.86, "%s, %s" %(proc, cat))

leg1 = ROOT.TLegend(0.53, 0.7, 0.83, 0.84)
leg1.SetTextFont(42)
leg1.SetTextSize(0.035)
leg1.SetFillColor(0)
leg1.SetLineColor(0)
leg1.AddEntry(xframe.findObject("120"), "PDF-120 GeV ", "l")
leg1.AddEntry(xframe.findObject("125"), "PDF-125 GeV ", "l")
leg1.AddEntry(xframe.findObject("130"), "PDF-130 GeV ", "l")
leg1.Draw()

outdir = "./plots/Interpolation/{}".format(year)
if not os.path.exists(outdir):
    os.makedirs(outdir)

CMS_lumi(c, 4, 11, "", year, True, "Work-in-progress", "H #rightarrow #gamma*#gamma #rightarrow ee#gamma", "")

c.SaveAs("{}/intp_{}_{}.pdf".format(outdir, cat, proc))
c.Close()

print("")
