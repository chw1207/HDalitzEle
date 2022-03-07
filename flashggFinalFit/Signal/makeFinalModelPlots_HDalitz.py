import ROOT
import sys
import os
import json
import numpy as np
from argparse import ArgumentParser
from collections import OrderedDict as od
from tools.CMS_lumi import CMS_lumi
from tools.sigmaEff import sigmaEff


def get_parser():
    parser = ArgumentParser(description = "Script to built the final signal model")
    parser.add_argument("-cat",  "--category",  help = "category of the signal", type = str)
    parser.add_argument("-proc", "--process",   help = "process of the signal",  type = str)
    parser.add_argument("-dBs",  "--dataBins",  default = (170 - 105), help = "number of bins for data", type = int)
    parser.add_argument("-m",    "--masspoint", default = 125, help = "mass point", type = int)

    return parser


def getWs(fwname, masspoint, cat, process):
    fw = ROOT.TFile.Open(fwname, "READ")
    if not fw:
        sys.exit(1)
    inWS = ROOT.RooWorkspace()
    fw.GetObject("w", inWS)
    fw.Close()

    return inWS


def plotSignalModel(hists, eff_sigma_low, eff_sigma_up, eff_sigma, outName, cat, process, offset = 0.03):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    colorMap = {"2016":38, "2017":30, "2018":46}

    canv = ROOT.TCanvas("canv", "", 850, 800)
    canv.cd()
    canv.SetRightMargin(0.05)
    canv.SetTopMargin(0.07)
    canv.SetLeftMargin(0.13)
    canv.SetBottomMargin(0.12)

    h_axes = hists["data"].Clone()
    h_axes.Reset()
    h_axes.SetMaximum(hists["data"].GetMaximum() * 1.2)
    h_axes.SetMinimum(0.)
    h_axes.GetXaxis().SetRangeUser(105, 140)
    h_axes.SetTitle("")
    h_axes.GetXaxis().SetTickSize(0.03)
    h_axes.GetXaxis().SetTitleSize(0.04)
    h_axes.GetXaxis().SetLabelSize(0.04)
    h_axes.GetXaxis().SetLabelOffset(0.02)
    h_axes.GetXaxis().SetTitleOffset(1.4)
    h_axes.GetXaxis().SetTitle("M_{ee#gamma} [GeV]")
    # h_axes.GetYaxis().SetTitle("Signal shape / ({} GeV)".format((170 - 105)/Nbins))
    h_axes.GetYaxis().SetNdivisions(510)
    h_axes.GetYaxis().SetTickSize(0.03)
    h_axes.GetYaxis().SetTitleSize(0.04)
    h_axes.GetYaxis().SetLabelSize(0.04)
    h_axes.GetYaxis().SetTitleOffset(1.7)
    h_axes.Draw()

    # Extract effSigma
    effSigma = eff_sigma["all_data"]
    effSigma_low, effSigma_high = eff_sigma_low["all_data"], eff_sigma_up["all_data"]
    h_effSigma = hists["pdf"].Clone()
    h_effSigma.GetXaxis().SetRangeUser(effSigma_low, effSigma_high)
    print("eff_sigma = {:.2f} [GeV]".format(effSigma))

    # Set style effSigma
    h_effSigma.SetLineColor(15)
    h_effSigma.SetLineWidth(3)
    h_effSigma.SetFillStyle(1001)
    h_effSigma.SetFillColor(19)
    h_effSigma.Draw("Same Hist F")
    vline_effSigma_low = ROOT.TLine(effSigma_low, 0, effSigma_low, hists["pdf"].GetBinContent(hists["pdf"].FindBin(effSigma_low)))
    vline_effSigma_high = ROOT.TLine(effSigma_high, 0, effSigma_high, hists["pdf"].GetBinContent(hists["pdf"].FindBin(effSigma_high)))
    vline_effSigma_low.SetLineColor(15)
    vline_effSigma_high.SetLineColor(15)
    vline_effSigma_low.SetLineWidth(3)
    vline_effSigma_high.SetLineWidth(3)
    vline_effSigma_low.Draw("Same")
    vline_effSigma_high.Draw("Same")

    # Set style pdf
    hists["pdf"].SetLineColor(4)
    hists["pdf"].SetLineWidth(3)
    hists["pdf"].Draw("Same Hist C")
   
    for year in ["2016", "2017", "2018"]:
        hists["pdf_{}".format(year)].SetLineColor(colorMap[year])  
        hists["pdf_{}".format(year)].SetLineStyle(2)
        hists["pdf_{}".format(year)].SetLineWidth(3)
        hists["pdf_{}".format(year)].Draw("Same Hist C")

    # Set style: data
    hists["data"].SetMarkerStyle(25)
    hists["data"].SetMarkerColor(1)
    hists["data"].SetMarkerSize(1.5)
    hists["data"].SetLineColor(1)
    hists["data"].SetLineWidth(2)
    hists["data"].Draw("Same ep X0")

    # legend
    leg0 = ROOT.TLegend(0.14+offset, 0.68, 0.5+offset, 0.82)
    leg0.SetFillStyle(0)
    leg0.SetLineColorAlpha(0,0)
    leg0.SetTextSize(0.033)
    leg0.AddEntry(hists["data"], "Simulation", "ep")
    leg0.AddEntry(hists["pdf"], "#splitline{Parametric}{model}", "l")
    leg0.Draw("Same")

    leg1 = ROOT.TLegend(0.16+offset, 0.52, 0.4+offset, 0.68)
    leg1.SetFillStyle(0)
    leg1.SetLineColorAlpha(0,0)
    leg1.SetTextSize(0.029)
    for year in [2016, 2017, 2018]: 
        leg1.AddEntry(hists["pdf_{}".format(year)], "%i: #sigma_{eff} = %1.2f GeV" %(year, eff_sigma["{}".format(year)]), "l")
    leg1.Draw("Same")

    leg2 = ROOT.TLegend(0.14+offset, 0.42, 0.4+offset, 0.52)
    leg2.SetFillStyle(0)
    leg2.SetLineColorAlpha(0,0)
    leg2.SetTextSize(0.033)
    leg2.AddEntry(h_effSigma, "#sigma_{eff} = %1.2f GeV" %effSigma, "fl")
    leg2.Draw("Same")

    mode = ROOT.TLatex()
    mode.SetTextFont(42)
    mode.SetNDC()
    mode.SetTextSize(0.04)
    mode.DrawLatex(0.16+offset, 0.86, "H #rightarrow #gamma*#gamma #rightarrow ee#gamma")

    catProc = ROOT.TLatex()
    catProc.SetTextFont(42)
    catProc.SetNDC()
    catProc.SetTextSize(0.035)
    catProc.DrawLatex(0.55+offset, 0.86, "{}, {}".format(process, cat))

    CMS_lumi(canv, 5, 0, "", 2017, True, "Work-in-progress", "", "")
    canv.Update()
    canv.RedrawAxis()

    canv.SaveAs(outName)
    canv.Close()
    del canv

    # Write effSigma to file
    outDir = "./eff_sigma"
    if not os.path.exists(outDir):
        os.makedirs(outDir)
    es = {}
    es["combined"] = effSigma
    for year in [2016, 2017, 2018]: 
        es[year] = eff_sigma["{}".format(year)]
    print ("Effective sigma json file: {}/effSigma_{}_{}_{}.json".format(outDir, masspoint, cat, process))
    with open("{}/effSigma_{}_{}_{}.json".format(outDir, masspoint, cat, process), "w") as jf: 
        json.dump(es, jf)

#############################
##                         ##
##    Main script part     ##
##                         ##
#############################
def main():
    # container
    hists, norm, hpdfs = od(), od(), od()
    data = od()

    # build data histogram
    CMS_higgs_mass = ROOT.RooRealVar("CMS_higgs_mass", "CMS_higgs_mass", 105, 170, "GeV")
    CMS_higgs_mass.setRange("NormRange", 105, 170)
    CMS_higgs_mass.setMin(105) #! do not remove this line
    CMS_higgs_mass.setMax(170) #! do not remove this line
    hists["data"] = CMS_higgs_mass.createHistogram("h_data", ROOT.RooFit.Binning(dataBins)) # create a empty histogram for dataset 

    # datasets and eff_sigma
    ScaleNumber = 400 # use to create smooth histogram
    eff_sigma_down, eff_sigma_up, eff_sigma = od(), od(), od()
    vall = []
    for year in [2016, 2017, 2018]:
        inWS = getWs("./workspace/{}/workspace_HDalitz_sigMC_{}_{}_{}.root".format(year, masspoint, cat, proc), masspoint, cat, proc)
        
        data[year] = inWS.data("set")
        data[year].fillHistogram(hists["data"], ROOT.RooArgList(CMS_higgs_mass)) # fill the histogram for dataset
        v = []
        for i in range(data[year].numEntries()):
            v.append(data[year].get(i).getRealValue(CMS_higgs_mass.GetName()))
            vall.append(data[year].get(i).getRealValue(CMS_higgs_mass.GetName()))
        eff_sigma_down["{}".format(year)], eff_sigma_up["{}".format(year)], eff_sigma["{}".format(year)] = sigmaEff(v)

        # norm[year] = data[year].sumEntries() # used to sale the pdf
        pdf = inWS.pdf("SigPdf")
        hpdfs[year] = pdf.createHistogram("h_pdf_%s" %year, CMS_higgs_mass, ROOT.RooFit.Binning(dataBins * ScaleNumber))

    eff_sigma_down["all_data"], eff_sigma_up["all_data"], eff_sigma["all_data"] = sigmaEff(vall)

    # Sum pdf histograms
    for k, p in hpdfs.iteritems():
        if "pdf" not in hists: 
            hists["pdf"] = p.Clone("h_pdf")
            hists["pdf"].Reset()
        hists["pdf"] += p
    hists["pdf"].Scale(hists["data"].Integral() * ScaleNumber / hists["pdf"].Integral())

    # Per-year pdf histograms
    for year in [2016, 2017, 2018]:
        if "pdf_{}".format(year) not in hists:
            hists["pdf_{}".format(year)] = hists["pdf"].Clone()
            hists["pdf_{}".format(year)].Reset()
        for i, p in hpdfs.iteritems():
            if (year == i):
                hists["pdf_{}".format(year)] += p 

        hists["pdf_{}".format(year)].Scale(data[year].sumEntries() * ScaleNumber / hists["pdf_{}".format(year)].Integral())

    # save the final signal model for all three years
    out = "./plots/final"
    print("[INFO] Saving the ws file in: {}/FinalModel_HDalitz_sigMC_{}_{}_{}.root".format(out, int(masspoint), cat, proc))
    f = ROOT.TFile("{}/FinalModel_HDalitz_sigMC_{}_{}_{}.root".format(out, int(masspoint), cat, proc), "RECREATE")
    f.cd()
    hists["pdf"].Write()
    f.Close()

    # draw the final model
    outPlot = "./plots/final"
    if not os.path.exists(outPlot):
        os.makedirs(outPlot)
    outName = "{}/FinalModel_{}_{}_{}.pdf".format(outPlot, masspoint, cat, proc)
    plotSignalModel(hists, eff_sigma_down, eff_sigma_up, eff_sigma, outName, cat, proc)

if __name__ == "__main__" :

    # get the arguments
    parser = get_parser()
    args = parser.parse_args()
    proc, cat = args.process, args.category
    dataBins = args.dataBins
    masspoint = args.masspoint

    state = main()
    sys.exit(state)

