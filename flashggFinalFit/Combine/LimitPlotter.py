import os, sys
import ROOT
import json
import numpy as np
import TdrStyle as tdr
from argparse import ArgumentParser
from CMS_lumi import CMS_lumi 
from array import array
from commonTools_HDalitz import *

def getValuesFromFile(fname):
    f = ROOT.TFile.Open(fname)
    if not f:
        sys.exit(1)
    tree = f.Get("limit")

    res = []
    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        res.append(getattr(tree, "limit"))
    f.Close()
    
    return res

def limitVsM(masses):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    exp2SigLow, exp1SigLow, expMean, exp1SigHi, exp2SigHi = [], [], [], [], []
    xAxis, xErr = [], []

    for i, m in enumerate(masses):
        l = getValuesFromFile("../Datacard/electron/higgsCombine_combine_{}.AsymptoticLimits.mH{}.root".format(m, m))
        exp2SigLow.append(l[0])
        exp1SigLow.append(l[1])
        expMean.append(l[2])
        exp1SigHi.append(l[3])
        exp2SigHi.append(l[4])

        xAxis.append(m)
        xErr.append(0.5)

    nPoints = len(xAxis)
    zeros_Array = array("f", [0 for i in range(nPoints)])
    xAxis_Array = array("f", xAxis)
    xErr_Array = zeros_Array

    exp_Array = array("f", expMean)
    exp2SigLowErr_Array = array("f", [a-b for a, b in zip(expMean, exp2SigLow)])
    exp1SigLowErr_Array = array("f", [a-b for a, b in zip(expMean, exp1SigLow)])
    exp1SigHiErr_Array = array("f", [b-a for a, b in zip(expMean, exp1SigHi)])
    exp2SigHiErr_Array = array("f", [b-a for a, b in zip(expMean, exp2SigHi)])

    # exp_max = np.max(exp2SigHiErr_Array)

    c1 = ROOT.TCanvas("c1", "c1", 700, 700)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetBottomMargin(0.15)
    ROOT.gPad.SetLeftMargin(0.15)
    c1.cd()
    mg = ROOT.TMultiGraph()
    mg.SetTitle("")

    
    # print exp_Array
    expected = ROOT.TGraphAsymmErrors(
        nPoints, array("f", xAxis), array("f", expMean)
    )
    oneSigma = ROOT.TGraphAsymmErrors(
        nPoints, xAxis_Array, exp_Array, xErr_Array, xErr_Array, exp1SigLowErr_Array, exp1SigHiErr_Array
    )
    twoSigma = ROOT.TGraphAsymmErrors(
        nPoints, xAxis_Array, exp_Array, xErr_Array, xErr_Array, exp2SigLowErr_Array, exp2SigHiErr_Array
    )

    oneSigma.SetFillColor(ROOT.kGreen+1)
    twoSigma.SetFillColor(ROOT.kOrange)

    expected.SetMarkerColor(ROOT.kBlack)
    expected.SetLineColor(ROOT.kBlack)
    expected.SetLineWidth(3)
    expected.SetLineStyle(7)

    mg.Add(twoSigma, "")
    mg.Add(oneSigma, "")
    mg.Add(expected, "")

    mg.Draw("AL3")
    mg.GetXaxis().SetTitle("M_{H} [GeV]")
    # mg.SetMinimum(0)
    # mg.SetMaximum(15)
    max_element = max(exp2SigHi)
    mg.GetXaxis().SetLimits(119.5, 130.5)
    mg.GetYaxis().SetTitle("95% CL limit on #sigma/#sigma_{SM}")
    mg.GetXaxis().SetTitleOffset(1.1)
    mg.GetYaxis().SetTitleOffset(1)
    # mg.GetYaxis().SetRangeUser(0, max_element * 2)
    mg.GetYaxis().SetRangeUser(0, 12)
    mg.Draw("zE2")
    mg.Draw("AL3")

    leg = ROOT.TLegend(0.6, 0.72, 0.9, 0.89)
    # leg.SetNColumns(2)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(expected, "Expected", "l")
    leg.AddEntry(oneSigma, "68% expected", "f")
    leg.AddEntry(twoSigma, "95% expected", "f")
    leg.Draw()

    line = ROOT.TLine(119.5, 1, 130.5, 1)
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(2)
    line.Draw()

    CMS_lumi(c1, 5, 0, "137.1 fb^{-1}", 0, True, "Work-in-progress", "", "")
    c1.Print("./LimitVSMass.pdf")


def limitSummary():
    BinOrder = [
        "combine",
        "Resolved",
        "Merged2Gsf_tagged",
        "Merged2Gsf_Untagged",
        "Merged1Gsf_tagged",
        "Merged1Gsf_Untagged",
        
    ]

    nameMap = {
        "combine": "Combined",
        "Merged2Gsf_tagged": "Merged-2Gsf tagged",
        "Merged2Gsf_Untagged": "Merged-2Gsf Untagged",
        "Merged1Gsf_tagged": "Merged-1Gsf tagged",
        "Merged1Gsf_Untagged": "Merged-1Gsf Untagged",
        "Resolved": "Resolved"
    }
    # BinOrder.append("combine")
    # BinOrder.reverse()

    exp2SigLow, exp1SigLow, expMean, exp1SigHi, exp2SigHi, expMeanSigInj = [], [], [], [], [], []
    yAxis = []

    for i, binname in enumerate(BinOrder):
        l = getValuesFromFile("../Datacard/electron/higgsCombine_{}_125.AsymptoticLimits.mH125.root".format(binname))
        exp2SigLow.append(l[0])
        exp1SigLow.append(l[1])
        expMean.append(l[2])
        exp1SigHi.append(l[3])
        exp2SigHi.append(l[4])

        yAxis.append(float(i*1.7))

    nPoints  = len(yAxis)
    zeros_Array = array("f", [0 for i in range(nPoints)])
    yAxis_Array = array("f", yAxis)
    yErr_Array = array("f", [0.15 for i in range(nPoints)])

    exp_Array = array("f", expMean)
    exp2SigLowErr_Array = array("f", [a-b for a,b in zip(expMean,exp2SigLow)])
    exp1SigLowErr_Array = array("f", [a-b for a,b in zip(expMean,exp1SigLow)])
    exp1SigHiErr_Array  = array("f", [b-a for a,b in zip(expMean,exp1SigHi)])
    exp2SigHiErr_Array  = array("f", [b-a for a,b in zip(expMean,exp2SigHi)])

    xMax = 500
    xMin = 0.5
    gPadRightMargin = 0.04
    gPadTopMargin = 0.06
    gPadBorromMargin = 0.12
    gPadLeftMargin = 0.32
    
    c1 = ROOT.TCanvas("c1", "c1", 850, 800)
    ROOT.gPad.SetRightMargin(gPadRightMargin)
    ROOT.gPad.SetTopMargin(gPadTopMargin)
    ROOT.gPad.SetBottomMargin(gPadBorromMargin)
    ROOT.gPad.SetLeftMargin(gPadLeftMargin)
    c1.cd()
    c1.SetLogx()
    mg = ROOT.TMultiGraph()
    mg.SetTitle("")

    expected = ROOT.TGraphAsymmErrors(nPoints,exp_Array,yAxis_Array,zeros_Array,zeros_Array,zeros_Array,zeros_Array)
    oneSigma = ROOT.TGraphAsymmErrors(nPoints,exp_Array,yAxis_Array,exp1SigLowErr_Array,exp1SigHiErr_Array,yErr_Array,yErr_Array)
    twoSigma = ROOT.TGraphAsymmErrors(nPoints,exp_Array,yAxis_Array,exp2SigLowErr_Array,exp2SigHiErr_Array,yErr_Array,yErr_Array)

    oneSigma.SetFillColor(ROOT.kGreen+1)
    twoSigma.SetFillColor(ROOT.kOrange)

    expected.SetMarkerStyle(20)
    expected.SetMarkerSize(1.5)
    expected.SetMarkerColor(ROOT.kBlack)
    expected.SetLineColor(ROOT.kBlack)
    # expected.SetLineWidth(5)
    # expected.SetLineStyle(2)

    mg.Add(twoSigma, "2")
    mg.Add(oneSigma, "2")
    mg.Add(expected, "e1p")

    mg.Draw("ap")
    mg.SetMinimum(-1)
    mg.SetMaximum(nPoints+7)
    mg.GetYaxis().SetTickLength(0)
    # mg.GetXaxis().SetLimits(1, 20) 
    mg.GetYaxis().SetLabelOffset(999)
    # mg.GetXaxis().SetLabelOffset()
    mg.GetXaxis().SetLabelFont(42)
    mg.GetXaxis().SetTitleSize(0.04)
    mg.GetXaxis().SetTitleOffset(1.3)
    mg.GetXaxis().SetLimits(xMin, xMax+50)  
    mg.GetXaxis().SetTitle("95% CL upper limit on #sigma/#sigma_{SM}")
    mg.GetXaxis().SetTitleOffset(1.5)

    l_1 = ROOT.TLine(1, -1, 1, 13)
    l_1.SetLineColor(ROOT.kRed+1)
    l_1.SetLineWidth(2)
    l_1.Draw()

    leg = ROOT.TLegend(0.65, 0.77, 0.92, 0.94*(1-gPadTopMargin))
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    leg.AddEntry(expected,"Expected", "p")
    leg.AddEntry(oneSigma,"68% expected ", "f")
    leg.AddEntry(twoSigma,"95% expected", "f")
    leg.Draw()

    lat = ROOT.TLatex()
    # Here we add the names of all categories  lat = TLatex()
    lat.SetTextColor(ROOT.kBlack)
    lat.SetTextSize(0.025)
    lat.SetTextAlign(32)
    for i, o in enumerate(BinOrder):
        lat.SetTextSize(0.025)
        lat.SetTextColor(ROOT.kBlack)
        lat.SetTextSize(0.030)
        if o != "combine":
            lat.SetTextFont(42)
            lat.DrawLatex(xMin-0.07, yAxis_Array[i], nameMap[o])
        else:
            lat.SetTextFont(72)
            lat.DrawLatex(xMin-0.07, yAxis_Array[i], "Combined")


    lat2 = ROOT.TLatex()
    lat.SetTextColor(ROOT.kBlack)
    lat.SetTextSize(0.03)
    lat.SetTextAlign(32)
    for i, exp in enumerate(exp_Array):
        # lat2.SetTextSize(0.025)
        lat2.SetTextColor(ROOT.kBlack)
        lat2.SetTextSize(0.04)
        lat2.SetTextFont(42)
        lat.DrawLatex(xMax*0.8, yAxis_Array[i], "{}".format(round(exp, 2)))

    CMS_lumi(c1, 5, 1, "137.1 fb^{-1}", 2017, "", "Work-in-progress", "H #rightarrow #gamma*#gamma #rightarrow ee#gamma", "")
    c1.Update()


    lat1 = ROOT.TLatex()
    lat1.SetTextColor(ROOT.kBlack)
    lat1.SetTextSize(0.030)
    lat1.SetTextFont(52)
    lat1.DrawLatex(1.5, 13.2, "Work-in-progress")
    lat1.Draw()


    c1.Print("./LimitSummary.pdf")












ROOT.gROOT.SetBatch()
tdr.setTDRStyle()
masses = [120+i for i in range(11)]

limitVsM(masses)
limitSummary()
