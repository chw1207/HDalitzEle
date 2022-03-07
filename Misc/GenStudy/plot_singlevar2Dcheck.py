#!/usr/bin/env python

import argparse
import ROOT
import os
import sys
import json

# import TdrStyle
import CMS_lumi
from plotUtil import gethist 

Lumi_2017 = 41.5
Lumi_2016 = 35.9
Lumi_2018 = 59.7
Lumi_FullRun2 = 35.9 + 41.5 + 59.7
refXS = 48.58 + 3.782 + 1.373 + 0.8839
BR_eeg = 8.07E-5
BR_mmg = 3.83E-5

ROOT.gROOT.SetBatch()

# The only way to make the tdrstyle work...
# Reference: https://jiafulow.github.io/blog/2017/08/22/cms-plot-styles/
ROOT.gROOT.LoadMacro(
    "/afs/cern.ch/work/h/hajheng/private/HDalitz/interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")


parser = argparse.ArgumentParser(
    description='Plot single variable for checking')
parser.add_argument('-v', '--verb', dest='verb', action='store_true', default=False,
                    help='Do more verbosity printout.')
parser.add_argument('-d', '--plotdir', dest="plotdir", type=str, default='./plots/SingleVar_check',
                    help='Directory to store plots')
opt = parser.parse_args()

if not opt.verb:
    ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = 1001;')

os.makedirs(opt.plotdir, exist_ok=True)


def Draw_2DHistSimp(hist, drawopt, logz, XaxisName, YaxisName, outname):
    ROOT.gStyle.SetPalette(ROOT.kBlueGreenYellow)
    TickSize = 0.03
    AxisTitleSize = 0.05
    AxisLabelSize = 0.04

    hist.Scale(1./hist.Integral(-1,-1,-1,-1))

    c = ROOT.TCanvas('c', '', 700, 600)
    c.cd()
    if logz:
        c.SetLogz()
    ROOT.gPad.SetRightMargin(0.13)
    ROOT.gPad.SetTopMargin(0.05)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)

    hist.GetXaxis().SetTitle(XaxisName)
    hist.GetYaxis().SetTitle(YaxisName)
    hist.GetXaxis().SetTickSize(TickSize)
    hist.GetYaxis().SetTickSize(TickSize)
    hist.GetXaxis().SetTitleSize(AxisTitleSize)
    hist.GetYaxis().SetTitleSize(AxisTitleSize)
    hist.GetXaxis().SetLabelSize(AxisLabelSize)
    hist.GetYaxis().SetLabelSize(AxisLabelSize)
    hist.GetXaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleOffset(1.3)
    hist.GetZaxis().SetLabelSize(AxisLabelSize)
    hist.SetContour(2000)
    hist.Draw(drawopt)
    c.RedrawAxis()

    c.SaveAs('{}/{}.png'.format(opt.plotdir, outname))
    c.SaveAs('{}/{}.pdf'.format(opt.plotdir, outname))


if __name__ == "__main__":
    if opt.verb:
        print("This is the __main__ part")

    with open('singlevar2Dcheck_param.json', 'r') as fp:
        myTags = json.load(fp)

        for tagName, par in myTags.items():
            print(tagName, par)

            hM = gethist('ele', tagName)
            Draw_2DHistSimp(hM, 'COLZ', True, par['xtitle'], par['ytitle'], par['plotname'])