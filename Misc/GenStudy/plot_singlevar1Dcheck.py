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


def Draw_1DHistSimp(hist, leg, drawopt, drawxoverflow, drawxunderflow, logy, logy_axisscale, legsty, XaxisName, yaxisunit, outname):
    TickSize = 0.03
    AxisTitleSize = 0.05
    AxisLabelSize = 0.05

    ROOT.gStyle.SetOptStat(0)

    hist.StatOverflows()
    hist.Scale(1./hist.Integral(-1,-1))
    ymaxval = hist.GetMaximum()
    yminval = hist.GetMinimum(0)
    binwidth = hist.GetXaxis().GetBinWidth(1)

    yaxtitletext = ''
    if yaxisunit == '':
        yaxtitletext = 'Normalized events / {:g}'.format(binwidth)
    else:
        yaxtitletext = 'Normalized events / {:g} {}'.format(binwidth, yaxisunit)

    c = ROOT.TCanvas('c', '', 700, 600)
    c.cd()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    if logy:
        c.SetLogy()
        hist.GetYaxis().SetRangeUser(yminval*0.95, ymaxval*logy_axisscale)
    else:
        hist.GetYaxis().SetRangeUser(0, ymaxval*1.15)

    if drawxoverflow:
        hist.GetXaxis().SetRange(1, hist.GetNbinsX() + 1)
    if drawxunderflow:
        hist.GetXaxis().SetRange(0, hist.GetNbinsX())
    if drawxoverflow and drawxunderflow:
        hist.GetXaxis().SetRange(0, hist.GetNbinsX() + 1)
    
    hist.SetLineWidth(3)
    hist.SetLineColor(ROOT.kBlack)
    hist.GetXaxis().SetTitle(XaxisName)
    hist.GetYaxis().SetTitle(yaxtitletext)
    hist.GetXaxis().SetTickSize(TickSize)
    hist.GetXaxis().SetTitleSize(AxisTitleSize)
    hist.GetXaxis().SetLabelSize(AxisLabelSize)
    hist.GetYaxis().SetTickSize(TickSize)
    hist.GetYaxis().SetTitleSize(AxisTitleSize)
    hist.GetYaxis().SetLabelSize(AxisLabelSize)
    hist.GetXaxis().SetTitleOffset(1.3)
    hist.GetYaxis().SetTitleOffset(1.4)
    hist.Draw(drawopt)
    # CMS_lumi.CMS_lumi(c, 4, 11, '41.5 fb^{-1}', 2017, True, 'Simulation', '', '') # CMS Simulation inside the frame
    # CMS_lumi.CMS_lumi(c, 5, 0, '', 2017, True, 'Simulation','', '')  # CMS Simulation outside the frame
    c.Modified()
    c.RedrawAxis()

    # l = ROOT.TLegend(0.7, 0.83, 0.85, 0.87)
    # l.SetTextSize(0.04)
    # l.AddEntry(hist, leg, legsty)
    # l.SetFillColor(0)  # Set the background to be white
    # l.SetLineColor(0)
    # l.SetFillStyle(0)
    # l.Draw("same")

    c.SaveAs('{}/{}.png'.format(opt.plotdir, outname))
    c.SaveAs('{}/{}.pdf'.format(opt.plotdir, outname))


if __name__ == "__main__":
    if opt.verb:
        print("This is the __main__ part")

    with open('singlevar1Dcheck_param.json', 'r') as fp:
        myTags = json.load(fp)

        for tagName, par in myTags.items():
            print(tagName, par)

            hM = gethist('ele', tagName)

            setlogy = False
            if par["logy"] == 'true':
                setlogy = True

            drawoverflow = False
            if par["xoverflow"] == 'true':
                drawoverflow = True

            drawunderflow = False
            if par["xunderflow"] == 'true':
                drawunderflow = True

            Draw_1DHistSimp(hM, par["leg"], 'hist', drawoverflow, drawunderflow, setlogy, par["logy_axisscale"], 'l', par['xtitle'], par['yunit'], par['plotname'])