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
refXS = 48.58 + 3.782 + 1.373 + 0.8839
BR_eeg = 8.10E-5
BR_mmg = 3.90E-5

ROOT.gROOT.SetBatch()

# The only way to make the tdrstyle work...
# Reference: https://jiafulow.github.io/blog/2017/08/22/cms-plot-styles/
ROOT.gROOT.LoadMacro(
    "/afs/cern.ch/work/h/hajheng/private/HDalitz/interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")


parser = argparse.ArgumentParser(
    description='Plot reco-level variables of electron')
parser.add_argument('-v', '--verb', dest='verb', action='store_true', default=False,
                    help='Do more verbosity printout.')
parser.add_argument('-d', '--plotdir', dest="plotdir", type=str, default='./plots/EleRecoVar_1D',
                    help='Directory to store plots')
opt = parser.parse_args()

if not opt.verb:
    ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = 1001;')

os.makedirs(opt.plotdir, exist_ok=True)


def Draw_1DHistComp(listhist, listleg, drawopt, drawxoverflow, drawxunderflow, logy, logy_axisscale, legsty, XaxisName, yaxisunit, outname):
    TickSize = 0.03
    AxisTitleSize = 0.05
    AxisLabelSize = 0.05

    ROOT.gStyle.SetOptStat(0)

    for hist in listhist:
        hist.StatOverflows()
        hist.Scale(1./hist.Integral(-1, -1))
        print('{} integral = {}'.format(hist.GetName(), hist.Integral(-1, -1)))

    ymaxval = 0.
    yminval = 1E10
    for hist in listhist:
        if hist.GetMaximum() > ymaxval:
            ymaxval = hist.GetMaximum()

        if hist.GetMinimum(0) < yminval:
            yminval = hist.GetMinimum(0)

    binwidth = listhist[0].GetXaxis().GetBinWidth(1)

    yaxtitletext = ''
    if yaxisunit == '':
        yaxtitletext = 'Normalized events / {:g}'.format(binwidth)
    else:
        yaxtitletext = 'Normalized events / {:g} {}'.format(
            binwidth, yaxisunit)

    c = ROOT.TCanvas('c', '', 700, 600)
    c.cd()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    if logy:
        c.SetLogy()
        listhist[0].GetYaxis().SetRangeUser(
            yminval*0.95, ymaxval*logy_axisscale)
    else:
        listhist[0].GetYaxis().SetRangeUser(0, ymaxval*3)

    if drawxoverflow:
        listhist[0].GetXaxis().SetRange(1, listhist[0].GetNbinsX() + 1)
    if drawxunderflow:
        listhist[0].GetXaxis().SetRange(0, listhist[0].GetNbinsX())
    if drawxoverflow and drawxunderflow:
        listhist[0].GetXaxis().SetRange(0, listhist[0].GetNbinsX() + 1)
    listhist[0].GetXaxis().SetTitle(XaxisName)
    listhist[0].GetYaxis().SetTitle(yaxtitletext)
    listhist[0].GetXaxis().SetTickSize(TickSize)
    listhist[0].GetXaxis().SetTitleSize(AxisTitleSize)
    listhist[0].GetXaxis().SetLabelSize(AxisLabelSize)
    listhist[0].GetYaxis().SetTickSize(TickSize)
    listhist[0].GetYaxis().SetTitleSize(AxisTitleSize)
    listhist[0].GetYaxis().SetLabelSize(AxisLabelSize)
    listhist[0].GetXaxis().SetTitleOffset(1.3)
    listhist[0].GetYaxis().SetTitleOffset(1.4)
    listhist[0].Draw(drawopt)
    for i in range(1, len(listhist)):
        listhist[i].GetXaxis().SetRange(1, listhist[i].GetNbinsX() + 1)
        listhist[i].Draw('{} same'.format(drawopt))

    # CMS_lumi.CMS_lumi(c, 4, 11, '41.5 fb^{-1}', 2017, True, 'Simulation', '', '') # CMS Simulation inside the frame
    CMS_lumi.CMS_lumi(c, 5, 0, '', 2017, True, 'Simulation',
                      '', '')  # CMS Simulation outside the frame
    c.Modified()
    c.RedrawAxis()

    l = ROOT.TLegend(0.54, 0.73, 0.8, 0.89)
    l.SetTextSize(0.04)
    for ileg, leg in enumerate(listleg):
        l.AddEntry(listhist[ileg], leg, legsty)
    l.SetFillColor(0) #Set the background to be white
    l.SetLineColor(0)
    l.SetFillStyle(0)
    l.Draw("same")

    c.SaveAs('{}/{}.png'.format(opt.plotdir, outname))
    c.SaveAs('{}/{}.pdf'.format(opt.plotdir, outname))


if __name__ == "__main__":
    if opt.verb:
        print("This is the __main__ part")

    list_histname = ['Resolved', 'Merged2Gsf', 'Merged1MissingGsf']
    list_leg = ['Resolved', 'Merged-2Gsf', 'Merged-1MissingGsf']
    list_color = ['#4D4C4C', '#8298CD', '#C06463']

    with open('varParams.json', 'r') as fp:
        myTags = json.load(fp)

        for tagName, par in myTags.items():
            print(tagName, par)
            list_hist = []
            for ih, hname in enumerate(list_histname):
                tmphist = gethist('ele', '{}_{}'.format(tagName, hname))
                tmphist.SetDirectory(0)
                tmphist.SetLineWidth(3)
                tmphist.SetLineColor(ROOT.TColor.GetColor(list_color[ih]))
                list_hist.append(tmphist)

            setlogy = False
            if par["logy"] == 'true':
                setlogy = True

            drawoverflow = True
            if par["xoverflow"] == 'false':
                drawoverflow = False

            drawunderflow = False
            if par["xunderflow"] == 'true':
                drawunderflow = True

            Draw_1DHistComp(list_hist, list_leg, 'HIST', drawoverflow, drawunderflow, setlogy,
                            par["logy_axisscale"], 'l', par["xAxisLabel"], par["yAxisUnit"], par["plotname"])
