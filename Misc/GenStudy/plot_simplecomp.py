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
BR_eeg = 8.10E-5
BR_mmg = 3.90E-5

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
parser.add_argument('-d', '--plotdir', dest="plotdir", type=str, default='./plots/Simple_comparison',
                    help='Directory to store plots')
opt = parser.parse_args()

if not opt.verb:
    ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = 1001;')

os.makedirs(opt.plotdir, exist_ok=True)


def Draw_1DHistComp(listhist, listleg, listcolor, drawopt, drawxoverflow, drawxunderflow, logy, logy_axisscale, legsty, XaxisName, yaxisunit, outname):
    TickSize = 0.03
    AxisTitleSize = 0.05
    AxisLabelSize = 0.05

    ROOT.gStyle.SetOptStat(0)

    ymaxval = 0.
    yminval = 1E10
    for hist in listhist:
        hist.StatOverflows()
        hist.Scale(1./hist.Integral(-1,-1))
        if hist.GetMaximum() > ymaxval:
            ymaxval = hist.GetMaximum()

        if hist.GetMinimum(0) < yminval:
            yminval = hist.GetMinimum(0)
    
    binwidth = listhist[0].GetXaxis().GetBinWidth(1)

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
        listhist[0].GetYaxis().SetRangeUser(yminval*0.95, ymaxval*logy_axisscale)
    else:
        listhist[0].GetYaxis().SetRangeUser(0, ymaxval*1.15)

    if drawxoverflow:
        listhist[0].GetXaxis().SetRange(1, listhist[0].GetNbinsX() + 1)
    if drawxunderflow:
        listhist[0].GetXaxis().SetRange(0, listhist[0].GetNbinsX())
    if drawxoverflow and drawxunderflow:
        listhist[0].GetXaxis().SetRange(0, listhist[0].GetNbinsX() + 1)
    
    listhist[0].SetLineWidth(3)
    listhist[0].SetLineColor(ROOT.TColor.GetColor(listcolor[0]))
    listhist[0].GetXaxis().SetTitle(XaxisName)
    listhist[0].GetYaxis().SetTitle(yaxtitletext)
    listhist[0].GetXaxis().SetTickSize(TickSize)
    listhist[0].GetXaxis().SetTitleSize(AxisTitleSize)
    listhist[0].GetXaxis().SetLabelSize(AxisLabelSize)
    listhist[0].GetYaxis().SetTickSize(TickSize)
    listhist[0].GetYaxis().SetTitleSize(AxisTitleSize)
    listhist[0].GetYaxis().SetLabelSize(AxisLabelSize)
    listhist[0].GetXaxis().SetTitleOffset(1.3)
    listhist[0].GetYaxis().SetTitleOffset(1.5)
    listhist[0].Draw(drawopt)
    for i in range(1, len(listhist)):
        listhist[i].GetXaxis().SetRange(1, listhist[i].GetNbinsX() + 1)
        listhist[i].SetLineWidth(3)
        listhist[i].SetLineColor(ROOT.TColor.GetColor(listcolor[i]))
        listhist[i].Draw('{} same'.format(drawopt))
    # CMS_lumi.CMS_lumi(c, 4, 11, '41.5 fb^{-1}', 2017, True, 'Simulation', '', '') # CMS Simulation inside the frame
    # CMS_lumi.CMS_lumi(c, 5, 0, '', 2017, True, 'Simulation','', '')  # CMS Simulation outside the frame
    c.Modified()
    c.RedrawAxis()

    legpos_y1 = 0.87
    legpos_dy = 0.04
    l = ROOT.TLegend(0.67, legpos_y1 - legpos_dy * len(listhist), 0.85, legpos_y1)
    l.SetTextSize(0.04)
    for ihist, hist in enumerate(listhist):
        l.AddEntry(hist, listleg[ihist], legsty)
    l.SetFillColor(0)  # Set the background to be white
    l.SetLineColor(0)
    l.SetFillStyle(0)
    l.Draw("same")

    c.SaveAs('{}/{}.png'.format(opt.plotdir, outname))
    c.SaveAs('{}/{}.pdf'.format(opt.plotdir, outname))


if __name__ == "__main__":
    if opt.verb:
        print("This is the __main__ part")

    with open('SimpComp_param.json', 'r') as fp:
        myTags = json.load(fp)

        for tagName, par in myTags.items():
            # print(tagName, par)
            listleg = []
            listcolor = []
            listhist = []
            for hist in par['hist']:
                print (par['hist']['{}'.format(hist)])
                listhist.append(gethist('ele', par['hist']['{}'.format(hist)]))
            for ent in par['leg']:
                listleg.append(par['leg']['{}'.format(ent)])
            for col in par['Color']:
                listcolor.append(par['Color']['{}'.format(col)])

            setlogy = False
            if par["logy"] == 'true':
                setlogy = True

            drawoverflow = False
            if par["xoverflow"] == 'true':
                drawoverflow = True

            drawunderflow = False
            if par["xunderflow"] == 'true':
                drawunderflow = True

            # Draw_1DHistComp(listhist, listleg, listcolor, drawopt, drawxoverflow, drawxunderflow, logy, logy_axisscale, legsty, XaxisName, yaxisunit, outname):
            Draw_1DHistComp(listhist, listleg, listcolor, "HIST", drawoverflow, drawunderflow, setlogy, par["logy_axisscale"], 'l', par['xtitle'], par['yunit'], par['plotname'])