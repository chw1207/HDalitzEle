#!/usr/bin/env python

import argparse
from ROOT import gROOT, gStyle, gPad, TH1F, TCanvas, TLegend, TChain, TLorentzVector, TColor, THStack
import os
import sys
import json
from array import array
import numpy
import math

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

gROOT.SetBatch()

# The only way to make the tdrstyle work...
# Reference: https://jiafulow.github.io/blog/2017/08/22/cms-plot-styles/
gROOT.LoadMacro(
    '/afs/cern.ch/work/h/hajheng/private/HDalitz/interface/tdrstyle.C')
gROOT.ProcessLine('setTDRStyle()')


parser = argparse.ArgumentParser(
    description='Plot reco-level variables of electron')
parser.add_argument('-v', '--verb', dest='verb', action='store_true', default=False,
                    help='Do more verbosity printout.')
parser.add_argument('-d', '--plotdir', dest='plotdir', type=str, default='./plots/sigeff_HLT',
                    help='Directory to store plots')
opt = parser.parse_args()

if not opt.verb:
    gROOT.ProcessLine('gErrorIgnoreLevel = 1001')

os.makedirs(opt.plotdir, exist_ok=True)


def Drawhist(listhist, listleg, listlinesty, drawopt, drawxoverflow, drawxunderflow, logy, axisscale, legsty, XaxisName, yaxisunit, outname):
    lcolor = ['#1b262c', '#0f4c75', '#3282b8', '#900d0d', '#ce6262', '#6a2c70']
    TickSize = 0.03
    AxisTitleSize = 0.05
    AxisLabelSize = 0.05

    gStyle.SetOptStat(0)

    # for hist in listhist:
    #     hist.StatOverflows()
    #     hist.Scale(1./hist.Integral(-1, -1))
    #     print('{} integral = {}'.format(hist.GetName(), hist.Integral(-1, -1)))

    ymaxval = 0.
    yminval = 1E10
    for i, hist in enumerate(listhist):
        if listlinesty[i] == 1:
            hist.SetLineWidth(3)
        else:
            hist.SetLineWidth(2)
        hist.SetLineColor(TColor.GetColor(lcolor[i]))
        hist.SetLineStyle(listlinesty[i])
        if hist.GetMaximum() > ymaxval:
            ymaxval = hist.GetMaximum()

        if hist.GetMinimum(0) < yminval:
            yminval = hist.GetMinimum(0)

    binwidth = listhist[0].GetXaxis().GetBinWidth(1)

    # yaxtitletext = ''
    # if yaxisunit == '':
    #     yaxtitletext = 'Expected events / {:g}'.format(binwidth)
    # else:
    #     yaxtitletext = 'Expected events / {:g} {}'.format(
    #         binwidth, yaxisunit)
    yaxtitletext = 'Expected events'

    c = TCanvas('c', '', 700, 600)
    c.cd()
    c.SetLogx()
    gPad.SetRightMargin(0.05)
    gPad.SetTopMargin(0.08)
    gPad.SetLeftMargin(0.15)
    gPad.SetBottomMargin(0.15)
    if logy:
        c.SetLogy()
        listhist[0].GetYaxis().SetRangeUser(
            yminval*0.95, ymaxval*axisscale)
    else:
        listhist[0].GetYaxis().SetRangeUser(0, ymaxval*axisscale)

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
    listhist[0].GetYaxis().SetTitleOffset(1.3)
    listhist[0].GetYaxis().SetMoreLogLabels()
    listhist[0].Draw(drawopt)
    for i in range(1, len(listhist)):
        if drawxoverflow:
            listhist[i].GetXaxis().SetRange(1, listhist[i].GetNbinsX() + 1)
        listhist[i].Draw('{} same'.format(drawopt))

    CMS_lumi.CMS_lumi(c, 5, 0, '', 2017, True, 'Simulation',
                      '', '')  # CMS Simulation outside the frame
    c.Modified()
    c.RedrawAxis()

    l = TLegend(0.18, 0.65, 0.92, 0.89)
    l.SetNColumns(2)
    l.SetTextSize(0.03)
    for ileg, leg in enumerate(listleg):
        if ileg == 0:
            l.AddEntry(listhist[ileg], leg, legsty)
            l.AddEntry(0, '','')
        else:
            l.AddEntry(listhist[ileg], leg, legsty)
    l.SetFillColor(0)  # Set the background to be white
    l.SetLineColor(0)
    l.SetFillStyle(0)
    l.Draw('same')

    c.SaveAs('{}/{}.png'.format(opt.plotdir, outname))
    c.SaveAs('{}/{}.pdf'.format(opt.plotdir, outname))
    
def DrawHistStack(listhist, listleg, logy, axisscale, listlegsty, XaxisName, yaxisunit, outname):
    lcolor = ['#1b262c', '#6a2c70', '#4a266a','#7f4a88','#de95ba','#ffd9e8']
    TickSize = 0.03
    AxisTitleSize = 0.05
    AxisLabelSize = 0.05
    
    ymaxval = 0.
    yminval = 1E10
    hs = THStack('hs', 'hs')
    for ih, hist in enumerate(listhist):
        if ih < 2:
            hist.SetLineWidth(3)
            hist.SetLineColor(TColor.GetColor(lcolor[ih]))
        else:
            hist.SetLineWidth(0)
            hist.SetFillColor(TColor.GetColor(lcolor[ih]))
            hs.Add(hist)
            
        if hist.GetMaximum() > ymaxval:
            ymaxval = hist.GetMaximum()

        if hist.GetMinimum(0) < yminval:
            yminval = hist.GetMinimum(0)

    binwidth = listhist[0].GetXaxis().GetBinWidth(1)
    gStyle.SetOptStat(0)

    # yaxtitletext = ''
    # if yaxisunit == '':
    #     yaxtitletext = 'Expected events / {:g}'.format(binwidth)
    # else:
    #     yaxtitletext = 'Expected events / {:g} {}'.format(
    #         binwidth, yaxisunit)
    yaxtitletext = 'Expected events'
    
    c = TCanvas('c', '', 700, 600)
    c.cd()
    c.SetLogx()
    gPad.SetRightMargin(0.05)
    gPad.SetTopMargin(0.08)
    gPad.SetLeftMargin(0.15)
    gPad.SetBottomMargin(0.15)
    if logy:
        c.SetLogy()
        listhist[0].GetYaxis().SetRangeUser(yminval*0.95, ymaxval*axisscale)
        listhist[0].GetYaxis().SetMoreLogLabels()
    else:
        listhist[0].GetYaxis().SetRangeUser(0, ymaxval*axisscale)

    listhist[0].GetXaxis().SetTitle(XaxisName)
    listhist[0].GetYaxis().SetTitle(yaxtitletext)
    listhist[0].GetXaxis().SetTickSize(TickSize)
    listhist[0].GetXaxis().SetTitleSize(AxisTitleSize)
    listhist[0].GetXaxis().SetLabelSize(AxisLabelSize)
    listhist[0].GetYaxis().SetTickSize(TickSize)
    listhist[0].GetYaxis().SetTitleSize(AxisTitleSize)
    listhist[0].GetYaxis().SetLabelSize(AxisLabelSize)
    listhist[0].GetXaxis().SetTitleOffset(1.3)
    listhist[0].GetYaxis().SetTitleOffset(1.3)
    listhist[0].Draw('hist')
    hs.Draw('hist same')
    listhist[1].Draw('hist same')
    CMS_lumi.CMS_lumi(c, 5, 0, '', 2017, True, 'Simulation', '', '')  # CMS Simulation outside the frame
    c.Modified()
    c.RedrawAxis()
    
    l = TLegend(0.18, 0.65, 0.92, 0.89)
    l.SetNColumns(2)
    l.SetTextSize(0.03)
    for ileg, leg in enumerate(listleg):
        l.AddEntry(listhist[ileg], leg, listlegsty[ileg])
    l.SetFillColor(0)  # Set the background to be white
    l.SetLineColor(0)
    l.SetFillStyle(0)
    l.Draw('same')
    
    c.SaveAs('{}/{}.png'.format(opt.plotdir, outname))
    c.SaveAs('{}/{}.pdf'.format(opt.plotdir, outname))

if __name__ == '__main__':
    if opt.verb:
        print('This is the __main__ part')

    NBins_mgs = 35
    xbin_mgs = numpy.logspace(-3, math.log10(50), NBins_mgs+1) # https://numpy.org/doc/stable/reference/generated/numpy.logspace.html
    hM_mgs_tot = TH1F('hM_mgs_tot', 'hM_mgs_tot', NBins_mgs, xbin_mgs)
    hM_mgs_kinsel_diElePho = TH1F('hM_mgs_kinsel_diElePho', 'hM_mgs_kinsel_diElePho', NBins_mgs, xbin_mgs)
    hM_mgs_kinsel_Ele23Ele12 = TH1F('hM_mgs_kinsel_Ele23Ele12', 'hM_mgs_kinsel_Ele23Ele12', NBins_mgs, xbin_mgs)
    hM_mgs_diElePho = TH1F('hM_mgs_diElePho', 'hM_mgs_diElePho', NBins_mgs, xbin_mgs)
    hM_mgs_Ele23Ele12 = TH1F('hM_mgs_Ele23Ele12', 'hM_mgs_Ele23Ele12', NBins_mgs, xbin_mgs)
    hM_mgs_combHLT = TH1F('hM_mgs_combHLT', 'hM_mgs_combHLT', NBins_mgs, xbin_mgs)
    hM_mgs_combHLT_resolved = TH1F('hM_mgs_combHLT_resolved', 'hM_mgs_combHLT_resolved', NBins_mgs, xbin_mgs)
    hM_mgs_combHLT_merged2Gsf = TH1F('hM_mgs_combHLT_merged2Gsf', 'hM_mgs_combHLT_merged2Gsf', NBins_mgs, xbin_mgs)
    hM_mgs_combHLT_merged1Gsf = TH1F('hM_mgs_combHLT_merged1Gsf', 'hM_mgs_combHLT_merged1Gsf', NBins_mgs, xbin_mgs)
    hM_mgs_combHLT_NPR = TH1F('hM_mgs_combHLT_NPR', 'hM_mgs_combHLT_NPR', NBins_mgs, xbin_mgs)
    NBins_dRee = 35
    xbin_dRee = numpy.logspace(math.log10(0.00005), math.log10(3), NBins_dRee+1) # https://numpy.org/doc/stable/reference/generated/numpy.logspace.html
    hM_dRee_tot = TH1F('hM_dRee_tot', 'hM_dRee_tot', NBins_dRee, xbin_dRee)
    hM_dRee_kinsel_diElePho = TH1F('hM_dRee_kinsel_diElePho', 'hM_dRee_kinsel_diElePho', NBins_dRee, xbin_dRee)
    hM_dRee_kinsel_Ele23Ele12 = TH1F('hM_dRee_kinsel_Ele23Ele12', 'hM_dRee_kinsel_Ele23Ele12', NBins_dRee, xbin_dRee)
    hM_dRee_diElePho = TH1F('hM_dRee_diElePho', 'hM_dRee_diElePho', NBins_dRee, xbin_dRee)
    hM_dRee_Ele23Ele12 = TH1F('hM_dRee_Ele23Ele12', 'hM_dRee_Ele23Ele12', NBins_dRee, xbin_dRee)
    hM_dRee_combHLT = TH1F('hM_dRee_combHLT', 'hM_dRee_combHLT', NBins_dRee, xbin_dRee)
    hM_dRee_combHLT_resolved = TH1F('hM_dRee_combHLT_resolved', 'hM_dRee_combHLT_resolved', NBins_dRee, xbin_dRee)
    hM_dRee_combHLT_merged2Gsf = TH1F('hM_dRee_combHLT_merged2Gsf', 'hM_dRee_combHLT_merged2Gsf', NBins_dRee, xbin_dRee)
    hM_dRee_combHLT_merged1Gsf = TH1F('hM_dRee_combHLT_merged1Gsf', 'hM_dRee_combHLT_merged1Gsf', NBins_dRee, xbin_dRee)
    hM_dRee_combHLT_NPR = TH1F('hM_dRee_combHLT_NPR', 'hM_dRee_combHLT_NPR', NBins_dRee, xbin_dRee)

    chain = TChain('outTree')
    chain.Add('./minitree/Minitree_HDalitz_ggF_eeg_m125_2017_RECO.root')
    chain.Add('./minitree/Minitree_HDalitz_VBF_eeg_m125_2017_RECO.root')
    chain.Add('./minitree/Minitree_HDalitz_WH_eeg_m125_2017_RECO.root')
    chain.Add('./minitree/Minitree_HDalitz_ZH_eeg_m125_2017_RECO.root')
    for index, ev in enumerate(chain):
        if index == 0:
            print('Whether the Ele23_Ele12 is prescaled?', ((
                (ev.HLTEleMuXIsPrescaled >> 5) & 1) == 1), (((ev.HLTEleMuXIsPrescaled >> 40) & 1) == 1))

        gen_lep1 = TLorentzVector()
        gen_lep1.SetPtEtaPhiM(ev.mcPt_lep1, ev.mcEta_lep1,
                              ev.mcPhi_lep1, 0.000510999)
        gen_lep2 = TLorentzVector()
        gen_lep2.SetPtEtaPhiM(ev.mcPt_lep2, ev.mcEta_lep2,
                              ev.mcPhi_lep2, 0.000510999)
        gen_pho = TLorentzVector()
        gen_pho.SetPtEtaPhiM(ev.mcPt_pho, ev.mcEta_pho, ev.mcPhi_pho, 0.0)

        hM_mgs_tot.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
        hM_dRee_tot.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)

        # Basic kinematic selections for diEle27 and diPho30and22
        passKinSel_diElePho = (abs((gen_lep1 + gen_lep2).Eta() < 2.5) and
                               abs(gen_pho.Eta() < 2.5) and
                               (gen_lep1 + gen_lep2).Pt() > 125 * 0.3 and
                               gen_pho.Pt() > 125 * 0.3)
        passKinSel_Ele23Ele12 = (gen_lep1.Pt() > 25. and gen_lep2.Pt() > 15.)

        if passKinSel_diElePho == True:
            hM_mgs_kinsel_diElePho.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
            hM_dRee_kinsel_diElePho.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)
            # HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v ==> 13
            # HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v ==> 14
            if ((ev.HLTPho >> 13) & 1) == 1 or ((ev.HLTPho >> 14) & 1) == 1:
                hM_mgs_diElePho.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
                hM_dRee_diElePho.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)

        if passKinSel_Ele23Ele12 == True:
            hM_mgs_kinsel_Ele23Ele12.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
            hM_dRee_kinsel_Ele23Ele12.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)
        #  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v OR HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v
            if ((ev.HLTEleMuX >> 5) & 1) == 1 or ((ev.HLTEleMuX >> 40) & 1) == 1:
                hM_mgs_Ele23Ele12.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
                hM_dRee_Ele23Ele12.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)

        passdiElePho = passKinSel_diElePho == True and (((ev.HLTPho >> 13) & 1) == 1 or ((ev.HLTPho >> 14) & 1) == 1)
        passEle23Ele12 = passKinSel_Ele23Ele12 == True and (((ev.HLTEleMuX >> 5) & 1) == 1 or ((ev.HLTEleMuX >> 40) & 1) == 1)
        if passdiElePho or passEle23Ele12:
            hM_mgs_combHLT.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
            hM_dRee_combHLT.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)
            if ev.category == 1:
                hM_mgs_combHLT_resolved.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
                hM_dRee_combHLT_resolved.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)
            elif ev.category == 2:
                hM_mgs_combHLT_merged2Gsf.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
                hM_dRee_combHLT_merged2Gsf.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)
            elif ev.category == 3:
                hM_mgs_combHLT_merged1Gsf.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
                hM_dRee_combHLT_merged1Gsf.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)
            elif ev.category == 4:
                hM_mgs_combHLT_NPR.Fill((gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
                hM_dRee_combHLT_NPR.Fill(gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)
            else:
                print ('[WARNING] Category',ev.category,'is not defined. Check the code.\n')
            

    lhist_mgs = [hM_mgs_tot, hM_mgs_kinsel_diElePho, hM_mgs_diElePho,hM_mgs_kinsel_Ele23Ele12, hM_mgs_Ele23Ele12, hM_mgs_combHLT]
    lhist_dRee = [hM_dRee_tot, hM_dRee_kinsel_diElePho, hM_dRee_diElePho,hM_dRee_kinsel_Ele23Ele12, hM_dRee_Ele23Ele12, hM_dRee_combHLT]
    lleg = ['Before selections', '#splitline{Kinematic selections}{(diEle27&diPho)}', 'HLT_diEle27 OR HLT_diPho',
            '#splitline{Kinematic selections}{(Ele23Ele12)}', 'HLT_Ele23Ele12', 'All HLTs combined']
    llinesty = [1, 2, 7, 2, 7, 1]
    # Drawhist(listhist, listleg, drawopt, drawxoverflow, drawxunderflow, logy, axisscale, legsty, XaxisName, yaxisunit, outname)
    Drawhist(lhist_mgs, lleg, llinesty, 'HIST', False, False, False, 1.7, 'l', 'GEN-level m_{ee} (GeV)', 'GeV', 'ExpYield_Mgs_allHLT')
    Drawhist(lhist_dRee, lleg, llinesty, 'HIST', False, False, False, 1.7, 'l', 'GEN-level #DeltaR(e,e)', '', 'ExpYield_dRee_allHLT')
    
    lhist_mgs_stack = [hM_mgs_tot, hM_mgs_combHLT, hM_mgs_combHLT_resolved, hM_mgs_combHLT_merged2Gsf, hM_mgs_combHLT_merged1Gsf, hM_mgs_combHLT_NPR]
    lhist_dRee_stack = [hM_dRee_tot, hM_dRee_combHLT, hM_dRee_combHLT_resolved, hM_dRee_combHLT_merged2Gsf, hM_dRee_combHLT_merged1Gsf, hM_dRee_combHLT_NPR]
    lleg_stack = ['Before selections', 'After combined HLT', '#splitline{After combined HLT}{Resolved}', '#splitline{After combined HLT}{Merged2Gsf}', '#splitline{After combined HLT}{Merge1Gsf}', '#splitline{After combined HLT}{NPR}']
    llegsty_stack = ['l', 'l', 'f', 'f', 'f', 'f']
    # DrawHistStack(listhist, listleg, logy, axisscale, listlegsty, XaxisName, yaxisunit, outname):
    DrawHistStack(lhist_mgs_stack, lleg_stack, False, 1.0, llegsty_stack, 'GEN-level m_{ee} (GeV)', 'GeV', 'ExpYield_Mgs_allHLT_StackCat')
    DrawHistStack(lhist_dRee_stack, lleg_stack, False, 1.0, llegsty_stack, 'GEN-level #DeltaR(e,e)', '', 'ExpYield_dRee_allHLT_StackCat')