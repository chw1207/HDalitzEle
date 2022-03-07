#!/usr/bin/env python

import argparse
import ROOT
import os
import sys
import json
from array import array

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
    '/afs/cern.ch/work/h/hajheng/private/HDalitz/interface/tdrstyle.C')
ROOT.gROOT.ProcessLine('setTDRStyle()')


parser = argparse.ArgumentParser(
    description='Plot reco-level variables of electron')
parser.add_argument('-v', '--verb', dest='verb', action='store_true', default=False,
                    help='Do more verbosity printout.')
parser.add_argument('-d', '--plotdir', dest='plotdir', type=str, default='./plots/sigeff_HLT',
                    help='Directory to store plots')
opt = parser.parse_args()

if not opt.verb:
    ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = 1001')

os.makedirs(opt.plotdir, exist_ok=True)


def Draw_effcomp(vgr, vleg, xlim, ylim, legsty, XaxisName, YaxisName, outname):
    lcolor = ['#4D4C4C', '#8298CD', '#C06463']

    for i, gr in enumerate(vgr):
        gr.SetLineWidth(2)
        gr.SetLineColor(ROOT.TColor.GetColor(lcolor[i]))
        gr.SetMarkerStyle(20)
        gr.SetMarkerColor(ROOT.TColor.GetColor(lcolor[i]))
        gr.SetMarkerSize(1.3)

    TickSize = 0.03
    AxisTitleSize = 0.05
    AxisLabelSize = 0.05

    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas('c', '', 700, 600)
    c.cd()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetLeftMargin(0.15)
    ROOT.gPad.SetBottomMargin(0.15)
    vgr[0].GetXaxis().SetLimits(xlim[0], xlim[1])
    vgr[0].SetMinimum(ylim[0])
    vgr[0].SetMaximum(ylim[1])
    vgr[0].GetXaxis().SetTitle(XaxisName)
    vgr[0].GetYaxis().SetTitle(YaxisName)
    vgr[0].GetXaxis().SetTickSize(TickSize)
    vgr[0].GetXaxis().SetTitleSize(AxisTitleSize)
    vgr[0].GetXaxis().SetLabelSize(AxisLabelSize)
    vgr[0].GetYaxis().SetTickSize(TickSize)
    vgr[0].GetYaxis().SetTitleSize(AxisTitleSize)
    vgr[0].GetYaxis().SetLabelSize(AxisLabelSize)
    vgr[0].GetXaxis().SetTitleOffset(1.3)
    vgr[0].GetYaxis().SetTitleOffset(1.4)
    vgr[0].Draw('APE1')
    for i, gr in enumerate(vgr):
        if i == 0:
            continue
        else:
            vgr[i].Draw('PE1 same')

    CMS_lumi.CMS_lumi(c, 5, 0, '', 2017, True, 'Simulation', '', '')
    c.RedrawAxis()
    c.Modified()

    l = ROOT.TLegend(0.6, 0.2, 0.93, 0.4)
    l.SetTextSize(0.035)
    for i, gr in enumerate(vgr):
        l.AddEntry(gr, vleg[i], legsty)
    l.SetFillColor(0)  # Set the background to be white
    l.SetLineColor(1)
    l.Draw('same')

    c.SaveAs('{}.png'.format(outname))
    c.SaveAs('{}.pdf'.format(outname))


if __name__ == '__main__':
    if opt.verb:
        print('This is the __main__ part')

    mgs_edge = [0., 1., 2., 3., 4., 5., 6., 7., 8.,
                9., 10., 15., 20., 30., 40., 50., 60.]
    dRee_edge = [0, 0.03, 0.06, 0.09, 0.12, 0.16, 0.22,
                 0.3, 0.4, 0.55, 0.7, 0.85, 1, 1.5, 2.0, 3.0]
    NBins_mgs = len(mgs_edge) - 1
    NBins_dRee = len(dRee_edge) - 1

    hM_mgs_tot = ROOT.TH1F('hM_mgs_tot', 'hM_mgs_tot',
                           NBins_mgs, array('d', mgs_edge))
    hM_dRee_tot = ROOT.TH1F('hM_dRee_tot', 'hM_dRee_tot',
                            NBins_dRee, array('d', dRee_edge))
    hM_mgs_kinsel_diElePho = ROOT.TH1F(
        'hM_mgs_kinsel_diElePho', 'hM_mgs_kinsel_diElePho', NBins_mgs, array('d', mgs_edge))
    hM_dRee_kinsel_diElePho = ROOT.TH1F(
        'hM_dRee_kinsel_diElePho', 'hM_dRee_kinsel_diElePho', NBins_dRee, array('d', dRee_edge))
    hM_mgs_kinsel_Ele23Ele12 = ROOT.TH1F(
        'hM_mgs_kinsel_Ele23Ele12', 'hM_mgs_kinsel_Ele23Ele12', NBins_mgs, array('d', mgs_edge))
    hM_dRee_kinsel_Ele23Ele12 = ROOT.TH1F(
        'hM_dRee_kinsel_Ele23Ele12', 'hM_dRee_kinsel_Ele23Ele12', NBins_dRee, array('d', dRee_edge))
    hM_mgs_HLTdiEle27 = ROOT.TH1F(
        'hM_mgs_HLTdiEle27', 'hM_mgs_HLTdiEle27', NBins_mgs, array('d', mgs_edge))
    hM_dRee_HLTdiEle27 = ROOT.TH1F(
        'hM_dRee_HLTdiEle27', 'hM_dRee_HLTdiEle27', NBins_dRee, array('d', dRee_edge))
    hM_mgs_HLTdiPho30and22 = ROOT.TH1F(
        'hM_mgs_HLTdiPho30and22', 'hM_mgs_HLTdiPho30and22', NBins_mgs, array('d', mgs_edge))
    hM_dRee_HLTdiPho30and22 = ROOT.TH1F(
        'hM_dRee_HLTdiPho30and22', 'hM_dRee_HLTdiPho30and22', NBins_dRee, array('d', dRee_edge))
    hM_mgs_HLTEle23Ele12 = ROOT.TH1F(
        'hM_mgs_HLTEle23Ele12', 'hM_mgs_HLTEle23Ele12', NBins_mgs, array('d', mgs_edge))
    hM_dRee_HLTEle23Ele12 = ROOT.TH1F(
        'hM_dRee_HLTEle23Ele12', 'hM_dRee_HLTEle23Ele12', NBins_dRee, array('d', dRee_edge))

    chain = ROOT.TChain('outTree')
    chain.Add('./minitree/Minitree_HDalitz_ggF_eeg_m125_2017_RECO.root')
    chain.Add('./minitree/Minitree_HDalitz_VBF_eeg_m125_2017_RECO.root')
    chain.Add('./minitree/Minitree_HDalitz_WH_eeg_m125_2017_RECO.root')
    chain.Add('./minitree/Minitree_HDalitz_ZH_eeg_m125_2017_RECO.root')
    for index, ev in enumerate(chain):
        if index == 0:
            print('Whether the Ele23_Ele12 is prescaled?', ((
                (ev.HLTEleMuXIsPrescaled >> 5) & 1) == 1), (((ev.HLTEleMuXIsPrescaled >> 40) & 1) == 1))

        gen_lep1 = ROOT.TLorentzVector()
        gen_lep1.SetPtEtaPhiM(ev.mcPt_lep1, ev.mcEta_lep1,
                              ev.mcPhi_lep1, 0.000510999)
        gen_lep2 = ROOT.TLorentzVector()
        gen_lep2.SetPtEtaPhiM(ev.mcPt_lep2, ev.mcEta_lep2,
                              ev.mcPhi_lep2, 0.000510999)
        gen_pho = ROOT.TLorentzVector()
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
            hM_mgs_kinsel_diElePho.Fill(
                (gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
            hM_dRee_kinsel_diElePho.Fill(
                gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)

            # HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v
            if ((ev.HLTPho >> 13) & 1) == 1:
                hM_mgs_HLTdiEle27.Fill(
                    (gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
                hM_dRee_HLTdiEle27.Fill(gen_lep1.DeltaR(
                    gen_lep2), ev.mcwei * ev.genwei)

            # HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v
            if ((ev.HLTPho >> 14) & 1) == 1:
                hM_mgs_HLTdiPho30and22.Fill(
                    (gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
                hM_dRee_HLTdiPho30and22.Fill(
                    gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)

        if passKinSel_Ele23Ele12 == True:
            hM_mgs_kinsel_Ele23Ele12.Fill(
                (gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
            hM_dRee_kinsel_Ele23Ele12.Fill(
                gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)

        #  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v OR HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v
            if ((ev.HLTEleMuX >> 5) & 1) == 1 or ((ev.HLTEleMuX >> 40) & 1) == 1:
                hM_mgs_HLTEle23Ele12.Fill(
                    (gen_lep1 + gen_lep2).M(), ev.mcwei * ev.genwei)
                hM_dRee_HLTEle23Ele12.Fill(
                    gen_lep1.DeltaR(gen_lep2), ev.mcwei * ev.genwei)

    eff_mgs_HLTdiEle27 = ROOT.TGraphAsymmErrors(
        hM_mgs_HLTdiEle27, hM_mgs_kinsel_diElePho)
    eff_mgs_HLTdiPho30and22 = ROOT.TGraphAsymmErrors(
        hM_mgs_HLTdiPho30and22, hM_mgs_kinsel_diElePho)
    eff_mgs_HLTEle23Ele12 = ROOT.TGraphAsymmErrors(
        hM_mgs_HLTEle23Ele12, hM_mgs_kinsel_Ele23Ele12)
    eff_dRee_HLTdiEle27 = ROOT.TGraphAsymmErrors(
        hM_dRee_HLTdiEle27, hM_dRee_kinsel_diElePho)
    eff_dRee_HLTdiPho30and22 = ROOT.TGraphAsymmErrors(
        hM_dRee_HLTdiPho30and22, hM_dRee_kinsel_diElePho)
    eff_dRee_HLTEle23Ele12 = ROOT.TGraphAsymmErrors(
        hM_dRee_HLTEle23Ele12, hM_dRee_kinsel_Ele23Ele12)

    l_eff_mgs = [eff_mgs_HLTdiEle27,
                 eff_mgs_HLTdiPho30and22, eff_mgs_HLTEle23Ele12]
    l_eff_dRee = [eff_dRee_HLTdiEle27,
                  eff_dRee_HLTdiPho30and22, eff_dRee_HLTEle23Ele12]
    llegend = ['diEle27', 'diPho30and22', 'Ele23Ele12']

    Draw_effcomp(l_eff_mgs, llegend, [0, 60], [
                 0, 1], 'PLEX0', 'm_{#gamma*} (GeV)', 'Efficiency', './plots/sigeff_HLT/HLT_eigeff_mgs')
    Draw_effcomp(l_eff_dRee, llegend, [0, 1], [
                 0, 1], 'PLEX0', '#DeltaR(e,e)', 'Efficiency', './plots/sigeff_HLT/HLT_eigeff_dRee')

    print('Expected yield before kinematic cuts = ', hM_mgs_tot.Integral(-1, -1), '\nExpected yield after kinematic cuts (for HLT_diEle and HLT_diPho) = ', hM_mgs_kinsel_diElePho.Integral(-1, -1), '\nExpected yield after HLT_diEle27 = ', hM_mgs_HLTdiEle27.Integral(-1, -1),
          '\nExpected yield after HLT_diPho30and22 = ', hM_mgs_HLTdiPho30and22.Integral(-1, -1), '\nExpected yield after kinematic cuts (for HLT_Ele23Ele12) = ', hM_mgs_kinsel_Ele23Ele12.Integral(-1, -1), '\nExpected yield after HLT_Ele23Ele12 = ', hM_mgs_HLTEle23Ele12.Integral(-1, -1))
