import ROOT
import os
import time
from CMS_lumi import CMS_lumi
from glob import glob
import numpy as np
ROOT.gROOT.LoadMacro("../interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")


def Draw_effcomp(vgr, vleg, xlim, ylim, legsty, XaxisName, YaxisName, outname):
    # lcolor = ["#202020", "#2666CF", "#CC4054", "#019267"]
    lcolor = [6, 4, 2, 1]
    makerstle = [24, 23, 21, 20]

    for i, gr in enumerate(vgr):
        gr.SetLineWidth(2)
        # gr.SetLineColor(ROOT.TColor.GetColor(lcolor[i]))
        gr.SetLineColor(lcolor[i])
        gr.SetMarkerStyle(makerstle[i])
        # gr.SetMarkerColor(ROOT.TColor.GetColor(lcolor[i]))
        gr.SetMarkerColor(lcolor[i])
        gr.SetMarkerSize(1.3)

    TickSize = 0.03
    AxisTitleSize = 0.05
    AxisLabelSize = 0.05

    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas("c", "", 700, 600)
    c.cd()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetBottomMargin(0.14)
    vgr[0].GetXaxis().SetLimits(xlim[0], xlim[1])
    vgr[0].SetMinimum(ylim[0])
    vgr[0].SetMaximum(ylim[1])
    vgr[0].GetXaxis().SetNdivisions(504)
    vgr[0].GetXaxis().SetTitle(XaxisName)
    vgr[0].GetYaxis().SetTitle(YaxisName)
    vgr[0].GetXaxis().SetTickSize(TickSize)
    vgr[0].GetXaxis().SetTitleSize(AxisTitleSize)
    vgr[0].GetXaxis().SetLabelSize(AxisLabelSize)
    vgr[0].GetYaxis().SetTickSize(TickSize)
    vgr[0].GetYaxis().SetTitleSize(AxisTitleSize)
    vgr[0].GetYaxis().SetLabelSize(AxisLabelSize)
    vgr[0].GetXaxis().SetTitleOffset(1.3)
    vgr[0].GetYaxis().SetTitleOffset(1.3)
    vgr[0].Draw("APE1")
    for i, gr in enumerate(vgr):
        if i == 0:
            continue
        else:
            vgr[i].Draw("PE1 same")
            
    CMS_lumi(c, 5, 0, "", 2017, True, "Simulation", "", "")
    c.RedrawAxis()
    c.Modified()

    l = ROOT.TLegend(0.17, 0.73, 0.4, 0.88)
    l.SetTextSize(0.028)
    for i, gr in enumerate(vgr):
        l.AddEntry(gr, vleg[i], legsty)
    l.SetTextFont(42)
    l.SetFillColor(0)  # Set the background to be white
    l.SetLineColor(0)
    l.Draw("same")

    c.SaveAs("{}.png".format(outname))
    c.SaveAs("{}.pdf".format(outname))

    

if __name__ == "__main__" :
    start_time = time.time()
    
    # PyROOT does not display any graphics(root "-b" option)
    ROOT.gROOT.SetBatch()
    
    ROOT.EnableImplicitMT(10)
    infile_sigs = glob("../miniTree/UL2018/miniTree_HDalitz_*_eeg_*.root")
    df = ROOT.RDataFrame("miniTree", infile_sigs).Define("weight", "mcwei*genwei*puwei")\
             .Define("dRee", "ROOT::VecOps::DeltaR(GenEle_Lead.Eta(), GenEle_subLead.Eta(), GenEle_Lead.Phi(), GenEle_subLead.Phi())")\
             .Filter("diGenEle.Pt() > 70 && GenPho_Lead.Pt() > 70 && diGenEle.M() < 10 && GenEle_Lead.Pt() > 15 && GenEle_subLead.Pt() > 15", "kin cuts")
    report = df.Report()
    
    h = df.Histo1D(("h1", " ", 15, 0, 0.3), "dRee", "weight")
    h1 = df.Filter("elePresel_Lead == 1 && phoPresel_Lead == 1").Histo1D(("h1", " ", 15, 0, 0.3), "dRee", "weight") # HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90_v
    # h2 = df.Filter("((HLTPho >> 13) & 1) == 1").Histo1D(("h2", " ", 15, 0, 0.3), "dRee", "weight") # HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v
    h3 = df.Filter("((HLTPho >> 22) & 1) == 1").Histo1D(("h3", " ", 15, 0, 0.3), "dRee", "weight") # HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v
    h4 = df.Filter("((HLTPho >> 13) & 1) == 1 || (elePresel_Lead == 1 && phoPresel_Lead == 1) ").Histo1D(("h4", " ", 15, 0, 0.3), "dRee", "weight") # HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v or hgg
    
    err1 = ROOT.TGraphAsymmErrors(h1.GetPtr(), h.GetPtr(), "cl=0.683 b(1,1) mode")
    # err2 = ROOT.TGraphAsymmErrors(h2.GetPtr(), h.GetPtr(), "cl=0.683 b(1,1) mode")
    err3 = ROOT.TGraphAsymmErrors(h3.GetPtr(), h.GetPtr(), "cl=0.683 b(1,1) mode")
    err4 = ROOT.TGraphAsymmErrors(h4.GetPtr(), h.GetPtr(), "cl=0.683 b(1,1) mode")
    
    l_eff_dRee = [
        err4, 
        err1, 
        # err2,
        err3
    ]
    llegend = [
        "Passes any of two triggers below",
        "H #rightarrow #gamma #gamma preselection", 
        # "HLT_DoublePhoton60", 
        "HLT_DiEle27_WPTightCaloOnly_L1DoubleEG_v",
    ]
    Draw_effcomp(l_eff_dRee, llegend, [0, 0.32], [
                 0, 1.25], "AP", "#DeltaR(e,e)", "Trigger Efficiency", "../plots/HLT_sigeff_dRee"
    )
    report.Print()
    
    # HLTEleMuX
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))