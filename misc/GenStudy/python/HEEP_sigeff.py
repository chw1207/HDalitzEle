import ROOT
import os
import time
from CMS_lumi import CMS_lumi
from glob import glob
import numpy as np
ROOT.gROOT.LoadMacro("../interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")


def Draw_effcomp(vgr, vleg, xlim, ylim, legsty, XaxisName, YaxisName, outname):
    lcolor = ["#202020", "#2666CF", "#CC4054", "#019267"]

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

    c = ROOT.TCanvas("c", "", 700, 600)
    c.cd()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.08)
    ROOT.gPad.SetLeftMargin(0.13)
    ROOT.gPad.SetBottomMargin(0.14)
    vgr[0].GetXaxis().SetLimits(xlim[0], xlim[1])
    vgr[0].SetMinimum(ylim[0])
    vgr[0].SetMaximum(ylim[1])
    # vgr[0].GetXaxis().SetNdivisions(504)
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
            
    CMS_lumi(c, 5, 10, "", 2017, True, "Simulation", "H#rightarrow#gamma*#gamma#rightarrowee#gamma", "")
    c.RedrawAxis()
    c.Modified()

    ltx = ROOT.TLatex()
    ltx.SetNDC()
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.045)
    ltx.DrawLatex(0.17, 0.72, "|#eta^{#gamma*}| < 2.5, |#eta^{#gamma}| < 2.5, p^{#gamma*, #gamma}_{T} > 125#times0.3, m_{ee} < 5 GeV")
    
    # abs(diGenEle.Eta()) < 2.5 && abs(GenPho_Lead.Eta()) < 2.5 && diGenEle.Pt() > 125*0.3 && GenPho_Lead.Pt() > 125*0.3 && diGenEle.M() < 5
    
    # l = ROOT.TLegend(0.14, 0.7, 0.4, 0.9)
    # l.SetTextSize(0.025)
    # for i, gr in enumerate(vgr):
    #     l.AddEntry(gr, vleg[i], legsty)
    # l.SetTextFont(42)
    # l.SetFillColor(0)  # Set the background to be white
    # l.SetLineColor(0)
    # l.Draw("same")

    c.SaveAs("{}.png".format(outname))
    c.SaveAs("{}.pdf".format(outname))

    

if __name__ == "__main__" :
    start_time = time.time()
    
    # PyROOT does not display any graphics(root "-b" option)
    ROOT.gROOT.SetBatch()
    
    ROOT.EnableImplicitMT(10)
    infile_sigs = glob("../miniTree/*/miniTree_HDalitz_*_eeg_*.root")
    df = ROOT.RDataFrame("miniTree", infile_sigs).Define("weight", "mcwei*puwei")\
             .Define("dRee", "ROOT::VecOps::DeltaR(GenEle_Lead.Eta(), GenEle_subLead.Eta(), GenEle_Lead.Phi(), GenEle_subLead.Phi())")\
             .Filter("abs(diGenEle.Eta()) < 2.5 && abs(GenPho_Lead.Eta()) < 2.5 && diGenEle.Pt() > 125*0.3 && GenPho_Lead.Pt() > 125*0.3 && diGenEle.M() < 5", "kin cuts")
    report = df.Report()
    
    # bin = [0, 0.03, 0.06, 0.09, 0.12, 0.16, 0.2, 0.24,
    #              0.3, 0.4, 0.5, 0.7,  1]
    # h = df.Histo1D(("h1", " ", len(bin) - 1, np.array(bin, dtype=np.float64)), "dRee", "weight")
    # h1 = df.Filter("((eleIDbit_Lead >> 4) & 1) == 1").Histo1D(("h1", " ", len(bin) - 1, np.array(bin, dtype=np.float64)), "dRee", "weight") 
    
    h = df.Histo1D(("h1", " ", 50, 0, 0.06), "dRee", "weight")
    h1 = df.Filter("((eleIDbit_Lead >> 4) & 1) == 1").Histo1D(("h1", " ", 50, 0, 0.06), "dRee", "weight") 
    
    err1 = ROOT.TGraphAsymmErrors(h1.GetPtr(), h.GetPtr(), "cl=0.683 b(1,1) mode")
    
    l_eff_dRee = [err1]
    llegend = [
        ""
    ]
    Draw_effcomp(l_eff_dRee, llegend, [0, 0.06], [
                 0, 1.22], "PLEX0", "#DeltaR(e,e)", "HEEP ID Efficiency", "../plots/HEEP_sigeff_dRee"
    )
    report.Print()
    
    # HLTEleMuX
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))