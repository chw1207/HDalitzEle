import ROOT
import os
import time
from CMS_lumi import CMS_lumi
from sigmaEff import sigmaEff
from glob import glob
import numpy as np


def Draw1DHist(c, vh, vc, vcl , xaxis="x-axis", yaxis="x-axis", option="hist", Log=True):
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.07)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetBottomMargin(0.15)
    if (Log == True):
        c.SetLogy()

    # Set the axis style
    if (Log == True):
        ymax = vh[0].GetBinContent(vh[0].GetMaximumBin()) * 10
        ymin = 1E-1
    else: 
        ymax = vh[0].GetBinContent(vh[0].GetMaximumBin()) * 1.5
        ymin = 0
    vh[0].SetMarkerStyle(20)
    vh[0].SetMarkerSize(1.2)
    vh[0].GetXaxis().SetTitle(xaxis)
    vh[0].GetXaxis().SetMoreLogLabels()
    # vh[0].GetXaxis().SetTickSize(0.02)
    vh[0].GetXaxis().SetTitleSize(0.05)
    vh[0].GetXaxis().SetLabelSize(0.045)
    vh[0].GetXaxis().SetLabelOffset(0.02)
    vh[0].GetXaxis().SetTitleOffset(1.4)
    vh[0].GetYaxis().SetTitle(yaxis)
    vh[0].GetYaxis().SetRangeUser(ymin, ymax)
    # vh[0].GetYaxis().SetNdivisions(506)
    # vh[0].GetYaxis().SetTickSize(0.02)
    vh[0].GetYaxis().SetTitleSize(0.05)
    vh[0].GetYaxis().SetLabelSize(0.045)
    # vh[0].GetYaxis().SetLabelOffset(0.02)
    vh[0].GetYaxis().SetTitleOffset(1.4)

    # Set the color style and draw option
    for i, h in enumerate(vh):
        h.SetLineColor(ROOT.TColor.GetColor(vcl[i]))
        h.SetLineWidth(4)
        if (vc[i] != None):
            h.SetFillColor(ROOT.TColor.GetColor(vc[i]))
        if (i == 0):
            h.Draw(option)
        else:
            h.Draw("%s same" %(option))


def main(_track):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    h1 = rdf_sigs.Histo1D(("h1", " ", 7, -0.5, 6.5), _track, "weight")
    h2 = rdf_gjet.Histo1D(("h2", " ", 7, -0.5, 6.5), _track, "weight")
    s1 = h1.GetBinContent(1)/np.sqrt(h2.GetBinContent(1))
    s2 = (h1.GetBinContent(1)+h1.GetBinContent(2))/np.sqrt(h2.GetBinContent(1)+h2.GetBinContent(2))
    print("significance for missing hit = 0: {}".format(s1))
    print("significance for missing hit = 1: {}".format(s2))
    
    h1.Scale(1/h1.Integral())
    h2.Scale(1/h2.Integral())
    c1 = ROOT.TCanvas("c1", "", 800, 800)
    c1.cd()
    
    if _track == "eleTrkMissHits_Lead":
        xName = "Missing hits(main Gsf track)"
    else:
        xName = "Missing hits(attached Gsf track)"
    Draw1DHist(c1, [h1, h2], [None, None], ["#E69D45", "#0061a8"], xaxis=xName, yaxis="A.U.", option="hist", Log=False)

    leg = ROOT.TLegend(0.55, 0.7, 0.9, 0.9)
    leg.SetHeader("Merged-2Gsf")
    leg.SetTextSize(0.04)
    leg.AddEntry(h1.GetPtr(), "Signal", "f")
    leg.AddEntry(h2.GetPtr(), "#gamma converts to e", "f")
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.Draw()
    
    ltx = ROOT.TLatex()
    ltx.SetNDC()
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.035)
    ltx.DrawLatex(0.5, 0.5, "#frac{S}{#sqrt{B}}(missing hits = 0) = %.4f" %(s1))
    ltx.DrawLatex(0.5, 0.4, "#frac{S}{#sqrt{B}}(missing hits = 1) = %.4f" %(s2))

    c1.RedrawAxis()
    CMS_lumi(c1, 5, 10, "", 2017, True, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    
    os.makedirs("../plots/GenStudy", exist_ok=True)
    outName = "../plots/GenStudy/missHits_choice_{}.pdf".format(_track)
    c1.Print(outName)
    c1.Close()

if __name__ == "__main__" :
    start_time = time.time()
    
    # PyROOT does not display any graphics(root "-b" option)
    ROOT.gROOT.SetBatch()
    
    ROOT.EnableImplicitMT(10)
    infile_sigs = glob("../miniTree/*/miniTree_HDalitz_*_eeg_*.root")
    rdf_sigs = ROOT.RDataFrame("miniTree", infile_sigs).Define("weight", "mcwei * genwei")\
                   .Filter("elePresel_Lead == 1 && category == 2 && nGsfMatchToReco_Lead >= 2")
    
    infile_gjet = glob("../miniTree/*/miniTree_GJets_*.root")
    rdf_gjet = ROOT.RDataFrame("miniTree", infile_gjet).Define("weight", "mcwei * genwei")\
                   .Filter("elePresel_Lead == 1 && nGsfMatchToReco_Lead >= 2")
    
    for track in ["eleTrkMissHits_Lead", "eleSubTrkMissHits_Lead"]:    
        main(track)

    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))