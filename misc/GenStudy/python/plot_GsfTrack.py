import ROOT
import os
import time
from CMS_lumi import CMS_lumi
from sigmaEff import sigmaEff
from glob import glob
import numpy as np

# script to determine the way to pick second gsf track associated to the electron cluster.
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
        h.SetLineWidth(3)
        if (vc[i] != None):
            h.SetFillColor(ROOT.TColor.GetColor(vc[i]))
        if (i == 0):
            h.Draw(option)
        else:
            h.Draw("%s same" %(option))
            
def main():
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    h1 = rdf.Filter("category == 2").Histo1D(("h1", " ", 80, 0, 2), "meeRatio", "weight")
    h2 = rdf.Filter("category == 2").Histo1D(("h2", " ", 80, 0, 2), "meeRatioPtMax", "weight")

    c1 = ROOT.TCanvas("c1", "", 800, 800)
    c1.cd()
    Draw1DHist(c1, [h1, h2], [None, None], ["#E69D45", "#0061a8"], xaxis="M^{gg}_{Reco} / M^{ee}_{Gen}", yaxis="Events", option="hist", Log=True)
    CMS_lumi(c1, 5, 0, "", 2017, True, "Simulation", "", "")

    leg = ROOT.TLegend(0.2, 0.75, 0.8, 0.9)
    leg.SetTextSize(0.04)
    leg.AddEntry(h1.GetPtr(), "Two smallest dR tracks, #sigma_{eff} = %.2f" %sigma_dr, "f")
    leg.AddEntry(h2.GetPtr(), "Two highest pT tracks, #sigma_{eff} = %.2f" %sigma_pt, "f")
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.Draw()

    c1.RedrawAxis()
    
    os.makedirs("../plots/GenStudy", exist_ok=True)
    outName = "../plots/GenStudy/GsfTrack_choice.pdf"
    c1.Print(outName)
    c1.Close()
    

if __name__ == "__main__" :
    start_time = time.time()
    
    # PyROOT does not display any graphics(root "-b" option)
    ROOT.gROOT.SetBatch()
    
    infile_list = glob("../miniTree/*/miniTree_HDalitz_ggF_eeg_*.root")
    ROOT.EnableImplicitMT(20)
    rdf = ROOT.RDataFrame("miniTree", infile_list).Define("weight", "mcwei * genwei")
    arr = rdf.Filter("category == 2").AsNumpy(columns = ["meeRatio", "meeRatioPtMax"])

    xmin_dr, xmax_dr, sigma_dr = sigmaEff(np.array(arr["meeRatio"]), threshold = 0.683)
    print("smallest dR: eff sigma = {}".format(sigma_dr))

    xmin_pt, xmax_pt, sigma_pt = sigmaEff(np.array(arr["meeRatioPtMax"]), threshold = 0.683)
    print("highest pT: eff sigma = {}".format(sigma_pt))
    
    main()
    
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))