import ROOT
import sys, os
import json
import numpy as np
from plugins.CMS_lumi import CMS_lumi

# script to check the signal and zg have similar distributions among the training features


def Draw1DHist(c, vh, vc=[None, None], vcl=["#E69D45", "#0061a8"], xaxis="x-axis", yaxis="x-axis", option="hist", y_axisscale=1.3, Log=False):
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.07)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetBottomMargin(0.15)
    if (Log == True):
        c.SetLogy()

    # Set the axis style
    if (Log == True):
        ymax = vh[0].GetBinContent(vh[0].GetMaximumBin()) * y_axisscale
        ymin = 1E-3
    else:
        ymax = vh[0].GetBinContent(vh[0].GetMaximumBin()) * y_axisscale
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


def main(pT_do, pT_up, region):
    # selections of zg
    sel_Zg = "&&".join([
        "convVtxRadius_lep1 < {}".format(rconv),
        "phoCalibEt_lep1 < {} && phoCalibEt_lep1 > {}".format(pT_up, pT_do),
        "isHggPho_lep1 == 1",
        "isEBPho_lep1 == 1"
    ])
    rdf_Zg_tmp = rdf_Zg.Filter(sel_Zg).Define("weight", "mcwei * puwei")

    # selections of signal
    sel_Si = "&&".join([
        "eleCalibPt_lep1 < {} && eleCalibPt_lep1 > {}".format(pT_up, pT_do),
        "elePresel_lep1 == 1",
        "abs(eleSCEta_lep1) < 1.4442",
        "category == 2"
    ])
    rdf_Si_tmp = rdf_Si.Filter(sel_Si).Define("weight", "mcwei * puwei")

    # Draw the distributions
    with open("FeaturePlotParam_MergedIDEB2Gsf.json", "r") as fp:
        myTags = json.load(fp)

        for tagName, par in myTags.items():
            setlogy = True if par["logy"] == "true" else False

            bininfo = par["bininfo"]
            h1 = rdf_Zg_tmp.Histo1D(("h1_{}".format(tagName), " ", bininfo[0], bininfo[1], bininfo[2]), "{}_lep1".format(tagName), "weight").GetPtr()
            h2 = rdf_Si_tmp.Histo1D(("h2_{}".format(tagName), " ", bininfo[0], bininfo[1], bininfo[2]), "{}_lep1".format(tagName), "weight").GetPtr()
            h1.Scale(1./h1.Integral())
            h2.Scale(1./h2.Integral())

            ROOT.gStyle.SetPadTickX(1)
            ROOT.gStyle.SetPadTickY(1)
            ROOT.gStyle.SetOptStat(0)
            c1 = ROOT.TCanvas("c1", "", 700, 600)
            c1.cd()

            Draw1DHist(c1, [h1, h2], xaxis=par["xAxisLabel"], yaxis="arb. unit", y_axisscale=par["y_axisscale"] , Log=setlogy)
            CMS_lumi(c1, 4, 0, "", 2017, True, "Work-in-progress", "", "")

            leg = ROOT.TLegend(0.2, 0.75, 0.55, 0.9)
            leg.SetTextSize(0.04)
            leg.AddEntry(h1, "Z#gamma with r_{conv} < %.1f cm" %rconv, "l")
            leg.AddEntry(h2, "H #rightarrow #gamma*#gamma #rightarrow ee#gamma", "l")
            leg.SetFillColor(0)
            leg.SetLineColor(0)
            leg.Draw()

            ltx1 = ROOT.TLatex()
            ltx1.SetNDC()
            ltx1.SetTextFont(42)
            ltx1.SetTextSize(0.04)
            ltx1.DrawLatex(0.65, 0.8, "P_{T}: %d to %d GeV" %(pT_do, pT_up))

            c1.RedrawAxis()

            outRes = "./plots/{}".format(region)
            if not os.path.exists(outRes):
                os.makedirs(outRes)
            c1.Print("{}/FeatureComparison_{}_{}.pdf".format(outRes, tagName, region))
            c1.Print("{}/FeatureComparison_{}_{}.png".format(outRes, tagName, region))
            c1.Close()


if __name__ == "__main__":
    rdf_Zg = ROOT.RDataFrame("outTree", "./miniTree/2017/miniTree_ZGToLLG_2017.root")
    rdf_Si = ROOT.RDataFrame("outTree", "../GENstudy/minitree/2017/Minitree_HDalitz_*_m125_*.root")

    rconv = 16
    pTbin = [15., 25., 40., 65, 150.]
    reg = ["pT15to25", "pT25to40", "pT40to65", "pT65to150"]
    for i in range(len(pTbin) - 1):
        main(pTbin[i], pTbin[i+1], reg[i])
        print(" ")