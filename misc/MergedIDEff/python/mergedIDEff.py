import ROOT
import os
import sys
import numpy as np
from argparse import ArgumentParser
from CMS_lumi import CMS_lumi
from colorPrint import *


def get_parser():
    parser = ArgumentParser(description="python script to make 1D efficiency plots")
    parser.add_argument("-e", "--era", help="era to run", default="combined", type=str)
    return parser



def Make1DEff(binArray, eff_Da, eff_MC, SFerr, option, outName):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetOptStat(0)

    c1 = ROOT.TCanvas("c1", "", 800, 800)
    c1.cd()

    #===============================================#
    #    1st Pad's Setting(Draw the Efficiency)     #
    #===============================================#
    pad1 = ROOT.TPad("pad1", " ", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.05)
    pad1.SetTopMargin(0.09)
    pad1.SetRightMargin(0.14)
    pad1.SetLeftMargin(0.12)
    pad1.SetBottomMargin(0.03)
    pad1.Draw()
    pad1.cd()

    # ("global title; x-axis title; y-axis title")
    eff_MC.SetTitle(";;Efficiency")
    eff_MC.SetMarkerColor(ROOT.TColor.GetColor("#C72C41"))
    eff_MC.SetMarkerSize(1.4)
    eff_MC.SetMarkerStyle(20)
    eff_MC.SetLineColor(ROOT.TColor.GetColor("#C72C41"))
    eff_MC.SetLineWidth(2)
    eff_MC.Draw("AP")

    eff_Da.SetMarkerColor(ROOT.TColor.GetColor("#202020"))
    eff_Da.SetMarkerSize(1.4)
    eff_Da.SetMarkerStyle(20)
    eff_Da.SetLineColor(ROOT.TColor.GetColor("#202020"))
    eff_Da.SetLineWidth(2)
    eff_Da.Draw("P SAME")

    pad1.Update()
    efferr = eff_MC.GetPaintedGraph()
    efferr.GetYaxis().SetRangeUser(0., 1.2)
    efferr.GetYaxis().SetTickSize(0.03)
    efferr.GetYaxis().SetTitleSize(0.06)
    efferr.GetYaxis().SetLabelSize(0.06)
    efferr.GetYaxis().SetTitleOffset(0.92)
    font = efferr.GetYaxis().GetLabelFont()

    efferr.GetXaxis().SetRangeUser(binArray[0], binArray[-1])
    efferr.GetXaxis().SetTickSize(0.03)
    efferr.GetXaxis().SetTitleSize(0.06)
    efferr.GetXaxis().SetLabelSize(0.06)
    efferr.GetXaxis().SetLabelOffset(0.1)
    efferr.GetXaxis().SetTitleOffset(1)

    pad1.Update()
    scal = 4 if option == "EB" else 6
    htotal_Da = eff_Da.GetCopyTotalHisto()
    rightmax = htotal_Da.GetBinContent(htotal_Da.GetMaximumBin()) * scal
    scale = 1.2/rightmax
    htotal_Da.SetLineWidth(2)
    htotal_Da.SetLineColor(ROOT.TColor.GetColor("#202020"))
    htotal_Da.Scale(scale)
    htotal_Da.Draw("hist same")

    hpass_Da = eff_Da.GetCopyPassedHisto()
    hpass_Da.SetLineStyle(7)
    hpass_Da.SetLineWidth(2)
    hpass_Da.SetLineColor(ROOT.TColor.GetColor("#202020"))
    hpass_Da.Scale(scale)
    hpass_Da.Draw("hist same")

    htotal_MC = eff_MC.GetCopyTotalHisto()
    htotal_MC.SetLineWidth(2)
    htotal_MC.SetLineColor(ROOT.TColor.GetColor("#C72C41"))
    htotal_MC.Scale(scale)
    htotal_MC.Draw("hist same")

    hpass_MC = eff_MC.GetCopyPassedHisto()
    hpass_MC.SetLineStyle(7)
    hpass_MC.SetLineWidth(2)
    hpass_MC.SetLineColor(ROOT.TColor.GetColor("#C72C41"))
    hpass_MC.Scale(scale)
    hpass_MC.Draw("hist same")

    axis = ROOT.TGaxis(pad1.GetUxmax(), 0, pad1.GetUxmax(), 1.2, 0, rightmax, 510, "+L")
    axis.SetLineColor(ROOT.TColor.GetColor("#202020"))
    axis.SetLabelColor(ROOT.TColor.GetColor("#202020"))
    axis.SetTitle("Events")
    axis.SetTitleFont(font)
    axis.SetTitleSize(0.057)
    axis.SetTitleOffset(1.15)
    axis.SetLabelFont(font)
    axis.SetLabelSize(0.057)
    axis.SetTickSize(0.03)
    axis.Draw()

    pad1.Update()

    ltx = ROOT.TLatex()
    ltx.SetNDC()
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.06)
    ltx.DrawLatex(0.42, 0.34, "Z #rightarrow #mu#mu#gamma, r_{conv}< 40cm")

    ext = "ECAL Barrel" if option == "EB" else "ECAL Endcap"
    ltx1 = ROOT.TLatex()
    ltx1.SetNDC()
    ltx1.SetTextFont(42)
    ltx1.SetTextSize(0.06)
    ltx1.DrawLatex(0.17, 0.8, ext)

    leg = ROOT.TLegend(0.4, 0.14, 0.77, 0.3)
    leg.SetTextSize(0.06)
    leg.AddEntry(eff_MC.GetPaintedGraph(), "Z#gamma", "LE1P")
    leg.AddEntry(eff_Da.GetPaintedGraph(), "Data", "LE1P")
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.Draw()

    if args.era == "combined":
        CMS_lumi(pad1, 5, 0, "%d fb^{-1}" % (Luminosity[args.era]), 2017, True, "Preliminary", "", "")
    else:
        CMS_lumi(pad1, 5, 0, "%.1f fb^{-1}" % (Luminosity[args.era]), 2017, True, "Preliminary", "", "")
    c1.cd()

    #===============================================#
    #   2nd Pad's Setting(Draw the Scale Factors)   #
    #===============================================#
    pad2 = ROOT.TPad("pad2", "", 0, 0, 1, 0.3)
    pad2.SetGridy()

    pad2.SetRightMargin(0.14)
    pad2.SetLeftMargin(0.12)
    pad2.SetTopMargin(0.06)
    pad2.SetBottomMargin(0.4)

    pad2.Draw()
    pad2.cd()

    SFerr.SetName("")
    SFerr.SetTitle("")
    SFerr.GetXaxis().SetTitle("p^{#gamma}_{T} [GeV]")
    SFerr.GetYaxis().SetTitle("Scale Factor")
    SFerr.GetYaxis().SetRangeUser(0.7, 1.3)

    SFerr.SetMarkerColor(ROOT.TColor.GetColor("#202020"))
    SFerr.SetMarkerSize(1.4)
    SFerr.SetMarkerStyle(20)
    SFerr.SetLineColor(ROOT.TColor.GetColor("#202020"))
    SFerr.SetLineWidth(2)

    SFerr.GetXaxis().SetRangeUser(binArray[0], binArray[-1])
    SFerr.GetXaxis().SetTickSize(0.03 * (7 / 3.))
    SFerr.GetXaxis().SetTitleSize(0.16)
    SFerr.GetXaxis().SetTitleOffset(1.15)
    SFerr.GetXaxis().SetLabelSize(0.06 * (7 / 3.))
    SFerr.GetXaxis().SetLabelOffset(0.05)
    SFerr.GetYaxis().SetTitleSize(0.13)
    SFerr.GetYaxis().SetTitleOffset(0.18 * (7 / 3.))
    SFerr.GetYaxis().SetLabelSize(0.06 * (7 / 3.))
    SFerr.GetYaxis().SetNdivisions(502)
    SFerr.Draw("AP")

    c1.Print("{}.pdf".format(outName))
    c1.Print("{}.png".format(outName))
    c1.Close()


# Error Propagation for Scale Factor (SF = Eff_A / Eff_B)
# * Some notes:
# *  1. When the quantity depending on other measurements, the error
# *     of this quantity comes from other measurements.
# *  2. In other words, the uncertainties of other measurements propagate to
# *     the uncertainty of this quantity.
# *  3. For example, the uncertainty of the scale factor comes from
# *     the uncertainty of the effs. of data and MC.
def ErrProp(f, sigmaA, A, sigmaB, B):
    error = f * np.sqrt((sigmaA / A) * (sigmaA / A) + (sigmaB / B) * (sigmaB / B))
    return error


def main(option):
    # selections
    region = "isEBPho_Lead == 1" if option == "EB" else "isEEPho_Lead == 1"
    rconv = 40

    # merged ID wp
    M2EBWP = 0.4263
    M2EEWP = 0.4205
    # M2EBWP = 0
    # M2EEWP = 0
    # M1EBWP = 0.4176
    # M1EEWP = 0.4663

    sel_den = "&&".join([
        "phoCalibEt_Lead > 25",
        "fsrmu.Pt() > 10.",
        "convVtxRadius_Lead < {}".format(rconv),
        "isHggPho_Lead == 1",
        "nGsfMatchToReco_Lead >= 2",
        region
    ])

    if option == "EB":
        sel_num = sel_den + "&& eleClass_Lead == 0 && eleXGBID_Lead > {}".format(M2EBWP)
    else:
        sel_num = sel_den + "&& eleClass_Lead == 0 && eleXGBID_Lead > {}".format(M2EEWP)

    binArray = np.asarray([25, 35, 50, 150], "d")
    print("Bin width used in efficiency calculation: {}".format(binArray))
    print("photon with radius < {}cm are used.".format(rconv))

    #===============================================#
    #               Efficiency for MC               #
    #===============================================#
    # Do not apply genweight!!!! (negative weight may makes passed > total)
    print(color.GREEN + "Efficiency calculation in MC {}".format(option) + color.END)
    h_num_MC = rdf_MC.Filter(sel_num).Histo1D(("h_num_MC", " ", len(binArray) - 1, binArray), "phoCalibEt_Lead", "wei").GetPtr()
    h_den_MC = rdf_MC.Filter(sel_den).Histo1D(("h_den_MC", " ", len(binArray) - 1, binArray), "phoCalibEt_Lead", "wei").GetPtr()
    print("[INFO] Passed entries: {}".format([round(h_num_MC.GetBinContent(i + 1), 4) for i in range(len(binArray) - 1)]))
    print("[INFO] Total entries: {}".format([round(h_den_MC.GetBinContent(i + 1), 4) for i in range(len(binArray) - 1)]))

    efferr_MC_up, efferr_MC_do = [], []
    if (ROOT.TEfficiency.CheckConsistency(h_num_MC, h_den_MC)):
        eff_MC = ROOT.TEfficiency(h_num_MC, h_den_MC)
        eff_MC.SetStatisticOption(ROOT.TEfficiency.kBUniform)
        eff_MC.SetConfidenceLevel(0.683)
        eff_MC.SetPosteriorMode(1)
        print("Bayesian confidence interval with uniform prior beta(1, 1) and posterior mode is used")
        for i in range(len(binArray) - 1):
            efferr_MC_up.append(eff_MC.GetEfficiencyErrorUp(i + 1))
            efferr_MC_do.append(eff_MC.GetEfficiencyErrorLow(i + 1))
        print("[INFO] Efficiency: {}".format([round(eff_MC.GetEfficiency(i + 1), 4) for i in range(len(binArray) - 1)]))
        print("[INFO] Efficiency err up: {}".format(np.round(efferr_MC_up, 4)))
        print("[INFO] Efficiency err do: {}".format(np.round(efferr_MC_do, 4)))

    #===============================================#
    #             Efficiency for Data               #
    #===============================================#
    print(color.GREEN + "Efficiency calculation in Data {}".format(option) + color.END)
    h_num_Da = rdf_Da.Filter(sel_num).Histo1D(("h_num_Da", " ", len(binArray) - 1, binArray), "phoCalibEt_Lead").GetPtr()
    h_den_Da = rdf_Da.Filter(sel_den).Histo1D(("h_den_Da", " ", len(binArray) - 1, binArray), "phoCalibEt_Lead").GetPtr()
    print("[INFO] Passed entries: {}".format([h_num_Da.GetBinContent(i + 1) for i in range(len(binArray) - 1)]))
    print("[INFO] Total entries: {}".format([h_den_Da.GetBinContent(i + 1) for i in range(len(binArray) - 1)]))

    efferr_Da_up, efferr_Da_do = [], []
    if (ROOT.TEfficiency.CheckConsistency(h_num_Da, h_den_Da)):
        eff_Da = ROOT.TEfficiency(h_num_Da, h_den_Da)
        eff_Da.SetStatisticOption(ROOT.TEfficiency.kBUniform)
        eff_Da.SetConfidenceLevel(0.683)
        eff_Da.SetPosteriorMode(1)
        print("Bayesian confidence interval with uniform prior beta(1, 1) and posterior mode is used")
        for i in range(len(binArray) - 1):
            efferr_Da_up.append(eff_Da.GetEfficiencyErrorUp(i + 1))
            efferr_Da_do.append(eff_Da.GetEfficiencyErrorLow(i + 1))

        print("[INFO] Efficiency: {}".format([round(eff_Da.GetEfficiency(i + 1), 4) for i in range(len(binArray) - 1)]))
        print("[INFO] Efficiency err up: {}".format(np.round(efferr_Da_up, 4)))
        print("[INFO] Efficiency err do: {}".format(np.round(efferr_Da_do, 4)))

    #===============================================#
    #                  SF calculation               #
    #===============================================#
    # https://root.cern.ch/doc/master/TGraphAsymmErrors_8cxx_source.html#l01129
    # get efficiency histogram of data

    print(color.GREEN + "Scale factors calculation in {}".format(option) + color.END)
    heff_Da = eff_Da.GetCopyPassedHisto()
    heff_Da.Divide(eff_Da.GetCopyTotalHisto())

    # get efficiency histogram of MC
    heff_MC = eff_MC.GetCopyPassedHisto()
    heff_MC.Divide(eff_MC.GetCopyTotalHisto())

    # get SF histogram
    hSF = heff_Da.Clone()
    hSF.Divide(heff_MC)

    x, y, xerr, yerr = [], [], [], []
    SFNbins = hSF.GetNbinsX()
    for i in range(SFNbins):
        Da_error = np.sqrt(0.5 * (np.power(efferr_Da_up[i], 2) + np.power(efferr_Da_do[i], 2)))
        MC_error = np.sqrt(0.5 * (np.power(efferr_MC_up[i], 2) + np.power(efferr_MC_do[i], 2)))
        x.append(hSF.GetBinCenter(i + 1))
        y.append(hSF.GetBinContent(i + 1))
        xerr.append(abs(hSF.GetBinCenter(i + 1) - hSF.GetBinLowEdge(i + 1)))
        yerr.append(ErrProp(hSF.GetBinContent(i + 1), Da_error, heff_Da.GetBinContent(i + 1), MC_error, heff_MC.GetBinContent(i + 1)))

    SFerr = ROOT.TGraphErrors(SFNbins, np.array(x), np.array(y), np.array(xerr), np.array(yerr))
    print("[INFO] Scale factors: {}".format(np.round(y, 4)))
    print("[INFO] Scale factors err: {}".format(np.round(yerr, 4)))

    # Draw the results
    directory = "../plots/1DEff/{}".format(args.era)
    if not os.path.exists(directory):
        os.makedirs(directory)
    Make1DEff(binArray, eff_Da, eff_MC, SFerr, option, "{}/MergedID_1DEff_{}_{}".format(directory, args.era, option))


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True) # turn off GUI

    parser = get_parser()
    args = parser.parse_args()

    fData = {
        "UL2016preVFP": "../miniTree/UL2016preVFP/miniTree_Data_UL2016preVFP.root",
        "UL2016postVFP": "../miniTree/UL2016postVFP/miniTree_Data_UL2016postVFP.root",
        "UL2016":[
            "../miniTree/UL2016preVFP/miniTree_Data_UL2016preVFP.root",
            "../miniTree/UL2016postVFP/miniTree_Data_UL2016postVFP.root"
        ],
        "UL2017": "../miniTree/UL2017/miniTree_Data_UL2017.root",
        "UL2018": "../miniTree/UL2018/miniTree_Data_UL2018.root"
    }

    fMC = {
        "UL2016preVFP": "../miniTree/UL2016preVFP/miniTree_ZGToLLG_UL2016preVFP.root",
        "UL2016postVFP": "../miniTree/UL2016postVFP/miniTree_ZGToLLG_UL2016postVFP.root",
        "UL2016": [
            "../miniTree/UL2016preVFP/miniTree_ZGToLLG_UL2016preVFP.root",
            "../miniTree/UL2016postVFP/miniTree_ZGToLLG_UL2016postVFP.root"
        ],
        "UL2017": "../miniTree/UL2017/miniTree_ZGToLLG_UL2017.root",
        "UL2018": "../miniTree/UL2018/miniTree_ZGToLLG_UL2018.root"
    }

    Luminosity = {
        "UL2016preVFP": 19.5,
        "UL2016postVFP": 16.8,
        "2016": 36.3,
        "UL2017": 41.5,
        "UL2018": 59.8,
        "combined": 138
    }

    ROOT.EnableImplicitMT(20)
    if (args.era != "combined"):
        rdf_MC = ROOT.RDataFrame("miniTree", fMC[args.era])
        rdf_Da = ROOT.RDataFrame("miniTree", fData[args.era])
    else:
        del fMC["UL2016"], fData["UL2016"]
        rdf_MC = ROOT.RDataFrame("miniTree", list(fMC.values()))
        rdf_Da = ROOT.RDataFrame("miniTree", list(fData.values()))

    for i in ["EB", "EE"]:
        main(i)
        print("")
