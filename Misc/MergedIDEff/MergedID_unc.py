import ROOT
import os
import sys
import numpy as np
import pandas as pd
from glob import glob
from argparse import ArgumentParser
from plugins.colorPrint import *
from plugins.CMS_lumi import CMS_lumi


def find_files(indir):
    Infile_list = sorted(glob(indir))
    f = ROOT.std.vector("string")()
    for i in Infile_list:
        f.push_back(i)
    return f


def Draw2DSFs(df_EB, df_EE):
    hSF_2D = ROOT.TH2D("hSF_2D", "", 3, np.array([0, 1.4442, 1.566, 2.5]), len(binArray)-1, binArray)
    for i in range(3):
        if i == 0:
            df = df_EB
        else:
            df = df_EE
        for j in range(len(binArray)-1):
            if i == 1:
                hSF_2D.SetBinContent(i+1, j+1, 0)
                hSF_2D.SetBinError(i+1, j+1, 0)
            else:
                hSF_2D.SetBinContent(i+1, j+1, df["nominal"][j])
                hSF_2D.SetBinError(i+1, j+1, df["tot_unc"][j])

    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    c1 = ROOT.TCanvas("c1", "", 700, 700)
    c1.cd()
    ROOT.gStyle.SetPaintTextFormat("1.3f")
    c1.SetLogy()
    c1.SetRightMargin(0.15)
    c1.SetBottomMargin(0.12)
    c1.SetTopMargin(0.06)
    c1.SetLeftMargin(0.15)

    hSF_2D.GetXaxis().SetTitle("|#eta^{SC}|")
    hSF_2D.GetXaxis().SetMoreLogLabels()
    hSF_2D.GetXaxis().SetTickSize(0.02)
    hSF_2D.GetXaxis().SetTitleSize(0.05)
    hSF_2D.GetXaxis().SetLabelSize(0.05)
    hSF_2D.GetXaxis().SetTitleOffset(1.1)

    hSF_2D.GetYaxis().SetTitle("p^{#gamma}_{T} [GeV]")
    hSF_2D.GetYaxis().SetMoreLogLabels()
    hSF_2D.GetYaxis().SetNdivisions(302)
    hSF_2D.GetYaxis().SetTickSize(0.02)
    hSF_2D.GetYaxis().SetTitleSize(0.05)
    hSF_2D.GetYaxis().SetLabelSize(0.05)
    hSF_2D.GetYaxis().SetTitleOffset(1.4)

    hSF_2D.GetZaxis().SetLabelSize(0.05)
    hSF_2D.GetZaxis().SetRangeUser(0.6, 1.2)
    hSF_2D.SetMarkerSize(2)
    hSF_2D.Draw("COLZ texte")
    c1.RedrawAxis()

    CMS_lumi(c1, 5, 0, "137.2 fb^{-1}", 2017, True, "Work-in-progress", "", "")

    directory = "./plots/2DSFs"
    if not os.path.exists(directory):
        os.makedirs(directory)

    c1.Print("{}/2DSFs_combined.png".format(directory))
    c1.Print("{}/2DSFs_combined.pdf".format(directory))
    c1.Close()


def Draw2Dunc(df_EB, df_EE, Type):
    hunc_2D = ROOT.TH2D("hunc_2D", "", 3, np.array([0, 1.4442, 1.566, 2.5]), len(binArray)-1, binArray)
    for i in range(3):
        if i == 0:
            df = df_EB
        else:
            df = df_EE
        for j in range(len(binArray)-1):
            if i == 1:
                hunc_2D.SetBinContent(i+1, j+1, 0)
            else:
                hunc_2D.SetBinContent(i+1, j+1, df[Type][j])

    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    c1 = ROOT.TCanvas("c1", "", 700, 700)
    c1.cd()
    ROOT.gStyle.SetPaintTextFormat("1.3f")
    c1.SetLogy()
    c1.SetRightMargin(0.15)
    c1.SetBottomMargin(0.12)
    c1.SetTopMargin(0.06)
    c1.SetLeftMargin(0.15)

    hunc_2D.GetXaxis().SetTitle("|#eta^{SC}|")
    hunc_2D.GetXaxis().SetMoreLogLabels()
    hunc_2D.GetXaxis().SetTickSize(0.02)
    hunc_2D.GetXaxis().SetTitleSize(0.05)
    hunc_2D.GetXaxis().SetLabelSize(0.05)
    hunc_2D.GetXaxis().SetTitleOffset(1.1)

    hunc_2D.GetYaxis().SetTitle("p^{#gamma}_{T} [GeV]")
    hunc_2D.GetYaxis().SetMoreLogLabels()
    hunc_2D.GetYaxis().SetNdivisions(302)
    hunc_2D.GetYaxis().SetTickSize(0.02)
    hunc_2D.GetYaxis().SetTitleSize(0.05)
    hunc_2D.GetYaxis().SetLabelSize(0.05)
    hunc_2D.GetYaxis().SetTitleOffset(1.4)

    hunc_2D.GetZaxis().SetLabelSize(0.05)
    hunc_2D.GetZaxis().SetRangeUser(0, 0.1)
    hunc_2D.SetMarkerSize(2)
    hunc_2D.Draw("COLZ text")
    c1.RedrawAxis()

    CMS_lumi(c1, 5, 0, "137.2 fb^{-1}", 2017, True, "Work-in-progress", "", "")

    directory = "./plots/2Dunc"
    if not os.path.exists(directory):
        os.makedirs(directory)

    c1.Print("{}/2Dunc_{}_combined.png".format(directory, Type))
    c1.Print("{}/2Dunc_{}_combined.pdf".format(directory, Type))
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


def get_sfs(varied, region):
    # selections
    reg = "isEBPho_lep1 == 1" if region == "EB" else "isEEPho_lep1 == 1"
    rconv = 40.
    if varied == "rconv_up":
        rconv = 30.
    if varied == "rconv_do":
        rconv = 50.
    sel_den = "&&".join([
        "convVtxRadius_lep1 < {}".format(rconv),
        "isHggPho_lep1 == 1",
        "nGsfMatchToReco_lep1 > 1",
        reg
    ])
    sel_num = sel_den + "&& eleClass_lep1 == 0"

    rdf_MC = rdf_Zg
    if (varied == "bkg"):
        rdf_MC = rdf_bkg
    if (varied == "alt"):
        rdf_MC = rdf_DY

    if (varied == "pu_up"):
        rdf_MC = rdf_MC.Define("wei", "puwei_up * mcwei")
    elif (varied == "pu_do"):
        rdf_MC = rdf_MC.Define("wei", "puwei_down * mcwei")
    else:
        rdf_MC = rdf_MC.Define("wei", "puwei * mcwei")

    h_num_MC = rdf_MC.Filter(sel_num).Histo1D(("h_num_MC", " ", len(binArray)-1, binArray), "phoCalibEt_lep1", "wei").GetPtr()
    h_den_MC = rdf_MC.Filter(sel_den).Histo1D(("h_den_MC", " ", len(binArray)-1, binArray), "phoCalibEt_lep1", "wei").GetPtr()
    eff_MC = ROOT.TEfficiency(h_num_MC, h_den_MC)
    eff_MC.SetStatisticOption(ROOT.TEfficiency.kBUniform)
    eff_MC.SetConfidenceLevel(0.683)
    eff_MC.SetPosteriorMode(1)

    h_num_Da = rdf_Da.Filter(sel_num).Histo1D(("h_num_Da", " ", len(binArray)-1, binArray), "phoCalibEt_lep1").GetPtr()
    h_den_Da = rdf_Da.Filter(sel_den).Histo1D(("h_den_Da", " ", len(binArray)-1, binArray), "phoCalibEt_lep1").GetPtr()
    eff_Da = ROOT.TEfficiency(h_num_Da, h_den_Da)
    eff_Da.SetStatisticOption(ROOT.TEfficiency.kBUniform)
    eff_Da.SetConfidenceLevel(0.683)
    eff_Da.SetPosteriorMode(1)

    heff_Da = eff_Da.GetCopyPassedHisto()
    heff_Da.Divide(eff_Da.GetCopyTotalHisto())
    heff_MC = eff_MC.GetCopyPassedHisto()
    heff_MC.Divide(eff_MC.GetCopyTotalHisto())
    hSF = heff_Da.Clone()
    hSF.Divide(heff_MC)

    sf, sf_err = [], []
    SFNbins = hSF.GetNbinsX()
    for i in range(SFNbins):
        Da_error = np.sqrt(0.5 * (np.power(eff_Da.GetEfficiencyErrorUp(i+1), 2) + np.power(eff_Da.GetEfficiencyErrorLow(i + 1), 2)))
        MC_error = np.sqrt(0.5 * (np.power(eff_MC.GetEfficiencyErrorUp(i+1), 2) + np.power(eff_MC.GetEfficiencyErrorLow(i + 1), 2)))
        sf.append(hSF.GetBinContent(i + 1))
        sf_err.append(ErrProp(hSF.GetBinContent(i + 1), Da_error, heff_Da.GetBinContent(i + 1), MC_error, heff_MC.GetBinContent(i + 1)))

    return sf, sf_err


if __name__ == "__main__":

    fData = ROOT.std.vector("string")((
        "./miniTree/2016_preVFP/miniTree_Data_2016_preVFP.root",
        "./miniTree/2016_postVFP/miniTree_Data_2016_postVFP.root",
        "./miniTree/2017/miniTree_Data_2017.root",
        "./miniTree/2018/miniTree_Data_2018.root"
    ))
    fZg = find_files("./miniTree/*/miniTree_ZGToLLG_*.root")
    fbkg = find_files("./miniTree/*/miniTree_TTJets_*.root")
    fbkg.insert(fbkg.end(), fZg.begin(), fZg.end())
    fDY = find_files("./miniTree/*/miniTree_DYJets_*.root")

    ROOT.EnableImplicitMT()
    rdf_Da = ROOT.RDataFrame("outTree", fData)
    rdf_Zg = ROOT.RDataFrame("outTree", fZg)
    rdf_bkg = ROOT.RDataFrame("outTree", fbkg)
    rdf_DY = ROOT.RDataFrame("outTree", fDY)

    binArray = np.array([25., 31., 55, 150.])
    sysType = ["nominal", "rconv_up", "rconv_do", "pu_up", "pu_do", "bkg", "alt"]
    df_list = []
    for region in ["EB", "EE"]:
        SFDict = {}
        for varied in sysType:
            SF, SFerr = get_sfs(varied, region)
            SFDict[varied] = SF
            if varied == "nominal":
                SFDict["nominal_staunc"] = SFerr

        df = pd.DataFrame.from_dict(SFDict)
        df["rconv_difup"] = (df["rconv_up"] - df["nominal"]).abs() / df["nominal"]
        df["rconv_difdo"] = (df["rconv_do"] - df["nominal"]).abs() / df["nominal"]
        df["rconv_sysunc"] = df[["rconv_difup", "rconv_difdo"]].max(axis = 1)

        df["pu_difup"] = (df["pu_up"] - df["nominal"]).abs() / df["nominal"]
        df["pu_difdo"] = (df["pu_do"] - df["nominal"]).abs() / df["nominal"]
        df["pu_sysunc"] = df[["pu_difup", "pu_difdo"]].max(axis = 1)

        df["bkg_sysunc"] = (df["bkg"] - df["nominal"]).abs() / df["nominal"]
        df["alt_sysunc"] = (df["alt"] - df["nominal"]).abs() / df["nominal"]

        df["tot_unc"] = np.sqrt(df["nominal_staunc"]**2 + df["rconv_sysunc"]**2 + df["pu_sysunc"]**2 + df["bkg_sysunc"]**2 + df["alt_sysunc"]**2)

        df_list.append(df)

    Draw2DSFs(df_list[0], df_list[1])

    for t in ["nominal_staunc", "rconv_sysunc", "pu_sysunc", "bkg_sysunc", "alt_sysunc"]:
        Draw2Dunc(df_list[0], df_list[1], t)
