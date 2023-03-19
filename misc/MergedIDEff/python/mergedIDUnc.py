import ROOT
import os
import sys
import numpy as np
import pandas as pd
from colorPrint import *
from CMS_lumi import CMS_lumi


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

    CMS_lumi(c1, 5, 0, "138 fb^{-1}", 2017, True, "Preliminary", "", "")

    directory = "../plots/2DSFs"
    if not os.path.exists(directory):
        os.makedirs(directory)

    if not os.path.exists("../data"):
        os.makedirs("../data")
    fout = ROOT.TFile("../data/MergedIDSFs_combined.root", "RECREATE")
    fout.cd()
    hSF_2D.Write()
    fout.Close()

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

    CMS_lumi(c1, 5, 0, "138 fb^{-1}", 2017, True, "Preliminary", "", "")

    directory = "../plots/2Dunc"
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
    reg = "isEBPho_Lead == 1" if region == "EB" else "isEEPho_Lead == 1"
    rconv = 40.
    if varied == "rconv_up":
        rconv = 30.
    if varied == "rconv_do":
        rconv = 50.
    sel_den = "&&".join([
        "fsrmu.Pt() > 10.",
        "convVtxRadius_Lead < {}".format(rconv),
        "isHggPho_Lead == 1",
        "nGsfMatchToReco_Lead > 1",
        reg
    ])

    # merged ID wp
    M2EBWP = 0.4469
    M2EEWP = 0.4552
    # M1EBWP = 0.4451
    # M1EEWP = 0.4188
    if region == "EB":
        sel_num = sel_den + "&& eleClass_Lead == 0 && eleXGBID_Lead > {}".format(M2EBWP)
    else:
        sel_num = sel_den + "&& eleClass_Lead == 0 && eleXGBID_Lead > {}".format(M2EEWP)

    rdf_MC = rdf_Zg
    if (varied == "bkg"):
        rdf_MC = rdf_bkg
    if (varied == "alt"):
        rdf_MC = rdf_DY

    if (varied == "pu_up"):
        rdf_MC = rdf_MC.Define("weight", "puwei_up * mcwei")
    elif (varied == "pu_do"):
        rdf_MC = rdf_MC.Define("weight", "puwei_down * mcwei")
    else:
        rdf_MC = rdf_MC.Define("weight", "puwei * mcwei")

    h_num_MC = rdf_MC.Filter(sel_num).Histo1D(("h_num_MC", " ", len(binArray)-1, binArray), "phoCalibEt_Lead", "weight").GetPtr()
    h_den_MC = rdf_MC.Filter(sel_den).Histo1D(("h_den_MC", " ", len(binArray)-1, binArray), "phoCalibEt_Lead", "weight").GetPtr()
    eff_MC = ROOT.TEfficiency(h_num_MC, h_den_MC)
    eff_MC.SetStatisticOption(ROOT.TEfficiency.kBUniform)
    eff_MC.SetConfidenceLevel(0.683)
    eff_MC.SetPosteriorMode(1)

    h_num_Da = rdf_Da.Filter(sel_num).Histo1D(("h_num_Da", " ", len(binArray)-1, binArray), "phoCalibEt_Lead").GetPtr()
    h_den_Da = rdf_Da.Filter(sel_den).Histo1D(("h_den_Da", " ", len(binArray)-1, binArray), "phoCalibEt_Lead").GetPtr()
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

    eff_da, eff_da_err = [], []
    eff_mc, eff_mc_err = [], []
    sf, sf_err = [], []
    SFNbins = hSF.GetNbinsX()
    for i in range(SFNbins):
        Da       = eff_Da.GetEfficiency(i+1)
        MC       = eff_MC.GetEfficiency(i+1)
        Da_error = np.sqrt(0.5 * (np.power(eff_Da.GetEfficiencyErrorUp(i+1), 2) + np.power(eff_Da.GetEfficiencyErrorLow(i + 1), 2)))
        MC_error = np.sqrt(0.5 * (np.power(eff_MC.GetEfficiencyErrorUp(i+1), 2) + np.power(eff_MC.GetEfficiencyErrorLow(i + 1), 2)))

        eff_da.append(Da)
        eff_mc.append(MC)
        eff_da_err.append(Da_error)
        eff_mc_err.append(MC_error)
        sf.append(hSF.GetBinContent(i + 1))
        sf_err.append(ErrProp(hSF.GetBinContent(i + 1), Da_error, heff_Da.GetBinContent(i + 1), MC_error, heff_MC.GetBinContent(i + 1)))

    print("[INFO] Eff. data        : {}".format(np.round(eff_da, 4)))
    print("[INFO] Eff. data err    : {}".format(np.round(eff_da_err, 4)))
    print("[INFO] Eff. mc          : {}".format(np.round(eff_mc, 4)))
    print("[INFO] Eff. mc err      : {}".format(np.round(eff_mc_err, 4)))
    print("[INFO] Scale factors    : {}".format(np.round(sf, 4)))
    print("[INFO] Scale factors err: {}".format(np.round(sf_err, 4)))

    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetOptStat(0)
    c1 = ROOT.TCanvas("c1", "", 800, 800)
    c1.cd()

    eff_MC.SetTitle(";p_{T} [GeV];Efficiency")
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

    c1.Update()
    efferr = eff_MC.GetPaintedGraph()
    efferr.GetYaxis().SetRangeUser(0., 1.2)

    directory = "../plots/1Deffunc"
    if not os.path.exists(directory):
        os.makedirs(directory)
    c1.Print("{}/1Deffunc_{}_{}.pdf".format(directory, varied, region))
    c1.Close()

    return sf, sf_err


if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True) # turn off GUI

    fData = ROOT.std.vector("string")((
        "../miniTree/UL2016preVFP/miniTree_Data_UL2016preVFP.root",
        "../miniTree/UL2016postVFP/miniTree_Data_UL2016postVFP.root",
        "../miniTree/UL2017/miniTree_Data_UL2017.root",
        "../miniTree/UL2018/miniTree_Data_UL2018.root"
    ))
    fZg = ROOT.std.vector("string")((
        "../miniTree/UL2016preVFP/miniTree_ZGToLLG_UL2016preVFP.root",
        "../miniTree/UL2016postVFP/miniTree_ZGToLLG_UL2016postVFP.root",
        "../miniTree/UL2017/miniTree_ZGToLLG_UL2017.root",
        "../miniTree/UL2018/miniTree_ZGToLLG_UL2018.root"
    ))
    fbkg = ROOT.std.vector("string")((
        "../miniTree/UL2016preVFP/miniTree_TTJets_UL2016preVFP.root",
        "../miniTree/UL2016postVFP/miniTree_TTJets_UL2016postVFP.root",
        "../miniTree/UL2017/miniTree_TTJets_UL2017.root",
        "../miniTree/UL2018/miniTree_TTJets_UL2018.root"
    ))
    fbkg.insert(fbkg.end(), fZg.begin(), fZg.end())
    fDY = ROOT.std.vector("string")((
        "../miniTree/UL2016preVFP/miniTree_DYJets_UL2016preVFP.root",
        "../miniTree/UL2016postVFP/miniTree_DYJets_UL2016postVFP.root",
        "../miniTree/UL2017/miniTree_DYJets_UL2017.root",
        "../miniTree/UL2018/miniTree_DYJets_UL2018.root"
    ))

    rdf_Da = ROOT.RDataFrame("miniTree", fData)
    rdf_Zg = ROOT.RDataFrame("miniTree", fZg)
    rdf_bkg = ROOT.RDataFrame("miniTree", fbkg)
    rdf_DY = ROOT.RDataFrame("miniTree", fDY)

    binArray = np.asarray([25, 35, 50, 150], "d")
    sysType = ["nominal", "rconv_up", "rconv_do", "pu_up", "pu_do", "bkg", "alt"]
    df_list = []
    for region in ["EB", "EE"]:
        SFDict = {}
        for varied in sysType:
            print(color.GREEN+"--------{} for {}--------".format(region, varied)+color.END)
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

    print(color.GREEN+"--------Draw SFs--------"+color.END)
    Draw2DSFs(df_list[0], df_list[1])

    for t in ["nominal_staunc", "rconv_sysunc", "pu_sysunc", "bkg_sysunc", "alt_sysunc"]:
        Draw2Dunc(df_list[0], df_list[1], t)
