# script to compare different energy reconstruction for merged electron in CMS
import sys, os
import ROOT
import numpy as np
import pickle as pkl
import xgboost as xgb
import pandas as pd
from glob import glob
from CMS_lumi import CMS_lumi
from colorPrint import *
from ROOT import RooFit
ROOT.gROOT.LoadMacro("../interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")

def fit_resp(hist, pt_do, pt_up, method):
    # define fiiting pdf
    x      = ROOT.RooRealVar("x", " ", 0, 2)
    mu      = ROOT.RooRealVar("mu",   " ",  hist.GetBinCenter(hist.GetMaximumBin()), 0.95, 1.05)
    sigma   = ROOT.RooRealVar("sigma",  " ", hist.GetStdDev(), 0.0, 0.5)
    alpha1  = ROOT.RooRealVar("alpha1", " ", 1,   1e-1,   3)
    alpha2  = ROOT.RooRealVar("alpha2", " ", 1,   1e-1,   3)
    n1      = ROOT.RooRealVar("n1",     " ", 80,    1,     200)
    n2      = ROOT.RooRealVar("n2",     " ", 80,    1,     200)
    pdf     = ROOT.RooCrystalBall("dcb", "dcb", x, mu, sigma, alpha1, n1, alpha2, n2)

    dh = ROOT.RooDataHist("dh", "dh", x, ROOT.RooFit.Import(hist))
    
    # perform the fitting
    x.setRange("fitRange", 0.65, 1.35)
    for i in range(10):
        results = pdf.fitTo(
            dh, 
            ROOT.RooFit.Minimizer("Minuit2", "minimize"),
            ROOT.RooFit.Range("fitRange"),
            ROOT.RooFit.Save(1), ROOT.RooFit.Strategy(2),
            ROOT.RooFit.SumW2Error(1), 
            ROOT.RooFit.PrintLevel(1)
        )
        if results.status() == 0:
            break
    results.Print("V")

    # visualization
    canv = ROOT.TCanvas("canv", "", 800, 650)
    canv.cd()
    canv.SetRightMargin(0.04)
    canv.SetTopMargin(0.075)
    canv.SetLeftMargin(0.14)
    canv.SetBottomMargin(0.14)

    # nbins = 60 if iBE == 0 else 400
    xframe = x.frame(0.65, 1.35)
    dh.plotOn(
        xframe,
        RooFit.Name("hist"),
        RooFit.MarkerStyle(ROOT.kFullCircle), RooFit.MarkerSize(1.2),
        RooFit.LineWidth(2)
    )
    pdf.plotOn(
        xframe,
        RooFit.Name("pdf"),
        RooFit.LineColor(ROOT.kRed - 4), RooFit.LineWidth(2)
    )
    dh.plotOn(
        xframe,
        RooFit.Name("hist"),
        RooFit.MarkerStyle(ROOT.kFullCircle), RooFit.MarkerSize(1.2),
        RooFit.LineWidth(2)
    )
    xframe.SetTitle("")
    xframe.GetXaxis().SetTitleSize(0.047)
    xframe.GetXaxis().SetLabelSize(0.046)
    xframe.GetXaxis().SetTitleOffset(1.24)
    xframe.GetXaxis().SetLabelOffset(0.012)
    xframe.GetXaxis().SetTitle("p^{reco}_{T} / p^{true}_{T}")
    # xframe.GetYaxis().SetTitle("Events / 0.005")
    xframe.GetYaxis().SetTitleSize(0.047)
    xframe.GetYaxis().SetLabelSize(0.046)
    xframe.GetYaxis().SetTitleOffset(1.4)
    xframe.SetMaximum(xframe.GetMaximum() * 1.5)
    xframe.SetMinimum(0)
    xframe.Draw()
    
    leg = ROOT.TLegend(0.17, 0.64, 0.45, 0.76)
    leg.SetTextFont(42)
    leg.SetTextSize(0.04)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.AddEntry(xframe.findObject("hist"), "Simulation", "LE1P")
    leg.AddEntry(xframe.findObject("pdf"), "DCB", "l")
    leg.Draw("same")

    var1 = results.floatParsFinal().find("mu")
    var2 = results.floatParsFinal().find("sigma")
    ltx = ROOT.TLatex()
    ltx.SetNDC()
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.045)
    ltx.DrawLatex(0.58, 0.86, " #mu = %.4f #pm %.4f" %(var1.getValV(), var1.getError()))
    ltx.DrawLatex(0.58, 0.81, " #sigma = %.4f #pm %.4f" %(var2.getValV(), var2.getError()))
    
    region = "ECAL Barrrel" if iBE == 0 else "ECAL Endcap"
    ltx.DrawLatex(0.18, 0.58, "%d < p^{true}_{T} < %d GeV" %(pt_do, pt_up))
    ltx.DrawLatex(0.18, 0.5, region)
    
    CMS_lumi(canv, 5, 10, "", 2017, True, "Simulation Preliminary", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    canv.Update()
    canv.RedrawAxis()

    reg_ext = "EB" if iBE == 0 else "EE"
    directory = "../plots/validation/resp_fit_HGGCheck/{}".format(reg_ext)
    if not os.path.exists(directory):
        os.makedirs(directory)
    canv.Print("{}/respFit_pt_{}_{}To{}GeV.pdf".format(directory, method, int(pt_do), int(pt_up)))
    canv.Close()
    
    return results


def fill_tgraph(h, ext):
    x_mu, y_mu, xerr_mu, yerr_mu = [], [], [], []
    x_width, y_width, xerr_width, yerr_width = [], [], [], []
    for i in range(h.GetNbinsX()):
        pt_do = h.GetXaxis().GetBinLowEdge(i+1) 
        pt_up = h.GetXaxis().GetBinUpEdge(i+1)
        print(color.GREEN + "[INFO] Perform the fitting in range [%d, %d]" %(int(pt_do), int(pt_up)) + color.END)
        histName = "hist1D_{}To{}".format(pt_do, pt_up) 
        hist = h.ProjectionY(histName, i+1, i+1).Clone()
        fitRes = fit_resp(hist, pt_do, pt_up, ext)
        
        x_mu.append(pt_do+(pt_up-pt_do)/2)
        y_mu.append(fitRes.floatParsFinal().find("mu").getValV())
        xerr_mu.append((pt_up-pt_do)/2)
        yerr_mu.append(fitRes.floatParsFinal().find("mu").getError())
        
        sigma_ave = fitRes.floatParsFinal().find("sigma").getValV()
        sigma_ave_err = fitRes.floatParsFinal().find("sigma").getError()
        x_width.append(pt_do+(pt_up-pt_do)/2)
        y_width.append(sigma_ave)
        xerr_width.append((pt_up-pt_do)/2)
        yerr_width.append(sigma_ave_err)
        
    mu_err = ROOT.TGraphErrors(h.GetNbinsX(), np.array(x_mu), np.array(y_mu), np.array(xerr_mu), np.array(yerr_mu)) 
    width_err = ROOT.TGraphErrors(h.GetNbinsX(), np.array(x_width), np.array(y_width), np.array(xerr_width), np.array(yerr_width)) 
    return mu_err, width_err


def main():
    # input features of dedicated regression
    features = [
        "rho",
        "nVtx",
        "eleSCEta_Lead",
        "eleSCPhi_Lead",
        "eleSCRawEn_Lead",
        "eleCalibPt_Lead",

        "eledEtaAtVtx_Lead",
        "eledPhiAtVtx_Lead",
        "elePtError_Lead",
        "eleHoverE_Lead",
        "eleEoverP_Lead",
        "eleEoverPout_Lead",
        "eleEoverPInv_Lead",

        "eleSCEtaWidth_Lead",
        "eleSCPhiWidth_Lead",
        "eleSigmaIEtaIEtaFull5x5_Lead",
        "eleSigmaIPhiIPhiFull5x5_Lead",
        "eleR9Full5x5_Lead",
        "eleBrem_Lead",

        "gsfPtSum_Lead",
        "gsfPtRatio_Lead",
        "diTrkPt_Lead",
        "gsfDeltaR_Lead"
    ]
    if iBE == 1:
        features.append("eleESEnToRawE_Lead")
        
    # preselection
    cut_base = "elePresel_Lead == 1 && category == 2"
    cut_region = "fabs(eleSCEta_Lead) < 1.479"
    if iBE == 1:
        cut_region = "fabs(eleSCEta_Lead) >= 1.479 && fabs(eleSCEta_Lead) < 2.5"
    
    # load data to pandas dataframe
    ROOT.EnableImplicitMT(20)
    infiles = glob("../../GenStudy/miniTree/*/miniTree_HDalitz_*_eeg_*.root")
    data = ROOT.RDataFrame("miniTree", infiles)\
               .Define("target",         "diGenEle.Pt()/elePt_Lead")\
               .Define("diTrkPt_Lead",   "diTrk.Pt()")\
               .Define("eleSCPt_Lead",   "eleSCEn_Lead/cosh(eleEta_Lead)")\
               .Define("eleECALPt_Lead", "eleEcalEn_Lead/cosh(eleEta_Lead)")\
               .Define("genPt",     "diGenEle.Pt()")\
               .Filter("{} && {}".format(cut_base, cut_region))\
               .AsNumpy(columns=features+["instwei", "target", "mcwei", "genwei", "puwei", "genPt", "elePt_Lead", "eleSCPt_Lead", "eleECALPt_Lead"])
    df = pd.DataFrame(data)
        
    # load the dedicated regression model and do the prediction
    region = "EB" if iBE == 0 else "EE"
    models_path = "/data4/chenghan/external/RegressionFinal/XGBRegression_NoRobustScaling_EGMRegTarget_{}/XGB_Regression.txt".format(region)
    reg = xgb.Booster()
    reg.load_model(models_path)
    preds = reg.predict(xgb.DMatrix(df.loc[:, features].values))
    df["eleHDALRegPt_Lead"] = preds * df["elePt_Lead"]
    
    # response of different reco methods 
    df["wei"]            = df["mcwei"] * df["genwei"] * df["puwei"]
    df["resp_eleHDALPt"] = df["eleHDALRegPt_Lead"]/df["genPt"]
    df["resp_elePt"]     = df["elePt_Lead"]/df["genPt"]
    df["resp_diTrkPt"]   = df["diTrkPt_Lead"]/df["genPt"]
    df["resp_eleSCPt"]   = df["eleSCPt_Lead"]/df["genPt"]
    df["resp_eleECALPt"] = df["eleECALPt_Lead"]/df["genPt"]
    weight = df["wei"].to_numpy()
    resp_eleHDALPt = df["resp_eleHDALPt"].to_numpy()
    resp_elePt     = df["resp_elePt"].to_numpy()
    # resp_diTrkPt   = df["resp_diTrkPt"].to_numpy()
    resp_eleSCPt   = df["resp_eleSCPt"].to_numpy()
    resp_eleECALPt = df["resp_eleECALPt"].to_numpy()
    genPt          = df["genPt"].to_numpy()
    
    bins = 280 if iBE == 0 else 350
    hist2D_eleHDALPt = ROOT.TH2D("hist2D_eleHDALPt", "", 24, 25, 125, bins, 0, 2)
    hist2D_elePt     = ROOT.TH2D("hist2D_elePt", "", 24, 25, 125, bins, 0, 2)   
    # hist2D_diTrkPt   = ROOT.TH2D("hist2D_diTrkPt", "", 13, 25, 125, 100, 0, 5)
    hist2D_eleSCPt   = ROOT.TH2D("hist2D_eleSCPt", "", 24, 25, 125, bins, 0, 2)   
    hist2D_eleECALPt = ROOT.TH2D("hist2D_eleECALPt", "", 24, 25, 125, bins, 0, 2)
    for i in range(len(resp_eleHDALPt)):
        hist2D_eleHDALPt.Fill(genPt[i], resp_eleHDALPt[i], weight[i])
        hist2D_elePt.Fill(genPt[i], resp_elePt[i], weight[i])
        # hist2D_diTrkPt.Fill(genPt[i], resp_diTrkPt[i], weight[i])
        hist2D_eleSCPt.Fill(genPt[i], resp_eleSCPt[i], weight[i])
        hist2D_eleECALPt.Fill(genPt[i], resp_eleECALPt[i], weight[i])
    
    mu_err_eleHDALPt, width_err_eleHDALPt = fill_tgraph(hist2D_eleHDALPt, "eleHDALPt")
    mu_err_elePt, width_err_elePt = fill_tgraph(hist2D_elePt, "elePt")
    # mu_err_diTrkPt, width_err_diTrkPt = fill_tgraph(hist2D_diTrkPt, "diTrkPt")
    mu_err_eleSCPt, width_err_eleSCPt = fill_tgraph(hist2D_eleSCPt, "eleSCPt")
    mu_err_eleECALPt, width_err_eleECALPt = fill_tgraph(hist2D_eleECALPt, "eleECALPt")
    
    c1 = ROOT.TCanvas("c1", "", 950, 800)
    c1.cd()
    # c1 = ROOT.TPad("c1", " ", 0, 0.3, 1, 1.0)
    c1.SetBottomMargin(0.12)
    c1.SetTopMargin(0.065)
    c1.SetRightMargin(0.04)
    c1.SetLeftMargin(0.15)
    
    mu_err_eleHDALPt.GetXaxis().SetTitle("")
    mu_err_eleHDALPt.GetYaxis().SetTitle("#LT p^{Reco}_{T} / p^{True}_{T} #GT")
    mu_err_eleHDALPt.GetYaxis().SetRangeUser(0.96, 1.02)
    mu_err_eleHDALPt.GetYaxis().SetNdivisions(506)
    mu_err_eleHDALPt.GetYaxis().SetTickSize(0.03)
    mu_err_eleHDALPt.GetYaxis().SetTitleSize(0.04)
    mu_err_eleHDALPt.GetYaxis().SetLabelSize(0.04)
    mu_err_eleHDALPt.GetYaxis().SetTitleOffset(1.7)

    mu_err_eleHDALPt.GetXaxis().SetRangeUser(18, 120)
    mu_err_eleHDALPt.GetXaxis().SetTickSize(0.03)
    mu_err_eleHDALPt.GetXaxis().SetTitleSize(0.04)
    mu_err_eleHDALPt.GetXaxis().SetLabelSize(0.04)
    mu_err_eleHDALPt.GetXaxis().SetTitleOffset(1.25)
    mu_err_eleHDALPt.GetXaxis().SetTitle("p^{True}_{T} (GeV)")
    mu_err_eleHDALPt.SetMarkerColor(ROOT.TColor.GetColor("#CC4054"))
    mu_err_eleHDALPt.SetMarkerSize(1.4)
    mu_err_eleHDALPt.SetMarkerStyle(20)
    mu_err_eleHDALPt.SetLineColor(ROOT.TColor.GetColor("#CC4054"))
    mu_err_eleHDALPt.SetLineWidth(2)
    mu_err_eleHDALPt.Draw("AP")

    mu_err_elePt.SetMarkerColor(ROOT.TColor.GetColor("#202020"))
    mu_err_elePt.SetMarkerSize(1.4)
    mu_err_elePt.SetMarkerStyle(20)
    mu_err_elePt.SetLineColor(ROOT.TColor.GetColor("#202020"))
    mu_err_elePt.SetLineWidth(2)
    mu_err_elePt.Draw("P same")
    
    mu_err_eleSCPt.SetMarkerColor(ROOT.TColor.GetColor("#019267"))
    mu_err_eleSCPt.SetMarkerSize(1.4)
    mu_err_eleSCPt.SetMarkerStyle(20)
    mu_err_eleSCPt.SetLineColor(ROOT.TColor.GetColor("#019267"))
    mu_err_eleSCPt.SetLineWidth(2)
    mu_err_eleSCPt.Draw("P same")
    
    mu_err_eleECALPt.SetMarkerColor(ROOT.TColor.GetColor("#2666CF"))
    mu_err_eleECALPt.SetMarkerSize(1.4)
    mu_err_eleECALPt.SetMarkerStyle(20)
    mu_err_eleECALPt.SetLineColor(ROOT.TColor.GetColor("#2666CF"))
    mu_err_eleECALPt.SetLineWidth(2)
    mu_err_eleECALPt.Draw("P same")
    
    # mu_err_diTrkPt.SetMarkerColor(ROOT.TColor.GetColor("#F0B86E"))
    # mu_err_diTrkPt.SetMarkerSize(1.4)
    # mu_err_diTrkPt.SetMarkerStyle(20)
    # mu_err_diTrkPt.SetLineColor(ROOT.TColor.GetColor("#F0B86E"))
    # mu_err_diTrkPt.SetLineWidth(2)
    # mu_err_diTrkPt.Draw("P same")
    
    # mu_err_eleHDALPt.Draw("P same")
    
    region = "ECAL Barrel" if iBE == 0 else "ECAL Endcap"
    leg = ROOT.TLegend(0.5, 0.15, 0.88, 0.45)
    leg.SetHeader(region)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.AddEntry(mu_err_eleHDALPt, "merged electron regression", "LE1P")
    leg.AddEntry(mu_err_elePt, "electron momentum", "LE1P")
    leg.AddEntry(mu_err_eleSCPt, "electron SC energy", "LE1P")
    leg.AddEntry(mu_err_eleECALPt, "electron ECAL energy", "LE1P")
    leg.Draw("same")
     
    l = ROOT.TLine(18, 1, 120, 1)
    l.SetLineStyle(2)
    l.SetLineColor(1)
    l.Draw()
    
    CMS_lumi(c1, 5, 10, "", 2017, True, "Simulation Preliminary", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    
    # c1.Print("test.pdf")
    # c1.Print("test.png")
    reg_ext = "EB" if iBE == 0 else "EE"
    directory = "../plots/validation/resp_mu_HGGCheck"
    if not os.path.exists(directory):
        os.makedirs(directory)

    c1.Print("{}/respMu_pt_{}.pdf".format(directory, reg_ext))
    c1.Close()
    

if __name__ == "__main__":
    iBE = int(sys.argv[1])
    ROOT.gROOT.SetBatch()
    ROOT.RooMsgService.instance().setSilentMode(True)
    
    main()