import sys, os
import ROOT
from ROOT import RooFit

import numpy as np
from CMS_lumi import CMS_lumi
from colorPrint import *
ROOT.gROOT.LoadMacro("../interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")


def fit_resp(hist, pt_do, pt_up, method):
    # define fiiting pdf
    # Cruijff is stable for asymmetrical distributions with long tails
    x      = ROOT.RooRealVar("x", " ", 0.8, 1.2)
    mu     = ROOT.RooRealVar("mu", " ", 1., 0.95, 1.05)
    sigmaL = ROOT.RooRealVar("sigmaL", " ", 0.002, 0.0, 0.5)
    sigmaR = ROOT.RooRealVar("sigmaR", " ", 0.002, 0.0, 0.5)
    alphaL = ROOT.RooRealVar("alphaL", " ", 0.1, 0, 2)
    alphaR = ROOT.RooRealVar("alphaR", " ", 0.1, 0, 2)
    pdf = ROOT.RooCruijff("cruijff", " ", x, mu, sigmaL, sigmaR, alphaL, alphaR)
    dh = ROOT.RooDataHist("dh", "dh", x, ROOT.RooFit.Import(hist))
    
    # perform the fitting
    x.setRange("fitRange", 0.8, 1.2)
    for i in range(10):
        results = pdf.fitTo(
            dh, 
            ROOT.RooFit.Minimizer("Minuit2", "minimize"),
            ROOT.RooFit.Range("fitRange"),
            ROOT.RooFit.Save(1), ROOT.RooFit.Strategy(1),
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

    nbins = 60 if iBE == 0 else 30
    xframe = x.frame(0.65, 1.2)
    dh.plotOn(
        xframe,
        RooFit.Name("hist"),
        RooFit.Binning(nbins),
        RooFit.MarkerStyle(ROOT.kFullCircle), RooFit.MarkerSize(1.2),
        RooFit.LineWidth(2)
    )
    pdf.plotOn(
        xframe,
        RooFit.Name("pdf"),
        RooFit.NormRange("fitRange"),
        RooFit.LineColor(ROOT.kRed - 4), RooFit.LineWidth(2)
    )
    dh.plotOn(
        xframe,
        RooFit.Name("hist"),
        RooFit.Binning(nbins),
        RooFit.MarkerStyle(ROOT.kFullCircle), RooFit.MarkerSize(1.2),
        RooFit.LineWidth(2)
    )
    xframe.SetTitle("")
    xframe.GetXaxis().SetTitleSize(0.047)
    xframe.GetXaxis().SetLabelSize(0.046)
    xframe.GetXaxis().SetTitleOffset(1.24)
    xframe.GetXaxis().SetLabelOffset(0.012)
    xframe.GetXaxis().SetTitle("p^{reco}_{T} / p^{true}_{T}")
    xframe.GetYaxis().SetTitle("Events / 0.005")
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
    leg.AddEntry(xframe.findObject("pdf"), "Cruijff", "l")
    leg.Draw("same")

    var1 = results.floatParsFinal().find("mu")
    var2 = results.floatParsFinal().find("sigmaL")
    var3 = results.floatParsFinal().find("sigmaR")
    ltx = ROOT.TLatex()
    ltx.SetNDC()
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.045)
    ltx.DrawLatex(0.58, 0.86, " #mu = %.4f #pm %.4f" %(var1.getValV(), var1.getError()))
    ltx.DrawLatex(0.58, 0.81, " #sigma_{L} = %.4f #pm %.4f" %(var2.getValV(), var2.getError()))
    ltx.DrawLatex(0.58, 0.76, " #sigma_{R} = %.4f #pm %.4f" %(var3.getValV(), var3.getError()))
    ltx.DrawLatex(0.58, 0.68, " #sigma_{ave} #equiv #frac{#sigma_{L} + #sigma_{R}}{2} = %.4f" %((var2.getValV()+var3.getValV())/2))
    
    region = "ECAL Barrrel" if iBE == 0 else "ECAL Endcap"
    ltx.DrawLatex(0.18, 0.58, "%d < p^{true}_{T} < %d GeV" %(pt_do, pt_up))
    ltx.DrawLatex(0.18, 0.5, region)
    
    CMS_lumi(canv, 5, 10, "", 2017, True, "Simulation Preliminary", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    canv.Update()
    canv.RedrawAxis()

    reg_ext = "EB" if iBE == 0 else "EE"
    directory = "../plots/validation/resp_fit/{}".format(reg_ext)
    if not os.path.exists(directory):
        os.makedirs(directory)
    canv.Print("{}/respFit_pt_{}_{}To{}GeV.pdf".format(directory, method, int(pt_do), int(pt_up)))
    canv.Close()
    
    return results


def plot_mu(mu_err_XGB, mu_err_EGM):
    c1 = ROOT.TCanvas("c1", "", 800, 800)
    c1.cd()

    pad1 = ROOT.TPad("pad1", " ", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.03)
    pad1.SetTopMargin(0.09)
    pad1.SetRightMargin(0.05)
    pad1.SetLeftMargin(0.14)
    pad1.Draw()
    pad1.cd()
    
    mu_err_XGB.GetXaxis().SetTitle("")
    mu_err_XGB.GetYaxis().SetTitle("#LT p^{reco}_{T} / p^{true}_{T} #GT")
    mu_err_XGB.GetYaxis().SetRangeUser(0.98, 1.02)
    mu_err_XGB.GetYaxis().SetNdivisions(506)
    mu_err_XGB.GetYaxis().SetTickSize(0.03)
    mu_err_XGB.GetYaxis().SetTitleSize(0.06)
    mu_err_XGB.GetYaxis().SetLabelSize(0.06)
    mu_err_XGB.GetYaxis().SetTitleOffset(1.15)

    mu_err_XGB.GetXaxis().SetRangeUser(25, 194)
    mu_err_XGB.GetXaxis().SetTickSize(0.03)
    mu_err_XGB.GetXaxis().SetTitleSize(0.06)
    mu_err_XGB.GetXaxis().SetLabelSize(0.05)
    mu_err_XGB.GetXaxis().SetLabelOffset(0.1)
    mu_err_XGB.GetXaxis().SetTitleOffset(1)
    mu_err_XGB.SetMarkerColor(ROOT.TColor.GetColor("#CC4054"))
    mu_err_XGB.SetMarkerSize(1.4)
    mu_err_XGB.SetMarkerStyle(20)
    mu_err_XGB.SetLineColor(ROOT.TColor.GetColor("#CC4054"))
    mu_err_XGB.SetLineWidth(2)
    mu_err_XGB.Draw("AP")

    mu_err_EGM.SetMarkerColor(ROOT.TColor.GetColor("#202020"))
    mu_err_EGM.SetMarkerSize(1.4)
    mu_err_EGM.SetMarkerStyle(20)
    mu_err_EGM.SetLineColor(ROOT.TColor.GetColor("#202020"))
    mu_err_EGM.SetLineWidth(2)
    
    l = ROOT.TLine(25, 1, 194, 1)
    l.SetLineStyle(2)
    l.SetLineColor(1)
    l.Draw()
    
    mu_err_EGM.Draw("P same")
    mu_err_XGB.Draw("P same")
    
    region = "ECAL Barrrel" if iBE == 0 else "ECAL Endcap"
    leg = ROOT.TLegend(0.62, 0.66, 0.88, 0.86)
    leg.SetHeader(region)
    leg.SetTextFont(42)
    leg.SetTextSize(0.05)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.AddEntry(mu_err_XGB, "XGB", "LE1P")
    leg.AddEntry(mu_err_EGM, "EGM", "LE1P")
    leg.Draw("same")
    
    CMS_lumi(pad1, 5, 10, "", 2017, True, "Simulation Preliminary", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    c1.cd()

    pad2 = ROOT.TPad("pad2", "", 0, 0, 1, 0.3)
    pad2.SetGridy()
    pad2.SetRightMargin(0.05)
    pad2.SetLeftMargin(0.14)
    pad2.SetTopMargin(0.06)
    pad2.SetBottomMargin(0.43)
    pad2.Draw()
    pad2.cd()
    
    nPoints = mu_err_XGB.GetN()
    ratio_err = ROOT.TGraphErrors(nPoints)
    for i in range(nPoints):
        xp = mu_err_XGB.GetPointX(i)
        yp = mu_err_XGB.GetPointY(i)/mu_err_EGM.GetPointY(i)
        xp_err = mu_err_XGB.GetErrorX(i)
        yp_err = yp * np.sqrt(pow(mu_err_XGB.GetErrorY(i)/mu_err_XGB.GetPointY(i), 2) + pow(mu_err_EGM.GetErrorY(i)/mu_err_EGM.GetPointY(i), 2))
        
        ratio_err.SetPoint(i, xp, yp)
        ratio_err.SetPointError(i, xp_err, yp_err)

    ratio_err.SetName("")
    ratio_err.SetTitle("")
    ratio_err.GetXaxis().SetTitle("p^{true}_{T} [GeV]")
    ratio_err.GetYaxis().SetTitle("Ratio")
    ratio_err.GetYaxis().SetRangeUser(0.97 , 1.03)

    ratio_err.SetMarkerColor(ROOT.TColor.GetColor("#202020"))
    ratio_err.SetMarkerSize(1.4)
    ratio_err.SetMarkerStyle(20)
    ratio_err.SetLineColor(ROOT.TColor.GetColor("#202020"))
    ratio_err.SetLineWidth(2)

    ratio_err.GetXaxis().SetRangeUser(25, 194)
    ratio_err.GetXaxis().SetTickSize(0.03 * (7/3.))
    ratio_err.GetXaxis().SetTitleSize(0.16)
    ratio_err.GetXaxis().SetTitleOffset(1.2)
    ratio_err.GetXaxis().SetLabelSize(0.06  * (7/3.))
    ratio_err.GetXaxis().SetLabelOffset(0.05)
    ratio_err.GetYaxis().SetTitleSize(0.13)
    ratio_err.GetYaxis().SetTitleOffset(0.22 * (7/3.))
    ratio_err.GetYaxis().SetLabelSize(0.06  * (7/3.))
    ratio_err.GetYaxis().SetNdivisions(502)
    ratio_err.Draw("AP")

    reg_ext = "EB" if iBE == 0 else "EE"
    directory = "../plots/validation/resp_mu"
    if not os.path.exists(directory):
        os.makedirs(directory)
    c1.Print("{}/respMu_pt_{}.pdf".format(directory, reg_ext))
    c1.Close()


def plot_width(width_err_XGB, width_err_EGM):
    c1 = ROOT.TCanvas("c1", "", 800, 800)
    c1.cd()

    pad1 = ROOT.TPad("pad1", " ", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.03)
    pad1.SetTopMargin(0.09)
    pad1.SetRightMargin(0.05)
    pad1.SetLeftMargin(0.14)
    pad1.Draw()
    pad1.cd()
    
    width_err_XGB.GetXaxis().SetTitle("")
    width_err_XGB.GetYaxis().SetTitle("#sigma^{ave}_{p^{reco}_{T} / p^{true}_{T}}")
    width_err_XGB.GetYaxis().SetRangeUser(0, 0.06)
    width_err_XGB.GetYaxis().SetNdivisions(506)
    width_err_XGB.GetYaxis().SetTickSize(0.03)
    width_err_XGB.GetYaxis().SetTitleSize(0.06)
    width_err_XGB.GetYaxis().SetLabelSize(0.06)
    width_err_XGB.GetYaxis().SetTitleOffset(1.15)

    width_err_XGB.GetXaxis().SetRangeUser(25, 194)
    width_err_XGB.GetXaxis().SetTickSize(0.03)
    width_err_XGB.GetXaxis().SetTitleSize(0.06)
    width_err_XGB.GetXaxis().SetLabelSize(0.05)
    width_err_XGB.GetXaxis().SetLabelOffset(0.1)
    width_err_XGB.GetXaxis().SetTitleOffset(1)
    width_err_XGB.SetMarkerColor(ROOT.TColor.GetColor("#CC4054"))
    width_err_XGB.SetMarkerSize(1.4)
    width_err_XGB.SetMarkerStyle(20)
    width_err_XGB.SetLineColor(ROOT.TColor.GetColor("#CC4054"))
    width_err_XGB.SetLineWidth(2)
    width_err_XGB.Draw("AP")

    width_err_EGM.SetMarkerColor(ROOT.TColor.GetColor("#202020"))
    width_err_EGM.SetMarkerSize(1.4)
    width_err_EGM.SetMarkerStyle(20)
    width_err_EGM.SetLineColor(ROOT.TColor.GetColor("#202020"))
    width_err_EGM.SetLineWidth(2)
    
    # l = ROOT.TLine(25, 1, 194, 1)
    # l.SetLineStyle(2)
    # l.SetLineColor(1)
    # l.Draw()
    
    width_err_EGM.Draw("P same")
    width_err_XGB.Draw("P same")
    
    region = "ECAL Barrrel" if iBE == 0 else "ECAL Endcap"
    leg = ROOT.TLegend(0.62, 0.66, 0.88, 0.86)
    leg.SetHeader(region)
    leg.SetTextFont(42)
    leg.SetTextSize(0.05)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.AddEntry(mu_err_XGB, "dedicated", "LE1P")
    leg.AddEntry(mu_err_EGM, "standard", "LE1P")
    leg.Draw("same")
    
    CMS_lumi(pad1, 5, 10, "", 2017, True, "Simulation Preliminary", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    c1.cd()

    pad2 = ROOT.TPad("pad2", "", 0, 0, 1, 0.3)
    pad2.SetGridy()
    pad2.SetRightMargin(0.05)
    pad2.SetLeftMargin(0.14)
    pad2.SetTopMargin(0.06)
    pad2.SetBottomMargin(0.43)
    pad2.Draw()
    pad2.cd()
    
    nPoints = width_err_XGB.GetN()
    ratio_err = ROOT.TGraphErrors(nPoints)
    for i in range(nPoints):
        xp = width_err_XGB.GetPointX(i)
        yp = width_err_XGB.GetPointY(i)/width_err_EGM.GetPointY(i)
        xp_err = width_err_XGB.GetErrorX(i)
        yp_err = yp * np.sqrt(pow(width_err_XGB.GetErrorY(i)/width_err_XGB.GetPointY(i), 2) + pow(width_err_EGM.GetErrorY(i)/width_err_EGM.GetPointY(i), 2))
        
        ratio_err.SetPoint(i, xp, yp)
        ratio_err.SetPointError(i, xp_err, yp_err)

    ratio_err.SetName("")
    ratio_err.SetTitle("")
    ratio_err.GetXaxis().SetTitle("p^{true}_{T} [GeV]")
    ratio_err.GetYaxis().SetTitle("Ratio")
    ratio_err.GetYaxis().SetRangeUser(0.6 , 1.4)

    ratio_err.SetMarkerColor(ROOT.TColor.GetColor("#202020"))
    ratio_err.SetMarkerSize(1.4)
    ratio_err.SetMarkerStyle(20)
    ratio_err.SetLineColor(ROOT.TColor.GetColor("#202020"))
    ratio_err.SetLineWidth(2)

    ratio_err.GetXaxis().SetRangeUser(25, 194)
    ratio_err.GetXaxis().SetTickSize(0.03 * (7/3.))
    ratio_err.GetXaxis().SetTitleSize(0.16)
    ratio_err.GetXaxis().SetTitleOffset(1.2)
    ratio_err.GetXaxis().SetLabelSize(0.06  * (7/3.))
    ratio_err.GetXaxis().SetLabelOffset(0.05)
    ratio_err.GetYaxis().SetTitleSize(0.13)
    ratio_err.GetYaxis().SetTitleOffset(0.22 * (7/3.))
    ratio_err.GetYaxis().SetLabelSize(0.06  * (7/3.))
    ratio_err.GetYaxis().SetNdivisions(502)
    ratio_err.Draw("AP")

    reg_ext = "EB" if iBE == 0 else "EE"
    directory = "../plots/validation/resp_width"
    if not os.path.exists(directory):
        os.makedirs(directory)
    c1.Print("{}/respWidth_pt_{}.pdf".format(directory, reg_ext))
    c1.Close()
    


if __name__ == "__main__":
    iBE = int(sys.argv[1])
    ROOT.gROOT.SetBatch()
    ROOT.RooMsgService.instance().setSilentMode(True)
    
    ROOT.EnableImplicitMT(20)
    df = ROOT.RDataFrame("regTree", "../reg_signal_XGB_{}_Pt25To200GeV.root".format(iBE))
    
    bins = 300 if iBE == 0 else 90
    hist2D_XGB = df.Histo2D(("hist2D_XGB", "", 13, 25, 194, bins, 0.8, 1.2), "genPt", "resp_XGB", "weight").GetPtr()
    hist2D_EGM = df.Histo2D(("hist2D_EGM", "", 13, 25, 194, bins, 0.8, 1.2), "genPt", "resp_EGMCalib", "weight").GetPtr()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    x_mu, y_mu, xerr_mu, yerr_mu = [], [], [], []
    x_width, y_width, xerr_width, yerr_width = [], [], [], []
    for i in range(hist2D_XGB.GetNbinsX()):
        pt_do = hist2D_XGB.GetXaxis().GetBinLowEdge(i+1) 
        pt_up = hist2D_XGB.GetXaxis().GetBinUpEdge(i+1)
        print(color.GREEN + "[INFO] Perform the fitting in range [%d, %d]" %(int(pt_do), int(pt_up)) + color.END)
        histName_XGB = "hist1D_XGB_{}To{}".format(pt_do, pt_up)
        hist_XGB = hist2D_XGB.ProjectionY(histName_XGB, i+1, i+1).Clone()
        fitRes = fit_resp(hist_XGB, pt_do, pt_up, "XGB")
        
        x_mu.append(pt_do+(pt_up-pt_do)/2)
        y_mu.append(fitRes.floatParsFinal().find("mu").getValV())
        xerr_mu.append((pt_up-pt_do)/2)
        yerr_mu.append(fitRes.floatParsFinal().find("mu").getError())
        
        sigma_ave = 0.5 * (fitRes.floatParsFinal().find("sigmaL").getValV() + fitRes.floatParsFinal().find("sigmaR").getValV())
        sigma_ave_err = np.sqrt(pow(fitRes.floatParsFinal().find("sigmaL").getError(), 2) + pow(fitRes.floatParsFinal().find("sigmaR").getError(), 2))
        x_width.append(pt_do+(pt_up-pt_do)/2)
        y_width.append(sigma_ave)
        xerr_width.append((pt_up-pt_do)/2)
        yerr_width.append(sigma_ave_err)
        
        print("")
    mu_err_XGB = ROOT.TGraphErrors(hist2D_XGB.GetNbinsX(), np.array(x_mu), np.array(y_mu), np.array(xerr_mu), np.array(yerr_mu)) 
    width_err_XGB = ROOT.TGraphErrors(hist2D_XGB.GetNbinsX(), np.array(x_width), np.array(y_width), np.array(xerr_width), np.array(yerr_width)) 
    width_err_XGB.Print()
    
    x_mu_egm, y_mu_egm, xerr_mu_egm, yerr_mu_egm = [], [], [], []
    x_width_egm, y_width_egm, xerr_width_egm, yerr_width_egm = [], [], [], []
    for i in range(hist2D_EGM.GetNbinsX()):
        pt_do = hist2D_EGM.GetXaxis().GetBinLowEdge(i+1) 
        pt_up = hist2D_EGM.GetXaxis().GetBinUpEdge(i+1)
        print(color.GREEN + "[INFO] Perform the fitting in range [%d, %d]" %(int(pt_do), int(pt_up)) + color.END)
        histName_EGM = "hist1D_EGM_{}To{}".format(pt_do, pt_up)
        hist_EGM = hist2D_EGM.ProjectionY(histName_EGM, i+1, i+1).Clone()
        fitRes = fit_resp(hist_EGM, pt_do, pt_up, "EGM")
        
        x_mu_egm.append(pt_do+(pt_up-pt_do)/2)
        y_mu_egm.append(fitRes.floatParsFinal().find("mu").getValV())
        xerr_mu_egm.append((pt_up-pt_do)/2)
        yerr_mu_egm.append(fitRes.floatParsFinal().find("mu").getError())
        
        sigma_ave = 0.5 * (fitRes.floatParsFinal().find("sigmaL").getValV() + fitRes.floatParsFinal().find("sigmaR").getValV())
        sigma_ave_err = np.sqrt(pow(fitRes.floatParsFinal().find("sigmaL").getError(), 2) + pow(fitRes.floatParsFinal().find("sigmaR").getError(), 2))
        x_width_egm.append(pt_do+(pt_up-pt_do)/2)
        y_width_egm.append(sigma_ave)
        xerr_width_egm.append((pt_up-pt_do)/2)
        yerr_width_egm.append(sigma_ave_err)
        
        print("")
    mu_err_EGM = ROOT.TGraphErrors(hist2D_EGM.GetNbinsX(), np.array(x_mu_egm), np.array(y_mu_egm), np.array(xerr_mu_egm), np.array(yerr_mu_egm)) 
    width_err_EGM = ROOT.TGraphErrors(hist2D_EGM.GetNbinsX(), np.array(x_width_egm), np.array(y_width_egm), np.array(xerr_width_egm), np.array(yerr_width_egm)) 
    width_err_EGM.Print()
    
    
    plot_mu(mu_err_XGB, mu_err_EGM)
    plot_width(width_err_XGB, width_err_EGM)    