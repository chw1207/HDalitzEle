import ROOT
import os, sys
from plugins.CMS_lumi import CMS_lumi
from plugins.colorPrint import *



def DrawFitting(var, dset, pdf, fitRes, region, outName, isMC):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    color = ROOT.kRed - 4 if isMC else ROOT.kAzure-3
    xframe = var.frame(-1, 1)
    dset.plotOn(xframe, ROOT.RooFit.Name("set"), ROOT.RooFit.Binning(100), ROOT.RooFit.MarkerStyle(ROOT.kFullCircle), ROOT.RooFit.MarkerSize(1.22), ROOT.RooFit.XErrorSize(10e-5))
    pdf.plotOn(xframe, ROOT.RooFit.Name("pdf"), ROOT.RooFit.LineColor(color), ROOT.RooFit.LineWidth(2))
    dset.plotOn(xframe, ROOT.RooFit.Name("set"), ROOT.RooFit.Binning(100), ROOT.RooFit.MarkerStyle(ROOT.kFullCircle), ROOT.RooFit.MarkerSize(1.22))
    # pdf.paramOn(xframe, ROOT.RooFit.Layout(0.7, 0.95, 0.95))

    xframe.SetTitle("")
    xframe.GetXaxis().SetTickSize(0.03)
    xframe.GetXaxis().SetTitleSize(0.045)
    xframe.GetXaxis().SetLabelSize(0.045)
    xframe.GetXaxis().SetLabelOffset(0.02)
    xframe.GetXaxis().SetTitleOffset(1.3)
    xframe.GetXaxis().SetTitle("s = E_{reco} / E_{kin} - 1")
    xframe.GetYaxis().SetTitle("Events / 0.01")
    xframe.GetYaxis().SetNdivisions(510)
    xframe.GetYaxis().SetTickSize(0.03)
    xframe.GetYaxis().SetTitleSize(0.045)
    xframe.GetYaxis().SetLabelSize(0.045)
    xframe.GetYaxis().SetTitleOffset(1.55)
    xframe.SetMaximum(xframe.GetMaximum() * 1.4)
    xframe.SetMinimum(0)

    c = ROOT.TCanvas("c", "", 700, 600)
    c.cd()
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.14)
    c.SetBottomMargin(0.13)
    xframe.Draw()

    ltx1 = ROOT.TLatex()
    ltx1.SetNDC()
    ltx1.SetTextFont(42)
    ltx1.SetTextSize(0.037)
    ltx1.DrawLatex(0.2, 0.86, region)

    leg1 = ROOT.TLegend(0.17, 0.72, 0.5, 0.84)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.037)
    leg1.SetFillColor(0)
    leg1.SetLineColor(0)

    name = "Simulation" if isMC else "Data"
    leg1.AddEntry(xframe.findObject("set"), name, "ep")
    # leg1.AddEntry(xframe.findObject("pdf"), "Breit-Wigner #otimes Crystal Ball", "l")
    leg1.AddEntry(xframe.findObject("pdf"), "Breit-Wigner #otimes Gaussian", "l")
    leg1.Draw("same")

    m = fitRes.floatParsFinal().find("mean")
    sig = fitRes.floatParsFinal().find("sigma_gau")
    gam = fitRes.floatParsFinal().find("width")
    ltx2 = ROOT.TLatex()
    ltx2.SetNDC()
    ltx2.SetTextFont(42)
    ltx2.SetTextSize(0.037)
    ltx2.DrawLatex(0.60, 0.86, "Mean = %.4f #pm %.4f" %(m.getVal(), m.getError()))
    ltx2.DrawLatex(0.60, 0.80, "#sigma = %.4f #pm %.4f" %(sig.getVal(), sig.getError()))
    ltx2.DrawLatex(0.60, 0.74, "#Gamma = %.4f #pm %.4f" %(gam.getVal(), gam.getError()))

    CMS_lumi(c, 5, 0, "137.2 fb^{-1}", 2017, True, "Work-in-progress", "", "")
    c.Update()
    c.RedrawAxis()

    outdir = "./plots/EnCalibFit"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    c.SaveAs("{}/{}.pdf".format(outdir, outName))
    c.SaveAs("{}/{}.png".format(outdir, outName))
    c.Close()
    print("")

def main(option, isMC):
    region = "isEBPho_lep1 == 1" if "EB" in option else "isEEPho_lep1 == 1"
    if isMC:
        rdf = rdf_Zg
    else:
        rdf = rdf_Da
    rdf = rdf.Filter(region)\
             .Define("Emes", "Z.M() * Z.M() - dimu.M() * dimu.M()")\
             .Define("Ekin", "91.188 * 91.188 - dimu.M() * dimu.M()")\
             .Define("scale", "(Emes / Ekin) - 1")

    # create RooDataset to fit
    s = ROOT.RooRealVar("scale", "scale", -1, 1)
    out = ["scale"] if isMC != True else ["scale", "puwei", "wei1", "wei2"]
    arr = rdf.AsNumpy(columns = out)

    s = ROOT.RooRealVar("scale", "scale", -1, 1)
    w = ROOT.RooRealVar("wei", "wei", 1.)
    aset = ROOT.RooArgSet(s, w)
    dset = ROOT.RooDataSet("dset", "dset", aset, "wei")
    for i in range(len(arr["scale"])):
        s.setVal(arr["scale"][i])
        if isMC:
            w.setVal(arr["wei2"][i])
        dset.add(aset, w.getVal())

    # if isMC:
    #     h = rdf.Histo1D(("h", " ", 80, -1, 1), "scale", "wei1").GetPtr()
    # else:
    #     h = rdf.Histo1D(("h", " ", 80, -1, 1), "scale").GetPtr()
    # dh = ROOT.RooDataHist("dh", "dh", s, ROOT.RooFit.Import(h))

    # fit function
    sigma_cb = ROOT.RooRealVar("sigma_cb", "sigma_cb", 0.1, 0, 0.2)
    mean = ROOT.RooRealVar("mean", "mean", 0.18, -0.2, 0.2)
    alpha = ROOT.RooRealVar("alpha", "alpha", 10, 1, 20)
    n = ROOT.RooRealVar("n", "n", 5, 1, 30)
    cball = ROOT.RooCBShape("cball", "cball", s, mean, sigma_cb, alpha, n)

    width = ROOT.RooRealVar("width", "width", 0.02, 0, 0.5)
    bw = ROOT.RooBreitWigner("bw", "bw", s, mean, width)

    sigma_gau = ROOT.RooRealVar("sigma_gau", "sigma_gau", 1, 0, 10)
    gauss = ROOT.RooGaussian("gauss", "gauss", s, mean, sigma_gau)

    # fitpdf = cball
    # fitpdf =  ROOT.RooFFTConvPdf("bxc", "bw (X) cball", s, bw, cball)
    fitpdf =  ROOT.RooFFTConvPdf("bxg", "bw (X) gauss", s, bw, gauss)

    # fitting
    s.setRange("fitRange", -1, 1)
    results = fitpdf.fitTo(dset, ROOT.RooFit.Save(ROOT.kTRUE), ROOT.RooFit.Range("fitRange"), ROOT.RooFit.Minimizer("Minuit", "minimize"), ROOT.RooFit.PrintLevel(-1))
    results.Print()
    fitpdf.Print()

    region = "ECAL Barrel" if "EB" in option else "ECAL Endcap"
    outName = "EnCalibFit_{}".format(option)
    DrawFitting(s, dset, fitpdf, results, region, outName, isMC)


if __name__ == "__main__":
    ROOT.EnableImplicitMT()
    # turn on the silent mode of RooFit
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    ROOT.RooMsgService.instance().setSilentMode(True)

    fData = [
        "./miniTree/2016_preVFP/miniTree_Data_2016_preVFP.root",
        "./miniTree/2016_postVFP/miniTree_Data_2016_postVFP.root",
        "./miniTree/2017/miniTree_Data_2017.root",
        "./miniTree/2018/miniTree_Data_2018.root"
    ]

    fZg = [
        "./miniTree/2016_preVFP/miniTree_ZGToLLG_2016_preVFP.root",
        "./miniTree/2016_postVFP/miniTree_ZGToLLG_2016_postVFP.root",
        "./miniTree/2017/miniTree_ZGToLLG_2017.root",
        "./miniTree/2018/miniTree_ZGToLLG_2018.root"
    ]

    sel_basic = "&&".join([
        "convVtxRadius_lep1 < 40",
        "isHggPho_lep1 == 1",
        "nGsfMatchToReco_lep1 > 1"
    ])

    rdf_Zg = ROOT.RDataFrame("outTree", fZg).Filter(sel_basic, "basic selection")
    rdf_Da = ROOT.RDataFrame("outTree", fData).Filter(sel_basic, "basic selection")

    for option in ["EB_Zg", "EB_Data", "EE_Zg", "EE_Data"]:
        isMC = True if "Data" not in option else False

        print(color.GREEN + "---------{}---------".format(option) + color.END)
        main(option, isMC)