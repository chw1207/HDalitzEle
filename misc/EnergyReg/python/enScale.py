import os, sys
import json
import ROOT
import numpy as np
from colorPrint import *
from CMS_lumi import CMS_lumi
from argparse import ArgumentParser

ROOT.gROOT.LoadMacro("../interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")

def get_parser():
    parser = ArgumentParser(
        description="macro to measure the energy scale differences between data and mc"
    )
    parser.add_argument(
        "-r", "--region",
        help="eta region [ EB | EE ], (default: EB)",
        default="EB",
        type=str
    )
    parser.add_argument(
        "-b", "--bin_array",
        help="Et bins of scale differences, (default: 25,35,50,150)",
        default="25,35,50,150",
        type=str
    )
    parser.add_argument(
        "-m", "--MaxTries",
        help="max tries when the fitting failed",
        default=10,
        type=int
    )
    parser.add_argument(
        "-op", "--out_dir_plot",
        help="directory to save the plots",
        default="../plots/energy_scale_diff",
        type=str
    )
    parser.add_argument(
        "-os", "--out_dir_scale",
        help="directory to save the scales root file (will be omitted if plotonly)",
        default="../data",
        type=str
    )
    parser.add_argument(
        "--unBinFit",
        action="store_true",
        help="un binned fit"
    )
    parser.add_argument(
        "--plotonly",
        action="store_true",
        help="do not save the scales to root file, plot only"
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print the detailed information of fitting"
    )
    return parser


def DrawFitting(var, dset, pdf, fitRes, isMC, ptcat, binwidth):
    cc = ROOT.kRed - 4 if isMC else ROOT.kAzure-3
    xframe = var.frame(-1, 1)

    dset.plotOn(
        xframe,
        ROOT.RooFit.Name("hist"),
        ROOT.RooFit.MarkerStyle(ROOT.kFullCircle), ROOT.RooFit.MarkerSize(1.2),
        ROOT.RooFit.LineWidth(2)
    )
    pdf.plotOn(
        xframe,
        ROOT.RooFit.Name("pdf"),
        ROOT.RooFit.LineColor(cc), ROOT.RooFit.LineWidth(2)
    )
    dset.plotOn(
        xframe,
        ROOT.RooFit.Name("hist"),
        ROOT.RooFit.MarkerStyle(ROOT.kFullCircle), ROOT.RooFit.MarkerSize(1.2),
        ROOT.RooFit.LineWidth(2)
    )
    if not args.unBinFit:
        n_param = fitRes.floatParsFinal().getSize()
        chi2 = xframe.chiSquare("pdf", "hist", n_param)
        print("[INFO] redued chi-square: {}".format(chi2))

    xframe.SetTitle("")
    xframe.GetXaxis().SetTitleSize(0.047)
    xframe.GetXaxis().SetLabelSize(0.046)
    xframe.GetXaxis().SetTitleOffset(1.24)
    xframe.GetXaxis().SetTitle("s = E_{reco} / E_{kin} - 1")
    xframe.GetYaxis().SetTitle("Events / {:.2f}".format((1. - (-1.))/binwidth))
    xframe.GetYaxis().SetTitleSize(0.047)
    xframe.GetYaxis().SetLabelSize(0.046)
    xframe.GetYaxis().SetTitleOffset(1.35)
    xframe.SetMaximum(xframe.GetMaximum() * 1.6)
    xframe.SetMinimum(0)

    c = ROOT.TCanvas("c", "", 800, 800)
    c.cd()
    # c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.14)
    # c.SetBottomMargin(0.13)
    xframe.Draw()

    ltx1 = ROOT.TLatex()
    ltx1.SetNDC()
    ltx1.SetTextFont(42)
    ltx1.SetTextSize(0.04)
    ltx1.DrawLatex(0.2, 0.86, "Z #rightarrow #mu#mu#gamma,  {}".format("r_{conv} < 40 cm"))
    ltx1.DrawLatex(0.2, 0.79, args.region)

    leg1 = ROOT.TLegend(0.17, 0.64, 0.5, 0.76)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.04)
    leg1.SetFillColor(0)
    leg1.SetLineColor(0)

    name = "Simulation" if isMC else "Data"
    leg1.AddEntry(xframe.findObject("hist"), name, "LE1P")
    leg1.AddEntry(xframe.findObject("pdf"), "DCB #otimes Breit-Wigner", "l")
    leg1.Draw("same")

    MeanPar = fitRes.floatParsFinal().find("mean")
    ltx = ROOT.TLatex()
    ltx.SetNDC()
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.04)
    ltx.DrawLatex(0.65, 0.86, " Mean = %.4f" %(MeanPar.getValV()))
    ltx.DrawLatex(0.65, 0.81, " #sigma = %.4f" %(fitRes.floatParsFinal().find("sigma").getValV()))
    ltx.DrawLatex(0.65, 0.76, " #alpha_{1} = %.4f" %(fitRes.floatParsFinal().find("alpha1").getValV()))
    ltx.DrawLatex(0.65, 0.71, " n_{1} = %.4f" %(fitRes.floatParsFinal().find("n1").getValV()))
    ltx.DrawLatex(0.65, 0.66, " #alpha_{2} = %.4f" %(fitRes.floatParsFinal().find("alpha2").getValV()))
    ltx.DrawLatex(0.65, 0.61, " n_{2} =  %.4f" %(fitRes.floatParsFinal().find("n2").getValV()))
    ltx.DrawLatex(0.65, 0.56, " #Gamma = %.2e" %(fitRes.floatParsFinal().find("width").getValV()))

    CMS_lumi(c, 5, 0, "138 fb^{-1}", 2017, True, "Preliminary", "", "")
    c.Update()
    c.RedrawAxis()

    ext = "MC" if isMC else "Data"
    directory = "{}/energy_scale_fit".format(args.out_dir_plot)
    if not os.path.exists(directory):
        os.makedirs(directory)
    c.Print("{}/EnScaleFit_ptcat{}_{}_{}.pdf".format(directory, ptcat, args.region, ext))
    c.Close()


# The scale corrections is defined as the relative shift between data and mc
# https://arxiv.org/pdf/1306.2016.pdf
# scale correction = 1 + deltaP,
# deltaP = (s_mc - s_data) -> scale difference
def plotScale():
    diff = np.asarray(sDict["{}_MC".format(args.region)], "d") - np.asarray(sDict["{}_Data".format(args.region)], "d")
    diff_err = np.sqrt(np.power(eDict["{}_Data".format(args.region)], 2) + np.power(eDict["{}_MC".format(args.region)], 2))

    x, y, ex, ey = [], [], [], []
    for i in range(len(bin_array)-1):
        binCenter = bin_array[i] + (bin_array[i+1] - bin_array[i])/2
        binContent = diff[i]
        binXErr = abs((bin_array[i+1] - bin_array[i])/2)
        binYErr = diff_err[i]

        x.append(binCenter)
        y.append(binContent)
        ex.append(binXErr)
        ey.append(binYErr)

    if (not args.plotonly):
        if not os.path.exists(args.out_dir_scale):
            os.makedirs(args.out_dir_scale)

        outName = "{}/EnScaleCorr_combined_{}.root".format(args.out_dir_scale, args.region)
        print("[INFO] Save scales correction in: {}".format(outName))

        fout = ROOT.TFile(outName, "RECREATE")
        fout.cd()
        hcorr = ROOT.TH1D("hscorr", "", len(bin_array)-1, np.asarray(bin_array, "d"))
        for i in range(len(y)):
            hcorr.SetBinContent(i+1, 1+y[i])
            hcorr.SetBinError(i+1, ey[i])
        hcorr.Write()
        fout.Close()

    print("[INFO] Energy scale correction:", flush=True)
    for i in range(len(diff)):
        print("scale diff = {:8.4f} +/- {:8.4f}".format(diff[i], diff_err[i]), flush=True)
    print("", flush=True)

    SFerr = ROOT.TGraphErrors(len(bin_array), np.array(x), np.array(y), np.array(ex), np.array(ey))
    c = ROOT.TCanvas("c", "", 800, 800)
    c.cd( )
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.17)
    c.SetBottomMargin(0.13)

    SFerr.SetTitle("")
    SFerr.GetXaxis().SetTickSize(0.03)
    SFerr.GetXaxis().SetTitleSize(0.045)
    SFerr.GetXaxis().SetLabelSize(0.045)
    SFerr.GetXaxis().SetLabelOffset(0.02)
    SFerr.GetXaxis().SetTitleOffset(1.3)
    SFerr.GetXaxis().SetRangeUser(bin_array[0]-0.0001, bin_array[-1]+0.0001)
    SFerr.GetXaxis().SetTitle("p^{corr}_{T} [GeV]")

    SFerr.GetYaxis().SetTitle("enrgy scale difference")
    SFerr.GetYaxis().SetRangeUser(-0.025, 0.015)
    SFerr.GetYaxis().SetTickSize(0.03)
    SFerr.GetYaxis().SetTitleSize(0.045)
    SFerr.GetYaxis().SetLabelSize(0.045)
    SFerr.GetYaxis().SetTitleOffset(1.92)
    SFerr.SetFillColor(ROOT.TColor.GetColor("#202020"))
    SFerr.SetLineColor(ROOT.TColor.GetColor("#202020"))
    SFerr.SetMarkerSize(1.5)
    SFerr.SetMarkerStyle(20)
    SFerr.SetLineWidth(2)
    SFerr.Draw("AP")

    extral = "ECAL Barrel" if "EB" in args.region else "ECAL Endcap"
    ltx = ROOT.TLatex()
    ltx.SetNDC()
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.05)
    ltx.DrawLatex(0.45, 0.28, "Z #rightarrow #mu#mu#gamma, {}".format("r_{conv} < 40 cm"))
    ltx.DrawLatex(0.45, 0.20, "{}".format(extral))

    l = ROOT.TLine(bin_array[0], 0, bin_array[-1], 0)
    l.SetLineStyle(7)
    l.SetLineWidth(2)
    l.SetLineColor(ROOT.kRed - 4)
    l.Draw()
    SFerr.Draw("P same")

    CMS_lumi(c, 5, 0, "138 fb^{-1}", 2017, True, "Preliminary", "", "")
    c.RedrawAxis()

    directory = "{}/energy_scale".format(args.out_dir_plot)
    if not os.path.exists(directory):
        os.makedirs(directory)
    c.Print("{}/EnScale_{}.pdf".format(directory, args.region))
    c.Close()


def main(isMC):
    rdf = rdf_Zg if isMC else rdf_Da
    rdf = rdf.Define("Emes",     "Z_reg.M() * Z_reg.M() - dimu.M() * dimu.M()")\
             .Define("Ekin",     "91.188 * 91.188 - dimu.M() * dimu.M()")\
             .Define("scale",    "(Emes / Ekin) - 1")

    # start to fit the scale distributions
    # ------------------------------------
    scaleList, errList = [], []
    min_float = sys.float_info.min # initial guess approches to 0
    # bins_mc = [90, 80, 60] if args.region == "EB" else [100, 85, 60]
    # bins_data = [60, 50, 40] if args.region == "EB" else [60, 50, 40]
    bins_mc = [85, 82, 60] if args.region == "EB" else [90, 80, 65]
    bins_data = [33, 28, 25] if args.region == "EB" else [37, 30, 22]
    for i in range(len(bin_array) - 1):
        rdf_pt = rdf.Filter("eleHDALRegPt_Lead >= {} && eleHDALRegPt_Lead < {}".format(bin_array[i], bin_array[i+1]))

        # count = rdf_pt.Count().GetValue()
        # bins = int(2 * count**(2/5.))
        if args.unBinFit:
            bins = 50
        else:
            bins = bins_mc[i] if isMC else bins_data[i]
            print("[INFO] {} bin is used to constrcut scale offset histogram".format(bins), flush=True)

        # create RooDataHist or DataSet to fit
        # ------------------------------------
        s       = ROOT.RooRealVar("scale", "scale", -1, 1)
        s.setBins(bins)
        if args.unBinFit:
            arr = rdf_pt.AsNumpy(columns = ["scale", "wei"]) if isMC else rdf_pt.AsNumpy(columns = ["scale"])
            wei       = ROOT.RooRealVar("weight", "weight", 1.)
            aset      = ROOT.RooArgSet(s, wei)
            dset      = ROOT.RooDataSet("dset", "dset", aset, "weight")
            for j in range(len(arr["scale"])):
                weight = arr["wei"][j] if isMC else 1.
                wei.setVal(weight)
                s.setVal(arr["scale"][j])
                dset.add(aset, weight)
            fitData = dset
            
        else:
            h = rdf_pt.Histo1D(("h", " ", bins, -1, 1), "scale", "wei").GetPtr() if isMC else rdf_pt.Histo1D(("h", " ", bins, -1, 1), "scale").GetPtr()
            dh      = ROOT.RooDataHist("dh", "dh", s, ROOT.RooFit.Import(h))
            fitData = dh

        # fit function
        # ------------
        if args.region == "EB" :
            mean    = ROOT.RooRealVar("mean",   " ", 0,     -0.05,  0.05)
            sigma   = ROOT.RooRealVar("sigma",  " ", 0.05,  1e-3,   1)
            alpha1  = ROOT.RooRealVar("alpha1", " ", 0.7,   1e-1,   3)
            alpha2  = ROOT.RooRealVar("alpha2", " ", 0.7,   1e-1,   3)
            n1      = ROOT.RooRealVar("n1",     " ", 80,    1,     200)
            n2      = ROOT.RooRealVar("n2",     " ", 80,    1,     200)
            width   = ROOT.RooRealVar("width",  " ", 0.02,  min_float,   2)
        else:
            mean    = ROOT.RooRealVar("mean",   " ", 0,     -0.05,  0.05)
            sigma   = ROOT.RooRealVar("sigma",  " ", 0.05,  1e-3,   1)
            alpha1  = ROOT.RooRealVar("alpha1", " ", 0.7,   1e-1,   3)
            alpha2  = ROOT.RooRealVar("alpha2", " ", 0.7,   1e-1,   3)
            n1      = ROOT.RooRealVar("n1",     " ", 80,    1,     200)
            n2      = ROOT.RooRealVar("n2",     " ", 80,    1,     200)
            width   = ROOT.RooRealVar("width",  " ", 0.02,  min_float,   2)         

        dcbPdf  = ROOT.RooCrystalBall("dcb", "dcb", s, mean, sigma, alpha1, n1, alpha2, n2)
        bw      = ROOT.RooBreitWigner("bw", "bw", s, mean, width)
        fitpdf  = ROOT.RooFFTConvPdf("dxb", "dcb (X) bw", s, dcbPdf, bw)

        # perform the fitting
        # -------------------
        s.setRange("fitRange", -0.8, 0.8)

        nll = fitpdf.createNLL(
            fitData,
            ROOT.RooFit.Range("fitRange"), ROOT.RooFit.SumW2Error(True)
        )
        minimizer = ROOT.RooMinimizer(nll)
        # minimizer.setStrategy(2)
        # minimizer.setEps(1e-3) # reduce the tolerance
        minimizer.setOffsetting(True)
        minimizer.optimizeConst(True)

        # algos = ["migrad", "minimize", "simplex", "fumili"] if args.region == "EB" else ["minimize", "migrad", "simplex", "fumili"]
        algos = ["migrad", "minimize", "simplex", "fumili"]
        for algo in algos:
            print("[INFO] Using {} algorithm to find the minimum".format(algo))

            status, status_hesse, status_minos = -999, -999, 0
            ntries, max_tries = 0, args.MaxTries
            while (status != 0 or status_hesse != 0 or status_minos != 0):
                if (ntries > 0):
                    print("[INFO] Fit failed! Refitting({}/{})...".format(ntries, max_tries))
                if (ntries >= max_tries):
                    break
                status = minimizer.minimize("Minuit", algo)
                status_hesse = minimizer.hesse()
                # status_minos = minimizer.minos()
                result = minimizer.save()
                ntries += 1
                
            if (not (status != 0 or status_hesse != 0 or status_minos != 0)):
                print(color.RED + "Fit succeed using {} algorithm!".format(algo) + color.END)
                break

        result.Print("V")
        print(result.status())
        print(result.covQual())

        # extract scale from the mean
        # ---------------------------
        mean_par = result.floatParsFinal().find("mean")
        scale_offset = mean_par.getValV()
        scale_err = mean_par.getError()

        # draw the fitting results
        # ------------------------
        if args.unBinFit:
            dh = dset.binnedClone()

        DrawFitting(s, dh, fitpdf, result, isMC, i, bins)
        print(color.CYAN + "[INFO] scale offset = {:8.4f} +/- {:8.4f}\n".format(scale_offset, scale_err) + color.END, flush=True)
        scaleList.append(scale_offset)
        errList.append(scale_err)

    return scaleList, errList


if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

    # turn on the silent mode of RooFit
    if (not args.verbose):
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
        ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
        # ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Minimization)
        # ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.NumIntegration)
        # ROOT.RooMsgService.instance().getStream(1).removeTopic(ROOT.RooFit.Eval)
        ROOT.RooMsgService.instance().setSilentMode(True)

    print(color.BLUE + color.BOLD + "Energy scale difference between data and mc" + color.END, flush=True)
    print("[INFO] Arguments: ", flush=True)
    print(json.dumps(vars(args), indent=4), flush=True)

    # specify input minitrees
    fDa = [
        "../miniTree/UL2016preVFP/miniTree_Data_UL2016preVFP.root",
        "../miniTree/UL2016postVFP/miniTree_Data_UL2016postVFP.root",
        "../miniTree/UL2017/miniTree_Data_UL2017.root",
        "../miniTree/UL2018/miniTree_Data_UL2018.root"
    ]
    fZg = [
        "../miniTree/UL2016preVFP/miniTree_ZGToLLG_UL2016preVFP.root",
        "../miniTree/UL2016postVFP/miniTree_ZGToLLG_UL2016postVFP.root",
        "../miniTree/UL2017/miniTree_ZGToLLG_UL2017.root",
        "../miniTree/UL2018/miniTree_ZGToLLG_UL2018.root"
    ]
    print("[INFO] Input trees of data: ", flush=True)
    print(json.dumps(fDa, indent=4), flush=True)
    print("[INFO] Input trees of mc: ", flush=True)
    print(json.dumps(fZg, indent=4), flush=True)

    # Et bins and pre-selections
    bin_array = np.asarray(args.bin_array.split(","), "d")
    print("[INFO] Corrections in bins of : {}\n".format(bin_array))
    cut = "isEBPho_Lead" if args.region == "EB" else "isEEPho_Lead"
    sel_basic = " && ".join([
        "Z_reg.M() > 80 && Z_reg.M() < 100",
        "eleHDALRegPt_Lead > 25",
        cut
    ])
    print("[INFO] Pre selections: {}".format(sel_basic), flush=True)

    # perform the energy scale correction
    # ROOT.EnableImplicitMT(10)
    rdf_Zg = ROOT.RDataFrame("miniTree", fZg).Filter(sel_basic)
    rdf_Da = ROOT.RDataFrame("miniTree", fDa).Filter(sel_basic)

    sDict, eDict = {}, {}
    options = ["EB_MC", "EB_Data"] if args.region == "EB" else ["EE_MC", "EE_Data"]
    for option in options:
        isMC = True if "Data" not in option else False

        print(color.GREEN + "Extract {} scale offset from fit".format(option) + color.END, flush=True)
        sList, eList = main(isMC)
        sDict[option] = sList
        eDict[option] = eList

    plotScale()
