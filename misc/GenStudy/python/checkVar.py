import ROOT
import os, sys
import time
from CMS_lumi import CMS_lumi
from sigmaEff import sigmaEff
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
ROOT.gROOT.LoadMacro("../interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")


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


# def main(_track):
#     ROOT.gStyle.SetPadTickX(1)
#     ROOT.gStyle.SetPadTickY(1)
#     ROOT.gStyle.SetOptStat(0)

#     h1 = rdf_sigs.Histo1D(("h1", " ", 7, -0.5, 6.5), _track, "weight")
#     h2 = rdf_gjet.Histo1D(("h2", " ", 7, -0.5, 6.5), _track, "weight")
#     s1 = h1.GetBinContent(1)/np.sqrt(h2.GetBinContent(1))
#     s2 = (h1.GetBinContent(1)+h1.GetBinContent(2))/np.sqrt(h2.GetBinContent(1)+h2.GetBinContent(2))
#     print("significance for missing hit = 0: {}".format(s1))
#     print("significance for missing hit = 1: {}".format(s2))
    
#     h1.Scale(1/h1.Integral())
#     h2.Scale(1/h2.Integral())
#     c1 = ROOT.TCanvas("c1", "", 800, 800)
#     c1.cd()
    
#     if _track == "eleTrkMissHits_Lead":
#         xName = "Missing hits(main Gsf track)"
#     else:
#         xName = "Missing hits(attached Gsf track)"
#     Draw1DHist(c1, [h1, h2], [None, None], ["#E69D45", "#0061a8"], xaxis=xName, yaxis="A.U.", option="hist", Log=False)

#     leg = ROOT.TLegend(0.55, 0.7, 0.9, 0.9)
#     leg.SetHeader("Merged-2Gsf")
#     leg.SetTextSize(0.04)
#     leg.AddEntry(h1.GetPtr(), "Signal", "f")
#     leg.AddEntry(h2.GetPtr(), "#gamma converts to e", "f")
#     leg.SetFillColor(0)
#     leg.SetLineColor(0)
#     leg.Draw()
    
#     ltx = ROOT.TLatex()
#     ltx.SetNDC()
#     ltx.SetTextFont(42)
#     ltx.SetTextSize(0.035)
#     ltx.DrawLatex(0.5, 0.5, "#frac{S}{#sqrt{B}}(missing hits = 0) = %.4f" %(s1))
#     ltx.DrawLatex(0.5, 0.4, "#frac{S}{#sqrt{B}}(missing hits = 1) = %.4f" %(s2))

#     c1.RedrawAxis()
#     CMS_lumi(c1, 5, 10, "", 2017, True, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    
#     os.makedirs("../plots/GenStudy", exist_ok=True)
#     outName = "../plots/GenStudy/missHits_choice_{}.pdf".format(_track)
#     c1.Print(outName)
#     c1.Close()

if __name__ == "__main__" :
    start_time = time.time()
    
    # PyROOT does not display any graphics(root "-b" option)
    ROOT.gROOT.SetBatch()
    
    iBE = int(sys.argv[1])
    BE_str = "EB" if iBE == 0 else "EE"
    acc_filter = "abs(eleSCEta_Lead) < 1.4442" if (iBE == 0) else "abs(eleSCEta_Lead) > 1.566 && abs(eleSCEta_Lead) < 2.5"
    
    ROOT.EnableImplicitMT(10)
    infile_sigs = glob("../miniTree/*/miniTree_HDalitz_*_eeg_125_*.root")
    rdf_sigs = ROOT.RDataFrame("miniTree", infile_sigs)\
                   .Filter("elePresel_Lead == 1 && category == 2 && nGsfMatchToReco_Lead >= 2 && eleCalibPt_Lead > 25 && {}".format(acc_filter))\
                   .Define("eleE2x5OverE5x5",   "eleE2x5Full5x5_Lead/eleE5x5Full5x5_Lead")
    
    infile_gjet = glob("../miniTree/*/miniTree_GJets_*.root")
    infile_qcd = glob("../miniTree/*/miniTree_QCD_*.root")
    infile_DYJets = glob("../miniTree/*/miniTree_DYJets_*.root")
    rdf_bkgs = ROOT.RDataFrame("miniTree", infile_gjet+infile_qcd+infile_DYJets)\
                   .Filter("elePresel_Lead == 1 && nGsfMatchToReco_Lead >= 2 && eleCalibPt_Lead > 25 && {}".format(acc_filter))\
                   .Define("eleE2x5OverE5x5",   "eleE2x5Full5x5_Lead/eleE5x5Full5x5_Lead")         
    
    xbins = [
        np.arange(0.005, 0.2, 0.005),
        np.arange(0.005, 0.2, 0.005),
        np.arange(1, 10, 1),
        np.arange(0.5, 1, 0.01)
    ]
        
    font = {
        "color":  "darkred",
        "weight": "normal",
    }
    
    # vars = ["D0", "Dz", "SIP", "E2x5OverE5x5"]
    # # vars = ["D0"]
    # for i in range(len(vars)):
    #     ams = []
    #     for j in xbins[i]:
    #         filter_str = ""
    #         if vars[i] == "D0":
    #             filter_str = "eleTrkD0_Lead < {} && eleSubTrkD0_Lead < {}".format(j, j)
    #         if vars[i] == "Dz":
    #             filter_str = "eleTrkDz_Lead < {} && eleSubTrkDz_Lead < {}".format(j, j)
    #         if vars[i] == "SIP":
    #             filter_str = "eleSIP_Lead < {}".format(j)
    #         if vars[i] == "E2x5OverE5x5":
    #             filter_str = "eleE2x5OverE5x5 > {}".format(j)
    #         s = rdf_sigs.Filter(filter_str).Sum("wei")
    #         b = rdf_bkgs.Filter(filter_str).Sum("wei")
            
    #         AMS = s.GetValue()/ROOT.TMath.Sqrt(b.GetValue())
    #         ams.append(AMS)
            
    #     xmax = xbins[i][np.argmax(ams)]
    #     ymax = np.max(ams)
        
    #     plt.tight_layout()
    #     plt.scatter(xbins[i], ams, color="lightgreen", edgecolor="green") 
    #     plt.title("{} {} optimization".format(vars[i], BE_str), **font)
    #     plt.ylabel("S/$\sqrt{B}$", **font)
    #     plt.xlabel("cut value", **font) 
    #     plt.savefig("../plots/{}_{}.pdf".format(vars[i], BE_str))        
    #     print("Best cut value for {} at {} is {}".format(vars[i], BE_str, xbins[i][np.argmax(ams)])) 
    #     plt.close("all")      
    
    
    
    
    
    
    # for var in [""]
    
    
    
    
    
    
    # ams = []
    # xbins = np.arange(0.005, 0.2, 0.005)
    # # xbins = np.arange(1, 10, 1)
    # # xbins = np.arange(0.5, 1, 0.01)
    # for i in xbins:
    #     # s = rdf_sigs.Filter("eleE1x5OverE5x5 > {}".format(i, i)).Sum("wei")
    #     # b = rdf_bkgs.Filter("eleE1x5OverE5x5 > {}".format(i, i)).Sum("wei")
    #     # s = rdf_sigs.Filter("eleSIP_Lead < {}".format(i)).Sum("wei")
    #     # b = rdf_bkgs.Filter("eleSIP_Lead < {}".format(i)).Sum("wei")
        
    #     # s = rdf_sigs.Filter("eleTrkD0_Lead < {} && eleSubTrkD0_Lead < {}".format(i, i)).Sum("wei")
    #     # b = rdf_bkgs.Filter("eleTrkD0_Lead < {} && eleSubTrkD0_Lead < {}".format(i, i)).Sum("wei")
        
    #     s = rdf_sigs.Filter("eleTrkDz_Lead < {} && eleSubTrkDz_Lead < {}".format(i, i)).Sum("wei")
    #     b = rdf_bkgs.Filter("eleTrkDz_Lead < {} && eleSubTrkDz_Lead < {}".format(i, i)).Sum("wei")
        
    #     # AMS = ROOT.TMath.Sqrt(2.*((s.GetValue()+b.GetValue())*ROOT.TMath.Log(1.+s.GetValue()/b.GetValue())-s.GetValue()))
    #     AMS = s.GetValue()/ROOT.TMath.Sqrt(b.GetValue())
    #     ams.append(AMS)
  
    # plt.scatter(xbins, ams, color="lightgreen", edgecolor="green") 
    # plt.ylabel("S/$\sqrt{B}$")
    # plt.xlabel("cut value")
    # # plt.savefig("../plots/E1x5OverE5x5.pdf")    
    # plt.savefig("../plots/dz.pdf")    
    
    # print(xbins[np.argmax(ams)])       
    
    
    # for track in ["eleTrkMissHits_Lead", "eleSubTrkMissHits_Lead"]:    
    #     main(track)
    
    # h_dxy_tk1 = rdf_sigs.Histo1D(("h_eleTrkD0_Lead", " ", 40, 0, 0.2), "eleTrkD0_Lead", "wei")
    # h_dxy_tk2 = rdf_sigs.Histo1D(("h_eleSubTrkD0_Lead", " ", 40, 0, 0.2), "eleTrkD0_Lead", "wei")
    # h_dz_tk1 = rdf_sigs.Histo1D(("h_eleTrkDz_Lead", " ", 40, 0, 0.2), "eleTrkDz_Lead", "wei")
    # h_dz_tk2 = rdf_sigs.Histo1D(("h_eleSubTrkDz_Lead", " ", 40, 0, 0.2), "eleTrkDz_Lead", "wei")
    # h_sip = rdf_sigs.Histo1D(("h_eleSIP_Lead", " ", 10, 0, 10), "eleSIP_Lead", "wei")
    
    
    h_eleE2x5OverE5x5_sig = rdf_sigs.Histo1D(("h_eleE2x5OverE5x5_sig", " ", 60, 0.7, 1), "eleE2x5OverE5x5", "wei")
    h_eleE2x5OverE5x5_bkg = rdf_bkgs.Histo1D(("h_eleE2x5OverE5x5_bkg", " ", 60, 0.7, 1), "eleE2x5OverE5x5", "wei")
    
    
    c = ROOT.TCanvas("c", "c", 800, 900)
    # print(h_eleE2x5OverE5x5_sig.Integral())
    
    h_eleE2x5OverE5x5_sig.Scale(h_eleE2x5OverE5x5_bkg.Integral()/h_eleE2x5OverE5x5_sig.Integral())
    Draw1DHist(c, [h_eleE2x5OverE5x5_sig.GetPtr(), h_eleE2x5OverE5x5_bkg.GetPtr()], [None, None], ["#E69D45", "#0061a8"] , xaxis="E_{2x5}/E_{5x5}", yaxis="Event", option="hist", Log=False)

    c.RedrawAxis()
    CMS_lumi(c, 5, 10, "", 2017, True, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    
    # os.makedirs("../plots/GenStudy", exist_ok=True)
    # outName = "../plots/GenStudy/missHits_choice_{}.pdf".format(_track)
    c.Print("testEE.png")
    c.Close()

    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))