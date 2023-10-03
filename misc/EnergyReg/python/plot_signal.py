import sys, os
import ROOT
from pprint import pprint
from glob import glob
from CMS_lumi import CMS_lumi

ROOT.gROOT.LoadMacro("../interface/tdrstyle.C")
ROOT.gROOT.ProcessLine("setTDRStyle();")
ROOT.gInterpreter.ProcessLine(""" #include "../interface/EnRegression_Signal.h" """)

def main(region):
    filter = "abs(eleSCEta_Lead) < 1.4442"
    if region == "EE":
        filter = "abs(eleSCEta_Lead) > 1.566 && abs(eleSCEta_Lead) < 2.5"
    
    # Merged-2Gsf
    ROOT.EnableImplicitMT(10)
    rdf = ROOT.RDataFrame("miniTree", sigfiles)\
              .Filter("elePresel_Lead == 1 && eleCalibPt_Lead > 25 && category == 2")\
              .Define("diTrkPt_Lead",       "(float) diTrk.Pt()")
    rdf_reg = ROOT.doEnRegression_XGB(ROOT.RDF.AsRNode(rdf))
        
    rdf_reg_re = rdf_reg.Filter(filter)\
                        .Define("RecoEleReg_Lead", "ROOT::Math::PtEtaPhiMVector v(eleHDALRegPt_Lead, RecoEle_Lead.Eta(), RecoEle_Lead.Phi(), RecoEle_Lead.M()); return v;")\
                        .Define("RecoEleOri_Lead", "ROOT::Math::PtEtaPhiMVector v(eleCalibPt_Lead, RecoEle_Lead.Eta(), RecoEle_Lead.Phi(), RecoEle_Lead.M()); return v;")\
                        .Define("higgsMass_NEW", "(RecoEleReg_Lead + RecoPho_Lead).M()")\
                        .Define("higgsMass_OLD", "(RecoEleOri_Lead + RecoPho_Lead).M()")
        
    bins = 60 if region == "EB" else 60
    h_new = rdf_reg_re.Histo1D(("h_new", "h_new", bins, 110, 140), "higgsMass_NEW", "wei").GetPtr()
    h = rdf_reg_re.Histo1D(("h", "h", bins, 110, 140), "higgsMass_OLD", "wei").GetPtr()
    
    h_new.Scale(1./h_new.Integral(-1, -1))
    h.Scale(1./h.Integral(-1, -1))

    c = ROOT.TCanvas("c", "", 800, 700)
    c.cd()
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.13)

    h.GetXaxis().SetTitle("M_{ee#gamma} [GeV]")
    h.GetYaxis().SetTitle("arb. unit")
    h.GetYaxis().SetRangeUser(0, h_new.GetBinContent(h_new.GetMaximumBin()) * 1.4)

    h.SetLineColor(ROOT.TColor.GetColor("#E16262"))
    h.SetLineWidth(3)

    h_new.SetLineColor(ROOT.TColor.GetColor("#3A9679"))
    h_new.SetLineWidth(3)

    h.Draw("hist")
    h_new.Draw("hist same")

    leg = ROOT.TLegend(0.5, 0.75, 0.89, 0.9)
    leg.SetTextSize(0.035)
    leg.AddEntry(h, "Before energy correction", "l")
    leg.AddEntry(h_new, "After energy correction", "l")
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.Draw()

    extral = "ECAL Barrel" if "EB" in region else "ECAL Endcap"
    ltx = ROOT.TLatex()
    ltx.SetNDC()
    ltx.SetTextFont(42)
    ltx.SetTextSize(0.04)
    ltx.DrawLatex(0.2, 0.73, extral)

    ltx.DrawLatex(0.6, 0.7, "M_{H} = 125 GeV")

    effsig_new = ROOT.effSigma(h_new)
    effsig = ROOT.effSigma(h)
    
    print("         :effective sigma before regression = {}".format(effsig))
    print("         :effective sigma after regression  = {}".format(effsig_new))
    ltx.DrawLatex(0.65, 0.65, "#sigma^{Before}_{eff} = %.3f GeV" %effsig)
    ltx.DrawLatex(0.65, 0.58, "#sigma^{After}_{eff} = %.3f GeV" %effsig_new)

    l = ROOT.TLine(125, 0, 125, h.GetBinContent(h.GetMaximumBin()) * 1.05)
    l.SetLineStyle(7)
    l.SetLineWidth(3)
    l.SetLineColor(ROOT.TColor.GetColor("#202020"))
    l.Draw()


    CMS_lumi(c, 5, 10, "138 fb^{-1}", 2017, True, "Simulation Preliminary", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")

    directory = "../plots/validation/{}".format(region)
    if not os.path.exists(directory):
        os.makedirs(directory)
    

    c.Print("{}/validation_reg_XGB_opt_{}.pdf".format(directory, region))
    
    c.Close()

if __name__ == "__main__":
    # setup the build directory
    ROOT.gROOT.SetBatch()
    
    sigfiles = glob("../../GenStudy/miniTree/*/miniTree_HDalitz_*_eeg_125_*.root") 
    pprint(sigfiles)
    
    print("ECAL Barrel Regression")
    main("EB")
    
    print("ECAL Endcap Regression")
    main("EE")