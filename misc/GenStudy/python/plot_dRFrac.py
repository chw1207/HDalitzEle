import ROOT
import os
from CMS_lumi import CMS_lumi
from glob import glob

# script to visualize the dielectron mass distribution in different dR(e, e) regions
# dR(e, e) < 0.4, 0.4 <= dR(e, e) < 1, dR(e, e) >= 1
def DrawFrac(h_Mee1, h_Mee2, h_Mee3, h_total, legs, name, h_Mee4=None):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas('c', '', 800, 800)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.07)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetBottomMargin(0.13)
    c.cd( )
    TickSize = 0.02
    AxisTitleSize = 0.05
    AxisLabelSize = 0.045

    efferr_Mee1, efferr_Mee2, efferr_Mee3 = ROOT.TGraphAsymmErrors(), ROOT.TGraphAsymmErrors(), ROOT.TGraphAsymmErrors()
    efferr_Mee1.BayesDivide(h_Mee1, h_total)
    efferr_Mee2.BayesDivide(h_Mee2, h_total)
    efferr_Mee3.BayesDivide(h_Mee3, h_total)
    if h_Mee4 != None:
        efferr_Mee4 = ROOT.TGraphAsymmErrors()
        efferr_Mee4.BayesDivide(h_Mee4, h_total)

    # fraction setting
    efferr_Mee1.GetXaxis().SetTitle("m_{ee} [GeV]")
    efferr_Mee1.GetYaxis().SetTitle("Fraction")
    efferr_Mee1.GetXaxis().SetRangeUser(0, 50)
    efferr_Mee1.GetYaxis().SetRangeUser(0, 1.5)
    efferr_Mee1.GetXaxis().SetTickSize(TickSize)
    efferr_Mee1.GetXaxis().SetTitleSize(AxisTitleSize)
    efferr_Mee1.GetXaxis().SetLabelSize(AxisLabelSize)
    efferr_Mee1.GetYaxis().SetTickSize(TickSize)
    efferr_Mee1.GetYaxis().SetTitleSize(AxisTitleSize)
    efferr_Mee1.GetYaxis().SetLabelSize(AxisLabelSize)
    efferr_Mee1.GetXaxis().SetTitleOffset(1.1)
    efferr_Mee1.GetYaxis().SetTitleOffset(1.3)

    colorN = ["#289672", "#E1701A", "#125D98", "#E76161"]
    efferr_Mee1.SetMarkerColor(ROOT.TColor.GetColor(colorN[0]))
    efferr_Mee1.SetMarkerSize(1.5)
    efferr_Mee1.SetMarkerStyle(20)
    efferr_Mee1.SetLineColor(ROOT.TColor.GetColor(colorN[0]))
    efferr_Mee1.SetLineWidth(2)
    efferr_Mee1.Draw("AP")

    efferr_Mee2.SetMarkerColor(ROOT.TColor.GetColor(colorN[1]))
    efferr_Mee2.SetMarkerSize(1.5)
    efferr_Mee2.SetMarkerStyle(20)
    efferr_Mee2.SetLineColor(ROOT.TColor.GetColor(colorN[1]))
    efferr_Mee2.SetLineWidth(2)
    efferr_Mee2.Draw("P same")

    efferr_Mee3.SetMarkerColor(ROOT.TColor.GetColor(colorN[2]))
    efferr_Mee3.SetMarkerSize(1.5)
    efferr_Mee3.SetMarkerStyle(20)
    efferr_Mee3.SetLineColor(ROOT.TColor.GetColor(colorN[2]))
    efferr_Mee3.SetLineWidth(2)
    efferr_Mee3.Draw("P same")
    
    if h_Mee4 != None:
        efferr_Mee4.SetMarkerColor(ROOT.TColor.GetColor(colorN[3]))
        efferr_Mee4.SetMarkerSize(1.5)
        efferr_Mee4.SetMarkerStyle(20)
        efferr_Mee4.SetLineColor(ROOT.TColor.GetColor(colorN[3]))
        efferr_Mee4.SetLineWidth(2)
        efferr_Mee4.Draw("P same")

    # legend
    l = ROOT.TLegend(0.38, 0.72, 0.6, 0.9)
    l.SetTextSize(0.04)
    # l.SetNColumns(2)
    l.AddEntry(efferr_Mee1, legs[0], "lep")
    l.AddEntry(efferr_Mee2, legs[1], "lep")
    l.AddEntry(efferr_Mee3, legs[2], "lep")
    if h_Mee4 != None:
        l.AddEntry(efferr_Mee4, legs[3], "lep")
    
    l.SetFillColor(0)
    l.SetLineColorAlpha(0, 0)
    l.SetFillStyle(0)
    l.Draw()

    c.RedrawAxis()
    CMS_lumi(c, 4, 10, "", 2017, True, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")

    c.SaveAs(name)
    c.Close()
    

def DrawStack(h_Mee1, h_Mee2, h_Mee3, legs, name, h_Mee4=None):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas("c", "c", 800, 800)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.14)
    c.SetBottomMargin(0.15)
    c.SetLogy()
    c.cd()

    fill_color = ["#9ad3bc", "#F8CE68", "#8ab6d6", "#E76161"]
    line_color = ["#202020", "#202020", "#202020", "#202020"]
    hs = ROOT.THStack("MeedR","")

    h_Mee1.SetLineColor(ROOT.TColor.GetColor(line_color[0]))
    # h_Mee1.SetLineWidth(2)
    h_Mee1.SetFillColor(ROOT.TColor.GetColor(fill_color[0]))
    hs.Add(h_Mee1)

    # h_Mee2.SetLineColor(ROOT.TColor.GetColor(line_color[1]))
    # # h_Mee2.SetLineWidth(2)
    # h_Mee2.SetFillColor(ROOT.TColor.GetColor(fill_color[1]))
    # hs.Add(h_Mee2)

    h_Mee3.SetLineColor(ROOT.TColor.GetColor(line_color[2]))
    # h_Mee3.SetLineWidth(2)
    h_Mee3.SetFillColor(ROOT.TColor.GetColor(fill_color[2]))
    hs.Add(h_Mee3)
    
    if h_Mee4 != None:
        h_Mee4.SetLineColor(ROOT.TColor.GetColor(line_color[3]))
        # h_Mee3.SetLineWidth(2)
        h_Mee4.SetFillColor(ROOT.TColor.GetColor(fill_color[3]))
        hs.Add(h_Mee4)
    
    hs.Draw("hist")

    ymax = h_Mee1.GetBinContent(h_Mee1.GetMaximumBin()) * 20
    hs.GetXaxis().SetTitle("m_{ee} (GeV)")
    hs.GetXaxis().SetTickSize(0.02)
    hs.GetXaxis().SetTitleSize(0.04)
    hs.GetXaxis().SetLabelSize(0.04)
    hs.GetXaxis().SetLabelOffset(0.02)
    hs.GetXaxis().SetTitleOffset(1.4)
    hs.GetYaxis().SetTitle("Events")
    hs.GetXaxis().SetRangeUser(0, 50)
    hs.GetYaxis().SetNdivisions(506)
    hs.SetMaximum(ymax)
    hs.GetYaxis().SetTickSize(0.02)
    hs.GetYaxis().SetTitleSize(0.04)
    hs.GetYaxis().SetLabelSize(0.04)
    c.Modified()

    legend = ROOT.TLegend(0.38, 0.76, 0.6, 0.9)
    legend.SetTextSize(0.04)
    legend.AddEntry(h_Mee1, legs[0], "f")
    # legend.AddEntry(h_Mee2, legs[1], "f")
    legend.AddEntry(h_Mee3, legs[2], "f")
    if h_Mee4 != None:
        legend.AddEntry(h_Mee4, legs[3], "f")
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.Draw()
    CMS_lumi(c, 5, 10, "", 2017, True, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    c.RedrawAxis()

    c.Print(name)
    c.Close()
    

def DrawStackV2(h_Mee1, h_Mee2, h_Mee3, legs, name, h_Mee4=None):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas("c", "c", 800, 800)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.14)
    c.SetBottomMargin(0.15)
    c.SetLogy()
    c.cd()

    fill_color = ["#9ad3bc", "#F8CE68", "#8ab6d6", "#E76161"]
    line_color = ["#202020", "#202020", "#202020", "#202020"]
    hs = ROOT.THStack("MeedR","")

    h_Mee1.SetLineColor(ROOT.TColor.GetColor(line_color[0]))
    # h_Mee1.SetLineWidth(2)
    h_Mee1.SetFillColor(ROOT.TColor.GetColor(fill_color[0]))
    hs.Add(h_Mee1)

    h_Mee2.SetLineColor(ROOT.TColor.GetColor(line_color[1]))
    # h_Mee2.SetLineWidth(2)
    h_Mee2.SetFillColor(ROOT.TColor.GetColor(fill_color[1]))
    hs.Add(h_Mee2)

    h_Mee3.SetLineColor(ROOT.TColor.GetColor(line_color[2]))
    # h_Mee3.SetLineWidth(2)
    h_Mee3.SetFillColor(ROOT.TColor.GetColor(fill_color[2]))
    hs.Add(h_Mee3)
    
    if h_Mee4 != None:
        h_Mee4.SetLineColor(ROOT.TColor.GetColor(line_color[3]))
        # h_Mee3.SetLineWidth(2)
        h_Mee4.SetFillColor(ROOT.TColor.GetColor(fill_color[3]))
        hs.Add(h_Mee4)
    
    hs.Draw("hist")

    ymax = h_Mee1.GetBinContent(h_Mee1.GetMaximumBin()) * 20
    hs.GetXaxis().SetTitle("#Delta#phi(e,e) (GeV)")
    hs.GetXaxis().SetTickSize(0.02)
    hs.GetXaxis().SetTitleSize(0.04)
    hs.GetXaxis().SetLabelSize(0.04)
    hs.GetXaxis().SetLabelOffset(0.02)
    hs.GetXaxis().SetTitleOffset(1.4)
    hs.GetYaxis().SetTitle("Events")
    # hs.GetXaxis().SetRangeUser(0, 50)
    hs.GetXaxis().SetNdivisions(508)
    hs.GetYaxis().SetNdivisions(506)
    hs.SetMaximum(ymax)
    hs.GetYaxis().SetTickSize(0.02)
    hs.GetYaxis().SetTitleSize(0.04)
    hs.GetYaxis().SetLabelSize(0.04)
    c.Modified()

    legend = ROOT.TLegend(0.38, 0.72, 0.6, 0.9)
    legend.SetTextSize(0.04)
    legend.AddEntry(h_Mee1, legs[0], "f")
    legend.AddEntry(h_Mee2, legs[1], "f")
    legend.AddEntry(h_Mee3, legs[2], "f")
    if h_Mee4 != None:
        legend.AddEntry(h_Mee4, legs[3], "f")
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.Draw()
    CMS_lumi(c, 5, 10, "", 2017, True, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    c.RedrawAxis()

    c.Print(name)
    c.Close()
    
    

if __name__ == "__main__" :
    # PyROOT does not display any graphics(root "-b" option)
    ROOT.gROOT.SetBatch()
    
    infile_list = glob("../miniTree/*/miniTree_HDalitz_*_eeg_*.root")
    
    ROOT.EnableImplicitMT(20)
    rdf = ROOT.RDataFrame("miniTree", infile_list)
    rdf = rdf.Define("mcdr", "ROOT::VecOps::DeltaR(mcEta_Lead, mcEta_subLead, mcPhi_Lead, mcPhi_subLead)")\
             .Define("mcdeta", "(mcEta_Lead - mcEta_subLead)")\
             .Define("mcdphi", "ROOT::VecOps::DeltaPhi(mcPhi_Lead, mcPhi_subLead)")\
             .Define("mc1", "ROOT::Math::PtEtaPhiMVector v(mcPt_Lead, mcEta_Lead, mcPhi_Lead, 0.000511); return v;")\
             .Define("mc2", "ROOT::Math::PtEtaPhiMVector v(mcPt_subLead, mcEta_subLead, mcPhi_subLead, 0.000511); return v;")\
             .Define("mll", "(mc1+mc2).M()")
             
    # h1 = rdf.Filter("mcdr < 0.1").Histo1D(("h_Mee1", " ", 25, 0, 50), "mll", "wei").GetPtr()
    # h2 = rdf.Filter("mcdr >= 0.1 && mcdr < 0.5").Histo1D(("h_Mee2", " ", 25, 0, 50), "mll", "wei").GetPtr()
    # h3 = rdf.Filter("mcdr >= 0.5 && mcdr < 1").Histo1D(("h_Mee3", " ", 25, 0, 50), "mll", "wei").GetPtr()
    # h4 = rdf.Filter("mcdr >= 1").Histo1D(("h_Mee4", " ", 25, 0, 50), "mll", "wei").GetPtr() 
    # legs = ["#DeltaR(e, e) < 0.1", "0.1 #leq #DeltaR(e, e) < 0.5", "0.5 #leq #DeltaR(e, e) < 1", "#DeltaR(e, e) \geq 1"]
    # name = "../plots/GenStudy/Mee_dR.pdf"
    # os.makedirs("../plots/GenStudy", exist_ok=True)
    # DrawStack(h1, h2, h3, legs, name, h4)
    
    # h_total = h1.Clone()
    # h_total.Add(h2, 1)
    # h_total.Add(h3, 1)
    # h_total.Add(h4, 1)
    # name = "../plots/GenStudy/Mee_frac.pdf"
    # os.makedirs("../plots/GenStudy", exist_ok=True)
    # DrawFrac(h1, h2, h3, h_total, legs, name, h4)
    
    h_m2 = rdf.Filter("category == 2").Histo1D(("h_m2", " ", 50, 0, 50), "mll", "wei").GetPtr()
    h_m1 = rdf.Filter("category == 3").Histo1D(("h_m1", " ", 50, 0, 50), "mll", "wei").GetPtr()
    h_re = rdf.Filter("category == 1").Histo1D(("h_re", " ", 50, 0, 50), "mll", "wei").GetPtr()
    legs = ["merging pair w 2^{nd} GSF track", "merging pair w/o 2^{nd} GSF track", "well-separated pair"]
    name = "../plots/GenStudy/Mee_cat.pdf"
    os.makedirs("../plots/GenStudy", exist_ok=True)
    DrawStack(h_m2, h_m1, h_re, legs, name)
    
    # h_total = h_m2.Clone()
    # h_total.Add(h_m1, 1)
    # h_total.Add(h_re, 1)
    # name = "../plots/GenStudy/Mee_frac_cat.pdf"
    # os.makedirs("../plots/GenStudy", exist_ok=True)
    # DrawFrac(h_m2, h_m1, h_re, h_total, legs, name)
    
    # h_m2_dr = rdf.Filter("category == 2").Histo1D(("h_m2", " ", 50, 0, 0.1), "mcdr", "wei").GetPtr()
    # h_m1_dr = rdf.Filter("category == 3").Histo1D(("h_m1", " ", 50, 0, 0.1), "mcdr", "wei").GetPtr()
    # h_re_dr = rdf.Filter("category == 1").Histo1D(("h_re", " ", 50, 0, 0.1), "mcdr", "wei").GetPtr()
    # legs = ["merging pair w 2^{nd} GSF track", "merging pair w/o 2^{nd} GSF track", "well-separated pair"]
    # name = "../plots/GenStudy/dr_cat.pdf"
    # os.makedirs("../plots/GenStudy", exist_ok=True)
    # DrawStackV2(h_m2_dr, h_m1_dr, h_re_dr, legs, name)
    
    # h_m2_deta = rdf.Filter("category == 2").Histo1D(("h_m2", " ", 50, -0.1, 0.1), "mcdeta", "wei").GetPtr()
    # h_m1_deta = rdf.Filter("category == 3").Histo1D(("h_m1", " ", 50, -0.1, 0.1), "mcdeta", "wei").GetPtr()
    # h_re_deta = rdf.Filter("category == 1").Histo1D(("h_re", " ", 50, -0.1, 0.1), "mcdeta", "wei").GetPtr()
    # legs = ["merging pair w 2^{nd} GSF track", "merging pair w/o 2^{nd} GSF track", "well-separated pair"]
    # name = "../plots/GenStudy/deta_cat.pdf"
    # os.makedirs("../plots/GenStudy", exist_ok=True)
    # DrawStackV2(h_m2_deta, h_m1_deta, h_re_deta, legs, name)
    
    # h_m2_dphi = rdf.Filter("category == 2").Histo1D(("h_m2", " ", 50, -0.1, 0.1), "mcdphi", "wei").GetPtr()
    # h_m1_dphi = rdf.Filter("category == 3").Histo1D(("h_m1", " ", 50, -0.1, 0.1), "mcdphi", "wei").GetPtr()
    # h_re_dphi = rdf.Filter("category == 1").Histo1D(("h_re", " ", 50, -0.1, 0.1), "mcdphi", "wei").GetPtr()
    # legs = ["merging pair w 2^{nd} GSF track", "merging pair w/o 2^{nd} GSF track", "well-separated pair"]
    # name = "../plots/GenStudy/dphi_cat.pdf"
    # os.makedirs("../plots/GenStudy", exist_ok=True)
    # DrawStackV2(h_m2_dphi, h_m1_dphi, h_re_dphi, legs, name)