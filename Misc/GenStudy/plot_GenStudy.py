import ROOT
import os
from interface.CMS_lumi import CMS_lumi


def DrawFrac(h_Mee1, h_Mee2, h_Mee3, h_total):
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

    # fraction setting
    efferr_Mee1.GetXaxis().SetTitle("M_{ee} [GeV]")
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

    colorN = ["#289672", "#E1701A", "#125D98"]
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

    # legend
    l = ROOT.TLegend(0.56, 0.72, 0.8, 0.9)
    l.SetTextSize(0.04)
    # l.SetNColumns(2)
    l.AddEntry(efferr_Mee1, "#DeltaR(e, e) < 0.4", "lep")
    l.AddEntry(efferr_Mee2, "0.4 #leq #DeltaR(e, e) < 1", "lep")
    l.AddEntry(efferr_Mee3, "#DeltaR(e, e) \geq 1", "lep")
    l.SetFillColor(0)
    l.SetLineColorAlpha(0, 0)
    l.SetFillStyle(0)
    l.Draw()

    CMS_lumi(c, 4, 10, "", 2017, True, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")

    os.makedirs("./plots/GenStudy", exist_ok = True)
    outName = "./plots/GenStudy/Mee_dR_frac.pdf"
    c.SaveAs(outName)
    c.Close()


def DrawStack(h_Mee1, h_Mee2, h_Mee3):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    c = ROOT.TCanvas("c","c",800,800)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.07)
    c.SetLeftMargin(0.14)
    c.SetBottomMargin(0.15)
    c.SetLogy()
    c.cd()

    fill_color = ["#9ad3bc", "#F8CE68", "#8ab6d6"]
    line_color = ["#202020", "#202020", "#202020"]
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
    
    hs.Draw("hist")

    ymax = h_Mee1.GetBinContent(h_Mee1.GetMaximumBin()) * 3
    hs.GetXaxis().SetTitle("M_{ee} [GeV]")
    hs.GetXaxis().SetTickSize(0.02)
    hs.GetXaxis().SetTitleSize(0.04)
    hs.GetXaxis().SetLabelSize(0.04)
    hs.GetXaxis().SetLabelOffset(0.02)
    hs.GetXaxis().SetTitleOffset(1.4)
    hs.GetYaxis().SetTitle("Events / 2.0 [GeV]")
    hs.GetXaxis().SetRangeUser(0, 50)
    hs.GetYaxis().SetNdivisions(506)
    hs.SetMaximum(ymax)
    hs.GetYaxis().SetTickSize(0.02)
    hs.GetYaxis().SetTitleSize(0.04)
    hs.GetYaxis().SetLabelSize(0.04)
    c.Modified()

    legend = ROOT.TLegend(0.58, 0.72, 0.8, 0.9)
    legend.SetTextSize(0.04)
    legend.AddEntry(h_Mee1, "#DeltaR(e, e) < 0.4", "f")
    legend.AddEntry(h_Mee2, "0.4 #leq #DeltaR(e, e) < 1", "f")
    legend.AddEntry(h_Mee3, "#DeltaR(e, e) \geq 1", "f")
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    legend.Draw()
    CMS_lumi(c, 4, 10, "", 2017, True, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "")
    c.RedrawAxis()

    os.makedirs("./plots/GenStudy", exist_ok = True)
    outName = "./plots/GenStudy/Mee_dR.pdf"
    c.SaveAs(outName)
    c.Close()


def main():
    h_Mee1 = ROOT.TH1F("h_Mee1", "", 25, 0, 50)
    h_Mee2 = ROOT.TH1F("h_Mee2", "", 25, 0, 50)
    h_Mee3 = ROOT.TH1F("h_Mee3", "", 25, 0, 50)

    for ev in range(tree.GetEntries()):
        tree.GetEntry(ev)

        # weight 
        mcwei       = tree.__getattr__("mcwei")
        genwei      = tree.__getattr__("genwei")

        # Gen variables
        mcPt_lep1       = tree.__getattr__("mcPt_lep1")
        mcEta_lep1      = tree.__getattr__("mcEta_lep1")
        mcPhi_lep1      = tree.__getattr__("mcPhi_lep1")
        mcPt_lep2       = tree.__getattr__("mcPt_lep2")
        mcEta_lep2      = tree.__getattr__("mcEta_lep2")
        mcPhi_lep2      = tree.__getattr__("mcPhi_lep2")

        mc_lep1, mc_lep2, mc_lep = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()
        mc_lep1.SetPtEtaPhiM(mcPt_lep1, mcEta_lep1, mcPhi_lep1, 0.000511)
        mc_lep2.SetPtEtaPhiM(mcPt_lep2, mcEta_lep2, mcPhi_lep2, 0.000511)
        mc_lep = mc_lep1 + mc_lep2

        w = mcwei * genwei
        if (mc_lep1.DeltaR(mc_lep2) < 0.4):
            h_Mee1.Fill(mc_lep.M(), w)
        elif ((mc_lep1.DeltaR(mc_lep2) >= 0.4) and (mc_lep1.DeltaR(mc_lep2) < 1)):
            h_Mee2.Fill(mc_lep.M(), w)
        elif (mc_lep1.DeltaR(mc_lep2) >= 1):
            h_Mee3.Fill(mc_lep.M(), w)

    # print(h_Mee1.GetEntries())
    # print(h_Mee2.GetEntries())
    # print(h_Mee3.GetEntries())
    DrawStack(h_Mee1, h_Mee2, h_Mee3)

    h_total = h_Mee1.Clone()
    h_total.Add(h_Mee2, 1)
    h_total.Add(h_Mee3, 1)
    DrawFrac(h_Mee1, h_Mee2, h_Mee3, h_total)


if __name__ == "__main__" :

    ROOT.gROOT.SetBatch() # PyROOT does not display any graphics(root "-b" option)
    infile_list = [
        "./minitree/2017/Minitree_HDalitz_ggF_eeg_m125_2017_RECO.root",
        "./minitree/2017/Minitree_HDalitz_VBF_eeg_m125_2017_RECO.root",
        "./minitree/2017/Minitree_HDalitz_WH_eeg_m125_2017_RECO.root",
        "./minitree/2017/Minitree_HDalitz_ZH_eeg_m125_2017_RECO.root"
    ]

    tree = ROOT.TChain("outTree")
    for i in infile_list:
        tree.Add(i)

    main()