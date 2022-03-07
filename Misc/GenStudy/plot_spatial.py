import sys
import uproot as uproot4
import numpy as np
import pandas as pd
import ROOT
import mplhep as hep
import matplotlib.pyplot as plt
from interface.Utilities import deltaPhi, deltaEta
from interface.CMS_lumi import CMS_lumi

def fill_2DHist(h, Xarray, Yarray, wei):
        [h.Fill(Xarray[i], Yarray[i], wei[i]) for i in Xarray.index]

def Draw2Dhist(hM):
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetOptStat(0)

    TickSize = 0.02
    TitleSize = 0.05
    LabelSize = TitleSize * 0.92

    c = ROOT.TCanvas('c', '', 700, 700)
    c.cd()
    ROOT.gPad.SetRightMargin(0.14)
    ROOT.gPad.SetBottomMargin(0.12)
    ROOT.gPad.SetTopMargin(0.06)
    ROOT.gPad.SetLeftMargin(0.14)
    
    hM.Scale(1. / hM.Integral(-1, -1, -1, -1))
    hM.GetXaxis().SetTitle("#Delta#eta")
    hM.GetXaxis().SetTickSize(TickSize)
    hM.GetXaxis().SetTitleSize(TitleSize)
    hM.GetXaxis().SetLabelSize(LabelSize)
    hM.GetXaxis().SetTitleOffset(1.1)
    hM.GetXaxis().SetNdivisions(508)

    hM.GetYaxis().SetTitle("#Delta#phi")
    hM.GetYaxis().SetTickSize(TickSize)
    hM.GetYaxis().SetTitleSize(TitleSize)
    hM.GetYaxis().SetLabelSize(LabelSize)
    hM.GetYaxis().SetTitleOffset(1.4)
    hM.GetYaxis().SetNdivisions(508)

    hM.GetZaxis().SetLabelSize(0.04)
    hM.GetZaxis().SetRangeUser(0, 0.5)
    hM.GetZaxis().SetNdivisions(510)
    hM.Draw("COLZ")

    CMS_lumi(c, 4, 0, "41.5 fb^{-1}", 2017, True, "Simulation", "", "")

    c.RedrawAxis()
    c.Print("./plot/spatial.pdf")
    c.Close()


def main():
    data = pd.DataFrame()
    
    branches = ["mcEta_lep1", "mcEta_lep2", "mcPhi_lep1", "mcPhi_lep2", "mcwei", "genwei"]
    for i in uproot4.iterate("./minitree/2017/Minitree_HDalitz_*_eeg_m125_2017_RECO.root:outTree", branches, library = "pd"):
        i["deltaEta"] = np.absolute(np.vectorize(deltaEta)(i.mcEta_lep1, i.mcEta_lep2))
        i["deltaPhi"] = np.absolute(np.vectorize(deltaPhi)(i.mcPhi_lep1, i.mcPhi_lep2))
        i["wei"] = i.mcwei * i.genwei
        data = data.append(i)
    # print(data)
    hM_spatial = ROOT.TH2F("hM_spatial"," ", 40, 0, 0.4, 40, 0, 0.4)
    fill_2DHist(hM_spatial, data["deltaEta"], data["deltaPhi"], data["wei"])
    Draw2Dhist(hM_spatial)

if __name__ == "__main__":
    state = main()
    sys.exit(state) 
