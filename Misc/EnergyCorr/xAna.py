import ROOT # version: 6.24/02
import os, sys
import time 
from plugins.colorPrint import *


ROOT.gInterpreter.ProcessLine("""
# include "./plugins/muon_utilities.h"
# include "./plugins/electron_utilities.h"
# include "./plugins/gen_utilities.h"
""")


class Analysis():
    #____________________________________________________________________________________
    def __init__(self, inputlist, ncpu):
        self.inputlist = inputlist
        self.do_gen = False

        ROOT.EnableImplicitMT(ncpu)
        self.df = ROOT.RDataFrame("ggNtuplizer/EventTree", self.inputlist)

        print(color.GREEN + "Loding done!" + color.END)
        print ("[INFO] using {} cores to do the analysis".format(ncpu))

    #____________________________________________________________________________________
    def run_muon_selection(self):
        # muon types are defined here: CMSSW/DataFormats/MuonReco/interface/Muon.h
        # namely: GlobalMuon     = 1<<1
        #         TrackerMuon    = 1<<2
        #         StandAloneMuon = 1<<3
        #         CaloMuon       = 1<<4
        #         PFMuon         = 1<<5
        #         RPCMuon        = 1<<6
        print(color.GREEN + "Run muon selection!" + color.END)
        self.df = (
            self.df
        
            # Define the muon type
            .Define("isGlobalMuon",      "MuTypeVector(muType, 1)")
            .Define("isTrackerMuon",     "MuTypeVector(muType, 2)")
            .Define("isPFMuon",          "MuTypeVector(muType, 5)")
            .Define("isTrkHighPtMuon",   "isTrackerMuon && muStations > 1 && muPixelHits > 0 && muTrkLayers > 5 && abs(muD0) < 0.2 && abs(muDz) < 0.5 && (muBestTrkPtError / muBestTrkPt < 0.3)")

            # Define the HZZ ID and Iso critera
            .Define("isHZZLooseMuon",    "muPt < 200 && muSIP < 4 && abs(muD0) < 0.5 && abs(muDz) < 1 && muBestTrkType != 2 && (isGlobalMuon ||(isTrackerMuon && muStations > 0))")
            .Define("isHZZTightMuon",    "muPt >= 200 && isHZZLooseMuon && (isPFMuon || isTrkHighPtMuon)")
            .Define("muIso03",           "NoneZeroIso03(muPFChIso03, muPFNeuIso03, muPFPhoIso03, muPFPUIso03)")
            .Define("isHZZIsoMuon",      "(muIso03/muPt) < 0.35")
            .Define("isGoodMuon",        "muPt > 4 && abs(muEta) < 2.4 && (isHZZLooseMuon || isHZZTightMuon) && isHZZIsoMuon")

            # Leave the events with good muons    
            .Filter("Sum(isGoodMuon) >= 2",                                         "HZZID selection")
            .Filter("muPt[isGoodMuon][0] > 20 && muPt[isGoodMuon][1] > 10",         "HLT threshold")
            .Filter("(muCharge[isGoodMuon][0] * muCharge[isGoodMuon][1]) < 0",      "Opposite charge")

            # Define the 4 vector and their variables
            .Define("mu1",               "TLorentzVector v; v.SetPtEtaPhiM(muPt[isGoodMuon][0], muEta[isGoodMuon][0], muPhi[isGoodMuon][0], 105.658*0.001); return v;")
            .Define("mu2",               "TLorentzVector v; v.SetPtEtaPhiM(muPt[isGoodMuon][1], muEta[isGoodMuon][1], muPhi[isGoodMuon][1], 105.658*0.001); return v;")
            .Define("Zmumu",             "mu1 + mu2")
            .Define("ZmumuMass",         "Zmumu.M()")
            .Define("ZmumuPt",           "Zmumu.Pt()")
            .Define("dR_mu1mu2",         "mu1.DeltaR(mu2)")

            .Filter("ZmumuMass > 60 && ZmumuMass < 120",                            "Z mass window cut")
        )

    #____________________________________________________________________________________
    def run_electron_selection(self):
        print(color.GREEN + "Run electron selection!" + color.END)
        self.df = (
            self.df
        
            # Define the good electron
            .Define("isM1Ele",          "eleClass == 1")
            .Define("isHggElectron",    "HggPreSelection(rhoAll, nEle, eleSCEta, nPho, phoSCEta, phoPFChIso, phoPFPhoIso, phoTrkIsoHollowConeDR03, phoR9Full5x5, phoCalibEt, phoSigmaIEtaIEtaFull5x5, phoHoverE)")
            .Define("isGoodEle_EB",     "isHggElectron && isM1Ele && eleCalibPt > 7. && abs(eleSCEta) < 1.4442 && eleEcalDrivenSeed == 1")
            .Define("isGoodEle_EE",     "isHggElectron && isM1Ele && eleCalibPt > 7. && abs(eleSCEta) > 1.566 && abs(eleSCEta) < 2.5 && eleEcalDrivenSeed == 1")
            .Define("isGoodEle",        "isGoodEle_EB || isGoodEle_EE")

            # Leave the events with good electron   
            .Filter("Sum(isGoodEle) >= 1",   "Electron selection")

            # Define the 4 vector and their variables
            .Define("ele1",             "TLorentzVector v; v.SetPtEtaPhiM(eleCalibPt[isGoodEle][0], eleEta[isGoodEle][0], elePhi[isGoodEle][0], 0.511*0.001); return v;")
            .Define("ele1Pt",           "ele1.Pt()")
            .Define("ele1ZmumuPt",      "ele1Pt/ZmumuPt")
            .Define("dR_ele1mu1",       "ele1.DeltaR(mu1)")
            .Define("dR_ele1mu2",       "ele1.DeltaR(mu2)")
            .Define("ele1SCEta",        "eleSCEta[isGoodEle][0]")
            .Define("ele1Charge",       "eleCharge[isGoodEle][0]")
        )

    #____________________________________________________________________________________
    def print_cut_flow(self):
        print(color.GREEN + "Cut-flow report: " + color.END)
        self.df.Report().Print() 

    #____________________________________________________________________________________
    def do_gen_match(self, purityCalc = False):
        print(color.GREEN + "Do gen matched!" + color.END)
        self.do_gen = True

        self.df = (
            self.df
            .Define("genInd_ele1",      "GenMatchInd(ele1, nMC, mcPt, mcEta, mcPhi, mcMass)")
            .Define("isHardEle",        "GenType(mcStatusFlag[genInd_ele1], 0) == 1")
            .Define("isPromptEle",      "GenType(mcStatusFlag[genInd_ele1], 1) == 1")
            .Define("isTrueEle",        "abs(mcPID[genInd_ele1]) == 11 && abs(mcMomPID[genInd_ele1]) == 23 && isHardEle && isPromptEle")
        )

        if (purityCalc == True):
            a = self.df.Stats("isTrueEle", "wei").GetW()
            b = self.df.Filter("isTrueEle == 1", "true").Stats("isTrueEle", "wei").GetW()
            purity = round(b*100/a, 2)
            print("[INFO] Gen-match purity = {}".format(purity))

    #____________________________________________________________________________________
    def save_tree(self, treeName, outname):
        
        finalVars = ["mu1", "mu2", "ele1", "Zmumu", "dR_ele1mu1", "dR_ele1mu2", "dR_mu1mu2", "ele1SCEta", "ele1Charge", "wei"]

        if self.do_gen == True:
            for v in ["isHardEle", "isPromptEle", "isTrueEle"]:
                finalVars.append(v)

        print(color.GREEN + "Save skimmed tree in(takes time to execute the event loop...):  " + color.END)
        print("[INFO] Yield = {}".format(self.df.Stats("ZmumuMass", "wei").GetW()))
        print(outname)

        self.df.Snapshot(treeName, outname, finalVars)



def main():
    f = [
        "/data4/chenghan/mc/V10_02_10_07/job_fall17_ZZ/skim.root",
        "/data4/chenghan/mc/V10_02_10_07/job_autumn18_ZZ/skim.root",
        "/data4/chenghan/mc/V10_02_10_07/job_summer16_ZZ//skim.root"
    ] 

    ncpus = ncpus = os.cpu_count() - 2
    Ana = Analysis(f, ncpus)

    Ana.run_muon_selection()
    Ana.run_electron_selection()
    Ana.do_gen_match(purityCalc = True)
    Ana.print_cut_flow()

    outdir = "./miniTree"
    os.makedirs(outdir, exist_ok = True)
    Ana.save_tree("outTree", "{}/miniTree_ZZ.root".format(outdir))


if __name__ == "__main__" :
    start_time = time.time()

    main()

    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))

        