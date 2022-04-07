import ROOT
import os, sys
import time
from glob import glob
from datetime import datetime
from argparse import ArgumentParser
from plugins.colorPrint import *
import plugins.SampleConfig as cfg


def get_parser():
    parser = ArgumentParser(description = "python script to produce mini trees")
    parser.add_argument(
        "-r", "--run", 
        help = "samples to run [ Data | ZGToLLG | TTJets ]", 
        type = str
    )
    parser.add_argument(
        "-e", "--era", 
        help = "era to run [ 2016_preVFP | 2016_postVFP | 2017 | 2018 ], (default = 2017)", 
        default = "2017",
        type = str
    )
    return parser


class Analysis():
    #____________________________________________________________________________________
    def __init__(self, _sample, _era):
        self.sample = _sample
        self.era = _era

        # make sure the input arguments are correct
        sample_list = ["Data", "ZGToLLG", "TTJets"]
        era_list = ["2016_preVFP", "2016_postVFP", "2017", "2018"]
        if (self.sample not in sample_list) or (self.era not in era_list):
            print("[ERROR] Please specify the correct sample or era list!")
            parser.print_help()
            sys.exit(1)

        # directory for the nominal minitree
        self.year = int(str(self.era)[:4])
        self.outdir = "./miniTree/{}".format(self.era)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # hadd target file 
        self.haddfileName = " "

    #____________________________________________________________________________________
    def runAna(self):
        isMC = False if self.sample == "Data" else True
        Era = "{}_{}".format(self.sample, self.era)
        inpath = cfg.MCSample[Era]["outpath"] if isMC else cfg.DataSample[Era]["outpath"]
        
        for i in range(len(inpath)):
            Run = ""
            if (isMC == True):
                if (len(cfg.MCSample[Era]["production"]) == 1):
                    Run = self.era
                else:
                    Run = cfg.MCSample[Era]["production"][i]
            else:
                Run = cfg.DataSample[Era]["run"][i]
            
            ROOT.rdfZmumug("{}/skim.root".format(inpath[i]), "{}/miniTree_{}_{}.root".format(self.outdir, self.sample, Run), self.year, self.era, isMC)
    
    #____________________________________________________________________________________
    # merge the files of data of each run to one single file
    def haddFiles(self):
        if self.sample != "Data":
            print("[WARN] hadd command is useless!")
            return False

        # make sure the hadd target file is not in the glob file list
        self.haddfileName = "{}/miniTree_{}_{}.root".format(self.outdir, self.sample, self.era)
        Filelist = sorted(glob("{}/miniTree_{}_*.root".format(self.outdir, self.sample)))
        if self.haddfileName in Filelist:
            Filelist.remove(self.haddfileName)

        filestr = ""
        for f in Filelist:
            filestr += " {}".format(f)

        os.system("hadd -f {}{}".format(self.haddfileName, filestr)) 


def main():
    ana = Analysis(args.run, args.era)
    ana.runAna()
    if args.run == "Data":
        ana.haddFiles()
    

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    
    ROOT.gROOT.SetBatch()
    ROOT.gSystem.AddIncludePath("-Iexternal")
    ROOT.gSystem.SetBuildDir("tmpdir", ROOT.kTRUE)
    ROOT.gROOT.ProcessLine(".L rdfZmumug.C+")
    
    start_time = time.time()  
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush = True)
    print(color.BLUE + "---Start to produce mini trees!---" + color.END, flush = True)
    
    main()
    
    print(color.BLUE + "---All done!---" + color.END, flush = True)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))
    