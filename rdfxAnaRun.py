import ROOT
import time
import os, sys
from glob import glob
from argparse import ArgumentParser
from datetime import datetime
from pluginsV2.colorPrint import *
import pluginsV2.SampleConfig as cfg


def get_parser():
    parser = ArgumentParser(description = "python script to run the rdfxAna")
    parser.add_argument(
        "-r", "--run", 
        help = "samples to run [ test | Data | HDalitz]", 
        type = str
    )
    parser.add_argument(
        "-y", "--year", 
        help = "year to run [ 2016 | 2017 | 2018 | all ], (default = 2017)", 
        default = "2017", 
        type = str
    )
    parser.add_argument(
        "-d", "--doSys", 
        help = "do systematics running or not [True | False]", 
        default = "False", 
        type = str,
    )
    parser.add_argument(
        "-n", "--NCPUs", 
        help = "number of cores", 
        default = -1, 
        type = int
    )
    return parser


class Analysis():
    #____________________________________________________________________________________
    def __init__(self, _sample, _era, _ncpu):
        self.sample = _sample
        self.era = _era
        self.ncpu = _ncpu

        # make sure the input arguments are correct
        sample_list = ["test", "Data", "HDalitz"]
        era_list = ["2016_preVFP", "2016_postVFP", "2017", "2018"]
        if (self.sample not in sample_list) or (self.era not in era_list):
            print("[ERROR] Please specify the correct sample or era list!")
            parser.print_help()
            sys.exit(1)

        # directory for the nominal minitree
        self.year = int(str(self.era)[:4])
        self.outdir = "./miniTree/{}".format(self.year)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

        # directory for the uncertainty minitree
        # use to estimate the systematic uncertainty
        self.sys_outdir = "./miniTree/uncertainty/{}".format(self.year)
        if not os.path.exists(self.sys_outdir):
            os.makedirs(self.sys_outdir)

        # hadd target file 
        self.haddfileName = " "

    #____________________________________________________________________________________
    def runAna(self, doSys):
        if ((self.sample not in ["HDalitz", "test"]) and (doSys == True)):
            print("[WARN] doSys option is usless to {}. It will be omitted!!!".format(self.sample))

        sysOption = ["UnPhoR9", "UnJECUp", "UnJECDo", "UnJERUp", "UnJERDo"]
        isMC = False if self.sample == "Data" else True
        
        if (self.sample != "test"):
            Era = "{}_{}".format(self.sample, self.era)
            inpath = cfg.MCSample[Era]["outpath"] if isMC else cfg.DataSample[Era]["outpath"]
            
            for i in range(len(inpath)):
                Run = cfg.MCSample[Era]["production"][i] if isMC else cfg.DataSample[Era]["run"][i]
                ROOT.rdfxAna("{}/skim.root".format(inpath[i]), "{}/miniTree_{}_{}.root".format(self.outdir, self.sample, Run), self.year, self.era, self.ncpu, isMC, "Nominal")

                if ((doSys == True) and (self.sample == "HDalitz")):
                    for sys in sysOption:
                        print("Process uncertainty miniTree: {}".format(sys), flush = True)
                        ROOT.rdfxAna("{}/skim.root".format(inpath[i]), "{}/miniTree_{}_{}.root".format(self.sys_outdir, self.sample, Run), self.year, self.era, self.ncpu, isMC, sys)

        else:
            ROOT.rdfxAna("/data4/chenghan/test/skim2.root", "./test/test.root".format(self.era), self.year, self.era, self.ncpu, isMC, "Nominal")
    
    #____________________________________________________________________________________
    # merge the files of data of each run to one single file
    def haddFiles(self):
        if self.sample != "Data":
            print("[WARN] hadd command is useless!")
            return False

        # make sure the hadd target file is not in the glob file list
        self.haddfileName = "{}/miniTree_{}_{}.root".format(self.outdir, self.sample, self.year)
        Filelist = sorted(glob("{}/miniTree_{}_*.root".format(self.outdir, self.sample)))
        if self.haddfileName in Filelist:
            Filelist.remove(self.haddfileName)

        filestr = ""
        for f in Filelist:
            filestr += " {}".format(f)

        os.system("hadd -f {}{}".format(self.haddfileName, filestr)) 


def str2bool(string):
    return string.lower() in ("yes", "true", "t", "1")


def main():
    doSysOption = True if (str2bool(args.doSys) == True) else False

    if (args.year == "2016"):
        for era in ["2016_preVFP", "2016_postVFP"]:
            ana = Analysis(args.run, era, args.NCPUs)
            ana.runAna(doSys = doSysOption)
            # hadd is only needed when both of "2016_preVFP" and "2016_postVFP" are done
            if ((args.run == "Data") and (era == "2016_postVFP")):
                ana.haddFiles()
    else:
        ana = Analysis(args.run, args.year, args.NCPUs)
        ana.runAna(doSys = doSysOption)
        if args.run == "Data":
            ana.haddFiles()


if __name__ == "__main__" :   
    parser = get_parser()
    args = parser.parse_args()

    ROOT.gROOT.SetBatch() # PyROOT does not display any graphics(root "-b" option)
    ROOT.gSystem.AddIncludePath("-Iexternal")
    ROOT.gSystem.SetBuildDir("tmpdir", ROOT.kTRUE)
    ROOT.gInterpreter.ProcessLine(".L rdfxAna.C++ ")
    print("", flush = True)
    
    start_time = time.time()
    
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush = True)
    print(color.BLUE + "---Start to process the skimmed ntuple!---" + color.END, flush = True)

    main()

    print(color.BLUE + "---All done!---" + color.END, flush = True)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))