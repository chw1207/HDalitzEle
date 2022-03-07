import os, sys, subprocess
from argparse import ArgumentParser
from commonTools_HDalitz import color, catNumMap, procToDatacardNameMap, years, massBaseList, massList
import time

def get_parser():
    parser = ArgumentParser(description = "Script to run all the signal scripts")
    parser.add_argument("-df", "--doFitting", help = "do the fitting", default = 0, type = int)
    parser.add_argument("-di", "--doInterpolation", help = "do the interpolation", default = 0, type = int)
    parser.add_argument("-mf", "--makeFinalPlot", help = "make model plot", default = 0, type = int)
    parser.add_argument("-t",  "--doTestOnly", help = "only do the test", default = 0, type = int)
    
    return parser

start_time = time.time()

parser = get_parser()
args = parser.parse_args()
cats = sorted(catNumMap.keys())
procs = sorted(procToDatacardNameMap.keys())
wsName = "w" #! It should be consistant with the name used in tree2ws_HDalitz.py

########################################################
##                                                    ##
##    Perform the fitting --> simpleFit_HDalitz.py    ##
##                                                    ##
########################################################
if ((args.doFitting == 1) and (args.doTestOnly == 0)):
    print(color.BLUE + "Perform the fitting --> simpleFit_HDalitz.py" + color.END)
    print(color.BLUE + "------------------------------------------------------------------------------------------------" + color.END)
    for m in massBaseList:
        for y in years:
            for p in procs:
                for c in cats:
                    print(color.GREEN + "Fitting: mass points: {}, year: {}, process: {}, category: {}".format(m, y, p, c) + color.END)
                    
                    dset_name = "set_{}_{}".format(m, c)
                    path = "../Trees2WS/WS/2017/signal_{}_m{}.root".format(p, m)

                    try:
                        output = subprocess.check_output("python simpleFit_HDalitz.py -if {} -iw {} -id {} -proc {} -y {}".format(path, wsName, dset_name, p, y), shell = True, stderr = subprocess.STDOUT)
                        print(output)
                    except subprocess.CalledProcessError as error:
                        print(color.RED + "[ERROR] " + error.output + color.END)
                        sys.exit(1)

if((args.doFitting == 1) and (args.doTestOnly == 1)):
    print(color.BLUE + "Perform the fitting --> simpleFit_HDalitz.py" + color.END)
    print(color.BLUE + "------------------------------------------------------------------------------------------------" + color.END)
    print(color.GREEN + "Only run one of the category(Merged1Gsf_HVBF) to debug!" + color.END)  
    try:
        output = subprocess.check_output("python simpleFit_HDalitz.py -if ../Trees2WS/WS/2017/signal_ZH_m120.root -iw w -id set_120_Merged1Gsf_HVBF -proc ZH -y 2017", shell = True, stderr = subprocess.STDOUT) 
        print(output)
    except subprocess.CalledProcessError as error:
        print(color.RED + "[ERROR] " + error.output + color.END)
        sys.exit(1)


#######################################################################
##                                                                   ##
##    Do the interpolation --> simpleFit_HDalitz_Interpolation.py    ##
##                                                                   ##
#######################################################################
if ((args.doInterpolation == 1) and (args.doTestOnly == 0)):
    print(color.BLUE + "Do the interpolation --> simpleFit_HDalitz_Interpolation.py" + color.END)
    print(color.BLUE + "-----------------------------------------------------------------------------------------------" + color.END)
    for y in years:
            for p in procs:
                for c in cats:
                    print(color.GREEN + "Interpolation: year: {}, process: {}, category: {}".format(y, p, c) + color.END)

                    fr_120 = "./FitResults/{}/fit_120_{}_{}.root".format(y, c, p)
                    fr_125 = "./FitResults/{}/fit_125_{}_{}.root".format(y, c, p)
                    fr_130 = "./FitResults/{}/fit_130_{}_{}.root".format(y, c, p)

                    fw_120 = "./workspace/{}/workspace_HDalitz_sigMC_120_{}_{}.root".format(y, c, p)
                    fw_125 = "./workspace/{}/workspace_HDalitz_sigMC_125_{}_{}.root".format(y, c, p)
                    fw_130 = "./workspace/{}/workspace_HDalitz_sigMC_130_{}_{}.root".format(y, c, p)

                    try:
                        output = subprocess.check_output("python simpleFit_HDalitz_Interpolation.py -c {} -y {} -p {} -fr120 {} -fr125 {} -fr130 {} -fw120 {} -fw125 {} -fw130 {}".format(c, y, p, fr_120, fr_125, fr_130, fw_120, fw_125, fw_130), shell = True, stderr = subprocess.STDOUT) 
                        print(output)
                    except subprocess.CalledProcessError as error:
                        print(color.RED + "[ERROR] " + error.output + color.END)
                        sys.exit(1)

if((args.doInterpolation == 1) and (args.doTestOnly == 1)):
    print(color.BLUE + "Do the interpolation --> simpleFit_HDalitz_Interpolation.py" + color.END)
    print(color.BLUE + "-----------------------------------------------------------------------------------------------" + color.END)
    print(color.GREEN + "Only run one of the category(Merged2Gsf_EBHR9) to debug!" + color.END)
    fr_120 = "./FitResults/2017/fit_120_Merged2Gsf_EBHR9_ggF.root"
    fr_125 = "./FitResults/2017/fit_125_Merged2Gsf_EBHR9_ggF.root"
    fr_130 = "./FitResults/2017/fit_130_Merged2Gsf_EBHR9_ggF.root"

    fw_120 = "./workspace/2017/workspace_HDalitz_sigMC_120_Merged2Gsf_EBHR9_ggF.root"
    fw_125 = "./workspace/2017/workspace_HDalitz_sigMC_125_Merged2Gsf_EBHR9_ggF.root"
    fw_130 = "./workspace/2017/workspace_HDalitz_sigMC_130_Merged2Gsf_EBHR9_ggF.root"

    try:
        output = subprocess.check_output("python simpleFit_HDalitz_Interpolation.py -c Merged2Gsf_EBHR9 -y 2017 -p ggF -fr120 {} -fr125 {} -fr130 {} -fw120 {} -fw125 {} -fw130 {}".format(fr_120, fr_125, fr_130, fw_120, fw_125, fw_130), shell = True, stderr = subprocess.STDOUT) 
        print(output)
    except subprocess.CalledProcessError as error:
        print(color.RED + "[ERROR] " + error.output + color.END)
        sys.exit(1)



#########################################################
##                                                     ##
##   make the model plots --> makeFinalModelPlots.py   ##
##                                                     ##
#########################################################
if ((args.makeFinalPlot == 1) and (args.doTestOnly == 0)):
    print(color.BLUE + "Make the model plot --> makeFinalModelPlots_HDalitz.py" + color.END)
    print(color.BLUE + "-----------------------------------------------------------------------------------------------" + color.END)
    for p in procs:
        for c in cats:
            print(color.GREEN + "Final Model: mass: 125GeV, process: {}, category: {}".format(p, c) + color.END)

            try:
                output = subprocess.check_output("python makeFinalModelPlots_HDalitz.py -cat {} -proc {}".format(c, p), shell = True, stderr = subprocess.STDOUT) 
                print(output)
            except subprocess.CalledProcessError as error:
                print(color.RED + "[ERROR] " + error.output + color.END)
                sys.exit(1)

if((args.makeFinalPlot == 1) and (args.doTestOnly == 1)):
    try:
        output = subprocess.check_output("python makeFinalModelPlots_HDalitz.py -cat Merged2Gsf_EBHR9 -proc ggF", shell = True, stderr = subprocess.STDOUT) 
        print(output)
    except subprocess.CalledProcessError as error:
        print(color.RED + "[ERROR] " + error.output + color.END)
        sys.exit(1)















seconds = time.time() - start_time
print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))