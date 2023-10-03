from glob import glob
from ROOT import gSystem

file_list = glob("../miniTree/*/miniTree_HDalitz_*_eeg_*.root")
cmd = "hadd -f ../reg_signal.root "
for f in file_list:
    cmd += f + " "
gSystem.Exec(cmd)
