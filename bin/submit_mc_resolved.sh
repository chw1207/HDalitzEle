#!/bin/bash

cd /home/chenghan/HDalitzEle/bin
source /opt/conda/etc/profile.d/conda.sh
conda activate hdalitzAna

echo "Starting job on " `date`
echo "Running on: `uname -a`"
echo "System software: `cat /etc/redhat-release`" 

./ResolvedAnalysis --config ../config_check/UL2016preVFP_SignalMC_Resolved.yaml
./ResolvedAnalysis --config ../config_check/UL2016postVFP_SignalMC_Resolved.yaml
./ResolvedAnalysis --config ../config_check/UL2017_SignalMC_Resolved.yaml
./ResolvedAnalysis --config ../config_check/UL2018_SignalMC_Resolved.yaml
