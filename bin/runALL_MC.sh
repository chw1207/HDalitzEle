#!/bin/sh

# ./MergedAnalysis --config ../config/UL2016preVFP_SignalMC_Merged.yaml   |& tee ../logger/UL2016preVFP_SignalMC_Merged.txt
# ./MergedAnalysis --config ../config/UL2016postVFP_SignalMC_Merged.yaml  |& tee ../logger/UL2016postVFP_SignalMC_Merged.txt
# ./MergedAnalysis --config ../config/UL2017_SignalMC_Merged.yaml         |& tee ../logger/UL2017_SignalMC_Merged.txt
# ./MergedAnalysis --config ../config/UL2018_SignalMC_Merged.yaml         |& tee ../logger/UL2018_SignalMC_Merged.txt

./ResolvedAnalysis --config ../config/UL2016preVFP_SignalMC_Resolved.yaml   |& tee ../logger/UL2016preVFP_SignalMC_Resolved.txt
./ResolvedAnalysis --config ../config/UL2016postVFP_SignalMC_Resolved.yaml  |& tee ../logger/UL2016postVFP_SignalMC_Resolved.txt
./ResolvedAnalysis --config ../config/UL2017_SignalMC_Resolved.yaml         |& tee ../logger/UL2017_SignalMC_Resolved.txt
./ResolvedAnalysis --config ../config/UL2018_SignalMC_Resolved.yaml         |& tee ../logger/UL2018_SignalMC_Resolved.txt