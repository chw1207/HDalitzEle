#!/bin/sh

# ./MergedAnalysis --config ../config_check/UL2016preVFP_SignalMC_Merged.yaml   |& tee ../logger/UL2016preVFP_SignalMC_Merged_check.txt
# ./MergedAnalysis --config ../config_check/UL2016postVFP_SignalMC_Merged.yaml  |& tee ../logger/UL2016postVFP_SignalMC_Merged_check.txt
# ./MergedAnalysis --config ../config_check/UL2017_SignalMC_Merged.yaml         |& tee ../logger/UL2017_SignalMC_Merged_check.txt
# ./MergedAnalysis --config ../config_check/UL2018_SignalMC_Merged.yaml         |& tee ../logger/UL2018_SignalMC_Merged_check.txt

./ResolvedAnalysis --config ../config_check/UL2016preVFP_SignalMC_Resolved.yaml   |& tee ../logger/UL2016preVFP_SignalMC_Resolved_check_WPLoose.txt
./ResolvedAnalysis --config ../config_check/UL2016postVFP_SignalMC_Resolved.yaml  |& tee ../logger/UL2016postVFP_SignalMC_Resolved_check_WPLoose.txt
./ResolvedAnalysis --config ../config_check/UL2017_SignalMC_Resolved.yaml         |& tee ../logger/UL2017_SignalMC_Resolved_check_WPLoose.txt
./ResolvedAnalysis --config ../config_check/UL2018_SignalMC_Resolved.yaml         |& tee ../logger/UL2018_SignalMC_Resolved_check_WPLoose.txt

# ./MergedAnalysis_HLTstudy --config ../config_check/UL2017_SignalMC_Merged_HLT.yaml         |& tee ../logger/UL2017_SignalMC_Merged_check_HLT.txt
# ./MergedAnalysis_HLTstudy --config ../config_check/UL2018_SignalMC_Merged_HLT.yaml         |& tee ../logger/UL2017_SignalMC_Merged_check_HLT.txt