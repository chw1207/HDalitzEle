#!/bin/sh

# ./MergedAnalysis --config ../config_check/UL2016preVFP_Data_Merged.yaml   |& tee ../logger/UL2016preVFP_Data_Merged_check_WPLoose.txt
# ./MergedAnalysis --config ../config_check/UL2016postVFP_Data_Merged.yaml  |& tee ../logger/UL2016postVFP_Data_Merged_check_WPLoose.txt
# ./MergedAnalysis --config ../config_check/UL2017_Data_Merged.yaml         |& tee ../logger/UL2017_Data_Merged_check_WPLoose.txt
# ./MergedAnalysis --config ../config_check/UL2018_Data_Merged.yaml         |& tee ../logger/UL2018_Data_Merged_check_WPLoose.txt

./ResolvedAnalysis --config ../config_check/UL2016preVFP_Data_Resolved.yaml   |& tee ../logger/UL2016preVFP_Data_Resolved_check_WPLoose.txt
./ResolvedAnalysis --config ../config_check/UL2016postVFP_Data_Resolved.yaml  |& tee ../logger/UL2016postVFP_Data_Resolved_check_WPLoose.txt
./ResolvedAnalysis --config ../config_check/UL2017_Data_Resolved.yaml         |& tee ../logger/UL2017_Data_Resolved_check_WPLoose.txt
./ResolvedAnalysis --config ../config_check/UL2018_Data_Resolved.yaml         |& tee ../logger/UL2018_Data_Resolved_check_WPLoose.txt


# ./MergedAnalysis_HLTstudy --config ../config_check/UL2017_Data_Merged_HLT.yaml         |& tee ../logger/UL2017_Data_Merged_check_WPLoose_HLT.txt