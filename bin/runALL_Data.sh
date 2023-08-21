#!/bin/sh

# ./MergedAnalysis --config ../config/UL2016preVFP_Data_Merged.yaml   |& tee ../logger/UL2016preVFP_Data_Merged.txt
# ./MergedAnalysis --config ../config/UL2016postVFP_Data_Merged.yaml  |& tee ../logger/UL2016postVFP_Data_Merged.txt
# ./MergedAnalysis --config ../config/UL2017_Data_Merged.yaml         |& tee ../logger/UL2017_Data_Merged.txt
# ./MergedAnalysis --config ../config/UL2018_Data_Merged.yaml         |& tee ../logger/UL2018_Data_Merged.txt

./ResolvedAnalysis --config ../config/UL2016preVFP_Data_Resolved.yaml   |& tee ../logger/UL2016preVFP_Data_Resolved.txt
./ResolvedAnalysis --config ../config/UL2016postVFP_Data_Resolved.yaml  |& tee ../logger/UL2016postVFP_Data_Resolved.txt
./ResolvedAnalysis --config ../config/UL2017_Data_Resolved.yaml         |& tee ../logger/UL2017_Data_Resolved.txt
./ResolvedAnalysis --config ../config/UL2018_Data_Resolved.yaml         |& tee ../logger/UL2018_Data_Resolved.txt