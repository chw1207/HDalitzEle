#!/bin/sh

./runAnalysis --config ../config/UL2016preVFP_SignalMC.yaml
./runAnalysis --config ../config/UL2016postVFP_SignalMC.yaml
./runAnalysis --config ../config/UL2017_SignalMC.yaml
./runAnalysis --config ../config/UL2018_SignalMC.yaml

./runAnalysis --config ../config/UL2016preVFP_Data.yaml
./runAnalysis --config ../config/UL2016postVFP_Data.yaml
./runAnalysis --config ../config/UL2017_Data.yaml
./runAnalysis --config ../config/UL2018_Data.yaml