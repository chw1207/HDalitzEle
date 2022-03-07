#!/bin/bash
cd /home/chenghan/CMSSW_10_2_13/src/flashggFinalFit/Background
eval `scramv1 runtime -sh`
$CMSSW_BASE/src/flashggFinalFit/Background/bin/makeBkgPlots -f UntaggedTag_0,UntaggedTag_1,UntaggedTag_2,UntaggedTag_3,UntaggedTag_4,VBFTag_0,VBFTag_1,VBFTag_2,VHHadronicTag,VHTightTag,VHLooseTag -b CMS-HGG_multipdf_auto.root -o outdir_auto/bkgPlots/BkgPlots_cat3.root -d outdir_auto/bkgPlots -c 3 -l "UntaggedTag_3" --sqrts 13  --intLumi 1.000000  --year 2016  --doBands --massStep 1.000 --nllTolerance 0.050 -L 100 -H 180 --higgsResolution 1.0 --isMultiPdf --useBinnedData
