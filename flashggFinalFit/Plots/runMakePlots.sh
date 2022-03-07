cats=(
    "Merged2Gsf_HVBF"   
    "Merged2Gsf_LVBF"   
    "Merged2Gsf_BST"   
    "Merged2Gsf_EBHR9"   
    "Merged2Gsf_EBLR9"   
    "Merged2Gsf_EE"  
    "Merged1Gsf_HVBF"   
    "Merged1Gsf_LVBF"   
    "Merged1Gsf_BST"   
    "Merged1Gsf_EBHR9"   
    "Merged1Gsf_EBLR9"   
    "Merged1Gsf_EE"  
    "Resolved"
)

cats_length=${#cats[@]}
for ((i=0; i<$cats_length; i=i+1))
do
    python makeSplusBModelPlot.py --inputWSFile ../Datacard/electron/datacard_heeg_runII_${cats[i]}_125.root --cats ${cats[i]}
done

# python makeSplusBModelPlot.py --inputWSFile ../Datacard/electron/datacard_heeg_runII_Merged2Gsf_EBHR9_125.root --cats Merged2Gsf_EBHR9 --doBands
# python makeToys.py --inputWSFile ../Datacard/electron/datacard_heeg_runII_Merged2Gsf_EBHR9_125.root --nToys 500 --dryRun
# combine ../Datacard/electron/datacard_heeg_runII_Merged2Gsf_EBHR9_125.root -m 125 -M GenerateOnly --saveWorkspace --expectSignal=0 -t -1 -s 1234 -n _Merged2Gsf_EBHR9_125_gen_step
# combine higgsCombine_Merged2Gsf_EBHR9_125_gen_step.GenerateOnly.mH125.1234.root -m 125 -M MultiDimFit -P r --floatOtherPOIs=1 --saveWorkspace --toysFrequentist --bypassFrequentistFit -t 500 -s 1234 -n _Merged2Gsf_EBHR9_125_fit_step --cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints --X-rtd MINIMIZER_multiMin_maskChannels=2
# combine higgsCombine_Merged2Gsf_EBHR9_125_fit_step.MultiDimFit.mH125.1234.root -m 125 --snapshotName MultiDimFit -M GenerateOnly --saveToys --toysFrequentist --bypassFrequentistFit -t -1 -n _Merged2Gsf_EBHR9_125_throw_step