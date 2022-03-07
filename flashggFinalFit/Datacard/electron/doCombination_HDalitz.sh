#!/bin/sh


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


# for test
# cardsName="datacard_heeg_runII_Merged2Gsf_EBHR9_125.txt"
# combine -M AsymptoticLimits $cardsName -m 125 -n _Merged2Gsf_EBHR9_125 --run=blind --cminDefaultMinimizerStrategy 0


cats_length=${#cats[@]}

masses=(120 121 122 123 124 125 126 127 128 129 130)
# masses=(125)
mass_length=${#masses[@]}

for ((i=0; i<$mass_length; i=i+1))
do
    echo "[INFO] Combine the cards of each category: mass = ${masses[i]} GeV"
    combineCards.py datacard_heeg_runII_Merged2Gsf_HVBF_${masses[i]}.txt datacard_heeg_runII_Merged2Gsf_LVBF_${masses[i]}.txt datacard_heeg_runII_Merged2Gsf_BST_${masses[i]}.txt datacard_heeg_runII_Merged2Gsf_EBHR9_${masses[i]}.txt datacard_heeg_runII_Merged2Gsf_EBLR9_${masses[i]}.txt datacard_heeg_runII_Merged2Gsf_EE_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_HVBF_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_LVBF_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_BST_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_EBHR9_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_EBLR9_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_EE_${masses[i]}.txt datacard_heeg_runII_Resolved_${masses[i]}.txt > datacard_heeg_runII_combine2_${masses[i]}.txt
    echo "[INFO] Compute the combine limit: mass = ${masses[i]} GeV"
    combine -M AsymptoticLimits datacard_heeg_runII_combine2_${masses[i]}.txt -m ${masses[i]} -n _combine_${masses[i]} --run=blind --cminDefaultMinimizerStrategy 0 
done

for ((i=0; i<$mass_length; i=i+1))
do
    echo "[INFO] Compute the combine significance: mass = ${masses[i]} GeV"
    combine -M Significance datacard_heeg_runII_combine2_${masses[i]}.txt -m ${masses[i]} -n _combine_${masses[i]}_expSignificance -t -1 --expectSignal=1 
done

masses=(125)
mass_length=${#masses[@]}
for ((i=0; i<$mass_length; i=i+1))
do
    echo "[INFO] Combine the cards of Merged2Gsf tagged category: mass = ${masses[i]} GeV"
    combineCards.py datacard_heeg_runII_Merged2Gsf_HVBF_${masses[i]}.txt datacard_heeg_runII_Merged2Gsf_LVBF_${masses[i]}.txt datacard_heeg_runII_Merged2Gsf_BST_${masses[i]}.txt  > datacard_heeg_runII_Merged2Gsf_tagged_${masses[i]}.txt
    echo "[INFO] Compute the combine limit: mass = ${masses[i]} GeV"
    combine -M AsymptoticLimits datacard_heeg_runII_Merged2Gsf_tagged_${masses[i]}.txt -m ${masses[i]} -n _Merged2Gsf_tagged_${masses[i]} --run=blind --cminDefaultMinimizerStrategy 0 

    echo "[INFO] Combine the cards of Merged2Gsf Untagged category: mass = ${masses[i]} GeV"
    combineCards.py datacard_heeg_runII_Merged2Gsf_EBHR9_${masses[i]}.txt datacard_heeg_runII_Merged2Gsf_EBLR9_${masses[i]}.txt datacard_heeg_runII_Merged2Gsf_EE_${masses[i]}.txt > datacard_heeg_runII_Merged2Gsf_Untagged_${masses[i]}.txt
    echo "[INFO] Compute the combine limit: mass = ${masses[i]} GeV"
    combine -M AsymptoticLimits datacard_heeg_runII_Merged2Gsf_Untagged_${masses[i]}.txt -m ${masses[i]} -n _Merged2Gsf_Untagged_${masses[i]} --run=blind --cminDefaultMinimizerStrategy 0

    echo "[INFO] Combine the cards of Merged1Gsf tagged category: mass = ${masses[i]} GeV"
    combineCards.py datacard_heeg_runII_Merged1Gsf_HVBF_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_LVBF_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_BST_${masses[i]}.txt  > datacard_heeg_runII_Merged1Gsf_tagged_${masses[i]}.txt
    echo "[INFO] Compute the combine limit: mass = ${masses[i]} GeV"
    combine -M AsymptoticLimits datacard_heeg_runII_Merged1Gsf_tagged_${masses[i]}.txt -m ${masses[i]} -n _Merged1Gsf_tagged_${masses[i]} --run=blind --cminDefaultMinimizerStrategy 0 

    echo "[INFO] Combine the cards of Merged1Gsf Untagged category: mass = ${masses[i]} GeV"
    combineCards.py datacard_heeg_runII_Merged1Gsf_EBHR9_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_EBLR9_${masses[i]}.txt datacard_heeg_runII_Merged1Gsf_EE_${masses[i]}.txt > datacard_heeg_runII_Merged1Gsf_Untagged_${masses[i]}.txt
    echo "[INFO] Compute the combine limit: mass = ${masses[i]} GeV"
    combine -M AsymptoticLimits datacard_heeg_runII_Merged1Gsf_Untagged_${masses[i]}.txt -m ${masses[i]} -n _Merged1Gsf_Untagged_${masses[i]} --run=blind --cminDefaultMinimizerStrategy 0
done

combine -M AsymptoticLimits datacard_heeg_runII_Resolved_125.txt -m 125 -n _Resolved_125 --run=blind --cminDefaultMinimizerStrategy 0





# use to make the S+B plots
for ((i=0; i<$cats_length; i=i+1))
do
    text2workspace.py datacard_heeg_runII_${cats[i]}_125.txt -o datacard_heeg_runII_${cats[i]}_125.root -m 125
done

text2workspace.py datacard_heeg_runII_combine_125.txt -o datacard_heeg_runII_combine_125.root -m 125