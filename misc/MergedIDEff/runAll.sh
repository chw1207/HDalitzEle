#!/bin/bash

echo ""
echo "Process ZGToLLG ....."
echo ""

python3 ./python/runAna.py ZGToLLG UL2016preVFP  |& tee ./logger/ZGToLLG_UL2016preVFP.txt
python3 ./python/runAna.py ZGToLLG UL2016postVFP |& tee ./logger/ZGToLLG_UL2016postVFP.txt
# python3 ./python/runAna.py ZGToLLG UL2017        |& tee ./logger/ZGToLLG_UL2017.txt
# python3 ./python/runAna.py ZGToLLG UL2018        |& tee ./logger/ZGToLLG_UL2018.txt

echo ""
echo "Process Data1Mu ....."
echo ""

python3 ./python/runAna.py Data1Mu UL2016preVFP  |& tee ./logger/Data1Mu_UL2016preVFP.txt
python3 ./python/runAna.py Data1Mu UL2016postVFP |& tee ./logger/Data1Mu_UL2016postVFP.txt
# python3 ./python/runAna.py Data1Mu UL2017        |& tee ./logger/Data1Mu_UL2017.txt
# python3 ./python/runAna.py Data1Mu UL2018        |& tee ./logger/Data1Mu_UL2018.txt

echo ""
echo "Process Data2Mu ....."
echo ""

python3 ./python/runAna.py Data2Mu UL2016preVFP  |& tee ./logger/Data2Mu_UL2016preVFP.txt
python3 ./python/runAna.py Data2Mu UL2016postVFP |& tee ./logger/Data2Mu_UL2016postVFP.txt
# python3 ./python/runAna.py Data2Mu UL2017        |& tee ./logger/Data2Mu_UL2017.txt
# python3 ./python/runAna.py Data2Mu UL2018        |& tee ./logger/Data2Mu_UL2018.txt

echo "Process DYJets ....."
echo ""

# python3 ./python/runAna.py DYJets UL2016preVFP   |& tee ./logger/DYJets_UL2016preVFP.txt
# python3 ./python/runAna.py DYJets UL2016postVFP  |& tee ./logger/DYJets_UL2016postVFP.txt
# python3 ./python/runAna.py DYJets UL2017         |& tee ./logger/DYJets_UL2017.txt
# python3 ./python/runAna.py DYJets UL2018         |& tee ./logger/DYJets_UL2018.txt

# echo "Process TTJets ....."
# echo ""

# python3 ./python/runAna.py TTJets UL2016preVFP   |& tee ./logger/TTJets_UL2016preVFP.txt
# python3 ./python/runAna.py TTJets UL2016postVFP  |& tee ./logger/TTJets_UL2016postVFP.txt
# python3 ./python/runAna.py TTJets UL2017         |& tee ./logger/TTJets_UL2017.txt
# python3 ./python/runAna.py TTJets UL2018         |& tee ./logger/TTJets_UL2018.txt

