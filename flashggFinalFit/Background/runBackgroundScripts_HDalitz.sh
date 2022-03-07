#!/bin/sh

make 
# echo " 2016 bkg...."
# nohup ./bin/fTest_HDalitz --infilename "~/CMSSW_10_2_13/src/flashggFinalFit/Trees2WS/WS/2016/data_obs_2016.root" --plotDir "./plots/fTest/2016" --multipdfPath "./MultiPdf/2016/Multipdf_2016.root" --year 2016 --luminosity 36.33 &> bkg_2016.txt &
# sleep 10

# echo " 2017 bkg...."
# nohup ./bin/fTest_HDalitz --infilename "~/CMSSW_10_2_13/src/flashggFinalFit/Trees2WS/WS/2017/data_obs_2017.root" --plotDir "./plots/fTest/2017" --multipdfPath "./MultiPdf/2017/Multipdf_2017.root" --year 2017 --luminosity 41.48 &> bkg_2017.txt &
# sleep 10

# echo " 2018 bkg...."
# nohup ./bin/fTest_HDalitz --infilename "~/CMSSW_10_2_13/src/flashggFinalFit/Trees2WS/WS/2018/data_obs_2018.root" --plotDir "./plots/fTest/2018" --multipdfPath "./MultiPdf/2018/Multipdf_2018.root" --year 2018 --luminosity 59.35 &> bkg_2018.txt &

./bin/fTest_HDalitz --infilename "../Trees2WS/WS/data_obs.root" --plotDir "./plots/fTest" --multipdfPath "./MultiPdf/Multipdf.root" --luminosity 137.1 --runFtestCheckWithToys 1


# nohup ./bin/fTest_HDalitz --infilename "~/CMSSW_10_2_13/src/flashggFinalFit/Trees2WS/WS/data_obs.root" --plotDir "./plots/fTest" --multipdfPath "./MultiPdf/Multipdf.root" --luminosity 137.1 --runFtestCheckWithToys 1 &> bkg_RunII.txt &