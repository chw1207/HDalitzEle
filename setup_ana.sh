#!/bin/sh

echo -e ""
echo -e "******************************************************"
echo -e "*       Search for Higgs Dalitz decay at CMS         *"
echo -e "*                  H -> g*g -> eeg                   *"
echo -e "*   Developed by Cheng-Han Wu, cheng-han.wu@cern.ch  *"
echo -e "******************************************************"
echo -e ""

wdir="${PWD}/bin"
dir="./build_ana"
if [ -d "$dir" ]; then
    rm -r $dir
fi
mkdir $dir
cd $dir

echo -e "1. Build the project for analysis: " $dir "\n"

echo -e "2. Construct the analysis project: HDalitzEle\n"
cmake ..

echo -e ""
echo -e "3. Compile the shared library in the project\n"
make -j10

echo -e ""
echo -e "4. Working directory: " $wdir "\n"