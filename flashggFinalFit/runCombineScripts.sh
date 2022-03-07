#!/bin/sh

Help() {
    # Display Help
    echo "Syntax: script Template [-t|s|b|m|c|a|h]"
    echo "options:" 
    echo -e "-t --tree2ws     Convert tree to ws" 
    echo -e "-s --signal      Create signal model"
    echo -e "-b --background  Create background model"
    echo -e "-m --makecards   Make datacards"
    echo -e "-c --combination Run Higgs combination plots"
    echo -e "-a --all         Run the full scripts(tree2ws->signal->bkg->datacards->combine)"
    echo -e "-h --help        Print this Help"
    echo
}

start=$(date +%s.%N)

OPTS=$(getopt -o tsbmcah --long tree2ws,signal,background,makecards,combination,all,help -n 'parse-options' -- "$@")


eval set -- "$OPTS"
echo "$OPTS"

if [ "$OPTS" = " --" ]; then
    echo -e "No argument has been passed to the script"
    Help
fi

TREE2WS=false
SIGNAL=false
BACKGROUND=false
MAKECARDS=false
COMBINATION=false
HELP=false

while true; do
    case "$1" in
    -t | --tree2ws)
        TREE2WS=true
        shift
        ;;
    -s | --signal)
        SIGNAL=true
        shift
        ;;
    -b | --background)
        BACKGROUND=true
        shift
        ;;
    -m | --makecards)
        MAKECARDS=true
        shift
        ;;
    -c | --combination)
        COMBINATION=true
        shift
        ;;
    -a | --all)
        TREE2WS=true
        SIGNAL=true
        BACKGROUND=true
        MAKECARDS=true
        COMBINATION=true
        shift
        ;;
    -h | --help)
        HELP=true
        shift
        shift
        ;;
    --)
        shift
        break
        ;;
    *) break ;;
    esac
done

if [ "$HELP" = "true" ]; then
    Help
fi

if [ "$TREE2WS" = "true" ]; then
    cd ./Trees2WS/
    echo -e "*********************************************************************************** \n[INFO] Start converting tree 2 ws\n change directory into ./Trees2WS/\n run the script (runTree2WS.sh)\n"
    sh runTree2WS.sh
    cd ../
    echo -e "[INFO] Generation of WS has done.\n Now in ${PWD}\n***********************************************************************************\n"
fi

if [ "$SIGNAL" = "true" ]; then

    echo -e "*********************************************************************************** \n[INFO] Start creating siganl model\n change directory into ./Signal/\n run the script (runSignalScripts.py)\n"
    cd ./Signal/
    python runSignalScripts.py --doFitting 1 --doInterpolation 1 --makeFinalPlot 1
    cd ../
    echo -e "[INFO] Generation of signal model has done.\n Now in ${PWD}\n***********************************************************************************\n"
fi

if [ "$BACKGROUND" = "true" ]; then
    echo -e "\n*********************************************************************************** \n[INFO] Start creating background model\n change directory into ./Background/\n run the script (runBackgroundScripts_HDalitz.sh)\n"
    cd ./Background/
    sh runBackgroundScripts_HDalitz.sh 
    cd ../
    echo -e "[INFO] Generation of Background model has done.\n Now in ${PWD}\n***********************************************************************************\n"
fi

if [ "$MAKECARDS" = "true" ]; then
    echo -e "\n*********************************************************************************** \n[INFO] Start creating datacards\n change directory into ./Datacard/\n run the script (run_makeYields.sh and makeDatacards_HDalitz.py) with option:\n"
    cd ./Datacard/
    sh run_makeYields.sh
    python makeDatacard_HDalitz.py
    cd ../
    echo -e "[INFO] Generation of datacards has done.\n Now in ${PWD}\n***********************************************************************************\n"
fi

if [ "$COMBINATION" = "true" ]; then
    echo -e "\n*********************************************************************************** \n[INFO] Start running Combination tool\n change directory into ./Datacard/electron \n"
    cd ./Datacard/electron # FIXEDME
    sh doCombination_HDalitz.sh
    cd ../../Combine
    python LimitPlotter.py
    cd ../
    echo -e "[INFO] Combination has done.\n Now in ${PWD}\n***********************************************************************************\n"
fi

end=$(date +%s.%N)
echo "Runtime = $(echo "$end - $start" | bc) secods"



