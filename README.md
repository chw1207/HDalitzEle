# Higgs Dalitz electron channel analysis

This repository has been tailored for conducting CMS Run2 H→γ*γ→eeγ analysis. The data input required for this framework is formatted in [ggNtuple](https://github.com/cmkuo/ggAnalysis/tree/106X). To process the data efficiently, the framework utilizes [RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html), enabling parallel data processing.

## Environment setup
This framework relies on two distinct environments. The `hdalitz2` environment comprises numerous Python packages such as **xgboost**, **optuna**, and **dask**. This environment is primarily utilized for machine learning purposes, specifically for ID and energy regression training within the programs found in the `misc` directory. 

On the other hand, the `hdalitzAna` environment is dedicated to event selections and is equipped with a smaller set of C++ packages like **ROOT** and **boost-cpp**.

```bash
git clone https://github.com/chw1207/HDalitzEle.git
cd HDalitzEle
conda create --name hdalitz2 --file requirements.txt
conda create --name hdalitzAna --file requirements_c.txt
```

After the installations of necessary dependencies, you can set up the environments by 
```bash
source env.sh
source setup.sh # build the framworks for the programs in misc directory (only need once)
```
and
```bash
source env_ana.sh
source setup_ana.sh # build the framworks for the main analysis (only need once)
```

Once the environments are established successfully, some executables - **MergedAnalysis**, **ResolvedAnalysis**, **drawKinematics** and **runSignificance** will become available within the  `bin` directory.

## Usage
### 1. Producing mini trees following event selections
The computations are carried out by two key executables.
- **MergedAnalysis** focuses on selecting a merged electron and a high-momentum photon.
- **ResolvedAnalysis** is responsible for selecting two separate electrons and a photon.

The **MergedAnalysis** and **ResolvedAnalysis** can be invoked with
```bash
./MergedAnalysis --config ../config_check/configuration_files files --range number [--skipSS]
# eg. ./MergedAnalysis --config ../config_check/UL2017_SignalMC_test_Merged.yaml --skipSS

./ResolvedAnalysis  --config ../config_check/configuration_files files --range number 
# eg. ./ResolvedAnalysis --config ../config_check/UL2017_SignalMC_test_Resolved.yaml
```
- `-c --config`: To specify the configuration file which contains global settings such as the N-tuple paths and paths to store the mini trees.
- `-r --range`: To decide the maximal number of events to process. If it is not specified, process all the events by default.
- `--skipSS`: If specified, program will not load the trees containing shower shape corrected variables and uses the default values for shower shape variables.

Submitting these tasks to a computing cluster with multiple nodes can be efficiently achieved using **HTCondor**. For instance, they can be seamlessly submitted to the ChiP server at NCU. 
```bash
cd bin
condor_submit submit_mc.sub 
```
----
### 2. Evaluating the performance
The **runSignificance** program computes the approximate median significance for each category. Here, the signal is estimated within the 2σ region, while the non-resonant background is estimated from a data side-band region, scaled to match the mass window of the signal.
```bash
./runSignificance --config ../config_check/ULfullRun2_plots_Merged.yaml
```

