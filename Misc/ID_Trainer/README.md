# HDalitz-trainer

### Original package: https://github.com/cms-egamma/ID-Trainer 

This trainer is used to train the dedicated merged electron ID in the <img src="https://render.githubusercontent.com/render/math?math=H\rightarrow\gamma^*\gamma\rightarrow ee\gamma"> analysis. It runs the xgboost on the multi-GPU node to increase the trainning speed and the avalible memory. (NCUHEP chip04 local server with GPU clusters)

- **Merged IDs:** Based on **the number of Gsf tracks(nGsfTracks)** matched to the electron, **Merged-2Gsf** and **Merged-1Gsf IDs** are developed to discriminate the merged electrons and the associated background.
  - **Merged-2Gsf ID**(nGsfTracks > 1)
  - **Merged-1Gsf ID**(nGsfTracks = 1)
  - Electrons will be didved into the **Merged-1/2Gsf**(merged <img src="https://render.githubusercontent.com/render/math?math=e">), **DYJets**(real single <img src="https://render.githubusercontent.com/render/math?math=e">) and **QCD**(fake <img src="https://render.githubusercontent.com/render/math?math=e">).

---
### Prepare training samples 
- PrepareData_MergedID.py

One could exucute the python scripts via
```bash
$ python3 PrepareData_MergedID.py
```
---

### Train the model
- Trainer-HDalitz-dask.py

One could exucute the python scripts via
```bash
$ python3 Trainer-HDalitz-dask.py configs/[Training config python file]
# eq. 
# python3 Trainer-HDalitz-dask.py configs/TrainConfig_Merged2GsfID_EB
```

- There are currently 4 training configuration files
    - TrainConfig_Merged2GsfID_EB.py
    - TrainConfig_Merged2GsfID_EE.py
    - TrainConfig_Merged1GsfID_EB.py
    - TrainConfig_Merged1GsfID_EE.py
---

### plot the results
The setting of the plots can be modified in the ```Tools/PlotTools.py```. The plotting style of the distributions of the features can be modified in the json files in Tools directory.
