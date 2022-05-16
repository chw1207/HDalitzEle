# ID-Trainer

### Original package: https://github.com/cms-egamma/ID-Trainer
This trainer is used to train the dedicated **merged electron ID**. Several features are added compared to the original package.
-  **Multi-GPUs training**: To increase the trainning speed and the avalible memory.
-  **Hyperparameter  optimization**: The hyperparameters can be tuned automatically based on a algorithm called [TPE](https://journals.sagepub.com/doi/pdf/10.1177/0020294020932347](https://proceedings.neurips.cc/paper/2011/file/86e8f7ab32cfd12577bc2619bc635690-Paper.pdf)) by the python package [hyperopt](https://github.com/hyperopt/hyperopt).

Based on the number of **Gsf tracks(nGsfTracks)** matched to the electron, **Merged-2Gsf** and **Merged-1Gsf IDs** are developed to discriminate the merged electrons and the associated background.
-  **Merged-2Gsf ID** (nGsfTracks > 1)
-  **Merged-1Gsf ID** (nGsfTracks = 1)
-  Associated background: real single electrons and jet fake electrons.
---

### Prepare training samples
Merged electrons are selected as signal from the <img src="https://render.githubusercontent.com/render/math?math=\gamma^*\rightarrow ee"> process. Real single electrons are provided by DYJetsToLL sample and jet fake electrons are provided by QCD sample. One could exucute the python scripts via
```bash
$ python3 PrepareData_MergedID.py
```
---

### Train the ID
Before performing the training, one could modify the training setting in the following training configuration files
- config/TrainConfig_Merged2GsfID_EB_opt.py
- config/TrainConfig_Merged2GsfID_EE_opt.py
- config/TrainConfig_Merged1GsfID_EB_opt.py
- config/TrainConfig_Merged1GsfID_EE_opt.py

One could further exucute the python scripts via
```bash
python Trainer-HDalitz-hyper.py config/[Training config python file]
# eq.
# python3 Trainer-HDalitz-hyper.py configs/TrainConfig_Merged2GsfID_EB_opt
```

---