#!/bin/sh

{
    echo -e "Merged-2Gsf ID training"
    python3 Trainer-HDalitz-hyper.py config/TrainConfig_Merged2GsfID_EB_opt
    python3 Trainer-HDalitz-hyper.py config/TrainConfig_Merged2GsfID_EE_opt

    echo -e "Merged-1Gsf ID training"
    python3 Trainer-HDalitz-hyper.py config/TrainConfig_Merged1GsfID_EB_opt
    python3 Trainer-HDalitz-hyper.py config/TrainConfig_Merged1GsfID_EE_opt
} || {
    set -e
}