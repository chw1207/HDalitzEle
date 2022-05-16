import pandas as pd
import pickle as pkl
from sys import exit
from typing import Tuple
from sklearn.model_selection import train_test_split
from Tools.plotUtils import color


def df_load(infile, cuts=None, extral_text="signal") -> pd.DataFrame:
    print("[INFO] Loading {} dataframes from: {}".format(extral_text, infile), flush=True)
    try:
        with open(infile, "rb") as f:
            if cuts is not None:
                print("[INFO] Applying the {} events/objects with cuts: {}".format(extral_text, cuts), flush=True)
                df = pkl.load(f).query(cuts)
            else:
                df = pkl.load(f)
        return df
    except FileNotFoundError:
        print("[ERROR] File does not exist :(", flush=True)
        print(">>> {}".format(infile), flush=True)
        exit(-1)


def df_split(Mdf, test_size=0.2, seed=42, isb=False, trueLable="EleType", Classes=[""]) -> Tuple[list, list, pd.DataFrame]:
    Mdf[trueLable] = 0
    if isb == True:
        Mdf.loc[Mdf.Class == "Signal", trueLable] = 1
        Mdf.loc[Mdf.Class == "Background", trueLable] = 0
    else:
        for i, k in enumerate(Classes):
            Mdf.loc[Mdf.Class == k, trueLable] = i

    index = Mdf.index
    TrainIndices, TestIndices = [], []
    for myclass in Classes:
        Indices = index[Mdf["Class"] == myclass].values.tolist()
        myclassTrainIndices, myclassTestIndices = train_test_split(Indices, test_size=test_size, random_state=seed, shuffle=True)
        TrainIndices = TrainIndices + myclassTrainIndices
        TestIndices = TestIndices + myclassTestIndices

    Mdf.loc[TrainIndices, "Dataset"] = "Train"
    Mdf.loc[TestIndices, "Dataset"] = "Test"

    Mdf.loc[TrainIndices, "TrainDataset"] = 1
    Mdf.loc[TestIndices, "TrainDataset"] = 0

    print(color.GREEN + "Reading classes:" + color.END, flush=True)
    print(Mdf.Class.unique().tolist(), flush=True)

    return TrainIndices, TestIndices, Mdf




