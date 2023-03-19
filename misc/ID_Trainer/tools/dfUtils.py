import pandas as pd
from fastparquet import ParquetFile
from sys import exit
from typing import Tuple
from sklearn.model_selection import train_test_split
from tools.plotUtils import color
from pprint import pprint

def df_load(infiles, columns=None, cuts=None, extral_text="signal") -> pd.DataFrame:
    print("[INFO] Read files:", flush=True)
    pprint(infiles)

    df = pd.concat([ParquetFile(parquet_file).to_pandas(columns=columns) for parquet_file in infiles])
    if cuts is not None:
        print("[INFO] Applying the {} events/objects with cuts: {}".format(extral_text, cuts), flush=True)
        df = df.query(cuts)

    return df


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




