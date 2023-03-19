import os, sys
import ROOT
import time

import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")

import pandas as pd
import dask.dataframe as dd
from glob import glob
from tqdm import tqdm
from tools.plotUtils import *
from pprint import pprint


# move data from RDataFrame to pandas DataFrame
def rdf2df(rdf, filter=None):
    branches_arr = [] # remove useless varaibles for training in tree
    for var in rdf.GetColumnNames():
        if ("mc" in str(var) and str(var) != "mcwei"): # generator variables
            continue
        if (str(var) == "isPVGood"): # bool type (not allowed by parquet file)
            continue
        if (rdf.GetColumnType(var) == "ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double32_t> >"): # LorentzVector variables
            continue
        if ("fCoordinates" in str(var)): # Elements of LorentzVector variables
            continue
        branches_arr.append(str(var))
    if filter == None:
        arr = rdf.AsNumpy(columns=branches_arr)
        df  = pd.DataFrame(arr)
    else:
        arr = rdf.Filter(filter).AsNumpy(columns=branches_arr)
        df  = pd.DataFrame(arr)

    return df


# split DataFrame into EB, EE regions and save to parquet files
def ddf2parquet(ddf, outpath):
    sel = ["(abs(eleSCEta) <= 1.479)", "(abs(eleSCEta) > 1.479 and abs(eleSCEta) <= 2.5)"]
    for j, reg in enumerate(["EB", "EE"]):
        ddf_reg = ddf.query(sel[j])

        out_path_reg = "{}-{}".format(outpath, reg)
        print("[INFO] Save the parquet files in {}".format(out_path_reg), flush=True)
        ddf_reg.to_parquet(
            "{}-{}".format(outpath, reg),
            engine="fastparquet",
            overwrite=True,
            write_metadata_file=False
        )


def main():
    ext = "FullRun2UL"
    samples = {
        "HDalitz": glob("../GenStudy/miniTree/*/miniTree_HDalitz_*_eeg_*.root"),
        "DYJets" : glob("../GenStudy/miniTree/*/miniTree_DYJets_*.root"),
        "QCD"    : glob("../GenStudy/miniTree/*/miniTree_QCD_*.root")
    }
    weights = ["rho", "nVtx", "nGoodVtx", "instwei", "procXS", "event"]

    ROOT.EnableImplicitMT(25)
    for proc, files in samples.items():
        print("[INFO] Find_files(): {} files are found for {}".format(len(files), proc), flush=True)
        if proc == "HDalitz": # H -> g*g -> e'g (e' is merged electron)
            rdf = ROOT.RDataFrame("miniTree", files)
            for cat in ["Merged-2Gsf", "Merged-1Gsf"]:
                print("[INFO] Preparing DataFrame for {}".format(cat), flush=True)

                # event % 2 == 0 for training, == 1 for analysis
                cat_cut = 2 if cat == "Merged-2Gsf" else 3
                filter = "category == {}".format(cat_cut)
                ddf = dd.from_pandas(rdf2df(rdf, filter=filter), npartitions=20)

                # extract the dataframe for leading electrons
                ddf_Lead   = ddf[[branch for branch in ddf.columns if ("_Lead" in branch) or (branch in weights)]]

                # rename the columns (remove Lead)
                ddf_Lead.columns = ddf_Lead.columns.str.replace("_Lead", "")

                # combine the dataframes
                # specify interleave_partitions=True to ignore order
                ddf_new = ddf_Lead.dropna()
                ddf_new["Class"] = cat # for multi-classification
                ddf2parquet(ddf_new, "./data/DataFrames-{}-{}".format(cat, ext))

        if proc == "DYJets": # Z -> ee (leading and trainling electrons)
            for i in range(len(files)): # memory issue, do not load all of the files at once
                print("[INFO] Preparing DataFrame for file-{}".format(i), flush=True)
                rdf = ROOT.RDataFrame("miniTree", files[i])
                ddf = dd.from_pandas(rdf2df(rdf, filter="category == 1"), npartitions=20)

                # extract the dataframes for leading and triling electrons
                ddf_Lead    = ddf[[branch for branch in ddf.columns if ("_Lead" in branch) or (branch in weights)]]
                ddf_subLead = ddf[[branch for branch in ddf.columns if ("_subLead" in branch) or (branch in weights)]]

                # rename the columns (remove Lead and subLead)
                ddf_Lead.columns    = ddf_Lead.columns.str.replace("_Lead", "")
                ddf_subLead.columns = ddf_subLead.columns.str.replace("_subLead", "")

                # combine the dataframes
                # specify interleave_partitions=True to ignore order
                ddf_new = dd.concat([ddf_Lead, ddf_subLead], interleave_partitions=True).dropna()
                ddf_new["Class"] = "DYJets" # for multi-classification
                ddf2parquet(ddf_new, "./data/DataFrames-DYJets-{}-{}".format(ext, i))

        if proc == "QCD":
            print("[INFO] Preparing DataFrame for QCD", flush=True)
            rdf = ROOT.RDataFrame("miniTree", files)
            ddf = dd.from_pandas(rdf2df(rdf), npartitions=20)

            # extract the dataframe for leading electrons
            ddf_Lead   = ddf[[branch for branch in ddf.columns if ("_Lead" in branch) or (branch in weights)]]

            # # rename the columns (remove Lead) and flatten the dataframe
            ddf_Lead.columns = ddf_Lead.columns.str.replace("_Lead", "")

            # specify interleave_partitions=True to ignore order
            ddf_new = ddf_Lead.dropna()
            ddf_new["Class"] = "QCD" # for multi-classification
            ddf2parquet(ddf_new, "./data/DataFrames-QCD-{}".format(ext))

        print("", flush=True)


if __name__ == "__main__":
    start_time = time.time()

    print(color.BLUE + "---Start to prepare the dataframe for training!---" + color.END, flush=True)
    os.makedirs("./data", exist_ok=True)
    main()
    print(color.BLUE + "---All done!---" + color.END, flush=True)

    seconds=time.time() - start_time
    print("[INFO] Time Taken:", time.strftime("%H:%M:%S", time.gmtime(seconds)), flush=True)