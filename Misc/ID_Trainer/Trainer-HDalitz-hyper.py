import sys, os
import time
import json as js
import pickle as pkl
from datetime import datetime
from importlib import import_module

# https://stackoverflow.com/a/15778297
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import pandas as pd
import numpy as np

import xgboost as xgb
import dask
from dask import array as da
from dask.distributed import Client
from dask_cuda import LocalCUDACluster
import hyperopt as hpt
from hyperopt import hp
from hyperopt.early_stop import no_progress_loss
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report
from sklearn.utils import shuffle

from Tools.plotUtils import *
from Tools.dfUtils import df_load, df_split
from Tools.hyperUtils import trial_count, NpEncoder
from Tools.ptEtaRwt import df_pteta_rwt, df_balance_rwt


# https://github.com/hyperopt/hyperopt/wiki/FMin#21-parameter-expressions
space = {
    "grow_policy": hp.choice("grow_policy", ["depthwise", "lossguide"]),
    "learning_rate": hp.loguniform("learning_rate", np.log(0.01), np.log(0.5)),
    "max_depth": hp.choice("max_depth", np.arange(5, 20)),
    "max_delta_step": hp.choice("max_delta_step", np.arange(5, 20)),
    "min_child_weight": hp.loguniform("min_child_weight", np.log(10), np.log(1000)),
    "subsample": hp.loguniform("subsample", np.log(0.1), np.log(1)),
    "min_split_loss": hp.uniform("min_split_loss", 0, 2),
    "colsample_bytree": hp.uniform("colsample_bytree", 0.1, 1)
}


@trial_count
def objective(params):
    print("Parameters: {}".format(params), flush=True)
    bst = xgb.dask.train(
        client, params, dtrain,
        num_boost_round=cfg.num_boost_round, early_stopping_rounds=cfg.early_stopping_rounds,
        evals=[(dtrain, "train"), (dtest, "test")], verbose_eval=False
    )
    loss = min(bst["history"]["test"][cfg.model_params["eval_metric"]])
    print("{}: {}".format(cfg.model_params["eval_metric"], loss), flush=True)
    print("", flush=True)

    return loss


def Training(params):
    bst = xgb.dask.train(
        client, params, dtrain,
        num_boost_round=cfg.num_boost_round, early_stopping_rounds=cfg.early_stopping_rounds,
        evals=[(dtrain, "train"), (dtest, "test")]
    )

    print(color.GREEN + "Save the training model" + color.END)
    print("[INFO] Save the model in {}".format("{}/{}_modelXGB.txt".format(mva_dir, MVA)), flush=True)
    bst["booster"].save_model("{}/{}_modelXGB.txt".format(mva_dir, MVA))

    print("[INFO] Save the model in {}".format("{}/{}_modelXGB.pkl".format(mva_dir, MVA)), flush=True)
    pkl.dump(bst["booster"], open("{}/{}_modelXGB.pkl".format(mva_dir, MVA), "wb"))

    return bst


# prepare the dask array for training
def PrepDataset(df, TrainIndices, TestIndices, features, cat, weight):
    # training dataset
    X_train = np.asarray(df.loc[TrainIndices, features])
    Y_train = np.asarray(df.loc[TrainIndices, cat])
    Wt_train = np.asarray(df.loc[TrainIndices, weight])

    # testing dataset
    X_test = np.asarray(df.loc[TestIndices, features])
    Y_test = np.asarray(df.loc[TestIndices, cat])
    Wt_test = np.asarray(df.loc[TestIndices, weight])

    return da.from_array(X_train, chunks=5000), da.from_array(Y_train, chunks=5000), da.from_array(Wt_train, chunks=5000),\
           da.from_array(X_test, chunks=5000), da.from_array(Y_test, chunks=5000), da.from_array(Wt_test, chunks=5000)


def main():
    #====================================================#
    #        load the training data and preprocess       #
    #====================================================#
    print(color.GREEN + "Loading the training dataframe!" + color.END, flush=True)
    sig_df = df_load(cfg.Pickle_signal, cuts=cfg.CommonCut, extral_text="signal")
    bkg_df = df_load(cfg.Pickle_bkg, cuts=cfg.CommonCut, extral_text="background")
    if isBinary == True:
        sig_df["Class"], bkg_df["Class"]  = "Signal", "Background"
    df = shuffle(pd.concat([sig_df, bkg_df], ignore_index=True, sort=False), random_state=cfg.RandomState)

    # change the dtype
    cols = df.select_dtypes(include=[object]).columns.drop("Class")
    df[cols] = df[cols].astype(np.float32)
    if cfg.uniformwt == True:
        df["instwei"] = 1. #! No relative xs weight

    # split df to taining and testing sample
    TrainIndices, TestIndices, df = df_split(
        df, test_size=cfg.testsize, seed=cfg.RandomState,
        isb=isBinary, trueLable=cat, Classes=cfg.Classes
    )

    #====================================================#
    #       Draw the training feature distributions      #
    #====================================================#
    print(color.GREEN + "Draw the training feature plots" + color.END, flush=True)
    with open("./Tools/{}".format(cfg.featureplotparam_json), "r") as fp:
        myTags = js.load(fp)
        for tagName, par in myTags.items():
            setlogy = True if par["logy"] == "true" else False
            if tagName not in cfg.features:
                continue
            plot_featureHist(
                df, tagName, len(cfg.Classes), cfg.Classes,
                par["bininfo"], par["xAxisLabel"],
                histcolor=cfg.ClassColors,
                logy=setlogy, y_axisscale=par["y_axisscale"],
                outname="{}/FeaturePlot/TrainFeature_{}".format(cfg.OutputDirName, tagName),
                wei="instwei",
                category=cat
            )

    #====================================================#
    #                 add training weight                #
    #====================================================#
    if cfg.Reweighing == "Balanced":
        print(color.GREEN + "Balanced reweighting for training sample" + color.END, flush=True)
        df.loc[TrainIndices, weight] = df_balance_rwt(df.loc[TrainIndices], SumWeightCol="instwei", NewWeightCol=weight, Classes=cfg.Classes)

        print(color.GREEN + "Balanced reweighting for testing sample" + color.END, flush=True)
        df.loc[TestIndices, weight] = df_balance_rwt(df.loc[TestIndices], SumWeightCol="instwei", NewWeightCol=weight, Classes=cfg.Classes)
    elif cfg.Reweighing in cfg.Classes:
        print(color.GREEN + "pt-eta reweighting to {} for training sample".format(cfg.Reweighing) + color.END, flush=True)
        df.loc[TrainIndices, weight] = df_pteta_rwt(
            df.loc[TrainIndices], "Class",
            ptw=cfg.ptbins, etaw=cfg.etabins, pt=cfg.ptwtvar, eta=cfg.etawtvar,
            SumWeightCol="instwei", NewWeightCol=weight, cand=cfg.Reweighing, Classes=cfg.Classes
        )
        plot_ptEtaRwt(df, cfg.ptwtvar, cfg.etawtvar, cfg.ptbins, cfg.etabins, mva_dir)

        print(color.GREEN + "pt-eta reweighting to {} for testing sample".format(cfg.Reweighing) + color.END, flush=True)
        df.loc[TestIndices, weight] = df_pteta_rwt(
            df.loc[TestIndices], "Class",
            ptw=cfg.ptbins, etaw=cfg.etabins, pt=cfg.ptwtvar, eta=cfg.etawtvar,
            SumWeightCol="instwei", NewWeightCol=weight, cand=cfg.Reweighing, Classes=cfg.Classes
        )
    else:
        df[weight] = 1.
        print(color.YELLOW + "[WARNING] {} reweighting is not supported!".format(cfg.Reweighing) + color.END, flush=True)
        print(color.YELLOW + "[WARNING] weight is temporarily assigned as 1.!" + color.END, flush=True)


    #====================================================#
    #     change the dtype to dask array for training    #
    #====================================================#
    global X_train, Y_train, Wt_train, X_test, Y_test, Wt_test
    X_train, Y_train, Wt_train, X_test, Y_test, Wt_test = PrepDataset(df, TrainIndices, TestIndices, cfg.features, cat, weight)

    global dtrain, dtest
    dtrain = xgb.dask.DaskDMatrix(client, X_train, Y_train, weight= Wt_train)
    dtest = xgb.dask.DaskDMatrix(client, X_test, Y_test, weight=Wt_test)

    #====================================================#
    #              hyperparameter searching              #
    #====================================================#
    bestPar = {**cfg.model_params, **cfg.hyper_params}
    if cfg.opt == True:
        print(color.GREEN + "Perform the hyper parameter searching" + color.END, flush=True)
        trials = hpt.Trials()
        search_params = {**cfg.model_params, **space}

        best = hpt.fmin(
            objective, search_params, trials=trials,
            algo=hpt.tpe.suggest, max_evals=1000, rstate=np.random.RandomState(cfg.RandomState),
            early_stop_fn=no_progress_loss(iteration_stop_count=20, percent_increase=0.0001),
            show_progressbar=False
        )
        minLoss = min(trials.losses())
        bestPar = hpt.space_eval(search_params, best)
        print(color.GREEN + "Best hyperparameters with {} = {}: ".format(search_params["eval_metric"], minLoss) + color.END, flush=True)
        print(bestPar, flush=True)

        with open("{}/{}_best_params.json".format(mva_dir, MVA), "w") as f:
            js.dump(bestPar, f, indent=4, cls=NpEncoder)

        main_plot_history_mod(trials, do_show=False, outName="{}/{}_hyperTrails.pdf".format(mva_dir, MVA))

    #====================================================#
    #   train the model and draw the training results    #
    #====================================================#
    bst = Training(bestPar)
    y_train_pred = xgb.dask.predict(client, bst, dtrain)
    y_test_pred = xgb.dask.predict(client, bst, dtest)
    Y_train, y_train_pred, Wt_train = dask.compute(Y_train, y_train_pred, Wt_train)
    Y_test, y_test_pred, Wt_test = dask.compute(Y_test, y_test_pred, Wt_test)

    print(color.GREEN + "Draw the training results!" + color.END, flush=True)
    plot_importance_mod(bst, features=cfg.features, outName="{}/{}_Importance.pdf".format(mva_dir, MVA))
    plot_loss_history(bst, outName="{}/{}_LogLoss.pdf".format(mva_dir, MVA))
    plot_MVAscore(
        df, Y_train, Y_test, y_train_pred, y_test_pred,
        wt="instwei", Classes=cfg.Classes, ClassColors=cfg.ClassColors, isb=isBinary, MVA=MVA,
        outName="{}/{}_MVA.pdf".format(mva_dir, MVA)
    )
    plotwt_train, plotwt_test, WPCuts = plot_ROC(
        df, Y_train, Y_test, y_train_pred, y_test_pred,
        cat="EleType", wt="instwei", Classes=cfg.Classes, isb=isBinary, MVA="XGB",
        outName="{}/{}_ROC.pdf".format(mva_dir, MVA)
    )

    # add the predicted class to dataframe
    if isBinary == True:
        y_train_cat = np.array([1 if i > WPCuts else 0 for i in y_train_pred])
        y_test_cat = np.array([1 if i > WPCuts else 0 for i in y_test_pred])
    else:
        y_train_cat = y_train_pred.argmax(axis=1)
        y_test_cat = y_test_pred.argmax(axis=1)
    df.loc[TrainIndices, predcat] = y_train_cat
    df.loc[TestIndices, predcat] = y_test_cat

    # draw the confusion matrix
    cm_train = confusion_matrix(Y_train, y_train_cat, sample_weight=plotwt_train)
    plot_confusion_matrix(
        cm_train, cfg.Classes, normalize=True,
        title="Normalized confusion matrix(Train sample)",
        outName="{}/{}_Normalized_CM_Training.pdf".format(mva_dir, MVA)
    )
    cm_test = confusion_matrix(Y_test, y_test_cat, sample_weight=plotwt_test)
    plot_confusion_matrix(
        cm_test, cfg.Classes, normalize=True,
        title="Normalized confusion matrix(Test sample)",
        outName="{}/{}_Normalized_CM_Testing.pdf".format(mva_dir, MVA)
    )

    # draw the efficiency of ID
    plot_eff(
        df, cfg.ptbins, trueLable=cat, predLable=predcat,
        var="eleCalibPt", xaxisName="p^{e}_{T} [GeV]", wt="instwei", isb=isBinary,
        Classes=cfg.Classes, ClassColors=cfg.ClassColors,
        outName="{}/{}_eleCalibPt_eff.pdf".format(mva_dir, MVA)
    )

    #====================================================#
    #             print the training results             #
    #====================================================#
    acc = accuracy_score(Y_test, y_test_cat, sample_weight=plotwt_test)
    logloss = min(bst["history"]["test"][cfg.model_params["eval_metric"]])
    print(color.GREEN + "The training results for test sample" + color.END, flush=True)
    print("Expected log loss of the test sample: {:.4f}".format(logloss))
    print("Accuracy of the test sample: {:.4f}".format(acc))
    clf = classification_report(Y_test, y_test_cat, target_names=cfg.Classes, sample_weight=plotwt_test, digits=4)
    print(clf)

    if isBinary == True:
        results = {"logloss": logloss, "accuracy": acc, "threshold": round(np.float(WPCuts), 4)}
    else:
        results = {"logloss": logloss, "accuracy": acc}

    with open("{}/{}_results.json".format(mva_dir, MVA), "w") as f:
        js.dump(results, f, indent=4)


if __name__ == "__main__":
    start_time = time.time()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush=True)
    print(color.BLUE + "---Start to train the ID!---" + color.END, flush=True)

    # get the information from config file
    train_cfg = sys.argv[1]
    print("[INFO] Importing settings: " + train_cfg)
    cfg = import_module(train_cfg.replace("/", "."))

    # global varibales used in the training
    # True lable coulumn name = cat. Predicted lable coulumn name = predcat
    cat, predcat, weight, Clfname, MVA = "EleType", "ElePredType","NewWt", cfg.Clfname, cfg.MVA
    isBinary = True if len(cfg.Classes) == 2 else False

    # setup the directory to put the training results
    print("[INFO] Making output directory: " + cfg.OutputDirName)
    mva_dir = "{}/{}".format(cfg.OutputDirName, MVA)
    feature_dir = "{}/FeaturePlot".format(cfg.OutputDirName)
    if (os.path.exists(cfg.OutputDirName) and len(os.listdir(cfg.OutputDirName)) != 0):
        os.system("rm -rf " + cfg.OutputDirName)
    os.system("mkdir -p " + cfg.OutputDirName)
    os.system("mkdir -p " + mva_dir)
    os.system("mkdir -p " + cfg.OutputDirName+"/CodeANDConfig")
    os.system("mkdir -p " + feature_dir)
    os.system("cp "+train_cfg+".py ./"+ cfg.OutputDirName+"/CodeANDConfig/")
    os.system("cp Trainer-HDalitz-hyper.py ./"+ cfg.OutputDirName+"/CodeANDConfig/")

    print(color.GREEN + "Initialize the cuda cluster!" + color.END)
    with LocalCUDACluster(threads_per_worker=4) as cluster:
        with Client(cluster) as client:
            main()

    print(color.BLUE + "---All done!---" + color.END, flush=True)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))