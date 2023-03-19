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
from dask import dataframe as dd
from dask.distributed import Client
from dask_cuda import LocalCUDACluster
import hyperopt as hpt
from hyperopt import hp
from hyperopt.early_stop import no_progress_loss
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report
from sklearn.utils import shuffle

from tools.plotUtils import *
from tools.dfUtils import df_load, df_split
from tools.hyperUtils import trial_count, NpEncoder
from tools.ptEtaRwt import df_pteta_rwt, df_balance_rwt

import ROOT


# https://github.com/hyperopt/hyperopt/wiki/FMin#21-parameter-expressions
space = {
    "grow_policy":      hp.choice("grow_policy", ["depthwise", "lossguide"]),
    "learning_rate":    hp.loguniform("learning_rate", np.log(0.01), np.log(0.5)),
    "max_depth":        hp.choice("max_depth", np.arange(5, 15)),
    "max_delta_step":   hp.choice("max_delta_step", np.arange(5, 15)),
    "min_child_weight": hp.uniform("min_child_weight", 1, 1000),
    "subsample":        hp.uniform("subsample", 0.1, 1),
    "min_split_loss":   hp.uniform("min_split_loss", 0, 2),
    "colsample_bytree": hp.uniform("colsample_bytree", 0.1, 1),
    "reg_lambda":       hp.loguniform("reg_lambda",  np.log(1),  np.log(10))
}


@trial_count
def objective(params):
    print("Parameters: {}".format(params), flush=True)
    bst = xgb.dask.train(
        client, params, dtrain,
        num_boost_round=cfg.num_boost_round, early_stopping_rounds=cfg.early_stopping_rounds,
        evals=[(dtrain, "train"), (dtest, "test")], verbose_eval=False
    )

    min_idx    = np.argmin(bst["history"]["test"][cfg.model_params["eval_metric"]])
    loss_test  = bst["history"]["test"][cfg.model_params["eval_metric"]][min_idx]
    loss_train = bst["history"]["train"][cfg.model_params["eval_metric"]][min_idx]

    print("{} for test: {}".format(cfg.model_params["eval_metric"], loss_test), flush=True)
    print("{} for train: {}".format(cfg.model_params["eval_metric"], loss_train), flush=True)
    print("", flush=True)

    return loss_test if cfg.model_params["eval_metric"] != "auc" else -loss_test


def Training(params):
    # return the model with best iteration (min loss in testing set) instead of last iteration
    early_stop = xgb.callback.EarlyStopping(
        rounds=cfg.early_stopping_rounds,
        data_name="test",
        save_best=True
    )

    # start training
    bst = xgb.dask.train(
        client, params, dtrain,
        evals=[(dtrain, "train"), (dtest, "test")],
        num_boost_round=cfg.num_boost_round, callbacks=[early_stop]
    )

    print(color.GREEN + "Save the training model" + color.END)
    print("[INFO] Save the model in {}".format("{}/{}_modelXGB.txt".format(mva_dir, MVA)), flush=True)
    bst["booster"].save_model("{}/{}_modelXGB.txt".format(mva_dir, MVA))
    print("[INFO] Save the model in {}".format("{}/{}_modelXGB.pkl".format(mva_dir, MVA)), flush=True)
    pkl.dump(bst["booster"], open("{}/{}_modelXGB.pkl".format(mva_dir, MVA), "wb"))

    return bst


# convert pandas dataframe to tree for tmva nn training
def df2tree(df, dtype="train"):
    data = {key: df[key].values for key in df.columns}
    rdf = ROOT.RDF.MakeNumpyDataFrame(data)
    if dtype == "train":
        rdf.Snapshot("train_tree", "{}/train_tmva.root".format(mva_dir))
    else:
        opt = ROOT.RDF.RSnapshotOptions()
        opt.fMode = "update"
        rdf.Snapshot("test_tree", "{}/train_tmva.root".format(mva_dir), "", opt)


# prepare the dask dataframe for training
def PrepDataset(df, TrainIndices, TestIndices):
    # training dataset
    X_train_    = df.loc[TrainIndices, cfg.features]
    Y_train_    = df.loc[TrainIndices, cat]
    Wt_train_   = df.loc[TrainIndices, weight]
    if cfg.saveROOT:
        df2tree(df.loc[TrainIndices, cfg.features+[weight]+[cat]], dtype="train")

    # testing dataset
    X_test_     = df.loc[TestIndices, cfg.features]
    Y_test_     = df.loc[TestIndices, cat]
    Wt_test_    = df.loc[TestIndices, weight]
    if cfg.saveROOT:
        df2tree(df.loc[TestIndices, cfg.features+[weight]+[cat]], dtype="test")

    return dd.from_pandas(X_train_, npartitions=10), dd.from_pandas(Y_train_, npartitions=10), dd.from_pandas(Wt_train_, npartitions=10),\
           dd.from_pandas(X_test_, npartitions=10), dd.from_pandas(Y_test_, npartitions=10), dd.from_pandas(Wt_test_, npartitions=10)


def main():
    #====================================================#
    #        load the training data and preprocess       #
    #====================================================#
    # Note: Normaliztion and scalar is not needed for xgboost.
    #       I have checked if standard scalr and PCA transformation can imporve the performance.
    #       The performance becomes worse for xgboost.
    #       TMVA also provides these transformation for BDTG and NN training
    # reference: https://github.com/dmlc/xgboost/issues/357
    print(color.GREEN + "Loading the training dataframe!" + color.END, flush=True)
    sig_df = df_load(cfg.Parquet_signal, columns=cfg.columns, cuts=cfg.CommonCut, extral_text="signal")
    bkg_df = df_load(cfg.Parquet_bkg, columns=cfg.columns, cuts=cfg.CommonCut, extral_text="background")
    if isBinary == True:
        sig_df["Class"], bkg_df["Class"]  = "Signal", "Background"
    df = shuffle(pd.concat([sig_df, bkg_df], ignore_index=True, sort=False), random_state=cfg.RandomState)

    if cfg.uniformwt == True:
        df["instwei"] = 1.

    # split df to taining and testing sample
    TrainIndices, TestIndices, df = df_split(
        df, test_size=cfg.testsize, seed=cfg.RandomState,
        isb=isBinary, trueLable=cat, Classes=cfg.Classes
    )

    #====================================================#
    # Draw the training feature distributions and correl #
    #====================================================#
    print(color.GREEN + "Draw the training feature plots" + color.END, flush=True)
    with open("./tools/{}".format(cfg.featureplotparam_json), "r") as fp:
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
    plot_corr(df, cfg.features, Classes=cfg.Classes, MVA=MVA, outdir=mva_dir)
    #====================================================#
    #                 add training weight                #
    #====================================================#
    if cfg.Reweighing == "Balanced":
        print(color.GREEN + "Balanced reweighting for training sample" + color.END, flush=True)
        df.loc[TrainIndices, weight] = df_balance_rwt(
            df.loc[TrainIndices],
            SumWeightCol="instwei", NewWeightCol=weight,
            Classes=cfg.Classes
        )

        print(color.GREEN + "Balanced reweighting for testing sample" + color.END, flush=True)
        df.loc[TestIndices, weight] = df_balance_rwt(
            df.loc[TestIndices],
            SumWeightCol="instwei", NewWeightCol=weight,
            Classes=cfg.Classes
        )
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

    # plot_stat(df, outname="{}/TotalStat_TrainANDTest".format(mva_dir))

    #====================================================#
    #     change the dtype to dask array for training    #
    #====================================================#
    global X_train, Y_train, Wt_train, X_test, Y_test, Wt_test
    X_train, Y_train, Wt_train, X_test, Y_test, Wt_test = PrepDataset(df, TrainIndices, TestIndices)

    global dtrain, dtest
    dtrain = xgb.dask.DaskDMatrix(client, X_train, Y_train, weight=Wt_train)
    dtest = xgb.dask.DaskDMatrix(client, X_test, Y_test, weight=Wt_test)

    #====================================================#
    #              hyperparameter searching              #
    #====================================================#
    bestPar = {**cfg.model_params, **cfg.hyper_params}
    if cfg.opt == True:
        print(color.GREEN + "Perform the hyper parameter searching" + color.END, flush=True)
        
        trials = hpt.Trials()
        search_params = {**cfg.model_params, **space}
        
        # atpe, tpe
        best = hpt.fmin(
            objective, search_params, trials=trials,
            algo=hpt.tpe.suggest, max_evals=500, rstate=np.random.RandomState(cfg.RandomState),
            early_stop_fn=no_progress_loss(iteration_stop_count=15, percent_increase=0.0001),
            show_progressbar=False
        )
        minLoss = min(trials.losses())
        bestPar = hpt.space_eval(search_params, best)
        print(color.GREEN + "Best hyperparameters with {} = {}: ".format(search_params["eval_metric"], minLoss) + color.END, flush=True)
        print(bestPar, flush=True)

        with open("{}/{}_best_params.json".format(mva_dir, MVA), "w") as f:
            js.dump(bestPar, f, indent=4, cls=NpEncoder)

        main_plot_history_mod(trials, metrix=cfg.model_params["eval_metric"], do_show=False, outName="{}/{}_hyperTrails.pdf".format(mva_dir, MVA))

    #====================================================#
    #   train the model and draw the training results    #
    #====================================================#
    bst = Training(bestPar)

    y_train_pred = xgb.dask.predict(client, bst, dtrain).compute()
    y_test_pred = xgb.dask.predict(client, bst, dtest).compute()

    Y_train, Wt_train = Y_train.compute(), Wt_train.compute()
    Y_test, Wt_test = Y_test.compute(), Wt_test.compute()
    X_train = X_train.compute()

    print(color.GREEN + "Draw the training results!" + color.END, flush=True)
    plot_importance_mod(bst, features=cfg.features, outName="{}/{}_Importance.pdf".format(mva_dir, MVA))

    # https://www.kaggle.com/code/bextuychiev/model-explainability-with-shap-only-guide-u-need
    # expected_vale(base_vale) is the average model score based on the provided training set
    # expected_vale ~ bst["booster"].predict(dtrain).mean()
    # calculating shap values is time consuming -> may run out of GPU memeory for high statistics :(
    if cfg.doShapPlot:
        shap_values = xgb.dask.predict(client, bst, dtrain, pred_contribs=True).compute()
        shap_values = [np.array(shap_values[i], dtype=np.float32)[:, :-1] for i in range(len(cfg.Classes))]
        plot_shap(
            shap_values, X_train,
            cmapList=cfg.ClassColors, Classes=cfg.Classes, features=cfg.features,
            outName="{}/{}_shap_summary.pdf".format(mva_dir, MVA)
        )
    plot_loss_history(bst, metric=cfg.model_params["eval_metric"], outName="{}/{}_loss.pdf".format(mva_dir, MVA))

    # Note: EGM ID convenors said we should NOT include the "unphysical" weight when drawing ROC and MVA distribution.
    #       "unphysical" weight: such as balanced weight or pt-eta reweight.
    #       However, in order to compare with the TMVA NN results, i still include these weights to visualize the results.
    #       The other reason that why i still include "unphysical" weight is to make sure the training is overtrained or not from the auc of ROC.
    #       wt=weight -> wt="instwei" to remove the "unphysical" weight when visulization,
    #       and the returned weight array (plotwt_train, plotwt_test) will also be the instwei extracted from dataframe.
    plot_MVAscore(
        df, Y_train, Y_test, y_train_pred, y_test_pred,
        wt=weight, Classes=cfg.Classes, ClassColors=cfg.ClassColors, isb=isBinary, MVA=MVA,
        outName="{}/{}_MVA.pdf".format(mva_dir, MVA)
    )
    plotwt_train, plotwt_test, WPCuts = plot_ROC( # WPCuts is useless for multiclass classification
        df, Y_train, Y_test, y_train_pred, y_test_pred,
        cat="EleType", wt=weight, Classes=cfg.Classes, isb=isBinary, MVA="XGB",
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
    # Note: use instwei instead of including "unphysical" weight is to make sure the ID efficiency.
    #       Since we will not include "unphysical" weight when measureing the ID efficiency.
    #       Merged ID efficiency measurements can be found in the ../misc/MergedIDEff directory.
    plot_eff(
        df, cfg.ptbins, trueLable=cat, predLable=predcat,
        var=cfg.ptwtvar, xaxisName="p^{e}_{T} [GeV]", wt="instwei", isb=isBinary,
        Classes=cfg.Classes, ClassColors=cfg.ClassColors,
        outName="{}/{}_{}_eff.pdf".format(mva_dir, MVA, cfg.ptwtvar)
    )

    # ====================================================#
    #             print the training results             #
    # ====================================================#
    acc = accuracy_score(Y_test, y_test_cat, sample_weight=plotwt_test)
    logloss = min(bst["history"]["test"][cfg.model_params["eval_metric"]])
    print(color.GREEN + "The training results for test sample" + color.END, flush=True)
    print("Expected {} of the test sample: {:.4f}".format(cfg.model_params["eval_metric"], logloss))
    print("Accuracy of the test sample: {:.4f}".format(acc))
    clf = classification_report(Y_test, y_test_cat, target_names=cfg.Classes, sample_weight=plotwt_test, digits=4)
    print(clf)

    if isBinary == True:
        results = {cfg.model_params["eval_metric"]: logloss, "accuracy": acc, "threshold": round(np.float(WPCuts), 4)}
    else:
        results = {cfg.model_params["eval_metric"]: logloss, "accuracy": acc, "threshold": round(np.float(WPCuts), 4)}

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
    with LocalCUDACluster() as cluster: # device_memory_limit="auto", threads_per_worker=10
        with Client(cluster) as client: # , timeout=1200
            main()

    print(color.BLUE + "---All done!---" + color.END, flush=True)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))