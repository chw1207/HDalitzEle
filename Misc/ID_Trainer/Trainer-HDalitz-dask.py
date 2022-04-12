import sys, os
import subprocess
import json
import pickle
import pandas as pd
import numpy as np
import time
import matplotlib.pyplot as plt

import xgboost as xgb
import dask.dataframe as dd
import dask
from dask import array as da
from dask.distributed import Client
from dask_cuda import LocalCUDACluster

from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest
from sklearn.metrics import confusion_matrix, roc_curve, roc_auc_score, accuracy_score, auc, log_loss, classification_report
from Tools.PlotTools import *

# Some useful links
# *     1) Dask: https://medium.com/rapids-ai/a-new-official-dask-api-for-xgboost-e8b10f3d1eb7
# *     2) Weight: https://stackoverflow.com/questions/35983565/how-is-the-parameter-weight-dmatrix-used-in-the-gradient-boosting-procedure
# *     3) Feature Selection: https://amueller.github.io/aml/05-advanced-topics/12-feature-selection.html
# *     4) Samplings: https://stackoverflow.com/a/66125815

#! To-Do:
# *     1) hyper parameters: https://www.kaggle.com/yassinealouini/hyperopt-the-xgboost-model
# *     2) plots of bkg rejection and signal efficiency as fuction of electron pT


# Load the pkl files (type: pandas Dataframe)
def Load_df():
    print (color.BLUE + "---Loading background dataframes from: {}---".format(Conf.Pickle_bkg) + color.END)

    # load the signal dataframe
    with open(Conf.Pickle_bkg, "rb") as f:
        if Conf.CommonCut:
            print("Select the background events/objects with commen cuts: {}".format(Conf.CommonCut))
            Bkgdf = pickle.load(f).query(Conf.CommonCut)
        else:
            Bkgdf = pickle.load(f)

    print (color.BLUE + "---Loading signal dataframes from: {}---".format(Conf.Pickle_signal) + color.END)

    # load the background dataframe
    with open(Conf.Pickle_signal, "rb") as f:
        if Conf.CommonCut:
            print ("Select the signal events/objects with commen cuts: {}".format(Conf.CommonCut))
            Sigdf = pickle.load(f).query(Conf.CommonCut)
        else:
            Sigdf = pickle.load(f)

    df = pd.concat([Sigdf, Bkgdf], ignore_index = True, sort = False)

    # change the dtype
    cols = df.select_dtypes(include = [object]).columns.drop("Class")
    df[cols] = df[cols].astype(np.float32)

    print (color.BLUE + "---Dataframes loading done!---" + color.END)

    return df


# split the dataframe to train and test sample
def SplitTestTrain(df):
    df[cat] = 0
    for i, k in enumerate(Conf.Classes):
        df.loc[df.Class == k, cat] = i

    index = df.index
    TrainIndices, TestIndices = [], []
    for myclass in Conf.Classes:
        Indices = index[df["Class"] == myclass].values.tolist()
        myclassTrainIndices, myclassTestIndices = train_test_split(Indices, test_size = Conf.testsize, random_state = Conf.RandomState, shuffle = True)
        TrainIndices = TrainIndices + myclassTrainIndices
        TestIndices = TestIndices + myclassTestIndices

    df.loc[TrainIndices, "Dataset"] = "Train"
    df.loc[TestIndices, "Dataset"] = "Test"

    df.loc[TrainIndices, "TrainDataset"] = 1
    df.loc[TestIndices, "TrainDataset"] = 0

    prGreen("Reading classes:")
    print(df.Class.unique().tolist())

    return TrainIndices, TestIndices, df


# deal with pT-eta reweighting
def df_pteta_rwt(Mdf, label, returnOnlyPosWeights = 0, ptw = [10, 30, 40, 50, 200, 10000], etaw = [-1.5, -1.0, 1.0, 1.5], eta = "", pt = "", SumWeightCol = "wt", NewWeightCol = "NewWt", cand = "", Classes = [""]):
    Mdf["rwt"] = 1
    Mdf[NewWeightCol] = 1

    ptwt = [1.0] * len(ptw)
    etawt = [1.0] * len(etaw)

    for k in range(len(etaw)):
        if k == len(etaw) - 1:
            continue
        for i in range(len(ptw)):
            if i == len(ptw) - 1:
                continue
            for target in Classes:
                if target != cand:
                    targetSum = Mdf.loc[(Mdf[pt] < ptw[i+1]) & (Mdf[pt] > ptw[i]) & (Mdf[eta] < etaw[k+1]) & (Mdf[eta] > etaw[k]) & (Mdf[label] == target), SumWeightCol].sum()
                    candSum = Mdf.loc[(Mdf[pt] < ptw[i+1]) & (Mdf[pt] >ptw[i]) & (Mdf[eta] < etaw[k+1]) & (Mdf[eta] > etaw[k]) & (Mdf[label] == cand), SumWeightCol].sum()

                    if candSum > 0 and targetSum > 0:
                        ptwt[i] = candSum / (targetSum)
                    else:
                        ptwt[i] = 0

                    Mdf.loc[(Mdf[pt] < ptw[i+1]) & (Mdf[pt] > ptw[i]) & (Mdf[eta] < etaw[k+1]) & (Mdf[eta] > etaw[k]) & (Mdf[label] == cand), "rwt"] = 1.0
                    Mdf.loc[(Mdf[pt] < ptw[i+1]) & (Mdf[pt] > ptw[i]) & (Mdf[eta] < etaw[k+1]) & (Mdf[eta] > etaw[k]) & (Mdf[label] == target), "rwt"] = ptwt[i]

    Mdf.loc[:, NewWeightCol] = Mdf.loc[:, "rwt"] * Mdf.loc[:, SumWeightCol]

    for justclass in Classes:
        Sum = Mdf.loc[Mdf[label] == justclass, NewWeightCol].sum()
        print("Number of events in {} after pt-eta reweighting = {:.2f}".format(justclass, Sum))

    return Mdf[NewWeightCol]


# unbalanced multi-classification: https://datascience.stackexchange.com/a/49067
def AddWeight(df):
    if Conf.Reweighing == "Balanced":
        df[weight] = 1
        sum_w, wei = [1.0 for i in range(len(Conf.Classes))], [1.0 for i in range(len(Conf.Classes))]
        for i, k in enumerate(Conf.Classes):
            sum_w[i] = df["instwei"][df.Class == k].sum()
        minimum = min(sum_w)
        for i, k in enumerate(Conf.Classes):
            wei[i] = minimum / sum_w[i]
            df.loc[df.Class == k, "rwt"] = wei[i]
            df.loc[:, weight] = df.loc[:, "rwt"] * df.loc[:, "instwei"]
            print("Class = %s, n = %.2f, balanced weight = %.2e" %(k, sum_w[i], wei[i]))

    elif Conf.Reweighing in Conf.Classes:
        # pt-eta reweighting for training set
        df.loc[TrainIndices, weight] = df_pteta_rwt(df.loc[TrainIndices], "Class", ptw = Conf.ptbins, etaw = Conf.etabins, pt = Conf.ptwtvar, eta = Conf.etawtvar, SumWeightCol = "instwei", NewWeightCol = weight, cand = Conf.Reweighing, Classes = Conf.Classes)
        plotPtEtaRwt(df, Conf.ptwtvar, Conf.etawtvar, Conf.ptbins, Conf.etabins, Conf.OutputDirName+"/"+MVA)

        # pt-eta reweighting for testing set
        df.loc[TestIndices, weight] = df_pteta_rwt(df.loc[TestIndices], "Class", ptw = Conf.ptbins, etaw = Conf.etabins, pt = Conf.ptwtvar, eta = Conf.etawtvar, SumWeightCol = "instwei", NewWeightCol = weight, cand = Conf.Reweighing, Classes = Conf.Classes)

    else:
        df[weight] = 1
        print("This reweighting method is not yet available!")
        sys.exit(1)

    return df


# prepare the dask array for training
# Note: use scaler after train_test_split, otherwise data leakage will happen
def PrepDataset(df, TrainIndices, TestIndices, features, cat, weight):
    X_train = np.asarray(df.loc[TrainIndices, features])
    Y_train = np.asarray(df.loc[TrainIndices, cat])
    Wt_train = np.asarray(df.loc[TrainIndices, weight])

    X_test = np.asarray(df.loc[TestIndices, features])
    Y_test = np.asarray(df.loc[TestIndices, cat])
    Wt_test = np.asarray(df.loc[TestIndices, weight])

    return da.from_array(X_train, chunks = 5000), da.from_array(Y_train, chunks = 5000), da.from_array(Wt_train, chunks = 5000), da.from_array(X_test, chunks = 5000), da.from_array(Y_test, chunks = 5000), da.from_array(Wt_test, chunks = 5000)


# trainging model
# Reference: https://github.com/dmlc/xgboost/blob/master/demo/dask/gpu_training.py
def Training(client: Client):
    dtrain = xgb.dask.DaskDMatrix(client, X_train, Y_train, weight = Wt_train)
    dtest = xgb.dask.DaskDMatrix(client, X_test, Y_test, weight = Wt_test)

    prGreen("Performing XGB model training!")
    bst = xgb.dask.train(
        client, Conf.param, dtrain,
        num_boost_round = Conf.num_boost_round, early_stopping_rounds = Conf.early_stopping_rounds,
        evals = [(dtrain, "train"), (dtest, "test")]
    )

    prGreen("Save the training model")
    bst["booster"].save_model(Conf.OutputDirName + "/" + MVA + "/" + MVA + "_" + "modelXGB.json")
    print("Save the model in {}".format(Conf.OutputDirName + "/" + MVA + "/" + MVA + "_" + "modelXGB.json"))

    pickle.dump(bst["booster"], open(Conf.OutputDirName + "/" + MVA + "/" + MVA + "_" + "modelXGB.pkl", "wb"))
    print("Save the model in {}".format(Conf.OutputDirName + "/" + MVA + "/" + MVA + "_" + "modelXGB.pkl"))

    y_train_pred = xgb.dask.predict(client, bst, dtrain)
    y_test_pred = xgb.dask.predict(client, bst, dtest)

    return bst, y_train_pred, y_test_pred


# copy from the TensorFlow to_categoriacal function
# Reference: https://github.com/keras-team/keras/blob/v2.7.0/keras/utils/np_utils.py#L21-L74
def ToCategorical(y, num_classes = None, dtype = "float32"):
    y = np.array(y, dtype='int')
    input_shape = y.shape
    if input_shape and input_shape[-1] == 1 and len(input_shape) > 1:
        input_shape = tuple(input_shape[:-1])
    y = y.ravel()
    if not num_classes:
        num_classes = np.max(y) + 1
    n = y.shape[0]
    categorical = np.zeros((n, num_classes), dtype=dtype)
    categorical[np.arange(n), y] = 1
    output_shape = input_shape + (num_classes,)
    categorical = np.reshape(categorical, output_shape)

    return categorical


if __name__ == "__main__":
    start_time = time.time()

    TrainConfig = sys.argv[1]

    prGreen("Importing settings from "+ TrainConfig.replace("/", "."))
    importConfig = TrainConfig.replace("/", ".")
    exec("import "+importConfig+" as Conf")
    cat, weight, label, Clfname, MVA = "EleType", "NewWt", Conf.Classes, Conf.Clfname, Conf.MVA

    print (color.BOLD + color.BLUE + "---Making output directory: {}---".format(Conf.OutputDirName)+ color.END)
    if (os.path.exists(Conf.OutputDirName) and len(os.listdir(Conf.OutputDirName)) != 0):
        os.system("rm -rf " + Conf.OutputDirName)
    os.system("mkdir -p " + Conf.OutputDirName)
    os.system("mkdir -p " + Conf.OutputDirName+"/"+MVA)
    os.system("mkdir -p " + Conf.OutputDirName+"/CodeANDConfig")
    os.system("mkdir -p " + Conf.OutputDirName+"/FeaturePlot")
    os.system("cp "+TrainConfig+".py ./"+ Conf.OutputDirName+"/CodeANDConfig/")
    os.system("cp Trainer-HDalitz-dask.py ./"+ Conf.OutputDirName+"/CodeANDConfig/")

    # pre-processing the dataframe
    df = Load_df()
    TrainIndices, TestIndices, df = SplitTestTrain(df)
    df = AddWeight(df)
    prGreen("Draw the statistics plots")
    DrawStatistic(df, Conf.OutputDirName)

    prGreen("Draw the training feature plots")
    os.makedirs("{}/FeaturePlot/".format(Conf.OutputDirName), exist_ok = True)
    with open("./Tools/{}".format(Conf.featureplotparam_json), "r") as fp:
        myTags = json.load(fp)

        for tagName, par in myTags.items():
            setlogy = True if par["logy"] == "true" else False

            if tagName not in Conf.features:
                continue

            DrawFeatureHist_multi(
                df, tagName, len(Conf.Classes), Conf.Classes,
                par["bininfo"], par["xAxisLabel"], par["yAxisUnit"],
                histcolor = Conf.ClassColors,
                logy = setlogy, y_axisscale = par["y_axisscale"],
                outname = "{}/FeaturePlot/TrainFeature_{}".format(Conf.OutputDirName, tagName),
                wei = "instwei",
                category = cat
            )

    # change the dtype to dask array for training
    X_train, Y_train, Wt_train, X_test, Y_test, Wt_test = PrepDataset(df, TrainIndices, TestIndices, Conf.features, cat, weight)

    # start training
    # n_workers: all GPUs(default), threads_per_worker: # of threads to be used for each worker
    prGreen("Initialize the cuda cluster......")
    with LocalCUDACluster() as cluster:
        with Client(cluster) as client:
            bst, y_train_pred, y_test_pred = Training(client)

            # compute all of the dask array
            Y_train, y_train_pred, Wt_train = dask.compute(Y_train, y_train_pred, Wt_train)
            Y_test, y_test_pred, Wt_test = dask.compute(Y_test, y_test_pred, Wt_test)

            # print the training results
            logloss = log_loss(Y_test, y_test_pred, sample_weight = Wt_test)
            acc = accuracy_score(Y_test, y_test_pred.argmax(axis = 1), sample_weight = Wt_test)
            AUC = roc_auc_score(Y_test, y_test_pred, sample_weight = Wt_test, multi_class = "ovr", average = "weighted")
            prGreen("The training results for test sample......")
            print("Expected log loss of the test sample: {:.4f}".format(logloss))
            print("Accuracy of the test sample: {:.4f}".format(acc))
            print("The average roc auc score of test sample is {:.4f}".format(AUC))
            print(classification_report(Y_test, y_test_pred.argmax(axis = 1), target_names = Conf.Classes, sample_weight = Wt_test))

            results = {"logloss": logloss, "accuracy": acc, "auc": AUC}
            results_path = Conf.OutputDirName+"/"+MVA+"/"+MVA+"_"+"results.json"
            with open(results_path, "w") as f:
                json.dump(results, f)

            # visualize the results
            prGreen("Draw the training results......")
            #! 1) feature importance
            bst["booster"].feature_names = Conf.features
            importance = bst["booster"].get_score(importance_type = "gain")
            for key in importance.keys():
                importance[key] = round(importance[key], 2)
            max_importance = max(list(importance.values()))

            fig, ax = plt.subplots(1, 1, figsize = (6, 7))
            xgb.plot_importance(importance, ax = ax, importance_type = "gain", height = 0.7, grid = False, title = None, xlabel = None, ylabel = None, xlim = (0, max_importance * 1.2), color = "#5893D4")
            pltSty(ax, xName = "Importance", yName = "Features", yAuto = False)
            fig.savefig(Conf.OutputDirName+"/"+MVA+"/"+MVA+"_"+"Importance.pdf", bbox_inches = "tight")
            print("Save fig in {}".format(Conf.OutputDirName+"/"+MVA+"/"+MVA+"_"+"Importance.pdf"))
            plt.close("all")

            #! 2) confusion metrices
            cm_train = confusion_matrix(Y_train, y_train_pred.argmax(axis = 1), sample_weight = Wt_train)
            plot_confusion_matrix(
                cm_train, Conf.Classes, normalize = True,
                title = "Normalized confusion matrix(Train sample)",
                outName = Conf.OutputDirName+"/"+MVA+"/"+MVA+"_Normalized_CM_Training"+".pdf"
            )
            cm_test = confusion_matrix(Y_test, y_test_pred.argmax(axis = 1), sample_weight = Wt_test)
            plot_confusion_matrix(
                cm_test, Conf.Classes, normalize = True,
                title = "Normalized confusion matrix(Test sample)",
                outName = Conf.OutputDirName+"/"+MVA+"/"+MVA+"_Normalized_CM_Testing"+".pdf"
            )

            #! 3) history of epochs
            eval_metric = Conf.param["eval_metric"]
            fig, axes = plt.subplots(1, 1, figsize = (6, 6))
            axes.plot(range(1, len(bst["history"]["train"][eval_metric]) + 1), bst["history"]["train"][eval_metric], linewidth = 3, linestyle = "-", color = "#184d47", label = "Training loss")
            axes.plot(range(1, len(bst["history"]["test"][eval_metric]) + 1), bst["history"]["test"][eval_metric], linewidth = 3, linestyle = "--", color = "#e2703a", label = "Testing loss")
            axes.legend(loc = "best", edgecolor= "none", fontsize = 15, title_fontsize = 15)
            pltSty(axes, xName = "Epochs", yName = "Log loss")
            fig.savefig(Conf.OutputDirName+"/"+MVA+"/"+MVA+"_"+"LogLoss.pdf", bbox_inches = "tight")
            print("Save fig in %s" %(Conf.OutputDirName+"/"+MVA+"/"+MVA+"_"+"LogLoss.pdf"))
            plt.close("all")

            #! 4) MVA score and ROC
            # * ROC threshold: -> threshold seems rebundent for multi-classification ?
            # *     1) https://towardsdatascience.com/demystifying-roc-curves-df809474529a
            # *     2) https://machinelearningmastery.com/threshold-moving-for-imbalanced-classification/
            Y_train_categorical = ToCategorical(Y_train, num_classes = len(Conf.Classes))
            Y_test_categorical  = ToCategorical(Y_test, num_classes = len(Conf.Classes))
            n_classes = len(Conf.Classes)
            figMVA, axesMVA = plt.subplots(1, n_classes, figsize = (n_classes*8.5, 7))
            fig, axes = plt.subplots(1, n_classes, figsize = (n_classes*8.5, 7))
            for i in range(n_classes):
                axMVA = axesMVA[i]
                ax = axes[i]
                for k in range(n_classes):
                    axMVA.hist(
                        y_test_pred[:, i][Y_test_categorical[:, k] == 1],
                        bins = np.linspace(0, 1, num = 40), label = Conf.Classes[k]+": Test",
                        weights = Wt_test[Y_test_categorical[:, k] == 1]/np.sum(Wt_test[Y_test_categorical[:, k] == 1]),
                        histtype = "step", linewidth = 2, color = Conf.ClassColors[k]
                    )
                    axMVA.hist(
                        y_train_pred[:, i][Y_train_categorical[:, k]==1],
                        bins = np.linspace(0, 1, num = 40), label = Conf.Classes[k]+": Train",
                        weights = Wt_train[Y_train_categorical[:, k] == 1]/np.sum(Wt_train[Y_train_categorical[:, k] == 1]),
                        histtype = "stepfilled", alpha = 0.5, linewidth = 2, color = Conf.ClassColors[k]
                    )

                pltSty(axMVA, xName = "{}_pred".format(MVA), yName = "Normalized entries", TitleSize = 16, LabelSize = 17, TickSize = 17, MajTickLength = 10, MinTickLength = 6)
                axMVA.legend(title = "{} vs Rest".format(Conf.Classes[i]), loc = "upper center", fontsize = 15, ncol = 2, edgecolor = "none", title_fontsize = 17)
                axMVA.set_yscale("log")
                axMVA.set_ylim([1E-5, 800])
                axMVA.set_xlim([-0.05, 1.05])

                # tpr: true positive rate = signal efficiency
                # fpr: false positive rate = background efficiency
                fpr, tpr, th = roc_curve(Y_test_categorical[:, i], y_test_pred[:, i], sample_weight = Wt_test)
                fpr_tr, tpr_tr, th_tr = roc_curve(Y_train_categorical[:, i], y_train_pred[:, i], sample_weight = Wt_train)

                bkgrej = np.array([1 - x for x in fpr])
                bkgrej_tr = np.array([1 - x for x in fpr_tr])

                # calculate the g-mean for each threshold
                # if (i == 0):
                #     print("Assume the first class is signal!")
                #     gmeans = np.sqrt(tpr_tr * (1 - fpr_tr))
                #     ix = np.argmax(gmeans)
                #     print("Best Threshold = %.3f, G-Mean = %.3f" %(th_tr[ix], gmeans[ix]))

                roc_auc = auc(fpr, tpr)
                roc_auc_tr = auc(fpr_tr, tpr_tr)

                ax.plot(tpr_tr, bkgrej_tr, label = "XGB Training auc = %.1f%s" %(roc_auc_tr*100, "%"),linewidth = 3, linestyle = "-", color = "#064635", zorder = 1)
                ax.plot(tpr, bkgrej, label = "XGB Testing auc = %.1f%s"% (roc_auc*100, "%"), linewidth = 3, linestyle = "--", color = "#e2703a", zorder = 2)
                # if (i == 0):
                #     ax.scatter(tpr_tr[ix], bkgrej_tr[ix], marker = "o", color = "#202020", label = "Threshold with G-Mean = %.3f" %(gmeans[ix]), s = 100, zorder = 3)
                pltSty(ax, xName = "Signal efficiency", yName = "Background rejection", TitleSize = 16, LabelSize = 17, TickSize = 17, MajTickLength = 10, MinTickLength = 6)

                ax.legend(title = "{} vs Rest".format(Conf.Classes[i]), loc = "best", edgecolor = "none", facecolor = "none", fontsize = 17, title_fontsize = 17)

            fig.savefig(Conf.OutputDirName+"/"+MVA+"/"+MVA+"_"+"ROC.pdf", bbox_inches = "tight")
            print("Save fig in {}".format(Conf.OutputDirName+"/"+MVA+"/"+MVA+"_"+"ROC.pdf"))

            figMVA.savefig(Conf.OutputDirName+"/"+MVA+"/"+MVA+"_"+"MVA.pdf", bbox_inches = "tight")
            print("Save fig in {}".format(Conf.OutputDirName+"/"+MVA+"/"+MVA+"_"+"MVA.pdf"))

    print(color.BOLD + color.BLUE + "---All done---!" + color.END)
    seconds = time.time() - start_time
    print("Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))