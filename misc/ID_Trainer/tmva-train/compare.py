import ROOT
import pandas as pd
import numpy as np
from sklearn.metrics import roc_curve, auc, confusion_matrix
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator, AutoMinorLocator
import itertools
import xgboost as xgb

def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          title = "Confusion matrix",
                          cmap = plt.cm.Blues,
                          outName = ""):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype("float") / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print("Confusion matrix, without normalization")

    plt.figure(figsize = (8, 6))
    plt.imshow(cm, interpolation = "nearest", cmap = cmap)
    plt.title(title, fontsize = 15)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize = 10)
    if normalize:
        plt.clim(0, 1)
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, fontsize = 13)
    plt.yticks(tick_marks, classes, fontsize = 13)

    fmt = ".3f" if normalize else "d"
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment="center", fontsize=13,
                 color="white" if cm[i, j] > thresh else "black")

    plt.ylabel("True class", fontsize=13, loc="top")
    plt.xlabel("Predicted class", fontsize=13, loc="right")
    plt.savefig("{}".format(outName), bbox_inches="tight")
    print("[INFO] Save fig in {}".format(outName), flush=True)
    plt.close("all")


def ToCategorical(y, num_classes=None, dtype="float64"):
    y = np.array(y, dtype="int")
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


def pltSty(ax, xName = "x-axis", yName = "y-axis", TitleSize = 15, LabelSize = 15, TickSize = 13, MajTickLength = 7, MinTickLength = 4, yAuto = True):
    ax.set_xlabel(xName, fontsize = LabelSize, loc = "right")
    ax.set_ylabel(yName, fontsize = LabelSize, loc = "top")
    ax.text(1, 1, "(13 TeV)", horizontalalignment = "right", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize)
    ax.text(0, 1, "CMS", horizontalalignment = "left", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize * 1.3, fontweight = "bold")
    ax.text(TitleSize * 0.009 + 0.02, 1, "Preliminary", horizontalalignment = "left", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize)

    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if (yAuto):
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction = "in", length = MajTickLength, labelsize = TickSize, top = True, right = True)
    ax.tick_params(direction = "in", length = MinTickLength, which = "minor", labelsize = TickSize, top = True, right = True)


def plot_roc(y_cat, y_pred, wt, y_cat_bdtg, y_pred_bdtg, wt_bdtg, y_cat_xgb, y_pred_xgb, wt_xgb, outName):
    # fig, axes = plt.subplots(1, len(Classes), figsize=(len(Classes)*8.5, 7))
    
    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    for i in range(1):
        # ax = axes[i]

        # tpr: true positive rate = signal efficiency
        # fpr: false positive rate = background efficiency
        fpr, tpr, th = roc_curve(
            y_cat[:, i], y_pred[:, i],
            sample_weight=wt
        )
        bkgrej = (1 - fpr)
        roc_auc = auc(fpr, tpr)
        
        fpr_bdtg, tpr_bdtg, th_bdtg = roc_curve(
            y_cat_bdtg[:, i], y_pred_bdtg[:, i],
            sample_weight=wt_bdtg
        )
        bkgrej_bdtg = (1 - fpr_bdtg)
        roc_auc_bdtg = auc(fpr_bdtg, tpr_bdtg)

        fpr_xgb, tpr_xgb, th = roc_curve(
            y_cat_xgb[:, i], y_pred_xgb[:, i],
            sample_weight=wt_xgb
        )
        bkgrej_xgb = (1 - fpr_xgb)
        roc_auc_xgb = auc(fpr_xgb, tpr_xgb)

        # ax.set_yscale("log")
        ax.plot(tpr_xgb, bkgrej_xgb, label="XGBoost testing auc = %.1f%s"%(roc_auc_xgb*100, "%"), linewidth=2.5, linestyle="-", color="#C70039", zorder=1)
        ax.plot(tpr, bkgrej, label="TMVA DNN testing auc = %.1f%s"%(roc_auc*100, "%"), linewidth=2.5, linestyle="-", color="#00A88F", zorder=2)
        ax.plot(tpr_bdtg, bkgrej_bdtg, label="TMVA BDTG testing auc = %.1f%s"%(roc_auc_bdtg*100, "%"), linewidth=2.5, linestyle="-", color="#202020", zorder=2)
        

        pltSty(ax, xName="Signal efficiency", yName="Background rejection", TitleSize=16, LabelSize=17, TickSize=17, MajTickLength=10, MinTickLength=6)
        ax.legend(title="{} vs Rest".format(Classes[i]), loc="best", edgecolor="none", facecolor="none", fontsize=17, title_fontsize=17)

    fig.savefig(outName, bbox_inches="tight")
    plt.close("all")

def main(tree):
    ROOT.EnableImplicitMT(20)
    features = [
        "rho",
        "eleSCEta",
        "eleSCRawEn",

        "eledEtaAtVtx",
        "eledPhiAtVtx",
        "elePtError",
        "eleHoverE",
        "eleEoverP",
        "eleEoverPout",
        "eleEoverPInv",

        "eleSCEtaWidth",
        "eleSCPhiWidth",
        "eleSigmaIEtaIEtaFull5x5",
        "eleSigmaIPhiIPhiFull5x5",
        "eleR9Full5x5",
        "eleBrem",

        "elePFChIso",
        "elePFPhoIso",
        "elePFNeuIso",

        "gsfPtRatio",
        "gsfDeltaR",
        "gsfRelPtRatio",

        # "eleIDMVAIso"
    ]
    rdf = ROOT.RDataFrame(tree, "train_results_EE.root")
    npy = rdf.AsNumpy(columns = ["classID", "weight", "DNN.sig0", "DNN.bkg1", "DNN.bkg2", "BDTG.sig0", "BDTG.bkg1", "BDTG.bkg2"] + features)
    df = pd.DataFrame(npy)

    y = npy["classID"]
    w = npy["weight"]
    y_pred_DNN = np.vstack([npy[var] for var in ["DNN.sig0", "DNN.bkg1", "DNN.bkg2"]]).T
    y_pred_BDTG = np.vstack([npy[var] for var in ["BDTG.sig0", "BDTG.bkg1", "BDTG.bkg2"]]).T

    if tree == "dataset/TrainTree":
        ext1 = "Train"
        ext2 = "Training"
    else:
        ext1 = "Test"
        ext2 = "Testing"

    # cm = confusion_matrix(y, y_pred.argmax(axis=1), sample_weight=w, normalize="true")
    # plot_confusion_matrix(
    #     cm, ["Merged-2Gsf", "DYJets", "QCD"], normalize=True,
    #     title="Normalized confusion matrix({} sample)".format(ext1),
    #     outName = "DNN_Normalized_CM_{}.pdf".format(ext2)
    # )

    # for XGB
    bst = xgb.Booster()
    bst.load_model("/data4/chenghan/external/MergedID/Output_Merged2GsfID_hyperTune_FullRun2ULWPPtWeiFinal_EE/XGB/XGB_modelXGB.txt")
    dmatrix = xgb.DMatrix(df.loc[:, features].values)
    y_pred_xgb = bst.predict(dmatrix)

    # for DNN
    y_cat = ToCategorical(y, num_classes=len(Classes))
    plot_roc(y_cat, y_pred_DNN, w, y_cat, y_pred_BDTG, w, y_cat, y_pred_xgb, w, outName="MVA_Methods_ROC_EE.pdf")



if __name__ == "__main__":
    Classes = ["Merged-2Gsf", "DYJets", "QCD"]
    for tree in [
        # "dataset/TrainTree",
        "dataset/TestTree"
    ]:
        main(tree)