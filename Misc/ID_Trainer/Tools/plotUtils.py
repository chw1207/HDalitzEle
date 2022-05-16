import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import ROOT
from Tools.CMS_lumi import CMS_lumi
from matplotlib.ticker import AutoLocator, AutoMinorLocator
from xgboost import plot_importance
from sklearn.metrics import roc_curve, auc


class color:
    PURPLE = "\033[95m"
    CYAN = "\033[96m"
    DARKCYAN = "\033[36m"
    BLUE = "\033[94m"
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"
    END = "\033[0m"


def pltSty(ax, xName = "x-axis", yName = "y-axis", TitleSize = 15, LabelSize = 15, TickSize = 13, MajTickLength = 7, MinTickLength = 4, yAuto = True):
    ax.set_xlabel(xName, fontsize = LabelSize, loc = "right")
    ax.set_ylabel(yName, fontsize = LabelSize, loc = "top")
    ax.text(1, 1, "(13 TeV)", horizontalalignment = "right", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize)
    ax.text(0, 1, "CMS", horizontalalignment = "left", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize * 1.3, fontweight = "bold")
    ax.text(TitleSize * 0.009 + 0.02, 1, "Work-in-progress", horizontalalignment = "left", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize)

    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if (yAuto):
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction = "in", length = MajTickLength, labelsize = TickSize, top = True, right = True)
    ax.tick_params(direction = "in", length = MinTickLength, which = "minor", labelsize = TickSize, top = True, right = True)


# Reference: https://stackoverflow.com/a/29643643
def hex2rgb(h, a):
    h = h.lstrip('#')
    rgb = tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))
    return (rgb[0], rgb[1], rgb[2], a)


def plot_featureHist(df, feature, num_class, Class, bininfo, XaxisName, histcolor=[""], logy=False, y_axisscale=1.3, outname="testfeatureplot", wei="", category="EleType"):

    plt.figure(figsize=(6, 6))
    binContentMax_list, binContentMin_list = [], []
    for i in range(num_class):
        binContent, binEdget, patches = plt.hist(
            df[df[category] == i][feature],
            range = (bininfo[1], bininfo[2]) ,
            bins = bininfo[0],
            ls = "-", lw = 2, fc = hex2rgb(histcolor[i], 0), ec = histcolor[i],
            histtype = "stepfilled",
            weights = df[df[category] == i][wei],
            density = True,
            label = Class[i]
        )
        # print(binContent, Class[i])
        binContentMax_list.append(np.max(binContent))
        binContentMin_list.append(np.min(binContent[binContent != 0]))

    ymaxval = np.max(binContentMax_list)
    yminval = np.min(binContentMin_list) if np.min(binContentMin_list) < 1E5 else 1E5

    ax = plt.gca()
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(yminval * 0.9, ymaxval * y_axisscale)
    else:
        ax.set_ylim(0, ymaxval * y_axisscale)
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(bininfo[1], bininfo[2])
    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.set_ylabel("A.U.", fontsize=15, loc="top")
    ax.set_xlabel(XaxisName, fontsize=15, loc="right")

    # lumi information text
    ax.text(0, 1, "CMS", horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes, fontsize=15, fontweight="bold")
    ax.text(1, 1, "$(13TeV,\ 2017)$", horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, fontsize=13)

    ax.legend(loc="best", edgecolor= "none", fontsize=13)
    ax.tick_params( # major tick
        direction="in", length=8,
        labelsize=13, top=True, right=True
    )
    ax.tick_params( # minor tick
        direction="in", length=4, which="minor",
        labelsize=13, top=True, right=True
    )

    ax.legend(loc="best", edgecolor="none", fontsize=13)
    plt.tight_layout()
    plt.draw()
    plt.savefig("{}.pdf".format(outname))
    print("[INFO] Save fig in {}.pdf".format(outname), flush=True)
    plt.close("all")


def plot_ptEtaRwt(df, ptName, etaName, ptBin, etaBin, outdir):
    fig, ax = plt.subplots(1, 2, figsize = (10, 5))
    for i,group_df in df[df["Dataset"] == "Train"].groupby("Class"):
        group_df[ptName].hist(histtype="step", bins=ptBin, alpha=0.7,label=i, ax=ax[0], density=False, ls="-", weights=group_df["instwei"], linewidth=2)
        ax[0].set_title("$p_T$ before reweighting")
        ax[0].legend()
        ax[0].set_xscale("log")
        group_df[ptName].hist(histtype="step", bins=ptBin, alpha=0.7,label=i, ax=ax[1], density=False, ls="-", weights=group_df["NewWt"], linewidth=2)
        ax[1].set_title("$p_T$ after reweighting")
        ax[1].legend()
        ax[1].set_xscale("log")
    fig.savefig(outdir+"/pT_rwt.pdf")
    print("[INFO] Save fig in {}".format(outdir+"/pT_rwt.pdf"), flush=True)

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    for i,group_df in df[df["Dataset"] == "Train"].groupby("Class"):
        group_df[etaName].hist(
            histtype="step",
            bins=etaBin,
            alpha=0.7, label=i, ax=ax[0], density=False, ls="-", weights=group_df["instwei"], linewidth=2
        )
        ax[0].set_title("$\eta$ before reweighting")
        ax[0].legend()
        group_df[etaName].hist(
            histtype="step",
            bins=etaBin,
            alpha=0.7, label=i, ax=ax[1], density=False, ls="-", weights=group_df["NewWt"], linewidth=2
        )
        ax[1].set_title("$\eta$ after reweighting")
        ax[1].legend()
    fig.savefig(outdir+"/eta_rwt.pdf")
    print("[INFO] Save fig in {}".format(outdir+"/eta_rwt.pdf"), flush=True)
    plt.close("all")


# for the hyperparameter tuning
def main_plot_history_mod(trials, do_show=True, status_colors=None, title="", outName="test_trials.png"):
    # self is an Experiment
    if status_colors is None:
        status_colors = "b"

    best_err = trials.average_best_error()
    print("avg best error:", round(best_err, 4))
    plt.axhline(best_err, c="g", zorder = 1, label = "min. loss = {}".format(round(best_err, 4)))

    Ys = trials.losses()
    plt.scatter(range(len(Ys)), Ys, c=status_colors, zorder = 2, label = "loss")

    plt.title(title)
    if do_show:
        plt.show()

    ax = plt.gca()

    up, do = max(Ys), min(Ys)
    ax.set_ylim([do*0.97, up*1.09])
    pltSty(ax, xName = "Trials", yName = "Log Loss")
    plt.legend(loc = "best", facecolor = "none", fontsize = 13)
    plt.tight_layout()

    print("[INFO] Save fig in {}".format(outName), flush=True)
    plt.savefig(outName)
    plt.close("all")


def plot_importance_mod(bst, features=[""], outName="test_importance.pdf"):
    bst["booster"].feature_names = features
    importance = bst["booster"].get_score(importance_type = "gain")
    for key in importance.keys():
        importance[key] = round(importance[key], 2)
    max_importance = max(list(importance.values()))

    fig, ax = plt.subplots(1, 1, figsize = (6, 7))
    plot_importance(
        importance,
        ax = ax, importance_type = "gain",
        height = 0.7, grid = False, title = None,
        xlabel = None, ylabel = None, xlim = (0, max_importance * 1.2), color = "#5893D4"
    )
    pltSty(ax, xName = "Importance", yName = "Features", yAuto = False)
    fig.savefig(outName, bbox_inches = "tight")
    print("[INFO] Save fig in {}".format(outName), flush=True)
    plt.close("all")


def plot_loss_history(bst, metric="mlogloss", outName="test_importance.pdf"):
    eval_metric = metric
    fig, axes = plt.subplots(1, 1, figsize = (6, 6))
    axes.plot(range(1, len(bst["history"]["train"][eval_metric]) + 1), bst["history"]["train"][eval_metric], linewidth=3, linestyle="-", color="#184d47", label="Training loss")
    axes.plot(range(1, len(bst["history"]["test"][eval_metric]) + 1), bst["history"]["test"][eval_metric], linewidth=3, linestyle="--", color="#e2703a", label="Testing loss")
    axes.legend(loc="best", edgecolor="none", fontsize=15, title_fontsize=15)
    pltSty(axes, xName="Epochs", yName="Log loss")
    fig.savefig(outName, bbox_inches="tight")
    print("[INFO] Save fig in %s" %(outName), flush=True)
    plt.close("all")


# copy from the TensorFlow to_categoriacal function
# Reference: https://github.com/keras-team/keras/blob/v2.7.0/keras/utils/np_utils.py#L21-L74
def ToCategorical(y, num_classes = None, dtype = "float64"):
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


def plot_MVAscore(
    df, Y_train, Y_test, y_train_pred, y_test_pred,
    wt="instwei", Classes=[""], ClassColors=[""], isb=False, MVA="XGB", outName="test_score.pdf"
):
    if isb == True:
        df[(df["TrainDataset"] == 1)][MVA+"_pred"] = y_train_pred
        df[(df["TrainDataset"] == 0)][MVA+"_pred"] = y_test_pred

        figMVA, axesMVA = plt.subplots(1, 1, figsize = (6, 6))
        for i, j in enumerate(Classes):
            axesMVA.hist(
                df[(df["TrainDataset"] == 0) & (df.Class == j)][MVA+"_pred"],
                bins=np.linspace(0, 1, num=40), label=Classes[i]+": Test",
                weights=df[(df["TrainDataset"] == 0) & (df.Class == j)][wt]/np.sum(df[(df["TrainDataset"] == 0) & (df.Class == j)][wt]),
                histtype="step", linewidth=2, color=ClassColors[i]
            )
            axesMVA.hist(
                df[(df["TrainDataset"] == 1) & (df.Class == j)][MVA+"_pred"],
                bins=np.linspace(0, 1, num=40), label=Classes[i]+": Train",
                weights=df[(df["TrainDataset"] == 1) & (df.Class == j)][wt]/np.sum(df[(df["TrainDataset"] == 1) & (df.Class == j)][wt]),
                histtype="stepfilled", alpha=0.5, linewidth=2, color=ClassColors[i]
            )
        pltSty(axesMVA, xName="{}_pred".format(MVA), yName="Normalized entries", TitleSize=16, LabelSize=15, TickSize=14, MajTickLength=10, MinTickLength=6)
        axesMVA.legend(loc="upper center", ncol=2, fontsize=13, edgecolor="none")
        axesMVA.set_yscale("log")
        axesMVA.set_ylim([1E-4, 10])
        axesMVA.set_xlim([-0.05, 1.05])

    else:
        Y_train_categorical = ToCategorical(Y_train, num_classes=len(Classes))
        Y_test_categorical  = ToCategorical(Y_test, num_classes=len(Classes))
        plotwt_train = np.array(df[(df["TrainDataset"] == 1)][wt])
        plotwt_test = np.array(df[(df["TrainDataset"] == 0)][wt])

        figMVA, axesMVA = plt.subplots(1, len(Classes), figsize=(len(Classes)*8.5, 7))
        for i in range(len(Classes)):
            axMVA = axesMVA[i]
            for k in range(len(Classes)):
                axMVA.hist(
                    y_test_pred[:, i][Y_test_categorical[:, k] == 1],
                    bins=np.linspace(0, 1, num=40), label=Classes[k]+": Test",
                    weights=plotwt_test[Y_test_categorical[:, k] == 1]/np.sum(plotwt_test[Y_test_categorical[:, k] == 1]),
                    histtype="step", linewidth=2, color=ClassColors[k]
                )
                axMVA.hist(
                    y_train_pred[:, i][Y_train_categorical[:, k]==1],
                    bins=np.linspace(0, 1, num=40), label=Classes[k]+": Train",
                    weights=plotwt_train[Y_train_categorical[:, k] == 1]/np.sum(plotwt_train[Y_train_categorical[:, k] == 1]),
                    histtype="stepfilled", alpha=0.5, linewidth=2, color=ClassColors[k]
                )

            pltSty(axMVA, xName="{}_pred".format(MVA), yName="Normalized entries", TitleSize=16, LabelSize=17, TickSize=17, MajTickLength=10, MinTickLength=6)
            axMVA.legend(title="{} vs Rest".format(Classes[i]), loc="upper center", fontsize=15, ncol=2, edgecolor="none", title_fontsize=17)
            axMVA.set_yscale("log")
            axMVA.set_ylim([1E-4, 100])
            axMVA.set_xlim([-0.05, 1.05])

    figMVA.savefig(outName, bbox_inches="tight")
    print("[INFO] Save fig in {}".format(outName), flush=True)
    plt.close("all")


def plot_ROC(
    df, Y_train, Y_test, y_train_pred, y_test_pred, cat="EleType",
    wt="instwei", Classes=[""], isb=False, MVA="XGB", outName="test_roc.pdf"
):
    plotwt_train = np.array(df[(df["TrainDataset"] == 1)][wt])
    plotwt_test = np.array(df[(df["TrainDataset"] == 0)][wt])
    WPCuts = 0.5 # useless to multi-classification
    if isb == True:
        df[(df["TrainDataset"] == 1)][MVA+"_pred"] = y_train_pred
        df[(df["TrainDataset"] == 0)][MVA+"_pred"] = y_test_pred

        fig, axes = plt.subplots(1, 1, figsize=(6, 6))

        # tpr: true positive rate = signal efficiency
        # fpr: false positive rate = background efficiency
        fpr, tpr, th = roc_curve(
            df[(df["TrainDataset"] == 0)][cat], y_test_pred,
            sample_weight=plotwt_test
        )
        fpr_tr, tpr_tr, th_tr = roc_curve(
            df[(df["TrainDataset"] == 1)][cat], y_train_pred,
            sample_weight=plotwt_train
        )

        bkgrej = (1 - fpr)
        bkgrej_tr = (1 - fpr_tr)

        gmeans = np.sqrt(tpr_tr * bkgrej_tr)
        ix = np.argmax(gmeans)
        WPCuts = th_tr[ix]
        print("Best Threshold = %.3f, G-Mean = %.3f" %(WPCuts, gmeans[ix]), flush=True)

        roc_auc = auc(fpr, tpr)
        roc_auc_tr = auc(fpr_tr, tpr_tr)

        axes.plot(tpr_tr, bkgrej_tr, label="XGB Training auc = %.1f%s" %(roc_auc_tr*100, "%"), linewidth=3, linestyle="-", color="#064635", zorder=1)
        axes.plot(tpr, bkgrej, label="XGB Testing auc = %.1f%s"% (roc_auc*100, "%"), linewidth=3, linestyle="--", color="#e2703a", zorder=2)
        axes.scatter(tpr_tr[ix], bkgrej_tr[ix], marker="o", color="#202020", label="Threshold with G-Mean = %.3f" %(gmeans[ix]), s=90, zorder=3)

        pltSty(axes, xName="Signal efficiency", yName="Background rejection", TitleSize=15, LabelSize=15, TickSize=15, MajTickLength=7, MinTickLength=4)
        axes.legend(loc="best", edgecolor="none", facecolor="none", fontsize=13)

    else:
        Y_train_categorical = ToCategorical(Y_train, num_classes=len(Classes))
        Y_test_categorical  = ToCategorical(Y_test, num_classes=len(Classes))

        fig, axes = plt.subplots(1, len(Classes), figsize=(len(Classes)*8.5, 7))
        for i in range(len(Classes)):
            ax = axes[i]

            # tpr: true positive rate = signal efficiency
            # fpr: false positive rate = background efficiency
            fpr, tpr, th = roc_curve(
                Y_test_categorical[:, i], y_test_pred[:, i],
                sample_weight=plotwt_test
            )
            fpr_tr, tpr_tr, th_tr = roc_curve(
                Y_train_categorical[:, i], y_train_pred[:, i],
                sample_weight=plotwt_train
            )

            bkgrej = (1 - fpr)
            bkgrej_tr = (1 - fpr_tr)

            roc_auc = auc(fpr, tpr)
            roc_auc_tr = auc(fpr_tr, tpr_tr)

            ax.plot(tpr_tr, bkgrej_tr, label="XGB Training auc = %.1f%s" %(roc_auc_tr*100, "%"),linewidth=3, linestyle="-", color="#064635", zorder=1)
            ax.plot(tpr, bkgrej, label="XGB Testing auc = %.1f%s"% (roc_auc*100, "%"), linewidth=3, linestyle="--", color="#e2703a", zorder=2)

            pltSty(ax, xName="Signal efficiency", yName="Background rejection", TitleSize=16, LabelSize=17, TickSize=17, MajTickLength=10, MinTickLength=6)
            ax.legend(title="{} vs Rest".format(Classes[i]), loc="best", edgecolor="none", facecolor="none", fontsize=17, title_fontsize=17)

    fig.savefig(outName, bbox_inches="tight")
    print("[INFO] Save fig in {}".format(outName), flush=True)
    plt.close("all")

    return plotwt_train, plotwt_test, WPCuts



## ----- for multi-classification -----##
# Reference:
# https://github.com/javaidnabi31/Multi-class-with-imbalanced-dataset-classification/blob/master/20-news-group-classification.ipynb
import itertools
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


def fillHist1D(hist, val, wei):
    for i in range(len(val)):
        hist.Fill(val[i], wei[i])
    return hist


def plot_eff(df, bin, trueLable="EleType", predLable="ElePredType", var="eleCalibPt", xaxisName="p^{e}_{T} [GeV]", wt="instwei", isb=False, Classes=[""], ClassColors=[""], outName="test_eff.png"):
    eff_bin = bin[:-1] if "Pt" in var else bin
    sigNum = 0 if isb == False else 1
    if isb == False:
        print("Assume the first class is the signal!")
    hdenDict, hnumDict, efferrDict = {}, {}, {}
    for i in range(len(Classes)):
        arr_den = df[df[trueLable] == i][var].to_numpy()
        arr_num = df[(df[trueLable] == i) & (df[predLable] == sigNum)][var].to_numpy()
        arr_den_wei = df[df[trueLable] == i][wt].to_numpy()
        arr_num_wei = df[(df[trueLable] == i) & (df[predLable] == sigNum)][wt].to_numpy()

        hdenDict[Classes[i]] = ROOT.TH1F("hden{}".format(i), " ", len(eff_bin)-1, np.asarray(eff_bin, "d"))
        fillHist1D(hdenDict[Classes[i]], arr_den, arr_den_wei)

        hnumDict[Classes[i]] = ROOT.TH1F("hnum{}".format(i), " ", len(eff_bin)-1, np.asarray(eff_bin, "d"))
        fillHist1D(hnumDict[Classes[i]], arr_num, arr_num_wei)

        efferrDict[Classes[i]] = ROOT.TGraphAsymmErrors()
        efferrDict[Classes[i]].BayesDivide(hnumDict[Classes[i]], hdenDict[Classes[i]])

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)

    c = ROOT.TCanvas("c", "", 800, 800)
    c.cd()
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.07)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetBottomMargin(0.13)
    TickSize = 0.02
    AxisTitleSize = 0.05
    AxisLabelSize = 0.045

    efferrDict[Classes[0]].GetXaxis().SetTitle(xaxisName)
    efferrDict[Classes[0]].GetYaxis().SetTitle("Efficiency")
    efferrDict[Classes[0]].GetXaxis().SetRangeUser(eff_bin[0], eff_bin[-1])
    efferrDict[Classes[0]].GetYaxis().SetRangeUser(0, 1.2)
    efferrDict[Classes[0]].GetXaxis().SetTickSize(TickSize)
    efferrDict[Classes[0]].GetXaxis().SetTitleSize(AxisTitleSize)
    efferrDict[Classes[0]].GetXaxis().SetLabelSize(AxisLabelSize)
    efferrDict[Classes[0]].GetYaxis().SetTickSize(TickSize)
    efferrDict[Classes[0]].GetYaxis().SetTitleSize(AxisTitleSize)
    efferrDict[Classes[0]].GetYaxis().SetLabelSize(AxisLabelSize)
    efferrDict[Classes[0]].GetXaxis().SetTitleOffset(1.1)
    efferrDict[Classes[0]].GetYaxis().SetTitleOffset(1.2)

    efferrDict[Classes[0]].SetMarkerColor(ROOT.TColor.GetColor(ClassColors[0]))
    efferrDict[Classes[0]].SetMarkerSize(1.5)
    efferrDict[Classes[0]].SetMarkerStyle(20)
    efferrDict[Classes[0]].SetLineColor(ROOT.TColor.GetColor(ClassColors[0]))
    efferrDict[Classes[0]].SetLineWidth(2)
    efferrDict[Classes[0]].Draw("AP")

    for j in range(1, len(Classes)):
        efferrDict[Classes[j]].SetMarkerColor(ROOT.TColor.GetColor(ClassColors[j]))
        efferrDict[Classes[j]].SetMarkerSize(1.5)
        efferrDict[Classes[j]].SetMarkerStyle(20)
        efferrDict[Classes[j]].SetLineColor(ROOT.TColor.GetColor(ClassColors[j]))
        efferrDict[Classes[j]].SetLineWidth(2)
        efferrDict[Classes[j]].Draw("P same")

    l = ROOT.TLegend(0.17, 0.82, 0.25*len(Classes), 0.88)
    l.SetTextSize(0.04)
    l.SetNColumns(len(Classes))
    for key, value in efferrDict.items():
        l.AddEntry(value, key, "lep")
    l.SetFillColor(0)
    l.SetLineColorAlpha(0, 0)
    l.SetFillStyle(0)
    l.Draw()

    CMS_lumi(c, 4, 0, "", 2017, True, "Work-in-progress", "", "")

    c.SaveAs(outName)
    c.Close()
