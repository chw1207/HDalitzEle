import sys, os
import time
from datetime import datetime
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
import numpy as np
import pandas as pd
import xgboost as xgb
import ROOT
import optuna
import json
import pickle as pkl
from contextlib import contextmanager
from sklearn.utils import shuffle
from sklearn.preprocessing import RobustScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error as MSE
from sklearn.metrics import median_absolute_error as MAE
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator, AutoMinorLocator
from colorPrint import *


def get_parser():
    parser = ArgumentParser(description="Script for merged electron energy regression")
    parser.add_argument("-r",  "--region", help="Which region to run. Options: [EB, EE]", default="", type=str)
    parser.add_argument("-e",  "--ext", help="external text pass to plot name and model name", default="", type=str)
    parser.add_argument("-t",  "--target", help="which target are used [reg or calib]", default="calib", type=str)
    parser.add_argument("--opt", help="perform parameter optimization", default=False, action="store_true")
    parser.add_argument("--robust", help="use RobustScaler or not", default=False, action="store_true")

    return parser


# setup the environment variables and mkdir to store the results       
def setup_env():
    global pkg_loc, file_loc, reg_loc, out_loc
    pkg_loc   = ROOT.gSystem.Getenv("HDalitzEle_LOC") 
    if (pkg_loc == ""):
        raise ValueError("HDalitzEle_LOC does not exist")
    reg_loc   = "{}/misc/EnergyReg".format(pkg_loc)
    file_loc  = "{}/python".format(reg_loc)
    interface = " #include \"{}/interface/EnRegression_Signal.h\" ".format(reg_loc)
    ROOT.gInterpreter.ProcessLine(interface) # effSigma
    
    out_loc = "{}/reg_results/XGBRegression{}_{}".format(reg_loc, args.ext, args.region)
    print("[INFO] Making the directory to save the results", flush=True)
    print("     - {}".format(out_loc), flush=True)
    os.makedirs(out_loc, exist_ok=True)


# function to cd to the working directory
# https://stackoverflow.com/a/24176022
@contextmanager
def cd(newdir):
    print("[INFO] Change to working directory", flush=True)
    print("     - {}".format(newdir), flush=True)
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
        

def data_preparing(signal_filename):
    
    global features    
    features = [
        "rho",
        "nVtx",
        "eleSCEta_Lead",
        "eleSCPhi_Lead",
        "eleSCRawEn_Lead",
        "eleCalibPt_Lead",

        "eledEtaAtVtx_Lead",
        "eledPhiAtVtx_Lead",
        "elePtError_Lead",
        "eleHoverE_Lead",
        "eleEoverP_Lead",
        "eleEoverPout_Lead",
        "eleEoverPInv_Lead",

        "eleSCEtaWidth_Lead",
        "eleSCPhiWidth_Lead",
        "eleSigmaIEtaIEtaFull5x5_Lead",
        "eleSigmaIPhiIPhiFull5x5_Lead",
        "eleR9Full5x5_Lead",
        "eleBrem_Lead",

        "gsfPtSum_Lead",
        "gsfPtRatio_Lead",
        "diTrkPt",
        "gsfDeltaR_Lead"
    ]
    if args.region == "EE":
        features.append("eleESEnToRawE_Lead")
    
    cut_base = "elePresel_Lead == 1 && category == 2 && eleCalibPt_Lead > 25"
    cut_region = "fabs(eleSCEta_Lead) < 1.479"
    if args.region == "EE":
        cut_region = "fabs(eleSCEta_Lead) >= 1.479"
    cut_target = "target < 1.2 && target > 0.8" 
    
    target_bra = "diGenEle.Pt()/eleCalibPt_Lead" if args.target == "calib" else "diGenEle.Pt()/elePt_Lead"
    print("[INFO] Reading the input file:", flush=True)
    print("     - {}".format(signal_filename), flush=True)
    print("[INFO] Training target is set as {}".format(target_bra), flush=True)
    print("[INFO] Applying the pre-selections:", flush=True)
    print("     - {} && {} && {}".format(cut_base, cut_region, cut_target), flush=True)
    with open("{}/features.pkl".format(out_loc), "wb") as f:
        print("[INFO] Save the training features in:", flush=True)
        print("     - {}/features.pkl".format(out_loc), flush=True)
        pkl.dump(features, f)
    
    ROOT.EnableImplicitMT(20) # comment it to ensure we have the same training and testing set each time
    data_sig = ROOT.RDataFrame("miniTree", signal_filename)\
                   .Define("target",            target_bra)\
                   .Define("diTrkPt",           "diTrk.Pt()")\
                   .Define("RatioRL",           "(eleEright_Lead - eleEleft_Lead)/(eleEright_Lead + eleEleft_Lead)")\
                   .Define("RatioBT",           "(eleEtop_Lead - eleEbottom_Lead)/(eleEtop_Lead + eleEbottom_Lead)")\
                   .Define("regwei",            "instwei * puwei")\
                   .Filter("{} && {} && {}".format(cut_base, cut_region, cut_target))\
                   .AsNumpy(columns=features+["regwei", "target"])
                   
    x = np.vstack([data_sig[var] for var in features]).T
    y = data_sig["target"]
    w = data_sig["regwei"]
    x, y, w = shuffle(x, y, w, random_state=123456)
    x_train_, x_test_, y_train_, y_test_, w_train_, w_test_ = train_test_split(x, y, w, test_size=0.2, random_state=42)
    
    if args.robust:
        print("[INFO] Applying the Robust scaling", flush=True)
        scaler = RobustScaler()
        x_train_ = scaler.fit_transform(x_train_)
        x_test_ = scaler.transform(x_test_)
    
        # extract the attributes from scaler and then save them to root file
        scaler_path_root = "{}/RobustScaler.root".format(out_loc)
        center_values = scaler.center_
        scale_values = scaler.scale_
        data = {"center": center_values, "scale": scale_values}
        df = ROOT.RDF.MakeNumpyDataFrame(data) # berfore root cersion 6.26, start from 6.28 using FromNumpy instead
        df.Snapshot("RobustScaler", scaler_path_root)
        print("[INFO] Save scaler in {}".format(scaler_path_root), flush=True)

        # save the scaler in normal way
        scaler_path_pkl = "{}/RobustScaler.pkl".format(out_loc)
        with open(scaler_path_pkl, "wb") as f:
            pkl.dump(scaler, f)
            print("[INFO] Save scaler in {}".format(scaler_path_pkl), flush=True)
  
    return x_train_, x_test_, y_train_, y_test_, w_train_, w_test_


# obj to be minimized by TPE algorithm
def objective(trial):
    params = {
        "max_depth":        trial.suggest_int("max_depth", 5, 18),
        "learning_rate":    trial.suggest_float("learning_rate", 0.01, 0.05),
        "subsample":        trial.suggest_float("subsample", 0.3, 0.99),
        "min_child_weight": trial.suggest_float("min_child_weight", 50, 300)
    }
    search_params = {**model_params, **params}
    
    model = xgb.XGBRegressor(**search_params)
    model.fit(
        x_train, y_train, 
        sample_weight=w_train, 
        eval_set=[(x_train, y_train), (x_test, y_test)], 
        sample_weight_eval_set=[w_train, w_test], 
        verbose=False, early_stopping_rounds=100
    )
    
    results = model.evals_result()
    return min(results["validation_1"]["mphe"])


# plotting style
def pltSty(ax, xName="x-axis", yName="y-axis", TitleSize=15, LabelSize=15, TickSize=13, yAuto=True):
    ax.set_xlabel(xName, fontsize=LabelSize, loc="right")
    ax.set_ylabel(yName, fontsize=LabelSize, loc="top")
    ax.text(1, 1, "(13 TeV)", horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, fontsize=TitleSize)
    ax.text(0, 1, "CMS", horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes, fontsize=TitleSize * 1.2, fontweight="bold")
    ax.text(TitleSize * 0.009 + 0.025, 1, "Simulation Preliminary", horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes, fontsize=TitleSize)

    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if (yAuto):
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction="in", length=7, labelsize=TickSize, top=True, right=True)
    ax.tick_params(direction="in", length=4, which="minor", labelsize=TickSize, top=True, right=True)
    
    
# https://github.com/scikit-learn/scikit-learn/blob/8c9c1f27b/sklearn/metrics/_plot/regression.py#L10     
# PredictionErrorDisplay          
def plot_deviation(y_true, y_pred, outName):
    max_value = max(np.max(y_true), np.max(y_pred))
    min_value = min(np.min(y_true), np.min(y_pred))
    plt.plot([min_value, max_value], [0, 0], linestyle="--", alpha=0.7, color="black")
    
    my_cmap = plt.cm.jet
    my_cmap.set_under("w", 1)
    plt.hist2d(y_true, y_pred - y_true, bins=(50, 50), range=[[0.8, 1.2],[-0.5, 0.5]], cmap=my_cmap, vmin=1)
    plt.colorbar() 
    pltSty(plt.gca(), xName="$\mathrm{target_{True}}$", yName="$\mathrm{target_{MVA} - target_{True}}$", TitleSize=11, LabelSize=11, TickSize=11)
    plt.savefig(outName, bbox_inches="tight")
    print("[INFO] Save the plot in {}".format(outName), flush=True)
    plt.close("all")
    

# relation between target_EGM target_MVA
def plot_relation(y_true, y_pred, outName):
    max_value = max(np.max(y_true), np.max(y_pred))
    min_value = min(np.min(y_true), np.min(y_pred))
    plt.plot([min_value, max_value], [min_value, max_value], linestyle="--", alpha=0.7, color="black")
    
    plt.hist2d(y_true, y_pred, bins=(50, 50), range=[[0.95, 1.05],[0.95, 1.05]], cmap=plt.cm.jet)
    plt.colorbar()
    pltSty(plt.gca(), xName="$\mathrm{target_{True}}$", yName="$\mathrm{target_{MVA}}$", TitleSize=11, LabelSize=11, TickSize=11)
    plt.savefig(outName, bbox_inches="tight")
    print("[INFO] Save the plot in {}".format(outName), flush=True)
    plt.close("all")
    

# target distribution
def plot_fvalue(y_test_, w_test_, w_train_, y_pred_test_, y_pred_train_, outName):
    h1 = ROOT.TH1D("h_egm", "", 80, 0.8, 1.2)
    h1.FillN(len(y_test_), np.array(y_test_, dtype=np.float64), np.array(w_test_, dtype=np.float64))
    effsig1 = ROOT.effSigma(h1)
    print("[INFO] Effective sigma of target", flush=True)
    print("     - EGM calibration: {} GeV".format(effsig1), flush=True)
    
    h2 = ROOT.TH1D("h_pred_test", "", 80, 0.8, 1.2)
    h2.FillN(len(y_pred_test_), np.array(y_pred_test_, dtype=np.float64), np.array(w_test_, dtype=np.float64))
    effsig2 = ROOT.effSigma(h2)
    print("     - XGB regression test: {} GeV".format(effsig2), flush=True)
    
    h3 = ROOT.TH1D("h_pred_train", "", 80, 0.8, 1.2)
    h3.FillN(len(y_pred_train_), np.array(y_pred_train_, dtype=np.float64), np.array(w_train_, dtype=np.float64))
    effsig3 = ROOT.effSigma(h3)
    print("     - XGB regression train: {} GeV".format(effsig3), flush=True)
    
    plt.hist(
        y_test_, bins=np.linspace(0.8, 1.2, num=80), 
        label="EGM($\mathrm{\sigma_{eff} = %.3f GeV}$)"%(effsig1), 
        histtype="step", linewidth=2, color="#202020", density=True
    )
    plt.hist(
        y_pred_train_, bins=np.linspace(0.8, 1.2, num=80), 
        label="XGB-train($\mathrm{\sigma_{eff} = %.3f GeV}$)"%(effsig3), 
        histtype="stepfilled", linewidth=2, color="#E16262", alpha=0.5, density=True
    )
    n, bins, patches = plt.hist(
        y_pred_test_, bins=np.linspace(0.8, 1.2, num=80), 
        label="XGB-test($\mathrm{\sigma_{eff} = %.3f GeV}$)"%(effsig2), 
        histtype="step", linewidth=2, color="#E16262", density=True
    )
    plt.legend(loc="best", fontsize=13, edgecolor="none")
    plt.ylim([0, max(n)*1.5])
    pltSty(plt.gca(), xName="$\mathrm{target}$", yName="$arb. unit$", TitleSize=11, LabelSize=11, TickSize=11)
    plt.savefig(outName, bbox_inches="tight")
    print("[INFO] Save the plot in {}".format(outName), flush=True)
    plt.close("all")


# feature importance
def plot_importance(bst, outName):
    bst.feature_names = [f.replace("_Lead", "") for f in features]
    importance = bst.get_score(importance_type="gain")
    max_importance = max(list(importance.values()))
    fig, ax = plt.subplots(1, 1, figsize=(6, 7))
    xgb.plot_importance(
        bst,
        ax=ax, importance_type="gain",
        height=0.7, grid=False, title=None, show_values=False,
        xlabel=None, ylabel=None, xlim=(0, max_importance * 1.2), color="#5893D4"
    )
    pltSty(ax, xName="Average Gain", yName="Features", yAuto=False)
    fig.savefig(outName, bbox_inches="tight")
    print("[INFO] Save plot in {}".format(outName), flush=True)
    plt.close("all")
    
        
if __name__ == "__main__":   
    start_time = time.time()
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(color.RED+"Execution date and time = {}".format(dt_string)+color.END, flush=True)
    print(color.BLUE + "---Start to train the regression!---" + color.END, flush=True)
         
    parser = get_parser()
    args = parser.parse_args()
    setup_env()
    
    with cd(file_loc):
        os.system("cp regression.py {}".format(out_loc))
        x_train, x_test, y_train, y_test, w_train, w_test = data_preparing("{}/reg_signal.root".format(reg_loc))
        
        # setup the training parameters
        model_params = {
            "objective":    "reg:pseudohubererror",
            "eval_metric":  "mphe",
            "tree_method":  "gpu_hist",
            "random_state": 123456,
            "grow_policy":  "lossguide",
            "n_estimators": 4000,
            "missing":      -9999.
        }
        opt_params = {
            "learning_rate": 0.012,
            "max_depth": 17,
            "subsample": 0.8,
            "min_child_weight": 60
        }
        train_params = {**model_params, **opt_params}
        
        # optimization
        if args.opt:
            print("[INFO] Performing the hyperparameter optimization", flush=True)
            study = optuna.create_study(
                direction="minimize", study_name="regression", 
                sampler=optuna.samplers.TPESampler(seed=123456),
                pruner=optuna.pruners.MedianPruner()
            )
            study.optimize(objective, n_trials=100)
            
            print(color.GREEN+"Best parameters: {}".format(study.best_params)+color.END, flush=True)
            train_params = {**model_params, **study.best_params}
        
        with open("{}/XGB_params.json".format(out_loc), "w") as f:
            print("[INFO] Save the paramters in {}".format(f.name), flush=True)
            json.dump(train_params, f, indent=4)
            
        # start training
        # https://heartbeat.comet.ml/5-regression-loss-functions-all-machine-learners-should-know-4fb140e9d4b0
        reg = xgb.XGBRegressor(**train_params)
        reg.fit(
            x_train, y_train, sample_weight=w_train, 
            eval_set=[(x_train, y_train), (x_test, y_test)], 
            sample_weight_eval_set=[w_train, w_test],
            early_stopping_rounds=100
        )
        
        bst = reg.get_booster()
        bst.save_model("{}/XGB_Regression.txt".format(out_loc))
        print("[INFO] Save the model in {}".format("{}/XGB_Regression.txt".format(out_loc)), flush=True)
        
        # Predict the model
        y_pred_train = reg.predict(x_train)
        y_pred_test  = reg.predict(x_test)

        # MSE Computation
        print("[INFO] Statistical results:")
        mse_train = MSE(y_train, y_pred_train, sample_weight=w_train)
        print("     - Training MSE : %f" %(mse_train), flush=True)
        mse_test = MSE(y_test, y_pred_test, sample_weight=w_test)
        print("     - Testing  MSE : %f" %(mse_test), flush=True)
        mae_train = MAE(y_train, y_pred_train, sample_weight=w_train)
        print("     - Training MAE : %f" %(mae_train), flush=True)
        mae_test = MAE(y_test, y_pred_test, sample_weight=w_test)
        print("     - Testing  MAE : %f" %(mae_test), flush=True)

        val_results = reg.evals_result()
        results = {
            "Training MPHE": min(val_results["validation_0"]["mphe"]),
            "Testing MPHE": min(val_results["validation_1"]["mphe"]), 
            "Training MSE": mse_train,
            "Testing MSE": mse_test, 
            "Training MAE": mae_train,
            "Testing MAE": mae_test, 
        }
        with open("{}/XGB_results.json".format(out_loc), "w") as f:
            print("[INFO] Save the statistical results in {}".format(f.name), flush=True)
            json.dump(results, f, indent=4)
            
        # visualize the result
        plot_importance(bst, "{}/Importancee_XGB.pdf".format(out_loc))
        plot_relation(y_test, y_pred_test, "{}/relation_XGB.pdf".format(out_loc))
        plot_deviation(y_test, y_pred_test, "{}/deviation_XGB.pdf".format(out_loc))
        plot_fvalue(y_test, w_test, w_train, y_pred_test, y_pred_train, "{}/fvalue_XGB.pdf".format(out_loc))
        
    print(color.BLUE + "---All done!---" + color.END)
    seconds = time.time() - start_time
    print("Total Time Taken: {}".format(time.strftime("%H:%M:%S",time.gmtime(seconds))))
    
  
    
    
    
