#===============================================#
#              Merged ID features               #
#              Merged ID models                 #
#===============================================#
features_m1 = [
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

    # "gsfPtRatio",
    # "gsfDeltaR",
    "gsfRelPtRatio"
]

features_m2 = [
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
    "gsfRelPtRatio"
]

features = {"M1": features_m1, "M2": features_m2}
models = {
    "M1EB": "../ID_Trainer/Results/Output_Merged1GsfID_dask_Fall17_EB_CheckALL/XGB/XGB_modelXGB.pkl",
    "M2EB": "../ID_Trainer/Results/Output_Merged2GsfID_dask_Fall17_EB_CheckALL/XGB/XGB_modelXGB.pkl",
    "M1EE": "../ID_Trainer/Results/Output_Merged1GsfID_dask_Fall17_EE_CheckALL/XGB/XGB_modelXGB.pkl",
    "M2EE": "../ID_Trainer/Results/Output_Merged2GsfID_dask_Fall17_EE_CheckALL/XGB/XGB_modelXGB.pkl"
}
models_b = {
    "M1EB": "../ID_Trainer/Results/Output_Merged1GsfID_dask_Fall17_EB_CheckALL_Binary/XGB/XGB_modelXGB.pkl",
    "M2EB": "../ID_Trainer/Results/Output_Merged2GsfID_dask_Fall17_EB_CheckALL_Binary/XGB/XGB_modelXGB.pkl",
    "M1EE": "../ID_Trainer/Results/Output_Merged1GsfID_dask_Fall17_EE_CheckALL_Binary/XGB/XGB_modelXGB.pkl",
    "M2EE": "../ID_Trainer/Results/Output_Merged2GsfID_dask_Fall17_EE_CheckALL_Binary/XGB/XGB_modelXGB.pkl"
}

#===============================================#
#                MC samples list                #
#            Signal: ggF, VBF, WH, ZH           #
#            background: DYJets, GJets, QCD     #
#===============================================#
# [unit]: luminosity(fb), xs(invfb)
# Luminosity information can be found in: /data3/ggNtuples/V10_06_20_00/LumiInfo
# xs information:
#*      1) Zg: https://cms-gen-dev.cern.ch/xsdb/?searchQuery=DAS=ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8
#*      2) TT: https://cms-gen-dev.cern.ch/xsdb/?searchQuery=DAS=TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8

MCSample = {
    "ZGToLLG_2016_preVFP": {
        "production": ["ZGToLLG"],
        "lumi": [19.28],
        "xs": [51.1 * 1000.],
        "path": ["/data6/ggNtuples/V10_06_30_00/job_UL16_Zg_aMCatNLO_preVFP",],
        "outpath": ["/data4/chenghan/mc/V10_06_30_00/job_UL16_Zg_aMCatNLO_preVFP",]
    },
    "ZGToLLG_2016_postVFP": {
        "production": ["ZGToLLG"],
        "lumi": [16.64],
        "xs": [51.1 * 1000.],
        "path": ["/data6/ggNtuples/V10_06_30_00/job_UL16_Zg_aMCatNLO_postVFP",],
        "outpath": ["/data4/chenghan/mc/V10_06_30_00/job_UL16_Zg_aMCatNLO_postVFP",]
    },
    "ZGToLLG_2017": {
        "production": ["ZGToLLG"],
        "lumi": [41.48],
        "xs": [51.1 * 1000.],
        "path": ["/data6/ggNtuples/V10_06_30_00/job_UL17_Zg_aMCatNLO",],
        "outpath": ["/data4/chenghan/mc/V10_06_30_00/job_UL17_Zg_aMCatNLO",]
    },
    "ZGToLLG_2018": {
        "production": ["ZGToLLG"],
        "lumi": [59.82],
        "xs": [51.1 * 1000.],
        "path": ["/data6/ggNtuples/V10_06_30_00/job_UL18_Zg_aMCatNLO",],
        "outpath": ["/data4/chenghan/mc/V10_06_30_00/job_UL18_Zg_aMCatNLO",]
    },
    "TTJets_2017": {
        "production": ["TTJets"],
        "lumi": [41.48],
        "xs": [831.76 * 1000.],
        "path": ["/data6/ggNtuples/V10_02_10_08/job_UL17_TT_aMCatNLO",],
        "outpath": ["/data4/chenghan/mc/V10_02_10_08/job_UL17_TT_aMCatNLO",]
    },


    # Fall17 samples
    # https://cms-gen-dev.cern.ch/xsdb/?searchQuery=DAS=DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8
    "DYJets_2017": {
        "production": ["DYJets"],
        "lumi": [41.525],
        "xs": [6529 * 1000.],
        "path": [
            "/data6/ggNtuples/V10_02_10_08/job_fall17_newPU_DYJetsToLL_m50_aMCatNLO*/"
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_08/job_fall17_newPU_DYJetsToLL_m50_aMCatNLO/"
        ]
    },
}

#===============================================#
#              Data samples list                #
#  Run: 2016BCDEF1(pre_VFP), 2016F2GH(post_VFP) #
#          Run: 2017BCDEF, 2018ABCD             #
#===============================================#
DataSample = {
    "Data_2016_preVFP": {
        "run": ["2016B", "2016C", "2016D", "2016E", "2016F1"],
        "lumi": [4.671835672, 2.621295973, 4.285851496, 4.065974560, 2.717344923],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2016B_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2016C_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2016D_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2016E_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2016F1_UL/"
        ],
        "outpath": [
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2016B_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2016C_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2016D_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2016E_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2016F1_UL/"
        ]
    },
    "Data_2016_postVFP": {
        "run": ["2016F2", "2016G", "2016H"],
        "lumi": [0.418120616, 7.652808345, 8.739883563],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2016F2_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2016G_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2016H_UL/"
        ],
        "outpath": [
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2016F2_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2016G_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2016H_UL/"
        ]
    },
    "Data_2017": {
        "run": ["2017B", "2017C", "2017D", "2017E", "2017F"],
        "lumi": [4.803363110, 9.572498264, 4.247682053, 9.313642402, 13.539221759],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2017B_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2017C_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2017D_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2017E_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2017F_UL/"
        ],
        "outpath": [
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2017B_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2017C_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2017D_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2017E_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2017F_UL/"
        ]
    },
    "Data_2018": {
        "run": ["2018A", "2018B", "2018C", "2018D"],
        "lumi": [14.026724850, 7.060617355, 6.894770971, 31.834115905],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2018A_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2018B_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2018C_UL/",
            "/data3/ggNtuples/V10_06_30_00/job_DoubleMu_Run2018D_UL/"
        ],
        "outpath": [
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2018A_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2018B_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2018C_UL/",
            "/data4/chenghan/data/V10_06_30_00/job_DoubleMu_Run2018D_UL/"
        ]
    }
}
