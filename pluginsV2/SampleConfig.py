
#===============================================#
#              Merged ID features               #
#              Merged ID models                 #
#===============================================#
features_m1 = [
    # With GJets ID
    # "rho",
    # "eleSCEta",
    # "eleSCRawEn",

    # "eleD0",
    # "eleDz",
    # "eleSIP",
    # "elePtError",
    # "eleSCPhiWidth",
    # "eleEoverP",
    # "eleEoverPout",
    # "eleEoverPInv",
    # "eledEtaAtVtx",
    # "eledPhiAtVtx",
    # "eleSigmaIEtaIEtaFull5x5",
    # "elePFChIso",
    # "elePFPhoIso",
    # "eleR9Full5x5",
    # "eleTrkdxy",
    # "eleMissHits",
    # "gsfRelPtRatio"


    # Without GJets ID
    "rho",
    "eleSCEta",
    "eleSCRawEn",
    
    "eleD0",
    "eleDz",
    "eleSIP",
    "elePtError",
    "eleSCEtaWidth",
    "eleHoverE",
    "eleEoverP",
    "eleEoverPInv",
    "eledEtaAtVtx",
    "eledPhiAtVtx",
    "eleSigmaIEtaIEtaFull5x5",
    "elePFChIso",
    "elePFPhoIso",
    "elePFNeuIso",
    "eleR9Full5x5",
    "eleTrkdxy",
    "gsfRelPtRatio"
]

features_m2 = [
    # With GJets ID
    # "rho",
    # "eleSCEta",
    # "eleSCRawEn",

    # "eleD0",
    # "eleDz",
    # "elePtError",
    # "eleHoverE",
    # "eleEoverP",
    # "eleEoverPInv",
    # "eledPhiAtVtx",
    # "eleSigmaIEtaIEtaFull5x5",
    # "elePFChIso",
    # "elePFPhoIso",
    # "eleR9Full5x5",
    # "eleTrkdxy",
    # "eleMissHits",
    # "eleConvVeto",
    # "gsfPtRatio",
    # "gsfDeltaR",
    # "gsfRelPtRatio"


    # Without GJets ID
    "rho",
    "eleSCEta",
    "eleSCRawEn",

    "elePtError",
    "eleSCPhiWidth",
    "eleHoverE",
    "eleEoverP",
    "eleEoverPout",
    "eleEoverPInv",
    "eleBrem",
    "eledEtaAtVtx",
    "eledPhiAtVtx",
    "eleSigmaIEtaIEtaFull5x5",
    "elePFChIso",
    "elePFPhoIso",
    "elePFNeuIso",
    "eleR9Full5x5",
    "gsfPtRatio",
    "gsfDeltaR",
    "gsfRelPtRatio"
]

features = {"M1": features_m1, "M2": features_m2}
models = {
    # With GJets ID
    # "M1EB": "./Misc/ID_Trainer/Results/Output_Merged1GsfID_EB_2017_dask/XGB/XGB_modelXGB.pkl",
    # "M2EB": "./Misc/ID_Trainer/Results/Output_Merged2GsfID_EB_2017_dask/XGB/XGB_modelXGB.pkl",
    # "M1EE": "./Misc/ID_Trainer/Results/Output_Merged1GsfID_EE_2017_dask/XGB/XGB_modelXGB.pkl",
    # "M2EE": "./Misc/ID_Trainer/Results/Output_Merged2GsfID_EE_2017_dask/XGB/XGB_modelXGB.pkl"


    # Without GJets ID
    "M1EB": "./Misc/ID_Trainer/Results/Output_Merged1GsfID_EB_2017_daskNoGJets/XGB/XGB_modelXGB.pkl",
    "M2EB": "./Misc/ID_Trainer/Results/Output_Merged2GsfID_EB_2017_daskNoGJets/XGB/XGB_modelXGB.pkl",
    "M1EE": "./Misc/ID_Trainer/Results/Output_Merged1GsfID_EE_2017_daskNoGJets/XGB/XGB_modelXGB.pkl",
    "M2EE": "./Misc/ID_Trainer/Results/Output_Merged2GsfID_EE_2017_daskNoGJets/XGB/XGB_modelXGB.pkl"
}


#===============================================#
#                MC samples list                #
#            Signal: ggF, VBF, WH, ZH           #
#            background: DYJets, GJets, QCD     #
#===============================================#
# Luminosity informaation can be found in: /data3/ggNtuples/V10_06_20_00/LumiInfo
# [unit]: luminosity(invfb), xs(fb)

#! [FixedME]: xs, lumi and path need to modify once have UL signal sample!
MCSample = {

    "HDalitz_2017": {
        "production": [
            "ggF_m125", "VBF_m125", "WH_m125", "ZH_m125", 
            "ggF_m120", "VBF_m120", "WH_m120", "ZH_m120",
            "ggF_m130", "VBF_m130", "WH_m130", "ZH_m130"
        ],
        "lumi": [41.525],
        "xs": [
            48.58 * 1000 * 8.10E-5,
            3.782 * 1000 * 8.10E-5,
            1.373 * 1000 * 8.10E-5,
            0.8839 * 1000 * 8.10E-5,

            52.22 * 1000 * 7.88E-5,
            3.935 * 1000 * 7.88E-5,
            1.565 * 1000 * 7.88E-5,
            0.9939 * 1000 * 7.88E-5,

            45.31 * 1000 * 8.02E-5,
            3.637 * 1000 * 8.02E-5,
            1.209 * 1000 * 8.02E-5,
            0.4539 * 1000 * 8.02E-5
        ],
        "path": [
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_m125/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_VBF_m125/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_WH_m125/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_ZH_m125/",

            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_m120/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_VBF_m120/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_WH_m120/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_ZH_m120/",

            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_m130/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_VBF_m130/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_WH_m130/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Dalitz_eeg_ZH_m130/",
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_m125/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_VBF_m125/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_WH_m125/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_ZH_m125/",

            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_m120/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_VBF_m120/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_WH_m120/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_ZH_m120/",

            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_m130/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_VBF_m130/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_WH_m130/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Dalitz_eeg_ZH_m130/"
        ]
    },

    "DYJetsToLL_2017": {
        "production": ["DYJetsToLL"],
        "lumi": [41.525],
        "xs": [6077.22 * 1000.],
        "path": [
            "/data6/ggNtuples/V10_02_10_07/job_fall17_DYJetsToLL_m50_aMCatNLO*/"
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_DYJetsToLL_m50_aMCatNLO/"
        ]
    },

    "GJets_2017": {
        "production": ["gjet_pt15to6000", "gjet_pt20_MGG_40to80", "gjet_pt20to40_MGG_80toInf", "gjet_pt40_MGG_80toInf"],
        "lumi": [41.525],
        "xs": [
            290100 * 1000.,
            154500 * 0.0239 * 1000.,
            137751 * 0.0029 * 1000., 
            16792 * 0.0558 * 1000.
        ],
        "path": [
            "/data6/ggNtuples/V10_02_10_07/job_fall17_gjet_pt15to6000/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_gjet_pt20_MGG_40to80/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_gjet_pt20to40_MGG_80toInf/",
            "/data6/ggNtuples/V10_02_10_07/job_fall17_gjet_pt40_MGG_80toInf/"
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_gjet_pt15to6000/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_gjet_pt20_MGG_40to80/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_gjet_pt20to40_MGG_80toInf/",
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_gjet_pt40_MGG_80toInf/"
        ]
    }
}


#===============================================#
#              Data samples list                #
#  Run: 2016BCDEF1(pre_VFP), 2016F2GH(post_VFP) #
#          Run: 2017BCDEF, 2018ABCD             #
#===============================================#
DataSample = {
    "Data_2016_preVFP": {
        "run": ["2016B", "2016C", "2016D", "2016E", "2016F1"],
        "lumi": [
            4.671835672, 
            2.621295973, 
            4.285851496, 
            4.065974560, 
            2.717344923
        ],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2016B_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2016C_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2016D_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2016E_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2016F1_UL/"
        ],
        "outpath": [
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2016B_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2016C_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2016D_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2016E_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2016F1_UL/"
        ]
    },

    "Data_2016_postVFP": {
        "run": ["2016F2", "2016G", "2016H"],
        "lumi": [
            0.418120616, 
            7.652808345, 
            8.739883563
        ],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2016F2_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2016G_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2016H_UL/"
        ],
        "outpath": [
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2016F2_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2016G_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2016H_UL/"
        ]
    },

    "Data_2017": {
        "run": ["2017B", "2017C", "2017D", "2017E", "2017F"],
        "lumi": [
            4.803363110, 
            9.572498264, 
            4.247682053, 
            9.313642402, 
            13.539221759
        ],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2017B_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2017C_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2017D_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2017E_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_DoubleEG_Run2017F_UL/"
        ],
        "outpath": [
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2017B_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2017C_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2017D_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2017E_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_DoubleEG_Run2017F_UL/"
        ]
    },

    "Data_2018": {
        "run": ["2018A", "2018B", "2018C", "2018D"],
        "lumi": [
            14.026724850, 
            7.060617355, 
            6.894770971, 
            31.834115905
        ],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_20_00/job_EGamma_Run2018A_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_EGamma_Run2018B_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_EGamma_Run2018C_UL/",
            "/data3/ggNtuples/V10_06_20_00/job_EGamma_Run2018D_UL/"
        ],
        "outpath": [
            "/data4/chenghan/data/V10_06_20_00/job_EGamma_Run2018A_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_EGamma_Run2018B_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_EGamma_Run2018C_UL/",
            "/data4/chenghan/data/V10_06_20_00/job_EGamma_Run2018D_UL/"
        ]
    }
}
