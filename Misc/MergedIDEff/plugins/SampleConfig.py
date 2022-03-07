
#===============================================#
#              Merged ID features               #
#              Merged ID models                 #
#===============================================#
features_m1 = [
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
    "M1EB": "../ID_Trainer/Results/Output_Merged1GsfID_EB_2017_daskNoGJets/XGB/XGB_modelXGB.pkl",
    "M2EB": "../ID_Trainer/Results/Output_Merged2GsfID_EB_2017_daskNoGJets/XGB/XGB_modelXGB.pkl",
    "M1EE": "../ID_Trainer/Results/Output_Merged1GsfID_EE_2017_daskNoGJets/XGB/XGB_modelXGB.pkl",
    "M2EE": "../ID_Trainer/Results/Output_Merged2GsfID_EE_2017_daskNoGJets/XGB/XGB_modelXGB.pkl"
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
    # reference: https://cms.cern.ch/iCMS/jsp/db_notes/noteInfo.jsp?cmsnoteid=CMS%20AN-2019/139
    "ZZ_2017": {
        "production": ["ZZ"],
        "lumi": [41.525],
        "xs": [1.212 * 1000.],
        "path": [
            "/data6/ggNtuples/V10_02_10_07/job_fall17_ZZ/",
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_ZZ/",
        ]
    },

    "ZZ_2018": {
        "production": ["ZZ"],
        "lumi": [59.7],
        "xs": [1.212 * 1000.],
        "path": [
            "/data6/ggNtuples/V10_02_10_07/job_autumn18_ZZ/",
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_07/job_autumn18_ZZ/",
        ]
    },

    "ZZ_2016": {
        "production": ["ZZ"],
        "lumi": [35.9],
        "xs": [1.212 * 1000.],
        "path": [
            "/data6/ggNtuples/V10_02_10_07/job_summer16_ZZ/",
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_07/job_summer16_ZZ/",
        ]
    },

    "Zg_2017": {
        "production": ["Zg"],
        "lumi": [41.525],
        "xs": [106.17*1000.],
        "path": [
            "/data6/ggNtuples/V10_02_10_07/job_fall17_Zg_lowMll_lowGPt_aMCatNLO/",
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_07/job_fall17_Zg_lowMll_lowGPt_aMCatNLO/",
        ]
    },

    "Zg_2016": {
        "production": ["Zg"],
        "lumi": [35.9],
        "xs": [106.17*1000.],
        "path": [
            "/data6/ggNtuples/V10_02_10_07/job_summer16_Zg_lowMll_lowGPt_aMCatNLO/",
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_07/job_summer16_Zg_lowMll_lowGPt_aMCatNLO/",
        ]
    },

    "Zg_2018": {
        "production": ["Zg"],
        "lumi": [59.7],
        "xs": [106.17*1000.],
        "path": [
            "/data6/ggNtuples/V10_02_10_07/job_autumn18_Zg_lowMll_lowGPt_aMCatNLO/",
        ],
        "outpath": [
            "/data4/chenghan/mc/V10_02_10_07/job_autumn18_Zg_lowMll_lowGPt_aMCatNLO/",
        ]
    },
}


#===============================================#
#              Data samples list                #
#  Run: 2016BCDEF1(pre_VFP), 2016F2GH(post_VFP) #
#          Run: 2017BCDEF, 2018ABCD             #
#===============================================#
DataSample = {

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
}
