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
    "ZGToLLG_UL2016preVFP": {
        "production": ["ZGToLLG"],
        "lumi": [19.52],
        "xs": [51.1 * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL16_Zg_aMCatNLO_preVFP/*.root"]
    },
    "ZGToLLG_UL2016postVFP": {
        "production": ["ZGToLLG"],
        "lumi": [16.81],
        "xs": [51.1 * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL16_Zg_aMCatNLO_postVFP/*.root"]
    },
    "ZGToLLG_UL2017": {
        "production": ["ZGToLLG"],
        "lumi": [41.48],
        "xs": [51.1 * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL17_Zg_aMCatNLO/*.root"]
    },
    "ZGToLLG_UL2018": {
        "production": ["ZGToLLG"],
        "lumi": [59.82],
        "xs": [51.1 * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL18_Zg_aMCatNLO/*.root"]
    },


    "TTJets_UL2016preVFP": {
        "production": ["TTJets"],
        "lumi": [19.52],
        "xs": [831.76 * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL16_TT_aMCatNLO_preVFP/*.root"]
    },
    "TTJets_UL2016postVFP": {
        "production": ["TTJets"],
        "lumi": [16.81],
        "xs": [831.76 * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL16_TT_aMCatNLO_postVFP/*.root"]
    },
    "TTJets_UL2017": {
        "production": ["TTJets"],
        "lumi": [41.48],
        "xs": [831.76 * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL17_TT_aMCatNLO/*.root"]
    },
    "TTJets_UL2018": {
        "production": ["TTJets"],
        "lumi": [59.82],
        "xs": [831.76 * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL18_TT_aMCatNLO/*.root"]
    },


    "DYJets_UL2016preVFP": {
        "production": ["DYJets"],
        "lumi": [19.52],
        "xs": [6404. * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL16_DYJetsToLL_m50_aMCatNLO_preVFP/*.root"]
    },
    "DYJets_UL2016postVFP": {
        "production": ["DYJets"],
        "lumi": [16.81],
        "xs": [6404. * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL16_DYJetsToLL_m50_aMCatNLO_postVFP/*.root"]
    },
    "DYJets_UL2017": {
        "production": ["DYJets"],
        "lumi": [41.48],
        "xs": [6435. * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL17_DYJetsToLL_m50_aMCatNLO/*.root"]
    },
    "DYJets_UL2018": {
        "production": ["DYJets"],
        "lumi": [59.82],
        "xs": [6529. * 1000.],
        "path": ["/data5/ggNtuples/V10_06_30_00/job_UL18_DYJetsToLL_m50_aMCatNLO/*.root"]
    },
}


#===============================================#
#              Data samples list                #
#  Run: 2016BCDEF1(pre_VFP), 2016F2GH(post_VFP) #
#          Run: 2017BCDEF, 2018ABCD             #
#===============================================#
DataSample = {
    "Data_UL2016preVFP": {
        "run": ["2016B", "2016C", "2016D", "2016E", "2016F1"],
        "lumi": [5.824235614, 2.621295973, 4.285851496, 4.065974639, 2.717344923],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2016B_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2016C_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2016D_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2016E_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2016F1_UL/*.root"
        ]
    },
    "Data_UL2016postVFP": {
        "run": ["2016F2", "2016G", "2016H"],
        "lumi": [0.418120616, 7.652808366, 8.739883636],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2016F2_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2016G_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2016H_UL/*.root"
        ]
    },
    "Data_UL2017": {
        "run": ["2017D", "2017E", "2017F"],
        "lumi": [4.247682053, 9.273176759, 13.538886094],
        "xs": [1.],
        "path": [
            # "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2017B_UL/*.root",
            # "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2017C_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2017D_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2017E_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2017F_UL/*.root"
        ]
    },
    "Data_UL2018": {
        "run": ["2018A", "2018B", "2018C", "2018D"],
        "lumi": [14.027047499, 7.055203602, 6.873175375, 31.834889477],
        "xs": [1.],
        "path": [
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2018A_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2018B_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2018C_UL/*.root",
            "/data3/ggNtuples/V10_06_30_00/job_SingleMu_Run2018D_UL/*.root"
        ]
    }
}