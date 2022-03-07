# Python file to store systematics: for Higgs Dalitz analysis

# THEORY SYSTEMATICS: The values are extracted from the origional hgg package json file
# ! Things to check:
#   1) the value of branching fraction
#   2) should we include alphaS or not?
theory_systematics = [
    {"name":"BR_higgs_dalitz",    "title":"BR_higgs_dalitz",    "type":"constant",  "prior":"lnN",  "yearsDependence":0,   "value":"0.94/1.06"}, 
    {"name":"pdf_Higgs_ggH",      "title":"pdf_Higgs_ggH",      "type":"constant",  "prior":"lnN",  "yearsDependence":0,   "value":"ggH:1.019"},
    {"name":"pdf_Higgs_qqbar",    "title":"pdf_Higgs_qqbar",    "type":"constant",  "prior":"lnN",  "yearsDependence":0,   "value":"qqH:1.021,WH:1.017,ZH:1.013"},
    {"name":"QCDscale_ggH",       "title":"QCDscale_ggH",       "type":"constant",  "prior":"lnN",  "yearsDependence":0,   "value":"ggH:1.047/0.931"},
    {"name":"QCDscale_qqH",       "title":"QCDscale_qqH",       "type":"constant",  "prior":"lnN",  "yearsDependence":0,   "value":"qqH:1.004/0.997"},
    {"name":"QCDscale_VH",        "title":"QCDscale_VH",        "type":"constant",  "prior":"lnN",  "yearsDependence":0,   "value":"WH:1.005/0.993,ZH:1.038/0.969"},
]

# EXPERIMENTAL SYSTEMATICS
# * yearsDependence = 0 : no dependence
# * yearsDependence = 1 : dependence and uncorrelated
# * yearsDependence = 2 : dependence and correlated
# ! Things to check:
#   1) estimate the experimaental uncertainties
experimental_systematics = [
    # Integrated luminosity 
    {"name":"lumi_13TeV",         "title":"lumi_13TeV",         "type":"constant",  "prior":"lnN",  "yearsDependence":1,   "value":"2016:1.022,2017:1.020,2018:1.015"},
    {"name":"lumi_13TeV_calib",   "title":"lumi_13TeV_calib",   "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.009,2017:1.008,2018:1.002"},
    {"name":"lumi_13TeV_deflect", "title":"lumi_13TeV_deflect", "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.004,2017:1.004,2018:1.000"},
    {"name":"lumi_13TeV_dynamic", "title":"lumi_13TeV_dynamic", "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.005,2017:1.005,2018:1.000"},
    {"name":"lumi_13TeV_sat",     "title":"lumi_13TeV_sat",     "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.000,2017:1.004,2018:1.001"},
    {"name":"lumi_13TeV_sc",      "title":"lumi_13TeV_sc",      "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.000,2017:1.003,2018:1.002"},
    {"name":"lumi_13TeV_xy",      "title":"lumi_13TeV_xy",      "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.009,2017:1.008,2018:1.02"},

    # ! FIXEDME: use 3% to estimate the following experimental uncertainties
    # Todo: HLT eff & Merged ID eff
    {"name":"CMS_HLTeff_13TeV",   "title":"CMS_HLTeff_13TeV",    "type":"constant",  "prior":"lnN",  "yearsDependence":1,   "value":"2016:1.030,2017:1.030,2018:1.030"},
    {"name":"CMS_IDeff_e_13TeV",  "title":"CMS_IDeff_e_13TeV",   "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.030,2017:1.030,2018:1.030"},
    {"name":"CMS_UE_13TeV",       "title":"CMS_UE_13TeV",        "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.030,2017:1.030,2018:1.030"},
    {"name":"CMS_PS_13TeV",       "title":"CMS_PS_13TeV",        "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.030,2017:1.030,2018:1.030"},
    
    # L1 prefiring weights, underlying event, parton showering, pileup reweighting, ID, HLT
    {"name":"CMS_L1_13TeV",       "title":"CMS_L1_13TeV",        "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.030,2017:1.030,2018:1.030"},
    {"name":"CMS_PU_13TeV",       "title":"CMS_PU_13TeV",        "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.030,2017:1.030,2018:1.030"},
    {"name":"CMS_IDeff_g_13TeV",  "title":"CMS_IDeff_g_13TeV",   "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.030,2017:1.030,2018:1.030"},
    
    # R9 reweighting: only for EBHR9, LR9 categories
    {"name":"CMS_R9_13TeV",       "title":"CMS_R9_13TeV",        "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.030,2017:1.030,2018:1.030"},
    
    # Jet energy scale/resolution: only for VBF categories
    {"name":"CMS_JER_13TeV",      "title":"CMS_JER_13TeV",       "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.030,2017:1.030,2018:1.030"},
    {"name":"CMS_JEC_13TeV",      "title":"CMS_JEC_13TeV",       "type":"constant",  "prior":"lnN",  "yearsDependence":2,   "value":"2016:1.030,2017:1.030,2018:1.030"},
]