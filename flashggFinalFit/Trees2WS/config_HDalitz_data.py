# Input config file for running trees2ws_HDalitz_data.py

Info = { 

    "MergeYears": True,
    
    # Input root files which contain TTree
    "inputTreeFiles": [
        "~/Analysis/Dalitz/electron/miniTree/2016/miniTree_Data_2016.root",
        "~/Analysis/Dalitz/electron/miniTree/2017/miniTree_Data_2017.root",
        "~/Analysis/Dalitz/electron/miniTree/2018/miniTree_Data_2018.root",
    ],

    # Name of the input tree(miniTree produced by xAna)
    "inputTreeName":  "outTree",

    # Output root files which contain WS
    "outputWSFiles": [
        "./WS/2016/data_obs_2016.root",
        "./WS/2017/data_obs_2017.root",
        "./WS/2018/data_obs_2018.root",
    ],

    # used to put the ws for merging all 3 years if MergeYears is true.
    "outputWSFile": "./WS/data_obs.root",
    
    # mass point
    "MassPoint": 125, # useless for data

    # mass cut 
    "massBounds": [105, 170],
    
    # Vars in the minitree(they corresponds to the WSVars)
    # ! For "TreeVars" and "WSVars", assuming the first var is mass and the second one is the weight
    # ! The weight of data is assigned as 1. in xAna
    "TreeVars": ["CMS_higgs_mass", "weight"],

    # Vars to be added to work space
    # ! Don't change "WSVars"
    "WSVars": ["CMS_higgs_mass", "weight"], 

    # category branch in the minitree 
    #   1. Assuming the sequence of the number stored in the branch corresponds to the "cats" 
    #   2. 1 -> "Merge2Gsf_HVBF", 2 -> "Merge2Gsf_LVBF" ....
    "catBranch": "category", 

    # Analysis categories
    "cats": [
        "Merged2Gsf_HVBF", 
        "Merged2Gsf_LVBF", 
        "Merged2Gsf_BST", 
        "Merged2Gsf_EBHR9", 
        "Merged2Gsf_EBLR9", 
        "Merged2Gsf_EE",
        "Merged1Gsf_HVBF", 
        "Merged1Gsf_LVBF", 
        "Merged1Gsf_BST", 
        "Merged1Gsf_EBHR9", 
        "Merged1Gsf_EBLR9", 
        "Merged1Gsf_EE",
        "Resolved"
    ]
}