from commonTools_HDalitz import massBaseList, years


# Input config file for running trees2ws_HDalitz.py

Info = {}

for mass in massBaseList:
    Info[mass] = {}
    for year in years:
        Info[mass][year] = {

            # Input root files which contain TTree
            "inputTreeFiles": [
                "~/Analysis/Dalitz/electron/miniTree/{}/miniTree_HDalitz_ggF_m{}.root".format(year, mass),
                "~/Analysis/Dalitz/electron/miniTree/{}/miniTree_HDalitz_VBF_m{}.root".format(year, mass),
                "~/Analysis/Dalitz/electron/miniTree/{}/miniTree_HDalitz_WH_m{}.root".format(year, mass),
                "~/Analysis/Dalitz/electron/miniTree/{}/miniTree_HDalitz_ZH_m{}.root".format(year, mass),
            ],

            # Name of the input tree(miniTree produced by xAna)
            "inputTreeName":  "outTree",

            # Output root files which contain WS
            "outputWSFile": [
                "./WS/{}/signal_ggF_m{}.root".format(year, mass),
                "./WS/{}/signal_VBF_m{}.root".format(year, mass),
                "./WS/{}/signal_WH_m{}.root".format(year, mass),
                "./WS/{}/signal_ZH_m{}.root".format(year, mass),
            ],
            

            # mass point
            "MassPoint": mass,

            # mass cut 
            "massBounds": [105, 170],
            
            # Vars in the minitree(they corresponds to the WSVars)
            # ! For "TreeVars" and "WSVars", assuming the first var is mass and the second one is the weight
            "TreeVars": ["CMS_higgs_mass", "weight"],

            # Vars to be added to work space
            # ! Don't change "WSVars"
            "WSVars": ["CMS_higgs_mass", "weight"], 

            # Variables to add to sytematic
            "systematicsVars": ["CMS_higgs_mass", "weight"],
            
            # List of systematics:
#             "sysTreeName":["unTreePhoR9"]
# {{"UnPhoR9", "unTreePhoR9"}, {"UnJERUp", "unTreeJERUp"}, {"UnJERDo", "unTreeJERDo"}, {"UnJECUp", "unTreeJECUp"}, {"UnJECDo", "unTreeJECDo"}};

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
