import os, sys, re, json
import ROOT
from collections import OrderedDict as od


# sd = "systematics dataframe"


# Add column to dataFrame with default value for constant systematics:
# eg.
#   1) "name":"QCDscale_VH", "value":"WH:1.005/0.993,ZH:1.038/0.969" -> {WH:1.005/0.993, ZH:1.038/0.969} (proc as key) -> store in valueDict
#   2) "name":"BR_higgs_dalitz", "value":"0.94/1.06" -> {BR_higgs_dalitz:0.94/1.06} (name as key) -> store in onevalueDict
def addConstantSyst(sd, _syst):
    # extract the values map 
    valueDict_list = _syst["value"].split(",")
    valueDict, onevalueDict = od(), od()
    for i in valueDict_list:
        key_plus_value = i.split(":")
        if len(key_plus_value) > 1:
            valueDict[key_plus_value[0]] = key_plus_value[1]
        else:
            onevalueDict[_syst["name"]] = key_plus_value[0]

    if _syst["yearsDependence"] == 0:
        sd[_syst["name"]] = "-" # initial value
        for k, v in valueDict.iteritems():
            sd.loc[(sd["type"] == "sig")&(sd["proc"].str.contains(k)), _syst["name"]] = v
        for k, v in onevalueDict.iteritems():
            sd.loc[(sd["type"] == "sig"), k] = v
    
    elif _syst["yearsDependence"] == 1:
        for year, v in valueDict.iteritems():
            column_name = "{}_{}".format(_syst["name"], year)
            sd[column_name] = "-"

            # JEC and JER only for VBF categories
            if (_syst["name"] == "CMS_JEC_13TeV" or _syst["name"] == "CMS_JER_13TeV"):
                sd.loc[(sd["type"] == "sig")&(sd["cat"].str.contains("VBF"))&(year == sd["year"]), column_name] = v

            # R9 reweighting only for R9 related categories
            elif (_syst["name"] == "CMS_R9_13TeV"):
                sd.loc[(sd["type"] == "sig")&(sd["cat"].str.contains("R9"))&(year == sd["year"]), column_name] = v
                
            else:
                sd.loc[(sd["type"] == "sig")&(year == sd["year"]), column_name] = v

    else:
        sd[_syst["name"]] = "-"
        for year, v in valueDict.iteritems():

            # JEC and JER only for VBF categories
            if (_syst["name"] == "CMS_JEC_13TeV" or _syst["name"] == "CMS_JER_13TeV"):
                sd.loc[(sd["type"] == "sig")&(sd["cat"].str.contains("VBF")), _syst["name"]] = v
            
            # R9 reweighting only for R9 related categories
            elif (_syst["name"] == "CMS_R9_13TeV"):
                sd.loc[(sd["type"] == "sig")&(sd["cat"].str.contains("R9")), _syst["name"]] = v

            else:
                sd.loc[(sd["type"] == "sig")&(sd["year"] == year), _syst["name"]] = v

    return sd
            

