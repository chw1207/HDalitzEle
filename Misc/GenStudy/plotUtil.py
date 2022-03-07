import ROOT as rt

refXS = 48.58 + 3.782 + 1.373 + 0.8839

def gethist(channel, histname):
    chan = ""
    if channel == "ele":
        chan = "eeg"
    elif channel == "mu":
        chan = "mmg"
    else:
        print(
            '[WARNING] {} is neither [ele] nor [mu]. Use default [eeg].'.format(channel, ))
        chan = "eeg"

    l_prod = ['ggF', 'VBF', 'WH', 'ZH']
    # l_XSprod = [48.58, 3.782, 1.373, 0.8839]  # unit: pb
    hist = rt.TH1F()

    for iprod, prod in enumerate(l_prod):
        fin = rt.TFile(
            './minitree/Minitree_HDalitz_{}_{}_m125_2017_RECO.root'.format(prod, chan), "READ")
        if iprod == 0:
            hist = fin.Get(histname)
            hist.SetDirectory(0)
            # hist.Scale(l_XSprod[iprod] / refXS)
        else:
            hist_tmp = fin.Get(histname)
            hist_tmp.SetDirectory(0)
            # hist_tmp.Scale(l_XSprod[iprod] / refXS)
            hist.Add(hist_tmp)

        fin.Close()

    return hist

def gethist_single(fname, histname):
    hist = rt.TH1F()
    
    fin = rt.TFile(fname, "READ")
    hist = fin.Get(histname)
    hist.SetDirectory(0)
    fin.Close()
    return hist