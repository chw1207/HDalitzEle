Float_t Lumi_2017 = 4.793 + 9.631 + 4.248 + 9.314 + 13.539;
Float_t Lumi_2016 = 35.917;
Float_t Lumi_2018 = 59.725;
Float_t Lumi_FullRun2 = Lumi_2017 + Lumi_2016 + Lumi_2018; // invfb

Float_t mcwei(TreeReader &data, Float_t XS, Float_t lumi, Int_t aMCatNLO, Long64_t &TotEv)
{
    Float_t mcwei = 1.;
    TotEv = 0;
    if (!data.HasMC())
    {
        return 1.;
    }
    else
    {
        for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++)
        {
            data.GetEntry(ev);
            float genWeight = data.GetFloat("genWeight");
            Int_t nVtx = data.GetInt("nVtx");

            // if (nVtx <= 1)
            //     continue;

            if (genWeight > 0)
                TotEv++;
            else
                TotEv--;
        }

        if (aMCatNLO == 1)
            mcwei = (TotEv != 0) ? XS * lumi / TotEv : 1;
        else
            mcwei = XS * lumi / (double)data.GetEntriesFast();
        
        return mcwei;
    }
}

Float_t genlepmass(Int_t lepID)
{
    if (abs(lepID) == 11)
        return 0.000510999;
    else if (abs(lepID) == 13)
        return 0.105658;
    else if (abs(lepID) == 15)
        return 1.77686;
    else
    {
        cout << "[WARNING] " << lepID << " is not lepton! Please correct it!" << endl;
        return -999.;
    }
}

std::map<string, Float_t> XS_HDalitz() // fb
{
    std::map<string, Float_t> XSmap;
    XSmap["ggF_125GeV"] = 48.58 * 1000 * 8.10E-5;
    XSmap["VBF_125GeV"] = 3.782 * 1000 * 8.10E-5;
    XSmap["WH_125GeV"] = 1.373 * 1000 * 8.10E-5;
    XSmap["ZH_125GeV"] = 0.8839 * 1000 * 8.10E-5;

    XSmap["ggF_120GeV"] = 52.22 * 1000 * 7.88E-5;
    XSmap["VBF_120GeV"] = 3.935 * 1000 * 7.88E-5;
    XSmap["WH_120GeV"] = 1.565 * 1000 * 7.88E-5;
    XSmap["ZH_120GeV"] = 0.9939 * 1000 * 7.88E-5;

    XSmap["ggF_130GeV"] = 45.31 * 1000 * 8.02E-5;
    XSmap["VBF_130GeV"] = 3.637 * 1000 * 8.02E-5;
    XSmap["WH_130GeV"] = 1.209 * 1000 * 8.02E-5;
    XSmap["ZH_130GeV"] = 0.4539 * 1000 * 8.02E-5;
    return XSmap;
}

std::map<string, Float_t> XS_GJets() // fb
{
    std::map<string, Float_t> XSmap;
    // From McM
    XSmap["pt15to6000"] = 290100 * 1000.;
    XSmap["pt20_MGG_40to80"] = 154500 * 0.0239 * 1000.;
    XSmap["pt20to40_MGG_80toInf"] = 137751 * 0.0029 * 1000.;
    XSmap["pt40_MGG_80toInf"] = 16792 * 0.0558 * 1000.;
    return XSmap;
}

// From EXO-20-005 AN2019_237
std::map<string, double> XS_QCD() // fb
{
    std::map<string, double> XSmap;
    XSmap["HT100to200"] = 27990000 * 1000.;
    XSmap["HT200to300"] = 1735000 * 1000.;
    XSmap["HT300to500"] = 366800 * 1000.;
    XSmap["HT500to700"] = 29370 * 1000.;
    XSmap["HT700to1000"] = 6524 * 1000.;
    XSmap["HT1000to1500"] = 1064 * 1000.;
    XSmap["HT1500to2000"] = 121.5 * 1000.;
    XSmap["HT2000toInf"] = 25.42 * 1000.;
    return XSmap;
}

Float_t luminosity(const char *runyear)
{
    if (strcmp(runyear, "2016") == 0)
        return Lumi_2016;
    else if (strcmp(runyear, "2017") == 0)
        return Lumi_2017;
    else if (strcmp(runyear, "2018") == 0)
        return Lumi_2018;
    else if (strcmp(runyear, "FullRun2") == 0)
        return Lumi_FullRun2;
    else 
    {
        cout << "[WARNING] runyear " << runyear << " is not correct. Please check! 1 fbinv is used for now." << endl;
        return 1.;
    }
}
