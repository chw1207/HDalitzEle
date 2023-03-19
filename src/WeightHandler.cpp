#include <string>
#include <mutex>
#include <vector>
#include <map>
#include "WeightHandler.h"
#include "PUWeightCalculator.h"
#include "Generator.h"
#include "Utilities.h"
#include "TString.h"


ROOT::RDF::RNode wei::ApplyWeight(ROOT::RDF::RNode df, const YAML::Node cfg){
    // auto filters = df.GetFilterNames();
    // if ((int)filters.size() > 0)
    //     throw std::runtime_error("This function should be applid in the begining(i.e. no event is filterd out)");

    // //! [NOTE] these filters only needed for signal MC samples!!
    // auto ff = df.Define("LHEMEE",           gen::CalcLHEMee,    {"nLHE", "lhePID", "lhePx", "lhePy", "lhePz"})
    //             .Define("HasPhoIntConv",    gen::IsPhoIntConv,  {"nMC", "mcPID", "mcMomPID", "mcStatusFlag"})
    //             .Filter("event % 2 == 1", "split smaples") // == 1 for analysis, == 0 for ID training
    //             .Filter("LHEMEE < 60",    "remove HZg")
    //             .Filter("HasPhoIntConv != 1", "remove phoIntCov");
    // auto ff = df;
    // calculate the mc weight
    // printf(" [+] Weight application:\n");
    // auto pos = ff.Filter("genWeight > 0",   "positive event").Count();
    // auto neg = ff.Filter("genWeight <= 0",  "negative event").Count();
    // const int totalev = pos.GetValue() - neg.GetValue();
    // const float xs = cfg["cross_section"].as<std::vector<float>>().at(iset);
    // const float lumi = cfg["luminosity"].as<float>();
    // const float mcwei = ((xs * lumi) / totalev);
    // printf("     - XS = %f fb, Luminosity = %f fb^-1, Nevents = %d \n", xs, lumi, totalev);
    // printf("     - MC weight (XS * Luminosity / Nevents): %f \n", mcwei);

    // pu rewei calculation
    auto pufiles = cfg["external_files"]["pileup_rewei"].as<std::map<std::string, std::string>>();
    std::vector<std::shared_ptr<PUWeightCalculator>> puCalc = {
        std::make_shared<PUWeightCalculator>(),
        std::make_shared<PUWeightCalculator>(),
        std::make_shared<PUWeightCalculator>()
    };
    puCalc[0]->Init(pufiles["nominal"].c_str());
    puCalc[1]->Init(pufiles["up"].c_str());
    puCalc[2]->Init(pufiles["down"].c_str());

    auto get_pu = [puCalc](const int run, const ROOT::RVec<float>& puTrue){
        float puwei = puCalc[0]->GetWeight(run, puTrue[0]);
        return puwei;
    };
    auto get_pu_up = [puCalc](const int run, const ROOT::RVec<float>& puTrue){
        float puwei_up = puCalc[1]->GetWeight(run, puTrue[0]);
        return puwei_up;
    };
    auto get_pu_do = [puCalc](const int run, const ROOT::RVec<float>& puTrue){
        float puwei_down = puCalc[2]->GetWeight(run, puTrue[0]);
        return puwei_down;
    };

    auto nf = df.Define("genwei",           "if (genWeight > 0) return 1.; else return -1.;")
                .Define("puwei",            get_pu,             {"run", "puTrue"})
                .Define("puwei_up",         get_pu_up,          {"run", "puTrue"})
                .Define("puwei_down",       get_pu_do,          {"run", "puTrue"});

    return nf;
}


ROOT::RDF::RNode wei::GetScaleFactors(ROOT::RDF::RNode df, const YAML::Node cfg){

    const YAML::Node sfs = cfg["external_files"]["scale_factors"];

    // reconstruction scale factors
    auto fNameRECOEle = sfs["RECOEle"]["file"].as<std::string>();
    auto hNameRECOEle = sfs["RECOEle"]["hist"].as<std::string>();
    std::unique_ptr<TFile> fRECOEle(TFile::Open(fNameRECOEle.c_str(), "READ"));
    if (!fRECOEle || fRECOEle->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", fNameRECOEle.c_str()));
    std::shared_ptr<TH2F> hRECOEle((TH2F*) fRECOEle->Get(hNameRECOEle.c_str()));
    if (!hRECOEle)
        throw std::runtime_error(Form("TFile::Get() failed: %s", hNameRECOEle.c_str()));
    hRECOEle->SetDirectory(nullptr);
    auto get_RECOEleSF = [&, hRECOEle](const float eleSCEta, const float elePt){
        float SF = hRECOEle->GetBinContent(hRECOEle->FindBin(eleSCEta, elePt));
        return SF;
    };
    auto get_RECOEleSFErr = [&, hRECOEle](const float eleSCEta, const float elePt){
        float SF_err = hRECOEle->GetBinError(hRECOEle->FindBin(eleSCEta, elePt));
        return SF_err;
    };

    // Fall 17 electron MVA ID SFs (for resolved electrons)
    auto fNameFall17EleID = sfs["Fall17EleID"]["file"].as<std::string>();
    auto hNameFall17EleID = sfs["Fall17EleID"]["hist"].as<std::string>();
    std::unique_ptr<TFile> fFall17EleID(TFile::Open(fNameFall17EleID.c_str(), "READ"));
    if (!fFall17EleID || fFall17EleID->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", fNameFall17EleID.c_str()));
    std::shared_ptr<TH2F> hFall17EleID((TH2F*) fFall17EleID->Get(hNameFall17EleID.c_str()));
    if (!hFall17EleID)
        throw std::runtime_error(Form("TFile::Get() failed:  %s", hNameFall17EleID.c_str()));
    hFall17EleID->SetDirectory(nullptr);
    auto get_Fall17EleIDSF = [&, hFall17EleID](const float eleSCEta, const float elePt){
        float SF = hFall17EleID->GetBinContent(hFall17EleID->FindBin(eleSCEta, elePt));
        return SF;
    };
    auto get_Fall17EleIDSFErr = [&, hFall17EleID](const float eleSCEta, const float elePt){
        float SF_err = hFall17EleID->GetBinError(hFall17EleID->FindBin(eleSCEta, elePt));
        return SF_err;
    };

    // Fall 17 photon MVA ID SFs
    auto fNameFall17PhoID = sfs["Fall17PhoID"]["file"].as<std::string>();
    auto hNameFall17PhoID = sfs["Fall17PhoID"]["hist"].as<std::string>();
    std::unique_ptr<TFile> fFall17PhoID(TFile::Open(fNameFall17PhoID.c_str(), "READ"));
    if (!fFall17PhoID || fFall17PhoID->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", fNameFall17PhoID.c_str()));
    std::shared_ptr<TH2F> hFall17PhoID((TH2F*) fFall17PhoID->Get(hNameFall17PhoID.c_str()));
    if (!hFall17PhoID)
        throw std::runtime_error(Form("TFile::Get() failed: %s", hNameFall17PhoID.c_str()));
    hFall17PhoID->SetDirectory(nullptr);
    auto get_Fall17PhoIDSF = [&, hFall17PhoID](const float phoSCEta, const float phoPt){
        float SF = hFall17PhoID->GetBinContent(hFall17PhoID->FindBin(phoSCEta, phoPt));
        return SF;
    };
    auto get_Fall17PhoIDSFErr = [&, hFall17PhoID](const float phoSCEta, const float phoPt){
        float SF_err = hFall17PhoID->GetBinError(hFall17PhoID->FindBin(phoSCEta, phoPt));
        return SF_err;
    };

    // CSEV SFs (conversion safe electron veto)
    //! NOTE: SFs histogram is TH1F instead of TH2F for CSEV
    auto fNameCSEV = sfs["CSEV"]["file"].as<std::string>();
    auto hNameCSEV = sfs["CSEV"]["hist"].as<std::string>();
    std::unique_ptr<TFile> fCSEV(TFile::Open(fNameCSEV.c_str(), "READ"));
    if (!fCSEV || fCSEV->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", fNameCSEV.c_str()));
    std::shared_ptr<TH1F> hCSEV((TH1F*) fCSEV->Get(hNameCSEV.c_str()));
    if (!hCSEV)
        throw std::runtime_error(Form("TFile::Get() failed: %s", hNameCSEV.c_str()));
    hCSEV->SetDirectory(nullptr);
    auto get_CSEVSF = [&, hCSEV](const float phoSCEta, const float phoR9){
        float SF = 1.;
        if (fabs(phoSCEta) < 1.4446 && phoR9 > 0.96)
            SF = hCSEV->GetBinContent(2);
        else if (fabs(phoSCEta) < 1.4446 && phoR9 < 0.96)
            SF = hCSEV->GetBinContent(3);
        else if (fabs(phoSCEta) > 1.566 && phoR9 > 0.96)
            SF = hCSEV->GetBinContent(5);
        else if (fabs(phoSCEta) > 1.566 && phoR9 < 0.96)
            SF = hCSEV->GetBinContent(6);
        return SF;
    };
    auto get_CSEVSFErr = [&, hCSEV](const float phoSCEta, const float phoR9){
        float SF_err = 1.;
        if (fabs(phoSCEta) < 1.4446 && phoR9 > 0.96)
            SF_err = hCSEV->GetBinError(2);
        else if (fabs(phoSCEta) < 1.4446 && phoR9 < 0.96)
            SF_err = hCSEV->GetBinError(3);
        else if (fabs(phoSCEta) > 1.566 && phoR9 > 0.96)
            SF_err = hCSEV->GetBinError(5);
        else if (fabs(phoSCEta) > 1.566 && phoR9 < 0.96)
            SF_err = hCSEV->GetBinError(6);
        return SF_err;
    };

    // H -> gg preselection SFs
    auto etabins = sfs["HggPresel"]["absEta"].as<std::vector<float>>();
    auto r9bins  = sfs["HggPresel"]["fullR9"].as<std::vector<std::vector<float>>>();
    auto values  = sfs["HggPresel"]["values"].as<std::vector<std::vector<float>>>();
    auto uncert  = sfs["HggPresel"]["uncert"].as<std::vector<std::vector<float>>>();
    auto get_HggPreselALL = [&, etabins, r9bins, values, uncert](const float sceta, const float full5x5r9){
        std::vector<float> hggsf(2);
        int etabin = utils::FindBins(etabins, fabs(sceta));
        if (etabin == -1)
            throw std::runtime_error(Form("SCEta = %f, out of range of HggPresel absEta bins", fabs(sceta)));

        int r9bin  = utils::FindBins(r9bins.at(etabin), full5x5r9);
        if (r9bin == -1)
            throw std::runtime_error(Form("full5x5r9 = %f, out of range of HggPresel fullR9 bins", full5x5r9));

        hggsf[0] = values.at(etabin).at(r9bin);
        hggsf[1] = uncert.at(etabin).at(r9bin);
        return hggsf;
    };

    auto ff = df.Define("RECOEleSF_Lead",               get_RECOEleSF,                {"eleSCEta_Lead", "elePt_Lead"})
                .Define("RECOEleSFErr_Lead",            get_RECOEleSFErr,             {"eleSCEta_Lead", "elePt_Lead"})
                .Define("RECOEleSFUp_Lead",             "RECOEleSF_Lead + RECOEleSFErr_Lead")
                .Define("RECOEleSFDo_Lead",             "RECOEleSF_Lead - RECOEleSFErr_Lead")

                .Define("RECOEleSF_subLead",            get_RECOEleSF,                {"eleSCEta_subLead", "elePt_subLead"})
                .Define("RECOEleSFErr_subLead",         get_RECOEleSFErr,             {"eleSCEta_subLead", "elePt_subLead"})
                .Define("RECOEleSFUp_subLead",          "RECOEleSF_subLead + RECOEleSFErr_subLead")
                .Define("RECOEleSFDo_subLead",          "RECOEleSF_subLead - RECOEleSFErr_subLead")

                .Define("Fall17EleIDSF_Lead",           get_Fall17EleIDSF,            {"eleSCEta_Lead", "elePt_Lead"})
                .Define("Fall17EleIDSFErr_Lead",        get_Fall17EleIDSFErr,         {"eleSCEta_Lead", "elePt_Lead"})
                .Define("Fall17EleIDSFUp_Lead",         "Fall17EleIDSF_Lead + Fall17EleIDSFErr_Lead")
                .Define("Fall17EleIDSFDo_Lead",         "Fall17EleIDSF_Lead - Fall17EleIDSFErr_Lead")

                // Hgg SFs are only for Lead electrons (merged electrons)
                .Define("HggPreselEleALL",              get_HggPreselALL,             {"eleSCEta_Lead", "eleR9Full5x5_Lead"})
                .Define("HggPreselEleSF_Lead",          "HggPreselEleALL[0]")
                .Define("HggPreselEleSFErr_Lead",       "HggPreselEleALL[1]")
                .Define("HggPreselEleSFUp_Lead",        "HggPreselEleSF_Lead + HggPreselEleSFErr_Lead")
                .Define("HggPreselEleSFDo_Lead",        "HggPreselEleSF_Lead - HggPreselEleSFErr_Lead")

                .Define("Fall17EleIDSF_subLead",        get_Fall17EleIDSF,            {"eleSCEta_subLead", "elePt_subLead"})
                .Define("Fall17EleIDSFErr_subLead",     get_Fall17EleIDSFErr,         {"eleSCEta_subLead", "elePt_subLead"})
                .Define("Fall17EleIDSFUp_subLead",      "Fall17EleIDSF_subLead + Fall17EleIDSFErr_subLead")
                .Define("Fall17EleIDSFDo_subLead",      "Fall17EleIDSF_subLead - Fall17EleIDSFErr_subLead")

                .Define("Fall17PhoIDSF_Lead",           get_Fall17PhoIDSF,            {"phoSCEta_Lead", "phoEt_Lead"})
                .Define("Fall17PhoIDSFErr_Lead",        get_Fall17PhoIDSFErr,         {"phoSCEta_Lead", "phoEt_Lead"})
                .Define("Fall17PhoIDSFUp_Lead",         "Fall17PhoIDSF_Lead + Fall17PhoIDSFErr_Lead")
                .Define("Fall17PhoIDSFDo_Lead",         "Fall17PhoIDSF_Lead - Fall17PhoIDSFErr_Lead")

                .Define("CSEVSF_Lead",                  get_CSEVSF,                   {"phoSCEta_Lead", "phoR9_Lead"})
                .Define("CSEVSFErr_Lead",               get_CSEVSFErr,                {"phoSCEta_Lead", "phoR9_Lead"})
                .Define("CSEVSFUp_Lead",                "CSEVSF_Lead + CSEVSFErr_Lead")
                .Define("CSEVSFDo_Lead",                "CSEVSF_Lead - CSEVSFErr_Lead")

                .Define("HggPreselPhoALL",              get_HggPreselALL,             {"phoSCEta_Lead", "phoR9Full5x5_Lead"})
                .Define("HggPreselPhoSF_Lead",          "HggPreselPhoALL[0]")
                .Define("HggPreselPhoSFErr_Lead",       "HggPreselPhoALL[1]")
                .Define("HggPreselPhoSFUp_Lead",        "HggPreselPhoSF_Lead + HggPreselPhoSFErr_Lead")
                .Define("HggPreselPhoSFDo_Lead",        "HggPreselPhoSF_Lead - HggPreselPhoSFErr_Lead");

    //! NOTE: To be estimated, temporarily assign as 0.9 and with 5% uncertainties
    auto nf = ff.Define("MergedEleIDSF_Lead",           "(float) 0.9")
                .Define("MergedEleIDSFErr_Lead",        "(float) 0.05")
                .Define("MergedEleIDSFUp_Lead",         "MergedEleIDSF_Lead + MergedEleIDSFErr_Lead")
                .Define("MergedEleIDSFDo_Lead",         "MergedEleIDSF_Lead - MergedEleIDSFErr_Lead")

                .Define("ConvVetoEleSF_Lead",           "(float) 0.9")
                .Define("ConvVetoEleSFErr_Lead",        "(float) 0.05")
                .Define("ConvVetoEleSFUp_Lead",         "ConvVetoEleSF_Lead + ConvVetoEleSFErr_Lead")
                .Define("ConvVetoEleSFDo_Lead",         "ConvVetoEleSF_Lead - ConvVetoEleSFErr_Lead")

                .Define("MissHitsTrkSF_Lead",           "(float) 0.9")
                .Define("MissHitsTrkSFErr_Lead",        "(float) 0.05")
                .Define("MissHitsTrkSFUp_Lead",         "MissHitsTrkSF_Lead + MissHitsTrkSFErr_Lead")
                .Define("MissHitsTrkSFDo_Lead",         "MissHitsTrkSF_Lead - MissHitsTrkSFErr_Lead")

                .Define("MissHitsSubTrkSF_Lead",        "(float) 0.9")
                .Define("MissHitsSubTrkSFErr_Lead",     "(float) 0.05")
                .Define("MissHitsSubTrkSFUp_Lead",      "MissHitsSubTrkSF_Lead + MissHitsSubTrkSFErr_Lead")
                .Define("MissHitsSubTrkSFDo_Lead",      "MissHitsSubTrkSF_Lead - MissHitsSubTrkSFErr_Lead")

                .Define("DiPhoHLTSF",                   "(float) 0.9")
                .Define("DiPhoHLTSFErr",                "(float) 0.05")
                .Define("DiPhoHLTSFUp",                 "DiPhoHLTSF + DiPhoHLTSFErr")
                .Define("DiPhoHLTSFDo",                 "DiPhoHLTSF - DiPhoHLTSFErr")

                .Define("DiEleHLTSF",                   "(float) 0.9")
                .Define("DiEleHLTSFErr",                "(float) 0.05")
                .Define("DiEleHLTSFUp",                 "DiEleHLTSF + DiEleHLTSFErr")
                .Define("DiEleHLTSFDo",                 "DiEleHLTSF - DiEleHLTSFErr")

                .Define("HLTSF",                        "if (isRe) return DiEleHLTSF; else return DiPhoHLTSF")
                .Define("HLTSFUp",                      "if (isRe) return DiEleHLTSFUp; else return DiPhoHLTSFUp")
                .Define("HLTSFDo",                      "if (isRe) return DiEleHLTSFDo; else return DiPhoHLTSFDo");
    return nf;
}


ROOT::RDF::RNode wei::GetFinalWeights(ROOT::RDF::RNode df){

    const char* merged2g_nominal = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_puweiUp = "mcwei * genwei * puwei_up   * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_puweiDo = "mcwei * genwei * puwei_down * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_l1pfUp  = "mcwei * genwei * puwei      * L1ECALPrefireUp   * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_l1pfDo  = "mcwei * genwei * puwei      * L1ECALPrefireDown * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_hltUp   = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSFUp * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_hltDo   = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSFDo * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_recoEUp = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSFUp_Lead * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_recoEDo = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSFDo_Lead * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_eleIDUp = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSFUp_Lead * MergedEleIDSFUp_Lead * ConvVetoEleSFUp_Lead * MissHitsTrkSFUp_Lead * MissHitsSubTrkSFUp_Lead * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_eleIDDo = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSFDo_Lead * MergedEleIDSFDo_Lead * ConvVetoEleSFDo_Lead * MissHitsTrkSFDo_Lead * MissHitsSubTrkSFDo_Lead * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged2g_phoIDUp = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSFUp_Lead * Fall17PhoIDSFUp_Lead * CSEVSF_Lead";
    const char* merged2g_phoIDDo = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSFDo_Lead * Fall17PhoIDSFDo_Lead * CSEVSF_Lead";
    const char* merged2g_csevUp  = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSFUp_Lead";
    const char* merged2g_csevDo  = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * MissHitsSubTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSFDo_Lead";

    // no MissHitsSubTrkSF comapred with merged2g
    const char* merged1g_nominal = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_puweiUp = "mcwei * genwei * puwei_up   * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_puweiDo = "mcwei * genwei * puwei_down * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_l1pfUp  = "mcwei * genwei * puwei      * L1ECALPrefireUp   * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_l1pfDo  = "mcwei * genwei * puwei      * L1ECALPrefireDown * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_hltUp   = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSFUp * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_hltDo   = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSFDo * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_recoEUp = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSFUp_Lead * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_recoEDo = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSFDo_Lead * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_eleIDUp = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSFUp_Lead * MergedEleIDSFUp_Lead * ConvVetoEleSFUp_Lead * MissHitsTrkSFUp_Lead * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_eleIDDo = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSFDo_Lead * MergedEleIDSFDo_Lead * ConvVetoEleSFDo_Lead * MissHitsTrkSFDo_Lead * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* merged1g_phoIDUp = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSFUp_Lead * Fall17PhoIDSFUp_Lead * CSEVSF_Lead";
    const char* merged1g_phoIDDo = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSFDo_Lead * Fall17PhoIDSFDo_Lead * CSEVSF_Lead";
    const char* merged1g_csevUp  = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSFUp_Lead";
    const char* merged1g_csevDo  = "mcwei * genwei * puwei      * L1ECALPrefire     * DiPhoHLTSF   * RECOEleSF_Lead   * HggPreselEleSF_Lead   * MergedEleIDSF_Lead   * ConvVetoEleSF_Lead   * MissHitsTrkSF_Lead   * HggPreselPhoSF_Lead   * Fall17PhoIDSF_Lead   * CSEVSFDo_Lead";

    const char* resolved_nominal = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_puweiUp = "mcwei * genwei * puwei_up   * L1ECALPrefire     * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_puweiDo = "mcwei * genwei * puwei_down * L1ECALPrefire     * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_l1pfUp  = "mcwei * genwei * puwei      * L1ECALPrefireUp   * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_l1pfDo  = "mcwei * genwei * puwei      * L1ECALPrefireDown * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_hltUp   = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSFUp * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_hltDo   = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSFDo * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_recoEUp = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSF   * RECOEleSFUp_Lead * RECOEleSFUp_subLead * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_recoEDo = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSF   * RECOEleSFDo_Lead * RECOEleSFDo_subLead * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_eleIDUp = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSFUp_Lead * Fall17EleIDSFUp_subLead * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_eleIDDo = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSFDo_Lead * Fall17EleIDSFDo_subLead * Fall17PhoIDSF_Lead   * CSEVSF_Lead";
    const char* resolved_phoIDUp = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSFUp_Lead * CSEVSF_Lead";
    const char* resolved_phoIDDo = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSFDo_Lead * CSEVSF_Lead";
    const char* resolved_csevUp  = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSFUp_Lead";
    const char* resolved_csevDo  = "mcwei * genwei * puwei      * L1ECALPrefire     * DiEleHLTSF   * RECOEleSF_Lead   * RECOEleSF_subLead   * Fall17EleIDSF_Lead   * Fall17EleIDSF_subLead   * Fall17PhoIDSF_Lead   * CSEVSFDo_Lead";

    auto nf = df.Define("weight",           Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_nominal,    merged1g_nominal,    resolved_nominal))
                .Define("weight_puweiUp",   Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_puweiUp,    merged1g_puweiUp,    resolved_puweiUp))
                .Define("weight_puweiDo",   Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_puweiDo,    merged1g_puweiDo,    resolved_puweiDo))
                .Define("weight_l1pfUp",    Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_l1pfUp,     merged1g_l1pfUp,     resolved_l1pfUp))
                .Define("weight_l1pfDo",    Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_l1pfDo,     merged1g_l1pfDo,     resolved_l1pfDo))
                .Define("weight_hltUp",     Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_hltUp,      merged1g_hltUp,      resolved_hltUp))
                .Define("weight_hltDo",     Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_hltDo,      merged1g_hltDo,      resolved_hltDo))
                .Define("weight_recoEUp",   Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_recoEUp,    merged1g_recoEUp,    resolved_recoEUp))
                .Define("weight_recoEDo",   Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_recoEDo,    merged1g_recoEDo,    resolved_recoEDo))
                .Define("weight_eleIDUp",   Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_eleIDUp,    merged1g_eleIDUp,    resolved_eleIDUp))
                .Define("weight_eleIDDo",   Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_eleIDDo,    merged1g_eleIDDo,    resolved_eleIDDo))
                .Define("weight_phoIDUp",   Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_phoIDUp,    merged1g_phoIDUp,    resolved_phoIDUp))
                .Define("weight_phoIDDo",   Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_phoIDDo,    merged1g_phoIDDo,    resolved_phoIDDo))
                .Define("weight_csevUp",    Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_csevUp,     merged1g_csevUp,     resolved_csevUp))
                .Define("weight_csevDo",    Form("if (isM2) return %s; else if (isM1) return %s; else return %s;", merged2g_csevDo,     merged1g_csevDo,     resolved_csevDo));

    return nf;
}