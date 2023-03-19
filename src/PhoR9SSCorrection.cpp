#include <vector>
#include <mutex>
#include "PhoR9SSCorrection.h"
#include "TString.h"
#include "ROOT/RVec.hxx"
#include "TMVA/Reader.h"
#include "TFormula.h"

ROOT::RDF::RNode doSSCorrections(ROOT::RDF::RNode df, std::string corrColumn, std::string era, std::map<std::string, std::string> weightMap){

    // decalre MVA variables
    float phoEt_, phoPhi_, phoR9Full5x5_;
    float phoSCEtaWidth_, phoSCPhiWidth_, rho_;
    float phoSCEta_;
    float sieieFull5x5_, sieipFull5x5_, s4Full5x5_;

    // initialize TMVA Reader
    std::vector<std::shared_ptr<TMVA::Reader>> bdt = {
        std::make_shared<TMVA::Reader>("!Color:Silent"),
        std::make_shared<TMVA::Reader>("!Color:Silent")
    };

    // https://github.com/cms-analysis/flashgg/tree/dev_legacy_runII/MetaData/data/MetaConditions
    // /eos/cms/store/group/phys_higgs/cmshgg/flashgg-data/Taggers/data/PhoIdInputsCorrections
    for (int i = 0; i < 2; i++){
        std::string region = (i == 0) ? "EB" : "EE";
        bdt[i]->AddVariable("f0", &phoEt_);
        bdt[i]->AddVariable("f1", &phoSCEta_);
        bdt[i]->AddVariable("f2", &phoPhi_);
        bdt[i]->AddVariable("f3", &rho_);
        bdt[i]->AddVariable("f4", &sieipFull5x5_);
        bdt[i]->AddVariable("f5", &s4Full5x5_);
        bdt[i]->AddVariable("f6", &phoR9Full5x5_);
        bdt[i]->AddVariable("f7", &phoSCPhiWidth_);
        bdt[i]->AddVariable("f8", &sieieFull5x5_);
        bdt[i]->AddVariable("f9", &phoSCEtaWidth_);
        bdt[i]->BookMVA("BDT", weightMap[region].c_str());
    }

    // initialize correction formula
    std::vector<std::shared_ptr<TFormula>> fR9 = {
        std::make_shared<TFormula>(nullptr),
        std::make_shared<TFormula>(nullptr)
    };

    if (era == "UL2016preVFP"){
        fR9[0] = std::make_shared<TFormula>("", "x[0]*0.0054010312397019256+0.00034956078302361693");
        fR9[1] = std::make_shared<TFormula>("", "x[0]*0.007378665956671193-0.0010303157348496295");
    }
    else if (era == "UL2016postVFP"){
        fR9[0] = std::make_shared<TFormula>("", "x[0]*0.005245345997489048+0.0004184281464778561");
        fR9[1] = std::make_shared<TFormula>("", "x[0]*0.0077575458552175125-0.0012002614286774627");
    }
    else if (era == "UL2017") {
        fR9[0] = std::make_shared<TFormula>("", "x[0]*0.008903793308924074-0.0031396514172706835");
        fR9[1] = std::make_shared<TFormula>("", "x[0]*0.01370282259444966-0.007243889920611868");
    }
    else if (era == "UL2018") {
        fR9[0] = std::make_shared<TFormula>("", "x[0]*0.00974529926595541-0.0032455444688871404");
        fR9[1] = std::make_shared<TFormula>("", "x[0]*0.015158424850746033-0.008216661026002492");
    }
    else
        throw std::runtime_error(Form("No suitable formula to perform ss correction for era: %s", era.c_str()));

    // std::mutex to ensure thread safty
    auto pred = [&, predm = std::make_shared<std::mutex>(), bdt, fR9, era](
        // const ROOT::RVec<int>& nPho, 
        const int phoIdx1, // only perform the correction on the selected photon
        const ROOT::RVec<float>& phoEt,
        const ROOT::RVec<float>& phoPhi,
        const ROOT::RVec<float>& phoSCEta,
        const ROOT::RVec<float>& phoSCEtaWidth,
        const ROOT::RVec<float>& phoSCPhiWidth,
        const ROOT::RVec<float>& phoSigmaIEtaIEtaFull5x5,
        const ROOT::RVec<float>& phoSigmaIEtaIPhiFull5x5,
        const ROOT::RVec<float>& phoR9Full5x5,
        const ROOT::RVec<float>& phoE2x2Full5x5,
        const ROOT::RVec<float>& phoE5x5Full5x5,
        const float rho
    ){
        int iBE = (fabs(phoSCEta[phoIdx1]) < 1.479) ? 0 : 1;

        // set MVA variables
        phoEt_ = phoEt[phoIdx1];
        phoSCEta_ = phoSCEta[phoIdx1];
        phoPhi_ = phoPhi[phoIdx1];
        rho_ = rho;
        phoSCPhiWidth_ = phoSCPhiWidth[phoIdx1];
        sieipFull5x5_ = phoSigmaIEtaIPhiFull5x5[phoIdx1];
        s4Full5x5_ = phoE2x2Full5x5[phoIdx1] / phoE5x5Full5x5[phoIdx1];
        phoR9Full5x5_ = phoR9Full5x5[phoIdx1];
        sieieFull5x5_ = phoSigmaIEtaIEtaFull5x5[phoIdx1];
        phoSCEtaWidth_ = phoSCEtaWidth[phoIdx1];

        std::lock_guard<std::mutex> plock(*predm); // lock the mutex at construction, releases it at destruction
        float reg = bdt[iBE]->EvaluateRegression(0, "BDT");

        float corrR9 = phoR9Full5x5[phoIdx1] + fR9[iBE]->Eval(reg);

        return corrR9;

        // ROOT::RVec<float> corrPhoR9;
        // corrPhoR9.clear();
        // for (int i = 0; i < nPho; i++){
        //     int iBE = (fabs(phoSCEta[i]) < 1.479) ? 0 : 1;

        //     // set MVA variables
        //     phoEt_ = phoEt[i];
        //     phoSCEta_ = phoSCEta[i];
        //     phoPhi_ = phoPhi[i];
        //     rho_ = rho;
        //     phoSCPhiWidth_ = phoSCPhiWidth[i];
        //     sieipFull5x5_ = phoSigmaIEtaIPhiFull5x5[i];
        //     s4Full5x5_ = phoE2x2Full5x5[i] / phoE5x5Full5x5[i];
        //     phoR9Full5x5_ = phoR9Full5x5[i];
        //     sieieFull5x5_ = phoSigmaIEtaIEtaFull5x5[i];
        //     phoSCEtaWidth_ = phoSCEtaWidth[i];

        //     // lock the mutex at construction, releases it at destruction
        //     // evaluate the correction
        //     std::lock_guard<std::mutex> plock(*predm); 
        //     float reg = bdt[iBE]->EvaluateRegression(0, "BDT");
        //     float corrR9 = phoR9Full5x5[i] + fR9[iBE]->Eval(reg);

        //     corrPhoR9.push_back(corrR9);
        // }
        // return corrPhoR9;
    };

    auto nf = df.Define(corrColumn, pred, {
        "phoIdx1",
        "phoEt",
        "phoPhi",
        "phoSCEta",
        "phoSCEtaWidth",
        "phoSCPhiWidth",
        "phoSigmaIEtaIEtaFull5x5",
        "phoSigmaIEtaIPhiFull5x5",
        "phoR9Full5x5",
        "phoE2x2Full5x5",
        "phoE5x5Full5x5",
        "rho",
    });

    return nf;
}