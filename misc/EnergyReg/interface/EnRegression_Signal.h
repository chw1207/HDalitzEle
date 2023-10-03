#include <vector>
#include <mutex>
#include "TString.h"
#include "ROOT/RVec.hxx"
#include "TMVA/Reader.h"
#include "TFormula.h"

// R__LOAD_LIBRARY($HDalitzEle_LOC/lib/libHDalitzEle.so)
// R__ADD_INCLUDE_PATH($HDalitzEle_LOC/include)
// #include "XGBReader.h"



// ROOT::RDF::RNode doEnRegression_BDTG(ROOT::RDF::RNode df){
//     // decalre MVA variables
//     float rho_, nVtx_;
//     float eleSCEta_, eleSCPhi_;
//     float eleSCRawEn_, eleCalibPt_;
//     float eledEtaAtVtx_, eledPhiAtVtx_;
//     float elePtError_, eleHoverE_, eleEoverP_, eleEoverPout_, eleEoverPInv_;
//     float eleSCEtaWidth_, eleSCPhiWidth_, eleSigmaIEtaIEtaFull5x5_, eleSigmaIPhiIPhiFull5x5_;
//     float eleR9Full5x5_, eleBrem_;
//     float diTrkPt_, gsfDeltaR_;
//     float eleESEnToRawE_;

//     // initialize TMVA Reader
//     std::vector<std::string> region = {"EB", "EE"};
//     std::vector<std::shared_ptr<TMVA::Reader>> model = {
//         std::make_shared<TMVA::Reader>("!Color:Silent"),
//         std::make_shared<TMVA::Reader>("!Color:Silent")
//     };

//     std::vector<std::string> paths = {
//         "/data4/chenghan/external/RegressionV2/dataset_reg_EB_diGenPtNoShower/weights/TMVARegression_BDTG.weights.xml",
//         "/data4/chenghan/external/RegressionV2/dataset_reg_EE_diGenPtNoShower/weights/TMVARegression_BDTG.weights.xml"
//     };
//     for (size_t i = 0; i < region.size(); i++){
//         std::string path = paths[i];

//         model[i]->AddVariable("rho",                          &rho_);
//         model[i]->AddVariable("nVtx",                         &nVtx_);
//         model[i]->AddVariable("eleSCEta_Lead",                &eleSCEta_);
//         model[i]->AddVariable("eleSCPhi_Lead",                &eleSCPhi_);
//         model[i]->AddVariable("eleSCRawEn_Lead",              &eleSCRawEn_);
//         model[i]->AddVariable("eleCalibPt_Lead",              &eleCalibPt_);   
//         model[i]->AddVariable("eledEtaAtVtx_Lead",            &eledEtaAtVtx_);
//         model[i]->AddVariable("eledPhiAtVtx_Lead",            &eledPhiAtVtx_);
//         model[i]->AddVariable("elePtError_Lead",              &elePtError_);
//         model[i]->AddVariable("eleHoverE_Lead",               &eleHoverE_);
//         model[i]->AddVariable("eleEoverP_Lead",               &eleEoverP_);
//         model[i]->AddVariable("eleEoverPout_Lead",            &eleEoverPout_);
//         model[i]->AddVariable("eleEoverPInv_Lead",            &eleEoverPInv_);
//         model[i]->AddVariable("eleSCEtaWidth_Lead",           &eleSCEtaWidth_);
//         model[i]->AddVariable("eleSCPhiWidth_Lead",           &eleSCPhiWidth_);
//         model[i]->AddVariable("eleSigmaIEtaIEtaFull5x5_Lead", &eleSigmaIEtaIEtaFull5x5_);
//         model[i]->AddVariable("eleSigmaIPhiIPhiFull5x5_Lead", &eleSigmaIPhiIPhiFull5x5_);
//         model[i]->AddVariable("eleR9Full5x5_Lead",            &eleR9Full5x5_);
//         model[i]->AddVariable("eleBrem_Lead",                 &eleBrem_);
//         model[i]->AddVariable("diTrk.Pt()",                   &diTrkPt_);
//         model[i]->AddVariable("gsfDeltaR_Lead",               &gsfDeltaR_);
//         if (region[i] == "EE")
//             model[i]->AddVariable("eleESEnToRawE_Lead",       &eleESEnToRawE_);

//         model[i]->BookMVA("BDTG", path.c_str());
//     }

//     auto regression = [&, regm = std::make_shared<std::mutex>(), model](
//         const float rho,
//         const int nVtx,
//         const float eleSCEta,
//         const float eleSCPhi,
//         const float eleSCRawEn,
//         const float eleCalibPt,
//         const float eledEtaAtVtx,
//         const float eledPhiAtVtx,
//         const float elePtError,
//         const float eleHoverE,
//         const float eleEoverP,
//         const float eleEoverPout,
//         const float eleEoverPInv,
//         const float eleSCEtaWidth,
//         const float eleSCPhiWidth,
//         const float eleSigmaIEtaIEtaFull5x5,
//         const float eleSigmaIPhiIPhiFull5x5,
//         const float eleR9Full5x5,
//         const float eleBrem,
//         const float diTrkPt,
//         const double gsfDeltaR,
//         const float eleESEnToRawE
//     ){
//         int iBE = (fabs(eleSCEta) < 1.479) ? 0 : 1;

//         rho_ = rho;
//         nVtx_ = (float) nVtx;
//         eleSCEta_ = eleSCEta;
//         eleSCPhi_ = eleSCPhi;
//         eleSCRawEn_ =  eleSCRawEn;
//         eleCalibPt_ =  eleCalibPt;
//         eledEtaAtVtx_ = eledEtaAtVtx;
//         eledPhiAtVtx_ = eledPhiAtVtx;
//         elePtError_ = elePtError;
//         eleHoverE_ = eleHoverE;
//         eleEoverP_ = eleEoverP;
//         eleEoverPout_ = eleEoverPout;
//         eleEoverPInv_ = eleEoverPInv;
//         eleSCEtaWidth_ = eleSCEtaWidth;
//         eleSCPhiWidth_ = eleSCPhiWidth;
//         eleSigmaIEtaIEtaFull5x5_ = eleSigmaIEtaIEtaFull5x5;
//         eleSigmaIPhiIPhiFull5x5_ = eleSigmaIPhiIPhiFull5x5;
//         eleR9Full5x5_ = eleR9Full5x5;
//         eleBrem_ = eleBrem;
//         diTrkPt_ = diTrkPt;
//         gsfDeltaR_ = (float) gsfDeltaR;
//         eleESEnToRawE_ = eleESEnToRawE;

//         // lock the mutex at construction, releases it at destruction
//         // evaluate the correction
//         std::lock_guard<std::mutex> rlock(*regm);
//         float reg = model[iBE]->EvaluateRegression(0, "BDTG");
//         float corrEt = eleCalibPt * reg;

//         return corrEt;
//     };

//     auto nf = df.Define("eleHDALRegPt_Lead", regression, {
//         "rho",
//         "nVtx",
//         "eleSCEta_Lead",
//         "eleSCPhi_Lead",
//         "eleSCRawEn_Lead",
//         "eleCalibPt_Lead",
//         "eledEtaAtVtx_Lead",
//         "eledPhiAtVtx_Lead",
//         "elePtError_Lead",              
//         "eleHoverE_Lead",               
//         "eleEoverP_Lead",               
//         "eleEoverPout_Lead",            
//         "eleEoverPInv_Lead",            
//         "eleSCEtaWidth_Lead",           
//         "eleSCPhiWidth_Lead",           
//         "eleSigmaIEtaIEtaFull5x5_Lead", 
//         "eleSigmaIPhiIPhiFull5x5_Lead", 
//         "eleR9Full5x5_Lead",            
//         "eleBrem_Lead",                
//         "diTrkPt_Lead",
//         "gsfDeltaR_Lead",
//         "eleESEnToRawE_Lead"         
//     });

//     return nf;
// }


// ROOT::RDF::RNode doEnRegression_BDTG_Trk(ROOT::RDF::RNode df){
//     // decalre MVA variables
//     float nVtx_;
//     float eleSCRawEn_, eleCalibPt_;
//     float eledEtaAtVtx_, eledPhiAtVtx_;
//     float elePtError_, eleHoverE_, eleEoverP_, eleEoverPout_, eleEoverPInv_;
//     float eleSCEtaWidth_, eleSCPhiWidth_, eleSigmaIEtaIEtaFull5x5_, eleSigmaIPhiIPhiFull5x5_;
//     float eleR9Full5x5_, eleBrem_;
//     float diTrkPt_, gsfDeltaR_;
//     float eleESEnToRawE_;

//     // initialize TMVA Reader
//     std::vector<std::string> region = {"EB", "EE"};
//     std::vector<std::shared_ptr<TMVA::Reader>> model = {
//         std::make_shared<TMVA::Reader>("!Color:Silent"),
//         std::make_shared<TMVA::Reader>("!Color:Silent")
//     };
//     for (size_t i = 0; i < region.size(); i++){
//         std::string path = Form("/data4/chenghan/external/Regression/dataset_reg_%s_trk/weights/TMVARegression_BDTG.weights.xml", region[i].c_str());

//         model[i]->AddVariable("nVtx",                         &nVtx_);
//         model[i]->AddVariable("eleSCRawEn_Lead",              &eleSCRawEn_);
//         model[i]->AddVariable("eleCalibPt_Lead",              &eleCalibPt_);   
//         model[i]->AddVariable("eledEtaAtVtx_Lead",            &eledEtaAtVtx_);
//         model[i]->AddVariable("eledPhiAtVtx_Lead",            &eledPhiAtVtx_);
//         model[i]->AddVariable("elePtError_Lead",              &elePtError_);
//         model[i]->AddVariable("eleHoverE_Lead",               &eleHoverE_);
//         model[i]->AddVariable("eleEoverP_Lead",               &eleEoverP_);
//         model[i]->AddVariable("eleEoverPout_Lead",            &eleEoverPout_);
//         model[i]->AddVariable("eleEoverPInv_Lead",            &eleEoverPInv_);
//         model[i]->AddVariable("eleSCEtaWidth_Lead",           &eleSCEtaWidth_);
//         model[i]->AddVariable("eleSCPhiWidth_Lead",           &eleSCPhiWidth_);
//         model[i]->AddVariable("eleSigmaIEtaIEtaFull5x5_Lead", &eleSigmaIEtaIEtaFull5x5_);
//         model[i]->AddVariable("eleSigmaIPhiIPhiFull5x5_Lead", &eleSigmaIPhiIPhiFull5x5_);
//         model[i]->AddVariable("eleR9Full5x5_Lead",            &eleR9Full5x5_);
//         model[i]->AddVariable("eleBrem_Lead",                 &eleBrem_);
//         model[i]->AddVariable("diTrk.Pt()",                   &diTrkPt_);
//         model[i]->AddVariable("gsfDeltaR_Lead",               &gsfDeltaR_);
//         if (region[i] == "EE")
//             model[i]->AddVariable("(eleESEnP1_Lead + eleESEnP2_Lead)/eleSCRawEn_Lead", &eleESEnToRawE_);

//         model[i]->BookMVA("BDTG", path.c_str());
//     }

//     auto regression = [&, regm = std::make_shared<std::mutex>(), model](
//         const float eleSCEta,
//         const int nVtx,
//         const float eleSCRawEn,
//         const float eleCalibPt,
//         const float elePt,
//         const float eledEtaAtVtx,
//         const float eledPhiAtVtx,
//         const float elePtError,
//         const float eleHoverE,
//         const float eleEoverP,
//         const float eleEoverPout,
//         const float eleEoverPInv,
//         const float eleSCEtaWidth,
//         const float eleSCPhiWidth,
//         const float eleSigmaIEtaIEtaFull5x5,
//         const float eleSigmaIPhiIPhiFull5x5,
//         const float eleR9Full5x5,
//         const float eleBrem,
//         const float diTrkPt,
//         const double gsfDeltaR,
//         const float eleESEnToRawE
//     ){
//         int iBE = (fabs(eleSCEta) < 1.479) ? 0 : 1;

//         nVtx_ = (float) nVtx;
//         eleSCRawEn_ =  eleSCRawEn;
//         eleCalibPt_ =  eleCalibPt;
//         eledEtaAtVtx_ = eledEtaAtVtx;
//         eledPhiAtVtx_ = eledPhiAtVtx;
//         elePtError_ = elePtError;
//         eleHoverE_ = eleHoverE;
//         eleEoverP_ = eleEoverP;
//         eleEoverPout_ = eleEoverPout;
//         eleEoverPInv_ = eleEoverPInv;
//         eleSCEtaWidth_ = eleSCEtaWidth;
//         eleSCPhiWidth_ = eleSCPhiWidth;
//         eleSigmaIEtaIEtaFull5x5_ = eleSigmaIEtaIEtaFull5x5;
//         eleSigmaIPhiIPhiFull5x5_ = eleSigmaIPhiIPhiFull5x5;
//         eleR9Full5x5_ = eleR9Full5x5;
//         eleBrem_ = eleBrem;
//         diTrkPt_ = diTrkPt;
//         gsfDeltaR_ = (float) gsfDeltaR;
//         eleESEnToRawE_ = eleESEnToRawE;

//         // lock the mutex at construction, releases it at destruction
//         // evaluate the correction
//         std::lock_guard<std::mutex> rlock(*regm);
//         float reg = model[iBE]->EvaluateRegression(0, "BDTG");
//         float corrEt = diTrkPt * reg;

//         return corrEt;
//     };

//     auto nf = df.Define("eleHDALRegPt_Lead", regression, {
//         "eleSCEta_Lead",

//         "nVtx",
//         "eleSCRawEn_Lead",
//         "eleCalibPt_Lead",
//         "eledEtaAtVtx_Lead",
//         "eledPhiAtVtx_Lead",
//         "elePtError_Lead",              
//         "eleHoverE_Lead",               
//         "eleEoverP_Lead",               
//         "eleEoverPout_Lead",            
//         "eleEoverPInv_Lead",            
//         "eleSCEtaWidth_Lead",           
//         "eleSCPhiWidth_Lead",           
//         "eleSigmaIEtaIEtaFull5x5_Lead", 
//         "eleSigmaIPhiIPhiFull5x5_Lead", 
//         "eleR9Full5x5_Lead",            
//         "eleBrem_Lead",                
//         "diTrkPt_Lead",
//         "gsfDeltaR_Lead",
//         "eleESEnToRawE_Lead"         
//     });

//     return nf;
// }


// ROOT::RDF::RNode doEnRegression_XGB(ROOT::RDF::RNode df){
//     std::vector<std::shared_ptr<XGBReader>> model = {
//         std::make_shared<XGBReader>("../models/MergedEn_Regression_XGB_EB_Pt25To150GeV_StdDev.txt"),
//         std::make_shared<XGBReader>("../models/MergedEn_Regression_XGB_EE_Pt25To150GeV_StdDev.txt")
//     };
    
//     auto regression = [&, model](
//         const float rho,
//         const int nVtx,
//         const float eleSCEta,
//         const float eleSCPhi,
//         const float eleSCRawEn,
//         const float eleCalibPt,
//         // const float elePt,
//         const float eledEtaAtVtx,
//         const float eledPhiAtVtx,
//         const float elePtError,
//         const float eleHoverE,
//         const float eleEoverP,
//         const float eleEoverPout,
//         const float eleEoverPInv,
//         const float eleSCEtaWidth,
//         const float eleSCPhiWidth,
//         const float eleSigmaIEtaIEtaFull5x5,
//         const float eleSigmaIPhiIPhiFull5x5,
//         const float eleR9Full5x5,
//         const float eleBrem,
//         const double gsfPtSum,
//         const double gsfPtRatio,
//         const float diTrkPt,
//         const double gsfDeltaR,
//         // const float eleEmax,
//         // const float eleE2nd,
//         // const float RatioRL,
//         // const float RatioBT,
//         const float eleESEnToRawE
//     ){
//         int iBE = (fabs(eleSCEta) < 1.479) ? 0 : 1;

//         std::vector<float> features = {
//             rho,
//             (float) nVtx,
//             eleSCEta,
//             eleSCPhi,
//             eleSCRawEn,
//             eleCalibPt,
//             // elePt,
//             eledEtaAtVtx,
//             eledPhiAtVtx,
//             elePtError,
//             eleHoverE,
//             eleEoverP,
//             eleEoverPout,
//             eleEoverPInv,
//             eleSCEtaWidth,
//             eleSCPhiWidth,
//             eleSigmaIEtaIEtaFull5x5,
//             eleSigmaIPhiIPhiFull5x5,
//             eleR9Full5x5,
//             eleBrem,
//             (float) gsfPtSum,
//             (float) gsfPtRatio,
//             diTrkPt,
//             (float) gsfDeltaR
//         };
//         if (iBE == 1)
//             features.push_back(eleESEnToRawE);
        
//         float reg = model[iBE]->Compute(features)[0];
//         // float corrEt = elePt * reg;
//         float corrEt = eleCalibPt * reg;

//         return corrEt;
//     };

//     auto nf = df.Define("eleHDALRegPt_Lead", regression, {
//         "rho",
//         "nVtx",
//         "eleSCEta_Lead",
//         "eleSCPhi_Lead",
//         "eleSCRawEn_Lead",
//         // "elePt_Lead",
//         "eleCalibPt_Lead",
//         "eledEtaAtVtx_Lead",
//         "eledPhiAtVtx_Lead",
//         "elePtError_Lead",              
//         "eleHoverE_Lead",               
//         "eleEoverP_Lead",               
//         "eleEoverPout_Lead",            
//         "eleEoverPInv_Lead",            
//         "eleSCEtaWidth_Lead",           
//         "eleSCPhiWidth_Lead",           
//         "eleSigmaIEtaIEtaFull5x5_Lead", 
//         "eleSigmaIPhiIPhiFull5x5_Lead", 
//         "eleR9Full5x5_Lead",            
//         "eleBrem_Lead", 
//         "gsfPtSum_Lead",
//         "gsfPtRatio_Lead",               
//         "diTrkPt_Lead",
//         "gsfDeltaR_Lead",
//         // "eleEmax_Lead",
//         // "eleE2nd_Lead",
//         // "RatioRL",
//         // "RatioBT",
//         "eleESEnToRawE_Lead"         
//     });

//     return nf;
// }


ROOT::RDF::RNode doEnRegression_BDTG(ROOT::RDF::RNode df){
    // decalre MVA variables
    float rho_, nVtx_;
    float eleSCEta_, eleSCPhi_;
    float eleSCRawEn_, eleCalibPt_;
    float eledEtaAtVtx_, eledPhiAtVtx_;
    float elePtError_, eleHoverE_, eleEoverP_, eleEoverPout_, eleEoverPInv_;
    float eleSCEtaWidth_, eleSCPhiWidth_, eleSigmaIEtaIEtaFull5x5_, eleSigmaIPhiIPhiFull5x5_;
    float eleR9Full5x5_, eleBrem_;
    float diTrkPt_, gsfDeltaR_;
    float eleESEnToRawE_;

    // initialize TMVA Reader
    std::vector<std::string> region = {"EB", "EE"};
    std::vector<std::shared_ptr<TMVA::Reader>> model = {
        std::make_shared<TMVA::Reader>("!Color:Silent"),
        std::make_shared<TMVA::Reader>("!Color:Silent")
    };

    std::vector<std::string> paths = {
        "/data4/chenghan/external/RegressionV2/dataset_reg_EB_diGenPtNoShower/weights/TMVARegression_BDTG.weights.xml",
        "/data4/chenghan/external/RegressionV2/dataset_reg_EE_diGenPtNoShower/weights/TMVARegression_BDTG.weights.xml"
    };
    for (size_t i = 0; i < region.size(); i++){
        std::string path = paths[i];

        model[i]->AddVariable("rho",                          &rho_);
        model[i]->AddVariable("nVtx",                         &nVtx_);
        model[i]->AddVariable("eleSCEta_Lead",                &eleSCEta_);
        model[i]->AddVariable("eleSCPhi_Lead",                &eleSCPhi_);
        model[i]->AddVariable("eleSCRawEn_Lead",              &eleSCRawEn_);
        model[i]->AddVariable("eleCalibPt_Lead",              &eleCalibPt_);   
        model[i]->AddVariable("eledEtaAtVtx_Lead",            &eledEtaAtVtx_);
        model[i]->AddVariable("eledPhiAtVtx_Lead",            &eledPhiAtVtx_);
        model[i]->AddVariable("elePtError_Lead",              &elePtError_);
        model[i]->AddVariable("eleHoverE_Lead",               &eleHoverE_);
        model[i]->AddVariable("eleEoverP_Lead",               &eleEoverP_);
        model[i]->AddVariable("eleEoverPout_Lead",            &eleEoverPout_);
        model[i]->AddVariable("eleEoverPInv_Lead",            &eleEoverPInv_);
        model[i]->AddVariable("eleSCEtaWidth_Lead",           &eleSCEtaWidth_);
        model[i]->AddVariable("eleSCPhiWidth_Lead",           &eleSCPhiWidth_);
        model[i]->AddVariable("eleSigmaIEtaIEtaFull5x5_Lead", &eleSigmaIEtaIEtaFull5x5_);
        model[i]->AddVariable("eleSigmaIPhiIPhiFull5x5_Lead", &eleSigmaIPhiIPhiFull5x5_);
        model[i]->AddVariable("eleR9Full5x5_Lead",            &eleR9Full5x5_);
        model[i]->AddVariable("eleBrem_Lead",                 &eleBrem_);
        model[i]->AddVariable("diTrk.Pt()",                   &diTrkPt_);
        model[i]->AddVariable("gsfDeltaR_Lead",               &gsfDeltaR_);
        if (region[i] == "EE")
            model[i]->AddVariable("eleESEnToRawE_Lead",       &eleESEnToRawE_);

        model[i]->BookMVA("BDTG", path.c_str());
    }

    auto regression = [&, regm = std::make_shared<std::mutex>(), model](
        const float rho,
        const int nVtx,
        const float eleSCEta,
        const float eleSCPhi,
        const float eleSCRawEn,
        const float eleCalibPt,
        const float eledEtaAtVtx,
        const float eledPhiAtVtx,
        const float elePtError,
        const float eleHoverE,
        const float eleEoverP,
        const float eleEoverPout,
        const float eleEoverPInv,
        const float eleSCEtaWidth,
        const float eleSCPhiWidth,
        const float eleSigmaIEtaIEtaFull5x5,
        const float eleSigmaIPhiIPhiFull5x5,
        const float eleR9Full5x5,
        const float eleBrem,
        const float diTrkPt,
        const double gsfDeltaR,
        const float eleESEnToRawE
    ){
        int iBE = (fabs(eleSCEta) < 1.479) ? 0 : 1;

        rho_ = rho;
        nVtx_ = (float) nVtx;
        eleSCEta_ = eleSCEta;
        eleSCPhi_ = eleSCPhi;
        eleSCRawEn_ =  eleSCRawEn;
        eleCalibPt_ =  eleCalibPt;
        eledEtaAtVtx_ = eledEtaAtVtx;
        eledPhiAtVtx_ = eledPhiAtVtx;
        elePtError_ = elePtError;
        eleHoverE_ = eleHoverE;
        eleEoverP_ = eleEoverP;
        eleEoverPout_ = eleEoverPout;
        eleEoverPInv_ = eleEoverPInv;
        eleSCEtaWidth_ = eleSCEtaWidth;
        eleSCPhiWidth_ = eleSCPhiWidth;
        eleSigmaIEtaIEtaFull5x5_ = eleSigmaIEtaIEtaFull5x5;
        eleSigmaIPhiIPhiFull5x5_ = eleSigmaIPhiIPhiFull5x5;
        eleR9Full5x5_ = eleR9Full5x5;
        eleBrem_ = eleBrem;
        diTrkPt_ = diTrkPt;
        gsfDeltaR_ = (float) gsfDeltaR;
        eleESEnToRawE_ = eleESEnToRawE;

        // lock the mutex at construction, releases it at destruction
        // evaluate the correction
        std::lock_guard<std::mutex> rlock(*regm);
        float reg = model[iBE]->EvaluateRegression(0, "BDTG");
        float corrEt = eleCalibPt * reg;

        return corrEt;
    };

    auto nf = df.Define("eleHDALRegPt_Lead", regression, {
        "rho",
        "nVtx",
        "eleSCEta_Lead",
        "eleSCPhi_Lead",
        "eleSCRawEn_Lead",
        "eleCalibPt_Lead",
        "eledEtaAtVtx_Lead",
        "eledPhiAtVtx_Lead",
        "elePtError_Lead",              
        "eleHoverE_Lead",               
        "eleEoverP_Lead",               
        "eleEoverPout_Lead",            
        "eleEoverPInv_Lead",            
        "eleSCEtaWidth_Lead",           
        "eleSCPhiWidth_Lead",           
        "eleSigmaIEtaIEtaFull5x5_Lead", 
        "eleSigmaIPhiIPhiFull5x5_Lead", 
        "eleR9Full5x5_Lead",            
        "eleBrem_Lead",                
        "diTrkPt_Lead",
        "gsfDeltaR_Lead",
        "eleESEnToRawE_Lead"         
    });

    return nf;
}


float effSigma(TH1D* hist, double quantile = TMath::Erf(1.0/sqrt(2.0))){
    TAxis* xaxis = hist->GetXaxis();
    int nb = xaxis->GetNbins();
    if(nb < 10) {
        cout << "effsigma: Not a valid histo(too few). nbins = " << nb << endl;
        return 0.;
    }

    float bwid = xaxis->GetBinWidth(1);
    if(bwid == 0) {
        cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
        return 0.;
    }

    float xmax = xaxis->GetXmax();
    float xmin = xaxis->GetXmin();
    float ave = hist->GetMean();
    // float ave = hist->GetBinCenter(hist->GetMaximumBin()); // modified by CH 
    float rms = hist->GetRMS();
    float total = 0.;
    for(int i = 0; i < nb+2; i++) {
        total += hist->GetBinContent(i);
    }
    int ierr = 0;
    int ismin = 999;

    // Scan around window of mean: window RMS/binWidth (cannot be bigger than 0.1*number of bins)
    // Determine minimum width of distribution which holds 0.693 of total
    float rlim = quantile*total;
    int nrms = rms/bwid; // Set scan size to +/- rms
    if(nrms > nb/10) 
        nrms = nb/10; // Could be tuned...

    float widmin = 9999999.;
    for(int iscan = -nrms; iscan < nrms+1; iscan++){ // Scan window centre
        int ibm = (ave - xmin)/bwid + 1 + iscan; // Find bin idx in scan: iscan from mean
        float x = (ibm - 0.5)*bwid + xmin; // 0.5 for bin centre
        float xj = x;
        float xk = x;
        int jbm = ibm;
        int kbm = ibm;

        // Define counter for yield in bins: stop when counter > rlim
        float bin = hist->GetBinContent(ibm);
        total = bin;
        for(int j = 1; j < nb; j++){
            if(jbm < nb) {
                jbm++;
                xj += bwid;
                bin = hist->GetBinContent(jbm);
                total += bin;
                if(total > rlim) 
                    break;
            }
            else 
                ierr = 1;

            if(kbm > 0) {
                kbm--;
                xk -= bwid;
                bin = hist->GetBinContent(kbm);
                total += bin;
                if(total > rlim) 
                    break;
            }
            else 
                ierr = 2;
        }

        // Calculate fractional width in bin takes above limt (assume linear)
        float dxf = (total-rlim)*bwid/bin;
        float wid = (xj-xk+bwid-dxf)*0.5; // Total width: half of peak

        if(wid < widmin) {
            widmin = wid;
            ismin = iscan;
        }
    }

    if(ismin == nrms || ismin == -nrms) 
        ierr=3;
    if(ierr != 0) 
        cout << "effsigma: Error of type " << ierr << endl;
    return widmin;
}