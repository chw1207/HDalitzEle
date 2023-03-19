#include "MergedEnCorrector.h"
#include <mutex>
#include "ROOT/RVec.hxx"
#include "TMVA/Reader.h"
#include "XGBReader.h"

ROOT::RDF::RNode doHDALRegression(ROOT::RDF::RNode df, std::vector<std::string> path){
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
    for (size_t i = 0; i < region.size(); i++){
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

        model[i]->BookMVA("BDTG", path[i].c_str());
    }

    auto regression = [&, regm = std::make_shared<std::mutex>(), model](
        const int nEle,
        const float rho,
        const int nVtx,
        const ROOT::RVec<float>& eleSCEta,
        const ROOT::RVec<float>& eleSCPhi,
        const ROOT::RVec<float>& eleSCRawEn,
        const ROOT::RVec<float>& eleCalibPt,
        const ROOT::RVec<float>& eledEtaAtVtx,
        const ROOT::RVec<float>& eledPhiAtVtx,
        const ROOT::RVec<float>& elePtError,
        const ROOT::RVec<float>& eleHoverE,
        const ROOT::RVec<float>& eleEoverP,
        const ROOT::RVec<float>& eleEoverPout,
        const ROOT::RVec<float>& eleEoverPInv,
        const ROOT::RVec<float>& eleSCEtaWidth,
        const ROOT::RVec<float>& eleSCPhiWidth,
        const ROOT::RVec<float>& eleSigmaIEtaIEtaFull5x5,
        const ROOT::RVec<float>& eleSigmaIPhiIPhiFull5x5,
        const ROOT::RVec<float>& eleR9Full5x5,
        const ROOT::RVec<float>& eleBrem,
        const ROOT::RVec<float>& diTrkPt,
        const ROOT::RVec<float>& gsfDeltaR,
        const ROOT::RVec<float>& eleESEnToRawE
    ){
        ROOT::RVec<float> corrPt(nEle);
        for (int i = 0; i < nEle; i++){
            int iBE = (fabs(eleSCEta[i]) < 1.479) ? 0 : 1;

            rho_ = rho;
            nVtx_ = (float) nVtx;
            eleSCEta_ = eleSCEta[i];
            eleSCPhi_ = eleSCPhi[i];
            eleSCRawEn_ =  eleSCRawEn[i];
            eleCalibPt_ =  eleCalibPt[i];
            eledEtaAtVtx_ = eledEtaAtVtx[i];
            eledPhiAtVtx_ = eledPhiAtVtx[i];
            elePtError_ = elePtError[i];
            eleHoverE_ = eleHoverE[i];
            eleEoverP_ = eleEoverP[i];
            eleEoverPout_ = eleEoverPout[i];
            eleEoverPInv_ = eleEoverPInv[i];
            eleSCEtaWidth_ = eleSCEtaWidth[i];
            eleSCPhiWidth_ = eleSCPhiWidth[i];
            eleSigmaIEtaIEtaFull5x5_ = eleSigmaIEtaIEtaFull5x5[i];
            eleSigmaIPhiIPhiFull5x5_ = eleSigmaIPhiIPhiFull5x5[i];
            eleR9Full5x5_ = eleR9Full5x5[i];
            eleBrem_ = eleBrem[i];
            diTrkPt_ = diTrkPt[i];
            gsfDeltaR_ = gsfDeltaR[i];
            eleESEnToRawE_ = eleESEnToRawE[i];

            // lock the mutex at construction, releases it at destruction
            // evaluate the correction
            std::lock_guard<std::mutex> rlock(*regm);
            corrPt[i] = eleCalibPt[i] * model[iBE]->EvaluateRegression(0, "BDTG");
        }
        return corrPt;
    };

    auto nf = df.Define("eleHDALRegPt", regression, {
        "nEle",
        "rho",
        "nVtx",
        "eleSCEta",
        "eleSCPhi",
        "eleSCRawEn",
        "eleCalibPt",
        "eledEtaAtVtx",
        "eledPhiAtVtx",
        "elePtError",              
        "eleHoverE",               
        "eleEoverP",               
        "eleEoverPout",            
        "eleEoverPInv",            
        "eleSCEtaWidth",           
        "eleSCPhiWidth",           
        "eleSigmaIEtaIEtaFull5x5", 
        "eleSigmaIPhiIPhiFull5x5", 
        "eleR9Full5x5",            
        "eleBrem",                
        "diTrkPt",
        "gsfDeltaR",
        "eleESEnToRawE"         
    });

    return nf;
}


ROOT::RDF::RNode doHDALRegressionXGB(ROOT::RDF::RNode df, std::vector<std::string> path){
    std::vector<std::shared_ptr<XGBReader>> model = {
        std::make_shared<XGBReader>(path[0]),
        std::make_shared<XGBReader>(path[1])
    };

    auto regression = [&, model](
        const int nEle,
        const float rho,
        const int nVtx,
        const ROOT::RVec<float>& eleSCEta,
        const ROOT::RVec<float>& eleSCPhi,
        const ROOT::RVec<float>& eleSCRawEn,
        const ROOT::RVec<float>& eleCalibPt,
        const ROOT::RVec<float>& eledEtaAtVtx,
        const ROOT::RVec<float>& eledPhiAtVtx,
        const ROOT::RVec<float>& elePtError,
        const ROOT::RVec<float>& eleHoverE,
        const ROOT::RVec<float>& eleEoverP,
        const ROOT::RVec<float>& eleEoverPout,
        const ROOT::RVec<float>& eleEoverPInv,
        const ROOT::RVec<float>& eleSCEtaWidth,
        const ROOT::RVec<float>& eleSCPhiWidth,
        const ROOT::RVec<float>& eleSigmaIEtaIEtaFull5x5,
        const ROOT::RVec<float>& eleSigmaIPhiIPhiFull5x5,
        const ROOT::RVec<float>& eleR9Full5x5,
        const ROOT::RVec<float>& eleBrem,
        const ROOT::RVec<float>& gsfPtSum,
        const ROOT::RVec<float>& gsfPtRatio,
        const ROOT::RVec<float>& diTrkPt,
        const ROOT::RVec<float>& gsfDeltaR,
        const ROOT::RVec<float>& eleESEnToRawE
    ){
        ROOT::RVec<float> corrPt(nEle);
        for (int i = 0; i < nEle; i++){
            int iBE = (fabs(eleSCEta[i]) < 1.479) ? 0 : 1;

            std::vector<float> features = {
                rho,
                (float) nVtx,
                eleSCEta[i],
                eleSCPhi[i],
                eleSCRawEn[i],
                eleCalibPt[i],
                eledEtaAtVtx[i],
                eledPhiAtVtx[i],
                elePtError[i],
                eleHoverE[i],
                eleEoverP[i],
                eleEoverPout[i],
                eleEoverPInv[i],
                eleSCEtaWidth[i],
                eleSCPhiWidth[i],
                eleSigmaIEtaIEtaFull5x5[i],
                eleSigmaIPhiIPhiFull5x5[i],
                eleR9Full5x5[i],
                eleBrem[i],
                gsfPtSum[i],
                gsfPtRatio[i],
                diTrkPt[i],
                gsfDeltaR[i]
            };
            if (iBE == 1)
                features.push_back(eleESEnToRawE[i]);
        
            corrPt[i] = eleCalibPt[i] * model[iBE]->Compute(features)[0];
        }
        return corrPt;
    };

    auto nf = df.Define("eleHDALRegPt", regression, {
        "nEle",
        "rho",
        "nVtx",
        "eleSCEta",
        "eleSCPhi",
        "eleSCRawEn",
        "eleCalibPt",
        "eledEtaAtVtx",
        "eledPhiAtVtx",
        "elePtError",              
        "eleHoverE",               
        "eleEoverP",               
        "eleEoverPout",            
        "eleEoverPInv",            
        "eleSCEtaWidth",           
        "eleSCPhiWidth",           
        "eleSigmaIEtaIEtaFull5x5", 
        "eleSigmaIPhiIPhiFull5x5", 
        "eleR9Full5x5",            
        "eleBrem",   
        "gsfPtSum",
        "gsfPtRatio",              
        "diTrkPt",
        "gsfDeltaR",
        "eleESEnToRawE"         
    });

    return nf;
}


// ROOT::RDF::RNode doHDALRegressionV2(ROOT::RDF::RNode df, std::vector<std::string> path){
    
//     int nslot = df.GetNSlots();

//     // decalre MVA variables for each slot
//     std::vector<float> rho_(nslot);
//     std::vector<float> nVtx_(nslot);
//     std::vector<float> eleSCEta_(nslot);
//     std::vector<float> eleSCPhi_(nslot);
//     std::vector<float> eleSCRawEn_(nslot);
//     std::vector<float> eleCalibPt_(nslot);
//     std::vector<float> eledEtaAtVtx_(nslot);
//     std::vector<float> eledPhiAtVtx_(nslot);
//     std::vector<float> elePtError_(nslot);
//     std::vector<float> eleHoverE_(nslot);
//     std::vector<float> eleEoverP_(nslot);
//     std::vector<float> eleEoverPout_(nslot);
//     std::vector<float> eleEoverPInv_(nslot);
//     std::vector<float> eleSCEtaWidth_(nslot);
//     std::vector<float> eleSCPhiWidth_(nslot);
//     std::vector<float> eleSigmaIEtaIEtaFull5x5_(nslot);
//     std::vector<float> eleSigmaIPhiIPhiFull5x5_(nslot);
//     std::vector<float> eleR9Full5x5_(nslot);
//     std::vector<float> eleBrem_(nslot);
//     std::vector<float> diTrkPt_(nslot);
//     std::vector<float> gsfDeltaR_(nslot);
//     std::vector<float> eleESEnToRawE_(nslot);

//     // initialize TMVA Reader
//     std::vector<std::vector<TMVA::Reader*>> Readers(2, std::vector<TMVA::Reader*>(nslot));
//     for (int i = 0 ; i < 2; i++){
//         for (int slot = 0; slot < nslot; slot++){
//             Readers[i][slot] = new TMVA::Reader("!Color:Silent");
//             Readers[i][slot]->AddVariable("rho",                          &rho_[slot]);
//             Readers[i][slot]->AddVariable("nVtx",                         &nVtx_[slot]);
//             Readers[i][slot]->AddVariable("eleSCEta_Lead",                &eleSCEta_[slot]);
//             Readers[i][slot]->AddVariable("eleSCPhi_Lead",                &eleSCPhi_[slot]);
//             Readers[i][slot]->AddVariable("eleSCRawEn_Lead",              &eleSCRawEn_[slot]);
//             Readers[i][slot]->AddVariable("eleCalibPt_Lead",              &eleCalibPt_[slot]);   
//             Readers[i][slot]->AddVariable("eledEtaAtVtx_Lead",            &eledEtaAtVtx_[slot]);
//             Readers[i][slot]->AddVariable("eledPhiAtVtx_Lead",            &eledPhiAtVtx_[slot]);
//             Readers[i][slot]->AddVariable("elePtError_Lead",              &elePtError_[slot]);
//             Readers[i][slot]->AddVariable("eleHoverE_Lead",               &eleHoverE_[slot]);
//             Readers[i][slot]->AddVariable("eleEoverP_Lead",               &eleEoverP_[slot]);
//             Readers[i][slot]->AddVariable("eleEoverPout_Lead",            &eleEoverPout_[slot]);
//             Readers[i][slot]->AddVariable("eleEoverPInv_Lead",            &eleEoverPInv_[slot]);
//             Readers[i][slot]->AddVariable("eleSCEtaWidth_Lead",           &eleSCEtaWidth_[slot]);
//             Readers[i][slot]->AddVariable("eleSCPhiWidth_Lead",           &eleSCPhiWidth_[slot]);
//             Readers[i][slot]->AddVariable("eleSigmaIEtaIEtaFull5x5_Lead", &eleSigmaIEtaIEtaFull5x5_[slot]);
//             Readers[i][slot]->AddVariable("eleSigmaIPhiIPhiFull5x5_Lead", &eleSigmaIPhiIPhiFull5x5_[slot]);
//             Readers[i][slot]->AddVariable("eleR9Full5x5_Lead",            &eleR9Full5x5_[slot]);
//             Readers[i][slot]->AddVariable("eleBrem_Lead",                 &eleBrem_[slot]);
//             Readers[i][slot]->AddVariable("diTrk.Pt()",                   &diTrkPt_[slot]);
//             Readers[i][slot]->AddVariable("gsfDeltaR_Lead",               &gsfDeltaR_[slot]);
//             if (i == 1) // EE
//                 Readers[i][slot]->AddVariable("eleESEnToRawE_Lead",       &eleESEnToRawE_[slot]);

//             Readers[i][slot]->BookMVA("BDTG", path[i].c_str());
//         }
//     }

//     auto regression = [&, Readers](
//         unsigned int slot,
//         const int nEle,
//         const float rho,
//         const int nVtx,
//         const ROOT::RVec<float>& eleSCEta,
//         const ROOT::RVec<float>& eleSCPhi,
//         const ROOT::RVec<float>& eleSCRawEn,
//         const ROOT::RVec<float>& eleCalibPt,
//         const ROOT::RVec<float>& eledEtaAtVtx,
//         const ROOT::RVec<float>& eledPhiAtVtx,
//         const ROOT::RVec<float>& elePtError,
//         const ROOT::RVec<float>& eleHoverE,
//         const ROOT::RVec<float>& eleEoverP,
//         const ROOT::RVec<float>& eleEoverPout,
//         const ROOT::RVec<float>& eleEoverPInv,
//         const ROOT::RVec<float>& eleSCEtaWidth,
//         const ROOT::RVec<float>& eleSCPhiWidth,
//         const ROOT::RVec<float>& eleSigmaIEtaIEtaFull5x5,
//         const ROOT::RVec<float>& eleSigmaIPhiIPhiFull5x5,
//         const ROOT::RVec<float>& eleR9Full5x5,
//         const ROOT::RVec<float>& eleBrem,
//         const ROOT::RVec<float>& diTrkPt,
//         const ROOT::RVec<float>& gsfDeltaR,
//         const ROOT::RVec<float>& eleESEnToRawE
//     ){
//         ROOT::RVec<float> corrPt(nEle);
//         for (int i = 0; i < nEle; i++){
//             int iBE = (fabs(eleSCEta[i]) < 1.479) ? 0 : 1;

//             rho_[slot] = rho;
//             nVtx_[slot] = (float) nVtx;
//             eleSCEta_[slot] = eleSCEta[i];
//             eleSCPhi_[slot] = eleSCPhi[i];
//             eleSCRawEn_[slot] =  eleSCRawEn[i];
//             eleCalibPt_[slot] =  eleCalibPt[i];
//             eledEtaAtVtx_[slot] = eledEtaAtVtx[i];
//             eledPhiAtVtx_[slot] = eledPhiAtVtx[i];
//             elePtError_[slot] = elePtError[i];
//             eleHoverE_[slot] = eleHoverE[i];
//             eleEoverP_[slot] = eleEoverP[i];
//             eleEoverPout_[slot] = eleEoverPout[i];
//             eleEoverPInv_[slot] = eleEoverPInv[i];
//             eleSCEtaWidth_[slot] = eleSCEtaWidth[i];
//             eleSCPhiWidth_[slot] = eleSCPhiWidth[i];
//             eleSigmaIEtaIEtaFull5x5_[slot] = eleSigmaIEtaIEtaFull5x5[i];
//             eleSigmaIPhiIPhiFull5x5_[slot] = eleSigmaIPhiIPhiFull5x5[i];
//             eleR9Full5x5_[slot] = eleR9Full5x5[i];
//             eleBrem_[slot] = eleBrem[i];
//             diTrkPt_[slot] = diTrkPt[i];
//             gsfDeltaR_[slot] = gsfDeltaR[i];
//             eleESEnToRawE_[slot] = eleESEnToRawE[i];

//             corrPt[i] = eleCalibPt[i] * Readers[iBE][slot]->EvaluateRegression(0, "BDTG");
//         }
//         return corrPt;
//     };

//     auto nf = df.DefineSlot("eleHDALRegPt", regression, {
//         "nEle",
//         "rho",
//         "nVtx",
//         "eleSCEta",
//         "eleSCPhi",
//         "eleSCRawEn",
//         "eleCalibPt",
//         "eledEtaAtVtx",
//         "eledPhiAtVtx",
//         "elePtError",              
//         "eleHoverE",               
//         "eleEoverP",               
//         "eleEoverPout",            
//         "eleEoverPInv",            
//         "eleSCEtaWidth",           
//         "eleSCPhiWidth",           
//         "eleSigmaIEtaIEtaFull5x5", 
//         "eleSigmaIPhiIPhiFull5x5", 
//         "eleR9Full5x5",            
//         "eleBrem",                
//         "diTrkPt",
//         "gsfDeltaR",
//         "eleESEnToRawE"         
//     });

//     return nf;
// }