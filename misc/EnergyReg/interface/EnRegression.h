#include <vector>
#include "TString.h"
#include "ROOT/RVec.hxx"
#include "TFormula.h"
#include "XGBReader.h"
#include "MyRobustScaler.h"


ROOT::RDF::RNode doHDALRegressionXGB(ROOT::RDF::RNode df, std::vector<std::string> path, std::vector<std::string> spath){
    std::vector<std::shared_ptr<XGBReader>> model = {
        std::make_shared<XGBReader>(path[0]),
        std::make_shared<XGBReader>(path[1])
    };

    std::vector<std::shared_ptr<MyRobustScaler>> scaler = {
        std::make_shared<MyRobustScaler>(spath[0]),
        std::make_shared<MyRobustScaler>(spath[1])
    };

    auto regression = [&, model, scaler](
        const int nEle,
        const ROOT::RVec<float>& elePt,

        // features for EB
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
        // const ROOT::RVec<float>& eleEmax,
        // const ROOT::RVec<float>& eleE2nd,
        // const ROOT::RVec<float>& RatioRL,
        // const ROOT::RVec<float>& RatioBT,

        // additional feature for EE
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
                // eleEmax[i],
                // eleE2nd[i],
                // RatioRL[i],
                // RatioBT[i]
            };
            if (iBE == 1)
                features.push_back(eleESEnToRawE[i]);

            std::vector<float> features_scale = scaler[iBE]->Transform(features);
            // corrPt[i] = eleCalibPt[i] * model[iBE]->Compute(features_scale)[0];
            corrPt[i] = elePt[i] * model[iBE]->Compute(features_scale)[0];
        }
        return corrPt;
    };

    auto nf = df.Define("eleHDALRegPt", regression, {
        "nEle",
        "elePt",

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
        // "eleEmax",
        // "eleE2nd",
        // "RatioRL",
        // "RatioBT",

        "eleESEnToRawE"         
    });

    return nf;
}