#include <vector>
#include <map>
#include <algorithm>
#include "MergedIDPred.h"
#include "XGBReader.h"
#include "ROOT/RVec.hxx"


ROOT::RDF::RNode MergedIDPred(ROOT::RDF::RNode df, std::string predColumn, std::map<std::string, std::string> modelMap){

    auto rM2EB = std::make_shared<XGBReader>(modelMap["M2EB"]);
    auto rM2EE = std::make_shared<XGBReader>(modelMap["M2EE"]);
    auto rM1EB = std::make_shared<XGBReader>(modelMap["M1EB"]);
    auto rM1EE = std::make_shared<XGBReader>(modelMap["M1EE"]);

    auto pred = [rM2EB, rM2EE, rM1EB, rM1EE](
        const int nEle,
        const ROOT::RVec<int>& nGsfMatchToReco,
        const float rho,
        const ROOT::RVec<float>& eleSCEta,
        const ROOT::RVec<float>& eleSCRawEn,
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
        const ROOT::RVec<float>& elePFChIso,
        const ROOT::RVec<float>& elePFPhoIso,
        const ROOT::RVec<float>& elePFNeuIso,
        const ROOT::RVec<float>& gsfPtRatio,
        const ROOT::RVec<float>& gsfDeltaR,
        const ROOT::RVec<float>& gsfRelPtRatio
    ){
        std::map<int, ROOT::RVec<float>> vmap;
        vmap.clear();

        ROOT::RVec<float> v;
        v.clear();

        ROOT::RVec<float> vs;
        vs.clear();

        for (int i = 0; i < nEle; i++){
            std::vector<float> features = {
                rho,
                eleSCEta[i],
                eleSCRawEn[i],
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
                elePFChIso[i],
                elePFPhoIso[i],
                elePFNeuIso[i],
                gsfPtRatio[i],
                gsfDeltaR[i],
                gsfRelPtRatio[i]
            };
            if (nGsfMatchToReco[i] == 1){
                features = {
                    rho,
                    eleSCEta[i],
                    eleSCRawEn[i],
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
                    elePFChIso[i],
                    elePFPhoIso[i],
                    elePFNeuIso[i],
                    // gsfPtRatio[i],
                    // gsfDeltaR[i],
                    gsfRelPtRatio[i]
                };
            }

            std::vector<float> scores;
            scores.clear();
            bool M2EB = nGsfMatchToReco[i] > 1 && fabs(eleSCEta[i]) < 1.479;
            bool M2EE = nGsfMatchToReco[i] > 1 && fabs(eleSCEta[i]) >= 1.479;
            bool M1EB = nGsfMatchToReco[i] == 1 && fabs(eleSCEta[i]) < 1.479;
            bool M1EE = nGsfMatchToReco[i] == 1 && fabs(eleSCEta[i]) >= 1.479;
            if (M2EB)
                scores = rM2EB->Compute(features);
            if (M2EE)
                scores = rM2EE->Compute(features);
            if (M1EB)
                scores = rM1EB->Compute(features);
            if (M1EE)
                scores = rM1EE->Compute(features);

            // choose the one with highest score
            v.push_back(std::max_element(scores.begin(), scores.end()) - scores.begin());
            vs.push_back(scores[0]);
        }

        vmap[0] = v;
        vmap[1] = vs;

        return vmap;
    };

    auto nf = df.Define("eleClassScore", pred, {
        "nEle",
        "nGsfMatchToReco",

        // training features
        "rho",
        "eleSCEta",
        "eleSCRawEn",

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

        "elePFChIso",
        "elePFPhoIso",
        "elePFNeuIso",

        "gsfPtRatio",
        "gsfDeltaR",
        "gsfRelPtRatio"
    }).Define(predColumn, "eleClassScore[0]").Define("eleXGBID", "eleClassScore[1]");

    return nf;
}