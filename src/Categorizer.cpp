#include "Categorizer.h"


ROOT::RVec<float> cat::dRVector(
    const ROOT::Math::PtEtaPhiMVector lep,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& phi
){
    ROOT::RVec<float> vec(eta.size());
    for (size_t i = 0; i < eta.size(); i++) {
        const float dr = ROOT::VecOps::DeltaR((float) lep.Eta(), eta[i], (float) lep.Phi(), phi[i]);
        vec[i] = dr;
    }
    return vec;
}


ROOT::RVec<int> cat::CutBasedJet(
    const int nJet,
    const ROOT::RVec<float>& jetEta,
    const ROOT::RVec<float>& jetNHF,
    const ROOT::RVec<float>& jetNEF,
    const ROOT::RVec<int>& jetNNP,
    const ROOT::RVec<int>& jetNCH,
    const ROOT::RVec<float>& jetCHF,
    const ROOT::RVec<float>& jetCEF,
    const ROOT::RVec<float>& jetMUF,
    std::string era
){
    ROOT::RVec<int> v(nJet);
    for (int i = 0; i < nJet; i++){
        bool cutID = false;

        if (era.find("2016") != std::string::npos){
            const bool eta1 = fabs(jetEta[i]) <= 2.4;
            const bool eta2 = fabs(jetEta[i]) > 2.4 && fabs(jetEta[i]) <= 2.7;
            const bool eta3 = fabs(jetEta[i]) > 2.7 && fabs(jetEta[i]) <= 3.;
            const bool eta4 = fabs(jetEta[i]) > 3. && fabs(jetEta[i]) <= 5.;

            int numConst = jetNCH[i] + jetNNP[i]; // Number of Constituents
            if (eta1)
                cutID = jetNHF[i] < 0.9 && jetNEF[i] < 0.9 && numConst > 1 /*&& jetMUF[i] < 0.8*/ && jetCHF[i] > 0. && jetNCH[i] > 0 /*&& jetCEF[i] < 0.8*/;
            else if (eta2)
                cutID = jetNHF[i] < 0.9 && jetNEF[i] < 0.9;
            else if (eta3)
                cutID = jetNHF[i] < 0.9 && jetNEF[i] < 0.99 && jetNEF[i] > 0. && jetNNP[i] > 1;
            else if (eta4)
                cutID = jetNHF[i] > 0.2 && jetNEF[i] < 0.9 && jetNNP[i] > 10;
        }
        else{
            const bool eta1 = fabs(jetEta[i]) <= 2.6;
            const bool eta2 = fabs(jetEta[i]) > 2.6 && fabs(jetEta[i]) <= 2.7;
            const bool eta3 = fabs(jetEta[i]) > 2.7 && fabs(jetEta[i]) <= 3.;
            const bool eta4 = fabs(jetEta[i]) > 3. && fabs(jetEta[i]) <= 5.;

            int numConst = jetNCH[i] + jetNNP[i]; // Number of Constituents
            if (eta1)
                cutID = jetNHF[i] < 0.9 && jetNEF[i] < 0.9 && numConst > 1 && jetCHF[i] > 0. && jetNCH[i] > 0 /*&& jetMUF[i] < 0.8 && jetCEF[i] < 0.8*/;
            else if (eta2)
                cutID = jetNHF[i] < 0.9 && jetNEF[i] < 0.99 && jetNCH[i] > 0 /*&& jetMUF[i] < 0.8 && jetCEF[i] < 0.8*/;
            else if (eta3)
                cutID = jetNEF[i] <  0.99 && jetNEF[i] > 0.01 && jetNNP[i] > 1;
            else if (eta4)
                cutID = jetNHF[i] > 0.2 && jetNEF[i] < 0.9 && jetNCH[i] > 10;
        }

        int pass = (cutID) ? 1 : 0;
        v[i] = pass;
    }

    return v;
}


ROOT::RVec<int> cat::VBFtag(
    const ROOT::RVec<int>& isGoodJet,
    const ROOT::RVec<float>& jetPt,
    const ROOT::RVec<float>& jetEta,
    const ROOT::RVec<float>& jetPhi,
    const ROOT::RVec<float>& jetEn,
    const ROOT::Math::PtEtaPhiMVector H
){
    int tag_number = 3;
    int jet1_ind = -1, jet2_ind = -1;
    ROOT::RVec<int> vbf{tag_number, jet1_ind, jet2_ind};

    if(jetPt.size() > 1){
        const auto c = ROOT::VecOps::Combinations(jetPt, 2);
        const auto make_p4 = [&](size_t idx) {
            ROOT::Math::PtEtaPhiEVector v(jetPt[idx], jetEta[idx], jetPhi[idx], jetEn[idx]);
            return v;
        };

        for (size_t i = 0; i < c[0].size(); i++){
            const auto i1 = c[0][i];
            const auto i2 = c[1][i];
            if(isGoodJet[i1] == 1 && isGoodJet[i2] == 1){
                const auto jet1 = make_p4(i1);
                const auto jet2 = make_p4(i2);
                const auto dijet = jet1 + jet2;
                const float zepen = H.Eta() - ((jetEta[i1] + jetEta[i2]) * 0.5);

                if (fabs(jet1.Eta() - jet2.Eta()) < 3.5)
                    continue;
                if (fabs(ROOT::VecOps::DeltaPhi(dijet.Phi(), H.Phi())) < 2.4)
                    continue;
                if (fabs(zepen) > 2.5)
                    continue;

                if (dijet.M() > 500.){
                    tag_number = 1;
                    jet1_ind = i1;
                    jet2_ind = i2;

                    vbf = {tag_number, jet1_ind, jet2_ind};
                    return vbf;
                }
                else if (dijet.M() <= 500. && dijet.M() > 360){
                    tag_number = 2;
                    jet1_ind = i1;
                    jet2_ind = i2;
                    vbf = {tag_number, jet1_ind, jet2_ind};
                    return vbf;
                }
            }
        }
    }

    return vbf;
}


int cat::makeCat(
    bool isM2,
    bool isM1,
    bool isRe,
    bool isHVbf,
    bool isLVbf,
    bool isBst,
    bool isEE,
    bool isEBHR9,
    bool isEBLR9
){
    int cat = 0;
    if (isM2){
        if (isHVbf == 1) cat = 1;
        else if (isLVbf == 1) cat = 2;
        else if (isBst == 1) cat = 3;
        else if (isEBHR9 == 1) cat = 4;
        else if (isEBLR9 == 1) cat = 5;
        else if (isEE == 1) cat = 6;
        else
            throw std::runtime_error("No proper category for M2");
    }
    else if (isM1 == 1){
        if (isHVbf == 1) cat = 7;
        else if (isLVbf == 1) cat = 8;
        else if (isBst == 1) cat = 9;
        else if (isEBHR9 == 1) cat = 10;
        else if (isEBLR9 == 1) cat = 11;
        else if (isEE == 1) cat = 12;
        else
            throw std::runtime_error("No proper category for M1");
    }
    else if (isRe == 1)
        cat = 13;
    else
        throw std::runtime_error("No proper category for Resolved");

    return cat;
}