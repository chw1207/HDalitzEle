#include "ROOT/RVec.hxx"


ROOT::RVec<float> dRVector(TLorentzVector lep, ROOT::RVec<float>& pt, ROOT::RVec<float>& eta, ROOT::RVec<float>& phi, ROOT::RVec<float>& e){
    ROOT::RVec<float> vec;
    vec.clear();
    
    for (auto i = 0; i < pt.size(); i++) {
        TLorentzVector v;
        v.SetPtEtaPhiE(pt[i], eta[i], phi[i], e[i]);
        
        const float dr = lep.DeltaR(v);
        vec.push_back(dr);
    }
    return vec;
}


// Reference:
// [1] https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
bool JetCutBasedID(int ijet, ROOT::RVec<float> jetEta, ROOT::RVec<float> jetNHF, ROOT::RVec<float> jetNEF, ROOT::RVec<int> jetNNP, ROOT::RVec<int> jetNCH, ROOT::RVec<float> jetCHF, ROOT::RVec<float> jetCEF, int year){
    
    bool isGoodJet = false;

    int numConst = jetNCH[ijet] + jetNNP[ijet];
    if (year == 2016){
        if (fabs(jetEta[ijet]) <= 2.4){
            if ((jetNHF[ijet] < 0.9) && (jetNEF[ijet] < 0.9) && (numConst > 1) && (jetCHF[ijet] > 0.) && (jetNCH[ijet] > 0)){
               isGoodJet = true; 
            }
        }
        else if ((fabs(jetEta[ijet]) > 2.4) && (fabs(jetEta[ijet]) <= 2.7)){
            if ((jetNHF[ijet] < 0.9) && (jetNEF[ijet] < 0.99)){
               isGoodJet = true; 
            }
        }
        else if ((fabs(jetEta[ijet]) > 2.7) && (fabs(jetEta[ijet]) <= 3.)){
            if ((jetNHF[ijet] < 0.9) && (jetNEF[ijet] < 0.99) && (jetNEF[ijet] > 0.) && (jetNNP[ijet] > 1)){
               isGoodJet = true; 
            }
        }
        else if ((fabs(jetEta[ijet]) > 3) && (fabs(jetEta[ijet]) <= 5.)){
            if ((jetNHF[ijet] > 0.2) && (jetNEF[ijet] < 0.9) && (jetNNP[ijet] > 10)){
               isGoodJet = true; 
            }
        }
    }

    else if (year == 2017 || year == 2018){
        if (fabs(jetEta[ijet]) <= 2.6){
            if ((jetNHF[ijet] < 0.9) && (jetNEF[ijet] < 0.9) && (numConst > 1) && (jetCHF[ijet] > 0.) && (jetNCH[ijet] > 0)){
               isGoodJet = true; 
            }
        }
        else if ((fabs(jetEta[ijet]) > 2.6) && (fabs(jetEta[ijet]) <= 2.7)){
            if ((jetNHF[ijet] < 0.9) && (jetNEF[ijet] < 0.99) && (jetNCH[ijet] > 0)){
               isGoodJet = true; 
            }
        }
        else if ((fabs(jetEta[ijet]) > 2.7) && (fabs(jetEta[ijet]) <= 3.)){
            if (jetNHF[ijet] <  0.9999){
               isGoodJet = true; 
            }
        }
        else if ((fabs(jetEta[ijet]) > 3) && (fabs(jetEta[ijet]) <= 5.)){
            if ((jetNEF[ijet] < 0.9) && (jetNNP[ijet] > 2)){
               isGoodJet = true; 
            }
        }
    }

    else {
        printf("ERROR: No specific year!");
        exit(-1);
    }

    return isGoodJet;
}


ROOT::RVec<int> CutBasedJet(int nJet, ROOT::RVec<float>& jetEta, ROOT::RVec<float>& jetNHF, ROOT::RVec<float>& jetNEF, ROOT::RVec<int>& jetNNP, ROOT::RVec<int>& jetNCH, ROOT::RVec<float>& jetCHF, ROOT::RVec<float>& jetCEF, int year){
    ROOT::RVec<int> v;
    v.clear();
    
    for (int i = 0; i < nJet; i++){
        const bool isGoodJet = JetCutBasedID(
            i,
            jetEta, jetNHF, jetNEF, jetNNP, jetNCH, jetCHF, jetCEF,
            year
        ); 
        
        if (isGoodJet) v.push_back(1);
        else v.push_back(0);
    }
    
    return v;
}

// VBFtag: [tag_number, jet1_ind, jet2_ind]
// tag_number: 1->isHVBF, 2->isLVBF, 3->isNone
// https://github.com/root-project/opendata-benchmarks/blob/master/tasks/8/rdataframe_compiled.cxx#L9-L47
ROOT::RVec<int> VBFtag(ROOT::RVec<int>&isGoodJet, ROOT::RVec<float>&jetPt, ROOT::RVec<float>&jetEta, ROOT::RVec<float>&jetPhi, ROOT::RVec<float>&jetEn, TLorentzVector H){
    int tag_number = 3;
    int jet1_ind = -1, jet2_ind = -1;
    ROOT::RVec<int> vbf{tag_number, jet1_ind, jet2_ind};
    
    if(jetPt.size() > 1) {
        const auto c = ROOT::VecOps::Combinations(jetPt, 2);
        const auto make_p4 = [&](std::size_t idx) {
            TLorentzVector v;
            v.SetPtEtaPhiE(jetPt[idx], jetEta[idx], jetPhi[idx], jetEn[idx]);
            return v;
        };
    
        for (auto i = 0; i < c[0].size(); i++){
            const auto i1 = c[0][i];
            const auto i2 = c[1][i];
            if(isGoodJet[i1] == 1 && isGoodJet[i2] == 1){
                const auto jet1 = make_p4(i1);
                const auto jet2 = make_p4(i2);
                const auto dijet = jet1 + jet2;
                const float zepen = H.Eta() - ((jetEta[i1] + jetEta[i2]) * 0.5);
                
                if (fabs(jet1.Eta() - jet2.Eta()) < 3.5) continue;
                if (fabs(dijet.DeltaPhi(H)) < 2.4) continue;
                if (fabs(zepen) > 2.5) continue;
                
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