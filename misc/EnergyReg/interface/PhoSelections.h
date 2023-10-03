#include <cmath>
#include <limits>
#include <algorithm>
#include "ROOT/RVec.hxx"
#include "TVector3.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"


// FSR selection and Hgg preselection
// https://arxiv.org/pdf/2012.06888.pdf
// https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/Taggers/python/flashggPreselectedDiPhotons_cfi.py#L72
// https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/MicroAOD/plugins/MuMuGammaProducer.cc


ROOT::RVec<int> FSRSelec(
    const ROOT::Math::PtEtaPhiMVector& mu1,
    const ROOT::Math::PtEtaPhiMVector& mu2,
    const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& phoP4
){
    ROOT::RVec<int> v(phoP4.size());
    for (size_t i = 0; i < phoP4.size(); i++){
        int isfsr = 0;
        
        //reject the event that neither of the muons emitted the photon (dR < 0.8)
        const bool cut1 = ROOT::VecOps::DeltaR(mu1.Eta(), phoP4[i].Eta(), mu1.Phi(), phoP4[i].Phi()) < 0.8
                            || ROOT::VecOps::DeltaR(mu2.Eta(), phoP4[i].Eta(), mu2.Phi(), phoP4[i].Phi()) < 0.8;

        // remove this cut to increase statistics
        // prevent photon picked up the track from one of the muons 
        const bool cut2 = ROOT::VecOps::DeltaR(mu1.Eta(), phoP4[i].Eta(), mu1.Phi(), phoP4[i].Phi()) > 0.1
                            && ROOT::VecOps::DeltaR(mu2.Eta(), phoP4[i].Eta(), mu2.Phi(), phoP4[i].Phi()) > 0.1;

        //reject isr phoP4[i]tons (muug of isr is higher)
        const bool cut3 = ((mu1 + mu2).M() + (mu1 + mu2 + phoP4[i]).M()) < 180;

        //three body mass must be close to Z mass
        const bool cut4 = (mu1 + mu2 + phoP4[i]).M() > 80 && (mu1 + mu2 + phoP4[i]).M() < 100;

        if (cut1 && cut2 && cut3 && cut4)
            isfsr = 1;

        v[i] = isfsr;
    }

    return v;
}


int GetZPho(
    const ROOT::Math::PtEtaPhiMVector& mu1,
    const ROOT::Math::PtEtaPhiMVector& mu2,
    const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& phoP4,
    const ROOT::RVec<int> &isGoodPho
){
    int finipho = -1;
    float mindiff = std::numeric_limits<float>::max();
    for (size_t i = 0; i < phoP4.size(); i++){
        if (isGoodPho[i] == 1){
            float Muug = (mu1 + mu2 + phoP4[i]).M();
            float diff = abs(91.188 - Muug); //close to Z mass

            if (mindiff > diff){
                mindiff = diff;
                finipho = i;
            }
        }
    }
    return finipho;
}


// find which bin a value fall into
int FindBins(const std::vector<float> bin, const float var){
    auto lower = std::lower_bound(bin.begin(), bin.end(), var);
    int position = std::distance(bin.begin(), lower) - 1;
    if (position < 0 || position > (int) bin.size()-2) // underflow or overflow
        return (int) -1;
    return position;
}


ROOT::RVec<int> HggPreSelection(
    const ROOT::Math::PtEtaPhiMVector& mu1_old,
    const ROOT::Math::PtEtaPhiMVector& mu2_old,
    const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& phoP4,
    const float rhoAll,
    const int nPho,
    const ROOT::RVec<float>& phoSCEta,
    const ROOT::RVec<float>& phoPFChIso,
    const ROOT::RVec<float>& phoPFPhoIso,
    const ROOT::RVec<float>& phoTrkIsoHollowConeDR03,
    const ROOT::RVec<float>& phoR9Full5x5,
    const ROOT::RVec<float>& phoEt,
    const ROOT::RVec<float>& phoSigmaIEtaIEtaFull5x5,
    const ROOT::RVec<float>& phoHoverE
){
    ROOT::RVec<int> vec;
    vec.clear();

    const std::vector<float> bins = {0., 0.9, 1.5, 2, 2.2, 3};
    const std::vector<float> reff = {0.16544, 0.16544, 0.13212, 0.13212, 0.13212};

    for (int i = 0; i < nPho; i++){
        // to mitigate PU effect
        const int effbin = FindBins(bins, fabs(phoSCEta[i]));
        const float phoEffArea = (effbin != -1) ? reff[effbin] : 0.;
        const float phoPFPhoIso_corr = std::max(phoPFPhoIso[i] - rhoAll * phoEffArea, (float) 0.);

        const bool isEB = (fabs(phoSCEta[i]) < 1.4442);
        const bool isEE = (fabs(phoSCEta[i]) > 1.566 && fabs(phoSCEta[i]) < 2.5);

        // the muon track pT is not summed if the muon track from the dimuon candidate is within the isolation photon
        // to increase the statistic
        float phoTrkIsoHollowConeDR03_mucorr = phoTrkIsoHollowConeDR03[i];
        const float deltR1 = ROOT::VecOps::DeltaR(mu1_old.Eta(), phoP4[i].Eta(), mu1_old.Phi(), phoP4[i].Phi());
        const float deltR2 = ROOT::VecOps::DeltaR(mu2_old.Eta(), phoP4[i].Eta(), mu2_old.Phi(), phoP4[i].Phi());
        if (deltR1 < 0.3 && phoTrkIsoHollowConeDR03[i] > 0.99 * mu1_old.Pt()){
            phoTrkIsoHollowConeDR03_mucorr -= mu1_old.Pt();
        }
        if (deltR2 < 0.3 && phoTrkIsoHollowConeDR03[i] > 0.99 * mu2_old.Pt()){
            phoTrkIsoHollowConeDR03_mucorr -= mu2_old.Pt();
        }

        const bool isHR9_EB = isEB && phoR9Full5x5[i] > 0.85;
        const bool isLR9_EB = isEB && phoR9Full5x5[i] <= 0.85 && phoR9Full5x5[i] > 0.5 && phoPFPhoIso_corr < 4 && phoTrkIsoHollowConeDR03_mucorr < 6 && phoSigmaIEtaIEtaFull5x5[i] < 0.015;
        const bool isHR9_EE = isEE && phoR9Full5x5[i] > 0.9;
        const bool isLR9_EE = isEE && phoR9Full5x5[i] <= 0.9 && phoR9Full5x5[i] > 0.8 && phoPFPhoIso_corr < 4 && phoTrkIsoHollowConeDR03_mucorr < 6 && phoSigmaIEtaIEtaFull5x5[i] < 0.035;

        // cuts here mimic the miniAOD photon cuts and the non-category based trigger cuts
        const bool isAOD = phoHoverE[i] < 0.08 && (phoR9Full5x5[i] > 0.8 || phoPFChIso[i] < 20. || (phoPFChIso[i]/phoEt[i] < 0.3));
        const bool base_pT_cut = phoP4[i].Pt() > 25.;

        // Hgg preselection without passElectronVeto
        const int isHgg = (base_pT_cut && isAOD && (isHR9_EB || isLR9_EB || isHR9_EE || isLR9_EE)) ? 1 : 0;
        vec.push_back(isHgg);
    }

    return vec;
}


int GetMatchEle(const float phoSCEta, const float phoSCPhi, const ROOT::RVec<float>& eleSCEta, const ROOT::RVec<float>& eleSCPhi){
    int idx = -1;
    for (size_t i = 0; i < eleSCEta.size(); i++){
        if ((phoSCEta == eleSCEta[i]) && (phoSCPhi == eleSCPhi[i])){
            idx = i;
            break;
        }
    }
    return idx;
}


// matching of conversion vertex
int GetConv(const ROOT::Math::RhoEtaPhiVectorF& sc,
            const ROOT::RVec<ROOT::Math::XYZVectorF>& vtxs,
            const ROOT::RVec<ROOT::Math::XYZVectorF>& vtxPs,
            const ROOT::RVec<int>& isGoodVtx){
    
    int idx = -1;
    float mindR = std::numeric_limits<float>::max();
    for (size_t i = 0; i < vtxs.size(); i++){
        if (isGoodVtx[i] < 1)
            continue;
        
        ROOT::Math::XYZVectorF vtx2sc(sc.x() - vtxs[i].x(), sc.y() - vtxs[i].y(), sc.z() - vtxs[i].z());
        float dR = ROOT::VecOps::DeltaR(vtx2sc.Eta(), vtxPs[i].Eta(), vtx2sc.Phi(), vtxPs[i].Phi());
        if (dR < mindR){
            mindR = dR;
            idx = i;
        }
    }

    return (mindR < 0.1) ? idx : -1;
}



int GenMatchParticles(
    const ROOT::Math::PtEtaPhiMVector reco,
    const int nMC,
    const ROOT::RVec<float>& mcPt,
    const ROOT::RVec<float>& mcEta,
    const ROOT::RVec<float>& mcPhi,
    const ROOT::RVec<float>& mcMass
){
    float min = std::numeric_limits<float>::max();
    int idx = -1;
    for (int i = 0; i < nMC; i++){

        ROOT::Math::PtEtaPhiMVector gen(mcPt[i], mcEta[i], mcPhi[i], mcMass[i]);

        if (ROOT::VecOps::DeltaR((float)reco.Eta(), mcEta[i], (float)reco.Phi(), mcPhi[i]) > 0.1)
            continue;

        if (fabs((gen.Pt() / reco.Pt()) - 1.) < min){
            min = fabs((gen.Pt() / reco.Pt()) - 1.);
            idx = i;
        }
    }

    return idx;
}


