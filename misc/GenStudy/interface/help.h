#include <algorithm>
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;
using namespace ROOT::VecOps;

namespace Helper{
    // Note: mcStatusFlag
    //*     1) num = 0 -> from Hard Process Final State
    //*     2) num = 1 -> is Prompt Final State
    ROOT::RVec<int> GenType(ROOT::RVec<unsigned short>& mcStatusFlag, int num){
        ROOT::RVec<int> vec;
        vec.clear();

        for (size_t i = 0; i < mcStatusFlag.size(); i++){
            if ((mcStatusFlag[i] >> num & 1) == 1)
                vec.push_back(1);
            else
                vec.push_back(0);
        }

        return vec;
    }


    // Function to match reco to gen electron
    int RecoIdx(
        const ROOT::Math::PtEtaPhiMVector& gen,
        const ROOT::RVec<float>& pt,
        const ROOT::RVec<float>& eta,
        const ROOT::RVec<float>& phi,
        const float m
    ){
        float min = std::numeric_limits<float>::max();
        int idx = -1;
        for (size_t i = 0; i < pt.size(); i++){
            ROOT::Math::PtEtaPhiMVector reco(pt[i], eta[i], phi[i], m);
            if (DeltaR(gen.Eta(), reco.Eta(), gen.Phi(), reco.Phi()) > 0.1)
                continue;

            if (fabs((gen.Pt() / reco.Pt()) - 1.) < min){
                min = fabs((gen.Pt() / reco.Pt()) - 1.);
                idx = i;
            }
        }

        return idx;
    }


    // only used for xAnaGenQCD.C
    ROOT::RVec<int> getGenIdxVec(
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> reco,
        const ROOT::RVec<float>& mcPt,
        const ROOT::RVec<float>& mcEta,
        const ROOT::RVec<float>& mcPhi,
        const ROOT::RVec<float>& mcMass,
        const ROOT::RVec<int>& mcPID,
        const ROOT::RVec<int>& mcMomPID
    ){
        ROOT::RVec<int> idxVec;
        idxVec.clear();

        for (int i = 0; i < reco.size(); i++){
            float min = std::numeric_limits<float>::max();
            int idx = -1;
            for (int j = 0; j < mcPt.size(); j++){
                ROOT::Math::PtEtaPhiMVector gen(mcPt[j], mcEta[j], mcPhi[j], mcMass[j]);
                if (DeltaR(gen.Eta(), reco[i].Eta(), gen.Phi(), reco[i].Phi()) > 0.1)
                    continue;

                bool isgoodgen = (abs(mcPID[j]) == 11) && (mcMomPID[j] != 22) && (fabs(mcEta[j]) < 2.5);
                if (!isgoodgen)
                    continue;

                if (fabs((gen.Pt() / reco[i].Pt()) - 1.) < min){
                    min = fabs((gen.Pt() / reco[i].Pt()) - 1.);
                    idx = j;
                }
            }
            idxVec.push_back(idx);
        }

        return idxVec;
    }


    // only used for xAnaGenGJets.C
    ROOT::RVec<int> getGenPhoIdxVec(
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> reco,
        const ROOT::RVec<float>& mcPt,
        const ROOT::RVec<float>& mcEta,
        const ROOT::RVec<float>& mcPhi,
        const ROOT::RVec<float>& mcMass,
        const ROOT::RVec<int>& mcPID
    ){
        ROOT::RVec<int> idxVec;
        idxVec.clear();

        for (int i = 0; i < reco.size(); i++){
            float min = std::numeric_limits<float>::max();
            int idx = -1;
            for (int j = 0; j < mcPt.size(); j++){
                ROOT::Math::PtEtaPhiMVector gen(mcPt[j], mcEta[j], mcPhi[j], mcMass[j]);
                if (DeltaR(gen.Eta(), reco[i].Eta(), gen.Phi(), reco[i].Phi()) > 0.1)
                    continue;

                bool isgoodgen = (mcPID[j] == 22) && (fabs(mcEta[j]) < 2.5);
                if (!isgoodgen)
                    continue;

                if (fabs((gen.Pt() / reco[i].Pt()) - 1.) < min){
                    min = fabs((gen.Pt() / reco[i].Pt()) - 1.);
                    idx = j;
                }
            }
            idxVec.push_back(idx);
        }

        return idxVec;
    }

    // Reference:
    // https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/MicroAOD/plugins/LegacyVertexSelector.cc#L499-L579
    // https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/CommonTools/Egamma/src/ConversionTools.cc#L103-L127
    // https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/CommonTools/Egamma/interface/ConversionTools.h#L48-L52
    ROOT::RVec<int> FindConvVtxWrtEle(
        const int nEle,
        const ROOT::RVec<float>& eleSCEta,
        const ROOT::RVec<float>& eleSCPhi,
        const ROOT::RVec<float>& eleSCEn,
        const int nConv,
        const ROOT::RVec<int>& convNTrks,
        const ROOT::RVec<float>& convVtxX,
        const ROOT::RVec<float>& convVtxY,
        const ROOT::RVec<float>& convVtxZ,
        const ROOT::RVec<float>& convFitPairPX,
        const ROOT::RVec<float>& convFitPairPY,
        const ROOT::RVec<float>& convFitPairPZ,
        const ROOT::RVec<float>& convFitProb
    ){
        ROOT::RVec<int> matchedVtx;
        matchedVtx.clear();

        for (int iele = 0; iele < nEle; iele++){
            TVector3 SC;
            SC.SetPtEtaPhi(eleSCEn[iele]/cosh(fabs(eleSCEta[iele])), eleSCEta[iele], eleSCPhi[iele]);
        
            float mindR = std::numeric_limits<float>::max();
            int selected_conversion_index = -1;
        
            for (int iconv = 0; iconv < nConv; iconv++){
                if (convNTrks[iconv] != 2) //vertex validity
                    continue;

                const float pairPt = sqrt(convFitPairPX[iconv]*convFitPairPX[iconv] + convFitPairPY[iconv]*convFitPairPY[iconv]);
                if (pairPt < 10.)
                    continue;
                if (convFitProb[iconv] < 1e-6) //fit probability
                    continue;

                TVector3 VtxtoSC;
                VtxtoSC.SetXYZ(SC.x() - convVtxX[iconv], SC.y() - convVtxY[iconv], SC.z() - convVtxZ[iconv]);

                TVector3 RefPairMo;
                RefPairMo.SetXYZ(convFitPairPX[iconv], convFitPairPY[iconv], convFitPairPZ[iconv]);

                float dR = VtxtoSC.DeltaR(RefPairMo);
                if(dR < mindR){
                    mindR = dR;
                    selected_conversion_index = iconv;
                }
            }

            if (mindR < 0.1)
                matchedVtx.push_back(selected_conversion_index);
            else
                matchedVtx.push_back(-1);
        }
        
        return  matchedVtx;
    }


    int make_cat(bool Gen2Reco2, bool Gen2Reco1, float nGsfMatchToReco_Lead){
        int cat = 0;
        if (Gen2Reco2)
            cat = 1;
        else if (Gen2Reco1 && nGsfMatchToReco_Lead > 1)
            cat = 2;
        else if (Gen2Reco1 && nGsfMatchToReco_Lead == 1)
            cat = 3;
        else{
            throw std::runtime_error("[ERROR] No proper category");
        }
        return cat;
    }

}
