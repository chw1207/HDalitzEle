#include "ROOT/RVec.hxx"
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

        for (int i = 0; i < mcStatusFlag.size(); i++){
            if ((mcStatusFlag[i] >> num & 1) == 1) 
                vec.push_back(1);
            else
                vec.push_back(0);
        }

        return vec;
    }


    // Function to get the index vector sorted by pT
    // * Reference: https://root.cern/doc/master/vo006__IndexManipulation_8C.html
    ROOT::RVec<int> getIdx(ROOT::RVec<int>& isgood, ROOT::RVec<float>& pt){
        ROOT::RVec<int> idx_select = Nonzero(isgood);
        ROOT::RVec<int> idx_sort = Reverse(Argsort(pt));
        ROOT::RVec<int> idx = Intersect(idx_sort, idx_select);

        return idx;
    }


    // Function to match gen to reco electron 
    int RecoIdx(TLorentzVector gen, ROOT::RVec<float>& pt, ROOT::RVec<float>& eta, ROOT::RVec<float>& phi, float m){
        float min = 999.;
        int idx = -1;

        for (int i = 0; i < pt.size(); i++){
            TLorentzVector reco;
            reco.SetPtEtaPhiM(pt[i], eta[i], phi[i], m);

            if (gen.DeltaR(reco) > 0.1) continue;
            if (fabs((gen.Pt() / reco.Pt()) - 1.) < min){
                min = fabs((gen.Pt() / reco.Pt()) - 1.);
                idx = i;
            }
        }

        return idx;
    }


    // Reference: (only 2 track conv)
    // https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/MicroAOD/plugins/LegacyVertexSelector.cc#L499-L579
    // https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/CommonTools/Egamma/src/ConversionTools.cc#L103-L127
    // https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/CommonTools/Egamma/interface/ConversionTools.h#L48-L52
    int RecoConvMatch(float eleSCEta, float eleSCPhi, float eleSCEn, int nConv, ROOT::RVec<int>& convNTrks, ROOT::RVec<float>& convVtxX, ROOT::RVec<float>& convVtxY, ROOT::RVec<float>& convVtxZ, ROOT::RVec<float>& convFitPairPX, ROOT::RVec<float>& convFitPairPY, ROOT::RVec<float>& convFitPairPZ, ROOT::RVec<float>& convFitProb){
        TVector3 scpos;
        scpos.SetPtEtaPhi(eleSCEn/cosh(eleSCEta), eleSCEta, eleSCPhi);

        float mindR = 999;
        int selected_conversion_index = -1;
        for (int i = 0; i < nConv; i++){
            if (convNTrks[i] != 2) continue;

            float pairPt = sqrt(convFitPairPX[i]*convFitPairPX[i] + convFitPairPY[i]*convFitPairPY[i]);
            if (pairPt < 10.) continue;
            if (convFitProb[i] < 1e-6) continue;

            TVector3 VtxtoSC;
            VtxtoSC.SetXYZ(scpos.x() - convVtxX[i], scpos.y() - convVtxY[i], scpos.z() - convVtxZ[i]);

            TVector3 RefPairMo;
            RefPairMo.SetXYZ(convFitPairPX[i], convFitPairPY[i], convFitPairPZ[i]);

            float dR = VtxtoSC.DeltaR(RefPairMo);
            if(dR < mindR){
                mindR = dR;
                selected_conversion_index = i;
            }
        }

        if (mindR < 0.1)
            return selected_conversion_index;
        
        return -1;
    }


    int make_cat(bool Gen2Reco2, bool Gen2Reco1, float nGsfMatchToReco_lep1){
        int cat = 0;
        if (Gen2Reco2) 
            cat = 1;
        else if (Gen2Reco1 && nGsfMatchToReco_lep1 > 1)
            cat = 2;
        else if (Gen2Reco1 && nGsfMatchToReco_lep1 == 1)
            cat = 3;
        else{
            cout << "[ERROR] No proper category for Resolved" << endl;
            exit(-1);
        }
        return cat;
    }



}
