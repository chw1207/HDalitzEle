#include <utility> // std::pair
#include "ROOT/RVec.hxx"
#include "TLorentzVector.h"
#include "TChain.h"

using namespace ROOT::VecOps;
using namespace std;

// Function to generate photon internal conversion weight
// no good idea to calculate the photon internal conversion weight via the RDataFrame :(
float phoIntConv(vector<string> fileVec){
    TChain* tree = new TChain("ggNtuplizer/EventTree");
    for(auto i : fileVec){
        tree->Add(i.c_str());
    }
    tree->SetBranchStatus("*",              0);
    tree->SetBranchStatus("nMC",            1);
    tree->SetBranchStatus("mcPID",          1);
    tree->SetBranchStatus("mcMomPID",       1);
    tree->SetBranchStatus("mcStatusFlag",   1);

    int nMC;      
    vector<int>* mcPID = new vector<int>(); 
    vector<int>* mcMomPID = new vector<int>();
    vector<unsigned short>* mcStatusFlag = new vector<unsigned short>();
    tree->SetBranchAddress("nMC",           &nMC);
    tree->SetBranchAddress("mcPID",         &mcPID);
    tree->SetBranchAddress("mcMomPID",      &mcMomPID);
    tree->SetBranchAddress("mcStatusFlag",  &mcStatusFlag);

    float IntPho = 0.;
    int tev = tree->GetEntries();
    for (int ev = 0; ev < tev; ev++){
        tree->GetEntry(ev);

        for (int j = 0; j < nMC; j++){
            if ((*mcPID)[j] != 22) continue;
            if ((*mcMomPID)[j] != 25) continue;
            if (((((*mcStatusFlag)[j] >> 0) & 1) == 1) && ((((*mcStatusFlag)[j] >> 1) & 1) == 1)) IntPho += 1;
        }
    }
    
    delete tree;
    const float wei = IntPho / (float)tev;
    return wei;
}

// Function to check if this is a main gsf track 
ROOT::RVec<int> IsMainGSF(ROOT::RVec<float>& eleD0, ROOT::RVec<float>& gsfD0){
    ROOT::RVec<int> v;
    v.clear();

    for (int i = 0; i < gsfD0.size(); i++){
        // if the GSF track with D0 the same as electron, then this GSF track is the main GSF track of the electron
        if (count(eleD0.begin(), eleD0.end(), gsfD0[i]))
            v.push_back(1);
        else 
            v.push_back(0);
    }

    const int nMainGSF = Sum(v);
    if (nMainGSF != eleD0.size()){
        cout << "[ERROR]: There is electron without proper matched GSF track!!!!!!" << endl;
        exit(-1);
    }

    return v;
} 


// Find the ambiguous gsf tracks of a electron
//* Description:
//*     1) The collection of GSF tracks are split into several chunks associated with the different electrons(ambiguousGSF)
//*     2) Start from the main GSF track of an electron, I keep adding the rest GSF tracks to this electron's ambiguousGSF until meeting the main GSF track of the other electron.
//*     3) No quality cuts applied 
//*     4) Each GSF track in the reco::GsfTrack should all pass the dEta and dPhi requirements described in the CMS reconstruction paper, otherwise, it will not be stored in the reco::GsfTrack.
//*     5) https://iopscience.iop.org/article/10.1088/1748-0221/16/05/P05014
ROOT::RVec<ROOT::RVec<int>> TrackElectronMatching(ROOT::RVec<float>& eleD0, ROOT::RVec<float>& gsfD0, ROOT::RVec<int>& isMainGSF){
    ROOT::RVec<ROOT::RVec<int>> v_associateIdx;
    v_associateIdx.clear();
    
    const ROOT::RVec<int> mainGsfIdx = Nonzero(isMainGSF);
    for (int i = 0; i < eleD0.size(); i++){
        
        const float mainGSFD0 = eleD0[i];
        pair<int, int> results(-1, -1);
        for (int j = 0; j < mainGsfIdx.size(); j++){
            if (gsfD0[mainGsfIdx[j]] == mainGSFD0){
                results.first = j;
                results.second = mainGsfIdx[j];
                break;
            }
        }

        // the first element is always the index of main gsf track  
        ROOT::RVec<int> associateIdx; // the indecies of gsf tracks associated with the electron. 
        associateIdx.clear();
        if (results.second == mainGsfIdx.back()){ // the last main gsf track in this event
            for (int k = results.second; k < gsfD0.size(); k++){
                associateIdx.push_back(k);
            }
        }
        else{
            for (int k = results.second; k < mainGsfIdx[results.first + 1]; k++){
                associateIdx.push_back(k);
            }
        }

        v_associateIdx.push_back(associateIdx);
    }
    
    return v_associateIdx;
}

// Function to calculate the number of ambiguous gsf tracks
ROOT::RVec<int> CalcNGsfMatchToReco(int nEle, ROOT::RVec<ROOT::RVec<int>>& ambiguousGSF){
    ROOT::RVec<int> v;
    v.clear();

    for (int i = 0; i < nEle; i++){
        v.push_back(ambiguousGSF[i].size());
    }

    return v;
}


// Function to Find the main gsf track
ROOT::RVec<int> FindMainGSF(int nEle, ROOT::RVec<ROOT::RVec<int>>& ambiguousGSF){
    ROOT::RVec<int> v;
    v.clear();

    for (int i = 0; i < nEle; i++){
        v.push_back(ambiguousGSF[i][0]); 
    }

    return v;
}


// Function to Find the second gsf track (pT max)
ROOT::RVec<int> FindSubGSF_PtMax(int nEle, ROOT::RVec<ROOT::RVec<int>>& ambiguousGSF, ROOT::RVec<float>& gsfPt){
    ROOT::RVec<int> v;
    v.clear();

    for (int i = 0; i < nEle; i++){

        const int nGsfMatchToReco = ambiguousGSF[i].size();
        const int eleTrkIdx = ambiguousGSF[i][0];

        float tmp = -999.;
        int eleSubTrkIdx = -1;
        for (int j = 1; j < nGsfMatchToReco; j++){
            if (gsfPt[ambiguousGSF[i][j]] > tmp){
                tmp = gsfPt[ambiguousGSF[i][j]];
                eleSubTrkIdx = ambiguousGSF[i][j];
            }
        }
        
        v.push_back(eleSubTrkIdx);
    }

    return v;
}


// Function to Find the second gsf track (dR min)
ROOT::RVec<int> FindSubGSF_dRMin(int nEle, ROOT::RVec<ROOT::RVec<int>>& ambiguousGSF, ROOT::RVec<float>& gsfEta, ROOT::RVec<float>& gsfPhi){
    ROOT::RVec<int> v;
    v.clear();

    for (int i = 0; i < nEle; i++){

        const int nGsfMatchToReco = ambiguousGSF[i].size();
        const int eleTrkIdx = ambiguousGSF[i][0];

        float tmp = 999.;
        int eleSubTrkIdx = -1;
        for (int j = 1; j < nGsfMatchToReco; j++){
            const float dR = DeltaR(gsfEta[eleTrkIdx], gsfEta[ambiguousGSF[i][j]], gsfPhi[eleTrkIdx], gsfPhi[ambiguousGSF[i][j]]);
            if (dR < tmp){
                tmp = dR;
                eleSubTrkIdx = ambiguousGSF[i][j];
            }
        }
        
        v.push_back(eleSubTrkIdx);
    }

    return v;
}


// Function to match the gsf variable to the electron based on the track index
//* Description:
//*     1) eleTrkPt = MatchIdex(eleTrkIdx, gsfPt)
//*     2) eleSubTrkPt = MatchIdex(eleSubTrkIdx, gsfPt)
template <typename T>
ROOT::RVec<T> MatchIdex(ROOT::RVec<int>& eleTrkIdx, ROOT::RVec<T>& var){
    ROOT::RVec<T> v;
    v.clear();

    for (int i = 0; i < eleTrkIdx.size(); i++){
        v.push_back(var[eleTrkIdx[i]]);
    }

    return v;
}


// Function to construct RVec of TLorentzVector
ROOT::RVec<TLorentzVector> P4Vector(ROOT::RVec<float>& pt, ROOT::RVec<float>& eta, ROOT::RVec<float>& phi, float m){
    ROOT::RVec<TLorentzVector> vec;
    const auto size = pt.size();
    for (auto i = 0; i < size; ++i) {
        TLorentzVector p4;
        p4.SetPtEtaPhiM(pt.at(i), eta.at(i), phi.at(i), m);
        vec.push_back(p4);
    }
    return vec;
}


// Function to get the pt ratio between electron trak and sub-track (if nGsfMatchToReco > 1)
ROOT::RVec<float> GetTrkPtRatio(int nEle, ROOT::RVec<int>& nGsfMatchToReco, ROOT::RVec<TLorentzVector>& eleTrk1, ROOT::RVec<TLorentzVector>& eleTrk2){
    ROOT::RVec<float> v;
    v.clear();
    for (int iele = 0; iele < nEle; iele++){
        if (nGsfMatchToReco[iele] > 1) v.push_back(eleTrk2[iele].Pt()/eleTrk1[iele].Pt());
        else v.push_back(-1.);
    }

    return v;
}


// Function to get the deltaR between electron trak and sub-track (if nGsfMatchToReco > 1)
ROOT::RVec<float> GetTrkdR(int nEle, ROOT::RVec<int>& nGsfMatchToReco, ROOT::RVec<TLorentzVector>& eleTrk1, ROOT::RVec<TLorentzVector>& eleTrk2){
    ROOT::RVec<float> v;
    v.clear();
    for (int iele = 0; iele < nEle; iele++){
        if (nGsfMatchToReco[iele] > 1) v.push_back(eleTrk1[iele].DeltaR(eleTrk2[iele]));
        else v.push_back(-999.);
    }

    return v;
}


// Function to get the deltaR between electron trak and sub-track (if nGsfMatchToReco > 1)
ROOT::RVec<float> GetTrkRelPtRatio(int nEle, ROOT::RVec<float>& eleCalibPt, ROOT::RVec<int>& nGsfMatchToReco, ROOT::RVec<TLorentzVector>& eleTrk1, ROOT::RVec<TLorentzVector>& eleTrk2){
    ROOT::RVec<float> v;
    v.clear();
    for (int iele = 0; iele < nEle; iele++){
        if (nGsfMatchToReco[iele] == 0) 
            v.push_back(-999.);
        else if (nGsfMatchToReco[iele] == 1) 
            v.push_back(eleTrk1[iele].Pt()/eleCalibPt[iele]);
        else 
            v.push_back((eleTrk1[iele] + eleTrk2[iele]).Pt()/eleCalibPt[iele]);
    }

    return v;
}



//// Function to match the gsf tracks to electron (old version)
//* Description:
//*     1) eleTrkIdx = GetTrkIdx(int nEle, eleCalibEn, eleEta, elePhi, nGSFTrk, gsfPt, gsfEta, gsfPhi)[0]
//*     2) eleSubTrkIdx = GetTrkIdx(int nEle, eleCalibEn, eleEta, elePhi, nGSFTrk, gsfPt, gsfEta, gsfPhi)[1]
//*     3) nGsfMatchToReco = GetTrkIdx(int nEle, eleCalibEn, eleEta, elePhi, nGSFTrk, gsfPt, gsfEta, gsfPhi)[2][0]
//*     4) Matching criteria: https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/RecoEgamma/EgammaElectronProducers/plugins/GsfElectronProducer.cc#L188-L191
// ROOT::RVec<ROOT::RVec<int>> GetTrkIdx(int nEle, ROOT::RVec<float>& eleCalibEn, ROOT::RVec<float>& eleEta, ROOT::RVec<float>& elePhi, int nGSFTrk, ROOT::RVec<float>& gsfPt, ROOT::RVec<float>& gsfEta, ROOT::RVec<float>& gsfPhi){
//     ROOT::RVec<int> LeadTrk; // leading track
//     LeadTrk.clear();
    
//     ROOT::RVec<int> SubLeadTrk; // trailing track
//     SubLeadTrk.clear();

//     ROOT::RVec<int> vngsf; // nTrks
//     vngsf.clear();

//     for (int iele = 0; iele < nEle; iele++){
//         TLorentzVector Ele;
//         Ele.SetPtEtaPhiE(eleCalibEn[iele]/cosh(eleEta[iele]), eleEta[iele], elePhi[iele], eleCalibEn[iele]);

//         ROOT::RVec<int> idx;
//         idx.clear();
        
//         ROOT::RVec<float> pt;
//         pt.clear();
//         for (int igsf = 0; igsf < nGSFTrk; igsf++){
//             TLorentzVector Gsf;
//             Gsf.SetPtEtaPhiM(gsfPt[igsf], gsfEta[igsf], gsfPhi[igsf], 0.000511);

//             // why some gsf trcks have very high pT?
//             if (Gsf.Pt()/Ele.Pt() > 5.) continue; // remove gsf track with unreasonable high pT
//             if (fabs(Ele.Eta() - Gsf.Eta()) > 0.02) continue;
//             if (Ele.DeltaPhi(Gsf) > 0.15) continue;
            
//             idx.push_back(igsf);
//             pt.push_back(gsfPt[igsf]);
//         }

//         const int ngsf = idx.size();
//         vngsf.push_back(ngsf);

//         // sort the idx via pT
//         ROOT::RVec<int> sort_idx = idx;
//         if (ngsf > 1){
//             auto sortIndices = Reverse(Argsort(pt)); // descending sorting
//             sort_idx = Take(idx, sortIndices);
//         }
        
//         if (ngsf != 0) LeadTrk.push_back(sort_idx[0]);
//         else LeadTrk.push_back(-1);
        
//         if (ngsf > 1) SubLeadTrk.push_back(sort_idx[1]);
//         else SubLeadTrk.push_back(-1);
//     }

//     ROOT::RVec<ROOT::RVec<int>> TrkIdx = {LeadTrk, SubLeadTrk, vngsf};

//     return TrkIdx;
// }


// // Function to find the main gsf track to electron (old version)
// //* Description:
// //*     1) pt, eta, phi of main gsf of electron are not stored in the ggNtuple
// //*     2) match eleD0 (main gsf track dxy) and gsfD0 to find the main gsf index 
// ROOT::RVec<int> FindMainGSF(ROOT::RVec<float>& eleD0, ROOT::RVec<float>& gsfD0){
//     ROOT::RVec<int> maingsfIdx;
//     maingsfIdx.clear();

//     for (int i = 0; i < eleD0.size(); i++){
//         int match = -1;
        
//         for (int j = 0; j < gsfD0.size(); j++){
//             if (fabs(eleD0[i] - gsfD0[j]) < 0.0001){
//                 match = j;
//                 break;
//             }
//         }
        
//         if (match == -1) {
//             cout << "[ERROR]: No matched main gsf track!!!!!!" << endl;
//             exit(-1);
//         }

//         maingsfIdx.push_back(match);
//     }

//     return maingsfIdx;
// }


// // Function to find the sub gsf track to electron (old version)
// //* Description:
// //*     1) find the gsf trck near the main gsf track
// //*     2) within cone size 0.1
// ROOT::RVec<ROOT::RVec<int>> FindSubGSF(ROOT::RVec<int>& eleTrkIdx, ROOT::RVec<float>& gsfPt, ROOT::RVec<float>& gsfEta, ROOT::RVec<float>& gsfPhi){
//     ROOT::RVec<int> subLeadTrk; // trailing track
//     subLeadTrk.clear();
    
//     ROOT::RVec<int> vngsf; // nGsfMatchToReco
//     vngsf.clear();

//     for (int i = 0; i < eleTrkIdx.size(); i++){

//         TLorentzVector mainGSF;
//         mainGSF.SetPtEtaPhiM(gsfPt[eleTrkIdx[i]], gsfEta[eleTrkIdx[i]], gsfPhi[eleTrkIdx[i]], 0.000511);

//         ROOT::RVec<int> idx;
//         idx.clear();
//         ROOT::RVec<float> pt;
//         pt.clear();
//         for (int j = 0; j < gsfPt.size(); j++){
//             if (j == eleTrkIdx[i]) continue; // cannot be main gsf track
//             if (gsfPt[j] > gsfPt[eleTrkIdx[i]]) continue; // pt < main gsf track pt

//             TLorentzVector subGSF;
//             subGSF.SetPtEtaPhiM(gsfPt[j], gsfEta[j], gsfPhi[j], 0.000511);

//             if (mainGSF.DeltaR(subGSF) > 0.1) continue;
            
//             idx.push_back(j);
//             pt.push_back(gsfPt[j]);
//         }

//         const int ngsf = idx.size();
//         vngsf.push_back(ngsf+1); // number of matched sub gsf tracks + 1 main gsf track = total number of gsf track matched to electron

//         // sort the idx via pT
//         ROOT::RVec<int> sort_idx = idx;
//         if (ngsf > 1){
//             auto sortIndices = Reverse(Argsort(pt)); // descending sorting
//             sort_idx = Take(idx, sortIndices);
//         }

//         if (ngsf > 0) subLeadTrk.push_back(sort_idx[0]); // most energenic one among all trailing gsf tracks
//         else subLeadTrk.push_back(-1);

//     }

//     ROOT::RVec<ROOT::RVec<int>> TrkIdx = {subLeadTrk, vngsf};

//     return TrkIdx;
// }



// // Function to generat the Gsf track informations (old version)
// // Reference:
// // [1] Matching criteria: https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/RecoEgamma/EgammaElectronProducers/plugins/GsfElectronProducer.cc#L188-L191
// auto Gsf_fun(int nEle, int nGSFTrk, ROOT::RVec<float>& eleCalibEn, ROOT::RVec<float>& eleEta, ROOT::RVec<float>& elePhi, ROOT::RVec<float>& gsfPt, ROOT::RVec<float>& gsfEta, ROOT::RVec<float>& gsfPhi){

//     ROOT::RVec<float> vnGsfmatchToReco;
//     vnGsfmatchToReco.clear();

//     ROOT::RVec<float> vGsfPtRatio;
//     vGsfPtRatio.clear();

//     ROOT::RVec<float> vGsfDeltaR;
//     vGsfDeltaR.clear();

//     ROOT::RVec<float> vGsfRelPtRatio;
//     vGsfRelPtRatio.clear();

//     for(int i = 0; i < nEle; i++){
//         vector<int> accInd;
//         accInd.clear();
        
//         TLorentzVector Ele;
//         Ele.SetPtEtaPhiE(eleCalibEn[i]/cosh(eleEta[i]), eleEta[i], elePhi[i], eleCalibEn[i]);

//         for (int j = 0; j < nGSFTrk; j++){
//             TLorentzVector Gsf;
//             Gsf.SetPtEtaPhiM(gsfPt[j], gsfEta[j], gsfPhi[j], 0.000511);

//             // if (Gsf.Pt()/Ele.Pt() > 5.) continue; // remove gsf track with unreasonable high pT
//             if (fabs(Ele.Eta() - Gsf.Eta()) > 0.02) continue;
//             if (Ele.DeltaPhi(Gsf) > 0.15) continue;

//             accInd.push_back(j);
//         }

//         float ngsf = accInd.size();
//         vnGsfmatchToReco.push_back(ngsf);

//         float ratio = -1., dR = -999.; // ngsf = 0 || 1 (should be the same as the value put into training)
//         float Relratio = -999; // ngsf = 0
//         if (ngsf > 1){

//             // argsort the index 
//             for (int j = 0; j < (accInd.size() - 1); j++) {
//                 for (int x = 0; x < (accInd.size() - j - 1); x++) {
//                     if (gsfPt[accInd[x]] < gsfPt[accInd[x + 1]]) {
//                         int temp = accInd[x];
//                         accInd[x] = accInd[x + 1];
//                         accInd[x + 1] = temp;
//                     }
//                 }
//             }

//             ratio = gsfPt[accInd[1]]/gsfPt[accInd[0]];
            
//             TLorentzVector Gsf1, Gsf2, digsf;
//             Gsf1.SetPtEtaPhiM(gsfPt[accInd[0]], gsfEta[accInd[0]], gsfPhi[accInd[0]], 0.000511);
//             Gsf2.SetPtEtaPhiM(gsfPt[accInd[1]], gsfEta[accInd[1]], gsfPhi[accInd[1]], 0.000511);

//             digsf = Gsf1 + Gsf2;
//             dR = Gsf1.DeltaR(Gsf2);
//             Relratio = digsf.Pt()/Ele.Pt();
//         }

//         if (ngsf == 1){
//             TLorentzVector digsf;
//             Relratio = gsfPt[accInd[0]]/Ele.Pt();
//         }

//         vGsfPtRatio.push_back(ratio);
//         vGsfDeltaR.push_back(dR);
//         vGsfRelPtRatio.push_back(Relratio);
//     }

//     vector<ROOT::RVec<float>> gsf_vec;
//     gsf_vec.clear();
//     gsf_vec.push_back(vnGsfmatchToReco);
//     gsf_vec.push_back(vGsfPtRatio);
//     gsf_vec.push_back(vGsfDeltaR);
//     gsf_vec.push_back(vGsfRelPtRatio);
    
//     return gsf_vec;
// }
