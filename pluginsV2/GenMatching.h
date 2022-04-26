#include "ROOT/RVec.hxx"

using namespace ROOT::VecOps;

ROOT::RVec<int> GenMatch(TLorentzVector reco, int nMC, ROOT::RVec<float>& mcPt, ROOT::RVec<float>& mcEta, ROOT::RVec<float>& mcPhi, ROOT::RVec<float>& mcMass, ROOT::RVec<int>& mcPID, ROOT::RVec<int>& mcMomPID, ROOT::RVec<int>& mcGMomPID, ROOT::RVec<unsigned short> mcStatusFlag){
    ROOT::RVec<int> idx;
    idx.clear();

    TLorentzVector gen;
    for (int i = 0; i < nMC; i++){
        gen.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], mcMass[i]);

        if (gen.DeltaR(reco) > 0.1) continue;

        bool isDalEle = abs(mcPID[i]) == 11 && mcMomPID[i] == 25;
        bool isHZgEle = abs(mcPID[i]) == 11 && mcMomPID[i] == 23 && mcGMomPID[i] == 25;
        if (!isDalEle && !isHZgEle) continue;
        if (((mcStatusFlag[i] >> 0) & 1) != 1) continue; // from hard process
        if (((mcStatusFlag[i] >> 1) & 1) != 1) continue; // prompt final state

        idx.push_back(i);
    }

    return idx;
}