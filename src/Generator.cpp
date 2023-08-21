#include <vector>
#include <limits>
#include "Generator.h"

float gen::CalcLHEMee(
    const ROOT::RVec<int>& lhePID,
    const ROOT::RVec<float>& lhePx,
    const ROOT::RVec<float>& lhePy,
    const ROOT::RVec<float>& lhePz
){
    auto mask = abs(lhePID) == 11; // gen electron
    float mass = std::numeric_limits<float>::max();
    if (ROOT::VecOps::Nonzero(mask).size() == 2){
        auto LHEElePx = lhePx[mask];
        auto LHEElePy = lhePy[mask];
        auto LHEElePz = lhePz[mask];
        ROOT::Math::PxPyPzMVector ele1(LHEElePx[0], LHEElePy[0], LHEElePz[0], 0.000511);
        ROOT::Math::PxPyPzMVector ele2(LHEElePx[1], LHEElePy[1], LHEElePz[1], 0.000511);
        mass = (ele1 + ele2).M();
    }
    return mass;
}


bool gen::IsPhoIntConv(
    const ROOT::RVec<int>& mcPID,
    const ROOT::RVec<int>& mcMomPID,
    const ROOT::RVec<unsigned short> mcStatusFlag
){
    // remove internal conversion photon
    auto hardProc = ROOT::VecOps::Map(mcStatusFlag, [](unsigned short bit){return (bit >> 0) & 1;});
    auto isPrompt = ROOT::VecOps::Map(mcStatusFlag, [](unsigned short bit){return (bit >> 1) & 1;});
    auto isRealPh = mcPID == 22 && mcMomPID == 25 && hardProc && isPrompt;
    
    int realpho = ROOT::VecOps::Sum(isRealPh);
    bool isPhoIntConv = (realpho == 0) ? true : false; // no real photon

    return isPhoIntConv;
}


ROOT::RVec<int> gen::MatchedGenEle(
    const ROOT::Math::PtEtaPhiMVector reco,
    const int nMC,
    const ROOT::RVec<float>& mcEta,
    const ROOT::RVec<float>& mcPhi,
    const ROOT::RVec<int>& mcPID,
    const ROOT::RVec<int>& mcMomPID,
    const ROOT::RVec<int>& mcGMomPID,
    const ROOT::RVec<unsigned short> mcStatusFlag
){
    ROOT::RVec<int> idx;
    idx.clear();

    for (int i = 0; i < nMC; i++){
        const float dr = ROOT::VecOps::DeltaR(mcEta[i], (float)reco.Eta(), mcPhi[i], (float)reco.Phi());
        if (dr > 0.1)
            continue;

        bool isDalEle = abs(mcPID[i]) == 11 && mcMomPID[i] == 25;
        bool isHZgEle = abs(mcPID[i]) == 11 && mcMomPID[i] == 23 && mcGMomPID[i] == 25;
        if (!isDalEle && !isHZgEle)
            continue;

        if (((mcStatusFlag[i] >> 0) & 1) != 1) // from hard process
            continue;

        if (((mcStatusFlag[i] >> 1) & 1) != 1) // prompt final state
            continue;

        idx.push_back(i);
    }

    return idx;
}


int gen::MatchedGenPho(
    const ROOT::Math::PtEtaPhiMVector reco,
    const int nMC,
    const ROOT::RVec<float>& mcEta,
    const ROOT::RVec<float>& mcPhi,
    const ROOT::RVec<int>& mcPID,
    const ROOT::RVec<int>& mcMomPID,
    const ROOT::RVec<unsigned short> mcStatusFlag
){
    for (int i = 0; i < nMC; i++){
        const float dr = ROOT::VecOps::DeltaR(mcEta[i], (float)reco.Eta(), mcPhi[i], (float)reco.Phi());
        bool isDalPho = mcPID[i] == 22 && mcMomPID[i] == 25;
        bool isHardPro = ((mcStatusFlag[i] >> 0) & 1) == 1;
        bool isPrompt = ((mcStatusFlag[i] >> 1) & 1) == 1;
        
        if (dr < 0.1 && isDalPho && isHardPro && isPrompt)
            return i;
    }
    return -1;
}