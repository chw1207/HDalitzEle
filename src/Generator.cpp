#include <vector>
#include <limits>
#include "Generator.h"

float gen::CalcLHEMee(
    const int nLHE,
    const ROOT::RVec<int>& lhePID,
    const ROOT::RVec<float>& lhePx,
    const ROOT::RVec<float>& lhePy,
    const ROOT::RVec<float>& lhePz
){
    std::vector<float> eleIdx; // find lhe level electrons
    for (int i = 0; i < nLHE; i++){
        if (abs(lhePID[i]) != 11)
            continue;
        eleIdx.push_back(i);
    }

    float mass = std::numeric_limits<float>::max(); // assign a large value that will be filtered out
    if (eleIdx.size() > 1){
        ROOT::Math::PxPyPzMVector ele1(lhePx[eleIdx[0]], lhePy[eleIdx[0]], lhePz[eleIdx[0]], 0.000511);
        ROOT::Math::PxPyPzMVector ele2(lhePx[eleIdx[1]], lhePy[eleIdx[1]], lhePz[eleIdx[1]], 0.000511);
        mass = (ele1 + ele2).M();
    }
    return mass;
}


bool gen::IsPhoIntConv(
    const int nMC,
    const ROOT::RVec<int>& mcPID,
    const ROOT::RVec<int>& mcMomPID,
    const ROOT::RVec<unsigned short> mcStatusFlag
){
    int realpho = 0;
    for (int i = 0; i < nMC; i++){
        bool hardProc = (mcStatusFlag[i] >> 0 & 1) == 1;
        bool isPrompt = (mcStatusFlag[i] >> 1 & 1) == 1;
        bool isRealPho = mcPID[i] == 22 && mcMomPID[i] == 25 && hardProc && isPrompt;

        if (isRealPho)
            realpho += 1;
    }
    bool isPhoIntConv = (realpho == 0) ? true : false; // no real photon

    return isPhoIntConv;
}


ROOT::RVec<int> gen::GenMatch(
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