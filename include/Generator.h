#ifndef GENERATOR_H_
#define GENERATOR_H_

#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"

namespace gen{

    // calculate di-electron mass in lhe level
    float CalcLHEMee(
        const int nLHE,
        const ROOT::RVec<int>& lhePID,
        const ROOT::RVec<float>& lhePx,
        const ROOT::RVec<float>& lhePy,
        const ROOT::RVec<float>& lhePz
    );

    // Find the events that photon is internal conversion
    bool IsPhoIntConv(
        const int nMC,
        const ROOT::RVec<int>& mcPID,
        const ROOT::RVec<int>& mcMomPID,
        const ROOT::RVec<unsigned short> mcStatusFlag
    );

    // Find the signal electrons
    ROOT::RVec<int> GenMatch(
        const ROOT::Math::PtEtaPhiMVector reco,
        const int nMC,
        const ROOT::RVec<float>& mcEta,
        const ROOT::RVec<float>& mcPhi,
        const ROOT::RVec<int>& mcPID,
        const ROOT::RVec<int>& mcMomPID,
        const ROOT::RVec<int>& mcGMomPID,
        const ROOT::RVec<unsigned short> mcStatusFlag
    );
}

#endif