#ifndef CATEGORIZER_H_
#define CATEGORIZER_H_

#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"

namespace cat{
    ROOT::RVec<float> dRVector(
        const ROOT::Math::PtEtaPhiMVector lep,
        const ROOT::RVec<float>& eta,
        const ROOT::RVec<float>& phi
    );


    // AK4CHS jets: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
    ROOT::RVec<int> CutBasedJet(
        const int nJet,
        const ROOT::RVec<float>& jetEta,
        const ROOT::RVec<float>& jetNHF, // Neutral Hadron Energy Fraction
        const ROOT::RVec<float>& jetNEF, // Neutral EM Energy Fraction
        const ROOT::RVec<int>& jetNNP,   // Number of Neutral Particles
        const ROOT::RVec<int>& jetNCH,   // Charged Multplicity
        const ROOT::RVec<float>& jetCHF, // Charged Hadron Energy Fraction
        const ROOT::RVec<float>& jetCEF, // Charged EM Energy Fraction
        const ROOT::RVec<float>& jetMUF,  // Muon Fraction
        std::string era
    );


    // VBFtag: [tag_number, jet1_ind, jet2_ind]
    // tag_number: 1->isHVBF, 2->isLVBF, 3->isNone
    // https://github.com/root-project/opendata-benchmarks/blob/master/tasks/8/rdataframe_compiled.cxx#L9-L47
    ROOT::RVec<int> VBFtag(
        const ROOT::RVec<int>& isGoodJet,
        const ROOT::RVec<float>& jetPt,
        const ROOT::RVec<float>& jetEta,
        const ROOT::RVec<float>& jetPhi,
        const ROOT::RVec<float>& jetEn,
        const ROOT::Math::PtEtaPhiMVector H
    );

    int makeCat(
        bool isM2,
        bool isM1,
        bool isRe,
        bool isHVbf,
        bool isLVbf,
        bool isBst,
        bool isEE,
        bool isEBHR9,
        bool isEBLR9
    );
}

#endif