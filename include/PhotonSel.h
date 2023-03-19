#ifndef PHOTONSEL_H_
#define PHOTONSEL_H_

#include "ROOT/RVec.hxx"

namespace phoSel{
    // It requires the photon to pass H->gg preselection
    // H->gg preselection: https://arxiv.org/pdf/1804.02716.pdf
    // https://github.com/cms-analysis/flashgg/blob/1453740b1e4adc7184d5d8aa8a981bdb6b2e5f8e/Taggers/python/flashggPreselectedDiPhotons_cfi.py
    // https://github.com/andreypz/nwu-dalitz-analysis/blob/master/zgamma/egamma.C#L493
    ROOT::RVec<int> HggPresel(
        const int nPho,
        const float rhoAll,
        const ROOT::RVec<float>& phoSCEta,
        const ROOT::RVec<float>& phoPFChIso,
        const ROOT::RVec<float>& phoPFPhoIso,
        const ROOT::RVec<float>& phoTrkIsoHollowConeDR03,
        const ROOT::RVec<float>& phoR9Full5x5,
        // const ROOT::RVec<float>& phoCalibEt,
        const ROOT::RVec<float>& phoEt,
        const ROOT::RVec<float>& phoSigmaIEtaIEtaFull5x5,
        const ROOT::RVec<float>& phoHoverE
    );
}
#endif