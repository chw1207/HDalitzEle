#ifndef ELECTRONSEL_H_
#define ELECTRONSEL_H_

#include "yaml-cpp/yaml.h"
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"

namespace eleSel{
    //  Fall17V2 MVA ID selection
    // https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/RecoEgamma/ElectronIdentification/python/Identification/mvaElectronID_Fall17_iso_V2_cff.py#L51-L59
    double MVARawCutWP90(float pt, int eleCase);
    ROOT::RVec<int> Fall17V2ID(
        const int nEle,
        const ROOT::RVec<float>& elePt,
        const ROOT::RVec<float>& eleSCEta,
        const ROOT::RVec<float>& eleIDMVAIso
    );


    // It requires the photon which is matched to the electron by SC to pass H->gg preselection
    // H->gg preselection: https://arxiv.org/pdf/1804.02716.pdf
    // https://github.com/cms-analysis/flashgg/blob/dev_legacy_runII/Taggers/python/flashggPreselectedDiPhotons_cfi.py
    // https://github.com/andreypz/nwu-dalitz-analysis/blob/master/zgamma/egamma.C#L493
    ROOT::RVec<int> HggPresel(
        const int nEle,
        const ROOT::RVec<float>& eleSCEta,
        const ROOT::RVec<float>& eleSCPhi,
        const int nPho,
        const float rhoAll,
        const ROOT::RVec<float>& phoSCEta,
        const ROOT::RVec<float>& phoSCPhi,
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



