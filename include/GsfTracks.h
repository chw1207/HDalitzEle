#ifndef GSFTRACKS_H_
#define GSFTRACKS_H_

#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4D.h"

namespace gsf{
    // create the dataframe with tracks related columns using the following functions
    ROOT::RDF::RNode DefineGSFColumns(ROOT::RDF::RNode df);


    // lable the main gsf tracks in reco::GsfTrack.
    ROOT::RVec<int> IsMainGSF(
        const Long64_t event,
        const int nGSFTrk,
        const ROOT::RVec<float>& gsfD0,
        const ROOT::RVec<float>& gsfDz,
        const int nEle,
        const ROOT::RVec<float>& eleD0,
        const ROOT::RVec<float>& eleDz
    );


    // Find the ambiguous gsf tracks of a electron
    //* Description:
    //*     1) The collection of GSF tracks is split into several chunks associated with the different electrons(ambiguous GSF tracks)
    //*     2) Starting from the main GSF track of an electron, I keep adding the rest GSF tracks as the ambiguous GSF tracks of this electron until meeting the main GSF track of the other electron.
    //*     3) No quality cuts applied
    //*     4) Each GSF track in the reco::GsfTrack should all pass the dEta and dPhi requirements described in the CMS reconstruction paper, otherwise, it will not be stored in the reco::GsfTrack.
    //*     5) https://iopscience.iop.org/article/10.1088/1748-0221/16/05/P05014
    ROOT::RVec<ROOT::RVec<int>> TrkEleAssociation(
        const int nGSFTrk,
        const ROOT::RVec<float>& gsfD0,
        const ROOT::RVec<float>& gsfDz,
        const int nEle,
        const ROOT::RVec<float>& eleD0,
        const ROOT::RVec<float>& eleDz,
        const ROOT::RVec<int>& isMainGSF // branch added by the gsf::IsMainGSF
    );


    // Function to calculate the number of ambiguous gsf tracks
    ROOT::RVec<int> CalcNGsfMatchToReco(
        const int nEle,
        const ROOT::RVec<ROOT::RVec<int>>& ambGSF // branch added by the gsf::TrkEleAssociation
    );


    // Find the main gsf track Idx
    ROOT::RVec<int> FindMainGSF(
        const int nEle,
        const ROOT::RVec<ROOT::RVec<int>>& ambGSF // branch added by the gsf::TrkEleAssociation
    );


    // Function to Find the second gsf track Idx (pT max)
    ROOT::RVec<int> FindSubGSF_PtMax(
        const int nEle,
        const ROOT::RVec<ROOT::RVec<int>>& ambGSF, // branch added by the gsf::TrkEleAssociation
        const ROOT::RVec<float>& gsfPt,
        const ROOT::RVec<int>& gsfCharge
    );


    // Function to Find the second gsf track Idx (dR min)
    ROOT::RVec<int> FindSubGSF_dRMin(
        const int nEle,
        const ROOT::RVec<ROOT::RVec<int>>& ambGSF, // branch added by the gsf::TrkEleAssociation
        const ROOT::RVec<float>& gsfEta,
        const ROOT::RVec<float>& gsfPhi,
        const ROOT::RVec<int>& gsfCharge
    );


    //! NEW: Function to Find the second gsf track Idx (dR min with several cuts on IP)
    ROOT::RVec<int> FindSubGSF_dRMinWithCuts(
        const int nEle,
        const ROOT::RVec<float>& eleSCEta,
        const ROOT::RVec<ROOT::RVec<int>>& ambGSF, // branch added by the gsf::TrkEleAssociation
        const ROOT::RVec<float>& gsfEta,
        const ROOT::RVec<float>& gsfPhi,
        const ROOT::RVec<int>& gsfCharge,
        const ROOT::RVec<float>& gsfD0,
        const ROOT::RVec<float>& gsfDz,
        const ROOT::RVec<int>& gsfMissHits
    );


    // Function to match the gsf variable to the electron based on the track index
    //* Description:
    //*     1) eleTrkPt = MatchIdex(eleTrkIdx, gsfPt)
    //*     2) eleSubTrkPt = MatchIdex(eleSubTrkIdx, gsfPt)
    ROOT::RVec<float> MatchIndexF(
        const int nEle,
        const ROOT::RVec<int>& eleTrkIdx, // branch added by the gsf::FindMainGSF or gsf::FindSubGSF_dRMin
        const ROOT::RVec<float>& var
    );
    ROOT::RVec<int> MatchIndexI(
        const int nEle,
        const ROOT::RVec<int>& eleTrkIdx, // branch added by the gsf::FindMainGSF or gsf::FindSubGSF_dRMin
        const ROOT::RVec<int>& var
    );


    // Function to get the pt ratio between electron trak and sub-track (if nGsfMatchToReco > 1)
    ROOT::RVec<float> GetTrkPtRatio(
        const int nEle,
        const ROOT::RVec<int>& nGsfMatchToReco, // branch added by the gsf::CalcNGsfMatchToReco
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& eleTrk1,
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& eleTrk2
    );


    // Function to get the deltaR between electron trak and sub-track (if nGsfMatchToReco > 1)
    ROOT::RVec<float> GetTrkdR(
        const int nEle,
        const ROOT::RVec<int>& nGsfMatchToReco, // branch added by the gsf::CalcNGsfMatchToReco
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& eleTrk1,
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& eleTrk2
    );


    // Function to get the deltaR between electron trak and sub-track (if nGsfMatchToReco > 1)
    ROOT::RVec<float> GetTrkPtSum(
        const int nEle,
        const ROOT::RVec<int>& nGsfMatchToReco, // branch added by the gsf::CalcNGsfMatchToReco
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& eleTrk1,
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& eleTrk2
    );


    // Function to get the deltaR between electron trak and sub-track
    ROOT::RVec<float> GetTrkRelPtRatio(
        const int nEle,
        const ROOT::RVec<float>& eleSCRawEn,
        const ROOT::RVec<int>& nGsfMatchToReco,
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& eleTrk1,
        const ROOT::RVec<ROOT::Math::PtEtaPhiMVector>& eleTrk2
    );


    // Function to find the basic cluster that the main gsf track points to
    ROOT::RVec<int> FindBC(
        const ROOT::RVec<int> eleTrkIdx,
        const ROOT::RVec<float>& gsfEta,
        const ROOT::RVec<float>& gsfPhi,
        const int nBC,
        const ROOT::RVec<float>& bcEta,
        const ROOT::RVec<float>& bcPhi
    );


    // Function to find the basic cluster that the sub gsf track points to
    ROOT::RVec<int> FindSubBC(
        const ROOT::RVec<int>& eleSubTrkIdx,
        const ROOT::RVec<float>& gsfEta,
        const ROOT::RVec<float>& gsfPhi,
        const int nBC,
        const ROOT::RVec<float>& bcEta,
        const ROOT::RVec<float>& bcPhi,
        const ROOT::RVec<int>& eleBCIdx // branch added by the gsf::FindBC
    );
}
#endif
