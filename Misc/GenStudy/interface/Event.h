#include <TLorentzVector.h>
#include <TVector3.h>
#include <TString.h>

struct Event
{
    Float_t procXS;
    Float_t mcwei;
    Float_t genwei;
    ULong64_t HLTEleMuX;
    ULong64_t HLTEleMuXIsPrescaled;
    ULong64_t HLTPho;
    ULong64_t HLTPhoIsPrescaled;
    Float_t rho;
    Float_t rhoCentral;
    Int_t nVtx;
    Int_t nGoodVtx;
    Bool_t isPVGood;
    Long64_t totalEvents;
    Float_t instwei;

    Int_t category;
    Int_t GenReco_case;
    Int_t RecoGsf_case;
    Int_t GsfGen_case;

    Int_t iGenLep[2];
    Int_t iGenPho;
    Int_t iele[2];
    Int_t igsf[2];
    vector<Int_t> GsfIdxLep1; // Gsf tracks to which Reco1 matches
    vector<Int_t> GsfIdxLep2; // Gsf tracks to which Reco2 matches

    vector<Int_t> BCIdxLep1; // BC to which Reco1 matches
    vector<Int_t> BCIdxLep2; // BC to which Reco2 matches

    Int_t ConvIdxLep1; // Conversion to which Reco1 matches
    Int_t ConvIdxLep2; // Conversion to which Reco2 matches

};