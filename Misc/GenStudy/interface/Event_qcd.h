#include <TLorentzVector.h>
#include <TVector3.h>
#include <TString.h>

struct Event_qcd // per object
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

    Int_t iGen;
    Int_t iReco;
    vector<Int_t> GsfIdx;
    vector<Int_t> BCIdx;
};