#include "ROOT/RVec.hxx"

using namespace std;
using namespace ROOT::VecOps;

namespace hp{
    template <typename T>
    bool HasMC(T &df){
        vector<string> colNames = df.GetColumnNames();
        if (count(colNames.begin(), colNames.end(), "nMC")) return true;
        else return false;
    }

    // Function to get the index vector sorted by pT
    // * Reference: https://root.cern/doc/master/vo006__IndexManipulation_8C.html
    ROOT::RVec<int> getIdx(ROOT::RVec<int>& isgood, ROOT::RVec<float>& pt){
        ROOT::RVec<int> idx_select = Nonzero(isgood);
        ROOT::RVec<int> idx_sort = Reverse(Argsort(pt));
        ROOT::RVec<int> idx = Intersect(idx_sort, idx_select);

        return idx;
    }

    // Function to get the square root of quadrature sum
    float Mod(float a, float b, float c){
        return sqrt(a*a + b*b + c*c);
    }

    // Function to get the square root of quadrature sum
    float Mod(float a, float b){
        return sqrt(a*a + b*b);
    }

    // Function to construct TLorentzVector
    TLorentzVector P4Mass(float pt, float eta, float phi, float m){
        TLorentzVector v;
        v.SetPtEtaPhiM(pt, eta, phi, m);
        return v;
    }

    // Function to construct TLorentzVector
    TLorentzVector P4En(float pt, float eta, float phi, float e){
        TLorentzVector v;
        v.SetPtEtaPhiE(pt, eta, phi, e);
        return v;
    }
}