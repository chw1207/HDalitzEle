#include "ROOT/RVec.hxx"


ROOT::RVec<int> PassMuonID(const ROOT::RVec<int>& muIDbit, const int num){
    // different ID bit number corresponds to different muon ID
    // https://github.com/cmkuo/ggAnalysis/blob/106X/ggNtuplizer/plugins/ggNtuplizer_muons.cc#L162-L183

    ROOT::RVec<int> vec;
    vec.clear();

    for (int i = 0; i < muIDbit.size(); ++i) {
        if ((muIDbit[i] >> num & 1) == 1)
            vec.push_back(1);
        else
            vec.push_back(0);
    }

    return vec;
}

