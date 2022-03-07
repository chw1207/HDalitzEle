auto MuTypeVector(ROOT::RVec<int>& muType, int num){
    ROOT::RVec<int> vec;
    vec.clear();
    const auto size = muType.size();
    for (auto i = 0; i < size; ++i) {
        if ((muType[i] >> num & 1) == 1) 
            vec.push_back(1);
        else 
            vec.push_back(0);
    }
    return vec;
}

auto NoneZeroIso03(ROOT::RVec<float>& muPFChIso03, ROOT::RVec<float>& muPFNeuIso03, ROOT::RVec<float>& muPFPhoIso03, ROOT::RVec<float>& muPFPUIso03){
    ROOT::RVec<float> vec;
    vec.clear();
    const auto size = muPFChIso03.size();
    for (auto i = 0; i < size; ++i) {
            vec.push_back(muPFChIso03[i] + TMath::Max(0., muPFNeuIso03[i] + muPFPhoIso03[i] - 0.5*muPFPUIso03[i]));
    }
    return vec;
}