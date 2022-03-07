int GenMatchInd(TLorentzVector reco, int nMC, ROOT::RVec<float>& mcPt, ROOT::RVec<float>& mcEta, ROOT::RVec<float>& mcPhi, ROOT::RVec<float>& mcMass){
    
    float tmpRat = 999.;
    int tmpInd = -1;
    
    TLorentzVector gen;
    for (int i = 0; i < nMC; i++){
        gen.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], mcMass[i]);
        
        if (gen.DeltaR(reco) > 0.1) continue;
        if (fabs((gen.Pt() / reco.Pt()) - 1.) < tmpRat){
            tmpRat = fabs((gen.Pt() / reco.Pt()) - 1.);
            tmpInd = i;
            continue;
        }
        else continue;
    }
    
    return tmpInd;
}

int GenType(unsigned short mcStatusFlag, int num){
    int tmp = 0;    
    if ((mcStatusFlag >> num & 1) == 1) 
            tmp = 1;
    
    return tmp;
}