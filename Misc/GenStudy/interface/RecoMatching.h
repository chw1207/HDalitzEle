// #include "/afs/cern.ch/work/h/hajheng/private/HDalitz/interface/Utilities.h"

std::map<string, Int_t> genreco_match(TreeReader &data, vector<TLorentzVector> v_genobj, Int_t &genreco_case)
{
    if (v_genobj[0].Pt() < v_genobj[1].Pt())
        std::reverse(v_genobj.begin(), v_genobj.end());

    Int_t nEle = data.GetInt("nEle");
    float *elePt = data.GetPtrFloat("elePt");
    float *eleEta = data.GetPtrFloat("eleEta");
    float *elePhi = data.GetPtrFloat("elePhi");
    float *eleSCEta = data.GetPtrFloat("eleSCEta");
    float *eleSCPhi = data.GetPtrFloat("eleSCPhi");

    std::map<string, int> GenReco_idx;
    GenReco_idx["GenLep1"] = -1;
    GenReco_idx["GenLep2"] = -1;

    vector<Int_t> v_matchidx;
    v_matchidx.clear();

    for (size_t iobj = 0; iobj < v_genobj.size(); iobj++)
    {
        int temp_index = -1;
        float temp_ptratio = 999.;
        for (int i = 0; i < nEle; i++)
        {
            TLorentzVector tmpreco;
            tmpreco.SetPtEtaPhiM(elePt[i], eleEta[i], elePhi[i], 0.000510999);

            // printf("nEle=%d; iEle=%d; (elePt,eleEta,elePhi)=(%f,%f,%f); (eleSCEta,eleSCPhi)=(%f,%f) ; DeltaR=%f; elePt/genPt=%f \n", nEle, i, elePt[i], eleEta[i], elePhi[i], eleSCEta[i], eleSCPhi[i], tmpreco.DeltaR(v_genobj[iobj]), elePt[i] / v_genobj[iobj].Pt());

            if (tmpreco.DeltaR(v_genobj[iobj]) > 0.1)
                continue;

            if (fabs((elePt[i] / v_genobj[iobj].Pt()) - 1.) < temp_ptratio)
            {
                temp_ptratio = fabs((elePt[i] / v_genobj[iobj].Pt()) - 1.);
                temp_index = i;
                continue;
            }
            else
                continue;
        }
        v_matchidx.push_back(temp_index);
    }
    GenReco_idx["GenLep1"] = v_matchidx[0];
    GenReco_idx["GenLep2"] = v_matchidx[1];

    if (GenReco_idx["GenLep1"] != GenReco_idx["GenLep2"] && GenReco_idx["GenLep1"] != -1 && GenReco_idx["GenLep2"] != -1)
        genreco_case = 1; // 2GEN-2Reco
    if (GenReco_idx["GenLep1"] == GenReco_idx["GenLep2"] && GenReco_idx["GenLep1"] != -1 && GenReco_idx["GenLep2"] != -1)
        genreco_case = 2; // 2GEN-1Reco
    if (GenReco_idx["GenLep1"] != GenReco_idx["GenLep2"] && (GenReco_idx["GenLep1"] == -1 || GenReco_idx["GenLep2"] == -1))
        genreco_case = 3; // 2GEN-1Match1N
    if (GenReco_idx["GenLep1"] == GenReco_idx["GenLep2"] && (GenReco_idx["GenLep1"] == -1 && GenReco_idx["GenLep2"] == -1))
        genreco_case = 4; // 2GEN-2N

    return GenReco_idx;
}

std::map<string, vector<Int_t>> recogsf_match(TreeReader &data, std::map<string, int> GenMatchReco_idx, Int_t &recogsf_case)
{
    std::map<string, vector<Int_t>> RecoGsf_idx;
    RecoGsf_idx["GenLep1"] = {};
    RecoGsf_idx["GenLep2"] = {};

    Int_t nGSFTrk = data.GetInt("nGSFTrk");
    float *gsfPt = data.GetPtrFloat("gsfPt");
    float *gsfEta = data.GetPtrFloat("gsfEta");
    float *gsfPhi = data.GetPtrFloat("gsfPhi");
    Int_t nEle = data.GetInt("nEle");
    float *elePt = data.GetPtrFloat("eleCalibPt");
    float *eleEta = data.GetPtrFloat("eleEta");
    float *elePhi = data.GetPtrFloat("elePhi");

    for (std::map<string, int>::iterator it = GenMatchReco_idx.begin(); it != GenMatchReco_idx.end(); ++it)
    {
        if (it->second != -1)
        {
            TLorentzVector tmpreco;
            tmpreco.SetPtEtaPhiM(elePt[it->second], eleEta[it->second], elePhi[it->second], 0.000510999);

            vector<Int_t> v_matchedidx;
            v_matchedidx.clear();
            for (int i = 0; i < nGSFTrk; i++)
            {
                // printf("nGSFTrk=%d; igsf=%d; (gsfPt,gsfEta,fsgPhi)=(%f,%f,%f); (recoPt,recoEta,recoPhi)=(%f,%f,%f); DeltaR=%f; elePt/ObjPt=%f \n", nGSFTrk, i, gsfPt[i], gsfEta[i], gsfPhi[i], Obj.Pt(), Obj.Eta(), Obj.Phi(), deltaR(gsfEta[i], gsfPhi[i], Obj.Eta(), Obj.Phi()), gsfPt[i] / Obj.Pt());

                // if (deltaR(gsfEta[i], gsfPhi[i], tmpreco.Eta(), tmpreco.Phi()) > 0.1)
                //     continue;

                if (gsfPt[i]/tmpreco.Pt() > 5) continue; // remove the gsf track with unreaonable high Pt

                // Reference: https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/RecoEgamma/EgammaElectronProducers/plugins/GsfElectronProducer.cc#L188-L191
                
                if (fabs(deltaEta(gsfEta[i], tmpreco.Eta())) > 0.02 ||
                    fabs(deltaPhi(gsfPhi[i], tmpreco.Phi())) > 0.15)
                    continue;

                RecoGsf_idx[it->first].push_back(i);
            }
        }
    }

    // 2Reco2Gsf
    if (GenMatchReco_idx["GenLep1"] != -1 && GenMatchReco_idx["GenLep2"] != -1 && GenMatchReco_idx["GenLep1"] != GenMatchReco_idx["GenLep2"] &&
        RecoGsf_idx["GenLep1"].size() >= 1 && RecoGsf_idx["GenLep2"].size() >= 1)
        recogsf_case = 1;
    // 1Reco2Gsf
    if (GenMatchReco_idx["GenLep1"] != -1 && GenMatchReco_idx["GenLep2"] != -1 && GenMatchReco_idx["GenLep1"] == GenMatchReco_idx["GenLep2"] &&
        RecoGsf_idx["GenLep1"].size() >= 2)
        recogsf_case = 2;
    // 1Reco1Gsf (either from 2Gen1Reco or from 2Gen1Match1N)
    if (((GenMatchReco_idx["GenLep1"] != -1 && GenMatchReco_idx["GenLep2"] != -1) && GenMatchReco_idx["GenLep1"] == GenMatchReco_idx["GenLep2"] &&
         (RecoGsf_idx["GenLep1"].size() == 1)) ||
        ((GenMatchReco_idx["GenLep1"] == -1 || GenMatchReco_idx["GenLep2"] == -1) && GenMatchReco_idx["GenLep1"] != GenMatchReco_idx["GenLep2"] &&
         ((RecoGsf_idx["GenLep1"].size() >= 1 || RecoGsf_idx["GenLep2"].size() >= 1))))
        recogsf_case = 3;
    // 0Reco0Gsf
    if ((GenMatchReco_idx["GenLep1"] == -1 && GenMatchReco_idx["GenLep2"] == -1) &&
        (RecoGsf_idx["GenLep1"].size() == 0 && RecoGsf_idx["GenLep2"].size() == 0))
        recogsf_case = 4;

    return RecoGsf_idx;
}

std::map<string, Int_t> gsfgen_match(TreeReader &data, vector<TLorentzVector> v_genobj, std::map<string, vector<Int_t>> RecoMatchGsf_idx, Int_t &gsfgen_case)
{
    if (v_genobj[0].Pt() < v_genobj[1].Pt())
        std::reverse(v_genobj.begin(), v_genobj.end());

    Int_t nGSFTrk = data.GetInt("nGSFTrk");
    float *gsfPt = data.GetPtrFloat("gsfPt");
    float *gsfEta = data.GetPtrFloat("gsfEta");
    float *gsfPhi = data.GetPtrFloat("gsfPhi");

    std::map<string, int> GsfGen_idx;
    GsfGen_idx["GenLep1"] = -1;
    GsfGen_idx["GenLep2"] = -1;
    vector<string> maplabel{"GenLep1", "GenLep2"};

    vector<Int_t> v_matchgsfidx;
    v_matchgsfidx.clear();

    for (size_t iobj = 0; iobj < v_genobj.size(); iobj++)
    {
        int temp_index = -1;
        float temp_ptratio = 1E10;
        for (size_t i = 0; i < RecoMatchGsf_idx[maplabel[iobj]].size(); i++)
        {
            // printf("nGSFTrk=%d; igsf=%d; igenobj=%zu; (gsfPt,gsfEta,gsfPhi)=(%f,%f,%f); (genPt,genEta,genPhi)=(%f,%f,%f); DeltaR=%f; elePt/ObjPt=%f \n", nGSFTrk, i, iobj, gsfPt[i], gsfEta[i], gsfPhi[i], v_genobj[iobj].Pt(), v_genobj[iobj].Eta(), v_genobj[iobj].Phi(), deltaR(gsfEta[i], gsfPhi[i], v_genobj[iobj].Eta(), v_genobj[iobj].Phi()), gsfPt[i] / v_genobj[iobj].Pt());

            bool repeat = false;
            for (size_t ii = 0; ii < v_matchgsfidx.size(); ii++)
            {
                if (RecoMatchGsf_idx[maplabel[iobj]][i] == v_matchgsfidx[ii])
                    repeat = true;
            }
            if (repeat)
                continue;

            if (deltaR(gsfEta[RecoMatchGsf_idx[maplabel[iobj]][i]], gsfPhi[RecoMatchGsf_idx[maplabel[iobj]][i]], v_genobj[iobj].Eta(), v_genobj[iobj].Phi()) > 0.1)
                continue;

            if (fabs((gsfPt[RecoMatchGsf_idx[maplabel[iobj]][i]] / v_genobj[iobj].Pt()) - 1.) < temp_ptratio)
            {
                temp_ptratio = fabs((gsfPt[RecoMatchGsf_idx[maplabel[iobj]][i]] / v_genobj[iobj].Pt()) - 1.);
                temp_index = RecoMatchGsf_idx[maplabel[iobj]][i];
                continue;
            }
            else
                continue;
        }
        v_matchgsfidx.push_back(temp_index);
    }

    if (v_matchgsfidx[0] != v_matchgsfidx[1] && v_matchgsfidx[0] != -1 && v_matchgsfidx[1] != -1)
    {
        if (gsfPt[v_matchgsfidx[0]] < gsfPt[v_matchgsfidx[1]] &&
            deltaR(gsfEta[v_matchgsfidx[0]], gsfPhi[v_matchgsfidx[0]], gsfEta[v_matchgsfidx[1]], gsfPhi[v_matchgsfidx[1]]) < 0.1)
        {
            // When the pT of matched gsf of GenLep1 < that of GenLep2 AND the deltaR between the two gsf tracks is less than 0.1, swap the matched gsf tracks
            Int_t temp = v_matchgsfidx[0];
            v_matchgsfidx[0] = v_matchgsfidx[1];
            v_matchgsfidx[1] = temp;
        }
    }

    GsfGen_idx["GenLep1"] = v_matchgsfidx[0];
    GsfGen_idx["GenLep2"] = v_matchgsfidx[1];

    // 2Gsf2Gen
    if ((RecoMatchGsf_idx["GenLep1"].size() >= 1 && RecoMatchGsf_idx["GenLep2"].size() >= 1) &&
        (GsfGen_idx["GenLep1"] != GsfGen_idx["GenLep2"]) && GsfGen_idx["GenLep1"] != -1 && GsfGen_idx["GenLep2"] != -1)
        gsfgen_case = 1;
    // 1Gsf2Gen
    if ((RecoMatchGsf_idx["GenLep1"].size() >= 1 || RecoMatchGsf_idx["GenLep2"].size() >= 1) &&
        (GsfGen_idx["GenLep1"] != GsfGen_idx["GenLep2"]) && (GsfGen_idx["GenLep1"] == -1 || GsfGen_idx["GenLep2"] == -1))
        gsfgen_case = 2;
    // 0Gsf2Gen
    if (((RecoMatchGsf_idx["GenLep1"].size() == 0 && RecoMatchGsf_idx["GenLep2"].size() == 0) &&
         (GsfGen_idx["GenLep1"] == GsfGen_idx["GenLep2"]) && (GsfGen_idx["GenLep1"] == -1 && GsfGen_idx["GenLep2"] == -1)) ||
        ((RecoMatchGsf_idx["GenLep1"].size() > 0 || RecoMatchGsf_idx["GenLep2"].size() > 0) &&
         (GsfGen_idx["GenLep1"] == GsfGen_idx["GenLep2"]) && (GsfGen_idx["GenLep1"] == -1 && GsfGen_idx["GenLep2"] == -1)))
        gsfgen_case = 3;

    return GsfGen_idx;
}

// For debug
Int_t Category(Int_t genreco_case, Int_t recogsf_case, Int_t gsfgen_case, const char *&Category_str)
{
    Int_t category = -1;
    string GenReco_str = "None", RecoGsf_str = "None", GsfGen_str = "None", cat_str = "None";
    if (genreco_case == 1)
        GenReco_str = "2Gen2Reco";
    if (genreco_case == 2)
        GenReco_str = "2Gen1Reco";
    if (genreco_case == 3)
        GenReco_str = "2Gen1Match1N";
    if (genreco_case == 4)
        GenReco_str = "2Gen2N";

    if (recogsf_case == 1)
        RecoGsf_str = "2Reco2Gsf";
    if (recogsf_case == 2)
        RecoGsf_str = "1Reco2Gsf";
    if (recogsf_case == 3)
        RecoGsf_str = "1Reco1Gsf";
    if (recogsf_case == 4)
        RecoGsf_str = "0Reco0Gsf";

    if (gsfgen_case == 1)
        GsfGen_str = "2Gsf2Gen";
    if (gsfgen_case == 2)
        GsfGen_str = "1Gsf2Gen";
    if (gsfgen_case == 3)
        GsfGen_str = "0Gsf2Gen";

    if (genreco_case == 1 && recogsf_case == 1 && gsfgen_case == 1)
    {
        category = 1;
        cat_str = "Resolved";
    }
    else if (genreco_case == 2 && recogsf_case == 2 && gsfgen_case == 1)
    {
        category = 2;
        cat_str = "Merged-2Gsf";
    }
    //! Modified: 
    else if ((genreco_case == 2 && recogsf_case == 3 && gsfgen_case == 2)) 
            //   || (genreco_case == 2 && recogsf_case == 2 && gsfgen_case == 2))
    {
        category = 3;
        cat_str = "Merged-1 Missing track";
    }
    else
    {
        category = 4;
        cat_str = "Not properly reconstructed";
    }

    Category_str = Form("%s->%s->%s==>%s", GenReco_str.c_str(), RecoGsf_str.c_str(), GsfGen_str.c_str(), cat_str.c_str());

    return category;
}

void Category_count(Int_t genreco_case, Int_t recogsf_case, Int_t gsfgen_case, std::map<string, Int_t> &count_case, std::map<string, Int_t> &count_category)
{
    if (genreco_case == 1 && recogsf_case == 1 && gsfgen_case == 1)
    {
        count_case["2Gen2Reco->2Reco2Gsf->2Gsf2Gen"] += 1;
        count_category["Resolved"] += 1;
    }
    else if (genreco_case == 2 && recogsf_case == 2 && gsfgen_case == 1)
    {
        count_case["2Gen1Reco->1Reco2Gsf->2Gsf2Gen"] += 1;
        count_category["Merged2Gsf"] += 1;
    }
    else if ((genreco_case == 2 && recogsf_case == 3 && gsfgen_case == 2) ||
             (genreco_case == 2 && recogsf_case == 2 && gsfgen_case == 2))
    {
        count_category["Merged1MissingGsf"] += 1;
        count_case["2Gen1Reco->1Reco1Gsf->1Gsf2Gen"] += 1;
        // if (genreco_case == 2 && recogsf_case == 3 && gsfgen_case == 2)
        //     count_case["2Gen1Reco->1Reco1Gsf->1Gsf2Gen"] += 1;
        // if (genreco_case == 2 && recogsf_case == 2 && gsfgen_case == 2)
        //     count_case["2Gen1Reco->1Reco2Gsf->1Gsf2Gen"] += 1;
    }
    else
    {
        count_category["NPR"] += 1;
        count_case["Others"] += 1;
    }
}

void debuginfo(TreeReader &data, Long64_t ev, vector<TLorentzVector> v_genlepobj, std::map<string, Int_t> GenMatchReco_idx, std::map<string, vector<Int_t>> RecoMatchGsf_idx, std::map<string, Int_t> GenMatchGSF_idx, Int_t GenReco_case, Int_t RecoGsf_case, Int_t GsfGen_case)
{
    // GEN information
    int nMC = data.GetInt("nMC"); // MC
    float *mcPt = data.GetPtrFloat("mcPt");
    float *mcEta = data.GetPtrFloat("mcEta");
    float *mcPhi = data.GetPtrFloat("mcPhi");
    // RECO information
    Int_t nEle = data.GetInt("nEle");
    float *elePt = data.GetPtrFloat("elePt");
    float *eleEta = data.GetPtrFloat("eleEta");
    float *elePhi = data.GetPtrFloat("elePhi");
    float *eleSCEta = data.GetPtrFloat("eleSCEta");
    float *eleSCPhi = data.GetPtrFloat("eleSCPhi");
    // Gsf track information
    Int_t nGSFTrk = data.GetInt("nGSFTrk");
    float *gsfPt = data.GetPtrFloat("gsfPt");
    float *gsfEta = data.GetPtrFloat("gsfEta");
    float *gsfPhi = data.GetPtrFloat("gsfPhi");

    printf("ev=%lld\nGenlep1 (pT,eta,phi=%fGeV,%f,%f) matches to RECO %d; Genlep2 (pT,eta,phi=%fGeV,%f,%f) matches to RECO %d \n", ev, v_genlepobj[0].Pt(), v_genlepobj[0].Eta(), v_genlepobj[0].Phi(), GenMatchReco_idx["GenLep1"], v_genlepobj[1].Pt(), v_genlepobj[1].Eta(), v_genlepobj[1].Phi(), GenMatchReco_idx["GenLep2"]);
    printf("nEle=%d, nGSFTrk=%d, RECO %d matched gsf size = %zd; RECO %d matched gsf size = %zd \n", nEle, nGSFTrk, GenMatchReco_idx["GenLep1"], RecoMatchGsf_idx["GenLep1"].size(), GenMatchReco_idx["GenLep2"], RecoMatchGsf_idx["GenLep2"].size());

    if (GenReco_case == 1)
    {
        TLorentzVector reco1, reco2, reco1_SC, reco2_SC;
        reco1.SetPtEtaPhiM(elePt[GenMatchReco_idx["GenLep1"]], eleEta[GenMatchReco_idx["GenLep1"]], elePhi[GenMatchReco_idx["GenLep1"]], genlepmass(11));
        reco2.SetPtEtaPhiM(elePt[GenMatchReco_idx["GenLep2"]], eleEta[GenMatchReco_idx["GenLep2"]], elePhi[GenMatchReco_idx["GenLep2"]], genlepmass(11));
        reco1_SC.SetPtEtaPhiM(elePt[GenMatchReco_idx["GenLep1"]], eleSCEta[GenMatchReco_idx["GenLep1"]], eleSCPhi[GenMatchReco_idx["GenLep1"]], genlepmass(11));
        reco2_SC.SetPtEtaPhiM(elePt[GenMatchReco_idx["GenLep2"]], eleSCEta[GenMatchReco_idx["GenLep2"]], eleSCPhi[GenMatchReco_idx["GenLep2"]], genlepmass(11));
        printf("dR(genlep1,genlep2)=%f, dR(reco1,reco2)=%f, dR(reco1_SC,reco2_SC)=%f\n", v_genlepobj[0].DeltaR(v_genlepobj[1]), reco1.DeltaR(reco2), reco1_SC.DeltaR(reco2_SC));
    }
    Int_t igenlep = 0;
    for (std::map<string, vector<Int_t>>::iterator it = RecoMatchGsf_idx.begin(); it != RecoMatchGsf_idx.end(); ++it)
    {
        for (size_t i = 0; i < RecoMatchGsf_idx[it->first].size(); i++)
        {
            printf("%s's RECO matches to gsf (pT,eta,phi)=(%f,%f,%f), dR(%s,gsf)=%f\n", it->first.c_str(), gsfPt[RecoMatchGsf_idx[it->first][i]], gsfEta[RecoMatchGsf_idx[it->first][i]], gsfPhi[RecoMatchGsf_idx[it->first][i]], it->first.c_str(), deltaR(v_genlepobj[igenlep].Eta(), v_genlepobj[igenlep].Phi(), gsfEta[RecoMatchGsf_idx[it->first][i]], gsfPhi[RecoMatchGsf_idx[it->first][i]]));
        }
        igenlep++;
    }

    for (std::map<string, Int_t>::iterator it = GenMatchGSF_idx.begin(); it != GenMatchGSF_idx.end(); ++it)
    {
        printf("%s matches to gsf %d (pT,eta,phi)=(%f,%f,%f)\n", it->first.c_str(), GenMatchGSF_idx[it->first], gsfPt[GenMatchGSF_idx[it->first]], gsfEta[GenMatchGSF_idx[it->first]], gsfPhi[GenMatchGSF_idx[it->first]]);
    }

    const char *Category_str = "";
    Int_t category = Category(GenReco_case, RecoGsf_case, GsfGen_case, Category_str);

    printf("Category label: %s\n", Category_str);
    printf("----------------------\n");
}
