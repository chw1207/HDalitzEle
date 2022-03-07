using namespace std;

void GenSelection(TreeReader &data, vector<Int_t> &v_genlep, vector<Int_t> &v_genpho, Int_t lepID, const char *proc)
{
    int nMC = data.GetInt("nMC"); // MC
    int *mcPID = data.GetPtrInt("mcPID");
    int *mcMomPID = data.GetPtrInt("mcMomPID");
    float *mcMomMass = data.GetPtrFloat("mcMomMass");
    int *mcGMomPID = data.GetPtrInt("mcGMomPID");
    float *mcPt = data.GetPtrFloat("mcPt");
    float *mcEta = data.GetPtrFloat("mcEta");
    float *mcPhi = data.GetPtrFloat("mcPhi");
    UShort_t *mcStatusFlag = (UShort_t *)data.GetPtrShort("mcStatusFlag");

    for (int i = 0; i < nMC; i++)
    {
        Bool_t target_lep = false;
        Bool_t target_pho = false;

        if (strcmp(proc, "HDalitz") == 0)
        {
            target_lep = ((abs(mcPID[i]) == 13 || abs(mcPID[i]) == 11) &&
                          mcMomPID[i] == 25 &&
                          ((mcStatusFlag[i] >> 0) & 1) == 1 && ((mcStatusFlag[i] >> 1) & 1) == 1) &&
                         (fabs(mcEta[i]) < 2.5);
            target_pho = (mcPID[i] == 22 && mcMomPID[i] == 25 &&
                          ((mcStatusFlag[i] >> 0) & 1) == 1 && ((mcStatusFlag[i] >> 1) & 1) == 1) &&
                         (fabs(mcEta[i]) < 2.5);

            if (target_lep)
                v_genlep.push_back(i);

            if (target_pho)
                v_genpho.push_back(i);
        }
        else if (strcmp(proc, "DYJetsToLL") == 0)
        {
            target_lep = ((abs(mcPID[i]) == 11) && mcMomPID[i] == 23 &&
                          ((mcStatusFlag[i] >> 0) & 1) == 1 && ((mcStatusFlag[i] >> 1) & 1) == 1) &&
                         (fabs(mcEta[i]) < 2.5);

            if (target_lep)
                v_genlep.push_back(i);

            if (v_genpho.size() < 1)
                v_genpho.push_back(-1);
        }
        else if (strcmp(proc, "gjets") == 0)
        {
            target_lep = ((abs(mcPID[i]) == 11) && mcMomPID[i] == 22 && fabs(mcEta[i]) < 2.5);
            // if (abs(mcPID[i]) == 11  && mcMomPID[i] == 22 && fabs(mcEta[i]) < 2.5)
            // {
            //     TLorentzVector genele;
            //     genele.SetPtEtaPhiM(mcPt[i], mcEta[i], mcPhi[i], 0.000510999);
            //     target_lep = gjet_convgam(data, genele);
            // }

            if (target_lep)
                v_genlep.push_back(i);

            if (v_genpho.size() < 1)
                v_genpho.push_back(-1);
        }
        else
        {
            cout << "Process [" << proc << "] is not available. Please check!" << endl;
            gSystem->Exit(0);
        }
    }
}