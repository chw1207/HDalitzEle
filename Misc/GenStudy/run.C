{
    gSystem->AddIncludePath("-Iexternal");
    gSystem->SetBuildDir("tmpdir", kTRUE);
    gROOT->ProcessLine(".L RecoLevel.C++");

    system("mkdir -p ./minitree/ ./minitree/2016/  ./minitree/2017/ ./minitree/2018/");

    system("rm yield_prod.txt");
    system("rm mcwei.txt");

    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_fall17_Dalitz_eeg_m125/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_ggF_eeg_m125_2017_RECO.root", "2017", "HDalitz", "ggF", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_fall17_Dalitz_eeg_VBF_m125/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_VBF_eeg_m125_2017_RECO.root", "2017", "HDalitz", "VBF", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_fall17_Dalitz_eeg_ZH_m125/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_ZH_eeg_m125_2017_RECO.root", "2017", "HDalitz", "ZH", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_fall17_Dalitz_eeg_WH_m125/ggtree_mc_*.root", "./minitree/2017/Minitree_HDalitz_WH_eeg_m125_2017_RECO.root", "2017", "HDalitz", "WH", 11);

    // DY+Jets
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_summer16_DYJetsToLL_m50_aMCatNLO_ext2/ggtree_mc_*.root", "./minitree/2016/Minitree_DYJetsToLL_2016_RECO.root", "2016", "DYJetsToLL", "", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_fall17_DYJetsToLL_m50_aMCatNLO*/ggtree_mc_*.root", "./minitree/2017/Minitree_DYJetsToLL_2017_RECO.root", "2017", "DYJetsToLL", "", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_autumn18_DYJetsToLL_m50_aMCatNLO*/ggtree_mc_*.root", "./minitree/2018/Minitree_DYJetsToLL_2018_RECO.root", "2018", "DYJetsToLL", "", 11);

    // Gamma+Jets
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_summer16_gjet_pt15to6000/ggtree_mc_*.root", "./minitree/2016/Minitree_gjet_pt15to6000_2016_RECO.root", "2016", "gjets", "pt15to6000", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_summer16_gjet_pt20_MGG_40to80/ggtree_mc_*.root", "./minitree/2016/Minitree_gjet_pt20_MGG_40to80_2016_RECO.root", "2016", "gjets", "pt20_MGG_40to80", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_summer16_gjet_pt20to40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2016/Minitree_gjet_pt20to40_MGG_80toInf_2016_RECO.root", "2016", "gjets", "pt20to40_MGG_80toInf", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_summer16_gjet_pt40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2016/Minitree_gjet_pt40_MGG_80toInf_2016_RECO.root", "2016", "gjets", "pt40_MGG_80toInf", 11);

    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_fall17_gjet_pt15to6000/ggtree_mc_*.root", "./minitree/2017/Minitree_gjet_pt15to6000_2017_RECO.root", "2017", "gjets", "pt15to6000", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_fall17_gjet_pt20_MGG_40to80/ggtree_mc_*.root", "./minitree/2017/Minitree_gjet_pt20_MGG_40to80_2017_RECO.root", "2017", "gjets", "pt20_MGG_40to80", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_fall17_gjet_pt20to40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2017/Minitree_gjet_pt20to40_MGG_80toInf_2017_RECO.root", "2017", "gjets", "pt20to40_MGG_80toInf", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_fall17_gjet_pt40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2017/Minitree_gjet_pt40_MGG_80toInf_2017_RECO.root", "2017", "gjets", "pt40_MGG_80toInf", 11);

    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_autumn18_gjet_pt15to6000/ggtree_mc_*.root", "./minitree/2018/Minitree_gjet_pt15to6000_2018_RECO.root", "2018", "gjets", "pt15to6000", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_autumn18_gjet_pt20_MGG_40to80/ggtree_mc_*.root", "./minitree/2018/Minitree_gjet_pt20_MGG_40to80_2018_RECO.root", "2018", "gjets", "pt20_MGG_40to80", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_autumn18_gjet_pt20to40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2018/Minitree_gjet_pt20to40_MGG_80toInf_2018_RECO.root", "2018", "gjets", "pt20to40_MGG_80toInf", 11);
    RecoLevel("/data6/ggNtuples/V10_02_10_05/job_autumn18_gjet_pt40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2018/Minitree_gjet_pt40_MGG_80toInf_2018_RECO.root", "2018", "gjets", "pt40_MGG_80toInf", 11);
}
