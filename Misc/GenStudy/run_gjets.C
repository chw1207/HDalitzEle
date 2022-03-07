{
    gSystem->AddIncludePath("-Iexternal");
    gSystem->SetBuildDir("tmpdir", kTRUE);
    gROOT->ProcessLine(".L RecoLevel_gjets.C++");

    system("mkdir -p ./minitree/ ./minitree/2016/  ./minitree/2017/ ./minitree/2018/");

    // Gamma+Jets
    RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_summer16_gjet_pt15to6000/ggtree_mc_*.root", "./minitree/2016/Minitree_gjet_pt15to6000_2016_RECO.root", "2016", "gjets", "pt15to6000", 11);
    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_summer16_gjet_pt20_MGG_40to80/ggtree_mc_*.root", "./minitree/2016/Minitree_gjet_pt20_MGG_40to80_2016_RECO.root", "2016", "gjets", "pt20_MGG_40to80", 11);
    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_summer16_gjet_pt20to40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2016/Minitree_gjet_pt20to40_MGG_80toInf_2016_RECO.root", "2016", "gjets", "pt20to40_MGG_80toInf", 11);
    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_summer16_gjet_pt40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2016/Minitree_gjet_pt40_MGG_80toInf_2016_RECO.root", "2016", "gjets", "pt40_MGG_80toInf", 11);

    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_fall17_gjet_pt15to6000/ggtree_mc_*.root", "./minitree/2017/Minitree_gjet_pt15to6000_2017_RECO.root", "2017", "gjets", "pt15to6000", 11);
    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_fall17_gjet_pt20_MGG_40to80/ggtree_mc_*.root", "./minitree/2017/Minitree_gjet_pt20_MGG_40to80_2017_RECO.root", "2017", "gjets", "pt20_MGG_40to80", 11);
    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_fall17_gjet_pt20to40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2017/Minitree_gjet_pt20to40_MGG_80toInf_2017_RECO.root", "2017", "gjets", "pt20to40_MGG_80toInf", 11);
    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_fall17_gjet_pt40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2017/Minitree_gjet_pt40_MGG_80toInf_2017_RECO.root", "2017", "gjets", "pt40_MGG_80toInf", 11);

    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_autumn18_gjet_pt15to6000/ggtree_mc_*.root", "./minitree/2018/Minitree_gjet_pt15to6000_2018_RECO.root", "2018", "gjets", "pt15to6000", 11);
    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_autumn18_gjet_pt20_MGG_40to80/ggtree_mc_*.root", "./minitree/2018/Minitree_gjet_pt20_MGG_40to80_2018_RECO.root", "2018", "gjets", "pt20_MGG_40to80", 11);
    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_autumn18_gjet_pt20to40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2018/Minitree_gjet_pt20to40_MGG_80toInf_2018_RECO.root", "2018", "gjets", "pt20to40_MGG_80toInf", 11);
    // RecoLevel_gjets("/data6/ggNtuples/V10_02_10_05/job_autumn18_gjet_pt40_MGG_80toInf/ggtree_mc_*.root", "./minitree/2018/Minitree_gjet_pt40_MGG_80toInf_2018_RECO.root", "2018", "gjets", "pt40_MGG_80toInf", 11);



}
