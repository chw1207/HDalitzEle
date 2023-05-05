void ExecMacro(std::string infile, std::string outfile){
    gROOT->ProcessLine(Form(".x mergeTnpNtuples.C+(\"%s\", \"%s\")", infile.c_str(), outfile.c_str()));
}


void runMergeTnpNtuples(std::string leg){
    // Setup the compile directory and library
    gSystem->SetBuildDir("build", true);
    const char* loc = gSystem->Getenv("HDalitzEle_LOC");
    gSystem->AddIncludePath(Form("-I%s/include", loc));
    gSystem->Load(Form("%s/lib/libHDalitzEle.so", loc));
    
    const char* ntuple_input_dir = (leg == "unseed") ? "/data4/chenghan/tnpTuples-unseedDiPhoHLT" : "/data4/chenghan/tnpTuples-seedDiPhoHLT";
    // UL2016preVFP
    gSystem->Exec(Form("mkdir -p %s/UL2016preVFP_merged", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016preVFP_Run2016B/*.root", ntuple_input_dir),    Form("%s/UL2016preVFP_merged/Run2016B.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016preVFP_Run2016C/*.root", ntuple_input_dir),    Form("%s/UL2016preVFP_merged/Run2016C.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016preVFP_Run2016D/*.root", ntuple_input_dir),    Form("%s/UL2016preVFP_merged/Run2016D.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016preVFP_Run2016E/*.root", ntuple_input_dir),    Form("%s/UL2016preVFP_merged/Run2016E.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016preVFP_Run2016F/*.root", ntuple_input_dir),    Form("%s/UL2016preVFP_merged/Run2016F.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016preVFP_DY_LO/*.root",    ntuple_input_dir),    Form("%s/UL2016preVFP_merged/DY_LO.root",    ntuple_input_dir));
    ExecMacro(Form("%s/UL2016preVFP_DY_NLO/*.root",   ntuple_input_dir),    Form("%s/UL2016preVFP_merged/DY_NLO.root",   ntuple_input_dir));

    // UL2016postVFP
    gSystem->Exec(Form("mkdir -p %s/UL2016postVFP_merged", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016postVFP_Run2016F/*.root", ntuple_input_dir),   Form("%s/UL2016postVFP_merged/Run2016F.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016postVFP_Run2016G/*.root", ntuple_input_dir),   Form("%s/UL2016postVFP_merged/Run2016G.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016postVFP_Run2016H/*.root", ntuple_input_dir),   Form("%s/UL2016postVFP_merged/Run2016H.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2016postVFP_DY_LO/*.root",    ntuple_input_dir),   Form("%s/UL2016postVFP_merged/DY_LO.root",    ntuple_input_dir));
    ExecMacro(Form("%s/UL2016postVFP_DY_NLO/*.root",   ntuple_input_dir),   Form("%s/UL2016postVFP_merged/DY_NLO.root",   ntuple_input_dir));

    // UL2017
    gSystem->Exec(Form("mkdir -p %s/UL2017_merged", ntuple_input_dir));
    ExecMacro(Form("%s/UL2017_Run2017B/*.root", ntuple_input_dir),   Form("%s/UL2017_merged/Run2017B.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2017_Run2017C/*.root", ntuple_input_dir),   Form("%s/UL2017_merged/Run2017C.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2017_Run2017D/*.root", ntuple_input_dir),   Form("%s/UL2017_merged/Run2017D.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2017_Run2017E/*.root", ntuple_input_dir),   Form("%s/UL2017_merged/Run2017E.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2017_Run2017F/*.root", ntuple_input_dir),   Form("%s/UL2017_merged/Run2017F.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2017_DY_LO/*.root",    ntuple_input_dir),   Form("%s/UL2017_merged/DY_LO.root",    ntuple_input_dir));
    ExecMacro(Form("%s/UL2017_DY_NLO/*.root",   ntuple_input_dir),   Form("%s/UL2017_merged/DY_NLO.root",   ntuple_input_dir));

    // UL2018
    gSystem->Exec(Form("mkdir -p %s/UL2018_merged", ntuple_input_dir));
    ExecMacro(Form("%s/UL2018_Run2018A/*.root", ntuple_input_dir),   Form("%s/UL2018_merged/Run2018A.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2018_Run2018B/*.root", ntuple_input_dir),   Form("%s/UL2018_merged/Run2018B.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2018_Run2018C/*.root", ntuple_input_dir),   Form("%s/UL2018_merged/Run2018C.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2018_Run2018D/*.root", ntuple_input_dir),   Form("%s/UL2018_merged/Run2018D.root", ntuple_input_dir));
    ExecMacro(Form("%s/UL2018_DY_LO/*.root",    ntuple_input_dir),   Form("%s/UL2018_merged/DY_LO.root",    ntuple_input_dir));
    ExecMacro(Form("%s/UL2018_DY_NLO/*.root",   ntuple_input_dir),   Form("%s/UL2018_merged/DY_NLO.root",   ntuple_input_dir));
}