#include <iostream>
#include <vector>
#include <string>
#include <random> 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <TH1.h>
#include <TFile.h>
#include <TRandom3.h>


class MyEnCalibrater{
public:
    MyEnCalibrater(){};

    ~MyEnCalibrater(){
        // if (hasScorr){
        //     delete hscorrEB;
        //     delete hscorrEE;
        // }
        if (hasScale){
            delete hscaleEB;
            delete hscaleEE;
        }
        if (hasSmear){
            delete hsmearEB;
            delete hsmearEE;
            delete rnd;
        }
    }


    // void InitScorr(std::vector<std::string> scorr_files){
    //     assert(scorr_files.size() == 2);

    //     hscorrEB = getHist(scorr_files[0], "hscorr");
    //     hscorrEE = getHist(scorr_files[1], "hscorr");
    //     hasScorr = true;
    // }


    void InitScale(std::vector<std::string> scale_files){
        assert(scale_files.size() == 2);

        hscaleEB = getHist(scale_files[0], "hscorr");
        hscaleEE = getHist(scale_files[1], "hscorr");
        hasScale = true;
    }


    void InitSmear(std::vector<std::string> smear_files){
        assert(smear_files.size() == 2);

        hsmearEB = getHist(smear_files[0], "hscorr");
        hsmearEE = getHist(smear_files[1], "hscorr");
        rnd = new TRandom3(rd());
        hasSmear = true;
    }


    // ROOT::RDF::RNode GetScorrRDF(ROOT::RDF::RNode df, std::string column, std::string column_eta, std::string column_corr){
    //     auto nf = df.Define(column_corr, [this](const float sceta, const float pt){
    //         bool isEB = fabs(sceta) < 1.4442;
    //         float SF = 1.;
    //         if (isEB){
    //             int bin = hscorrEB->GetXaxis()->FindBin(pt);
    //             int Nbin = hscorrEB->GetNbinsX();
    //             if (bin > Nbin)
    //                 bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
    //             if (bin < 1)
    //                 bin = 1;
    //             SF = hscorrEB->GetBinContent(bin);
    //         }
    //         else{
    //             int bin = hscorrEE->GetXaxis()->FindBin(pt);
    //             int Nbin = hscorrEE->GetNbinsX();
    //             if (bin > Nbin)
    //                 bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
    //             if (bin < 1)
    //                 bin = 1;
    //             SF = hscorrEE->GetBinContent(bin);
    //         }

    //         return (float) pt * SF;
    //     }, {column_eta, column});
        
    //     return nf;
    // }

    ROOT::RDF::RNode GetScaleRDF(ROOT::RDF::RNode df, bool isMC, std::string column, std::string column_eta, std::string column_corr){
        if (isMC)
            return df.Define(column_corr, column);
        
        auto nf = df.Define(column_corr, [this](const float sceta, const float pt){
            bool isEB = fabs(sceta) < 1.479;
            float SF = 1.;
            if (isEB){
                int bin = hscaleEB->GetXaxis()->FindBin(pt);
                int Nbin = hscaleEB->GetNbinsX();
                if (bin > Nbin)
                    bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (bin < 1)
                    bin = 1;
                SF = hscaleEB->GetBinContent(bin);
            }
            else{
                int bin = hscaleEE->GetXaxis()->FindBin(pt);
                int Nbin = hscaleEE->GetNbinsX();
                if (bin > Nbin)
                    bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (bin < 1)
                    bin = 1;
                SF = hscaleEE->GetBinContent(bin);
            }

            return (float) pt * SF;
        }, {column_eta, column});

        return nf;
    }

    ROOT::RDF::RNode GetSmearRDF(ROOT::RDF::RNode df, bool isMC, std::string column, std::string column_eta, std::string column_corr){
        if (!isMC)
            return df.Define(column_corr, column);
        
        auto nf = df.Define(column_corr, [this](const float sceta, const float pt){
            bool isEB = fabs(sceta) < 1.4442;
            float SF = 1.;
            if (isEB){
                int bin = hsmearEB->GetXaxis()->FindBin(pt);
                int Nbin = hsmearEB->GetNbinsX();
                if (bin > Nbin)
                    bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (bin < 1)
                    bin = 1;
                SF = rnd->Gaus(1, hsmearEB->GetBinContent(bin));
            }
            else{
                int bin = hsmearEE->GetXaxis()->FindBin(pt);
                int Nbin = hsmearEE->GetNbinsX();
                if (bin > Nbin)
                    bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (bin < 1)
                    bin = 1;
                SF = rnd->Gaus(1, hsmearEE->GetBinContent(bin));
            }

            return (float) pt * SF;
        }, {column_eta, column});

        return nf;
    }


private:
    TH1F* getHist(std::string histFile, std::string histName){
        TFile* fhist = new TFile(histFile.c_str(), "READ");
        assert(fhist);
        
        TH1F*  hist  = (TH1F*) fhist->Get(histName.c_str());
        assert(hist);

        hist->SetDirectory(nullptr);
        fhist->Close();
        delete fhist;
        
        return hist;
    }

    // for energy scale correction
    // bool hasScorr = false;
    // TH1F* hscorrEB = nullptr;
    // TH1F* hscorrEE = nullptr;
    
    // for residual correction of energy scale 
    bool hasScale = false;
    TH1F* hscaleEB = nullptr;
    TH1F* hscaleEE = nullptr;
    
    // for residual correction of energy resolution 
    std::random_device rd;
    bool hasSmear = false;
    TH1F* hsmearEB = nullptr;
    TH1F* hsmearEE = nullptr;
    TRandom3* rnd = nullptr;
};





















































































// // https://root-forum.cern.ch/t/adding-data-from-an-external-container-to-a-dataframe/46177/5
// // https://lbianch.github.io/programming/2016/12/19/root-histograms-uniqueptr.html
// ROOT::RDF::RNode EnConvCorrector(ROOT::RDF::RNode df){
//     auto fcorrEB = std::make_shared<TFile>("../data/EnConvCorr_combined_EB_nominal.root", "READ");
//     std::shared_ptr<TH1F> hcorrEB((TH1F*)fcorrEB->Get("hcorr"));
//     hcorrEB->SetDirectory(nullptr);

//     auto fcorrEE = std::make_shared<TFile>("../data/EnConvCorr_combined_EE_nominal.root", "READ");
//     std::shared_ptr<TH1F> hcorrEE((TH1F*)fcorrEE->Get("hcorr"));
//     hcorrEE->SetDirectory(nullptr);

//     auto nf = df.Define("phoConvCorrEt_Lead",
//         [hcorrEB, hcorrEE](
//             const float& phoSCEta,
//             const float& phoCalibEt
//         ){

//             bool isEB = abs(phoSCEta) < 1.4442;
//             float SF = 1.;
//             if (isEB){
//                 int bin = hcorrEB->GetXaxis()->FindBin(phoCalibEt);
//                 int Nbin = hcorrEB->GetNbinsX();
//                 if (bin > Nbin || bin < 1)
//                     bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin

//                 SF = hcorrEB->GetBinContent(bin);
//             }
//             else{
//                 int bin = hcorrEE->GetXaxis()->FindBin(phoCalibEt);
//                 int Nbin = hcorrEE->GetNbinsX();
//                 if (bin > Nbin || bin < 1)
//                     bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin

//                 SF = hcorrEE->GetBinContent(bin);
//             }

//             return (float) phoCalibEt * SF;
//         }, {"phoSCEta_Lead", "phoCalibEt_Lead"});

//     return nf;
// }


// ROOT::RDF::RNode EnScaleCorrector(ROOT::RDF::RNode df, bool isMC){
//     if (isMC)
//         return df.Define("phoScaleCorrEt_Lead", "phoConvCorrEt_Lead");

//     auto fcorrEB = std::make_shared<TFile>("../data/EnScaleCorr_combined_EB.root", "READ");
//     std::shared_ptr<TH1F> hcorrEB((TH1F*)fcorrEB->Get("hscorr"));
//     hcorrEB->SetDirectory(nullptr);

//     auto fcorrEE = std::make_shared<TFile>("../data/EnScaleCorr_combined_EE.root", "READ");
//     std::shared_ptr<TH1F> hcorrEE((TH1F*)fcorrEE->Get("hscorr"));
//     hcorrEE->SetDirectory(nullptr);

//     auto nf = df.Define("phoScaleCorrEt_Lead",
//         [hcorrEB, hcorrEE](
//             const float& phoSCEta,
//             const float& phoEt
//         ){

//             bool isEB = abs(phoSCEta) < 1.4442;
//             float SF = 1.;
//             if (isEB){
//                 int bin = hcorrEB->GetXaxis()->FindBin(phoEt);
//                 int Nbin = hcorrEB->GetNbinsX();
//                 if (bin > Nbin || bin < 1)
//                     bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin

//                 SF = 1 + hcorrEB->GetBinContent(bin);
//             }
//             else{
//                 int bin = hcorrEE->GetXaxis()->FindBin(phoEt);
//                 int Nbin = hcorrEE->GetNbinsX();
//                 if (bin > Nbin || bin < 1)
//                     bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin

//                 SF = 1 + hcorrEE->GetBinContent(bin);
//             }

//             return (float) phoEt * SF;
//         }, {"phoSCEta_Lead", "phoConvCorrEt_Lead"});

//     return nf;
// }


// // for Higgs dalitz process
// ROOT::RDF::RNode EnHDALCorrector(ROOT::RDF::RNode df){
//     auto fcorrEB = std::make_shared<TFile>("../data/EnConvCorr_combined_EB_nominal.root", "READ");
//     std::shared_ptr<TH1F> hcorrEB((TH1F*)fcorrEB->Get("hcorr"));
//     hcorrEB->SetDirectory(nullptr);

//     auto fcorrEE = std::make_shared<TFile>("../data/EnConvCorr_combined_EE_nominal.root", "READ");
//     std::shared_ptr<TH1F> hcorrEE((TH1F*)fcorrEE->Get("hcorr"));
//     hcorrEE->SetDirectory(nullptr);

//     auto nf = df.Define("eleHDALCorrEt_lep1",
//         [hcorrEB, hcorrEE](
//             const float& phoSCEta,
//             const float& phoCalibEt
//         ){

//             bool isEB = abs(phoSCEta) < 1.4442;
//             float SF = 1.;
//             if (isEB){
//                 int bin = hcorrEB->GetXaxis()->FindBin(phoCalibEt);
//                 int Nbin = hcorrEB->GetNbinsX();
//                 if (bin > Nbin || bin < 1)
//                     bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin

//                 SF = 1 + hcorrEB->GetBinContent(bin);
//             }
//             else{
//                 int bin = hcorrEE->GetXaxis()->FindBin(phoCalibEt);
//                 int Nbin = hcorrEE->GetNbinsX();
//                 if (bin > Nbin || bin < 1)
//                     bin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin

//                 SF = 1 + hcorrEE->GetBinContent(bin);
//             }

//             return (float) phoCalibEt * SF;
//         }, {"eleSCEta_lep1", "eleCalibEt_lep1"});

//     return nf;
// }