
#include <iostream>
#include "EnCalibrater.h"
#include "ROOT/RVec.hxx"
#include "TH1.h"
#include "TRandom3.h"
#include "TString.h"


void EnCalibrater::SetCorrFiles(std::vector<std::string> corr_files){
    if (corr_files.size() != 2)
        throw std::runtime_error("The size of input files for correction should be 2");

    _fcorrEB = corr_files[0];
    _fcorrEE = corr_files[1];

    // printf(" [+] Setup converted photon correction files:\n");
    // printf("    - EB: %s\n", _fcorrEB.c_str());
    // printf("    - EE: %s\n", _fcorrEE.c_str());
}


void EnCalibrater::SetScaleFiles(std::vector<std::string> scale_files){
    if (scale_files.size() != 2)
        throw std::runtime_error("The size of input files for scale correction should be 2");

    _fscaleEB = scale_files[0];
    _fscaleEE = scale_files[1];

    // printf(" [+] Setup residual correction files for scale:\n");
    // printf("    - EB: %s\n", _fscaleEB.c_str());
    // printf("    - EE: %s\n", _fscaleEE.c_str());
}


void EnCalibrater::SetSmearFiles(std::vector<std::string> smear_files){
    if (smear_files.size() != 2)
        throw std::runtime_error("The size of input files for smear correction should be 2");
    _fsmearEB = smear_files[0];
    _fsmearEE = smear_files[1];

    // printf(" [+] Setup residual correction files for smear:\n");
    // printf("    - EB: %s\n", _fsmearEB.c_str());
    // printf("    - EE: %s\n", _fsmearEE.c_str());
}


ROOT::RDF::RNode EnCalibrater::GetConvCorr(ROOT::RDF::RNode df, std::string column, std::string column_eta, std::string column_corr){
    std::unique_ptr<TFile> fcorrEB(TFile::Open(_fcorrEB.c_str(), "READ"));
    if (!fcorrEB || fcorrEB->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", _fcorrEB.c_str()));
    std::shared_ptr<TH1F>  hcorrEB((TH1F*)fcorrEB->Get("hcorr"));
    if (!hcorrEB)
        throw std::runtime_error("TFile::Get() failed: no histogram called hcorr for EB");
    hcorrEB->SetDirectory(nullptr);

    std::unique_ptr<TFile> fcorrEE(TFile::Open(_fcorrEE.c_str(), "READ"));
    if (!fcorrEE || fcorrEE->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", _fcorrEE.c_str()));
    std::shared_ptr<TH1F>  hcorrEE((TH1F*)fcorrEE->Get("hcorr"));
    if (!hcorrEE)
        throw std::runtime_error("TFile::Get() failed: no histogram called hcorr for EE");
    hcorrEE->SetDirectory(nullptr);

    // nominal
    auto corr = [hcorrEB, hcorrEE](
        const ROOT::RVec<float>& sceta,
        const ROOT::RVec<float>& pt
    ){
        ROOT::RVec<float> pt_corr(pt.size());
        for (size_t i = 0; i < pt.size(); i++){
            bool isEB = fabs(sceta[i]) < 1.4442;
            float SF = 1.;
            if (isEB){
                const int Nbin = hcorrEB->GetNbinsX();
                int ibin = hcorrEB->GetXaxis()->FindBin(pt[i]);
                if (ibin > Nbin)
                    ibin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (ibin < 1)
                    ibin = 1;
                SF     = hcorrEB->GetBinContent(ibin);
            }
            else{
                const int Nbin = hcorrEE->GetNbinsX();
                int ibin = hcorrEE->GetXaxis()->FindBin(pt[i]);
                if (ibin > Nbin)
                    ibin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (ibin < 1)
                    ibin = 1;
                SF     = hcorrEE->GetBinContent(ibin);
            }

            pt_corr[i]    = pt[i] * SF;
        }

        return pt_corr;
    };

    auto nf = df.Define(column_corr, corr, {column_eta, column});
    return nf;
}


ROOT::RDF::RNode EnCalibrater::GetScaleDif(ROOT::RDF::RNode df, std::string column, std::string column_eta, std::string column_corr, bool isMC){

    std::unique_ptr<TFile> fscaleEB(TFile::Open(_fscaleEB.c_str(), "READ"));
    if (!fscaleEB || fscaleEB->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", _fscaleEB.c_str()));
    std::shared_ptr<TH1F>  hscaleEB((TH1F*)fscaleEB->Get("hscorr"));
    if (!hscaleEB)
        throw std::runtime_error("TFile::Get() failed: no histogram called hscorr for EB");
    hscaleEB->SetDirectory(nullptr);

    std::unique_ptr<TFile> fscaleEE(TFile::Open(_fscaleEE.c_str(), "READ"));
    if (!fscaleEE || fscaleEE->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", _fscaleEE.c_str()));
    std::shared_ptr<TH1F>  hscaleEE((TH1F*)fscaleEE->Get("hscorr"));
    if (!hscaleEE)
        throw std::runtime_error("TFile::Get() failed: no histogram called hscorr for EE");
    hscaleEE->SetDirectory(nullptr);

    auto corr = [hscaleEB, hscaleEE](
        const ROOT::RVec<float>& sceta,
        const ROOT::RVec<float>& pt
    ){
        ROOT::RVec<float> pt_corr(pt.size());
        for (size_t i = 0; i < pt.size(); i++){
            bool isEB = fabs(sceta[i]) < 1.4442;
            float SF = 1.;
            if (isEB){
                const int Nbin = hscaleEB->GetNbinsX();
                int ibin = hscaleEB->GetXaxis()->FindBin(pt[i]);
                if (ibin > Nbin)
                    ibin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (ibin < 1)
                    ibin = 1;
                SF      = hscaleEB->GetBinContent(ibin);
            }
            else{
                const int Nbin = hscaleEE->GetNbinsX();
                int ibin = hscaleEE->GetXaxis()->FindBin(pt[i]);
                if (ibin > Nbin)
                    ibin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (ibin < 1)
                    ibin = 1;
                SF      = hscaleEE->GetBinContent(ibin);
            }
            pt_corr[i]    = pt[i] * SF;
        }
        return pt_corr;
    };

    auto corr_err = [hscaleEB, hscaleEE](
        const ROOT::RVec<float>& sceta,
        const ROOT::RVec<float>& pt
    ){
        ROOT::RVec<float> pt_corr_err(pt.size());
        for (size_t i = 0; i < pt.size(); i++){
            bool isEB = fabs(sceta[i]) < 1.4442;
            float SF_err = 0.;
            if (isEB){
                const int Nbin = hscaleEB->GetNbinsX();
                int ibin = hscaleEB->GetXaxis()->FindBin(pt[i]);
                if (ibin > Nbin)
                    ibin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (ibin < 1)
                    ibin = 1;
                SF_err = hscaleEB->GetBinError(ibin);
            }
            else{
                const int Nbin = hscaleEE->GetNbinsX();
                int ibin = hscaleEE->GetXaxis()->FindBin(pt[i]);
                if (ibin > Nbin)
                    ibin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (ibin < 1)
                    ibin = 1;
                SF_err = hscaleEE->GetBinError(ibin);
            }
            pt_corr_err[i] = SF_err;
        }
        return pt_corr_err;
    };

    if (isMC) {
        auto ff = df.Define(column_corr, column)
                    .Define("eleHDALScale", corr_err, {column_eta, column}); // for systematic variations
        return ff;
    }

    auto nf = df.Define(column_corr,    corr,     {column_eta, column});
    return nf;
}


ROOT::RDF::RNode EnCalibrater::GetSmearDif(ROOT::RDF::RNode df, std::string column, std::string column_eta, std::string column_corr, bool isMC){
    if (!isMC) // only for MC
        return df.Define(column_corr, column);

    std::unique_ptr<TFile> fsmearEB(TFile::Open(_fsmearEB.c_str(), "READ"));
    if (!fsmearEB || fsmearEB->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", _fsmearEB.c_str()));
    std::shared_ptr<TH1F>  hsmearEB((TH1F*)fsmearEB->Get("hscorr"));
    if (!hsmearEB)
        throw std::runtime_error("TFile::Get() failed: no histogram called hscorr for EB");
    hsmearEB->SetDirectory(nullptr);

    std::unique_ptr<TFile> fsmearEE(TFile::Open(_fsmearEE.c_str(), "READ"));
    if (!fsmearEE || fsmearEE->IsZombie())
        throw std::runtime_error(Form("TFile::Open failed: %s", _fsmearEE.c_str()));
    std::shared_ptr<TH1F>  hsmearEE((TH1F*)fsmearEE->Get("hscorr"));
    if (!hsmearEE)
        throw std::runtime_error("TFile::Get() failed: no histogram called hscorr for EE");
    hsmearEE->SetDirectory(nullptr);

    // gRandom->SetSeed(1234);
    auto corr = [hsmearEB, hsmearEE](
        const ROOT::RVec<float>& sceta,
        const ROOT::RVec<float>& pt
    ){
        ROOT::RVec<float> pt_corr(pt.size());
        ROOT::RVec<float> pt_corr_up(pt.size());
        ROOT::RVec<float> pt_corr_do(pt.size());
        for (size_t i = 0; i < pt.size(); i++){
            bool isEB = fabs(sceta[i]) < 1.4442;
            float SF = 1., SF_err = 0.;
            if (isEB){
                const int Nbin = hsmearEB->GetNbinsX();
                int ibin = hsmearEB->GetXaxis()->FindBin(pt[i]);
                if (ibin > Nbin)
                    ibin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (ibin < 1)
                    ibin = 1;
                SF      = hsmearEB->GetBinContent(ibin);
                SF_err  = hsmearEB->GetBinError(ibin);
            }
            else{
                const int Nbin = hsmearEE->GetNbinsX();
                int ibin = hsmearEE->GetXaxis()->FindBin(pt[i]);
                if (ibin > Nbin)
                    ibin = Nbin; // if pt > the last pt bin of hSF, use the sf of the last pt bin
                if (ibin < 1)
                    ibin = 1;
                SF      = hsmearEE->GetBinContent(ibin);
                SF_err  = hsmearEE->GetBinError(ibin);
            }

            pt_corr[i]      = pt[i] * gRandom->Gaus(1., SF);
            pt_corr_up[i]   = pt[i] * gRandom->Gaus(1., (SF + SF_err));
            pt_corr_do[i]   = pt[i] * gRandom->Gaus(1., (SF - SF_err));
        }

        ROOT::RVec<ROOT::RVec<float>> pt_corr_2d = {pt_corr, pt_corr_up, pt_corr_do};
        return pt_corr_2d;
    };

    auto nf = df.Define(Form("%s_2d", column_corr.c_str()),    corr,     {column_eta, column})
                .Define(column_corr,                           Form("%s_2d[0]", column_corr.c_str()));
                // .Define(Form("%s_up", column_corr.c_str()),    Form("%s_2d[1]", column_corr.c_str()))
                // .Define(Form("%s_do", column_corr.c_str()),    Form("%s_2d[2]", column_corr.c_str()));
    return nf;
}