
#include <iostream>
#include "EnCalibrater.h"
#include "ROOT/RVec.hxx"
#include "TH1.h"
#include "TString.h"
#include "TMath.h"


EnCalibrater::~EnCalibrater(){
    if (hasScale){
        delete hscaleEB;
        delete hscaleEE;
    }
    if (hasSmear){
        delete hsmearEB;
        delete hsmearEE;
    }
}


TH1F* EnCalibrater::getHist(std::string histFile, std::string histName){
    TFile* fhist = new TFile(histFile.c_str(), "READ");
    assert(fhist);
    
    TH1F*  hist  = (TH1F*) fhist->Get(histName.c_str());
    assert(hist);

    hist->SetDirectory(nullptr);
    fhist->Close();
    delete fhist;
    
    return hist;
}


void EnCalibrater::SetScaleFiles(std::vector<std::string> scale_files){
    if (scale_files.size() != 2)
        throw std::runtime_error("The size of input files for scale correction should be 2");

    hscaleEB = getHist(scale_files[0], "hscorr");
    hscaleEE = getHist(scale_files[1], "hscorr");
    hasScale = true;
}


void EnCalibrater::SetSmearFiles(std::vector<std::string> smear_files){
    if (smear_files.size() != 2)
        throw std::runtime_error("The size of input files for smear correction should be 2");
    
    hsmearEB = getHist(smear_files[0], "hscorr");
    hsmearEE = getHist(smear_files[1], "hscorr");
    // rnd = new TRandom3(rd());
    // hasSmear = true;
}


ROOT::RDF::RNode EnCalibrater::GetScaleRDF(ROOT::RDF::RNode df, bool isMC, std::string column, std::string column_eta, std::string column_corr){
    auto scale_err = [this](const ROOT::RVec<float>& eleSCEta,
                            const ROOT::RVec<float>& eleHDALRegPt){
                                ROOT::RVec<float> pterr(eleHDALRegPt.size());
                                for (size_t i = 0; i < eleHDALRegPt.size(); i++){
                                    bool isEB = fabs(eleSCEta[i]) < 1.479;
                                    float SFerr = 0.;
                                    if (isEB){
                                        int bin = TMath::Max(1, TMath::Min(hscaleEB->GetNbinsX(), hscaleEB->GetXaxis()->FindBin(eleHDALRegPt[i])));
                                        SFerr = hscaleEB->GetBinError(bin);
                                    }
                                    else{
                                        int bin = TMath::Max(1, TMath::Min(hscaleEE->GetNbinsX(), hscaleEE->GetXaxis()->FindBin(eleHDALRegPt[i])));;
                                        SFerr = hscaleEE->GetBinError(bin);
                                    }
                                    pterr[i] = SFerr;
                                }
                                return pterr;
                            };
            
    if (isMC){
        auto nf = df.Define(column_corr, column)
                    .Define("eleHDALScaleErr", scale_err, {column_eta, column}); // for systematic uncertainties
        return nf;
    }

    auto nf = df.Define(column_corr, [this](const ROOT::RVec<float>& eleSCEta,
                                            const ROOT::RVec<float>& eleHDALRegPt){
                                                ROOT::RVec<float> pt(eleHDALRegPt.size());
                                                for (size_t i = 0; i < pt.size(); i++){
                                                    bool isEB = fabs(eleSCEta[i]) < 1.479;
                                                    float SF = 1.;
                                                    if (isEB){
                                                        int bin = TMath::Max(1, TMath::Min(hscaleEB->GetNbinsX(), hscaleEB->GetXaxis()->FindBin(eleHDALRegPt[i])));
                                                        SF = hscaleEB->GetBinContent(bin);
                                                    }
                                                    else{
                                                        int bin = TMath::Max(1, TMath::Min(hscaleEE->GetNbinsX(), hscaleEE->GetXaxis()->FindBin(eleHDALRegPt[i])));
                                                        SF = hscaleEE->GetBinContent(bin);
                                                    }
                                                    pt[i] = eleHDALRegPt[i] * SF;
                                                }
                                                return pt;
                                            }, {column_eta, column});
    return nf;
}


ROOT::RDF::RNode EnCalibrater::GetSmearRDF(ROOT::RDF::RNode df, bool isMC, std::string column, std::string column_eta, std::string column_corr){
    if (!isMC) // only for MC
        return df.Define(column_corr, column);

    // https://root-forum.cern.ch/t/question-on-usage-of-trandom3-in-lambda-function-with-rdataframe-column-definition/38642/2
    // https://github.com/eguiraud/rdf-benchmarks/blob/main/src/without_io.cpp#L34-L36
    int nslot = TMath::Max(1, (int) df.GetNSlots());
    for (int i = 0; i < nslot; i++){
        rnds.emplace_back(std::make_shared<TRandom3>(rd()));
    }
    
    auto nf = df.DefineSlot("eleHDALSmear2D",[this](unsigned int slot,
                                                    const ROOT::RVec<float>& eleSCEta,
                                                    const ROOT::RVec<float>& eleHDALScalePt){
                                                        ROOT::RVec<float> pt(eleHDALScalePt.size());
                                                        ROOT::RVec<float> ptup(eleHDALScalePt.size());
                                                        ROOT::RVec<float> ptdo(eleHDALScalePt.size());
                                                        for (size_t i = 0; i < pt.size(); i++){
                                                            bool isEB = fabs(eleSCEta[i]) < 1.479;
                                                            float SF    = 0.;
                                                            float SFerr = 0.;
                                                            if (isEB){
                                                                int bin = TMath::Max(1, TMath::Min(hsmearEB->GetNbinsX(), hsmearEB->GetXaxis()->FindBin(eleHDALScalePt[i])));
                                                                SF      = hsmearEB->GetBinContent(bin);
                                                                SFerr   = hsmearEB->GetBinError(bin);
                                                            }
                                                            else{
                                                                int bin = TMath::Max(1, TMath::Min(hsmearEE->GetNbinsX(), hsmearEE->GetXaxis()->FindBin(eleHDALScalePt[i])));
                                                                SF      = hsmearEE->GetBinContent(bin);
                                                                SFerr   = hsmearEE->GetBinError(bin);
                                                            }
                                                            pt[i]   = eleHDALScalePt[i] * rnds[slot]->Gaus(1., SF);
                                                            ptup[i] = eleHDALScalePt[i] * rnds[slot]->Gaus(1., SF+SFerr);
                                                            ptdo[i] = eleHDALScalePt[i] * rnds[slot]->Gaus(1., SF-SFerr);
                                                        }
                                                        ROOT::RVec<ROOT::RVec<float>> pt_all = {pt, ptup, ptdo};
                                                        return pt_all;
                                                    }, {column_eta, column})
                .Define(column_corr,                        "eleHDALSmear2D[0]")
                .Define(Form("%sUp", column_corr.c_str()),  "eleHDALSmear2D[1]")
                .Define(Form("%sDo", column_corr.c_str()),  "eleHDALSmear2D[2]");
    return nf;
}