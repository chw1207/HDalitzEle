#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include "boost/program_options.hpp"
#include "boost/lexical_cast.hpp"

#include "TFile.h"
#include "TMath.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooMinuit.h"
#include "RooMinimizer.h"
#include "RooMsgService.h"
#include "RooDataHist.h"
#include "RooExtendPdf.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TMacro.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TArrow.h"
#include "TKey.h"

#include "RooCategory.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"

#include "../interface/PdfModelBuilder.h"
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFunc.h>
#include <iomanip>
#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"

#include "../../tdrStyle/tdrstyle.C"
#include "../../tdrStyle/CMS_lumi_mod.C"

#define RESET   "\033[0m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[93m"      

using namespace std;
using namespace RooFit;
using namespace boost;

namespace po = program_options;

bool BLIND = true;
int nBinsForMass;

string fileName_ = "~/CMSSW_10_2_13/src/flashggFinalFit/Trees2WS/WS/data_obs.root";
string plotDir_ = "plots/fTest";
string outfilename_;
// string multipdfDir_;
string multipdfPath_ = "./Multipdf/Multipdf.root";
string datfile_;
bool runFtestCheckWithToys;
string HDalitzCatsStr_ = "Merged2Gsf_HVBF,Merged2Gsf_LVBF,Merged2Gsf_BST,Merged2Gsf_EBHR9,Merged2Gsf_EBLR9,Merged2Gsf_EE,Merged1Gsf_HVBF,Merged1Gsf_LVBF,Merged1Gsf_BST,Merged1Gsf_EBHR9,Merged1Gsf_EBLR9,Merged1Gsf_EE,Resolved";
// string HDalitzCatsStr_ = "Merged1Gsf_HVBF,Merged1Gsf_LVBF,Merged1Gsf_BST,Merged1Gsf_EBHR9,Merged1Gsf_EBLR9,Merged1Gsf_EE";
vector<string> HDalitzCats_;
string HDalitzCatsNumStr_;
vector<int> HDalitzCatsNum_;
int year_ = 2017;
float lumi_ = 41.5;
float sqrts_ = 13;
float MassLow_ = 105.;
float MassHigh_ = 170.;
string trigname_;
int ncpu_ = 4;
bool verbosity_ = false;
TString extraText_ = "Work-in-progress";
TString lumitext_;
bool saveMultiPdf = true;
bool debug_ = true;
string inWS_ = "w"; // name of the work space
// bool isData_ = true;

TString procText_ = "H #rightarrow #gamma*#gamma #rightarrow ee#gamma";
string massName = "M_{ee#gamma} [GeV]"; // use to plot the x title //!changed
int ntoys_ = 500; //modified of Dalitz; formerly 500

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

TRandom3 *RandomGen = new TRandom3();
float blind_low = 120.;
float blind_high = 130.;

void OptionParser(int argc, char *argv[]){

    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,         h",         "Show help")
    ("infilename,   i",         po::value<string>(&fileName_)->default_value(fileName_.c_str()), "In file name")
    ("plotDir,     p",         po::value<string>(&plotDir_)->default_value(plotDir_.c_str()), "Put plots in this directory")
    // ("multipdfDir",             po::value<string>(&multipdfDir_)->default_value("Multipdf"), "Directory for MultiPdf model")
    ("multipdfPath",            po::value<string>(&multipdfPath_)->default_value(multipdfPath_.c_str()), "Path for MultiPdf model")
    ("inWS",                    po::value<string>(&inWS_)->default_value(inWS_.c_str()), "Name of the input WS")
    ("datfile,      d",         po::value<string>(&datfile_)->default_value("dat/fTest.dat"), "Write results to datfile for BiasStudy")
    ("runFtestCheckWithToys",   po::value<bool>(&runFtestCheckWithToys)->default_value(false), "When running the F-test, use toys to calculate p-value (and make plots)")
    ("HDalitzCats",             po::value<string>(&HDalitzCatsStr_)->default_value(HDalitzCatsStr_.c_str()), "Higgs Dalitz category")
    ("HDalitzCatsNum",          po::value<string>(&HDalitzCatsNumStr_)->default_value("1,1,1,1,1,1,1"), "Higgs Dalitz category number") //!FIXEDME
    ("year",                    po::value<int>(&year_)->default_value(year_), "Dataset year")
    ("luminosity",              po::value<float>(&lumi_)->default_value(lumi_), "Luminosity (Default: 122.7/fb)")
    ("sqrts",                   po::value<float>(&sqrts_)->default_value(sqrts_), "Center of energy [GeV]")
    ("MassLow",                 po::value<float>(&MassLow_)->default_value(MassLow_), "Lower bound of Mass")
    ("MassHigh",                po::value<float>(&MassHigh_)->default_value(MassHigh_), "Upper bound of Mass")
    ("proctext",                po::value<TString>(&procText_)->default_value(procText_), "Process")
    // ("nThreads,t",              po::value<int>(&ncpu_)->default_value(ncpu_), "Number of threads to be used for the fits (Default = 4)")
    ("verbosity",               po::value<bool>(&verbosity_)->default_value(verbosity_), "Roofit verbosity (Default = false)")
    ("blind",                   po::value<bool>(&BLIND)->default_value(BLIND), "Blind analysis");
    // ("isData",                  po::value<bool>(&isData_)->default_value(isData_), "is data or not");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")){
        cout << desc << endl;
        exit(1);
    }

    if (!verbosity_){
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(true);
        gErrorIgnoreLevel = kWarning;
    }

    split(HDalitzCats_, HDalitzCatsStr_, boost::is_any_of(","));
    vector<string> v_catnum;
    split(v_catnum, HDalitzCatsNumStr_, boost::is_any_of(","));
    for (vector<string>::iterator it = v_catnum.begin(); it != v_catnum.end(); it++){
        HDalitzCatsNum_.push_back(boost::lexical_cast<int>(*it));
    }

    nBinsForMass = (int)(MassHigh_ - MassLow_); //!FIXEDME
}

RooAbsPdf* getPdf(PdfModelBuilder &pdfsModel, string type, int order, const char* ext = " "){
  
    if (type == "Bernstein") 
        return pdfsModel.getBernstein(Form("%s_bern%d", ext, order), order); 
    else if (type == "Chebychev") 
        return pdfsModel.getChebychev(Form("%s_cheb%d", ext, order), order); 
    else if (type == "Exponential") 
        return pdfsModel.getExponentialSingle(Form("%s_exp%d", ext, order), order); 
    else if (type == "PowerLaw") 
        return pdfsModel.getPowerLawSingle(Form("%s_pow%d", ext, order), order); 
    else if (type == "Laurent") 
        return pdfsModel.getLaurentSeries(Form("%s_lau%d", ext, order), order); 
    else {
        cerr << "[ERROR] -- getPdf() -- type " << type << " not recognised." << endl;
        return NULL;
    }
} 


void runFit(RooAbsPdf *pdf, RooDataSet *data, double *NLL, int *stat_t, int MaxTries){
    int ntries = 0;
    RooArgSet *params_test = pdf->getParameters((const RooArgSet*)(0));
    //params_test->Print("v");
    int stat = 1;
    double minnll = 10e8;
    while ((stat != 0) && (ntries <= MaxTries)){
        // if (ntries >= MaxTries) break;
        
        RooFitResult *fitTest = pdf->fitTo(*data, RooFit::Save(1), RooFit::Minimizer("Minuit2", "minimize"), RooFit::SumW2Error(true)); //FIXME
        stat = fitTest->status();
        minnll = fitTest->minNll();

        if (stat != 0) params_test->assignValueOnly(fitTest->randomizePars());
        ntries++; 
    }

    *stat_t = stat;
    *NLL = minnll;
}


void transferMacros(TFile *inFile, TFile *outFile){
  
    TIter next(inFile->GetListOfKeys());
    TKey *key;
    while ((key = (TKey*)next())){
        if (string(key->ReadObj()->ClassName()) == "TMacro") {
            //cout << key->ReadObj()->ClassName() << " : " << key->GetName() << endl;
            TMacro *macro = (TMacro*)inFile->Get(key->GetName());
            outFile->cd();
            macro->Write();
        }
    }
}


double getProbabilityFtest(double chi2, int ndof, RooAbsPdf *pdfNull, RooAbsPdf *pdfTest, RooRealVar *mass, RooDataSet *data, string name){

    double prob_asym = TMath::Prob(chi2, ndof);
    if (!runFtestCheckWithToys)
        return prob_asym;

    int ndata = data->sumEntries();

    // fit the pdfs to the data and keep this fit Result (for randomizing)
    RooFitResult *fitNullData = pdfNull->fitTo(*data, RooFit::Save(1), RooFit::Strategy(1) ,RooFit::Minimizer("Minuit2","minimize"), RooFit::SumW2Error(true), RooFit::PrintLevel(-1)); //FIXME
    RooFitResult *fitTestData = pdfTest->fitTo(*data, RooFit::Save(1), RooFit::Strategy(1), RooFit::Minimizer("Minuit2", "minimize"), RooFit::SumW2Error(true), RooFit::PrintLevel(-1)); //FIXME

    // Ok we want to check the distribution in toys then
    // Step 1, cache the parameters of each pdf so as not to upset anything
    RooArgSet *params_null = pdfNull->getParameters((const RooArgSet *)(0));
    RooArgSet preParams_null;
    params_null->snapshot(preParams_null);
    RooArgSet *params_test = pdfTest->getParameters((const RooArgSet *)(0));
    RooArgSet preParams_test;
    params_test->snapshot(preParams_test);

    int ntoys = ntoys_;
    TCanvas *can = new TCanvas();
    can->SetLogy();
    TH1F toyhist(Form("toys_fTest_%s.pdf", pdfNull->GetName()), ";Chi2;", 60, -2, 10);
    TH1I toyhistStatN(Form("Status_%s.pdf", pdfNull->GetName()), ";FitStatus;", 8, -4, 4);
    TH1I toyhistStatT(Form("Status_%s.pdf", pdfTest->GetName()), ";FitStatus;", 8, -4, 4);

    TGraph *gChi2 = new TGraph();
    gChi2->SetLineColor(kGreen + 2);
    double w = toyhist.GetBinWidth(1);

    int ipoint = 0;

    for (int b = 0; b < toyhist.GetNbinsX(); b++){
        double x = toyhist.GetBinCenter(b + 1);
        if (x > 0){
            gChi2->SetPoint(ipoint, x, (ROOT::Math::chisquared_pdf(x, ndof)));
            ipoint++;
        }
    }
    int npass = 0;
    int nsuccesst = 0;
    mass->setBins(nBinsForMass);
    for (int itoy = 0; itoy < ntoys; itoy++){

        params_null->assignValueOnly(preParams_null);
        params_test->assignValueOnly(preParams_test);
        RooDataHist *binnedtoy = pdfNull->generateBinned(RooArgSet(*mass), ndata, 0, 1);

        int stat_n = 1;
        int stat_t = 1;
        int ntries = 0;
        double nllNull, nllTest;
        // Iterate on the fit
        int MaxTries = 2;
        while (stat_n != 0){
            if (ntries >= MaxTries)
                break;
            RooFitResult *fitNull = pdfNull->fitTo(
                *binnedtoy, RooFit::Save(1), RooFit::Strategy(1), RooFit::SumW2Error(kTRUE), //FIXME
                RooFit::Minimizer("Minuit2", "minimize"), RooFit::Minos(0), RooFit::Hesse(0), RooFit::PrintLevel(-1)
            );
            //,RooFit::Optimize(0));

            nllNull = fitNull->minNll();
            stat_n = fitNull->status();
            if (stat_n != 0)
                params_null->assignValueOnly(fitNullData->randomizePars());
            ntries++;
        }

        ntries = 0;
        while (stat_t != 0){
            if (ntries >= MaxTries)
                break;
            RooFitResult *fitTest = pdfTest->fitTo(
                *binnedtoy, RooFit::Save(1), RooFit::Strategy(1), RooFit::SumW2Error(true), //FIXME
                RooFit::Minimizer("Minuit2", "minimize"), RooFit::Minos(0), RooFit::Hesse(0), RooFit::PrintLevel(-1)
            );
            nllTest = fitTest->minNll();
            stat_t = fitTest->status();
            if (stat_t != 0)
                params_test->assignValueOnly(fitTestData->randomizePars());
            ntries++;
        }

        toyhistStatN.Fill(stat_n);
        toyhistStatT.Fill(stat_t);

        if (stat_t != 0 || stat_n != 0)
            continue;
        nsuccesst++;
        double chi2_t = 2 * (nllNull - nllTest);
        if (chi2_t >= chi2)
            npass++;
        toyhist.Fill(chi2_t);
    }

    double prob = 0;
    if (nsuccesst != 0)
        prob = (double)npass / nsuccesst;
    toyhist.Scale(1. / (w * toyhist.Integral()));
    toyhist.Draw();
    TArrow lData(chi2, toyhist.GetMaximum(), chi2, 0);
    lData.SetLineWidth(2);
    lData.Draw();
    gChi2->Draw("L");
    TLatex *lat = new TLatex();
    lat->SetNDC();
    lat->SetTextFont(42);
    lat->DrawLatex(0.1, 0.91, Form("Prob (asymptotic) = %.4f (%.4f)", prob, prob_asym));
    can->Print(name.c_str());

    TCanvas *stas = new TCanvas();
    toyhistStatN.SetLineColor(2);
    toyhistStatT.SetLineColor(1);
    TLegend *leg = new TLegend(0.2, 0.6, 0.4, 0.87);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->AddEntry(&toyhistStatN, "Null Hyp", "L");
    leg->AddEntry(&toyhistStatT, "Test Hyp", "L");
    toyhistStatN.Draw();
    toyhistStatT.Draw("same");
    leg->Draw();
    stas->Print(Form("%s_fitstatus.pdf", name.c_str()));
    //reassign params
    params_null->assignValueOnly(preParams_null);
    params_test->assignValueOnly(preParams_test);

    delete can;
    delete stas;
    delete gChi2;
    delete leg;
    delete lat;

    // Still return the asymptotic prob (usually its close to the toys one)
    return prob_asym;
}

double getGoodnessOfFit(RooRealVar *mass, RooAbsPdf *mpdf, RooDataSet *data, string name){
    
    mass->setRange("unblindReg_1", MassLow_, blind_low);
    mass->setRange("unblindReg_2", blind_high, MassHigh_);

    double prob;
    int ntoys = ntoys_;
    // Routine to calculate the goodness of fit. 
    
    name += "_gofTest.pdf";
    RooRealVar norm("norm", "norm", data->sumEntries(), 0, 10E6);
    //norm.removeRange();

    RooExtendPdf *pdf = new RooExtendPdf("ext", "ext", *mpdf, norm);

    // (Originally the Chi2 value is taken from the data in full range)
    // The question is whether this is correct
    //!FIXEDME: How to properly calculate Chi2 value with RooExtendPdf in sideband fit? 
    RooPlot *plot_chi2 = mass->frame();
    data->plotOn(plot_chi2, Binning(nBinsForMass), Name("data"));

    pdf->plotOn(plot_chi2, Name("pdf"));
    int np = pdf->getParameters(*data)->getSize();

    double chi2 = plot_chi2->chiSquare("pdf", "data", np);
    cout << "[INFO] Calculating GOF for pdf " << pdf->GetName() << ", using " << np << " fitted parameters" << endl;

    // The first thing is to check if the number of entries in any bin is < 5 
    // if so, we don't rely on asymptotic approximations
    if ((double)data->sumEntries()/nBinsForMass < 5 ){

        cout << "[INFO] Running toys for GOF test " << endl;
        // store pre-fit params 
        RooArgSet *params = pdf->getParameters(*data);
        RooArgSet preParams;
        params->snapshot(preParams);
        int ndata = data->sumEntries();
    
        int npass = 0;
        vector<double> toy_chi2;
        for (int itoy = 0 ; itoy < ntoys ; itoy++){
            //  cout << "[INFO] " <<Form("\t.. %.1f %% complete\r",100*float(itoy)/ntoys) << flush;
            params->assignValueOnly(preParams);
            int nToyEvents = RandomGen->Poisson(ndata);
            RooDataHist *binnedtoy = pdf->generateBinned(RooArgSet(*mass), nToyEvents, 0,1 );
            pdf->fitTo(*binnedtoy, RooFit::Minimizer("Minuit2", "minimize"), RooFit::Minos(0), RooFit::Hesse(0), RooFit::PrintLevel(-1), RooFit::Strategy(0), RooFit::SumW2Error(kTRUE)); //FIXME

            RooPlot *plot_t = mass->frame();
            binnedtoy->plotOn(plot_t);
            pdf->plotOn(plot_t); //,RooFit::NormRange("fitdata_1,fitdata_2"));

            double chi2_t = plot_t->chiSquare(np);
            if( chi2_t >= chi2) npass++;
            toy_chi2.push_back(chi2_t*(nBinsForMass - np));
            delete plot_t;
        }
        cout << "[INFO] complete" << endl;
        prob = (double)npass / ntoys;

        TCanvas *can = new TCanvas();
        double medianChi2 = toy_chi2[(int)(((float)ntoys)/2)];
        double rms = TMath::Sqrt(medianChi2);

        TH1F toyhist(Form("gofTest_%s.pdf", pdf->GetName()), ";Chi2;" , 50, medianChi2 - 5*rms, medianChi2 + 5*rms);
        for (vector<double>::iterator itx = toy_chi2.begin(); itx!=toy_chi2.end(); itx++){
            toyhist.Fill((*itx));
        }
        toyhist.Draw();

        TArrow lData(chi2*(nBinsForMass-np), toyhist.GetMaximum(), chi2*(nBinsForMass - np), 0);
        lData.SetLineWidth(2);
        lData.Draw();
        can->SaveAs(name.c_str());

        // back to best fit 	
        params->assignValueOnly(preParams);
    } 
    else {
        prob = TMath::Prob(chi2*(nBinsForMass - np), nBinsForMass - np);
    }
    
    cout << "[INFO] Chi2 in Observed =  " << chi2*(nBinsForMass - np) << endl;
    cout << "[INFO] p-value  =  " << prob << endl;
    delete pdf;
    return prob;
}


// Plot single fit
void plot_1(RooRealVar *mass, RooAbsPdf *pdf, RooDataSet *data, string name, int status, double *prob){
    
    // (Originally the Chi2 is taken from full range fit)
    RooPlot *plot_chi2 = mass->frame();
    data->plotOn(plot_chi2, Binning(nBinsForMass));
    pdf->plotOn(plot_chi2);
    
    int np = pdf->getParameters(*data)->getSize() + 1; //Because this pdf has no extend
    cout << "np = " << np << endl; 
    double chi2 = plot_chi2->chiSquare(np);

    *prob = getGoodnessOfFit(mass, pdf, data, name);
    RooPlot *plot = mass->frame();
    plot->GetXaxis()->SetTitle(massName.c_str());
    // plot->GetYaxis()->SetTitle(Form("Events / %.1f GeV", (float) (MassHigh_ - MassLow_)/nBinsForMass));

    mass->setRange("unblindReg_1", MassLow_, blind_low);
    mass->setRange("unblindReg_2", blind_high, MassHigh_);
    if (BLIND){
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_1"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_2"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), Invisible());
    }
    else
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), XErrorSize(0.0001));

    // Draw the single fit results
    

    plot->SetMaximum(plot->GetMaximum() * 1.3);
    if (BLIND)
        plot->SetMinimum(0.0001);
    plot->SetTitle("");

    TCanvas *canv = new TCanvas("canv", "", 700, 600);
    pdf->plotOn(plot);//,RooFit::NormRange("fitdata_1,fitdata_2"));
    bool paramoncanv = false;
    if (paramoncanv)
        pdf->paramOn(plot, RooFit::Layout(0.17, 0.93, 0.89), RooFit::Format("NEA", AutoPrecision(1)));
    if (BLIND) plot->SetMinimum(0.0001);
    plot->SetTitle("");
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.15);
    canv->cd();
    plot->Draw();

    CMS_lumi(canv, 5, 11, lumitext_, year_, true, extraText_, procText_, " ");
    canv->Update();
    canv->RedrawAxis();
    TLatex *lat = new TLatex();
    lat->SetNDC();
    lat->SetTextFont(36);
    lat->DrawLatex(0.60, 0.83, Form("#chi^{2} = %g", chi2 * (nBinsForMass - np)));
    lat->DrawLatex(0.60, 0.76, Form("Prob. = %g", *prob));
    lat->DrawLatex(0.60, 0.70, Form("Fit Status = %d", status));

    canv->SaveAs(Form("%s.pdf", name.c_str()));

    delete canv;
    delete lat;
}

void plot_2(RooRealVar *mass, map<string, RooAbsPdf *> pdfs, RooDataSet *data, string name, vector<string> HDalitzCats_, int cat, int bestFitPdf = -1){
    
    TCanvas *canv = new TCanvas("canv", "", 700, 600);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.15);
    canv->cd();
    TLegend *leg = new TLegend(0.5, 0.6, 0.88, 0.88);
    leg->SetFillColor(0);
    leg->SetLineColor(0);
    leg->SetTextSize(0.032);
    RooPlot *plot = mass->frame();
    plot->GetXaxis()->SetTitle(massName.c_str());
    plot->GetXaxis()->SetTitleOffset(1.4);
    // plot->GetYaxis()->SetTitle(Form("Events / %.1f GeV", (float) (MassHigh_ - MassLow_)/nBinsForMass));
    plot->GetYaxis()->SetTitleOffset(1.6);
    mass->setRange("unblindReg_1", MassLow_, blind_low);
    mass->setRange("unblindReg_2", blind_high, MassHigh_);
    if (BLIND){
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_1"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_2"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), Invisible());
    }
    else
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), XErrorSize(0.0001));

    TObject *datLeg = plot->getObject(int(plot->numItems() - 1));
    if (HDalitzCats_.size() > 0){
        leg->AddEntry(datLeg, Form("Data - %s", HDalitzCats_[cat].c_str()), "EP");
    }
    else{
        leg->AddEntry(datLeg, Form("Data - %d", cat), "EP");
    }

    int i = 0;
    int style = 1;
    for (map<string, RooAbsPdf *>::iterator it = pdfs.begin(); it != pdfs.end(); it++){
        it->second->plotOn(plot, LineColor(TColor::GetColorPalette((i) * (int)(256/(int)pdfs.size() - 1))), LineStyle(style)); //,RooFit::NormRange("fitdata_1,fitdata_2"));
        TObject *pdfLeg = plot->getObject(int(plot->numItems() - 1));
        std::string ext = "";
        if (bestFitPdf == i)
            ext = " (Best Fit Pdf) ";
        leg->AddEntry(pdfLeg, Form("%s%s", it->first.c_str(), ext.c_str()), "L");
        i++;
    }
    // plot dataset again to make the data points on top of the fit functions
    if (BLIND){
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_1"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_2"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), Invisible());
    }
    else
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), XErrorSize(0.0001));

    plot->SetTitle(Form(" %s", HDalitzCats_[cat].c_str()));
    if (BLIND)
        plot->SetMinimum(0.0001);
    plot->SetMaximum(plot->GetMaximum() * 1.5);
    plot->Draw();
    leg->Draw("same");
    CMS_lumi(canv, 5, 11, lumitext_, year_, true, extraText_, procText_, "");
    
    canv->SaveAs(Form("%s.pdf", name.c_str()));
    // canv->SaveAs(Form("%s.png", name.c_str()));
    delete canv;
}

int getBestFitFunction(RooMultiPdf *bkg, RooDataSet *data, RooCategory *cat, bool silent = false){
    
    double global_minNll = 1E10;
	int best_index = 0;
	int number_of_indeces = cat->numTypes();
		
	RooArgSet snap,clean;
	RooArgSet *params = bkg->getParameters((const RooArgSet*)0);
	params->remove(*cat);
	params->snapshot(snap);
	params->snapshot(clean);
	// if (!silent) {
	// 	params->Print("V");
	// }
	
	//bkg->setDirtyInhibit(1);
	//RooAbsReal *nllm = bkg->createNLL(*data);
	//RooMinimizer minim(*nllm);
	//minim.setStrategy(1);
	
	for (int id = 0; id < number_of_indeces; id++){		
		params->assignValueOnly(clean);
		cat->setIndex(id);

		//RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);

		// if (!silent) {
		// 	std::cout << "BEFORE  MAKING FIT" << std::endl;
		// 	params->Print("V");
		// 	std::cout << "-----------------------" << std::endl;		
		// }
		
		//minim.minimize("Minuit2","minimize");
		double minNll = 0; //(nllm->getVal())+bkg->getCorrection();
		int fitStatus = 1;		
		runFit(bkg->getCurrentPdf(), data, &minNll, &fitStatus, /*max iterations*/3);
		// Add the penalty

		minNll = minNll + bkg->getCorrection();

		if (!silent) {
			/*
			std::cout << "After Minimization ------------------  " <<std::endl;
			std::cout << bkg->getCurrentPdf()->GetName() << " " << minNll <<std::endl;
			bkg->Print("v");
			bkg->getCurrentPdf()->getParameters(*data)->Print("V");
			std::cout << " ------------------------------------  " << std::endl;
	
			params->Print("V");
			*/
			cout << "[INFO] AFTER FITTING" << endl;
			cout << "[INFO] Function was " << bkg->getCurrentPdf()->GetName() <<endl;
			cout << "[INFO] Correction Applied is " << bkg->getCorrection() <<endl;
			cout << "[INFO] NLL + c = " <<  minNll << endl;
			cout << "-----------------------" << endl;
		}
			
		if (minNll < global_minNll){
        	global_minNll = minNll;
			snap.assignValueOnly(*params);
        	best_index = id;
		}
	}
    
    cat->setIndex(best_index);
	params->assignValueOnly(snap);
	
	if (!silent) {
		cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<endl;
		//bkg->getCurrentPdf()->getParameters(*data)->Print("v");
	}
	return best_index;
}

vector<string> split_string(string instr, string delimiter = "_"){
    // Ref: https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = instr.find(delimiter, prev);
        if (pos == string::npos) pos = instr.length();
        string token = instr.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delimiter.length();
    }
    while (pos < instr.length() && prev < instr.length());
    return tokens;
}

string ReplaceAll(string str, const string& from, const string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}


void plot_3(RooRealVar *mass, RooMultiPdf *pdfs, RooCategory *catIndex, RooDataSet *data, string name, vector<string> HDalitzCats_, int cat, int bestFitPdf = -1){

    // int color[7] = {kBlue, kRed, kMagenta, kGreen + 1, kOrange + 7, kAzure + 10, kBlack};

    TLegend *leg = new TLegend(0.42, 0.57, 0.92, 0.86);
    leg->SetFillColor(0);
    leg->SetLineColor(1);
    leg->SetTextSize(0.047);
    leg->SetTextFont(42);

    RooPlot *plot = mass->frame();
    mass->setRange("unblindReg_1", MassLow_, blind_low);
    mass->setRange("unblindReg_2", blind_high, MassHigh_);
    if (BLIND)
    {
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_1"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_2"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), Invisible());
    }
    else
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), XErrorSize(0.0001));
    
    TCanvas *canv = new TCanvas("canv", "", 800, 800);
    // gPad->SetLeftMargin(0.15);
    // gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.1);
    canv->cd();

    ///start extra bit for ratio plot///
    RooHist *plotdata = (RooHist *)plot->getObject(plot->numItems() - 1);
    // bool doRatioPlot_ = 1;
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad1->SetBottomMargin(0.05);
    // pad1->SetTopMargin(0.01);
    pad1->SetRightMargin(0.05);
    pad1->SetLeftMargin(0.13);
    pad1->SetBottomMargin(0.03);

    pad2->SetRightMargin(0.05);
    pad2->SetLeftMargin(0.13);
    pad2->SetTopMargin(0.06);
    pad2->SetBottomMargin(0.4);

    // pad1->SetLeftMargin(0.15);
    // pad1->SetRightMargin(0.05);
    // pad1->SetBottomMargin(0.0);
    // pad1->SetTopMargin(0.11);
    // pad2->SetLeftMargin(0.15);
    // pad2->SetRightMargin(0.05);
    // pad2->SetTopMargin(0.0);
    // pad2->SetBottomMargin(0.35);
    pad1->Draw();
    pad2->Draw();
    pad1->cd();
    // enf extra bit for ratio plot///

    int currentIndex = catIndex->getIndex();
    TObject *datLeg = plot->getObject(int(plot->numItems() - 1));
    leg->AddEntry(datLeg, Form("Data - %s", HDalitzCats_[cat].c_str()), "EP");
    int style = 1;
    RooAbsPdf *pdf;
    RooCurve *nomBkgCurve;
    int bestcol = -1;
    for (int icat = 0; icat < catIndex->numTypes(); icat++)
    {
        catIndex->setIndex(icat);
        pdfs->getCurrentPdf()->fitTo(*data, RooFit::Minos(0), RooFit::Minimizer("Minuit2", "minimize"), RooFit::SumW2Error(kTRUE)); //FIXME
        pdfs->getCurrentPdf()->plotOn(plot, LineColor(TColor::GetColorPalette((icat) * (int)(256/catIndex->numTypes() - 1))), LineStyle(style));                                                      
        TObject *pdfLeg = plot->getObject(int(plot->numItems() - 1));
        std::string ext = "";
        if (bestFitPdf == icat)
        {
            ext = " (Best Fit Pdf) ";
            pdf = pdfs->getCurrentPdf();
            nomBkgCurve = (RooCurve *)plot->getObject(plot->numItems() - 1);
            bestcol = icat;
        }

        vector<string> Pdfname = split_string(pdfs->getCurrentPdf()->GetName());
        leg->AddEntry(pdfLeg, Form("%s%s", Pdfname[Pdfname.size()-1].c_str(), ext.c_str()), "L");
    }

    if (BLIND)
    {
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_1"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), CutRange("unblindReg_2"), XErrorSize(0.0001));
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), Invisible());
    }
    else
        data->plotOn(plot, Binning(MassHigh_ - MassLow_), XErrorSize(0.0001));

    plot->SetTitle(Form("Category %s", HDalitzCats_[cat].c_str()));
    // plot->GetYaxis()->SetTitle(Form("Events / %.1f GeV", (float) (MassHigh_ - MassLow_)/nBinsForMass));

    plot->GetXaxis()->SetLabelOffset(0.05);
    plot->GetXaxis()->SetTickSize(0.03);
    plot->GetYaxis()->SetTitleSize(0.055);
    plot->GetYaxis()->SetTickSize(0.03);
    plot->GetYaxis()->SetLabelSize(0.05);
    plot->GetYaxis()->SetTitleOffset(1.2);

    // plot->GetYaxis()->SetTitleOffset(1.3);
    // plot->GetYaxis()->SetTitleSize(0.05);
    // plot->GetYaxis()->SetLabelSize(0.05);
    if (BLIND)
        plot->SetMinimum(0.0001);
    plot->SetMaximum(plot->GetMaximum() * 1.7);
    plot->Draw();
    leg->Draw("same");

    CMS_lumi(pad1, 5, 11, lumitext_, year_, true, extraText_, procText_, "");
    pad1->Update();
    canv->cd();


    ///start extra bit for ratio plot///
    TH1D *hbplottmp = (TH1D *)pdf->createHistogram("hbplottmp", *mass, Binning(MassHigh_ - MassLow_, MassLow_, MassHigh_));
    hbplottmp->Scale(plotdata->Integral());
    hbplottmp->Draw("same");
    int npoints = plotdata->GetN();
    double xtmp, ytmp; //
    int point = 0;
    TGraphAsymmErrors *hdatasub = new TGraphAsymmErrors(npoints);
    //hdatasub->SetMarkerSize(defmarkersize);
    for (int ipoint = 0; ipoint < npoints; ++ipoint)
    {
        //double bkgval = hbplottmp->GetBinContent(ipoint+1);
        plotdata->GetPoint(ipoint, xtmp, ytmp);
        double bkgval = nomBkgCurve->interpolate(xtmp);
        if (BLIND)
        {
            if ((xtmp > blind_low) && (xtmp < blind_high))
                continue;
        }
        std::cout << "[INFO] plotdata->Integral() " << plotdata->Integral() << " ( bins " << npoints << ") hbkgplots[i]->Integral() " << hbplottmp->Integral() << " (bins " << hbplottmp->GetNbinsX() << std::endl;
        double errhi = plotdata->GetErrorYhigh(ipoint);
        double errlow = plotdata->GetErrorYlow(ipoint);

        //std::cout << "[INFO]  Channel " << name  << " errhi " << errhi << " errlow " << errlow  << std::endl;
        std::cout << "[INFO] Channel  " << name << " setting point " << point << " : xtmp " << xtmp << "  ytmp " << ytmp << " bkgval  " << bkgval << " ytmp-bkgval " << ytmp - bkgval << std::endl;
        bool drawZeroBins_ = 1;
        if (!drawZeroBins_)
            if (fabs(ytmp) < 1e-5)
                continue;
        hdatasub->SetPoint(point, xtmp, ytmp - bkgval);
        hdatasub->SetPointError(point, 0., 0., errlow, errhi);
        point++;
    }
    pad2->cd();
    TH1 *hdummy = new TH1D("hdummyweight", "", MassHigh_ - MassLow_, MassLow_, MassHigh_);
    hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum() + 1);
    hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum() - 1);
    
    hdummy->GetYaxis()->SetTitle("Data - best fit");
    hdummy->GetYaxis()->SetTitleOffset(0.5);
    hdummy->GetYaxis()->SetTitleSize(0.12);
    hdummy->GetYaxis()->SetLabelSize(0.115);
    hdummy->GetYaxis()->SetNdivisions(505);

    hdummy->GetXaxis()->SetTitle(massName.c_str());
    hdummy->GetXaxis()->SetTitleSize(0.15);
    hdummy->GetXaxis()->SetLabelSize(0.115);
    hdummy->GetXaxis()->SetTitleOffset(1.2);
    hdummy->GetXaxis()->SetTickSize(0.07);
    hdummy->GetXaxis()->SetLabelOffset(0.045);
    
    hdummy->Draw("HIST");

    TLine *line3 = new TLine(MassLow_, 0., MassHigh_, 0.);
    line3->SetLineColor(TColor::GetColorPalette((bestcol) * (int)(256/catIndex->numTypes() - 1)));
    //line3->SetLineStyle(kDashed);
    line3->SetLineWidth(5.0);
    line3->Draw();
    hdatasub->Draw("PESAME");
    // enf extra bit for ratio plot///
    canv->SaveAs(Form("%s.pdf", name.c_str()));
    // canv->SaveAs(Form("%s.png", name.c_str()));
    catIndex->setIndex(currentIndex);
    delete canv;
}


string getDirName(string s, string del = "/"){
    // split the string into vector<string>
    // eg. s = "./test/testc/testa.root"(full file name) -> svec = {., test, testc, testa.root}
    vector<string> svec;
    svec.clear();
    split(svec, s, boost::is_any_of("/"));

    // get the directory (eg. "./test/testc/")
    string sDir = "";
    for (size_t i = 0; i < (svec.size()-1); i++){
        sDir.append(svec[i].c_str());
        sDir.append("/");
    }
    return sDir;
}

//====================================================================//
//                                                                    // 
//                      Main fitting part!                            //
//                                                                    //
//====================================================================//
int main(int argc, char *argv[]){
    OptionParser(argc, argv);

    setTDRStyle();
    lumitext_.Form("%.1f fb^{-1}", lumi_);

    // Setup the output file and ws 
    printf("[INFO] SaveMultiPdf? %s\n", saveMultiPdf ? "true" : "false");
    TFile *outputfile;
    RooWorkspace *outputws;
    string multipdfDir_ = getDirName(multipdfPath_, "/");
    if (saveMultiPdf){
        system(Form("mkdir -p %s", multipdfDir_.c_str()));
        outputfile = new TFile(multipdfPath_.c_str(), "RECREATE");
        outputws = new RooWorkspace();
        outputws->SetName("multipdf");
    }

    // Open the infile 
    TFile *inFile = TFile::Open(fileName_.c_str(), "READ");
    RooWorkspace *inWS = (RooWorkspace *) inFile->Get(inWS_.c_str());
    cout << "[INFO] inWS open " << inWS->GetName() << endl;

    if (saveMultiPdf){
        transferMacros(inFile, outputfile);

        RooRealVar *intLumi_ = new RooRealVar("IntLumi", "", lumi_);
        RooRealVar *sqrts = new RooRealVar("SqrtS", "SqrtS", sqrts_);

        outputws->import(*intLumi_);
        outputws->import(*sqrts);

        if (verbosity_)
            cout << "[INFO] Got intLumi and sqrts " << intLumi_ << ", " << sqrts << endl;
    }

    vector<string> functionClasses;
	functionClasses.push_back("Bernstein");
	functionClasses.push_back("Exponential");
	functionClasses.push_back("PowerLaw");
	// functionClasses.push_back("Laurent");
	map<string,string> namingMap;
	namingMap.insert(pair<string, string>("Bernstein", "pol"));
	namingMap.insert(pair<string, string>("Exponential", "exp"));
	namingMap.insert(pair<string, string>("PowerLaw", "pow"));
	// namingMap.insert(pair<string, string>("Laurent", "lau"));

    // store results here
    system(Form("mkdir -p %s/single_fit", plotDir_.c_str()));
    FILE *resFile;
    resFile = fopen(Form("%s/fTestResults.txt", plotDir_.c_str()), "w");
    cout << Form("[INFO] FitResults file: %s/fTestResults.txt", plotDir_.c_str()) << endl;
    vector<map<string, int>> choices_vec;
    vector<map<string, vector<int>>> choices_envelope_vec;
    vector<map<string, RooAbsPdf *>> pdfs_vec;

    PdfModelBuilder pdfsModel;
    RooRealVar *mass = (RooRealVar *)inWS->var("CMS_higgs_mass");
    cout << "[INFO] Got mass from ws: " << mass->GetName() << endl;

    pdfsModel.setObsVar(mass);
    double upperEnvThreshold = 0.1; // upper threshold on delta(chi2) to include function in envelope (looser than truth function)

    fprintf(resFile, "Truth Model & d.o.f & $\\Delta NLL_{N+1}$ & $p(\\chi^{2}>\\chi^{2}_{(N\\rightarrow N+1)})$ \\\\\n");
    fprintf(resFile, "\\hline\n");

    // string ext = Form("%s_13TeV", trigname_.c_str());
    string ext = "DiPho_13TeV"; //! FIXEDME
    // for (size_t cat = 0; cat < HDalitzCats_.size(); cat++){
    //     cout << HDalitzCats_[cat].c_str() << endl;
    // }

    for (size_t cat = 0; cat < HDalitzCats_.size(); cat++){
        map<string, int> choices;
        map<string, std::vector<int>> choices_envelope;
        map<string, RooAbsPdf *> pdfs;
        map<string, RooAbsPdf *> allPdfs;
        string catname_ = HDalitzCats_[cat].c_str();

        RooDataSet *dataFull = (RooDataSet *)inWS->data(Form("data_obs_%s", catname_.c_str()));
        cout << "[INFO] Got RooDataSet: " << dataFull->GetName() << ", Sum of entries: " << dataFull->sumEntries() << endl;

        mass->setBins(nBinsForMass);

        RooDataSet *data;
        string thisdataBinned_name = Form("roohist_data_mass_%s", catname_.c_str());
        
        RooDataHist thisdataBinned(thisdataBinned_name.c_str(), "data", *mass, *dataFull);
        data = (RooDataSet*) &thisdataBinned;

        RooArgList storedPdfs("store");

        string catName_split = ReplaceAll(catname_.c_str(), string("_"), string(" "));
        fprintf(resFile, "\\multicolumn{4}{|c|}{\\textbf{%s}} \\\\\n", catName_split.c_str());//!FIXEDME
        fprintf(resFile, "\\hline\n");

        double MinimimNLLSoFar = 1e10;
        int simplebestFitPdfIndex = 0;

        // Standard F-Test to find the truth functions
        for (vector<string>::iterator funcType = functionClasses.begin(); funcType != functionClasses.end(); funcType++){
            double thisNll = 0., prevNll = 0., chi2 = 0., prob = 0.;
            int order = 1, prev_order = 0, cache_order = 0;

            RooAbsPdf *prev_pdf = NULL;
            RooAbsPdf *cache_pdf = NULL;
            std::vector<int> pdforders;

            int counter = 0;

            while (prob < 0.05 && order < 7){ 
                RooAbsPdf *bkgPdf = getPdf(pdfsModel, *funcType, order, Form("ftest_pdf_%s_cat%d_%s", catname_.c_str(), (int)cat, ext.c_str())); //!FIXEDME

                if (!bkgPdf) order++; // assume this order is not allowed
                else{
                    int fitStatus = 0;
                    // bkgPdf->Print(); //! FIXEDME
                    cout << GREEN << "[INFO] Perform the F-test on: " << " category " << HDalitzCats_[cat] << ", " << bkgPdf->GetName() << RESET << endl;
                    runFit(bkgPdf, data, &thisNll, &fitStatus, /*max iterations*/ 3); 

                    if (fitStatus != 0){
                        cout << YELLOW << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus << RESET << endl;
                    }
                    chi2 = 2. * (prevNll - thisNll);

                    if (chi2 < 0. && order > 1) chi2 = 0.;

                    if (prev_pdf != NULL){
                        prob = getProbabilityFtest(chi2, order - prev_order, prev_pdf, bkgPdf, mass, data, Form("%s/Ftest_from_%s%d_%s_cat%d_%s.pdf", plotDir_.c_str(), funcType->c_str(), order, catname_.c_str(), (int)cat, ext.c_str()));//!FIXEDME

                        cout << "[INFO] F-test Prob(chi2 > chi2(data)) == " << prob << " with ndf = " << (order - prev_order) << endl;
                    }
                    else prob = 0;

                    cout << "[INFO] Function: " << *funcType << ", order: " << order << ", prevNll: " << prevNll << ", thisNll: " << thisNll << ", chi2: " << chi2 << ", prob: " << prob << endl;

                    prevNll = thisNll;
                    cache_order = prev_order;
                    cache_pdf = prev_pdf;
                    prev_order = order;
                    prev_pdf = bkgPdf;
                    order++;
                }
                counter += 1;
            }

            fprintf(resFile, "%15s & %d & %5.2f & %5.2f \\\\\n", funcType->c_str(), cache_order + 1, chi2, prob);
            choices.insert(pair<string, int>(*funcType, cache_order));
            pdfs.insert(pair<string, RooAbsPdf *>(Form("%s%d", funcType->c_str(), cache_order), cache_pdf));

            int truthOrder = cache_order;
            
            // Now run loop to determine functions inside envelope
            if (saveMultiPdf){
                chi2 = 0.;
				thisNll = 0.;
				prevNll = 0.;
				prob = 0.;
				order = 1;
				prev_order = 0;
				cache_order = 0;

                cout << "[INFO] Determining Envelope Functions for Family: " << *funcType << ", cat: " << catname_ << endl;
				cout << "[INFO] Upper end Threshold for highest order function " << upperEnvThreshold <<endl;

                while (prob < upperEnvThreshold){
                    RooAbsPdf *bkgPdf = getPdf(pdfsModel, *funcType, order, Form("env_pdf_%s_cat%d_%s", catname_.c_str(), (int)cat, ext.c_str())); //!FIXEDME

                    if (!bkgPdf){
                        if (order > 6){ 
                            cout << " [WARNING] could not add ] " << endl; 
                            break ;
                        }
						order++;
                    }
                    else{
                        //RooFitResult *fitRes;
						int fitStatus = 0;
						runFit(bkgPdf, data, &thisNll, &fitStatus, /*max iterations*/3);//bkgPdf->fitTo(*data,Save(true),RooFit::Minimizer("Minuit2","minimize"));
						//thisNll = fitRes->minNll();
						if (fitStatus != 0) cout << "[WARNING] Warning -- Fit status for " << bkgPdf->GetName() << " at " << fitStatus <<endl;
						
                        double myNll = 2.*thisNll;
						chi2 = 2. * (prevNll - thisNll);
						
                        if (chi2 < 0. && order > 1) chi2 = 0.;
						prob = TMath::Prob(chi2, order - prev_order);

                        cout << "[INFO] Function: " << *funcType << ", order: " << order << ", prevNll: " << prevNll << ", thisNll: " << thisNll << ", chi2: " << chi2 << ", prob: " << prob << endl;

						prevNll = thisNll;
						cache_order = prev_order;
						cache_pdf = prev_pdf;

                        // Calculate goodness of fit for the thing to be included (will use toys for lowstats)!
                        double gofProb = 0;
                        plot_1(mass, bkgPdf, data, Form("%s/single_fit/%s%d_%s_cat%d_%s", plotDir_.c_str(), funcType->c_str(), order, catname_.c_str(), (int)cat, ext.c_str()), fitStatus, &gofProb); //!FIXEDME

                        // Looser requirements for the envelope
                        if ((prob < upperEnvThreshold)){ 
                            // Good looking fit or one of our regular truth functions
                            if (gofProb > 0.01 || order == truthOrder){ 
                                cout << "[INFO] Adding to Envelope " << bkgPdf->GetName() << " " << gofProb
                                     << " 2xNLL + c is " << myNll + bkgPdf->getVariables()->getSize() 
                                     << endl;
                                allPdfs.insert(pair<string, RooAbsPdf *>(Form("%s%d", funcType->c_str(), order), bkgPdf));
                                storedPdfs.add(*bkgPdf);
                                pdforders.push_back(order);

                                // Keep track but we shall redo this later
                                if ((myNll + bkgPdf->getVariables()->getSize()) < MinimimNLLSoFar){
                                    simplebestFitPdfIndex = storedPdfs.getSize() - 1;
                                    MinimimNLLSoFar = myNll + bkgPdf->getVariables()->getSize();
                                }
                            }
                        }
                        prev_order = order;
                        prev_pdf = bkgPdf;
                        order++;
                    }
                }
                fprintf(resFile, "%15s & %d & %5.2f & %5.2f \\\\\n", funcType->c_str(), cache_order + 1, chi2, prob);
                choices_envelope.insert(pair<string, std::vector<int>>(*funcType, pdforders));
            }
        }

        // End of standard F-test
        fprintf(resFile, "\\hline\n");
        choices_vec.push_back(choices);
        choices_envelope_vec.push_back(choices_envelope);
        pdfs_vec.push_back(pdfs);
        plot_2(mass, pdfs, data, Form("%s/truths_%s_cat%d_%s", plotDir_.c_str(), HDalitzCats_[cat].c_str(), (int)cat, ext.c_str()), HDalitzCats_, cat); //!FIXEDME

        if (saveMultiPdf){
            // Put selectedModels into a MultiPdf
            string catindexname = Form("pdfindex_%s_cat%d_%s", catname_.c_str(), (int)cat, ext.c_str()); //!FIXEDME
            string catnameNum = Form("cat%d", (int)cat); //!FIXEDME
            
            RooCategory catIndex(catindexname.c_str(), "c");
            RooMultiPdf *pdf = new RooMultiPdf(Form("CMS_higgs_%s_%s_bkgshape", catname_.c_str(), ext.c_str()), "all pdfs", catIndex, storedPdfs);
            RooRealVar nBackground(Form("CMS_higgs_%s_%s_bkgshape_norm", catname_.c_str(), ext.c_str()), "nbkg", data->sumEntries(), 0, 3 * data->sumEntries());

            //double check the best pdf!
			int bestFitPdfIndex = getBestFitFunction(pdf, data, &catIndex, !verbosity_);
			catIndex.setIndex(bestFitPdfIndex);
			cout << "// ------------------------------------------------------------------------- //" <<endl; 
			cout << "[INFO] Created MultiPdf " << pdf->GetName() << ", in Category " << cat << " with a total of " << catIndex.numTypes() << " pdfs"<< endl;
			storedPdfs.Print();
			cout << "[INFO] Best Fit Pdf = " << bestFitPdfIndex << ", " << storedPdfs.at(bestFitPdfIndex)->GetName() << endl;
			cout << "// ------------------------------------------------------------------------- //" <<endl;
			cout << "[INFO] Simple check of index "<< simplebestFitPdfIndex <<endl;

            mass->setBins(nBinsForMass);
            RooDataHist dataBinned(Form("roohist_data_mass_%s", catname_.c_str()), "data", *mass, *dataFull);

            // Save it (also a binned version of the dataset
            outputws->import(*pdf);
            outputws->import(nBackground);
            outputws->import(catIndex);
            outputws->import(dataBinned);
            outputws->import(*data);
            outputws->import(*dataFull, Rename(Form("data_obs_cat%d", (int)cat))); //!FIXEDME
            plot_3(mass, pdf, &catIndex, data, Form("%s/multipdf_%s_%s", plotDir_.c_str(), catname_.c_str(), ext.c_str()), HDalitzCats_, cat, bestFitPdfIndex); //!FIXEDME

            cout << pdf->getCurrentPdf()->GetName() <<endl;
        }
    }
    if (saveMultiPdf){
        outputfile->cd();
        outputws->Write();
        outputfile->Close();
    }

    FILE *dfile = fopen(datfile_.c_str(), "w");
    cout << "[RESULT] Recommended options" << endl;

    for (size_t cat = 0; cat < HDalitzCats_.size(); cat++){
        cout << "Cat [" << HDalitzCats_[cat].c_str() << ", " << HDalitzCatsNum_[cat] << "]" << endl;
        
        fprintf(dfile, "cat = [%s, %d]\n", HDalitzCats_[cat].c_str(), (int)cat);
        for (map<string, int>::iterator it = choices_vec[cat].begin(); it != choices_vec[cat].end(); it++){
            cout << "\t" << it->first << " - " << it->second << endl;
            fprintf(dfile, "truth = %s:%d:%s%d\n", it->first.c_str(), it->second, namingMap[it->first].c_str(), it->second);
        }
        for (map<string, vector<int>>::iterator it = choices_envelope_vec[cat].begin(); it != choices_envelope_vec[cat].end(); it++){
            vector<int> ords = it->second;
            for (vector<int>::iterator ordit = ords.begin(); ordit != ords.end(); ordit++){
                fprintf(dfile, "paul = %s:%d:%s%d\n", it->first.c_str(), *ordit, namingMap[it->first].c_str(), *ordit);
            }
        }
        fprintf(dfile, "\n");
    }
    inFile->Close();

    return 0;
}