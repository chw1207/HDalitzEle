// Profiles normalisation in mass bins over multiple pdfs
// Author: M Kenzie
//
// There are a few different algorithms written, the default
// and most robust is guessNew() but it is very slow!!!
// If you would like to write something smarter that would be great!!

#include "TF1.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TMath.h"
#include "TLine.h"
#include "TDirectory.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooAbsPdf.h"
#include "RooHist.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMacro.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "TROOT.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"
#include "RooCategory.h"
#include "../interface/WSTFileWrapper.h"

#include "boost/program_options.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/predicate.hpp"
#include "../interface/ProfileMultiplePdfs.h"

#include "HiggsAnalysis/CombinedLimit/interface/RooMultiPdf.h"
#include "HiggsAnalysis/CombinedLimit/interface/RooBernsteinFast.h"

#include <iostream>

#include "../../tdrStyle/tdrstyle.C"
#include "../../tdrStyle/CMS_lumi_mod.C"

#define RESET   "\033[0m"
#define GREEN   "\033[32m"
#define YELLOW  "\033[93m"
#define RED     "\033[91m"

#include <TSystem.h>

using namespace RooFit;
using namespace std;
using namespace boost;
namespace po = boost::program_options;

vector<string> HDalitzCats_;
string HDalitzCatsStr_ = "Merged2Gsf_HVBF,Merged2Gsf_LVBF,Merged2Gsf_BST,Merged2Gsf_EBHR9,Merged2Gsf_EBLR9,Merged2Gsf_EE,Merged1Gsf_HVBF,Merged1Gsf_LVBF,Merged1Gsf_BST,Merged1Gsf_EBHR9,Merged1Gsf_EBLR9,Merged1Gsf_EE,Resolved";

int nBinsForMass;
string bkgFileName = "./MultiPdf/Multipdf.root";
string sigFileName;
string outFileName;
string outDir;
int cat = 3;
string catLabel;
double massStep;
double nllTolerance;
bool doBands = true;
bool useBinnedData = true;
bool doSignal = false;
bool unblind = false;
bool verbose_ = false;
bool makeCrossCheckProfPlots=false;
int mhLow;
int mhHigh;
double mhvalue_;
double higgsResolution_;
string ext;
int sqrts; 
int year_ = 2016;
//int year_=2017;
float intLumi = 137.1;
string massName = "M_{ee#gamma} [GeV]"; // use to plot the x title //!changed

RooRealVar *intLumi_ = new RooRealVar("IntLumi","hacked int lumi", 1000.);

vector<int> unblind_bond = {120, 130};

void OptionParser(int argc, char *argv[]){
    po::options_description desc("Allowed options");
    desc.add_options()
    ("help,h",          "Show help")
    ("bkgFileName,b",   po::value<string>(&bkgFileName)->default_value(bkgFileName), 				    "Input file name")
	// ("sigFileName,s",   po::value<string>(&sigFileName), 								                "Input file name")
    ("outDir,d",        po::value<string>(&outDir)->default_value("BkgPlots"),						    "Output directory")
	("outFileName,o",   po::value<string>(&outFileName)->default_value("BkgPlots.root"),	            "Output file name")
    ("cat,c",           po::value<int>(&cat),											                "Category")
    ("massStep,m",      po::value<double>(&massStep)->default_value(1),				                    "Mass step for calculating bands")
    ("nllTolerance,n",  po::value<double>(&nllTolerance)->default_value(0.05),			                "Tolerance for nll calc in %")
	("mhLow,L",         po::value<int>(&mhLow)->default_value(105),						                "Starting point for scan")
	("mhHigh,H",        po::value<int>(&mhHigh)->default_value(145),						            "End point for scan")
    ("mhVal",           po::value<double>(&mhvalue_)->default_value(125.),			                    "Choose the MH for the plots")
    ("higgsResolution", po::value<double>(&higgsResolution_)->default_value(1.),		                "Starting point for scan")
    ("intLumi",         po::value<float>(&intLumi)->default_value(0.),						            "What intLumi in fb^{-1}")
    ("year",            po::value<int>(&year_)->default_value(2016),						            "Dataset year")
	("sqrts,S",         po::value<int>(&sqrts)->default_value(13),							            "Which centre of mass is this data from?")
    ("HDalitzCats",     po::value<string>(&HDalitzCatsStr_)->default_value(HDalitzCatsStr_.c_str()),    "Higgs Dalitz category")
    ("extratext",       po::value<string>(&ext)->default_value("DiPho_13TeV"),                         "Higgs Dalitz category")
    
    ("doBands",			"Do error bands")
    ("unblind",		    "un blind central mass region")
    ("useBinnedData",	"Data binned")
    ("verbose,v", 		"Verbose");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help")) { 
        cout << desc << endl; 
        exit(1); 
    }
	if (vm.count("doBands")) doBands = true;
	if (vm.count("makeCrossCheckProfPlots")) makeCrossCheckProfPlots = true;
	if (vm.count("unblind")) unblind = true;
	if (vm.count("useBinnedData")) useBinnedData = true;
	if (vm.count("sigFileName")) doSignal = true;
	if (vm.count("verbose")) verbose_ = true;


    if (!verbose_){
        RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
        RooMsgService::instance().setSilentMode(true);
        // gErrorIgnoreLevel = kWarning;
    }
    split(HDalitzCats_, HDalitzCatsStr_, boost::is_any_of(","));

    nBinsForMass = (int)(mhHigh - mhLow); //!FIXEDME
}

int getBestFitFunction(RooMultiPdf *bkg, RooAbsData *data, RooCategory *cat, bool silent = false){
    double global_minNll = 1E10;
	int best_index = 0;
	int number_of_indeces = cat->numTypes();

    RooArgSet snap,clean;
	RooArgSet *params = bkg->getParameters(*data);
	params->snapshot(snap);
	params->snapshot(clean);
	if (!silent) {
		cout << "[INFO] CLEAN SET OF PARAMETERS" << endl;
		params->Print("V");
		cout << "-----------------------" << endl;
	}

    //bkg->setDirtyInhibit(1);
	RooAbsReal *nllm = bkg->createNLL(*data);
	RooMinimizer minim(*nllm);
	minim.setStrategy(2);

    for (int id = 0;id < number_of_indeces; id++){		
		params->assignValueOnly(clean);
		cat->setIndex(id);

		//RooAbsReal *nllm = bkg->getCurrentPdf()->createNLL(*data);
		if (!silent) {
			cout << "[INFO] BEFORE FITTING" << endl;
			params->Print("V");
			cout << "-----------------------" << endl;
		}
		
		minim.minimize("Minuit2","simplex");
		double minNll = nllm->getVal()+bkg->getCorrection();
		if (!silent) {
			cout << "[INFO] After Minimization ------------------  " <<endl;
			cout << "[INFO] " << bkg->getCurrentPdf()->GetName() << " " << minNll <<endl;
			bkg->Print("v");
			bkg->getCurrentPdf()->getParameters(*data)->Print("V");
			cout << " ------------------------------------  " << endl;
	
			cout << "[INFO] AFTER FITTING" << endl;
			params->Print("V");
			cout << "-----------------------" << endl;
		}
			
		if (minNll < global_minNll){
        		global_minNll = minNll;
			snap.assignValueOnly(*params);
        		best_index=id;
		}
	}

    params->assignValueOnly(snap);
    cat->setIndex(best_index);
	
	if (!silent) {
		cout << "[INFO] Best fit Function -- " << bkg->getCurrentPdf()->GetName() << " " << cat->getIndex() <<endl;
		bkg->getCurrentPdf()->getParameters(*data)->Print("v");
	}
	return best_index;
}


pair<double, double> getNormTermNllAndRes(RooRealVar *mass, RooAbsData *data, RooMultiPdf *mpdf, RooCategory *mcat, double normVal = -1., double massRangeLow = -1., double massRangeHigh = -1.){
    double bestFitNll = 1.e8;
	double bestFitNorm;

    for (int pInd = 0; pInd < mpdf->getNumPdfs(); pInd++){
        mcat->setIndex(pInd);
		RooRealVar *normVar = new RooRealVar(Form("%snorm", mpdf->getCurrentPdf()->GetName()), " " , 0., 1.e6);
		RooExtendPdf *extPdf;
		RooAbsReal *nll;

        if (massRangeLow > -1. && massRangeHigh > -1.){
            mass->setRange("errRange", massRangeLow, massRangeHigh);
			extPdf = new RooExtendPdf(Form("%sext", mpdf->getCurrentPdf()->GetName()), " ", *(mpdf->getCurrentPdf()), *normVar, "errRange");
			nll = extPdf->createNLL(*data, Extended());
        }
        else {
			extPdf = new RooExtendPdf(Form("%sext", mpdf->getCurrentPdf()->GetName()), " ", *(mpdf->getCurrentPdf()), *normVar);
			nll = extPdf->createNLL(*data, Extended());
		}

        if (normVal > -1.){
			normVar->setConstant(false);
			normVar->setVal(normVal);
			normVar->setConstant(true);
		}

        RooMinimizer minim(*nll);
		minim.setStrategy(0);
		//minim.minimize("Minuit2","simplex");
		minim.migrad();
		double corrNll = nll->getVal() + mpdf->getCorrection();

        if (corrNll < bestFitNll){
			bestFitNll = corrNll;
			bestFitNorm = normVar->getVal();
		}
		if (normVal > -1.) normVar->setConstant(false);
    }

    return make_pair(2*bestFitNll, bestFitNorm);
}



double getNormTermNll(RooRealVar *mass, RooAbsData *data, RooMultiPdf *mpdf, RooCategory *mcat, double normVal = -1., double massRangeLow = -1., double massRangeHigh = -1.){
	pair<double, double> temp = getNormTermNllAndRes(mass, data, mpdf, mcat, normVal, massRangeLow, massRangeHigh);
	return temp.first;
}


double guessNew(RooRealVar *mass, RooMultiPdf *mpdf, RooCategory *mcat, RooAbsData *data, double bestPoint, double nllBest, double boundary, double massRangeLow, double massRangeHigh, double crossing, double tolerance){
    bool isLowSide;
	double guess, guessNll, lowPoint, highPoint;
	if (boundary > bestPoint) {
		isLowSide = false;
		lowPoint = bestPoint;
		highPoint = boundary;
	}
	else {
		isLowSide = true;
		lowPoint = boundary;
		highPoint = bestPoint;
	}

	//double prevDistanceFromTruth = 1.e6;
	double distanceFromTruth = 1.e6;
	int nIts = 0;

    while (TMath::Abs(distanceFromTruth / crossing) > tolerance){
        //prevDistanceFromTruth = distanceFromTruth;
		guess = lowPoint + (highPoint - lowPoint) / 2.;
		guessNll = getNormTermNll(mass, data, mpdf, mcat, guess, massRangeLow, massRangeHigh) - nllBest;
        distanceFromTruth = crossing - guessNll;

        if (verbose_) {
			cout << "[INFO] "<< Form("\t lP: %7.3f hP: %7.3f xg: %7.3f yg: %7.3f", lowPoint, highPoint, guess, guessNll) << endl;;
			cout << "[INFO] \t ----------- " << distanceFromTruth/crossing << " -------------" << endl;
		}

        // for low side. if nll val is lower than target move right point left. if nll val is higher than target move left point right
		// vice versa for high side
		if (isLowSide){
			if (guessNll > crossing) 
                lowPoint = guess;
			else 
                highPoint = guess;
		}
		else {
			if (guessNll > crossing) 
                highPoint = guess;
			else 
                lowPoint = guess;
		}
		nIts++;

        // because these are envelope nll curves this algorithm can get stuck in local minima
		// hacked get out is just to start trying again
		if (nIts > 20) {
			return guess;
			lowPoint = TMath::Max(0.,lowPoint - 20);
			highPoint += 20;
			nIts = 0;
			if (verbose_) 
                cout << "[INFO] RESET:" << endl;
			// if it's a ridicolous value it wont converge so return value of bestGuess
			if (TMath::Abs(guessNll) > 2e4) 
                return 0.; 
		}
    }

    return guess;
}


int main(int argc, char *argv[]){
    OptionParser(argc, argv);

    setTDRStyle();
    
    system(Form("mkdir -p %s", outDir.c_str()));
	if (makeCrossCheckProfPlots) 
        system(Form("mkdir -p %s/normProfs", outDir.c_str()));

    cout << GREEN << "----> " << "Open the background file and extract ws:" << RESET << endl;
    TFile *inFile = TFile::Open(bkgFileName.c_str());
    WSTFileWrapper * inWS = new WSTFileWrapper(bkgFileName, "multipdf");
    if (!inWS) {
		cout << RED << "[ERROR] "<< "Cant find the workspace" << RESET <<endl;
		exit(0);
	}

    RooRealVar *mass = (RooRealVar *)inWS->var("CMS_higgs_mass");
    cout << "[INFO] Got mass from ws: " << mass->GetName() << endl;

    string catname = HDalitzCats_[cat].c_str();
    TFile *outFile = TFile::Open(outFileName.c_str(), "RECREATE");
	RooWorkspace *outWS = new RooWorkspace("bkgplotws", "bkgplotws");

    RooAbsData *data = (RooDataSet*)inWS->data(Form("data_obs_cat%d", cat));
	if (useBinnedData) 
        data = (RooDataHist*)inWS->data(Form("roohist_data_mass_%s", catname.c_str()));

	RooMultiPdf *mpdf = (RooMultiPdf*)inWS->pdf(Form("CMS_higgs_%s_%s_bkgshape", catname.c_str(), ext.c_str())); 
	RooCategory *mcat = (RooCategory*)inWS->cat(Form("pdfindex_%s_cat%d_%s", catname.c_str(), cat, ext.c_str()));
    if (!mpdf || !mcat){
        cout << RED << "[ERROR] " << "Can't find multipdfs (" << Form("CMS_higgs_%s_%s_bkgshape", catname.c_str(), ext.c_str()) << ") or multicat (" << Form("pdfindex_%s_cat%d_%s", catname.c_str(), cat, ext.c_str()) << ")" << RESET << endl;
        exit(0);
    }

    cout << GREEN << "----> " << "Current PDF and data:" << RESET << endl;
	cout << "[INFO] " << "\t"; 
    mpdf->getCurrentPdf()->Print();
	cout << "[INFO] " << "\t"; 
    data->Print();

    // reset to best fit
	int bf = getBestFitFunction(mpdf, data, mcat, !verbose_);
	mcat->setIndex(bf);
    cout << GREEN << "----> " << "Best fit PDF and data:" << RESET << endl;
	cout<< "[INFO] " << "\t"; 
    mpdf->getCurrentPdf()->Print();
	cout<< "[INFO] " << "\t"; 
    data->Print();  
  
    // plot the data
	TLegend *leg = new TLegend(0.6, 0.6, 0.89, 0.89);
	leg->SetFillColor(0);
	leg->SetLineColor(0);
    RooPlot *plot = mass->frame();
	RooPlot *plotLC = mass->frame();
	plot->GetXaxis()->SetTitle(massName.c_str());
	plot->SetTitle("");
	data->plotOn(plot, Binning(nBinsForMass), Invisible());

	cout << GREEN << "----> " << "Plotting data and nominal curve" << RESET << endl;
    ///start extra bit for ratio plot///
    RooHist *plotdata = (RooHist*)plot->getObject(plot->numItems()-1);
    // enf extra bit for ratio plot///
	TObject *dataLeg = (TObject*)plot->getObject(plot->numItems()-1);
	mpdf->getCurrentPdf()->plotOn(plot,LineColor(kRed),LineWidth(2));
	RooCurve *nomBkgCurve = (RooCurve*)plot->getObject(plot->numItems()-1);

	leg->AddEntry(dataLeg, "Data", "LEP");
	leg->AddEntry(nomBkgCurve, "Bkg Fit", "L");

    // Bands
	TGraphAsymmErrors *oneSigmaBand = new TGraphAsymmErrors();
	TGraphAsymmErrors *oneSigmaBand_r = new TGraphAsymmErrors();
	oneSigmaBand->SetName(Form("onesigma_%s", catname.c_str()));
	oneSigmaBand_r->SetName(Form("onesigma_%s_r", catname.c_str()));
	TGraphAsymmErrors *twoSigmaBand = new TGraphAsymmErrors();
	TGraphAsymmErrors *twoSigmaBand_r = new TGraphAsymmErrors();
	twoSigmaBand->SetName(Form("twosigma_%s", catname.c_str()));
	twoSigmaBand_r->SetName(Form("twosigma_%s_r", catname.c_str()));

    cout<< "[INFO] " << "Plot has " << plot->GetXaxis()->GetNbins() << " bins" << endl;
    if (doBands){
        int p = 0;
        for (double m = double(mhLow); m < double(mhHigh) + massStep; m += massStep){
            double lowedge = m - 0.5;
			double upedge = m + 0.5;
			double center = m;
            double nomBkg = nomBkgCurve->interpolate(center);
            double nomBkg_perGeV = (nomBkgCurve->average(center - higgsResolution_, center + higgsResolution_));
            double nllBest = getNormTermNll(mass, data, mpdf, mcat, nomBkg, lowedge, upedge);

            // sensible range
			double lowRange = TMath::Max(0.,nomBkg - 3*TMath::Sqrt(nomBkg));
			double highRange = nomBkg + 3*TMath::Sqrt(nomBkg);
            // cout << "[FOR TABLE] ,"<< HDalitzCats_[cat].c_str() << "," << m << "," << nomBkg_perGeV << ", assuming resolution of " << higgsResolution_ << endl;
            // if (verbose_) 
            // cout<< "[INFO] " << "mass: " << center << " nomBkg: " << nomBkg << " lR: " << lowRange << " hR: " << highRange << endl;

            double errLow1Value, errHigh1Value, errLow2Value, errHigh2Value;
            
            // cant handle having 0 events
			if (nomBkg < 1.e-5) {
				errLow1Value = 0.;
				errLow2Value = 0.;
				
                if (verbose_) 
                    cout << "[INFO] "<< "errHigh1" << endl;
				errHigh1Value = guessNew(mass, mpdf, mcat, data, nomBkg, nllBest, highRange, lowedge, upedge, 1., nllTolerance);
				
                if (verbose_) 
                    cout << "[INFO] "<< "errHigh2" << endl;
				errHigh2Value = guessNew(mass, mpdf, mcat, data, nomBkg, nllBest, highRange, lowedge, upedge, 4., nllTolerance);
			}

            double errLow1 = nomBkg - errLow1Value;
			double errHigh1 = errHigh1Value - nomBkg;
			double errLow2 = nomBkg - errLow2Value;
			double errHigh2 = errHigh2Value - nomBkg;

			oneSigmaBand->SetPoint(p, center, nomBkg);
			oneSigmaBand_r->SetPoint(p, center, 0);
			twoSigmaBand->SetPoint(p, center, nomBkg);
			twoSigmaBand_r->SetPoint(p, center, 0);
			oneSigmaBand->SetPointError(p, 0., 0., errLow1, errHigh1);
			oneSigmaBand_r->SetPointError(p, 0., 0., errLow1, errHigh1);
			twoSigmaBand->SetPointError(p, 0., 0., errLow2, errHigh2);
			twoSigmaBand_r->SetPointError(p, 0., 0., errLow2, errHigh2);

			cout<< "[INFO] " << "mass: " << center << " nomBkg: " << nomBkg << " +/- 1 (2) sig -- +" << errHigh1 << "(" << errHigh2 << ")" << " - " << errLow1 << "(" << errLow2 << ")" << endl;

            p++;
        }
        cout << endl;
    }

    // outFile->cd();
    // oneSigmaBand->Write();
    // twoSigmaBand->Write();
    // nomBkgCurve->Write();
    // outWS->import(*mcat);
    // outWS->import(*mpdf);
    // outWS->import(*data);

    TCanvas *canv = new TCanvas("c", " ", 800, 800);
    ///start extra bit for ratio plot///
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

    pad1->Draw();
    pad2->Draw();
    pad1->cd();

    // enf extra bit for ratio plot///
    canv->SetTickx(); canv->SetTicky();
    //RooRealVar *lumi = (RooRealVar*)inWS->var("IntLumi");
    RooRealVar *lumi = intLumi_;
    plot->Draw();

    mass->setRange("unblind_up", unblind_bond[1], mhHigh);
    mass->setRange("unblind_down", mhLow, unblind_bond[0]);
    data->plotOn(plot, Binning(nBinsForMass), CutRange("unblind_down, unblind_up"));

    // if (!unblind) {
    //     mass->setRange("unblind_up", unblind_bond[1], mhHigh);
    //     mass->setRange("unblind_down", mhLow, unblind_bond[0]);
    //     data->plotOn(plot, Binning(nBinsForMass), CutRange("unblind_down, unblind_up"));
    // }
    // else {
    //     data->plotOn(plot,Binning(nBinsForMass));
    // }

    // if (doBands) {
    //     twoSigmaBand->SetLineColor(kYellow);
    //     twoSigmaBand->SetFillColor(kYellow);
    //     twoSigmaBand->SetMarkerColor(kYellow);
    //     twoSigmaBand->Draw("L3 SAME");
    //     oneSigmaBand->SetLineColor(kGreen);
    //     oneSigmaBand->SetFillColor(kGreen);
    //     oneSigmaBand->SetMarkerColor(kGreen);
    //     oneSigmaBand->Draw("L3 SAME");
    //     leg->AddEntry(oneSigmaBand,"#pm1#sigma","F");
    //     leg->AddEntry(twoSigmaBand,"#pm2#sigma","F");
    //     twoSigmaBand_r->SetLineColor(kYellow);
    //     twoSigmaBand_r->SetFillColor(kYellow);
    //     twoSigmaBand_r->SetMarkerColor(kYellow);
    //     oneSigmaBand_r->SetLineColor(kGreen);
    //     oneSigmaBand_r->SetFillColor(kGreen);
    //     oneSigmaBand_r->SetMarkerColor(kGreen);
    // }

    twoSigmaBand->SetLineColor(kYellow);
    twoSigmaBand->SetFillColor(kYellow);
    twoSigmaBand->SetMarkerColor(kYellow);
    twoSigmaBand->Draw("L3 SAME");
    oneSigmaBand->SetLineColor(kGreen);
    oneSigmaBand->SetFillColor(kGreen);
    oneSigmaBand->SetMarkerColor(kGreen);
    oneSigmaBand->Draw("L3 SAME");
    leg->AddEntry(oneSigmaBand,"#pm1#sigma","F");
    leg->AddEntry(twoSigmaBand,"#pm2#sigma","F");
    twoSigmaBand_r->SetLineColor(kYellow);
    twoSigmaBand_r->SetFillColor(kYellow);
    twoSigmaBand_r->SetMarkerColor(kYellow);
    oneSigmaBand_r->SetLineColor(kGreen);
    oneSigmaBand_r->SetFillColor(kGreen);
    oneSigmaBand_r->SetMarkerColor(kGreen);

    plot->Draw("same");
    leg->Draw("same");

    TLatex *latex = new TLatex();	
    latex->SetTextSize(0.045);
    latex->SetNDC();
    TLatex *cmslatex = new TLatex();
    cmslatex->SetTextSize(0.03);
    cmslatex->SetNDC();

    TString catLabel_humanReadable = catname;
    catLabel_humanReadable.ReplaceAll("_", " ");
    latex->DrawLatex(0.15,0.85,catLabel_humanReadable);
    outWS->import(*lumi,RecycleConflictNodes());

    // if (unblind) 
    //     plot->SetMinimum(0.0001);
    plot->GetYaxis()->SetTitleOffset(1.3);
    canv->Modified();
    CMS_lumi(pad1, 4, 0);
    canv->Update();
    canv->cd();

    // int npoints = plotdata->GetN();
    // double xtmp, ytmp;
    // int point = 0;
    // TGraphAsymmErrors *hdatasub = new TGraphAsymmErrors(npoints);
    // //hdatasub->SetMarkerSize(defmarkersize);
    // for (int ipoint = 0; ipoint < npoints; ipoint++) {
    //     //double bkgval = hbplottmp->GetBinContent(ipoint+1);
    //     plotdata->GetPoint(ipoint, xtmp, ytmp);
    //     double bkgval = nomBkgCurve->interpolate(xtmp);
    //     if (!unblind) {
    //         if ((xtmp > unblind_bond[0]) && ( xtmp < unblind_bond[1])) continue;
    //     }
    //     //std::cout << "[INFO] plotdata->Integral() " <<  plotdata->Integral() << " ( bins " << npoints  << ") hbkgplots[i]->Integral() " << hbplottmp->Integral() << " (bins " << hbplottmp->GetNbinsX() << std::endl;
    //     double errhi = plotdata->GetErrorYhigh(ipoint);
    //     double errlow = plotdata->GetErrorYlow(ipoint);
            
    //     //std::cout << "[INFO]  Channel " << name  << " errhi " << errhi << " errlow " << errlow  << std::endl;
    //     cout << "[INFO] Channel  " << " setting point " << point << " : xtmp " << xtmp << "  ytmp " << ytmp << " bkgval  " << bkgval << " ytmp-bkgval " << ytmp-bkgval << endl;
    //     bool drawZeroBins_ = 1;
    //     if (!drawZeroBins_) 
    //         if(fabs(ytmp) < 1e-5) continue; 
    //     hdatasub->SetPoint(point, xtmp, ytmp-bkgval);
    //     hdatasub->SetPointError(point, 0., 0., errlow, errhi);
    //     point++;
    // }

    pad2->cd();


    // TH1 *hdummy = new TH1D("hdummyweight", " ", nBinsForMass, mhLow, mhHigh);
    // hdummy->SetMaximum(hdatasub->GetHistogram()->GetMaximum()+1);
    // hdummy->SetMinimum(hdatasub->GetHistogram()->GetMinimum()-1);
    // hdummy->GetYaxis()->SetTitle("data - best fit PDF");
    // hdummy->GetYaxis()->SetTitleSize(0.12);
    // hdummy->GetXaxis()->SetTitle(massName.c_str());
    // hdummy->GetXaxis()->SetTitleSize(0.12);
    // hdummy->Draw("HIST");
    //     if (doBands) twoSigmaBand_r->Draw("L3 SAME");
    //     if (doBands) oneSigmaBand_r->Draw("L3 SAME");
    // hdummy->GetYaxis()->SetNdivisions(808);

    // TLine *line3 = new TLine(100,0.,180,0.);
    // line3->SetLineColor(kRed);
    // //line3->SetLineStyle(kDashed);
    // line3->SetLineWidth(4.0);
    // line3->Draw();
    // hdatasub->Draw("PESAME");
    canv->Print(Form("%s/bkgplot_%s.pdf",outDir.c_str(),catname.c_str()));
    // canv->Print(Form("%s/bkgplot_%s.png",outDir.c_str(),catname.c_str()));
    // canv->Print(Form("%s/bkgplot_%s.C",outDir.c_str(),catname.c_str()));
    // canv->SetName(Form("bkgplot_%s",catname.c_str()));
    // outFile->cd();
    // canv->Write();
    // outWS->Write();
    outFile->Close();

    inFile->Close();

    return 0; 
	
}