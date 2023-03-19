#include "PUWeightCalculator.h"
#include "TFile.h"
#include "TSystem.h"

PUWeightCalculator::~PUWeightCalculator(){
    // Frees occupied memory.

    for (std::map<Int_t, TH1*>::iterator it = fRunHistoMap.begin(); it != fRunHistoMap.end(); it++){
        TH1* histo = it->second;
        delete histo;
    }
}

void PUWeightCalculator::Init(const char* path){
    /* Reads in "event weight vs pileup" histograms.
    *
    * path = path to a weight file.
    */

    // current working directory
    TDirectory* wd = gDirectory;

    // open file with histograms
    TFile f(path);
    if (f.IsZombie())
        throw std::runtime_error("TFile::Open() failed");

    // cd back into previous current working directory
    if (wd) wd->cd();
    else gDirectory = 0;

    TIter inext(f.GetListOfKeys());

    // loop over objects in file
    while (TObject* obj = inext()) {
        TString name(obj->GetName());

        if (!name.BeginsWith("mcwei_run"))
            continue;

        // remove "mcwei_run"
        name.Replace(0, 9, "");

        // extract run number
        if (!name.IsDigit())
            throw std::runtime_error("unexpected name for object in file");
        Int_t run = name.Atoi();

        // load the histogram
        TH1* histo = dynamic_cast<TH1*>(f.Get(obj->GetName()));
        if (!histo)
            throw std::runtime_error("TFile::Get() or dynamic_cast<TH1*> failed");

        histo->SetDirectory(0);
        fRunHistoMap[run] = histo;
    }
}


double PUWeightCalculator::GetWeight(Int_t run, float puTrue){
    /* Evaluates and returns event weight.
    *
    * NOTE: if event weight could not be determined, zero is returned.
    *
    * run = simulated run number;
    * puTrue = in-time pileup value.
    */

    // search for a histogram corresponding to the run number "run"
    std::map<Int_t, TH1*>::const_iterator got = fRunHistoMap.find(run);

    // if not found
    if (got == fRunHistoMap.end())
        throw std::runtime_error(Form("no weight histogram for run number %i", run));

    TH1* histo = got->second;

    Int_t bin = histo->FindBin(puTrue);
    if (bin < 1 || bin > histo->GetNbinsX())
        return 1;

    double wei = histo->GetBinContent(bin);
    if (wei < 1e-6) return 1;

    return wei;
}