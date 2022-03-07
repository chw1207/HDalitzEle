/* Class which calculates event weights for MC due to difference in pileup
 * distributions between data and MC.
 */

#ifndef PUWEICALC_H
#define PUWEICALC_H

#include <map>
#include <TH1D.h>
#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

// prints a message and exits gracefully
#ifndef FATAL
#define FATAL(msg) do { fprintf(stderr, "FATAL: %s\n", msg); gSystem->Exit(1); } while (0)
#endif

class PUWeightCalculator {
public:
   virtual ~PUWeightCalculator();

   void Init(const char* path);
   double GetWeight(Int_t run, float puTrue);

protected:
   std::map<Int_t,TH1*> fRunHistoMap;  // {run number: weights} map
};

//______________________________________________________________________________
PUWeightCalculator::~PUWeightCalculator(){
   // Frees occupied memory.

   for (std::map<Int_t,TH1*>::iterator it = fRunHistoMap.begin(); it != fRunHistoMap.end(); it++) {
      TH1* histo = it->second;
      delete histo;
   }
}

//______________________________________________________________________________
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
      FATAL("TFile::Open() failed");

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
         FATAL("unexpected name for object in file");
      Int_t run = name.Atoi();

      // load the histogram
      TH1* histo = dynamic_cast<TH1*>(f.Get(obj->GetName()));
      if (!histo)
         FATAL("TFile::Get() or dynamic_cast<TH1*> failed");

      histo->SetDirectory(0);
      fRunHistoMap[run] = histo;
   }
}

//______________________________________________________________________________
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
      FATAL(Form("no weight histogram for run number %i", run));

   TH1* histo = got->second;

   Int_t bin = histo->FindBin(puTrue);
   if (bin < 1 || bin > histo->GetNbinsX())
      return 1;

   double wei = histo->GetBinContent(bin);
   if (wei < 1e-6) return 1;

   return wei;
}

#endif

// Added function by CH
std::string PUfile(int year, std::string choice){
   
   if (choice == "nominal"){
      if (year == 2016){
         return "../../external/puweights/94X/PU_histo_13TeV_2016_GoldenJSON_69200nb.root";
      }
      else if (year == 2017){
         return "../../external/puweights/94X/PU_histo_13TeV_2017_GoldenJSON_69200nb.root";
      }
      else if (year == 2018){
         return "../../external/puweights/102X/PU_histo_13TeV_2018_GoldenJSON_69200nb.root";
      }
      else{
         printf("ERROR: No specific year for PU file!\n");
         exit(-1);
      }
   }
   else if (choice == "up"){
      if (year == 2016){
         return "../../external/puweights/94X/PU_histo_13TeV_2016_GoldenJSON_72400nb.root";
      }
      else if (year == 2017){
         return "../../external/puweights/94X/PU_histo_13TeV_2017_GoldenJSON_72400nb.root";
      }
      else if (year == 2018){
         return "../../external/puweights/102X/PU_histo_13TeV_2018_GoldenJSON_72383nb.root";
      }
      else{
         printf("ERROR: No specific year for PU-up file!\n");
         exit(-1);
      }
   }
   else if (choice == "down"){
      if (year == 2016){
         return "../../external/puweights/94X/PU_histo_13TeV_2016_GoldenJSON_66000nb.root";
      }
      else if (year == 2017){
         return "../../external/puweights/94X/PU_histo_13TeV_2017_GoldenJSON_66000nb.root";
      }
      else if (year == 2018){
         return "../../external/puweights/102X/PU_histo_13TeV_2018_GoldenJSON_66016nb.root";
      }
      else{
         printf("ERROR: No specific year for PU-down file!\n");
         exit(-1);
      }
   }
   else {
      printf("ERROR: No specific choice for PU file!\n");
      exit(-1);
   }
}
