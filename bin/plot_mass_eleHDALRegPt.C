R__LOAD_LIBRARY($HDalitzEle_LOC/lib/libHDalitzEle.so)
R__ADD_INCLUDE_PATH($HDalitzEle_LOC/include)
#include "ROOT/RVec.hxx"
#include "ROOT/RDF/InterfaceUtils.hxx"
#include "ROOT/RDataFrame.hxx"
#include "Math/Vector4D.h"
#include "CMS_lumi_modification.h"
#include "tdrstyle.h"
using namespace std;


vector<const char*> cname = {"#E16262", "#202020", "#387c6d"};
vector<string> categories = {
    "Merged-2Gsf-HVBF",
    "Merged-2Gsf-LVBF",
    "Merged-2Gsf-Bst",
    "Merged-2Gsf-EBHR9",
    "Merged-2Gsf-EBLR9",
    "Merged-2Gsf-EE"
};


// Function to extract the sigma effective of a histogram
// https://github.com/cms-analysis/flashggFinalFit/blob/dev_fggfinalfits_lite/Background/test/plotweightedsigEd.cpp#L124
float effSigma(TH1* hist, double quantile = TMath::Erf(1.0/sqrt(2.0))){
    TAxis* xaxis = hist->GetXaxis();
    int nb = xaxis->GetNbins();
    if(nb < 10) {
        cout << "effsigma: Not a valid histo(too few). nbins = " << nb << endl;
        return 0.;
    }

    float bwid = xaxis->GetBinWidth(1);
    if(bwid == 0) {
        cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
        return 0.;
    }

    float xmax = xaxis->GetXmax();
    float xmin = xaxis->GetXmin();
    float ave = hist->GetMean();
    // float ave = hist->GetBinCenter(hist->GetMaximumBin()); // modified by CH 
    float rms = hist->GetRMS();
    float total = 0.;
    for(int i = 0; i < nb+2; i++) {
        total += hist->GetBinContent(i);
    }
    int ierr = 0;
    int ismin = 999;

    // Scan around window of mean: window RMS/binWidth (cannot be bigger than 0.1*number of bins)
    // Determine minimum width of distribution which holds 0.693 of total
    float rlim = quantile*total;
    int nrms = rms/bwid; // Set scan size to +/- rms
    if(nrms > nb/10) 
        nrms = nb/10; // Could be tuned...

    float widmin = 9999999.;
    for(int iscan = -nrms; iscan < nrms+1; iscan++){ // Scan window centre
        int ibm = (ave - xmin)/bwid + 1 + iscan; // Find bin idx in scan: iscan from mean
        float x = (ibm - 0.5)*bwid + xmin; // 0.5 for bin centre
        float xj = x;
        float xk = x;
        int jbm = ibm;
        int kbm = ibm;

        // Define counter for yield in bins: stop when counter > rlim
        float bin = hist->GetBinContent(ibm);
        total = bin;
        for(int j = 1; j < nb; j++){
            if(jbm < nb) {
                jbm++;
                xj += bwid;
                bin = hist->GetBinContent(jbm);
                total += bin;
                if(total > rlim) 
                    break;
            }
            else 
                ierr = 1;

            if(kbm > 0) {
                kbm--;
                xk -= bwid;
                bin = hist->GetBinContent(kbm);
                total += bin;
                if(total > rlim) 
                    break;
            }
            else 
                ierr = 2;
        }

        // Calculate fractional width in bin takes above limt (assume linear)
        float dxf = (total-rlim)*bwid/bin;
        float wid = (xj-xk+bwid-dxf)*0.5; // Total width: half of peak

        if(wid < widmin) {
            widmin = wid;
            ismin = iscan;
        }
    }

    if(ismin == nrms || ismin == -nrms) 
        ierr=3;
    if(ierr != 0) 
        cout << "effsigma: Error of type " << ierr << endl;
    
    return widmin;
}


void plot_mass_eleHDALRegPt(){

    vector<string> infiles_reg = {
        "/data4/chenghan/electron/miniTree_reg/UL2017/miniTree_HDalitz_ggF_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_reg/UL2017/miniTree_HDalitz_VBF_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_reg/UL2017/miniTree_HDalitz_WH_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_reg/UL2017/miniTree_HDalitz_ZH_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_reg/UL2017/miniTree_HDalitz_ttH_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_reg/UL2017/miniTree_HDalitz_bbH_eeg_125_UL2017.root"
    };
    auto df_reg = ROOT::RDataFrame("miniTree", infiles_reg).Define("plotwei", "mcwei*genwei*puwei");

    vector<string> infiles_noreg = {
        "/data4/chenghan/electron/miniTree_noreg/UL2017/miniTree_HDalitz_ggF_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_noreg/UL2017/miniTree_HDalitz_VBF_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_noreg/UL2017/miniTree_HDalitz_WH_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_noreg/UL2017/miniTree_HDalitz_ZH_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_noreg/UL2017/miniTree_HDalitz_ttH_eeg_125_UL2017.root",
        "/data4/chenghan/electron/miniTree_noreg/UL2017/miniTree_HDalitz_bbH_eeg_125_UL2017.root"
    };
    auto df_noreg = ROOT::RDataFrame("miniTree", infiles_noreg).Define("plotwei", "mcwei*genwei*puwei");
    
    // plot the hostograms
    setTDRStyle();
    for (size_t i = 0; i < categories.size(); i++){
        auto h_reg = df_reg.Filter(Form("category == %lu", i+1))
                           .Histo1D<float>({Form("h_reg_%lu", i+1), "", 90, 110, 140}, "CMS_higgs_mass", "plotwei");
        
        auto h_noreg = df_noreg.Filter(Form("category == %lu", i+1))
                               .Histo1D<float>({Form("h_noreg_%lu", i+1), "", 90, 110, 140}, "CMS_higgs_mass", "plotwei");

        cout << Form("Comparison for %s", categories[i].c_str()) << endl;
        auto sigma_noreg = effSigma((TH1*) h_noreg.GetPtr());
        auto scale_noreg = h_noreg->GetMean();
        cout << "      :Mean of Higgs mass: " << scale_noreg << endl;
        cout << "      :Effective sigma before regression: " << sigma_noreg << endl;

        auto sigma_reg = effSigma((TH1*) h_reg.GetPtr());
        auto scale_reg = h_reg->GetMean();
        cout << "      :Mean of Higgs mass: " << scale_reg << endl;
        cout << "      :Effective sigma after regression: " << sigma_reg << endl;

        h_reg->Scale(1./h_reg->Integral());
        h_noreg->Scale(1./h_noreg->Integral());
        
        TCanvas* canv = new TCanvas("c", "c", 800, 700);
        canv->cd();
        canv->SetRightMargin(0.05);
        canv->SetBottomMargin(0.14);
        canv->SetTopMargin(0.07);
        canv->SetLeftMargin(0.17);

        h_reg->GetXaxis()->SetTitle("M_{ee#gamma} [GeV]");
        h_reg->GetXaxis()->SetTickSize(0.02);
        h_reg->GetXaxis()->SetTitleSize(0.05);
        h_reg->GetXaxis()->SetLabelSize(0.05);
        h_reg->GetXaxis()->SetLabelOffset(0.02);
        h_reg->GetXaxis()->SetTitleOffset(1.3);
        h_reg->GetYaxis()->SetTitle("A.U.");
        h_reg->GetYaxis()->SetRangeUser(0, h_reg->GetBinContent(h_reg->GetMaximumBin()) * 1.5); 
        h_reg->GetYaxis()->SetTickSize(0.02);
        h_reg->GetYaxis()->SetTitleSize(0.05);
        h_reg->GetYaxis()->SetLabelSize(0.05);
        h_reg->GetYaxis()->SetTitleOffset(1.6);

        h_reg->SetLineColor(TColor::GetColor(cname[0]));
        h_reg->SetLineWidth(4);
        h_reg->Draw("hist same");

        h_noreg->SetLineColor(TColor::GetColor(cname[1]));
        h_noreg->SetLineWidth(4);
        h_noreg->Draw("hist same");
        h_reg->Draw("hist same");

        TLegend *leg = new TLegend(0.34, 0.73, 0.82, 0.91); //(x1,y1,x2,y2)
        leg->SetTextSize(0.04);
        leg->SetFillColor(0);
        leg->SetLineColor(0);
        leg->SetHeader(categories[i].c_str());
        leg->AddEntry(h_reg.GetPtr(),   Form("EGM calibration (#sigma_{eff} = %.2f GeV)", sigma_noreg), "l");
        leg->AddEntry(h_noreg.GetPtr(), Form("XGB regression (#sigma_{eff} = %.2f GeV)", sigma_reg), "l");
        leg->SetFillColor(0);
        leg->SetLineColor(0);
        leg->Draw();

        CMS_lumi(canv, 4, 10, "41.5 fb^{-1}", 2017, true, "Simulation", "H #rightarrow #gamma* #gamma #rightarrow ee#gamma", "");
    
        gSystem->Exec("mkdir -p ../plots/regression");
        canv->Print(Form("../plots/regression/regressionCmp_%s.pdf", categories[i].c_str()));
        canv->Close();

        cout << endl;
    }
}