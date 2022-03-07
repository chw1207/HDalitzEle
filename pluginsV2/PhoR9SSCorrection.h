#include "TMVA/Reader.h"
#include "TFormula.h"

// Referece: 
// Hgg shower shape corrections are applied to mitigate discrepancies between data and simulation
// [1] https://arxiv.org/pdf/1804.02716.pdf

// This script is the simplified version of Ming's PhotonSSCorrection.h 
// This function would only return the corrected photonR95x5

float doR9SSCorrections(float phoEt, float phoPhi, float phoSCEta, float phoSCEtaWidth, float phoSCPhiWidth, float phoSigmaIEtaIEtaFull5x5, float phoSigmaIEtaIPhiFull5x5, float phoR9Full5x5, float phoE2x2Full5x5, float phoE5x5Full5x5, float rho, int year){
        
    const char* dir = "./external/models/tmvaWeightFiles";

    float corrR9 = 0;

    static TMVA::Reader *tmvaReaderR9[2] = {NULL, NULL};
    static TFormula *fR9[2] = {NULL, NULL};

    // MVA variables
	static float phoEt_, phoPhi_, phoR9Full5x5_;
	static float phoSCEtaWidth_, phoSCPhiWidth_, rho_;
	static float phoSCEta_;
	static float sieieFull5x5, sieipFull5x5, s4Full5x5;

    // 0 = ECAL barrel or 1 = ECAL endcaps
	int iBE = (fabs(phoSCEta) < 1.479) ? 0 : 1;

    // set MVA variables
    phoEt_ = phoEt;
    phoSCEta_ = phoSCEta;
    phoPhi_ = phoPhi;
	rho_ = rho;
    phoSCPhiWidth_ = phoSCPhiWidth;
    sieipFull5x5 = phoSigmaIEtaIPhiFull5x5;
    s4Full5x5 = phoE2x2Full5x5 / phoE5x5Full5x5;
    phoR9Full5x5_ = phoR9Full5x5;
    sieieFull5x5 = phoSigmaIEtaIEtaFull5x5;
    phoSCEtaWidth_ = phoSCEtaWidth;

    if (!tmvaReaderR9[iBE])
    {
        tmvaReaderR9[iBE] = new TMVA::Reader("!Color:Silent");

        if (year == 2016)
        {
            tmvaReaderR9[iBE]->AddVariable("f0", &phoEt_);
            tmvaReaderR9[iBE]->AddVariable("f1", &phoSCEta_);
            tmvaReaderR9[iBE]->AddVariable("f2", &phoPhi_);
            tmvaReaderR9[iBE]->AddVariable("f3", &rho_);
            tmvaReaderR9[iBE]->AddVariable("f4", &phoSCPhiWidth_);
            tmvaReaderR9[iBE]->AddVariable("f5", &sieipFull5x5);
            tmvaReaderR9[iBE]->AddVariable("f6", &s4Full5x5);
            tmvaReaderR9[iBE]->AddVariable("f7", &phoR9Full5x5_);
            tmvaReaderR9[iBE]->AddVariable("f8", &sieieFull5x5);
            tmvaReaderR9[iBE]->AddVariable("f9", &phoSCEtaWidth_);
        }
        else
        {
            tmvaReaderR9[iBE]->AddVariable("f0", &phoEt_);
            tmvaReaderR9[iBE]->AddVariable("f1", &phoSCEta_);
            tmvaReaderR9[iBE]->AddVariable("f2", &phoPhi_);
            tmvaReaderR9[iBE]->AddVariable("f3", &rho_);
            tmvaReaderR9[iBE]->AddVariable("f4", &sieipFull5x5);
            tmvaReaderR9[iBE]->AddVariable("f5", &s4Full5x5);
            tmvaReaderR9[iBE]->AddVariable("f6", &phoR9Full5x5_);
            tmvaReaderR9[iBE]->AddVariable("f7", &phoSCPhiWidth_);
            tmvaReaderR9[iBE]->AddVariable("f8", &sieieFull5x5);
            tmvaReaderR9[iBE]->AddVariable("f9", &phoSCEtaWidth_);
        }

        if (year == 2016)
        {
            if (iBE == 0)
            {
                tmvaReaderR9[0]->BookMVA("BDT", Form("%s/2016/weights_finalRegressor_EB_r9.xml", dir));
                fR9[0] = new TFormula("", "x[0]*0.00915591682273384+0.00042335154161760036");
            }
            else
            {
                tmvaReaderR9[1]->BookMVA("BDT", Form("%s/2016/weights_finalRegressor_EE_r9.xml", dir));
                fR9[1] = new TFormula("", "x[0]*0.010253359666186235-0.0031073467744174854");
            }
        }
        else if (year == 2017)
        {
            if (iBE == 0)
            {
                tmvaReaderR9[0]->BookMVA("BDT", Form("%s/2017/weights_finalRegressor_EB_r9.xml", dir));
                fR9[0] = new TFormula("", "x[0]*0.009489495430770628-0.002028252255487417");
            }
            else
            {
                tmvaReaderR9[1]->BookMVA("BDT", Form("%s/2017/weights_finalRegressor_EE_r9.xml", dir));
                fR9[1] = new TFormula("", "x[0]*0.014832657700482338-0.007372621685483638");
            }
        }
        else if (year == 2018)
        {
            if (iBE == 0)
            {
                tmvaReaderR9[0]->BookMVA("BDT", Form("%s/2018/weights_finalRegressor_EB_r9.xml", dir));
                fR9[0] = new TFormula("", "x[0]*0.010041591228904162-0.0034957641959179053");
            }
            else
            {
                tmvaReaderR9[1]->BookMVA("BDT", Form("%s/2018/weights_finalRegressor_EE_r9.xml", dir));
                fR9[1] = new TFormula("", "x[0]*0.014926955853444557-0.00869450644995956");
            }
        }
    }

    corrR9 = phoR9Full5x5_ + fR9[iBE]->Eval(tmvaReaderR9[iBE]->EvaluateRegression(0, "BDT"));
    return corrR9;
}