#ifndef CMS_LUMI_MODIFICATION_H_
#define CMS_LUMI_MODIFICATION_H_

#include <iostream>
#include "TPad.h"
#include "TString.h"

void CMS_lumi(TPad *pad, int iPeriod, int iPosX, TString lumitext_, int year_, bool MakeExtraText_, TString extraText_, TString procText_, TString extraExtraText);

#endif