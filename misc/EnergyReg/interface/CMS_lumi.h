#include <iostream>
#include "TPad.h"
#include "TLine.h"
#include "TBox.h"
#include "TString.h"
#include "TLatex.h"

void CMS_lumi(TPad *pad, int iPeriod, int iPosX, TString lumitext_="35.9 fb^{-1}", int year_=2016, bool MakeExtraText_=true, TString extraText_="Preliminary", TString procText_="", TString extraExtraText="")
{
    TString cmsText = "CMS";
    float cmsTextFont = 61; // default is helvetic-bold

    bool writeExtraText = MakeExtraText_;
    TString extraText = extraText_;
    float extraTextFont = 52; // default is helvetica-italics

    // text sizes and text offsets with respect to the top frame
    // in unit of the top margin size
    float lumiTextSize = 0.6;
    float lumiTextOffset = 0.2;
    float cmsTextSize = 0.8;
    // float cmsTextOffset = 0.1; // only used in outOfFrame version

    float relPosX = 0.045;
    float relPosY = 0.035;
    float relExtraDY = 1.2;

    // ratio of "CMS" and extra text size
    float extraOverCmsTextSize = 0.76;

    // TString lumi_13TeV = "20.1 fb^{-1}";
    // TString lumi_8TeV = "19.7 fb^{-1}";
    // TString lumi_7TeV = "5.1 fb^{-1}";
    TString lumi_sqrtS = "";

    bool drawLogo = false;

    bool outOfFrame = false;
    if (iPosX / 10 == 0)
    {
        outOfFrame = true;
    }
    int alignY_ = 3;
    int alignX_ = 2;
    if (iPosX / 10 == 0)
        alignX_ = 1;
    if (iPosX == 0)
        alignX_ = 1;
    if (iPosX == 0)
        alignY_ = 1;
    if (iPosX / 10 == 1)
        alignX_ = 1;
    if (iPosX / 10 == 2)
        alignX_ = 2;
    if (iPosX / 10 == 3)
        alignX_ = 3;
    //if( iPosX == 0  ) relPosX = 0.12;
    int align_ = 10 * alignX_ + alignY_;

    float H = pad->GetWh();
    float W = pad->GetWw();
    float l = pad->GetLeftMargin();
    float t = pad->GetTopMargin();
    float r = pad->GetRightMargin();
    float b = pad->GetBottomMargin();
    //  float e = 0.025;

    pad->cd();

    TString yeartext;
    yeartext.Form("%d", year_);
    TString lumiText;
    if (iPeriod == 1)
    {
        // OBSOLETE
        lumiText += lumitext_;
        lumiText += " (7 TeV, ";
        lumiText += yeartext;
        lumiText += ")";
    }
    else if (iPeriod == 2)
    {
        // OBSOLETE
        lumiText += lumitext_;
        lumiText += " (8 TeV, ";
        lumiText += yeartext;
        lumiText += ")";
    }
    else if (iPeriod == 3)
    {
        // OBSOLETE
        lumiText = lumitext_;
        lumiText += " (8 TeV)";
        lumiText += " + ";
        lumiText += lumitext_;
        lumiText += " (7 TeV)";
    }
    else if (iPeriod == 4)
    {
        lumiText += lumitext_;
        lumiText += " (13 TeV, ";
        lumiText += yeartext;
        lumiText += ")";
    }
    else if (iPeriod == 5)
    {
        lumiText += lumitext_;
        lumiText += " (13 TeV)";
    }
    else if (iPeriod == 7)
    {
        if (outOfFrame)
            lumiText += "#scale[0.85]{";
        lumiText += lumitext_;
        lumiText += " (13 TeV)";
        lumiText += " + ";
        lumiText += lumitext_;
        lumiText += " (8 TeV)";
        lumiText += " + ";
        lumiText += lumitext_;
        lumiText += " (7 TeV)";
        if (outOfFrame)
            lumiText += "}";
    }
    else if (iPeriod == 12)
    {
        lumiText += "8 TeV";
    }
    else if (iPeriod == 0)
    {
        lumiText += lumi_sqrtS;
    }

    // std::cout << lumiText << std::endl;

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(kBlack);

    float extraTextSize = extraOverCmsTextSize * cmsTextSize * 0.9;

    latex.SetTextFont(42);
    latex.SetTextAlign(31);
    latex.SetTextSize(0.9 * lumiTextSize * t);
    latex.DrawLatex(1 - r, 1 - t + lumiTextOffset * t, lumiText);

    // std::cout << procText_ << std::endl;

    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.9 * lumiTextSize * t);
    latex.DrawLatex(l, 1 - t + lumiTextOffset * t, procText_);

    if (outOfFrame)
    {
        latex.SetTextFont(cmsTextFont);
        latex.SetTextAlign(11);
        latex.SetTextSize(0.9 * cmsTextSize * t);
        latex.DrawLatex(l, 1 - t + lumiTextOffset * t, cmsText);
    }

    pad->cd();

    float posX_ = 0;
    if (iPosX % 10 <= 1)
    {
        posX_ = l + relPosX * (1 - l - r);
    }
    else if (iPosX % 10 == 2)
    {
        posX_ = l + 0.5 * (1 - l - r);
    }
    else if (iPosX % 10 == 3)
    {
        posX_ = 1 - r - relPosX * (1 - l - r);
    }
    float posY_ = 1 - t - relPosY * (1 - t - b);
    if (!outOfFrame)
    {
        if (drawLogo)
        {
            posX_ = l + 0.045 * (1 - l - r) * W / H;
            posY_ = 1 - t - 0.045 * (1 - t - b);
            float xl_0 = posX_;
            float yl_0 = posY_ - 0.15;
            float xl_1 = posX_ + 0.15 * H / W;
            float yl_1 = posY_;
            //TASImage* CMS_logo = new TASImage("CMS-BW-label.png");
            TPad *pad_logo = new TPad("logo", "logo", xl_0, yl_0, xl_1, yl_1);
            pad_logo->Draw();
            pad_logo->cd();
            //CMS_logo->Draw("X");
            pad_logo->Modified();
            pad->cd();
            delete pad_logo;
        }
        else
        {
            latex.SetTextFont(cmsTextFont);
            latex.SetTextSize(0.9 * cmsTextSize * t);
            latex.SetTextAlign(align_);
            latex.DrawLatex(posX_, posY_, cmsText);
            if (writeExtraText)
            {
                latex.SetTextFont(extraTextFont);
                latex.SetTextAlign(align_);
                latex.SetTextSize(0.9 * extraTextSize * t);
                latex.DrawLatex(posX_, posY_ - relExtraDY * cmsTextSize * t, extraText + " " + extraExtraText);
            }
        }
    }
    else if (writeExtraText)
    {
        if (iPosX == 0)
        {
            posX_ = l + relPosX * (1 - l - r);
            posY_ = 1 - t + lumiTextOffset * t;
        }
        latex.SetTextFont(extraTextFont);
        latex.SetTextSize(0.9 * extraTextSize * t);
        latex.SetTextAlign(align_);
        latex.DrawLatex(posX_ + 0.08, posY_, extraText + " " + extraExtraText);
    }
}