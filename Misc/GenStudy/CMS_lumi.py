import ROOT as rt

# void CMS_lumi(TPad *pad, int iPeriod, int iPosX, TString lumitext_="35.9 fb^{-1}", int year_=2016, bool MakeExtraText_=true, TString extraText_="Preliminary", TString procText_="H#rightarrow#gamma*#gamma#rightarrow#mu#mu#gamma", TString extraExtraText="")


def CMS_lumi(pad,  iPeriod,  iPosX, lumitext_, year_, MakeExtraText_, extraText_, procText_, extraExtraText):
    cmsText = "CMS"
    cmsTextFont = 61

    # writeExtraText = True
    writeExtraText = MakeExtraText_
    extraText = extraText_
    extraTextFont = 52

    # text sizes and text offsets with respect to the top frame
    # in unit of the top margin size
    lumiTextSize = 0.65
    lumiTextOffset = 0.2

    cmsTextSize = 0.8
    cmsTextOffset = 0.08 # only used in outOfFrame version

    relPosX = 0.045
    relPosY = 0.035
    relExtraDY = 1.2

    # ratio of "CMS" and extra text size
    extraOverCmsTextSize = 0.76

    # lumi_13TeV = "35.9 fb^{-1}"
    # lumi_8TeV = "19.7 fb^{-1}"
    # lumi_7TeV = "5.1 fb^{-1}"
    lumi_sqrtS = ""

    drawLogo = False

    outOfFrame = False
    if(iPosX/10 == 0):
        outOfFrame = True
        
    print ('CMS (Preliminary) out of frame? -->', outOfFrame)

    alignY_ = 3
    alignX_ = 2
    if(iPosX/10 == 0):
        alignX_ = 1
    if(iPosX == 0):
        alignY_ = 1
    if(iPosX/10 == 1):
        alignX_ = 1
    if(iPosX/10 == 2):
        alignX_ = 2
    if(iPosX/10 == 3):
        alignX_ = 3
    align_ = 10*alignX_ + alignY_

    H = pad.GetWh()
    W = pad.GetWw()
    l = pad.GetLeftMargin()
    t = pad.GetTopMargin()
    r = pad.GetRightMargin()
    b = pad.GetBottomMargin()
    # e = 0.025

    pad.cd()

    yeartext = '{}'.format(year_)
    lumiText = ""
    if(iPeriod == 1):
        # OBSOLETE
        lumiText += lumitext_
        lumiText += " (7 TeV, "
        lumiText += yeartext
        lumiText += ")"
    elif (iPeriod == 2):
        # OBSOLETE
        lumiText += lumitext_
        lumiText += " (8 TeV, "
        lumiText += yeartext
        lumiText += ")"
    elif(iPeriod == 3):
        # OBSOLETE
        lumiText = lumitext_
        lumiText += " (8 TeV)"
        lumiText += " + "
        lumiText += lumitext_
        lumiText += " (7 TeV)"
    elif (iPeriod == 4):
        lumiText += lumitext_
        lumiText += " (13 TeV, "
        lumiText += yeartext
        lumiText += ")"
    elif (iPeriod == 5):
        lumiText += lumitext_
        lumiText += " (13 TeV)"
    elif (iPeriod == 7):
        if outOfFrame:
            lumiText += "#scale[0.85]{"
            lumiText += lumitext_
            lumiText += " (13 TeV)"
            lumiText += " + "
            lumiText += lumitext_
            lumiText += " (8 TeV)"
            lumiText += " + "
            lumiText += lumitext_
            lumiText += " (7 TeV)"
        if outOfFrame:
            lumiText += "}"
    elif (iPeriod == 12):
        lumiText += "8 TeV"
    elif (iPeriod == 0):
        lumiText += lumi_sqrtS

    print('lumiText = ', lumiText)

    latex = rt.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(rt.kBlack)

    extraTextSize = extraOverCmsTextSize*cmsTextSize

    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(0.9*lumiTextSize*t)
    latex.DrawLatex(1-r, 1-t+lumiTextOffset*t, lumiText)

    print('procText = ', procText_)

    if not outOfFrame:
        latex.SetTextFont(42)
        latex.SetTextAlign(11)
        latex.SetTextSize(0.9 * lumiTextSize * t)
        latex.DrawLatex(l, 1 - t + lumiTextOffset * t, procText_)

    if outOfFrame:
        latex.SetTextFont(cmsTextFont)
        latex.SetTextAlign(11)
        latex.SetTextSize(0.9*cmsTextSize*t)
        latex.DrawLatex(l, 1-t+lumiTextOffset*t, cmsText)

    pad.cd()

    posX_ = 0
    if(iPosX % 10 <= 1):
        posX_ = l + relPosX*(1-l-r)
    elif(iPosX % 10 == 2):
        posX_ = l + 0.5*(1-l-r)
    elif(iPosX % 10 == 3):
        posX_ = 1-r - relPosX*(1-l-r)

    posY_ = 1-t - relPosY*(1-t-b)

    if not outOfFrame:
        if drawLogo:
            posX_ = l + 0.045*(1-l-r)*W/H
            posY_ = 1-t - 0.045*(1-t-b)
            xl_0 = posX_
            yl_0 = posY_ - 0.15
            xl_1 = posX_ + 0.15*H/W
            yl_1 = posY_
            # CMS_logo = rt.TASImage("CMS-BW-label.png")
            pad_logo = rt.TPad("logo", "logo", xl_0, yl_0, xl_1, yl_1)
            pad_logo.Draw()
            pad_logo.cd()
            # CMS_logo.Draw("X")
            pad_logo.Modified()
            pad.cd()
        else:
            latex.SetTextFont(cmsTextFont)
            latex.SetTextSize(0.9 * cmsTextSize*t)
            latex.SetTextAlign(align_)
            latex.DrawLatex(posX_, posY_, cmsText)
            if writeExtraText:
                latex.SetTextFont(extraTextFont)
                latex.SetTextAlign(align_)
                latex.SetTextSize(0.9*extraTextSize*t)
                latex.DrawLatex(posX_, posY_ - relExtraDY *
                                cmsTextSize*t, extraText+" "+extraExtraText)
    elif writeExtraText:
        if iPosX == 0:
            posX_ = l + relPosX*(1-l-r)
            posY_ = 1-t+lumiTextOffset*t

        latex.SetTextFont(extraTextFont)
        latex.SetTextSize(0.9*extraTextSize*t)
        latex.SetTextAlign(align_)
        latex.DrawLatex(posX_+cmsTextOffset, posY_, extraText+" "+extraExtraText)

    pad.Update()
