import ROOT
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator, AutoMinorLocator
from sigmaEff import sigmaEff


def pltSty(xName = "x-axis", yName = "y-axis", TitleSize = 15, LabelSize = 15, TickSize = 13, MajTickLength = 7, MinTickLength = 4, yAuto = True):
    ax = plt.gca()
    ax.set_xlabel(xName, fontsize = LabelSize, loc = "right")
    ax.set_ylabel(yName, fontsize = LabelSize, loc = "top")
    ax.text(1, 1, "(13 TeV)", horizontalalignment = "right", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize)
    ax.text(0, 1, "CMS", horizontalalignment = "left", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize * 1.3, fontweight = "bold")
    ax.text(TitleSize * 0.009 + 0.02, 1, "Work-in-progress", horizontalalignment = "left", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize)

    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if (yAuto):
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction = "in", length = MajTickLength, labelsize = TickSize, top = True, right = True)
    ax.tick_params(direction = "in", length = MinTickLength, which = "minor", labelsize = TickSize, top = True, right = True)


# Reference: https://stackoverflow.com/a/29643643
def hex2rgb(h, a):
    h = h.lstrip('#')
    rgb = tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))
    return (rgb[0], rgb[1], rgb[2], a)

ROOT.EnableImplicitMT(10)
cat = ["Resolved", "Merged-2Gsf", "Merged-1Gsf"]
for i in range(len(cat)):
    df_UL = ROOT.RDataFrame("miniTree", "./minitree/UL2017/minitree_HDalitz_ggF_eeg_m125_UL2017.root").Filter("category == {}".format(i+1))
    df_Le = ROOT.RDataFrame("miniTree", "./minitree/test/minitree_HDalitz_ggF_eeg_m125_fall17.root").Filter("category == {}".format(i+1))

    mass_UL = df_UL.AsNumpy(columns=["higgsMass"])["higgsMass"]
    xmin, xmax, sigma_UL = sigmaEff(mass_UL)
    print("effective sigma for UL2017 = {}".format(sigma_UL))
    mass_Le = df_Le.AsNumpy(columns=["higgsMass"])["higgsMass"]
    xmin, xmax, sigma_Le = sigmaEff(mass_Le)
    print("effective sigma for Fall17 = {}".format(sigma_Le))

    binContent, binEdget, patches = plt.hist(
        mass_UL,
        range=(110, 140),
        bins=60,
        ls="-", lw=2, fc=hex2rgb("#3A9679", 0), ec="#3A9679",
        histtype="stepfilled",
        density=True,
        label="UL2017 ($\sigma_{eff}$ = %.2f GeV)" %(sigma_UL)
    )

    binContent, binEdget, patches = plt.hist(
        mass_Le,
        range=(110, 140),
        bins=60,
        ls="-", lw=2, fc=hex2rgb("#E16262", 0), ec="#E16262",
        histtype="stepfilled",
        density=True,
        label="Fall17 ($\sigma_{eff}$ = %.2f GeV)" %(sigma_Le)
    )

    ax = plt.gca()
    ax.set_ylim(0, max(binContent) * 1.3)
    plt.text(0.05, 0.9, "$H \\rightarrow \gamma^*\gamma \\rightarrow ee\gamma, ggF$", fontsize=13, transform=ax.transAxes)
    plt.text(0.05, 0.8, cat[i], fontsize=13, transform=ax.transAxes)
    pltSty(xName="$\mathrm{M_{\mu\mu\gamma}}$ GeV", yName="A.U.")

    dir = "./plots/Legacy"
    os.makedirs("./plots/Legacy", exist_ok=True)
    plt.legend(loc="best", edgecolor="none", fontsize=13)
    plt.tight_layout()
    plt.savefig("{}/{}.pdf".format(dir, cat[i]))
    plt.close("all")