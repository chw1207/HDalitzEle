# Plotting Tools

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import ROOT
from Tools.CMS_lumi import CMS_lumi
from matplotlib.ticker import AutoLocator, AutoMinorLocator
import os
import sys
import pickle
import seaborn as sns  
# import tdrstyle
# ROOT.gROOT.LoadMacro("tdrstyle.C")
# ROOT.gROOT.ProcessLine("setTDRStyle()")

# ROOT.gROOT.SetBatch()
# ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")
# tdrstyle.setTDRStyle()

class color:
   PURPLE = "\033[95m"
   CYAN = "\033[96m"
   DARKCYAN = "\033[36m"
   BLUE = "\033[94m"
   GREEN = "\033[92m"
   YELLOW = "\033[93m"
   RED = "\033[91m"
   BOLD = "\033[1m"
   UNDERLINE = "\033[4m"
   END = "\033[0m"

TickSize = 0.02
AxisTitleSize = 0.04
AxisLabelSize = 0.04
   
def prGreen(prt): print("\033[92m {}\033[00m" .format(prt))


def plot_mva(df, column, bins, logscale = False, ax = None, title = None, ls = 'dashed', alpha = 0.3, sample='', cat = "Matchlabel", Wt = "Wt"):
    hcolor = ['#e84545', '#0061a8']
    histtype = "bar"
    if sample == 'test':
        histtype = "step"
    if ax is None:
        ax = plt.gca()
    for name, group in df.groupby(cat):
        if name == 0:
            label = "background"
        else:
            label = "signal"
        group[column].hist(
            bins = bins, histtype = histtype, alpha = alpha, grid = False,
            label = label + ' ' + sample, ax = ax, density = False, ls = ls, 
            color = hcolor[name], weights = group[Wt]/group[Wt].sum(), linewidth = 2
        )
    ax.set_ylabel("Normalized entries", fontsize = 15, loc = 'top')
    ax.set_xlabel(column, fontsize = 15, loc = 'right')
    # ax.set_title("CMS", fontsize = 15, loc = 'left', fontweight = 'bold')
    ax.text(0, 1, "CMS", horizontalalignment = 'left', verticalalignment = 'bottom', transform=ax.transAxes, fontsize = 15, fontweight = 'bold')
    ax.text(0.13, 1, "$\it{Simulation}$", horizontalalignment = 'left', verticalalignment = 'bottom', transform = ax.transAxes, fontsize = 13)
    ax.text(1, 1, "41.7 $fb^{-1} (13TeV,\ 2017)$", horizontalalignment = 'right', verticalalignment = 'bottom', transform = ax.transAxes, fontsize = 13)
    if logscale:
        ax.set_yscale("log", nonpositive = 'clip')
    else:
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.legend(loc = 'upper center', edgecolor = 'none', fontsize = 13)
    ax.tick_params( # major tick 
        direction = 'in', length = 8,
        labelsize = 13, top = True, right = True
    )
    ax.tick_params( # minor tick 
        direction = 'in', length = 4, which = 'minor',
        labelsize = 13, top = True, right = True
    )


def plot_roc_curve(df, score_column, tpr_threshold=0, ax=None, color=None, linestyle='-', label=None, cat="Matchlabel", Wt="Wt"):
    from sklearn import metrics
    if ax is None:
        ax = plt.gca()
    if label is None:
        label = score_column
    
    # * tpr: true positive rate = signal efficiency 
    # * fpr: false positive rate = background efficiency
    fpr, tpr, thresholds = metrics.roc_curve(
        df[cat], df[score_column], sample_weight=df[Wt])
    mask = tpr > tpr_threshold
    fpr_masked, tpr_masked = fpr[mask], tpr[mask]
    bkgrej_mask = [1 - fpr for fpr in fpr_masked]
    auc = metrics.auc(fpr_masked, tpr_masked)
    label = label+' auc = '+str(round(auc*100, 1))+'%'
    ax.plot(tpr_masked, bkgrej_mask, label = label, color = color,
            linestyle = linestyle, linewidth = 3, alpha = 1.0)
    # ax.set_yscale("log")
    ax.legend(loc = 'best', edgecolor= 'none', fontsize = 13)
    # ax.tick_params(direction='in', which='both', bottom=True, top=True, left=True, right=True)
    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_major_locator(AutoLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params( # major tick 
    direction = 'in', length = 8,
    labelsize = 13, top = True, right = True
    )
    ax.tick_params( # minor tick 
        direction = 'in', length = 4, which = 'minor',
        labelsize = 13, top = True, right = True
    )
    return [auc, tpr, bkgrej_mask, thresholds]


def plot_single_roc_point(df, var='Fall17isoV1wpLoose', ax=None, marker='o', markersize=6, color="red", label='', cat="Matchlabel", Wt="Wt"):
    backgroundpass = df.loc[(df[var] == 1) & (df[cat] == 0), Wt].sum()
    backgroundrej = df.loc[(df[var] == 0) & (df[cat] == 0), Wt].sum()
    signalpass = df.loc[(df[var] == 1) & (df[cat] == 1), Wt].sum()
    signalrej = df.loc[(df[var] == 0) & (df[cat] == 1), Wt].sum()
    backgroundrej = backgroundrej/(backgroundpass+backgroundrej)
    signaleff = signalpass/(signalpass+signalrej)
    ax.plot([signaleff], [1-backgroundrej], marker=marker,
            color=color, markersize=markersize, label=label)
    ax.set_yscale("log")
    ax.legend(loc='best')
    ax.tick_params(direction='in', which='both', bottom=True,
                   top=True, left=True, right=True)


def pngtopdf(ListPattern=[], Save="mydoc.pdf"):
    import glob
    import PIL.Image
    L = []
    for List in ListPattern:
        L += [PIL.Image.open(f) for f in glob.glob(List)]
    for i, Li in enumerate(L):
        rgb = PIL.Image.new('RGB', Li.size, (255, 255, 255))
        rgb.paste(Li, mask=Li.split()[3])
        L[i] = rgb
    L[0].save(Save, "PDF", resolution=100.0,
              save_all=True, append_images=L[1:])

def MakeFeaturePlotsComb_sep(df_final, features, feature_bins, leppos, MVA="XGB_1", OutputDirName='Output', cat="EleType", label=["Background", "Signal"], weight="NewWt"):
    prGreen("Making Combined"+" dataset feature plots")
    hcolor = ['#0f4c75', '#903749']
    for m in range(len(features)):
        print('Feature: {}'.format(features[m-1]))
        fig, ax = plt.subplots(1, 1, figsize=(5, 4))
        for i, group_df in df_final[df_final['Dataset'] == "Train"].groupby(cat):
            group_df[features[m-1]].hist(histtype='stepfilled', bins=feature_bins[m-1], alpha=0.4, label=label[i]+"_Train",
                                         ax=ax, density=False, ls='-', color=hcolor[i], weights=group_df[weight]/group_df[weight].sum(), linewidth=4)
        for i, group_df in df_final[df_final['Dataset'] == "Test"].groupby(cat):
            group_df[features[m-1]].hist(histtype='step', bins=feature_bins[m-1], alpha=0.9, label=label[i]+"_Test",
                                         ax=ax, density=False, color=hcolor[i], weights=group_df[weight]/group_df[weight].sum(), linewidth=1.5)
            # df_new = pd.concat([group_df, df_new],ignore_index=True, sort=False)
        ax.legend(loc=leppos[m-1])
        ax.set_xlabel(features[m-1])
        ax.set_ylabel('Normalized entries')
        ax.set_yscale("log")
        ax.set_title(features[m-1])
        ax.set_xlim([feature_bins[m-1][0], feature_bins[m-1][-1]])
        ax.tick_params(direction='in', which='both', bottom=True,
                       top=True, left=True, right=True)
        plt.savefig(OutputDirName+"/"+MVA+"/"+MVA+"_featureplots_comb_" +
                    features[m-1]+".pdf", bbox_inches = 'tight')
        plt.close('all')

def DrawFeatureHist(df, feature, bininfo, XaxisName, yaxisunit, histcolor = ['#0061a8', '#e84545'], logy = True, y_axisscale = 1.5, drawxoverflow = True, drawxunderflow = False, outname = 'testfeatureplot', wei='' ):
    plt.figure(figsize=(6, 6))
    binContent_s, binEdge_s, patches_s = plt.hist(
        df[df['Type'] == "Signal"][feature],
        range = (bininfo[1], bininfo[2]) , bins = bininfo[0],
        ls = "-", lw = 2, fc = (0/255, 97/255, 168/255, 0.5), ec = histcolor[0],
        histtype = 'stepfilled',
        weights = df[df['Type'] == "Signal"][wei],
        density = True,
        label = "Signal"
    )
    binContent_b, binEdge_b, patches_b = plt.hist(
        df[df['Type'] == "Background"][feature],
        range = (bininfo[1], bininfo[2]) , bins = bininfo[0],
        ls = "-", lw = 2, fc = (232/255, 69/255, 69/255, 0.5), ec = histcolor[1],
        histtype = 'stepfilled',
        weights = df[df['Type'] == "Background"][wei],
        density = True,
        label = "Background"
    )
    
    binwidth = binEdge_s[1] - binEdge_s[0]
    ymaxval = 0.
    if np.max(binContent_b) > np.max(binContent_s):
        ymaxval = np.max(binContent_b)
    else:
        ymaxval = np.max(binContent_s)
    
    yminval = 1E10
    res = binContent_b[binContent_b != 0]
    if np.min(res) < yminval:
        yminval = np.min(res)

    ax = plt.gca()
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(yminval * 0.9, ymaxval * y_axisscale)
    else:
        ax.set_ylim(0, ymaxval * y_axisscale)
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(bininfo[1], bininfo[2])
    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.set_ylabel("A.U.", fontsize = 15, loc = 'top')
    ax.set_xlabel(XaxisName, fontsize = 15, loc = 'right')
    
    # lumi information text
    ax.text(0, 1, "CMS", horizontalalignment = 'left', verticalalignment = 'bottom', transform=ax.transAxes, fontsize = 15, fontweight = 'bold')
    ax.text(0.13, 1, "$\it{Simulation}$", horizontalalignment = 'left', verticalalignment = 'bottom', transform = ax.transAxes, fontsize = 13)
    ax.text(1, 1, "41.7 $fb^{-1} (13TeV,\ 2017)$", horizontalalignment = 'right', verticalalignment = 'bottom', transform = ax.transAxes, fontsize = 13)
    
    ax.legend(loc = "best", edgecolor= "none", fontsize = 13)
    ax.tick_params( # major tick 
        direction = 'in', length = 8,
        labelsize = 13, top = True, right = True
    )
    ax.tick_params( # minor tick 
        direction = 'in', length = 4, which = 'minor',
        labelsize = 13, top = True, right = True
    )

    ax.legend(loc = "best", edgecolor= "none", fontsize = 13)
    plt.tight_layout()
    plt.draw()
    plt.savefig("{}.pdf".format(outname))
    print("Save fig in {}.pdf".format(outname))
    plt.close()

# def fillhist(hist, val, wei):
#     # val is a pandas Series object
#     for i in val:
#         hist.Fill(i)

# def DrawFeature(df, feature, bininfo, XaxisName, yaxisunit, histcolor = ['#0061a8', '#e84545'], logy = True, y_axisscale = 1.5, drawxoverflow = True, drawxunderflow = False, outname = 'testfeatureplot'):
#     ROOT.gStyle.SetOptStat(0)
#     ROOT.gStyle.SetPadTickX(1)
#     ROOT.gStyle.SetPadTickY(1)
    
#     hist_sig = ROOT.TH1F('hist_sig_{}'.format(feature),'',bininfo[0],bininfo[1],bininfo[2])
#     hist_bkg = ROOT.TH1F('hist_bkg_{}'.format(feature),'',bininfo[0],bininfo[1],bininfo[2])
#     fillhist(hist_sig, df[df['Type'] == "Signal"][feature], df[df['Type'] == "Signal"]["mcwei"])
#     fillhist(hist_bkg, df[df['Type'] == "Background"][feature], df[df['Type'] == "Background"]["mcwei"])
#     listhist = [hist_sig, hist_bkg]
    
#     ymaxval = 0.
#     yminval = 1E10
#     for i, hist in enumerate(listhist):
#         hist.StatOverflows()
#         hist.SetLineWidth(3)
#         hist.SetLineColor(ROOT.TColor.GetColor(histcolor[i]))
#         hist.SetFillColorAlpha(ROOT.TColor.GetColor(histcolor[i]), 0.35)
#         hist.Scale(1./hist.Integral(-1,-1))
#         if hist.GetMaximum() > ymaxval:
#             ymaxval = hist.GetMaximum()

#         if hist.GetMinimum(0) < yminval:
#             yminval = hist.GetMinimum(0)
    
#     binwidth = listhist[0].GetXaxis().GetBinWidth(1)
    
#     yaxtitletext = ''
#     if yaxisunit == '':
#         yaxtitletext = 'Normalized events / {:g}'.format(binwidth)
#     else:
#         yaxtitletext = 'Normalized events / {:g} {}'.format(binwidth, yaxisunit)

#     c = ROOT.TCanvas('c', '', 800, 800)
#     c.cd( )
#     ROOT.gPad.SetRightMargin(0.05)
#     ROOT.gPad.SetTopMargin(0.07)
#     ROOT.gPad.SetLeftMargin(0.14)
#     ROOT.gPad.SetBottomMargin(0.15)
#     if logy:
#         c.SetLogy()
#         listhist[0].GetYaxis().SetRangeUser(yminval*0.90, ymaxval*y_axisscale)
#         # listhist[0].GetYaxis().SetMoreLogLabels()
#     else:
#         listhist[0].GetYaxis().SetRangeUser(0, ymaxval*y_axisscale)
    
#     for i, hist in enumerate(listhist):
#         if i==0:
#             if drawxoverflow:
#                 listhist[i].GetXaxis().SetRange(1, listhist[i].GetNbinsX() + 1)
#             if drawxunderflow:
#                 listhist[i].GetXaxis().SetRange(0, listhist[i].GetNbinsX())
#             if drawxoverflow and drawxunderflow:
#                 listhist[i].GetXaxis().SetRange(0, listhist[i].GetNbinsX() + 1)
#             listhist[i].GetXaxis().SetTitle(XaxisName)
#             listhist[i].GetYaxis().SetTitle(yaxtitletext)
#             listhist[i].GetXaxis().SetTickSize(TickSize)
#             listhist[i].GetXaxis().SetTitleSize(AxisTitleSize)
#             listhist[i].GetXaxis().SetLabelSize(AxisLabelSize)
#             listhist[i].GetYaxis().SetTickSize(TickSize)
#             listhist[i].GetYaxis().SetTitleSize(AxisTitleSize)
#             listhist[i].GetYaxis().SetLabelSize(AxisLabelSize)
#             listhist[i].GetXaxis().SetTitleOffset(1.3)
#             listhist[i].GetYaxis().SetTitleOffset(1.5)
#             listhist[i].Draw('hist')
#         else:
#             if drawxoverflow:
#                 listhist[i].GetXaxis().SetRange(1, listhist[i].GetNbinsX() + 1)
#             if drawxunderflow:
#                 listhist[i].GetXaxis().SetRange(0, listhist[i].GetNbinsX())
#             if drawxoverflow and drawxunderflow:
#                 listhist[i].GetXaxis().SetRange(0, listhist[i].GetNbinsX() + 1)
                
#             listhist[i].Draw('histsame')
            
#     legpos_y1 = 0.87
#     legpos_dy = 0.06
#     listleg = ['Signal', 'Background']
#     legsty = 'f'
#     l = ROOT.TLegend(0.67, legpos_y1 - legpos_dy * len(listhist), 0.85, legpos_y1)
#     l.SetTextSize(0.04)
#     for ihist, hist in enumerate(listhist):
#         l.AddEntry(hist, listleg[ihist], legsty)
#     l.SetFillColor(0)  # Set the background to be white
#     l.SetLineColor(0)
#     l.SetFillStyle(0)
#     l.Draw("same")
            
#     CMS_lumi(c, 4, 0, "41.5 fb^{-1}", 2017, True, "Simulation", "", "")

#     c.SaveAs('{}.pdf'.format(outname))


def plot_featureimp(feature_importance, feature_names, MVA='XGB_1', OutputDirName='Output'):
    df_featimp = pd.DataFrame(
        feature_importance, index=feature_names, columns=['importance'])
    df_featimp.sort_values('importance', inplace=True, ascending=True)

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    df_featimp['importance'].plot(kind='barh', ax=ax)
    ax.set_xlabel('Feature importance', fontsize=13)
    ax.set_ylabel('Feature', fontsize=13)
    ax.tick_params(direction='in', which='both', bottom=True,
                   top=False, left=True, right=True)
    ax.set_xlim(0, df_featimp['importance'].iloc[-1]*1.15)
    for i in ax.patches:
        # get_width pulls left or right; get_y pushes up or down
        ax.text(i.get_width()+0.002, i.get_y(),
                '{:.3g}'.format(i.get_width()), fontsize=9)

    plt.savefig(OutputDirName+'/'+MVA+'/'+MVA +
                '_featureimportance.pdf', bbox_inches='tight')
    plt.close('all')


def plot_FeatureImpSel(df, MVA='XGB', OutputDirName='Output'):
    titlesize = 15
    labelsize = 13
    secondcol = 'green'
    fig, ax = plt.subplots(1, 1, figsize=(12, 5))
    df['importance'].plot(kind='bar',ax=ax)

    ax.set_ylim([0, df['importance'].iloc[-1]*1.2])

    ax2 = ax.twinx()
    df['accuracy'].plot(ax=ax2, lw=3, color=secondcol, rot=90)
    for i in ax.patches:
        textxpos = i.get_x()+0.4*i.get_width()
        ax.text(textxpos, i.get_height()+0.01,'{:.3g}'.format(i.get_height()), rotation=90, ha = 'center', fontsize=11)
    
    ax.tick_params(direction='in', which='both', bottom=True, top=False, left=True, right=False)
    ax.tick_params(axis='x',labelsize=labelsize)
    ax.tick_params(axis='y',labelsize=labelsize)
    ax.set_ylabel('Importance',fontsize=titlesize)
    ax.set_xlabel('Feature',fontsize=titlesize)

    ax2.tick_params(direction='in', which='both', bottom=False, top=False, left=False, right=True)
    ax2.set_ylabel('Accuracy (%)',fontsize=titlesize) #  including the most up to the nth most important features
    ax2.spines["right"].set_edgecolor(secondcol)
    ax2.spines["top"].set_edgecolor(secondcol)
    ax2.tick_params(axis='y', colors=secondcol)
    ax2.tick_params(axis='y',labelsize=labelsize)
    ax2.yaxis.label.set_color(secondcol)

    ax3 = ax.twiny()
    ax3.set_xticks(df['number_of_feature'].tolist())
    ax3.set_xticklabels(df['number_of_feature'].tolist())
    ax3.set_xlim(ax.get_xlim()[0]+1, ax.get_xlim()[1]+1)
    ax3.spines["right"].set_edgecolor(secondcol)
    ax3.spines["top"].set_edgecolor(secondcol)
    ax3.set_xlabel('Number of features included',fontsize=titlesize)
    ax3.xaxis.label.set_color(secondcol)
    ax3.tick_params(axis='x', colors=secondcol)
    ax3.tick_params(direction='in', which='both', bottom=False, top=True, left=False, right=True)

    ax.invert_xaxis()

    plt.savefig(OutputDirName+'/'+MVA+'/'+MVA+'_featimp_accu.pdf', bbox_inches='tight')
    print("Save the fig in: {}".format(OutputDirName+'/'+MVA+'/'+MVA+'_featimp_accu.pdf'))
    plt.close('all')

def df2table(df, MVA='XGB', OutputDirName='Output'):
    # plot the table
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax = plt.gca() #get current axis
    ax.axis('off')
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    rcolors = np.full(len(df.index), '#bedcfa')
    ccolors = np.full(len(df.columns), '#bedcfa')
    Table = pd.plotting.table(ax, df, rowLoc='left', loc = 'center', rowColours = rcolors, colColours = ccolors)
    Table.auto_set_font_size(True)
    Table.set_fontsize(13)
    plt.tight_layout()
    plt.savefig(OutputDirName+'/'+MVA+'/'+MVA+'_featsel_dataframe.pdf', bbox_inches='tight')
    
# Reference: 
# [1]: https://www.kaggle.com/drazen/heatmap-with-sized-markers
def heatmap(x, y, **kwargs):
    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = [1]*len(x)

    if 'palette' in kwargs:
        palette = kwargs['palette']
        n_colors = len(palette)
    else:
        n_colors = 256 # Use 256 colors for the diverging color palette
        palette = sns.color_palette("Blues", n_colors) 

    if 'color_range' in kwargs:
        color_min, color_max = kwargs['color_range']
    else:
        color_min, color_max = min(color), max(color) # Range of values that will be mapped to the palette, i.e. min and max possible correlation

    def value_to_color(val):
        if color_min == color_max:
            return palette[-1]
        else:
            val_position = float((val - color_min)) / (color_max - color_min) # position of value in the input range, relative to the length of the input range
            val_position = min(max(val_position, 0), 1) # bound the position betwen 0 and 1
            ind = int(val_position * (n_colors - 1)) # target index in the color palette
            return palette[ind]

    if 'size' in kwargs:
        size = kwargs['size']
    else:
        size = [1]*len(x)

    if 'size_range' in kwargs:
        size_min, size_max = kwargs['size_range'][0], kwargs['size_range'][1]
    else:
        size_min, size_max = min(size), max(size)

    size_scale = kwargs.get('size_scale', 250)

    def value_to_size(val):
        if size_min == size_max:
            return 1 * size_scale
        else:
            val_position = (val - size_min) * 0.99 / (size_max - size_min) + 0.01 # position of value in the input range, relative to the length of the input range
            val_position = min(max(val_position, 0), 1) # bound the position betwen 0 and 1
            return val_position * size_scale
    if 'x_order' in kwargs: 
        x_names = [t for t in kwargs['x_order']]
    else:
        x_names = [t for t in sorted(set([v for v in x]))]
    x_to_num = {p[1]:p[0] for p in enumerate(x_names)}

    if 'y_order' in kwargs: 
        y_names = [t for t in kwargs['y_order']]
    else:
        y_names = [t for t in sorted(set([v for v in y]))]
    y_to_num = {p[1]:p[0] for p in enumerate(y_names)}

    plot_grid = plt.GridSpec(1, 15, hspace=0.2, wspace=0.1) # Setup a 1x10 grid
    ax = plt.subplot(plot_grid[:,:-1]) # Use the left 14/15ths of the grid for the main plot

    marker = kwargs.get('marker', 's')

    kwargs_pass_on = {k:v for k,v in kwargs.items() if k not in [
         'color', 'palette', 'color_range', 'size', 'size_range', 'size_scale', 'marker', 'x_order', 'y_order'
    ]}

    ax.scatter(
        x=[x_to_num[v] for v in x],
        y=[y_to_num[v] for v in y],
        marker=marker,
        s=[value_to_size(v) for v in size], 
        c=[value_to_color(v) for v in color],
        **kwargs_pass_on
    )
    ax.set_xticks([v for k,v in x_to_num.items()])
    ax.set_xticklabels([k for k in x_to_num], rotation=90)
    ax.set_yticks([v for k,v in y_to_num.items()])
    ax.set_yticklabels([k for k in y_to_num])

    ax.grid(False, 'major')
    ax.grid(True, 'minor')
    ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
    ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)

    ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
    ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    ax.set_facecolor('#F1F1F1')
    
    # print ([x_to_num[v] for v in x])
    xpos = [x_to_num[v] for v in x]
    ypos = [y_to_num[v] for v in y]
    for i, cval in enumerate(color.tolist()):
        ax.text(xpos[i], ypos[i], '{0:.2f}'.format(cval), horizontalalignment='center', verticalalignment='center', size=4, color='black')

    # Add color legend on the right side of the plot
    if color_min < color_max:
        ax = plt.subplot(plot_grid[:,-1]) # Use the rightmost column of the plot

        col_x = [0]*len(palette) # Fixed x coordinate for the bars
        bar_y=np.linspace(color_min, color_max, n_colors) # y coordinates for each of the n_colors bars

        bar_height = bar_y[1] - bar_y[0]
        ax.barh(
            y=bar_y,
            width=[5]*len(palette), # Make bars 5 units wide
            left=col_x, # Make bars start at 0
            height=bar_height,
            color=palette,
            linewidth=0
        )
        ax.set_xlim(1, 2) # Bars are going from 0 to 5, so lets crop the plot somewhere in the middle
        ax.grid(False) # Hide grid
        ax.set_facecolor('white') # Make background white
        ax.set_xticks([]) # Remove horizontal ticks
        ax.set_yticks(np.linspace(min(bar_y), max(bar_y), 3)) # Show vertical ticks for min, middle and max
        ax.yaxis.tick_right() # Show vertical ticks on the right 


def corrplot(data, figsize=(8,8), size_scale=250, marker='s', plotname='plot'):
    import seaborn as sns
    sns.set() 

    corr = pd.melt(data.reset_index(), id_vars='index')
    corr.columns = ['x', 'y', 'value']
    plt.figure(figsize=figsize)
    heatmap(
        corr['x'], corr['y'],
        color=corr['value'], color_range=[-1, 1],
        palette=sns.diverging_palette(20, 220, n=2000),
        size=corr['value'].abs(), size_range=[0,1],
        marker=marker,
        x_order=data.columns,
        y_order=data.columns[::-1],
        size_scale=size_scale
    )
    plt.savefig('{}.png'.format(plotname), bbox_inches='tight')
    plt.savefig('{}.pdf'.format(plotname), bbox_inches='tight')
    plt.close('all')

def fillhist(hist, val, wei):
    [hist.Fill(val[i], wei[i]) for i in val.index]

def plot_eff(df, choice, xaxixs, TrainModel, features, bin, importConfig, score_cut, wei, outdir):

    clf_model = pickle.load(open(TrainModel, 'rb'))
    X = df.loc[:,features] # extract the features columns
    df['XGB_score'] = clf_model.predict_proba(X)[:,1] # BDT score

    outcl = importConfig.replace("Tools.TrainConfig_", "")
    os.makedirs(outdir, exist_ok = True) # mkdir the ouput directory of the plots
    xaxisName = xaxixs
    outName = "{}/eff_{}.pdf".format(outdir, choice)
    
    header = ""
    if outcl == "ResolvedMergedClassifier_EB":
        header = "Resolved-Merged Classifier (EB)"
    elif outcl == "ResolvedMergedClassifier_EE":
        header = "Resolved-Merged Classifier (EE)"
    elif outcl == "Merged2GsfID_EB":
        header = "Merged-2Gsf ID (EB)"
    elif outcl == "Merged2GsfID_EE":
        header = "Merged-2Gsf ID (EE)"
    elif outcl == "Merged1GsfID_EB":
        header = "Merged-1Gsf ID (EB)"
    elif outcl == "Merged1GsfID_EE":
        header = "Merged-1Gsf ID (EE)"
    elif outcl == "MergedID_EB":
        header = "Merged ID (EB)"
    elif outcl == "MergedID_EE":
        header = "Merged ID (EE)"
    else:
        print("[ERROR] No specific clf name")
        sys.exit(0)

    df_sig = df[df.Type == "Signal"]
    df_bkg = df[df.Type == "Background"]
    
    print("[INFO] Fill the histograms")
    hM_sig = ROOT.TH1F("hM_sig","", len(bin)-1, np.asarray(bin, "d"))
    hM_sig_pred = ROOT.TH1F("hM_sig_pred","", len(bin)-1, np.asarray(bin, "d"))
    hM_bkg = ROOT.TH1F("hM_bkg","", len(bin)-1, np.asarray(bin, "d"))
    hM_bkg_pred = ROOT.TH1F("hM_bkg_pred","", len(bin)-1, np.asarray(bin, "d"))
    fillhist(hM_sig, df_sig[choice], df_sig[wei])
    fillhist(hM_sig_pred, df_sig[df_sig.XGB_score > score_cut][choice], df_sig[df_sig.XGB_score > score_cut][wei])
    fillhist(hM_bkg, df_bkg[choice], df_bkg[wei])
    fillhist(hM_bkg_pred, df_bkg[df_bkg.XGB_score > score_cut][choice], df_bkg[df_bkg.XGB_score > score_cut][wei])
    print("[INFO] Histograms are done")

    efferr_sig, efferr_bkg = ROOT.TGraphAsymmErrors(), ROOT.TGraphAsymmErrors()
    efferr_sig.BayesDivide(hM_sig_pred, hM_sig)
    efferr_bkg.BayesDivide(hM_bkg_pred, hM_bkg)
    
    # plot setting
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    c = ROOT.TCanvas('c', '', 800, 800)
    c.cd( )
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.07)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetBottomMargin(0.13)
    TickSize = 0.02
    AxisTitleSize = 0.05
    AxisLabelSize = 0.045
    colorN = ["#202020", "#ec0101"]

    # signal efficiency setting
    efferr_sig.GetXaxis().SetTitle(xaxisName)
    efferr_sig.GetYaxis().SetTitle("Efficiency")
    efferr_sig.GetXaxis().SetRangeUser(bin[0], bin[-1])
    efferr_sig.GetYaxis().SetRangeUser(0, 1.3)
    efferr_sig.GetXaxis().SetTickSize(TickSize)
    efferr_sig.GetXaxis().SetTitleSize(AxisTitleSize)
    efferr_sig.GetXaxis().SetLabelSize(AxisLabelSize)
    efferr_sig.GetYaxis().SetTickSize(TickSize)
    efferr_sig.GetYaxis().SetTitleSize(AxisTitleSize)
    efferr_sig.GetYaxis().SetLabelSize(AxisLabelSize)
    efferr_sig.GetXaxis().SetTitleOffset(1.1)
    efferr_sig.GetYaxis().SetTitleOffset(1.3)

    efferr_sig.SetMarkerColor(ROOT.TColor.GetColor(colorN[0]))
    efferr_sig.SetMarkerSize(1.5)
    efferr_sig.SetMarkerStyle(20)
    efferr_sig.SetLineColor(ROOT.TColor.GetColor(colorN[0]))
    efferr_sig.SetLineWidth(2)
    efferr_sig.Draw("AP")

    # background efficiency setting
    efferr_bkg.SetMarkerColor(ROOT.TColor.GetColor(colorN[1]))
    efferr_bkg.SetMarkerSize(1.5)
    efferr_bkg.SetMarkerStyle(20)
    efferr_bkg.SetLineColor(ROOT.TColor.GetColor(colorN[1]))
    efferr_bkg.SetLineWidth(2)
    efferr_bkg.Draw("P same")

    # legend
    l = ROOT.TLegend(0.17, 0.75, 0.7, 0.92)
    l.SetHeader(header)
    l.SetTextSize(0.042)
    l.SetNColumns(2)
    l.AddEntry(efferr_sig, "Signal", "lep")
    l.AddEntry(efferr_bkg, "Background", "lep")
    l.SetFillColor(0)
    l.SetLineColorAlpha(0, 0)
    l.SetFillStyle(0)
    l.Draw()

    CMS_lumi(c, 4, 0, "41.5 fb^{-1}", 2017, True, "Simulation", "", "")

    c.SaveAs(outName)
    c.Close()

def plot_ptetaKin(df, category, label, ptName, etaName, ptBin, etaBin, outdir):

    os.makedirs(outdir, exist_ok = True) # mkdir the ouput directory of the plots

    # pT before and after reweighting
    fig, ax = plt.subplots(1, 2, figsize = (10, 5))
    for i, group_df in df[df['Dataset'] == "Train"].groupby(category):
        group_df[ptName].hist(histtype = 'step', bins = ptBin, alpha = 0.7, label = label[i], ax = ax[0], density = True, ls = '-', weights = group_df["instwei"], linewidth = 2)
        ax[0].set_title("$p_T$ before reweighting")
        ax[0].set_xlabel("$p_T$ [GeV]", fontsize = 12, loc = 'right')
        ax[0].set_ylabel("A.U.", fontsize = 12, loc = 'top')
        ax[0].legend()
        group_df[ptName].hist(histtype = 'step', bins = ptBin, alpha = 0.7, label = label[i], ax = ax[1], density = True, ls = '-', weights = group_df["NewWt"], linewidth = 2)
        ax[1].set_title("$p_T$ after reweighting")
        ax[1].set_xlabel("$p_T$ [GeV]", fontsize = 12, loc = 'right')
        ax[1].set_ylabel("A.U.", fontsize = 12, loc = 'top')
        ax[1].legend()
    outName_pt = outdir + "/pT_rwt.pdf"
    fig.savefig(outName_pt, bbox_inches='tight')
    print("Save fig in %s" %(outName_pt))

    # eta before and after reweighting
    fig, ax = plt.subplots(1, 2, figsize = (10, 5))
    for i, group_df in df[df['Dataset'] == "Train"].groupby(category):
        group_df[etaName].hist(histtype = 'step', bins = etaBin, alpha = 0.7,label = label[i], ax = ax[0],  density = True, ls = '-', weights = group_df["instwei"], linewidth = 2)
        ax[0].set_title("$\eta$ before reweighting")
        ax[0].legend()
        ax[0].set_xlabel("$\eta$", fontsize = 12, loc = 'right')
        ax[0].set_ylabel("A.U.", fontsize = 12, loc = 'top')
        group_df[etaName].hist(histtype = 'step', bins = etaBin, alpha = 0.7,label = label[i], ax = ax[1], density = True, ls = '-', weights = group_df["NewWt"], linewidth = 2)
        ax[1].set_title("$\eta$ after reweighting")
        ax[1].set_xlabel("$\eta$", fontsize = 12, loc = 'right')
        ax[1].set_ylabel("A.U.", fontsize = 12, loc = 'top')
        ax[1].legend()
    outName_eta = outdir + "/eta_rwt.pdf"
    fig.savefig(outName_eta, bbox_inches='tight')
    print("Save fig in %s" %(outName_eta))
    
    plt.close('all')

## ----- for multi-classification -----##
# Reference: https://github.com/javaidnabi31/Multi-class-with-imbalanced-dataset-classification/blob/master/20-news-group-classification.ipynb
import itertools
def plot_confusion_matrix(cm, classes,
                          normalize=False,
                          title = 'Confusion matrix',
                          cmap = plt.cm.Blues,
                          outName = ""):
    """
    This function prints and plots the confusion matrix.
    Normalization can be applied by setting `normalize=True`.
    """
    if normalize:
        cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
        print("Normalized confusion matrix")
    else:
        print('Confusion matrix, without normalization')

    plt.figure(figsize = (8, 6))
    plt.imshow(cm, interpolation = 'nearest', cmap = cmap)
    plt.title(title, fontsize = 15)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize = 10)
    if normalize:
        plt.clim(0, 1)
    tick_marks = np.arange(len(classes))
    plt.xticks(tick_marks, classes, fontsize = 13)
    plt.yticks(tick_marks, classes, fontsize = 13)

    fmt = '.3f' if normalize else 'd'
    thresh = cm.max() / 2.
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        plt.text(j, i, format(cm[i, j], fmt),
                 horizontalalignment = "center", fontsize = 13, 
                 color = "white" if cm[i, j] > thresh else "black")

    plt.ylabel('True class', fontsize = 13, loc = 'top')
    plt.xlabel('Predicted class', fontsize = 13, loc = 'right')
    plt.savefig("{}".format(outName), bbox_inches = 'tight')
    print("Save fig in {}".format(outName))
    plt.close('all')

def fillhist_multi(hist, val):
    for i in val:
        hist.Fill(i)

def plot_eff_multi(df, choice, xaxixs, bin, true_cat, pred_cat, importConfig, outdir):
    ROOT.gROOT.SetBatch() # PyROOT does not display any graphics(root "-b" option)
    outcl = importConfig.replace("Tools.TrainConfig_", "")
    os.makedirs(outdir, exist_ok = True) # mkdir the ouput directory of the plots
    xaxisName = xaxixs
    outName = outdir+"/XGB_Eff_"+choice+".pdf"
    
    header = ""
    if outcl == "MergedID_EB":
        header = "Merged ID (EB)"
    elif outcl == "MergedID_EE":
        header = "Merged ID (EE)"
    else:
        print("[Warning] No specific clf name: use MergedID_EB as default (only use to determine the legend title)")
        header = "Merged ID (EB)"
    
    df_sig = df.query(true_cat+"==1 | "+true_cat+"==0")
    df_bkg = df.query(true_cat+"!=1 & "+true_cat+"!=0")
    
    df_sig_pred = df_sig.query(pred_cat+"==1 | "+pred_cat+"==0")
    df_bkg_pred = df_bkg.query(pred_cat+"==1 | "+pred_cat+"==0")

    hM_sig = ROOT.TH1F("hM_sig","", len(bin)-1, np.asarray(bin, "d"))
    hM_sig_pred = ROOT.TH1F("hM_sig_pred","", len(bin)-1, np.asarray(bin, "d"))
    hM_bkg = ROOT.TH1F("hM_bkg","", len(bin)-1, np.asarray(bin, "d"))
    hM_bkg_pred = ROOT.TH1F("hM_bkg_pred","", len(bin)-1, np.asarray(bin, "d"))
    fillhist_multi(hM_sig, df_sig[choice])
    fillhist_multi(hM_sig_pred, df_sig_pred[choice])

    fillhist_multi(hM_bkg, df_bkg[choice])
    fillhist_multi(hM_bkg_pred, df_bkg_pred[choice])

    efferr_sig, efferr_bkg = ROOT.TGraphAsymmErrors(), ROOT.TGraphAsymmErrors()
    efferr_sig.BayesDivide(hM_sig_pred, hM_sig)
    efferr_bkg.BayesDivide(hM_bkg_pred, hM_bkg)
    
    # plot setting
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    c = ROOT.TCanvas('c', '', 800, 800)
    c.cd( )
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetTopMargin(0.07)
    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetBottomMargin(0.13)
    TickSize = 0.02
    AxisTitleSize = 0.05
    AxisLabelSize = 0.045
    colorN = ["#202020", "#ec0101"]

    # signal efficiency setting
    efferr_sig.GetXaxis().SetTitle(xaxisName)
    efferr_sig.GetYaxis().SetTitle("Efficiency")
    efferr_sig.GetXaxis().SetRangeUser(bin[0], bin[-1])
    efferr_sig.GetYaxis().SetRangeUser(0, 1.5)
    efferr_sig.GetXaxis().SetTickSize(TickSize)
    efferr_sig.GetXaxis().SetTitleSize(AxisTitleSize)
    efferr_sig.GetXaxis().SetLabelSize(AxisLabelSize)
    efferr_sig.GetYaxis().SetTickSize(TickSize)
    efferr_sig.GetYaxis().SetTitleSize(AxisTitleSize)
    efferr_sig.GetYaxis().SetLabelSize(AxisLabelSize)
    efferr_sig.GetXaxis().SetTitleOffset(1.1)
    efferr_sig.GetYaxis().SetTitleOffset(1.3)

    efferr_sig.SetMarkerColor(ROOT.TColor.GetColor(colorN[0]))
    efferr_sig.SetMarkerSize(1.5)
    efferr_sig.SetMarkerStyle(20)
    efferr_sig.SetLineColor(ROOT.TColor.GetColor(colorN[0]))
    efferr_sig.SetLineWidth(2)
    efferr_sig.Draw("AP")

    # background efficiency setting
    efferr_bkg.SetMarkerColor(ROOT.TColor.GetColor(colorN[1]))
    efferr_bkg.SetMarkerSize(1.5)
    efferr_bkg.SetMarkerStyle(20)
    efferr_bkg.SetLineColor(ROOT.TColor.GetColor(colorN[1]))
    efferr_bkg.SetLineWidth(2)
    efferr_bkg.Draw("P same")

    # legend
    l = ROOT.TLegend(0.17, 0.75, 0.7, 0.92)
    l.SetHeader(header)
    l.SetTextSize(0.042)
    # l.SetNColumns(2)
    l.AddEntry(efferr_sig, "Signal (Merged-2Gsf + Merged-1Gsf)", "lep")
    l.AddEntry(efferr_bkg, "Background (GJets + DYJets)", "lep")
    l.SetFillColor(0)
    l.SetLineColorAlpha(0, 0)
    l.SetFillStyle(0)
    l.Draw()

    CMS_lumi(c, 4, 0, "", 2017, True, "Simulation", "", "")

    c.SaveAs(outName)
    c.Close()

# Reference: https://stackoverflow.com/a/29643643
def hex2rgb(h, a):
    h = h.lstrip('#')
    rgb = tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))
    return (rgb[0], rgb[1], rgb[2], a)

# def DrawFeatureHist(df, feature, bininfo, XaxisName, yaxisunit, histcolor = ['#0061a8', '#e84545'], logy = True, y_axisscale = 1.5, drawxoverflow = True, drawxunderflow = False, outname = 'testfeatureplot', wei='' ):
def DrawFeatureHist_multi(df, feature, num_class, Class, bininfo, XaxisName, yaxisunit, histcolor, logy = False, y_axisscale = 1.3, outname = "testfeatureplot", wei = "", category = "EleType"):

    plt.figure(figsize = (6, 6))
    binContentMax_list, binContentMin_list = [], []
    for i in range(num_class):
        binContent, binEdget, patches = plt.hist(
            df[df[category] == i][feature],
            range = (bininfo[1], bininfo[2]) , 
            bins = bininfo[0],
            ls = "-", lw = 2, fc = hex2rgb(histcolor[i], 0), ec = histcolor[i],
            histtype = "stepfilled",
            weights = df[df[category] == i][wei],
            density = True,
            label = Class[i]
        )
        # print(binContent, Class[i])
        binContentMax_list.append(np.max(binContent))
        binContentMin_list.append(np.min(binContent[binContent != 0]))

    ymaxval = np.max(binContentMax_list)
    yminval = np.min(binContentMin_list) if np.min(binContentMin_list) < 1E5 else 1E5

    ax = plt.gca()
    if logy:
        ax.set_yscale("log")
        ax.set_ylim(yminval * 0.9, ymaxval * y_axisscale)
    else:
        ax.set_ylim(0, ymaxval * y_axisscale)
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.set_xlim(bininfo[1], bininfo[2])
    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())

    ax.set_ylabel("A.U.", fontsize = 15, loc = "top")
    ax.set_xlabel(XaxisName, fontsize = 15, loc = "right")
    
    # lumi information text
    ax.text(0, 1, "CMS", horizontalalignment = "left", verticalalignment = "bottom", transform=ax.transAxes, fontsize = 15, fontweight = "bold")
    # ax.text(0.13, 1, "$\it{Simulation}$", horizontalalignment = "left", verticalalignment = "bottom", transform = ax.transAxes, fontsize = 13)
    ax.text(1, 1, "$(13TeV,\ 2017)$", horizontalalignment = "right", verticalalignment = "bottom", transform = ax.transAxes, fontsize = 13)
    
    ax.legend(loc = "best", edgecolor= "none", fontsize = 13)
    ax.tick_params( # major tick 
        direction = "in", length = 8,
        labelsize = 13, top = True, right = True
    )
    ax.tick_params( # minor tick 
        direction = "in", length = 4, which = "minor",
        labelsize = 13, top = True, right = True
    )

    ax.legend(loc = "best", edgecolor= "none", fontsize = 13)
    plt.tight_layout()
    plt.draw()
    plt.savefig("{}.pdf".format(outname))
    print("Save fig in {}.pdf".format(outname))
    plt.close()


def pltSty(ax, xName = "x-axis", yName = "y-axis", TitleSize = 15, LabelSize = 15, TickSize = 13, MajTickLength = 7, MinTickLength = 4, yAuto = True):
    ax.set_xlabel(xName, fontsize = LabelSize, loc = "right")
    ax.set_ylabel(yName, fontsize = LabelSize, loc = "top")
    ax.text(1, 1, "(13 TeV)", horizontalalignment = "right", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize)
    ax.text(0, 1, "CMS", horizontalalignment = "left", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize * 1.3, fontweight = "bold")
    ax.text(TitleSize * 0.009 + 0.02, 1, "Simulation", horizontalalignment = "left", verticalalignment = "bottom", transform = ax.transAxes, fontsize = TitleSize)

    ax.xaxis.set_major_locator(AutoLocator())
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    if (yAuto):
        ax.yaxis.set_major_locator(AutoLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(direction = "in", length = MajTickLength, labelsize = TickSize, top = True, right = True)
    ax.tick_params(direction = "in", length = MinTickLength, which = "minor", labelsize = TickSize, top = True, right = True)


def DrawStatistic(df, outdir):
    fig, axes = plt.subplots(1, 1, figsize = (10, 5))
    kplot = sns.countplot(x = "Class", data = df, ax = axes, hue = "Dataset", palette = ["#432371","#FAAE7B","black"])
    for p in kplot.patches:
        kplot.annotate(format(p.get_height(), "d"), (p.get_x() + p.get_width() / 2., p.get_height()), ha = "center", va = "center", xytext = (0, 5), textcoords = "offset points", size = 10)
    axes.set_title("Number of samples")
    #axes.set_yscale("log")
    plt.savefig(outdir+"/TotalStat_TrainANDTest.pdf", bbox_inches = "tight")
    print("Save fig in {}".format(outdir+"/TotalStat_TrainANDTest.pdf"))
    plt.close("all")

# def corre(df, features, Classes = [""], MVA = "", outdir):
#     for C in Classes:
#         fig, axes = plt.subplots(1, 1, figsize = (len(features)/2, len(features)/2))
#         cor = df.loc[(df["Class"] == str(C))][features].corr()
#         sns.heatmap(
#             cor, annot = True, fmt = ".2f", cmap = plt.cm.Blues, 
#             ax = axes, annot_kws = {"size": 25/np.sqrt(len(cor))},
#             vmin = -1, vmax = 1
#         )
#         axes.tick_params(axis = "x", labelsize = len(features)/2)
#         axes.tick_params(axis = "y", labelsize = len(features)/2)
#         #axes.set_xticklabels(MVA["features"],fontsize= len(MVA["features"]))
#         #axes.set_yticklabels(MVA["features"],fontsize= len(MVA["features"]))
#         plt.savefig(outdir+"/"+MVA+"/"+MVA+"_"+C+"_CORRELATION"+".pdf", bbox_inches = "tight")
#         print("Save fig in {}".format(outdir+"/"+MVA+"/"+MVA+"_"+C+"_CORRELATION"+".pdf"))
#         plt.close("all")

