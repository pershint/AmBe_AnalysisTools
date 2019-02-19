from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.optimize as spc
import pandas
import time

def GetAnalysisMask():
    util = RAT.DU.Utility.Get()
    util.LoadDBAndBeginRun()
    dcmask = RAT.GetDataCleaningWord('analysis_mask')
    return dcmask

def flatline(x,b):
    return b

def weighted_stdev(vals,valavg, uncs):
    weights = 1/(uncs**2)
    return np.sqrt((len(weights)*np.sum(weights*(vals-valavg)**2))/((len(weights)-1)*np.sum(weights)))

def SacVSVar_Data(rootfiles=[],vartype = "prompt", precuts=None,var="nhits",nbins=10,xmin=10.0,xmax=100.0,dcmask=None):
    '''
    Takes in a rootfile and returns a PandaFrame object that can be used
    for plotting in matplotlib.  Returns nhits vs. fractional sacrifice for
    Each cut in the given dcmask.
    '''
    if vartype=="prompt":
        vt = "p"
    elif vartype=="delayed":
        vt = "d"
    else:
        print("Choose either prompt or delayed variable type...")
        return
    if dcmask == None:
        print("No DC mask given.  Please give a DC mask of cuts to plot")
        return
    data = ROOT.TChain("CombinedOutput")
    for rf in rootfiles:
        data.Add(rf)
    fullmask = mb.get_dcwords()
    print(fullmask)
    plotmask = {}
    for cut in fullmask:
        if 'prescale' in fullmask[cut]:
            continue
        if 'waterblind' in fullmask[cut]:
            continue
        if mb.binary_bit_to_int(cut)&dcmask==mb.binary_bit_to_int(cut):
            print("WE ARE ADDING IN CUT: " + str(cut))
            print("THAT IS... CUT: " + str(fullmask[cut]))
            plotmask[cut] = fullmask[cut]
    plotmask[dcmask] = "total"
    #Now, make a new dictionary object: key is cutname, value is histogram
    if precuts is not None: 
        data.Draw("%s_%s>>h_allevents(%i,%f,%f)"% (var,vt,nbins,xmin,xmax),
                precuts,"goff")
    else: 
        data.Draw("%s_%s>>h_allevents(%i,%f,%f)"% (var,vt,nbins,xmin,xmax),"","goff")
    h_allevents = gDirectory.Get("h_allevents")
    h_allevents.Sumw2()
    cutnames = ()
    allcutsacs = {}
    for cut in plotmask:
        if cut!= dcmask:
            cutint = mb.binary_bit_to_int(cut)
        else:
            cutint = cut
        graphdict={}
        vardat, fs,fs_unc =(), (), () #pandas wants ntuples
        h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", nbins,xmin,xmax)
        h_cut_FracFlagged.Sumw2()
        if precuts is not None:
            data.Draw("%s_%s>>h_flagged(%i,%f,%f)" % (var,vt,nbins,xmin,xmax),
                    "%s&&((dcFlagged_%s&%i)!=%i)" % (precuts,vt,cutint,cutint),"goff")
        else:
            data.Draw("%s_%s>>h_flagged(%i,%f,%f)" % (var,vt,nbins,xmin,xmax),
                    "((dcFlagged_%s&%i)!=%i)" % (vt,cutint,cutint),"goff")
        h_flagged = gDirectory.Get("h_flagged")
        h_flagged.Sumw2()
        h_cut_FracFlagged.Divide(h_flagged,h_allevents,1.,1.,"b")
        for i in xrange(h_cut_FracFlagged.GetNbinsX()):
            if i==0:
                continue
            vardat =  vardat + (float(h_cut_FracFlagged.GetBinWidth(i)) + float(h_cut_FracFlagged.GetBinLowEdge(i)),)
            fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
            fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
        graphdict["vardat"] = vardat
        graphdict["fractional_sacrifice"] = fs
        graphdict["fs_uncertainty"] = fs_unc
        allcutsacs[plotmask[cut]]= pandas.DataFrame(data=graphdict)
        del h_cut_FracFlagged
        del h_flagged
    #graphdict has the sacrifice information for each cut. Now, let's plot it.
    meta = {"binwidth":((xmax-xmin)/float(nbins)),"variable":"%s_%s"%(var,vt),"precuts":precuts}
    return allcutsacs, meta

def GetTopSacs(allcutsacs,topnumber=5,dcmask=None):
    if len(allcutsacs) <= topnumber:
        print("You already only have five cuts.  Not combining any")
        return
    names = []
    sacrifices = []
    for cut in allcutsacs:
        sacrifices.append(np.average(allcutsacs[cut]["fractional_sacrifice"]))
        names.append(cut)
    sortednames = [x for _,x in sorted(zip(sacrifices,names))]
    topnames = sortednames[(len(sortednames)-(topnumber+1)):len(sortednames)]
    return topnames
     

def SacVSVar_Plot(allcutsacs,topcuts=None,metadata={"binwidth":10.0,"variable":"nhits"},fittotal=True):
    sns.set_style("whitegrid")
    xkcd_colors = ['black','slate blue', 'fluro green', 'brown', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'leaf',
            'aqua blue','vomit', 'red','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(allcutsacs)))
    if topcuts is None:
        for cut in allcutsacs:
            plt.errorbar(x=allcutsacs[cut].vardat, 
                    y=allcutsacs[cut].fractional_sacrifice,
                    yerr=allcutsacs[cut].fs_uncertainty,
                    linestyle='none', marker='+', label=cut, markersize=8,
                    elinewidth=2, capsize=0)
    else:
        fracsum,fracuncsum=[],[]
        for cut in allcutsacs:
            if cut in topcuts:
                plt.errorbar(x=allcutsacs[cut].vardat, y=allcutsacs[cut].fractional_sacrifice,
                        yerr=allcutsacs[cut].fs_uncertainty, linestyle='none',
                        marker='o', label=cut, capsize=0, elinewidth=2, markersize=5)
            else:
                fracsum.append(allcutsacs[cut].fractional_sacrifice)
                fracuncsum.append(allcutsacs[cut].fs_uncertainty)
        plt.errorbar(x=allcutsacs[cut].vardat,y=sum(fracsum),yerr=sum(fracuncsum),
                linestyle='none',marker='o', capsize=0, elinewidth=2, label='All other cuts', markersize=5)

    if fittotal is True:
        popt, pcov = spc.curve_fit(flatline, allcutsacs["total"].vardat, allcutsacs["total"].fractional_sacrifice,
                p0=[0.02], sigma=allcutsacs["total"].fs_uncertainty)
        #one standard deviation
        print("BEST FIT: " + str(popt))
        print("PCOVARIANCE: " + str(pcov))
        stdev = np.sqrt(np.diag(pcov))
        plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
        print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(weighted_stdev(allcutsacs["total"].fractional_sacrifice,
            float(popt[0]),allcutsacs["total"].fs_uncertainty)))
    plt.legend(loc=3)
    plt.yscale("log")
    plt.ylabel("Fractional sacrifice",fontsize=32)
    plt.xlabel(metadata["variable"],fontsize=32)
    plt.tick_params(labelsize=30)
    plt.title("Fractional sacrifice of dataset by data cleaning cuts",fontsize=34)
    plt.ion()
    plt.show()
