from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas
import time

def GetAnalysisMask():
    util = RAT.DU.Utility.Get()
    util.LoadDBAndBeginRun()
    dcmask = RAT.GetDataCleaningWord('analysis_mask')
    return dcmask

def PrepData_AmBebkg(ambefiles=[],bkgfiles=[],var="nhits",nbins=10,xmin=10.0,xmax=100.0,dcmask=None):
    '''
    Takes in a list of AmBe central root files and background files formatted
    to have prompt and delayed IBD candidates chosen. Preps the data as a
    PandasFrame for plotting. NOTE: AmBe root files are
    Assumed to have prompt and delayeds paired into IBD candidates with
    variables labeled like "nhits_p" and "nhits_d" for prompt and delayed'''

    bkgdata = ROOT.TChain("CombinedOutput")
    for rf in bkgfiles:
        bkgdata.Add(rf)
    ambedata = ROOT.TChain("CombinedOutput")
    for af in ambefiles:
        ambedata.Add(af)
    #Now, make a new dictionary object: key is cutname, value is histogram
    bkgdata.Draw("%s_p>>h_allbkg_p(%i,%f,%f)"% (var,nbins,xmin,xmax),
            "interevent_time<500000&&fitValid_p==1&&posr_p<5000","goff")
    h_allbkg_p = gDirectory.Get("h_allbkg_p")
    h_allbkg_p.Sumw2()
    bkgdata.Draw("%s_d>>h_allbkg_d(%i,%f,%f)"% (var,nbins,xmin,xmax),
            "interevent_time<500000&&fitValid_d==1&&posr_d<5000","goff")
    h_allbkg_d = gDirectory.Get("h_allbkg_d")
    h_allbkg_d.Sumw2()
    ambedata.Draw("%s_p>>h_allambe_p(%i,%f,%f)"% (var,nbins,xmin,xmax),
            "interevent_time<500000&&fitValid_p==1&&posr_p<5000","goff")
    h_allambe_p = gDirectory.Get("h_allambe_p")
    h_allambe_p.Sumw2()
    ambedata.Draw("%s_d>>h_allambe_d(%i,%f,%f)"% (var,nbins,xmin,xmax),
            "interevent_time<500000&&fitValid_d==1&&posr_d<5000","goff")
    h_allambe_d = gDirectory.Get("h_allambe_d")
    h_allambe_d.Sumw2()
    print("MADE ALL HISTOGRAMS")

    allclasssacs = {}
    labeldict = ["AmBe_Prompt","AmBe_Delayed","Bkg_Prompt","Bkg_Delayed"]
    for dat in labeldict:
        graphdict={}
        vardat, fs,fs_unc =(), (), () #pandas wants ntuples
        h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", nbins,xmin,xmax)
        h_cut_FracFlagged.Sumw2()
        if dat == "AmBe_Prompt":
            ambedata.Draw("%s_p>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "interevent_time<500000&&fitValid_p==1&&posr_p<5000&&((dcFlagged_p&%s)!=%s)"%(dcmask,dcmask),"goff")
            h_all = h_allambe_p
        if dat == "AmBe_Delayed":
            ambedata.Draw("%s_d>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "interevent_time<500000&&fitValid_d==1&&posr_d<5000&&((dcFlagged_d&%s)!=%s)"%(dcmask,dcmask),"goff")
            h_all = h_allambe_d
        if dat == "Bkg_Prompt":
            bkgdata.Draw("%s_p>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "interevent_time<500000&&fitValid_p==1&&posr_p<5000&&((dcFlagged_p&%s)!=%s)"%(dcmask,dcmask),"goff")
            h_all = h_allbkg_p
        if dat == "Bkg_Delayed":
            bkgdata.Draw("%s_d>>h_flagged(%i,%f,%f)" % (var,nbins,xmin,xmax),
                    "interevent_time<500000&&fitValid_d==1&&posr_d<5000&&((dcFlagged_d&%s)!=%s)"%(dcmask,dcmask),"goff")
            h_all = h_allbkg_d
        print("CUT: " + str(dat))
        h_flagged = gDirectory.Get("h_flagged")
        h_flagged.Sumw2()
        h_cut_FracFlagged.Divide(h_flagged,h_all,1.,1.,"b")
        for i in xrange(int(h_cut_FracFlagged.GetNbinsX())):
            vardat =  vardat + (float(h_cut_FracFlagged.GetBinWidth(i)) + float(h_cut_FracFlagged.GetBinLowEdge(i)),)
            fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
            fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
        graphdict["vardat"] = vardat
        graphdict["fractional_sacrifice"] = fs
        graphdict["fs_uncertainty"] = fs_unc
        allclasssacs[dat]= pandas.DataFrame(data=graphdict)
        del h_cut_FracFlagged
        del h_flagged
    #graphdict has the sacrifice information for each cut. Now, let's plot it.
    return allclasssacs,var 



def Plot_SacAmBebkg(allclasssacs,variable="nhits"):
    '''Takes in A dictionary of pandaframes output from PrepData_BkgAmBe
    and plots it with seaborn/matplotlib'''

    sns.set_style("whitegrid")
    xkcd_colors = ['slate blue', 'fluro green', 'dark green', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
            'aqua blue','vomit', 'black','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(allclasssacs)))
    for cut in allclasssacs:
        plt.errorbar(x=allclasssacs[cut].vardat, 
                y=allclasssacs[cut].fractional_sacrifice,
                yerr=allclasssacs[cut].fs_uncertainty,
                linestyle='none', marker='o', label=cut, markersize=8,
                elinewidth=3, capsize=0)
    plt.legend(loc=3)
    plt.yscale("log")
    plt.ylabel("Fractional sacrifice",fontsize=34)
    plt.xlabel(variable,fontsize=34)
    plt.tick_params(labelsize=32)
    plt.title("Fractional sacrifice of prompt/delayed candidates in\n"+\
            "AmBe central data & background data with Valid Fit and radius<5 m",fontsize=36)
    plt.ion()
    plt.show()
