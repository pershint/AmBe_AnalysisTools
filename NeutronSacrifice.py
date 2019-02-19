#This script is a demonstration of how to use the AmBeDCSacrifice class
#for processing prompt/delayed candidate data and evaluating DC sacrifices
#for both data and MC files.

import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
import lib.AmBeDCSacrifice as ab
import lib.EarlyLatePlots as el
import lib.SimpleFitter as sf
import lib.NeutronSacCalculator as nsc

sns.set_style("darkgrid")
sns.set_context("poster")

DATADIR = "./rootfiles/all"
PCONFIGFILE = "./config/config_prompt_default.json"
DCONFIGFILE = "./config/config_delayed_default.json"
BOTHCONFIGFILE = "./config/full_pair_config.json"

EARLYCONFIGFILE = "./config/early_pair_config.json"
LATECONFIGFILE = "./config/late_pair_config.json"

def GetRanges(dataframe):
    ranges = []
    xvalues = dataframe["vardat"]
    width = xvalues[1] - xvalues[0]
    for i in xrange(xvalues):
        thisrange = [xvalues[i]- width/2.,xvalues[i]+width/2.]
        ranges.append(thisrange)
    return ranges

if __name__=='__main__':
    print("LET US BEGIN")
    with open(PCONFIGFILE,"r") as f:
        pconfig = json.load(f)
    with open(DCONFIGFILE,"r") as f:
        dconfig = json.load(f)
    with open(BOTHCONFIGFILE,"r") as f:
        bconfig = json.load(f)
    with open(EARLYCONFIGFILE,"r") as f:
        econfig = json.load(f)
    with open(LATECONFIGFILE,"r") as f:
        lconfig = json.load(f)
    
    #First, get the rootfile names for all data and MC
    datafiles = glob.glob("%s/data/*.root"%(DATADIR))
    mcfiles = glob.glob("%s/mc/*.root"%(DATADIR))
    #Now, start up the sacrifice analyzer
    EarlySacAnalyze = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig,econfig)
    LateSacAnalyze = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig,lconfig)
    #Set the number of bins we want when plotting/analyzing
    EarlySacAnalyze.SetPromptBinNumber(6)
    EarlySacAnalyze.SetDelayedBinNumber(6)
    LateSacAnalyze.SetPromptBinNumber(6)
    LateSacAnalyze.SetDelayedBinNumber(6)
    #Now, analyze the sacrifice per bin.  Start with nhits as the variable
    #accepts nhits, energy, udotr, or posr3 for nice plotting right now
    earlysacs, earlymeta = EarlySacAnalyze.AnalyzeData(var='nhitsCleaned')
    latesacs, latemeta = LateSacAnalyze.AnalyzeData(var='nhitsCleaned')

    SacAnalyzer = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig,bconfig)
    #Now, analyze the sacrifice per bin.  Start with nhits as the variable
    #accepts nhits, energy, udotr, or posr3 for nice plotting right now
    ##defaultdat, defaultmeta = SacAnalyzer.AnalyzeData(var='nhitsCleaned')
    ##nhitranges = GetRanges(defaultdat) #FIXME: check values in vardat are RHS of bins
    nhitranges=[[4,8],[8,12]]
    for i in xrange(len(nhitranges)):
        earliest_time = econfig["interevent_time_low"]
        latest_time = lconfig["interevent_time_high"]
        lownhit = str(nhitranges[i][0])
        highnhit = str(nhitranges[i][1])
        IThist,histMeta = SacAnalyzer.DrawCleanHist(evtype="pair",var='interevent_time',
                                                    dattype="data",xmin=earliest_time,
                                                    xmax=latest_time,nbins=50,
                                                    addlROOTcuts="nhitsCleaned_d>%s&&nhitsCleaned_d<%s"%(lownhit,highnhit))
        #Neat.  Now, let's fit to this and plot it
        doubleexp = lambda x,A1,l1,A2,l2: A1*np.exp(-l1*x) + A2*np.exp(-l2*x)
        singleexp = lambda x,A,l: A*np.exp(-l*x)
        myfitter = sf.Fitter(datax=IThist["x"],datay=IThist["y"],
                          datasigma=IThist["y_unc"])
        myfitter.SetFitFunction(doubleexp, 4)
        initvars = [200., 1./20000., 15.,1./2.E7]
        popt, pcov = myfitter.RunFit(initvars)
        plt.errorbar(x=IThist["x"],y=IThist["y"],yerr=IThist["y_unc"],
                     linestyle='none',marker='o',markersize=5)
        bestfitline = doubleexp(IThist["x"],popt[0],popt[1],popt[2],popt[3])
        corry = singleexp(IThist["x"],popt[0],popt[1])
        uncorry = singleexp(IThist["x"],popt[2],popt[3])
        plt.plot(IThist["x"],bestfitline,label="best fit")
        plt.plot(IThist["x"],corry,label="Neutron capture fit")
        plt.plot(IThist["x"],uncorry,label="Backgrounds fit")
        plt.xlabel("Interevent Time (ns)")
        plt.ylabel("Events")
        plt.title("Data cleaned interevent time distribution of AmBe data")
        plt.show()
        #Neato.  Now, load these into the data cleaning sacrifice class
        earlywindow = [econfig["interevent_time_low"],
                       econfig["interevent_time_high"]]
        latewindow = [lconfig["interevent_time_low"],
                      lconfig["interevent_time_high"]]
        esac = earlysacs["delayed"]["data"]["total"].fractional_sacrifice[i]
        lsac = latesacs["delayed"]["data"]["total"].fractional_sacrifice[i]
        NSCalc = nsc.NeutronSacrificeCalculator(IThist["x"],corry, uncorry,
                earlywindow, latewindow, esac, lsac) 
