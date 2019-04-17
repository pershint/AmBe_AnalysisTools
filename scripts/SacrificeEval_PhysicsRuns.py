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

sns.set_style("darkgrid")
sns.set_context("poster")

DATADIR = "./rootfiles/all"
PCONFIGFILE = "./config/config_prompt_default.json"
DCONFIGFILE = "./config/config_delayed_default.json"
BOTHCONFIGFILE = "./config/full_pair_config.json"

if __name__=='__main__':
    print("LET US BEGIN")
    with open(PCONFIGFILE,"r") as f:
        pconfig = json.load(f)
    with open(DCONFIGFILE,"r") as f:
        dconfig = json.load(f)
    with open(BOTHCONFIGFILE,"r") as f:
        bconfig = json.load(f)
    
    #First, get the rootfile names for all data and MC
    datafiles = glob.glob("%s/data/*Physics*"%(DATADIR))
    mcfiles = glob.glob("%s/mc/*.root"%(DATADIR))
    #Now, start up the sacrifice analyzer
    SacAnalyzer = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig,bconfig)
    #Set the number of bins we want when plotting/analyzing
    SacAnalyzer.SetPromptBinNumber(29)
    SacAnalyzer.SetDelayedBinNumber(11)
    #Now, analyze the sacrifice per bin.  Start with nhits as the variable
    #accepts nhits, energy, udotr, or posr3 for nice plotting right now
    defaultdat, defaultmeta = SacAnalyzer.AnalyzeData(var='nhitsCleaned')
    #Now, lets try to plot out the prompt event's sacrifice for the variable analyze
    #We will make the plot for the MC files' prompt events
    SacAnalyzer.ShowSacrificePlot(evtype='prompt',dattype='data',fittotal=False)
    #SacAnalyzer.ShowSacrificePlot(evtype='prompt',dattype='MC',fittotal=False)
    SacAnalyzer.ShowSacrificePlot(evtype='delayed',dattype='data',fittotal=False)
    #SacAnalyzer.ShowSacrificePlot(evtype='delayed',dattype='MC',fittotal=False)
    #Lets do a data/MC comparison as well for the delayed events
    SacAnalyzer.DataMCSacCompare(evtype="prompt")
    #And let's look at the Data/MC sacrifice ratio
    SacAnalyzer.PlotDataMCSacRatio(evtype='prompt',fittotal=False)
    SacAnalyzer.PlotDataMCAccRatio(evtype='prompt',fittotal=False)
    SacAnalyzer.PlotDataMCAccRatio(evtype='delayed',fittotal=False)
    #Get a histogram of the clean distribution for the input variable
    CleanIThist,CleanIThistMeta = SacAnalyzer.DrawHist(evtype="prompt",var='interevent_time',
                                                     dattype="data",xmin=3000,
                                                     xmax=990000.,nbins=50,
                                                     addlROOTcuts="interevent_time>3000",
                                                     dctype="clean")
    PDirtyIThist,PDirtyhistMeta = SacAnalyzer.DrawHist(evtype="prompt",var='interevent_time',
                                                     dattype="data",xmin=3000,
                                                     xmax=990000.,nbins=50,
                                                     addlROOTcuts="interevent_time>3000",
                                                     dctype="dirty")
    
    DDirtyIThist,DDirtyhistMeta = SacAnalyzer.DrawHist(evtype="delayed",var='interevent_time',
                                                     dattype="data",xmin=3000,
                                                     xmax=990000.,nbins=50,
                                                     addlROOTcuts="interevent_time>3000",
                                                     dctype="dirty")
    DDirtynhithist,DDirtynhithistMeta = SacAnalyzer.DrawHist(evtype="delayed",var='nhitsCleaned',
                                                     dattype="data",xmin=4,
                                                     xmax=15,nbins=11,
                                                     addlROOTcuts="interevent_time>3000",
                                                     dctype="dirty")
    PDirtynhithist,PDirtynhithistMeta = SacAnalyzer.DrawHist(evtype="prompt",var='nhitsCleaned',
                                                     dattype="data",xmin=15,
                                                     xmax=45,nbins=15,
                                                     addlROOTcuts="interevent_time>3000",
                                                     dctype="dirty")
    DCleannhithist,DCleannhithistMeta = SacAnalyzer.DrawHist(evtype="delayed",var='nhitsCleaned',
                                                     dattype="data",xmin=4,
                                                     xmax=15,nbins=11,
                                                     addlROOTcuts="interevent_time>3000",
                                                     dctype="clean")
    PCleannhithist,PCleannhithistMeta = SacAnalyzer.DrawHist(evtype="prompt",var='nhitsCleaned',
                                                     dattype="data",xmin=15,
                                                     xmax=45,nbins=15,
                                                     addlROOTcuts="interevent_time>3000",
                                                     dctype="clean")
    #Neat.  Now, let's fit to our clean and dirty distributions
    doubleexp = lambda x,A1,l1,A2,l2: A1*np.exp(-l1*x) + A2*np.exp(-l2*x)
    singleexp = lambda x,A1,l1: A1*np.exp(-l1*x)
    mycleanfitter = sf.Fitter(datax=CleanIThist["x"],datay=CleanIThist["y"],
                      datasigma=CleanIThist["y_unc"])
    mycleanfitter.SetFitFunction(singleexp, 2)
    initvars = [200., 1./20000.]
    popt, pcov = mycleanfitter.RunFit(initvars)
    plt.errorbar(x=CleanIThist["x"],y=CleanIThist["y"],yerr=CleanIThist["y_unc"],
                 linestyle='none',marker='o',markersize=7)
    bestlin1 = popt[0]*np.exp(-popt[1]*CleanIThist["x"])
    plt.plot(CleanIThist["x"],bestlin1,label=r"Best fit, $\tau=%s \, ms$"%(str(np.round((1./(1E6*popt[1])),2))))
    plt.legend(loc=1)
    plt.xlabel("Interevent Time (ns)")
    plt.ylabel("Events")
    plt.title(("Interevent time distribution of AmBe Physics data \n" +
               "Delayed and Prompt Pass livetime, trigger, and data cleaning cuts"))
    plt.show()
    
    mydirtyfitter = sf.Fitter(datax=PDirtyIThist["x"],datay=PDirtyIThist["y"],
                      datasigma=PDirtyIThist["y_unc"])
    mydirtyfitter.SetFitFunction(singleexp, 2)
    initvars = [200., 1./20000.]
    popt, pcov = mydirtyfitter.RunFit(initvars)
    plt.errorbar(x=PDirtyIThist["x"],y=PDirtyIThist["y"],yerr=PDirtyIThist["y_unc"],
                 linestyle='none',marker='o',markersize=7)
    bestfitline = popt[0]*np.exp(-popt[1]*PDirtyIThist["x"])
    plt.plot(PDirtyIThist["x"],bestfitline,label=r"Best fit, $\tau=%s \, ms$"%(str(np.round((1./(1E6*popt[1])),2))))
    plt.legend(loc=1)
    plt.xlabel("Interevent Time (ns)")
    plt.ylabel("Events")
    plt.title(("Interevent time distribution of AmBe Physics data \n" +
               "Delayed pass cuts, Prompt only fails DC cuts"))
    plt.show()
    
    myddirtyfitter = sf.Fitter(datax=DDirtyIThist["x"],datay=DDirtyIThist["y"],
                      datasigma=DDirtyIThist["y_unc"])
    myddirtyfitter.SetFitFunction(singleexp, 2)
    initvars = [200., 1./20000.]
    popt, pcov = myddirtyfitter.RunFit(initvars)
    plt.errorbar(x=DDirtyIThist["x"],y=DDirtyIThist["y"],yerr=DDirtyIThist["y_unc"],
                 linestyle='none',marker='o',markersize=7)
    bestfitline = popt[0]*np.exp(-popt[1]*DDirtyIThist["x"])
    plt.plot(DDirtyIThist["x"],bestfitline,label=r"Best fit, $\tau=%s \, ms$"%(str(np.round((1./(1E6*popt[1])),2))))
    plt.legend(loc=1)
    plt.xlabel("Interevent Time (ns)")
    plt.ylabel("Events")
    plt.title(("Interevent time distribution of AmBe Physics data \n" +
               "Prompt pass cuts, Delayed only fails DC cuts"))
    plt.show()

    plt.errorbar(x=DDirtynhithist["x"],y=DDirtynhithist["y"],yerr=DDirtynhithist["y_unc"],
                 xerr=0.5,
                 linestyle='none',marker='o',markersize=7,color="r",
                 label="Delayed dist. (delayed dirty, prompt clean)")
    plt.errorbar(x=PDirtynhithist["x"],y=PDirtynhithist["y"],yerr=PDirtynhithist["y_unc"],
                 xerr=1.0,
                 linestyle='none',marker='o',markersize=7,color="purple",
                 label="Prompt dist. (prompt dirty, delayed clean)")
    plt.xlabel("nhitsCleaned")
    plt.ylabel("Events")
    plt.legend()
    plt.title(("nhitsCleaned distribution of AmBe Physics data \n" +
               "Distribution for prompt/delayed candidates failing DC"))
    plt.show()

    plt.errorbar(x=DCleannhithist["x"],y=DCleannhithist["y"],yerr=DCleannhithist["y_unc"],
                 xerr=0.5,
                 linestyle='none',marker='o',markersize=7,color="black",
                 label="Delayed dist (Prompt clean, delayed clean)")
    plt.errorbar(x=PCleannhithist["x"],y=PCleannhithist["y"],yerr=PCleannhithist["y_unc"],
                 xerr=1.0,
                 linestyle='none',marker='o',markersize=7,color="blue",
                 label="Prompt dist (Prompt clean, delayed clean)")
    plt.xlabel("nhitsCleaned")
    plt.ylabel("Events")
    plt.legend()
    plt.title(("nhitsCleaned distribution of AmBe Physics data for prompt/delayed candidates\n" +
               "Both pass trigger, livetime, and DC cuts"))
    plt.show()
