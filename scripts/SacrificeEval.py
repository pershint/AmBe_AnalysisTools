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
BOTHCONFIGFILE = "./config/early_pair_config.json"

if __name__=='__main__':
    print("LET US BEGIN")
    with open(PCONFIGFILE,"r") as f:
        pconfig = json.load(f)
    with open(DCONFIGFILE,"r") as f:
        dconfig = json.load(f)
    with open(BOTHCONFIGFILE,"r") as f:
        bconfig = json.load(f)
    
    #First, get the rootfile names for all data and MC
    datafiles = glob.glob("%s/data/*.root"%(DATADIR))
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
    SacAnalyzer.ShowSacrificePlot(evtype='prompt',dattype='data',fittotal=True)
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
    IThist,histMeta = SacAnalyzer.DrawCleanHist(evtype="pair",var='interevent_time',
                                                dattype="data",xmin=3000,
                                                xmax=999000.,nbins=50,
                                                addlROOTcuts="interevent_time>3000&&nhitsCleaned_d==14")
    #Neat.  Now, let's fit to this and plot it
    doubleexp = lambda x,A1,l1,A2,l2: A1*np.exp(-l1*x) + A2*np.exp(-l2*x)
    myfitter = sf.Fitter(datax=IThist["x"],datay=IThist["y"],
                      datasigma=IThist["y_unc"])
    myfitter.SetFitFunction(doubleexp, 4)
    initvars = [200., 1./20000., 15.,1./2.E7]
    popt, pcov = myfitter.RunFit(initvars)
    plt.errorbar(x=IThist["x"],y=IThist["y"],yerr=IThist["y_unc"],
                 linestyle='none',marker='o',markersize=5)
    bestfitline = doubleexp(IThist["x"],popt[0],popt[1],popt[2],popt[3])
    plt.plot(IThist["x"],bestfitline)
    plt.xlabel("Interevent Time (ns)")
    plt.ylabel("Events")
    plt.title("Data cleaned interevent time distribution of AmBe data")
    plt.show()
