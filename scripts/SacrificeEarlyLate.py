#This script is a demonstration of how to use the AmBeDCSacrifice class
#for processing prompt/delayed candidate data and evaluating DC sacrifices
#for both data and MC files.

import json
import glob
import numpy as np
import os,sys
import lib.AmBeDCSacrifice as ab
import lib.EarlyLatePlots as el

DATADIR = "./rootfiles/all"
PCONFIGFILE = "./config/config_prompt_default.json"
DCONFIGFILE = "./config/config_delayed_default.json"
EARLYCONFIGFILE = "./config/early_pair_config.json"
LATECONFIGFILE = "./config/late_pair_config.json"

if __name__=='__main__':
    print("LET US BEGIN")
    with open(PCONFIGFILE,"r") as f:
        pconfig = json.load(f)
    with open(DCONFIGFILE,"r") as f:
        dconfig = json.load(f)
    with open(EARLYCONFIGFILE,"r") as f:
        econfig = json.load(f)
    with open(LATECONFIGFILE,"r") as f:
        lconfig = json.load(f)
    #First, get the rootfile names for all data and MC
    datafiles = glob.glob("%s/data/*.root"%(DATADIR))
    mcfiles = glob.glob("%s/mc/*.root"%(DATADIR))
    #Now, start up the sacrifice analyzer
    EarlySacAnalyze = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig,econfig)
    #Set the number of bins we want when plotting/analyzing
    EarlySacAnalyze.SetPromptBinNumber(6)
    EarlySacAnalyze.SetDelayedBinNumber(6)

    LateSacAnalyze = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig,lconfig)
    #Set the number of bins we want when plotting/analyzing
    LateSacAnalyze.SetPromptBinNumber(6)
    LateSacAnalyze.SetDelayedBinNumber(6)

    #Now, analyze the sacrifice per bin.  Start with nhits as the variable
    #accepts nhits, energy, udotr, or posr3 for nice plotting right now
    earlysacs, earlymeta = EarlySacAnalyze.AnalyzeData(var='nhitsCleaned')
    latesacs, latemeta = LateSacAnalyze.AnalyzeData(var='nhitsCleaned')
    #Now, lets try to plot out the prompt event's sacrifice for the variable analyze
    #We will make the plot for the MC files' prompt events
    EarlySacAnalyze.ShowSacrificePlot(evtype='delayed',dattype='data',fittotal=True)
    LateSacAnalyze.ShowSacrificePlot(evtype='delayed',dattype='data',fittotal=True)

    #Let's try the new tool; this is a general sacrifice comparer that can
    #compare two data cleaning sacrifice dataframes
    elplotter = el.SacComparer(earlysacs, latesacs, earlymeta, latemeta)
    elplotter.PlotSacrifices(fitdiff=True,title=None,dattype='data', 
                             evtype='delayed',label1=r"interevent_time<500 $\mu$s",
                             label2=r"interevent_time>500 $\mu$s")
    elplotter.PlotSacrifices(fitdiff=True,title=None,dattype='data',
                             evtype='prompt',label1=r"interevent_time<500 $\mu$s",
                             label2=r"interevent_time>500 $\mu$s")
    elplotter.PlotDataSacWeightedDifference(fitdiff=True,title=None,evtype='delayed')
