#This script is a demonstration of how to use the AmBeDCSacrifice class
#for processing prompt/delayed candidate data and evaluating DC sacrifices
#for both data and MC files.

import json
import glob
import numpy as np
import os,sys
import plots.AmBeDCSacrifice as ab
import plots.EarlyLatePlots as el

DATADIR = "./rootfiles/all"
PCONFIGFILE = "./config/config_prompt_default.json"
DCONFIGFILE = "./config/config_delayed_default.json"
BOTHCONFIGFILE = "./config/early_pair_config.json"

PSYSCONFIGFILE = "./config/nhit_systematic/config_prompt_default.json"
DSYSCONFIGFILE = "./config/nhit_systematic/config_delayed_default.json"

if __name__=='__main__':
    print("LET US BEGIN")
    with open(PCONFIGFILE,"r") as f:
        pconfig = json.load(f)
    with open(DCONFIGFILE,"r") as f:
        dconfig = json.load(f)
    with open(BOTHCONFIGFILE,"r") as f:
        bconfig = json.load(f)
    
    #Different configuration that can be compared to default for systematic evaluation
    with open(PSYSCONFIGFILE,"r") as f:
        psysconfig = json.load(f)
    with open(DSYSCONFIGFILE,"r") as f:
        dsysconfig = json.load(f)
    
    #First, get the rootfile names for all data and MC
    datafiles = glob.glob("%s/data/*.root"%(DATADIR))
    mcfiles = glob.glob("%s/mc/*.root"%(DATADIR))
    #Now, start up the sacrifice analyzer
    SacAnalyzer = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig,bconfig)
    #Set the number of bins we want when plotting/analyzing
    SacAnalyzer.SetPromptBinNumber(14)
    SacAnalyzer.SetDelayedBinNumber(1)
    #Now, analyze the sacrifice per bin.  Start with nhits as the variable
    #accepts nhits, energy, udotr, or posr3 for nice plotting right now
    defaultdat, defaultmeta = SacAnalyzer.AnalyzeData(var='nhitsCleaned')
    #Now, lets try to plot out the prompt event's sacrifice for the variable analyze
    #We will make the plot for the MC files' prompt events
    SacAnalyzer.ShowSacrificePlot(evtype='prompt',dattype='data',fittotal=True)
    SacAnalyzer.ShowSacrificePlot(evtype='prompt',dattype='MC',fittotal=True)
    SacAnalyzer.ShowSacrificePlot(evtype='delayed',dattype='data',fittotal=True)
    SacAnalyzer.ShowSacrificePlot(evtype='delayed',dattype='MC',fittotal=False)
    #Lets do a data/MC comparison as well for the delayed events
    SacAnalyzer.DataMCSacCompare(evtype="prompt")
    #And let's look at the Data/MC sacrifice ratio
    SacAnalyzer.PlotDataMCSacRatio(evtype='prompt',fittotal=False)
    SacAnalyzer.PlotDataMCAccRatio(evtype='prompt',fittotal=False)
    SacAnalyzer.PlotDataMCAccRatio(evtype='delayed',fittotal=False)
