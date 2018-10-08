#This script is a demonstration of how to use the AmBeDCSacrifice class
#for processing prompt/delayed candidate data and evaluating DC sacrifices
#for both data and MC files.

import json
import glob
import numpy as np
import os,sys
import plots.AmBeDCSacrifice as ab

DATADIR = "./rootfiles"
PCONFIGFILE = "./config/config_prompt_default.json"
DCONFIGFILE = "./config/config_delayed_default.json"

if __name__=='__main__':
    print("LET US BEGIN")
    with open(PCONFIGFILE,"r") as f:
        pconfig = json.load(f)
    with open(DCONFIGFILE,"r") as f:
        dconfig = json.load(f)
    #First, get the rootfile names for all data and MC
    datafiles = glob.glob("%s/data/*.root"%(DATADIR))
    mcfiles = glob.glob("%s/mc/*.root"%(DATADIR))
    #Now, start up the sacrifice analyzer
    SacAnalyzer = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig)
    #Set the number of bins we want when plotting/analyzing
    SacAnalyzer.SetPromptBinNumber(4)
    SacAnalyzer.SetDelayedBinNumber(4)
    #Now, analyze the sacrifice per bin.  Start with nhits as the variable
    #accepts nhits, energy, udotr, or posr3 for nice plotting right now
    SacAnalyzer.AnalyzeData(var='nhits')
    #Now, lets try to plot out the prompt event's sacrifice for the variable analyze
    #We will make the plot for the MC files' prompt events
    SacAnalyzer.ShowSacrificePlot(evtype='prompt',dattype='data',fittotal=False)
    #Lets do a data/MC comparison as well for the delayed events
    SacAnalyzer.DataMCSacCompare(evtype="delayed")
    #And let's look at the Data/MC sacrifice ratio
    SacAnalyzer.PlotDataMCSacRatio(evtype='prompt',fittotal=False)
