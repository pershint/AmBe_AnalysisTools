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
    xvalues = dataframe["delayed"]["data"]["total"].vardat
    width = xvalues[1] - xvalues[0]
    for i in xrange(len(xvalues)):
        thisrange = [xvalues[i]- width/2.,xvalues[i]+width/2.]
        ranges.append(thisrange)
    return ranges,width

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
    earlysacs, earlymeta = EarlySacAnalyze.AnalyzeData(var='nhitsCleaned',
                                                       onlycuts=['total'])
    latesacs, latemeta = LateSacAnalyze.AnalyzeData(var='nhitsCleaned',
                                                    onlycuts=['total'])

    SacAnalyzer = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig,bconfig)
    SacAnalyzer.SetPromptBinNumber(6)
    SacAnalyzer.SetDelayedBinNumber(6)
    #Now, analyze the sacrifice per bin.  Start with nhits as the variable
    #accepts nhits, energy, udotr, or posr3 for nice plotting right now
    defaultdat, defaultmeta = SacAnalyzer.AnalyzeData(var='nhitsCleaned',
                                                      onlycuts=['total'])
    nhitranges,binwidth = GetRanges(defaultdat) #TODO: check values in vardat are center of bins
    print("LENGTH OF NHIT RANGES: " + str(len(nhitranges)))
    esac_vardat = earlysacs["delayed"]["data"]["total"].vardat
    print("THE EARLYSAC VARDAT LENGTH: ")
    print(len(np.array(esac_vardat)))
    #nhitranges = [[4,12]]
    esacs = []
    esac_uncs = []
    neutronsacs = []
    neutronsac_uncs = []
    chartimes = []
    chartime_uncs = []
    for i in xrange(len(nhitranges)):
        earliest_time = econfig["interevent_time_low"]
        latest_time = lconfig["interevent_time_high"]
        lownhit = int(nhitranges[i][0])
        highnhit = int(nhitranges[i][1])
        print("NHIT RANGE MIN: " + str(lownhit))
        print("NHIT RANGE MAX: " + str(highnhit))
        IThist,histMeta = SacAnalyzer.DrawCleanHist(evtype="pair",var='interevent_time',
                                                    dattype="data",xmin=earliest_time,
                                                    xmax=latest_time,nbins=100,
                                                    addlROOTcuts="nhitsCleaned_d>%i&&nhitsCleaned_d<=%i"%(int(lownhit),int(highnhit)))
        #Neat.  Now, let's fit to this and plot it
        doubleexp = lambda x,A1,l1,A2,l2: A1*np.exp(-l1*x) + A2*np.exp(-l2*x)
        singleexp = lambda x,A,l: A*np.exp(-l*x)
        myfitter = sf.Fitter(datax=IThist["x"],datay=IThist["y"],
                          datasigma=IThist["y_unc"])
        myfitter.SetFitFunction(doubleexp, 4)
        initvars = [2000., 1./20000., 15.,1./2.E8]
        popt, pcov = myfitter.RunFit(initvars)
        try:
            if pcov == np.inf:
                print("OH SHIIIIIIT FAILED FIT")
                esacs.append(-1)
                continue
        except ValueError:
            pass
        charlamb = popt[1]
        charlamb_unc = np.sqrt(pcov[1][1])
        chartime = 1./charlamb
        chartime_unc = (1./(charlamb**2))*charlamb_unc
        chartimes.append(chartime)
        chartime_uncs.append(chartime_unc)
        plt.errorbar(x=IThist["x"],y=IThist["y"],yerr=IThist["y_unc"],
                     linestyle='none',marker='o',markersize=5)
        bestfitline = doubleexp(IThist["x"],popt[0],popt[1],popt[2],popt[3])
        corry = singleexp(IThist["x"],popt[0],popt[1])
        uncorry = singleexp(IThist["x"],popt[2],popt[3])
        plt.plot(IThist["x"],bestfitline,label="best fit")
        plt.plot(IThist["x"],corry,label="Neutron capture fit")
        plt.plot(IThist["x"],uncorry,label="Backgrounds fit")
        plt.legend()
        plt.xlabel("Interevent Time (ns)")
        plt.ylabel("Events")
        plt.title("Data cleaned interevent time distribution of AmBe data \n" + \
                  "delayed nhitsCleaned range: (%i,%i]"%(lownhit,highnhit))
        plt.show()
        #Neato.  Now, load these into the data cleaning sacrifice class
        earlywindow = [econfig["interevent_time_low"],
                       econfig["interevent_time_high"]]
        latewindow = [lconfig["interevent_time_low"],
                      lconfig["interevent_time_high"]]
        esac = earlysacs["delayed"]["data"]["total"].fractional_sacrifice[i]
        esacs.append(float(esac))
        lsac = latesacs["delayed"]["data"]["total"].fractional_sacrifice[i]
        esac_unc = earlysacs["delayed"]["data"]["total"].fs_uncertainty[i]
        esac_uncs.append(float(esac_unc))
        lsac_unc = latesacs["delayed"]["data"]["total"].fs_uncertainty[i]
        NSCalc = nsc.NSacCalculator(IThist["x"],corry, uncorry, earlywindow,
                                    latewindow, esac,esac_unc, lsac,lsac_unc)
        NSac = NSCalc.CalculateNeutronSacrifice()
        print("DELAYED EARLY WINDOW SACRIFICE IS: %f"%(esac))
        print("DELAYED LATE WINDOW SACRIFICE IS: %f"%(lsac))
        print("NEUTRON SACRIFICE IS: %f"%(NSac))
        neutronsacs.append(NSac)
        neutronsac_uncs.append(0)

    #Remove failed fits from the plot
    failedfits = np.where(esacs==-1)[0]
    thevardat = np.array(earlysacs["delayed"]["data"]["total"].vardat)
    thevardat = np.delete(thevardat,failedfits)
    esacs = np.delete(esacs,failedfits)

    #Make the delayed sacrifice/neutron sacrifice comparison vs. nhitsCleaned
    esacs = np.array(esacs)
    esac_uncs = np.array(esac_uncs)
    neutronsacs = np.array(neutronsacs)
    print(esacs)
    print(esac_uncs)
    print(neutronsacs)
    plt.errorbar(x=thevardat, y=esacs, xerr=binwidth/2.,
                 yerr=esac_uncs,linestyle='none',marker='o',markersize=7,
                 label='early window delayed sacrifice')
    plt.errorbar(x=thevardat, y=neutronsacs, xerr=binwidth/2.,
                 linestyle='none',marker='o', markersize=7,label='calculated neutron sacrifice')
    plt.legend()
    plt.xlabel("nhitsCleaned ")
    plt.ylabel("fractional DC sacrifice")
    plt.title("Comparison of delayed candidate DC sacrifice in early window to \n "+\
              "calculated neutron DC sacrifice, AmBe internal data")
    plt.show()
    
    #Make the plot showing how the fitted correlated time varies
    plt.errorbar(x=thevardat, y=chartimes, xerr=binwidth/2.,
                 yerr=chartime_uncs,linestyle='none',marker='o',markersize=7,
                 color='r',alpha=0.7)
    plt.xlabel("nhitsCleaned")
    plt.ylabel("Fitted neutron capture time (ns)")
    plt.title("Variation in fitted neutron capture time across delayed \n" +\
            "nhitsCleaned values")
    plt.show()
