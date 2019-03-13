#This script is a demonstration of how to use the AmBeDCSacrifice class
#for processing prompt/delayed candidate data and evaluating DC sacrifices
#for both data and MC files.

import json
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys

import lib.playDarts as pd
import lib.AmBeDCSacrifice as ab
import lib.EarlyLatePlots as el
import lib.SimpleFitter as sf
import lib.NeutronSacCalculator as nsc
import lib.CorrelationTools as ct

sns.set_style("darkgrid")
sns.set_context("poster")

DATADIR = "./rootfiles/all"
PCONFIGFILE = "./config/config_prompt_default.json"
DCONFIGFILE = "./config/config_delayed_default.json"
BOTHCONFIGFILE = "./config/full_pair_config.json"

EARLYCONFIGFILE = "./config/early_pair_config.json"
LATECONFIGFILE = "./config/late_pair_config.json"

def GetRanges(dataframe,evtype,dattype):
    ranges = []
    xvalues = dataframe[evtype][dattype]["total"].vardat
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
    EarlySacAnalyze.SetPromptBinNumber(14)
    EarlySacAnalyze.SetDelayedBinNumber(11)
    LateSacAnalyze.SetPromptBinNumber(14)
    LateSacAnalyze.SetDelayedBinNumber(11)
    #Now, analyze the sacrifice per bin.  Start with nhits as the variable
    #accepts nhits, energy, udotr, or posr3 for nice plotting right now
    earlysacs, earlymeta = EarlySacAnalyze.AnalyzeData(var='nhitsCleaned',
                                                       onlycuts=['total'])
    latesacs, latemeta = LateSacAnalyze.AnalyzeData(var='nhitsCleaned',
                                                    onlycuts=['total'])

    SacAnalyzer = ab.AmBeSacrificeComparer(datafiles,mcfiles,pconfig,dconfig,bconfig)
    SacAnalyzer.SetPromptBinNumber(14)
    SacAnalyzer.SetDelayedBinNumber(11)
    #Now, analyze the sacrifice per bin.  Start with nhits as the variable
    #accepts nhits, energy, udotr, or posr3 for nice plotting right now
    defaultdat, defaultmeta = SacAnalyzer.AnalyzeData(var='nhitsCleaned',
                                                      onlycuts=['total'])
    esacs = {"prompt": [], "delayed": []}
    esac_uncs = {"prompt": [], "delayed": []}
    lsacs = {"prompt": [], "delayed": []}
    lsac_uncs = {"prompt": [], "delayed": []}
    calcsacs = {'prompt': [], 'delayed': []}
    calcsac_uncs = {'prompt':[], 'delayed': []}
    chartimes = {'prompt':[],'delayed':[]}
    chartime_uncs = {'prompt': [], 'delayed': []}
    pairsuff = {"prompt":"_p","delayed":"_d"}
    thevardat = {}
    nhitranges = {}
    binwidths = {}
    for key in pairsuff:
        nhitranges[key],binwidths[key]= GetRanges(defaultdat,key,"data")
        print("LENGTH OF NHIT RANGES: " + str(len(nhitranges)))
        for i in xrange(len(nhitranges[key])):
            earliest_time = econfig["interevent_time_low"]
            latest_time = lconfig["interevent_time_high"]
            lownhit = int(nhitranges[key][i][0])
            highnhit = int(nhitranges[key][i][1])
            print("NHIT RANGE MIN: " + str(lownhit))
            print("NHIT RANGE MAX: " + str(highnhit))
            #Get dirty interevent time distribution.  Fit estimates the fraction of 
            #correlated and uncorrelated events in the dataset
            IThist,histMeta = SacAnalyzer.DrawDirtyHist(evtype=key,var='interevent_time',
                                                        dattype="data",xmin=earliest_time,
                                                        xmax=latest_time,nbins=100,
                                                        addlROOTcuts="nhitsCleaned%s>=%i&&nhitsCleaned%s<%i"%(pairsuff[key],int(lownhit),pairsuff[key],int(highnhit)))
            #Neat.  Now, let's fit to this and plot it
            #doubleexp = lambda x,A1,l1,A2: A1*np.exp(-(l1)*x) + A2
            doubleexp = lambda x,A1,l1,A2,l2: A1*np.exp(-(l1)*x) + A2*np.exp(-(l2)*x) 

            varlabels = ["A1","l1","A2","l2"]
            singleexp = lambda x,A,l: A*np.exp(-l*x)
            myfitter = sf.Fitter(datax=IThist["x"],datay=IThist["y"],
                              datasigma=IThist["y_unc"])
            initvars = [2000., 1./20000., 100., 1E-7]
            myfitter.SetFitFunction(doubleexp, len(initvars))
            popt, pcov = myfitter.RunFit(initvars)
            try:
                if pcov == np.inf:
                    print("OH SHIIIIIIT FAILED FIT")
                    esacs[key].append(-1)
                    continue
            except ValueError:
                print("THE EFF WE GETTING A VALUE ERROR FOR")
                pass
            ErrorShooter = ct.CovarianceMatrix(cov_matrix=pcov,
                                               variables=varlabels)
            ErrorShooter.CholeskyDecompose()
            var_shots = ErrorShooter.ShootCorrVars(numsets=1000)
            for j,var in enumerate(varlabels):
                print("VAR %s unc. has mean shift of,stdev of: %s,%s"%(var,
                      str(np.mean(var_shots,axis=0)[j]),
                      str(np.std(var_shots,axis=0)[j])))
            charlamb = popt[1]
            charlamb_unc = np.sqrt(pcov[1][1])
            chartime = 1./charlamb
            chartime_unc = (1./(charlamb**2))*charlamb_unc
            chartimes[key].append(chartime)
            chartime_uncs[key].append(chartime_unc)
            #plt.errorbar(x=IThist["x"],y=IThist["y"],yerr=IThist["y_unc"],
            #             linestyle='none',marker='o',markersize=5)
            #bestfitline = doubleexp(IThist["x"],popt[0],popt[1],popt[2],popt[3])
            corry = singleexp(IThist["x"],popt[0],popt[1])
            uncorry = singleexp(IThist["x"],popt[2],popt[3])
            #plt.plot(IThist["x"],bestfitline,label="best fit")
            #plt.plot(IThist["x"],corry,label="Neutron capture fit")
            #plt.plot(IThist["x"],uncorry,label="Backgrounds fit")
            #plt.legend()
            #plt.xlabel("Interevent Time (ns)")
            #plt.ylabel("Events")
            #plt.title("Non-data cleaned interevent time distribution of AmBe data \n" + \
            #          "%s nhitsCleaned range: [%i,%i)"%(key,lownhit,highnhit))
            #plt.show()
            #Neato.  Now, load these into the data cleaning sacrifice class
            earlywindow = [econfig["interevent_time_low"],
                           econfig["interevent_time_high"]]
            latewindow = [lconfig["interevent_time_low"],
                          lconfig["interevent_time_high"]]
            esac = earlysacs[key]["data"]["total"].fractional_sacrifice[i]
            esacs[key].append(float(esac))
            lsac = latesacs[key]["data"]["total"].fractional_sacrifice[i]
            esac_unc = earlysacs[key]["data"]["total"].fs_uncertainty[i]
            esac_uncs[key].append(float(esac_unc))
            lsac_unc = latesacs[key]["data"]["total"].fs_uncertainty[i]
            thesets = 1000
            var_shots = ErrorShooter.ShootCorrVars(numsets=thesets)
            Ssacs = [] 
            for i in xrange(thesets):
                thisesac = pd.RandShoot(esac,esac_unc,1)
                thislsac = pd.RandShoot(lsac,lsac_unc,1)
                corry = singleexp(IThist["x"],popt[0]+var_shots[i][0],
                                  popt[1]+var_shots[i][1])
                uncorry = np.ones(len(IThist["x"])) * \
                          ((popt[2]+var_shots[i][2]))
                SCalc = nsc.NSacCalculator(IThist["x"],corry, uncorry, earlywindow,
                                            latewindow, thisesac, thislsac)
                Sac = SCalc.CalculateNeutronSacrifice()
                Ssacs.append(Sac)
            #plt.hist(np.array(NSsacs),bins=10)
            #plt.show()
            meanSsac = np.mean(Ssacs)
            stdevSsac = np.std(Ssacs)
            print("MEAN %s SACRIFICE AT THIS STEP: "%(key) + str(meanSsac))
            print("STD. DEV. %s SACRIFICE AT THIS STEP: "%(key) + str(stdevSsac))
            calcsacs[key].append(meanSsac)
            calcsac_uncs[key].append(stdevSsac)

        #Remove failed fits from the plot
        print("%s ESACS AFTER ALL LOOPS: "%(key) + str(esacs[key]))
        failedfits = np.where(np.array(esacs[key])==-1)[0]
        print("THE FAILED FITS: " + str(failedfits))
        thevardat[key] = np.array(earlysacs[key]["data"]["total"].vardat)
        thevardat[key] = np.delete(thevardat[key],failedfits)
        esacs[key] = np.delete(esacs[key],failedfits)

        #Make the delayed sacrifice/neutron sacrifice comparison vs. nhitsCleaned
        esacs[key] = np.array(esacs[key])
        esac_uncs[key] = np.array(esac_uncs[key])
        calcsacs[key] = np.array(calcsacs[key])
        calcsac_uncs[key] = np.array(calcsac_uncs[key])
        label='NO LABEL'
        if key=='prompt': label='4.4 MeV gamma'
        elif key=='delayed': label='neutron'
        plt.errorbar(x=thevardat[key], y=esacs[key], xerr=binwidths[key]/2.,
                     yerr=esac_uncs[key],linestyle='none',marker='o',markersize=7,
                     label='early window %s sacrifice'%(key))
        plt.errorbar(x=thevardat[key], y=calcsacs[key], xerr=binwidths[key]/2., yerr=calcsac_uncs[key],
                     linestyle='none',marker='o', markersize=7,label='calculated %s sacrifice'%(label))
        plt.xlabel("nhitsCleaned ")
        plt.ylabel("fractional DC sacrifice")
        plt.legend()
        plt.title("Comparison of %s candidate DC sacrifice in early window to \n "%(key)+\
                  "calculated %s DC sacrifice, AmBe internal data"%(label))
        plt.show()
    plt.errorbar(x=thevardat["prompt"], y=calcsacs["prompt"], xerr=binwidths["prompt"]/2.,
                 yerr=calcsac_uncs["prompt"],linestyle='none',marker='o',markersize=7,
                 label='Calculated 4.4 MeV gamma sacrifice')
    plt.errorbar(x=thevardat["delayed"], y=calcsacs["delayed"], xerr=binwidths["delayed"]/2., yerr=calcsac_uncs["delayed"],
                 linestyle='none',marker='o', markersize=7,label='calculated neutron sacrifice')
    plt.xlabel("nhitsCleaned ")
    plt.ylabel("fractional DC sacrifice")
    plt.legend()
    plt.title("Plot of calculated neutron and 4.4 MeV gamma sacrifices \n AmBe internal data")
    plt.show()
    
    #Make the plot showing how the fitted correlated time varies
    #plt.errorbar(x=thevardat[key], y=chartimes[key], xerr=binwidths[key]/2.,
    #             yerr=chartime_uncs[key],linestyle='none',marker='o',markersize=7,
    #             color='r',alpha=0.7)
    #plt.xlabel("nhitsCleaned")
    #plt.ylabel("Fitted neutron capture time (ns)")
    #plt.title("Variation in fitted neutron capture time across %s \n"%(key) +\
    #        "nhitsCleaned values")
    #plt.show()
    g_results = {}
    g_results['sacrifice'] = list(calcsacs['prompt'])
    g_results['sacrifice_uncertainty'] = list(calcsac_uncs['prompt'])
    thevardat['prompt'] = np.array(thevardat['prompt'])-0.5
    g_results['nhitsCleaned'] = list(thevardat['prompt'])
    with open("gamma_sacrifice.json","w") as f:
        json.dump(g_results,f,sort_keys=True,indent=4)
    p_results = {}
    p_results['sacrifice'] = list(esacs['prompt'])
    p_results['sacrifice_uncertainty'] = list(esac_uncs['prompt'])
    thevardat['prompt'] = np.array(thevardat['prompt'])-0.5
    p_results['nhitsCleaned'] = list(thevardat['prompt'])
    with open("prompt_sacrifice.json","w") as f:
        json.dump(p_results,f,sort_keys=True,indent=4)
    n_results = {}
    n_results['sacrifice'] = list(calcsacs['delayed'])
    n_results['sacrifice_uncertainty'] = list(calcsac_uncs['delayed'])
    thevardat['delayed'] = np.array(thevardat['delayed'])-0.5
    n_results['nhitsCleaned'] = list(thevardat['delayed'])
    with open("neutron_sacrifice.json","w") as f:
        json.dump(n_results,f,sort_keys=True,indent=4)
    d_results = {}
    d_results['sacrifice'] = list(esacs['delayed'])
    d_results['sacrifice_uncertainty'] = list(esac_uncs['delayed'])
    thevardat['delayed'] = np.array(thevardat['delayed'])-0.5
    d_results['nhitsCleaned'] = list(thevardat['delayed'])
    with open("delayed_sacrifice.json","w") as f:
        json.dump(d_results,f,sort_keys=True,indent=4)

