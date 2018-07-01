#Some functions are included here for grouping data from AmBe IBD root files as
#Prepared by find_prompt_delayed.py in this same directory.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ROOT
import maskbuilder as mb

def SacVSTimeWindow_Data(rootfiles,event="prompt",nmin=10.0,nmax=100.0,tmin=0.0, tmax=800000.0,tsteps=40000.0, fitvalid=False,dcmask=None):
    '''
    Takes in a rootfile and returns a PandaFrame object that can be used
    for plotting in matplotlib.  Returns nhits vs. fractional sacrifice for
    Each cut in the given dcmask.
    '''
    if event == "prompt":
        ev = "_p"
    elif event == "delayed":
        ev = "_d"
    else:
        print("Choose either the prompt or delayed event to consider.")
        return

    if dcmask == None:
        print("No DC mask given.  Please give a DC mask of cuts to plot")
        return

    data = ROOT.TChain("CombinedOutput")
    for f in rootfiles:
        data.Add(f)
    fullmask = mb.get_dcwords()
    plotmask = {}
    for cut in fullmask:
        if 'prescale' in fullmask[cut]:
            continue
        if 'waterblind' in fullmask[cut]:
            continue
        if mb.binary_bit_to_int(cut)&dcmask==mb.binary_bit_to_int(cut):
            print("ADDING IN..." + str(fullmask[cut]))
            plotmask[cut] = fullmask[cut]
    #Get the total cut mask
    fullmask = 0
    for cut in plotmask:
        cutint = mb.binary_bit_to_int(cut)
        fullmask+=int(cutint)
    #Now, make a new dictionary object: key is cutname, value is histogram
    nbins = nmin-nmax
    twindow, fs,fs_unc = (), (), () #pandas wants ntuples
    timearr = np.arange(tmin,tmax,tsteps)
    base = ""
    if fitvalid is True:
        base = base + "fitValid"+ev + "==1" + "&&"
    for i,time in enumerate(timearr):
        if i==0: continue
        twindow = twindow + (time,)
        ntot = float(data.GetEntries(base+"nhits%s>%i && nhits%s<%i && interevent_time<%d && interevent_time >%d" % (ev,nmin,ev,nmax,timearr[i],timearr[i-1])))
        nflag = float(data.GetEntries(base+"nhits%s>%i && nhits%s<%i && interevent_time<%d && interevent_time >%d && ((dcFlagged%s&%i)!=%i)" % (ev,nmin,ev,nmax,timearr[i],timearr[i-1],ev,fullmask,fullmask)))
        fs = fs + (nflag/ntot,)
        fs_unc = fs_unc + (np.sqrt((np.sqrt(nflag)/ntot)**2 +
            (nflag*np.sqrt(ntot)/(ntot**2))**2),)
  
    print(twindow,fs,fs_unc)
    dat = {"twindow_edge": twindow, "fractional_sacrifice":fs, 
            "fs_uncertainty":fs_unc}
    timedf = pd.DataFrame(data=dat)
    timedf.cutmask = dcmask
    timedf.timewidth = tsteps
    timedf.nmin = nmin
    timedf.nmax = nmax
    timedf.evtype = event
    data.GetEntry(0)
    return timedf

def SacVSTimeWindow_Plot(timedat):
    sns.set_style("whitegrid")
    xkcd_color=['slate blue']
    sns.set_palette(sns.xkcd_palette(xkcd_color))
    timewindow = timedat.twindow_edge[1] - timedat.twindow_edge[0]
    plt.errorbar(x=timedat.twindow_edge, y=timedat.fractional_sacrifice,
            yerr=timedat.fs_uncertainty, xerr=timewindow/2, linestyle='none',
            marker='.', markersize=10,capsize=0,elinewidth=2)
    plt.legend(loc=1)
    plt.ylabel("Fractional sacrifice", fontsize=22)
    plt.xlabel("Time window", fontsize=22)
    plt.tick_params(labelsize=20)
    plt.title("Fractional sacrifice of %s events due to data cleaning cuts\n" % (timedat.evtype)+\
            "mask %i used on nhit range [%i,%i]" % (timedat.cutmask,timedat.nmin,timedat.nmax),
            fontsize=24)
    plt.ion()
    plt.show()
