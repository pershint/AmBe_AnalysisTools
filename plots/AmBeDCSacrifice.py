from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.optimize as spc
import pandas
import time


class AmBeSacrificeComparer(object):
    def __init__(self,rootfiles_data,rootfiles_mc,cdict_prompt,cdict_delayed):
        self.cdict = {'prompt':cdict_prompt,'delayed':cdict_delayed}
        self.top_cuts = {'prompt': {'data': [], 'MC': []}, 'delayed': {'data':[], 'MC':[]}}
        self.rootfiles_data = rootfiles_data
        self.rootfiles_mc = rootfiles_mc
        self.data_precuts = None
        self.mc_precuts = None
        self.nbins = 11
        self.evtype = evtype
        self.analyze_range = {}
        self._SetAnalyzeRanges()

    def _SetAnalyzeRanges(self):
        '''Uses variable ranges defined in cut dictionary and sets the analysis range
        to these windows'''
        self.analyze_range['prompt'] = {}
        self.analyze_range['delayed'] = {}
        for key in self.analyze_range:
            self.analyze_range[key]['energy'] = [self.cdict[key]['E_low'],self.cdict[key]['E_high']]
            self.analyze_range[key]['udotr'] = [self.cdict[key]['udotr_low'],self.cdict[key]['udotr_high']]
            self.analyze_range[key]['beta14'] = [self.cdict[key]['cut2_b14_low'],self.cdict[key]['cut2_b14_high']]
            self.analyze_range[key]['itr'] = [self.cdict[key]['cut2_itr_low'],1.0]
            self.analyze_range[key]['posr'] = [self.cdict[key]['r_low'],self.cdict[key]['r_high']] 
            self.analyze_range[key]['posr3'] = [self.cdict[key]['r_low'],self.cdict[key]['r_high']] 
            self.analyze_range[key]['nhits'] = [self.cdict[key]['nhits_low'],self.cdict[key]['nhits_high']] 

    def SetEventType(self,evtype):
        '''Sets whether to evaluate the DC sacrifice for the prompt or delayed events in
        the given IBD candidates'''
        self.evtype=evtype

    def SetBinNumber(self, nbins):
        self.nbins = nbins + 1
    
    def _GetAnalysisMask(self):
        util = RAT.DU.Utility.Get()
        util.LoadDBAndBeginRun()
        dcmask = RAT.GetDataCleaningWord('analysis_mask')
        return dcmask
    
    def _flatline(self,x,b):
        return b
    
    def _weighted_stdev(self,vals,valavg, uncs):
        weights = 1/(uncs**2)
        return np.sqrt((len(weights)*np.sum(weights*(vals-valavg)**2))/((len(weights)-1)*np.sum(weights)))

    def _DefinePrecuts(self):
        '''Defines the preliminary cuts for the prompt and delayed events for both
        data and MC files'''
        self.precuts_data = {}
        self.precuts_mc={}
        for key in self.cdict:
            if key == "prompt":
                suf = "_p"
            elif key == "delayed":
                suf = "_d"
            if self.cdict[key]['r_high'] is not None:
                self.precuts_data[key].append("posr%s<%s"%(suf,str(self.cdict[key]['r_high'])))
                self.precuts_mc[key].append("posr%s<%s"%(suf,str(self.cdict[key]['r_high'])))
            if self.cdict[key]['r_low'] is not None:
                self.precuts_data[key].append("posr%s>%s"%(suf,str(self.cdict[key]['r_low']))
                self.precuts_mc[key].append("posr>"+str(self.cdict[key]['r_low']))
            if self.cdict[key]['nhits_high'] is not None:
                self.precuts_data[key].append('nhits%s<%s'%(suf,str(self.cdict[key]['nhits_high'])))
                self.precuts_mc[key].append('nhits%s<%s'%(suf,str(self.cdict[key]['nhits_high'])))
            if self.cdict[key]['nhits_low'] is not None:
                self.precuts_data[key].append('nhits%s>%s'%(suf,str(self.cdict[key]['nhits_low'])))
                self.precuts_mc[key].append('nhits%s>%s'%(suf,str(self.cdict[key]['nhits_low'])))
            if self.cdict[key]['udotr_high'] is not None:
                self.precuts_data[key].append('udotr%s<%s'%(suf,str(self.cdict[key]['udotr_high'])))
                self.precuts_mc[key].append('udotr%s<%s'%(suf,str(self.cdict[key]['udotr_high'])))
            if self.cdict[key]['udotr_low'] is not None:
                self.precuts_data[key].append('udotr%s>%s'%(suf,str(self.cdict[key]['udotr_low'])))
                self.precuts_mc[key].append('udotr%s>%s'%(suf,str(self.cdict[key]['udotr_low'])))
            if self.cdict[key]['E_high'] is not None:
                self.precuts_data[key].append('energy%s<%s'%(suf,str(self.cdict[key]['E_high'])))
                self.precuts_mc[key].append('energy%s<%s'%(suf,str(self.cdict[key]['E_high'])))
            if self.cdict[key]['E_low'] is not None:
                self.precuts_data[key].append('energy%s>%s'%(suf,str(self.cdict[key]['E_low'])))
                self.precuts_mc[key].append('energy%s>%s'%(suf,str(self.cdict[key]['E_low'])))
            if self.cdict[key]['Z_low'] is not None:
                self.precuts_data[key].append('posz%s>%s'%(suf,str(self.cdict[key]['Z_low']*10.0)))
                self.precuts_mc[key].append('posz%s>%s'%(suf,str(self.cdict[key]['Z_low']*10.0)))
            if self.cdict[key]['Z_high'] is not None:
                self.precuts_data[key].append('posz%s<%s'%(suf,str(self.cdict[key]['Z_high']*10.0)))
                self.precuts_mc[key].append('posz%s<%s'%(suf,str(self.cdict[key]['Z_high']*10.0)))
            if self.cdict[key]['fitValid'] is not None: 
                if self.cdict[key]['fitValid'] is True: 
                    self.precuts_data[key].append("fitValid%s==1"%(suf))
                    self.precuts_mc[key].append("fitValid%s==1"%(suf))
                else:
                    self.precuts_data[key].append("fitValid%s==0"%(suf))
                    self.precuts_mc[key].append("fitValid%s==0"%(suf))
            if self.cdict[key]['AVudotrCut'] is not None:
                if self.cdict[key]['AVudotrCut'] is True:
                    rcut = self.cdict[key]["AVudotrCut_rcut"]
                    cutparam = (rcut/6000.0)**3
                    self.precuts_data[key].append("(udotr%s > (1.0 - 12 *((posr3%s-%f)**2)))"%(suf,suf,cutparam))
                    self.precuts_mc[key].append("(udotr%s > (1.0 - 12 *((posr3%s-%f)**2)))"%(suf,suf,cutparam))
            if self.cdict[key]["sacpath_DCmask"] is not None:
                self.precuts_data[key].append("((dcFlagged%s&%s)==%s)" % (suf,self.cdict[key]["sacpath_DCmask"],\
                    self.cdict[key]["sacpath_DCmask"]))
                self.precuts_mc[key].append("((dcFlagged%s&%s)==%s)" % (suf,self.cdict[key]["sacpath_DCmask"],\
                    self.cdict[key]["sacpath_DCmask"]))
            
            self.precuts_data[key].append("((triggerWord&%s)==0)" % (self.cdict[key]["path_trigmask"]))
            self.precuts_data[key] = self._cutfuse(self.precuts_data[key], "&&")
            self.precuts_mc[key] = self._cutfuse(self.precuts_mc[key], "&&")

    def _cutfuse(self,stringlist,delim):
        outstring = ""
        for j,s in enumerate(stringlist):
            if j == 0:
                outstring = s
            else:
                outstring=outstring+delim+s
        return outstring

    def _GetTopSacs(self,topnumber=7,dattype=None):
        if dattype=None:
            print("Must input a dattype that is either 'data' or 'MC'")
            return None
        for key in self.cdict:
            dcmask = self.cdict[key]['cut1_sacDCmask']
            if len(self.sac_percut[key][dattype]) <= topnumber:
                print("You already only have five cuts.  Not combining any")
                return
            names = []
            sacrifices = []
            for cut in self.sac_percut[key][dattype]:
                sacrifices.append(np.average(self.sac_percut[key][dattype][cut]["fractional_sacrifice"]))
                names.append(cut)
            sortednames = [x for _,x in sorted(zip(sacrifices,names))]
            topnames = sortednames[(len(sortednames)-(topnumber+1)):len(sortednames)]
            self.top_cuts[key][dattype] = topnames
    
    def AnalyzeData(self, var="nhits",xmin=None,xmax=None):

        '''Analyzes the DC sacrifice for the data and Monte Carlo data fed to
        the class.  Returns the sacrifice as a function of the variable given'''
        havebounds = False for v in self.analyze_range:
            if (xmin is None or xmax is None) and v == var: 
                xmin,xmax = self.analyze_range[var][0],self.analyze_range[var][1]
                havebounds = True
        if not havebounds:
            if xmin is None or xmax is None:
                print("You do not have bounds defined for analyzing this variable"+\
                        "defined.  Check your config file.")
                sys.exit(0)
        varname = var
        if var=='posr3':
            xmin = (float(xmin)/6000.0)**3
            xmax = (float(xmax)/6000.0)**3
        data = ROOT.TChain("output")
        MC = ROOT.TChain("output") 
        for rf in self.rootfiles_data:
            data.Add(rf)
        for mf in self.rootfiles_mc:
            MC.Add(mf)

        for key in self.cdict: 
            if key=="prompt":
                suf = "_p"
            elif key=="delayed":
                suf = "_d"
            dcmask = self.cdict[key]['cut1_sacDCmask']
            if dcmask == None:
                print("No DC mask in cuts config.  Please give a DC mask of cuts to plot")
                sys.exit(0)
            fullmask = mb.get_dcwords()
            print(fullmask)
            plotmask = {}
            for cut in fullmask:
                if 'prescale' in fullmask[cut]:
                    continue
                if 'waterblind' in fullmask[cut]:
                    continue
                if mb.binary_bit_to_int(cut)&dcmask==mb.binary_bit_to_int(cut):
                    print("WE ARE ADDING IN CUT: " + str(cut))
                    print("THAT IS... CUT: " + str(fullmask[cut]))
                    plotmask[cut] = fullmask[cut]
            plotmask[dcmask] = "total"
                #Now, make a new dictionary object: key is cutname, value is histogram
            data.Draw("%s>>h_dallevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "%s" % (self.precuts_data[key]),"goff")
            MC.Draw("%s>>h_mallevents(%i,%f,%f)"% (var,self.nbins,xmin,xmax),
                    "%s" % (self.precuts_mc[key]),"goff")
            h_dallevents = gDirectory.Get("h_dallevents")
            h_dallevents.Sumw2()
            h_mallevents = gDirectory.Get("h_mallevents")
            h_mallevents.Sumw2()
                
            cutnames = ()
            self.sac_percut[key] = {}
            for cut in plotmask:
                if cut!= dcmask:
                    cutint = mb.binary_bit_to_int(cut)
                else:
                    cutint = cut
                h_cut_dFracFlagged = ROOT.TH1D("h_cut_dFracFlagged", "h_cut_dFracFlagged", self.nbins,xmin,xmax)
                h_cut_dFracFlagged.Sumw2()
                h_cut_mFracFlagged = ROOT.TH1D("h_cut_mFracFlagged", "h_cut_mFracFlagged", self.nbins,xmin,xmax)
                h_cut_mFracFlagged.Sumw2()
                if self.precuts_data[key] is not None:
                    data.Draw("%s>>h_dflagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                            "%s"%(self.precuts_data[key])+\
                                 "&&((dcFlagged%s&%i)!=%i)" % (suf,cutint,cutint),"goff")
                else:
                    data.Draw("%s>>h_dflagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                            "((dcFlagged%s&%i)!=%i)" % (suf,cutint,cutint),"goff")
                if self.precuts_mc[key] is not None:
                    MC.Draw("%s>>h_mflagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                            "%s"%(self.precuts_mc[key])+\
                                 "&&((dcFlagged%s&%i)!=%i)" % (suf,cutint,cutint),"goff")
                else:
                    MC.Draw("%s>>h_mflagged(%i,%f,%f)" % (var,self.nbins,xmin,xmax),
                            "((dcFlagged%s&%i)!=%i)" % (suf,cutint,cutint),"goff")
                h_dflagged = gDirectory.Get("h_dflagged")
                h_dflagged.Sumw2()
                h_mflagged = gDirectory.Get("h_mflagged")
                h_mflagged.Sumw2()
                h_cut_dFracFlagged.Divide(h_dflagged,h_dallevents,1.,1.,"b")
                h_cut_mFracFlagged.Divide(h_mflagged,h_mallevents,1.,1.,"b")
                data_sacdict = self._histToDict(h_cut_dFracFlagged)
                mc_sacdict = self._histToDict(h_cut_mFracFlagged)
                self.sac_percut[key]["data"][plotmask[cut]]= pandas.DataFrame(data=data_sacdict)
                self.sac_percut[key]["MC"][plotmask[cut]]= pandas.DataFrame(data=mc_sacdict)
                del h_cut_dFracFlagged
                del h_dflagged
                del h_cut_mFracFlagged
                del h_mflagged
            #graphdict has the sacrifice information for each cut. Now, let's plot it.
            self.sac_percut_metadata = {"binwidth":((xmax-xmin)/float(self.nbins)),"variable":varname}
            self._GetTopSacs()
            self._DeleteEmptyBins()
        print self.sac_percut

    def _histToDict(self,roothist):
        graphdict={}
        vardat, fs,fs_unc =(), (), () #pandas wants ntuples
        for i in xrange(int(roothist.GetNbinsX()+1)):
            if i==0:
                continue
            vardat =  vardat + ((float(roothist.GetBinWidth(i))/2.0) + float(roothist.GetBinLowEdge(i)),)
            fs = fs + (roothist.GetBinContent(i),)
            fs_unc = fs_unc + (roothist.GetBinError(i),)
        graphdict["vardat"] = vardat
        graphdict["fractional_sacrifice"] = fs
        graphdict["fs_uncertainty"] = fs_unc
        return graphdict

def SacVSVar_Data(rootfiles=[],evtype = "prompt", precuts=None,var="nhits",nbins=10,xmin=10.0,xmax=100.0,dcmask=None):
    '''
    Takes in a rootfile and returns a PandaFrame object that can be used
    for plotting in matplotlib.  Returns nhits vs. fractional sacrifice for
    Each cut in the given dcmask.
    '''
    if evtype=="prompt":
        vt = "p"
    elif evtype=="delayed":
        vt = "d"
    else:
        print("Choose either prompt or delayed variable type...")
        return
    if dcmask == None:
        print("No DC mask given.  Please give a DC mask of cuts to plot")
        return
    data = ROOT.TChain("CombinedOutput")
    for rf in rootfiles:
        data.Add(rf)
    fullmask = mb.get_dcwords()
    print(fullmask)
    plotmask = {}
    for cut in fullmask:
        if 'prescale' in fullmask[cut]:
            continue
        if 'waterblind' in fullmask[cut]:
            continue
        if mb.binary_bit_to_int(cut)&dcmask==mb.binary_bit_to_int(cut):
            print("WE ARE ADDING IN CUT: " + str(cut))
            print("THAT IS... CUT: " + str(fullmask[cut]))
            plotmask[cut] = fullmask[cut]
    plotmask[dcmask] = "total"
    #Now, make a new dictionary object: key is cutname, value is histogram
    if precuts is not None: 
        data.Draw("%s_%s>>h_allevents(%i,%f,%f)"% (var,vt,nbins,xmin,xmax),
                precuts,"goff")
    else: 
        data.Draw("%s_%s>>h_allevents(%i,%f,%f)"% (var,vt,nbins,xmin,xmax),"","goff")
    h_allevents = gDirectory.Get("h_allevents")
    h_allevents.Sumw2()
    cutnames = ()
    allcutsacs = {}
    for cut in plotmask:
        if cut!= dcmask:
            cutint = mb.binary_bit_to_int(cut)
        else:
            cutint = cut
        graphdict={}
        vardat, fs,fs_unc =(), (), () #pandas wants ntuples
        h_cut_FracFlagged = ROOT.TH1D("h_cut_FracFlagged", "h_cut_FracFlagged", nbins,xmin,xmax)
        h_cut_FracFlagged.Sumw2()
        if precuts is not None:
            data.Draw("%s_%s>>h_flagged(%i,%f,%f)" % (var,vt,nbins,xmin,xmax),
                    "%s&&((dcFlagged_%s&%i)!=%i)" % (precuts,vt,cutint,cutint),"goff")
        else:
            data.Draw("%s_%s>>h_flagged(%i,%f,%f)" % (var,vt,nbins,xmin,xmax),
                    "((dcFlagged_%s&%i)!=%i)" % (vt,cutint,cutint),"goff")
        h_flagged = gDirectory.Get("h_flagged")
        h_flagged.Sumw2()
        h_cut_FracFlagged.Divide(h_flagged,h_allevents,1.,1.,"b")
        for i in xrange(h_cut_FracFlagged.GetNbinsX()):
            vardat =  vardat + (float(h_cut_FracFlagged.GetBinWidth(i)) + float(h_cut_FracFlagged.GetBinLowEdge(i)),)
            fs = fs + (h_cut_FracFlagged.GetBinContent(i),)
            fs_unc = fs_unc + (h_cut_FracFlagged.GetBinError(i),)
        graphdict["vardat"] = vardat
        graphdict["fractional_sacrifice"] = fs
        graphdict["fs_uncertainty"] = fs_unc
        allcutsacs[plotmask[cut]]= pandas.DataFrame(data=graphdict)
        del h_cut_FracFlagged
        del h_flagged
    #graphdict has the sacrifice information for each cut. Now, let's plot it.
    meta = {"binwidth":((xmax-xmin)/float(nbins)),"variable":"%s_%s"%(var,vt),"precuts":precuts}
    return allcutsacs, meta

    #TODO: We need to bring this into the "multiple datatypes" and "prompt or delayed"
    #Future
    def ShowPlottedSacrifice(self,fittotal=True,title=None,savedir="."):
        sns.set_style("whitegrid")
        xkcd_colors = ['black','slate blue', 'fluro green', 'brown', 'blue',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'leaf',
                'aqua blue','vomit', 'red','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        fig.set_size_inches(18.5,11.5)
        ax = fig.add_subplot(1,1,1)
        if self.cut1_mask['dcmask_cutnames'] is None:
            for cut in self.sac_percut:
                plt.errorbar(x=self.sac_percut[cut].vardat, 
                        y=self.sac_percut[cut].fractional_sacrifice,
                        yerr=self.sac_percut[cut].fs_uncertainty,
                        linestyle='none', marker='+', label=cut, markersize=7,
                        elinewidth=2, capsize=0)
        else:
            fracsum,fracuncsum=[],[]
            for cut in self.sac_percut:
                if cut in self.cut1_mask['dcmask_cutnames']:
                    print(cut)
                    plt.errorbar(x=self.sac_percut[cut].vardat, y=self.sac_percut[cut].fractional_sacrifice,
                            yerr=self.sac_percut[cut].fs_uncertainty, linestyle='none',
                            marker='o', label=cut, capsize=0, elinewidth=2, markersize=8)
                else:
                    print("CUT GOING INTO ALL OTHER CUTS: " + str(cut))
                    fracsum.append(self.sac_percut[cut].fractional_sacrifice)
                    print("FRACSUM ARRAY: " + str(fracsum))
                    fracuncsum.append(self.sac_percut[cut].fs_uncertainty)
            plt.errorbar(x=self.sac_percut[cut].vardat,y=sum(fracsum),yerr=sum(fracuncsum),
                    linestyle='none',marker='o', capsize=0, elinewidth=2, label='All other cuts', markersize=8)
    
        if fittotal is True:
            popt, pcov = spc.curve_fit(self._flatline, self.sac_percut["total"].vardat, self.sac_percut["total"].fractional_sacrifice,
                    p0=[0.02], sigma=self.sac_percut["total"].fs_uncertainty)
            #one standard deviation
            try:
                stdev = np.sqrt(np.diag(pcov))
            except ValueError:
                print("Fit likely failed.  You may not have enough data"+\
                        "in this region to fit a straight line.")
                raise
            plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
        legend = plt.legend(loc=3,frameon=1)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        plt.yscale("log")
        plt.ylabel("Fractional sacrifice",fontsize=32)
      
        varindict = False
        for var in self.xlabel_dict:
            if self.sac_percut_metadata["variable"] == var:
                plt.xlabel(self.xlabel_dict[var],fontsize=32)
                varindict = True
            else:
                if varindict is False: 
                    plt.xlabel(self.sac_percut_metadata["variable"],fontsize=32)
        plt.tick_params(labelsize=30)
        if title is None:
            plt.title("Fractional sacrifice due to data cleaning cuts",fontsize=36)
        else: 
            plt.title(title,fontsize=36)
        variable = self.sac_percut_metadata["variable"]
        plt.savefig(savedir+"/DCSac_%s.pdf"%(variable))
        #plt.show()
        plt.close()




def SacVSVar_Plot(allcutsacs,topcuts=None,metadata={"binwidth":10.0,"variable":"nhits"},fittotal=True):
    sns.set_style("whitegrid")
    xkcd_colors = ['black', 'slate blue',  'fluro green', 'brown', 'blue',
            'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'leaf',
            'aqua blue','grey', 'red','twilight']
    sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(allcutsacs)))
    if topcuts is None:
        for cut in allcutsacs:
            plt.errorbar(x=allcutsacs[cut].vardat, 
                    y=allcutsacs[cut].fractional_sacrifice,
                    yerr=allcutsacs[cut].fs_uncertainty,
                    linestyle='none', marker='o', label=cut, markersize=8,
                    elinewidth=2, capsize=0)
    else:
        fracsum,fracuncsum=[],[]
        for cut in allcutsacs:
            if cut in topcuts:
                plt.errorbar(x=allcutsacs[cut].vardat, y=allcutsacs[cut].fractional_sacrifice,
                        yerr=allcutsacs[cut].fs_uncertainty, linestyle='none',
                        marker='o', label=cut, capsize=0, elinewidth=2, markersize=8)
            else:
                fracsum.append(allcutsacs[cut].fractional_sacrifice)
                fracuncsum.append(allcutsacs[cut].fs_uncertainty)
        plt.errorbar(x=allcutsacs[cut].vardat,y=sum(fracsum),yerr=sum(fracuncsum),
                linestyle='none',marker='o', capsize=0, elinewidth=2, label='All other cuts', markersize=8)

    if fittotal is True:
        popt, pcov = spc.curve_fit(flatline, allcutsacs["total"].vardat, allcutsacs["total"].fractional_sacrifice,
                p0=[0.02], sigma=allcutsacs["total"].fs_uncertainty)
        #one standard deviation
        print("BEST FIT: " + str(popt))
        print("PCOVARIANCE: " + str(pcov))
        stdev = np.sqrt(np.diag(pcov))
        plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
        print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(_weighted_stdev(allcutsacs["total"].fractional_sacrifice,
            float(popt[0]),allcutsacs["total"].fs_uncertainty)))
    legend = plt.legend(loc=3,frameon=1)
    frame = legend.get_frame()
    frame.set_facecolor("white")
    plt.yscale("log")
    plt.ylabel("Fractional sacrifice",fontsize=32)
    plt.xlabel(metadata["variable"],fontsize=32)
    plt.tick_params(labelsize=30)
    plt.title("Fractional sacrifice of dataset by data cleaning cuts",fontsize=34)
    plt.ion()
    plt.show()
