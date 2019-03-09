from ROOT import gDirectory
import ROOT
from rat import RAT
import numpy as np
import maskbuilder as mb
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as spc
import pandas
import time


class AmBeSacrificeComparer(object):
    def __init__(self,rootfiles_data,rootfiles_mc,cdict_prompt,cdict_delayed,cdict_pair):
        self.cdict = {'prompt':cdict_prompt,'delayed':cdict_delayed}
        self.cdict_pair = cdict_pair
        self.top_cuts = {'prompt': {'data': [], 'MC': []}, 'delayed': {'data':[], 'MC':[]}}
        self.rootfiles_data = rootfiles_data
        self.rootfiles_mc = rootfiles_mc
        print(self.rootfiles_data)
        self.precuts_data = {'prompt':[], 'delayed':[]}
        self.precuts_mc = {'prompt':[], 'delayed':[]}
        self.nbins = {'prompt': 11, 'delayed': 11}
        self.sac_percut = {}
        self.analyze_range = {}
        self._SetAnalyzeRanges()
        self._DefinePrecuts()
        
        self.xlabel_dict = {'energy': 'Energy (MeV)', 'udotr': r'$U . R$',
                'posr3': r'$(R/R_{AV})^{3}$','posr':'Radius (mm)'}
    
    def _SetAnalyzeRanges(self):
        '''Uses variable ranges defined in cut dictionary and sets the analysis range
        to these windows'''
        self.analyze_range['prompt'] = {}
        self.analyze_range['delayed'] = {}
        for key in self.analyze_range:
            self.analyze_range[key]['energy'] = [self.cdict[key]['E_low'],self.cdict[key]['E_high']]
            self.analyze_range[key]['udotr'] = [self.cdict[key]['udotr_low'],self.cdict[key]['udotr_high']]
            self.analyze_range[key]['beta14'] = [self.cdict[key]['b14_low'],self.cdict[key]['b14_high']]
            self.analyze_range[key]['itr'] = [self.cdict[key]['itr_low'],self.cdict[key]['itr_high']]
            self.analyze_range[key]['posr'] = [self.cdict[key]['r_low'],self.cdict[key]['r_high']] 
            self.analyze_range[key]['posr3'] = [self.cdict[key]['r_low'],self.cdict[key]['r_high']] 
            self.analyze_range[key]['nhits'] = [self.cdict[key]['nhits_low'],self.cdict[key]['nhits_high']]
            self.analyze_range[key]['nhitsCleaned'] = [self.cdict[key]['nhitsCleaned_low'],self.cdict[key]['nhitsCleaned_high']] 


    def SetPromptBinNumber(self, nbins_p):
        self.nbins['prompt'] = nbins_p + 1

    def SetDelayedBinNumber(self, nbins_d):
        self.nbins['delayed'] = nbins_d + 1

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
        for key in self.cdict:
            for entry in self.cdict_pair: #apply to both prompt and delayed
                if self.cdict_pair[entry] is not None and entry == "interevent_dist":
                    self.precuts_data[key].append("interevent_dist<%s"%(str(self.cdict_pair[entry])))
                    self.precuts_mc[key].append("interevent_dist<%s"%(str(self.cdict_pair[entry])))
                if self.cdict_pair[entry] is not None and entry == "interevent_time_high":
                    self.precuts_data[key].append("interevent_time<%s"%(str(self.cdict_pair[entry])))
                    self.precuts_mc[key].append("interevent_time<%s"%(str(self.cdict_pair[entry])))
                if self.cdict_pair[entry] is not None and entry == "interevent_time_low":
                    self.precuts_data[key].append("interevent_time>%s"%(str(self.cdict_pair[entry])))
                    self.precuts_mc[key].append("interevent_time>%s"%(str(self.cdict_pair[entry])))
            if key == "prompt":
                suf = "_p"
            elif key == "delayed":
                suf = "_d"
            if self.cdict[key]['b14_high'] is not None:
                self.precuts_data[key].append("beta14%s<%s"%(suf,str(self.cdict[key]['b14_high'])))
                self.precuts_mc[key].append("beta14%s<%s"%(suf,str(self.cdict[key]['b14_high'])))
            if self.cdict[key]['b14_low'] is not None:
                self.precuts_data[key].append("beta14%s>%s"%(suf,str(self.cdict[key]['b14_low'])))
                self.precuts_mc[key].append("beta14%s>%s"%(suf,str(self.cdict[key]['b14_low'])))
            if self.cdict[key]['itr_high'] is not None:
                self.precuts_data[key].append("itr%s<%s"%(suf,str(self.cdict[key]['itr_high'])))
                self.precuts_mc[key].append("itr%s<%s"%(suf,str(self.cdict[key]['itr_high'])))
            if self.cdict[key]['itr_low'] is not None:
                self.precuts_data[key].append("itr%s>%s"%(suf,str(self.cdict[key]['itr_low'])))
                self.precuts_mc[key].append("itr%s>%s"%(suf,str(self.cdict[key]['itr_low'])))
            if self.cdict[key]['r_high'] is not None:
                self.precuts_data[key].append("posr%s<%s"%(suf,str(self.cdict[key]['r_high'])))
                self.precuts_mc[key].append("posr%s<%s"%(suf,str(self.cdict[key]['r_high'])))
            if self.cdict[key]['r_low'] is not None:
                self.precuts_data[key].append("posr%s>%s"%(suf,str(self.cdict[key]['r_low'])))
                self.precuts_mc[key].append("posr%s>%s"%(suf,str(self.cdict[key]['r_low'])))
            if self.cdict[key]['nhits_high'] is not None:
                self.precuts_data[key].append('nhits%s<%s'%(suf,str(self.cdict[key]['nhits_high'])))
                self.precuts_mc[key].append('nhits%s<%s'%(suf,str(self.cdict[key]['nhits_high'])))
            if self.cdict[key]['nhits_low'] is not None:
                self.precuts_data[key].append('nhits%s>%s'%(suf,str(self.cdict[key]['nhits_low'])))
                self.precuts_mc[key].append('nhits%s>%s'%(suf,str(self.cdict[key]['nhits_low'])))
            if self.cdict[key]['nhitsCleaned_high'] is not None:
                self.precuts_data[key].append('nhitsCleaned%s<%s'%(suf,str(self.cdict[key]['nhitsCleaned_high'])))
                self.precuts_mc[key].append('nhitsCleaned%s<%s'%(suf,str(self.cdict[key]['nhitsCleaned_high'])))
            if self.cdict[key]['nhitsCleaned_low'] is not None:
                self.precuts_data[key].append('nhitsCleaned%s>%s'%(suf,str(self.cdict[key]['nhitsCleaned_low'])))
                self.precuts_mc[key].append('nhitsCleaned%s>%s'%(suf,str(self.cdict[key]['nhitsCleaned_low'])))
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
            if self.cdict[key]["path_trigmask"] is not None:
                self.precuts_data[key].append("((triggerWord%s&%s)==0)" % (suf,self.cdict[key]["path_trigmask"]))
            if self.cdict[key]["path_trigmask_not_exact"] is not None:
                self.precuts_data[key].append("(triggerWord%s!=%s)" % (suf,self.cdict[key]["path_trigmask_not_exact"]))
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

    def _GetTopSacs(self,topnumber=7):
        for key in self.cdict:
            dcmask = self.cdict[key]['sacDCmask']
            for dattype in ['data','MC']:
                if len(self.sac_percut[key][dattype]) <= topnumber:
                    print("You already only have five cuts in %s.  Not combining any"%(dattype))
                    return
                names = []
                sacrifices = []
                for cut in self.sac_percut[key][dattype]:
                    sacrifices.append(np.average(self.sac_percut[key][dattype][cut]["fractional_sacrifice"]))
                    names.append(cut)
                sortednames = [x for _,x in sorted(zip(sacrifices,names))]
                topnames = sortednames[(len(sortednames)-(topnumber+1)):len(sortednames)]
                self.top_cuts[key][dattype] = topnames


    def AnalyzeData(self, var="nhits",onlycuts=[]):
        self.sac_percut_metadata = {}
        '''Analyzes the DC sacrifice for the data and Monte Carlo data fed to
        the class.  Returns the sacrifice as a function of the variable given'''
        havebounds = {'prompt':False,'delayed':False}
        xminimum = {'prompt':None, 'delayed':None}
        xmaximum = {'prompt':None, 'delayed':None}
        for key in self.analyze_range:
            for v in self.analyze_range[key]:
                xminimum[key],xmaximum[key] = self.analyze_range[key][var][0],self.analyze_range[key][var][1]
                havebounds[key] = True
            if not havebounds[key]:
                if xminimum[key] is None or xmaximum[key] is None:
                    print("You do not have bounds defined for analyzing this variable"+\
                            "defined.  Check your config file (%s type not defined)."%(key))
                    sys.exit(0)
            varname = var
            if var=='posr3':
                xminimum[key] = (float(xminimum[key])/6000.0)**3
                xmaximum[key] = (float(xmaximum[key])/6000.0)**3
        
        #Load the data trees from data and MC rootfiles
        data = ROOT.TChain("CombinedOutput")
        MC = ROOT.TChain("CombinedOutput") 
        for rf in self.rootfiles_data:
            data.Add(rf)
        for mf in self.rootfiles_mc:
            MC.Add(mf)

        for key in self.cdict: 
            xmin = xminimum[key]
            xmax = xmaximum[key]
            nbins = self.nbins[key]
            if key=="prompt":
                suf = "_p"
            elif key=="delayed":
                suf = "_d"
            dcmask = self.cdict[key]['sacDCmask']
            if dcmask == None:
                print("No DC mask in cuts config.  Please give a DC mask of cuts to plot")
                sys.exit(0)
            fullmask = mb.get_dcwords()
            print(fullmask)
            plotmask = {}
            for cut in fullmask:
                if cut not in onlycuts and len(onlycuts)>0:
                    continue
                if 'prescale' in fullmask[cut]:
                    continue
                if 'waterblind' in fullmask[cut]:
                    continue
                if mb.binary_bit_to_int(cut)&dcmask==mb.binary_bit_to_int(cut):
                    print("WE ARE ADDING IN CUT: " + str(cut))
                    print("THAT IS... CUT: " + str(fullmask[cut]))
                    plotmask[cut] = fullmask[cut]
            if "total" in onlycuts and len(onlycuts) > 0:
                plotmask[dcmask] = "total"
            elif len(onlycuts)==0:
                plotmask[dcmask] = "total"
            data.Draw("%s%s>>h_dallevents(%i,%f,%f)"% (var,suf,nbins,xmin,xmax),
                    "%s" % (self.precuts_data[key]),"goff")
            print("DATA CUTSTRING: " + "%s" % (self.precuts_data[key]))
            numdata_precuts = data.GetEntries("%s" % (self.precuts_data[key]))
            self.sac_percut_metadata["numdataevts_precuts_%s"%(key)] = numdata_precuts
            MC.Draw("%s%s>>h_mallevents(%i,%f,%f)"% (var,suf,nbins,xmin,xmax),
                    "%s" % (self.precuts_mc[key]),"goff")
            print("DATA CUTSTRING: " + "%s" % (self.precuts_mc[key]))
            nummc_precuts = data.GetEntries("%s" % (self.precuts_mc[key]))
            self.sac_percut_metadata["nummcevts_precuts_%s"%(key)] = nummc_precuts
            print("TEST: ENTRIES IN DATA AFTER PRECUTS")
            print(data.GetEntries("%s"%(self.precuts_data[key]))) 
            h_dallevents = gDirectory.Get("h_dallevents")
            h_dallevents.Sumw2()
            h_mallevents = gDirectory.Get("h_mallevents")
            h_mallevents.Sumw2()
                
            cutnames = ()
            self.sac_percut[key] = {'data':{}, 'MC':{}}
            for cut in plotmask:
                if cut!= dcmask:
                    cutint = mb.binary_bit_to_int(cut)
                else:
                    cutint = cut
                h_cut_dFracFlagged = ROOT.TH1D("h_cut_dFracFlagged", "h_cut_dFracFlagged", nbins,xmin,xmax)
                h_cut_dFracFlagged.Sumw2()
                h_cut_mFracFlagged = ROOT.TH1D("h_cut_mFracFlagged", "h_cut_mFracFlagged", nbins,xmin,xmax)
                h_cut_mFracFlagged.Sumw2()
                if self.precuts_data[key] is not None:
                    data.Draw("%s%s>>h_dflagged(%i,%f,%f)" % (var,suf,nbins,xmin,xmax),
                            "%s"%(self.precuts_data[key])+\
                                 "&&((dcFlagged%s&%i)!=%i)" % (suf,cutint,cutint),"goff")
                    print("TEST: NUM EVENTS AFTER DC APPLIED")
                    print("THECUTS: %s"%("%s"%(self.precuts_data[key])+\
                                 "&&((dcFlagged%s&%i)!=%i)" % (suf,cutint,cutint)))
                    print(data.GetEntries("%s"%(self.precuts_data[key])+\
                                 "&&((dcFlagged%s&%i)!=%i)" % (suf,cutint,cutint)))
                else:
                    data.Draw("%s%s>>h_dflagged(%i,%f,%f)" % (var,suf,nbins,xmin,xmax),
                            "((dcFlagged%s&%i)!=%i)" % (suf,cutint,cutint),"goff")
                if self.precuts_mc[key] is not None:
                    MC.Draw("%s%s>>h_mflagged(%i,%f,%f)" % (var,suf,nbins,xmin,xmax),
                            "%s"%(self.precuts_mc[key])+\
                                 "&&((dcFlagged%s&%i)!=%i)" % (suf,cutint,cutint),"goff")
                else:
                    MC.Draw("%s%s>>h_mflagged(%i,%f,%f)" % (var,suf,nbins,xmin,xmax),
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
                self.sac_percut_metadata["binwidth_%s"%(key)] = (xmaximum[key]-xminimum[key])/float(self.nbins[key])
                del h_cut_dFracFlagged
                del h_dflagged
                del h_cut_mFracFlagged
                del h_mflagged
        #graphdict has the sacrifice information for each cut. Now, let's plot it.
        self.sac_percut_metadata["variable"] = varname 
        self._GetTopSacs()
        #self._DeleteEmptyBins()
        print self.sac_percut
        return self.sac_percut, self.sac_percut_metadata

    def DrawCleanHist(self, evtype="pair",dattype="data",var="interevent_time",xmin=0., xmax=5.E5,nbins=100,addlROOTcuts=""):
        histMeta = {"var":var, "evtype":evtype, "datatype":dattype, "xrange":[xmin,xmax],"nbins":nbins}
        
        #Load the data trees from data and MC rootfiles
        data = ROOT.TChain("CombinedOutput")
        MC = ROOT.TChain("CombinedOutput") 
        for rf in self.rootfiles_data:
            data.Add(rf)
        for mf in self.rootfiles_mc:
            MC.Add(mf)
        if evtype=="prompt":
            suf = "_p"
        elif evtype=="delayed":
            suf = "_d"
        else:
            suf = ""
        dcmask_prompt = self.cdict["prompt"]['sacDCmask']
        dcmask_delayed = self.cdict["delayed"]['sacDCmask']
        fullmask = mb.get_dcwords()
        print(fullmask)
        plotmask = {}
        plotmask["_p"] = dcmask_prompt
        plotmask["_d"] = dcmask_delayed 
        print("THEPLOTMASK: " + str(plotmask))
        precuts_data = self.precuts_data["prompt"] + "&&" + self.precuts_data["delayed"]
        precuts_mc = self.precuts_mc["prompt"] + "&&" + self.precuts_mc["delayed"]
        dcstring = ""
        for cut in plotmask:
            dcstring+="&&((dcFlagged%s&%i)==%i)" % (cut, plotmask[cut],plotmask[cut])

        #h_CleanDistHist = ROOT.TH1D("h_CleanDistHist", "h_CleanDistHist", nbins,xmin,xmax)
        #h_CleanDistHist.Sumw2()
        if dattype=="data":
            if addlROOTcuts:
                data.Draw("%s%s>>h_CleanDistHist(%i,%f,%f)" % (var,suf,nbins,xmin,xmax),
                          "%s%s&&%s"%(precuts_data,dcstring,addlROOTcuts),"goff")
            else:
                data.Draw("%s%s>>h_CleanDistHist(%i,%f,%f)" % (var,suf,nbins,xmin,xmax),
                          "%s%s"%(precuts_data,dcstring),"goff")
            print("CUTS USED: %s%s&&%s"%(precuts_data,dcstring,addlROOTcuts))
        elif dattype=="MC":
            if addlROOTcuts:
                mc.Draw("%s%s>>h_CleanDistHist(%i,%f,%f)" % (var,suf,nbins,xmin,xmax),
                        "%s%s&&%s"%(precuts_mc,dcstring,addlROOTcuts),"goff")
            else:
                mc.Draw("%s%s>>h_CleanDistHist(%i,%f,%f)" % (var,suf,nbins,xmin,xmax),
                        "%s%s"%(precuts_mc,dcstring),"goff")
            print("CUTS USED: %s%s&&%s"%(precuts_mc,dcstring,addlROOTcuts))

        h_CleanDistHist = gDirectory.Get("h_CleanDistHist")
        h_CleanDistHist.Sumw2()
        histdict = self._histToDictXY(h_CleanDistHist)
        histPD= pandas.DataFrame(data=histdict)
        del h_CleanDistHist
        #graphdict has the sacrifice information for each cut. Now, let's plot it.
        return histPD, histMeta

    def _histToDictXY(self,roothist):
        graphdict={}
        x, y,y_unc =(), (), () #pandas wants ntuples
        for i in xrange(int(roothist.GetNbinsX()+1)):
            if i==0:
                continue
            x =  x + ((float(roothist.GetBinWidth(i))/2.0) + float(roothist.GetBinLowEdge(i)),)
            y = y + (roothist.GetBinContent(i),)
            y_unc = y_unc + (roothist.GetBinError(i),)
        graphdict["x"] = x
        graphdict["y"] = y
        graphdict["y_unc"] = y_unc
        return graphdict

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

    #TODO: We need to bring this into the "multiple datatypes" and "prompt or delayed"
    #Future
    def ShowSacrificePlot(self,fittotal=True,title=None,evtype="prompt",dattype='data',savedir="."):
        '''Plots the data cleaning sacrifice as a function of the analyzed variable'''
        sns.set_style("whitegrid")
        xkcd_colors = ['black','slate blue', 'fluro green', 'red', 'blue',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'leaf',
                'aqua blue','vomit', 'brown','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        fig.set_size_inches(18.5,11.5)
        ax = fig.add_subplot(1,1,1)
        if self.top_cuts[evtype][dattype] is None:
            for cut in self.sac_percut[evtype][dattype]:
                thiscut_dict = self.sac_percut[evtype][dattype][cut]
                plt.errorbar(x=thiscut_dict.vardat, 
                        y=thiscut_dict.fractional_sacrifice,
                        yerr=thiscut_dict.fs_uncertainty,
                        linestyle='none', marker='o', label=cut, markersize=9,
                        elinewidth=2, capsize=0)
        else:
            fracsum,fracuncsum=[],[]
            for cut in self.sac_percut[evtype][dattype]:
                if cut in self.top_cuts[evtype][dattype]:
                    thiscut_dict = self.sac_percut[evtype][dattype][cut]
                    plt.errorbar(x=thiscut_dict.vardat, 
                            y=thiscut_dict.fractional_sacrifice,
                            yerr=thiscut_dict.fs_uncertainty,
                            linestyle='none', marker='o', label=cut, markersize=9,
                            elinewidth=2, capsize=0)
                else:
                    fracsum.append(self.sac_percut[evtype][dattype][cut].fractional_sacrifice)
                    fracuncsum.append(self.sac_percut[evtype][dattype][cut].fs_uncertainty)
            plt.errorbar(x=self.sac_percut[evtype][dattype][cut].vardat,y=sum(fracsum),yerr=sum(fracuncsum),
                    linestyle='none',marker='o', capsize=0, elinewidth=2, label='All other cuts', markersize=9)
    
        if fittotal is True:
            popt, pcov = spc.curve_fit(self._flatline, 
                    self.sac_percut[evtype][dattype]["total"].vardat, 
                    self.sac_percut[evtype][dattype]["total"].fractional_sacrifice,
                    p0=[0.02], 
                    sigma=self.sac_percut[evtype][dattype]["total"].fs_uncertainty)
                    #one standard deviation
            try:
                stdev = np.sqrt(np.diag(pcov))
            except ValueError:
                print("Fit to total sacrifice likely failed.  You may not have enough data"+\
                        "in this region to fit a straight line.")
                raise
            wstdev = self._weighted_stdev(self.sac_percut[evtype][dattype]["total"].fractional_sacrifice,
                    float(popt[0]),self.sac_percut[evtype][dattype]["total"].fs_uncertainty)
            plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f, fit unc = %f, dev. from fit = %f $' % (float(popt[0]), float(stdev[0]),float(wstdev)))
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
            plt.title("Fractional sacrifice due to data cleaning cuts \n"+ \
                    "datatype %s for event type %s"%(dattype,evtype),fontsize=36)
        else: 
            plt.title(title,fontsize=36)
        variable = self.sac_percut_metadata["variable"]
        plt.savefig(savedir+"/DCSac_%s_%s_%s.pdf"%(dattype,evtype,variable))
        plt.show()
        plt.close()

    def DataMCSacCompare(self,title=None,evtype="prompt",savedir="."):
        '''Plots the total sacrifice evaluated for the Data and MC'''
        sns.set_style("whitegrid")
        xkcd_colors = ['black','slate blue', 'fluro green', 'red', 'blue',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'leaf',
                'aqua blue','vomit', 'brown','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        fig.set_size_inches(18.5,11.5)
        ax = fig.add_subplot(1,1,1)
        for dattype in self.sac_percut[evtype]: 
            thisdat_sactotal = self.sac_percut[evtype][dattype]["total"]
            plt.errorbar(x=thisdat_sactotal.vardat, 
                    y=thisdat_sactotal.fractional_sacrifice,
                    yerr=thisdat_sactotal.fs_uncertainty,
                    linestyle='none', marker='o', label=dattype, markersize=9,
                    elinewidth=2, capsize=0)
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
            plt.title("Sacrifice due to data cleaning cuts \n" + \
                    "in Data and MC (event type %s)"%(evtype),fontsize=36)
        else: 
            plt.title(title,fontsize=36)
        variable = self.sac_percut_metadata["variable"]
        plt.savefig(savedir+"/DataMCSacrifices_%s_%s.pdf"%(evtype,variable))
        plt.show()
        plt.close()

    def PlotDataMCSacRatio(self,fittotal=True,title=None,evtype='prompt',savedir="."):
        '''Plots the ratio of the data to MC total sacrifice bin-by-bin, as well as the best fit
        if title is None, has a default title to go with'''
        if self.sac_percut is None:
            print("You must first analyze the input root data!")
            return
        sns.set_style("whitegrid")
        xkcd_colors = ['blue','black',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
                'aqua blue','vomit', 'black','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        fig.set_size_inches(18.5,11.5)
        ax = fig.add_subplot(1,1,1)
        #self.sac_percut = self._DeleteEmptyBins_all(self.sac_percut) 
        ratio = pandas.Series(self.sac_percut[evtype]['data']["total"].fractional_sacrifice /
                self.sac_percut[evtype]['MC']["total"].fractional_sacrifice) 
        ratio_unc = np.sqrt(self.sac_percut[evtype]["MC"]["total"].fs_uncertainty**2 + \
                self.sac_percut[evtype]["data"]["total"].fs_uncertainty**2)
        ratio_unc = pandas.Series(ratio_unc)
        plt.errorbar(x=self.sac_percut[evtype]['data']["total"].vardat, 
                y=ratio,
                yerr=ratio_unc,
                linestyle='none', marker='o', label="Data/MC Ratio", markersize=12,
                elinewidth=4, capsize=0)
        if fittotal is True:
            popt, pcov = spc.curve_fit(self._flatline, self.sac_percut[evtype]['data']['total'].vardat, ratio,
                    p0=[0.98], sigma=ratio_unc)
            #one standard deviation
            print("BEST FIT: " + str(popt))
            print("PCOVARIANCE: " + str(pcov))
            stdev = np.sqrt(np.diag(pcov))
            plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
            print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(self._weighted_stdev(ratio,
                float(popt[0]),ratio_unc)))
        legend = plt.legend(loc=3,frameon=1)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        #plt.yscale("log")
        plt.ylabel("Ratio",fontsize=34)
        varindict = False
        for var in self.xlabel_dict:
            if self.sac_percut_metadata["variable"] == var:
                plt.xlabel(self.xlabel_dict[var],fontsize=32)
                varindict = True
            else:
                if varindict is False: 
                    plt.xlabel(self.sac_percut_metadata["variable"],fontsize=32)
        if title is None:
            plt.title("Data/MC data cleaning sacrifice ratio \n"+ \
                    "event type %s"%(evtype),fontsize=36)
        else: 
            plt.title(title,fontsize=36)
        plt.tick_params(labelsize=32)
        variable = self.sac_percut_metadata["variable"]
        plt.savefig(savedir+"/DataMCRatio_%s_%s.pdf"%(evtype,variable))
        plt.show()
        plt.close()

    def PlotDataMCAccRatio(self,fittotal=True,title=None,evtype='prompt',savedir="."):
        '''Plots the ratio of the data to MC total sacrifice bin-by-bin, as well as the best fit
        if title is None, has a default title to go with'''
        if self.sac_percut is None:
            print("You must first analyze the input root data!")
            return
        sns.set_style("whitegrid")
        xkcd_colors = ['blue','black',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
                'aqua blue','vomit', 'black','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        fig.set_size_inches(18.5,11.5)
        ax = fig.add_subplot(1,1,1)
        #self.sac_percut = self._DeleteEmptyBins_all(self.sac_percut)
        data_sac = self.sac_percut[evtype]['data']["total"].fractional_sacrifice 
        mc_sac = self.sac_percut[evtype]['MC']["total"].fractional_sacrifice
        data_acc = np.ones(len(data_sac)) - data_sac
        mc_acc = np.ones(len(mc_sac)) - mc_sac
        ratio =  data_acc / mc_acc
        ratio_unc = np.sqrt(self.sac_percut[evtype]["MC"]["total"].fs_uncertainty**2 + \
                self.sac_percut[evtype]["data"]["total"].fs_uncertainty**2)
        ratio_unc = pandas.Series(ratio_unc)
        plt.errorbar(x=self.sac_percut[evtype]['data']["total"].vardat, 
                y=ratio,
                yerr=ratio_unc,
                linestyle='none', marker='o', label="Data/MC Acc.Ratio", markersize=12,
                elinewidth=5, capsize=0)
        if fittotal is True:
            popt, pcov = spc.curve_fit(self._flatline, self.sac_percut[evtype]['data']['total'].vardat, ratio,
                    p0=[0.98], sigma=ratio_unc)
            #one standard deviation
            print("BEST FIT: " + str(popt))
            print("PCOVARIANCE: " + str(pcov))
            stdev = np.sqrt(np.diag(pcov))
            plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
            print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(self._weighted_stdev(ratio,
                float(popt[0]),ratio_unc)))
        legend = plt.legend(loc=3,frameon=1)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        #plt.yscale("log")
        plt.ylabel("Ratio",fontsize=34)
        varindict = False
        for var in self.xlabel_dict:
            if self.sac_percut_metadata["variable"] == var:
                plt.xlabel(self.xlabel_dict[var],fontsize=32)
                varindict = True
            else:
                if varindict is False: 
                    plt.xlabel(self.sac_percut_metadata["variable"],fontsize=32)
        if title is None:
            plt.title("Data/MC data cleaning acceptance ratio \n"+ \
                    "event type %s"%(evtype),fontsize=36)
        else: 
            plt.title(title,fontsize=36)
        plt.tick_params(labelsize=32)
        variable = self.sac_percut_metadata["variable"]
        plt.savefig(savedir+"/DataMCRatio_%s_%s.pdf"%(evtype,variable))
        plt.show()
        plt.close()

