#Here are our fitter tools
import ROOT
import numpy as np

class FitEngine(object):
    def __init__(self, TH1=None, distribution=None,fit_range=[0.0,100.0]):
        self.histogram = TH1
        self.distribution = distribution
        if distribution is not None:
            self.initial_parameters = self.__dlib[distribution]
        else:
            self.initial_parameters = None
        self.current_fit = None
        self.fit_range = fit_range

    def DrawHistogram(self):
        if self.histogram is None:
            print("No histogram available to draw right now.")
        else:
            self.histogram.Draw()

    def DrawFit(self):
        if self.current_fit is None:
            print("No fit available to draw.  Run MakeFit() first.")
        else:
            self.histogram.Draw()
            self.current_fit.Draw("same")

    def setHistogramToFit(self,TH1):
        self.histogram = TH1

    def setFitRange(self,fitrange):
        if len(fitrange) is not 2:
            print("Please give an array with two elements [min,max]")
        if fitrange[0] >= fitrange[1]:
            print("Please choose a min less than the max value")
        self.fit_range = fitrange
    
    def TimeDiffFit(self, x,p):
        a = p[0]*p[1]*p[2]*p[3]*np.exp(-p[3]*x[0])
        b = (1/p[4])*p[5]*p[2] + p[6] + p[0]*p[2]
        alpha = ((1-p[4])/p[4])*p[2]*p[1] + p[7] + p[2]*p[8]*(1.0-p[1])
        return p[9]*(a + b*alpha*np.exp(-alpha*x[0]))

    def MakeFit(self):
        self.current_fit = None
        #Function takes in a time histogram and fits a polynomial to it. Returns
        #The polynomial function for random firing of values from
        TimeFit = ROOT.TF1('TimeFit', self.TimeDiffFit, self.fit_range[0],
                    self.fit_range[1], 10)
        efficiencyParamters = [0,1,5,8]
        for e in efficiencyParamters:
            TimeFit.SetParLimits(e,0.0,1.0)
            TimeFit.SetParameter(e,0.5)
        TimeFit.FixParameter(2,0.58*66.0E-09)
        #TimeFit.SetParameter(2,0.58*66.0E-9)
        #TimeFit.SetParLimits(2,0.58*52.80E-9,0.58*79.2E-09) #Conservative 20% unc.
        TimeFit.FixParameter(4,0.58)
        TimeFit.SetParLimits(3,0,1E-4)
        TimeFit.SetParameter(3,1E-5)
        TimeFit.SetParLimits(6,2.95E-9,2.97E-9)
        TimeFit.SetParameter(6,2.96E-9)
        TimeFit.SetParLimits(7,489.29E-9,489.53E-9)
        TimeFit.SetParameter(7,489.40E-9)
        TimeFit.SetParameter(9,10000000000000000.0)
        TimeFit.SetLineColor(2)
        TimeFit.SetLineWidth(4)
        TimeFit.SetLineStyle(2)
        print(TimeFit)
        ROOT.gStyle.SetOptFit(0157)
        self.histogram.Fit('TimeFit','Lq')
        self.current_fit = TimeFit
        return TimeFit
   
class ChainShooter(object):
    def __init__(self,datafiles=None,variable_list=None):
        print("NEED TO MAKE")

    def MakeTH1(tree,branch):
        #Return a TH1 histogram made of all data for a branch in the given
        #data files
        for f in datafiles:
            thefile = ROOT.TFile(f,"READ")
            t=thefile.Get(tree)
            for entry in xrange(t.GetEntries()):
                t.GetEntry(entry)
                #Fill hisgotram with GetProperty(t,branch)
