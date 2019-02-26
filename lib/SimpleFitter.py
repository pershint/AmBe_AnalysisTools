import numpy as np
import scipy.optimize as scp
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")
sns.set_context("poster")

class Fitter(object):
    '''
    General class for using SciPy's curvefit class to
    fit a given function to data.
    Parameters
    ----------------
    data : numpy array
    numpy array of data that will be fit to using the given function
    fitfunc : function that will be fit to the data.
    '''
    def __init__(self,datax=None,datay=None,datasigma=None,fitfunc=None,numvars=None):
        self.datax = datax
        self.datay = datay
        self.datasigma=datasigma
        self.f = fitfunc
        self.numfitvar = numvars

    def SetFitFunction(self,fitfunc,numvars):
        self.f = fitfunc
        self.numfitvar = numvars

    def SetDataX(self,dat):
        self.datax = dat

    def SetDataY(self,dat):
        self.datay = dat

    def _CheckInitialVariableLength(self,p0):
        if len(p0) != self.numfitvar:
            print("OH NO, YOUR INITIAL VARIABLE ARRAY LENGTH DOESN'T " +\
                  "MATCH YOUR FUNCTION'S NUMBER OF VARIABLES.")

    def SetDataUncertainties(self,datunc):
        self.datasigma = datunc
    
    def RunFit(self,initialvars):
        print("WERUN")
        popt, pcov = None, None
        self._CheckInitialVariableLength(initialvars)
        if self.datasigma is None:
            popt, pcov = scp.curve_fit(self.f,self.datax, self.datay, p0=initialvars)
        else:
            popt, pcov = scp.curve_fit(self.f,self.datax, self.datay, p0=initialvars,
                    sigma=self.datasigma)
        print("POPT: " + str(popt))
        print("PCOV: " + str(pcov))
        return popt, pcov

if __name__ == '__main__':
    doubleexp = lambda x,A1,l1,A2,l2: A1*np.exp(-l1*x) + A2*np.exp(-l2*x)
    testy = np.array([20., 12., 7., 3., 2., 1.4, 1.2, 1.1, 1.05, 1.003])
    testx = np.arange(0,len(testy))
    myfitter = Fitter()
    myfitter.SetDataX(testx)
    myfitter.SetDataY(testy)
    myfitter.SetFitFunction(doubleexp, 4)
    initvars = [15., 0.1, 3.,1.]
    popt, pcov = myfitter.RunFit(initvars)
    plt.plot(testx,testy,linestyle='none',marker='o',markersize=5)
    bestfitline = doubleexp(testx,popt[0],popt[1],popt[2],popt[3])
    plt.plot(testx,bestfitline)
    plt.xlabel("testx")
    plt.ylabel("testy")
    plt.title("Test title")
    plt.show()
