import numpy as np
import scipy.optimize as scp
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")
sns.set_context("poster")

class NSacCalculator(object):
    '''
    General class for using SciPy's curvefit class to
    fit a given function to data.
    Parameters
    ----------------
    data : numpy array
    numpy array of data that will be fit to using the given function
    fitfunc : function that will be fit to the data.
    '''
    def __init__(self,ithist_x,itcorr_y,ituncorr_y,itearlywindow,itlatewindow,
                 earlysac, earlysac_unc, latesac,latesac_unc):
        self.itx = ithist_x
        self.corry = itcorr_y
        self.uncorry = ituncorr_y
        self.it_ewin = itearlywindow
        self.it_lwin = itlatewindow
        self.es = earlysac
        self.es_unc = earlysac_unc
        self.ls = latesac
        self.ls_unc = latesac_unc

    def _SumInWindow(self,x,y,xmin,xmax):
        sumind = np.where(x<xmax)[0]
        sumind2 = np.where(x>xmin)[0]
        sum_indices = np.intersect1d(sumind,sumind2)
        return np.sum(y[sum_indices])

    def CalculateNeutronSacrifice(self):
        N_ec = self._SumInWindow(self.itx,self.corry,self.it_ewin[0],
                                self.it_ewin[1])
        N_lc = self._SumInWindow(self.itx,self.corry,self.it_lwin[0],
                                self.it_lwin[1])
        N_eu = self._SumInWindow(self.itx,self.uncorry,self.it_ewin[0],
                                self.it_ewin[1])
        N_lu = self._SumInWindow(self.itx,self.uncorry,self.it_lwin[0],
                                self.it_lwin[1])
        N_l = N_lc + N_lu
        N_e = N_ec + N_eu
        numerator = N_l*N_eu*self.ls - N_e*N_lu*self.es
        denom = N_eu*N_lc - N_lu*N_ec
        corr_sacrifice = numerator/denom
        return corr_sacrifice

    #WHAT WE NEED:
    #  - EVENTUALLY, we need to figure out the uncertainties here.  Let's
    #    Have a function that re-shoots all values based on their 
    #    uncertainties.  We'll have to give the fit function and parameters 
    #    to the class to do this though.


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
