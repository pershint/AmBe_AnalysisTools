import numpy as np
import scipy.optimize as scp
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")
sns.set_context("poster")

class NeutronSacrificeCalculator(object):
    '''
    General class for using SciPy's curvefit class to
    fit a given function to data.
    Parameters
    ----------------
    data : numpy array
    numpy array of data that will be fit to using the given function
    fitfunc : function that will be fit to the data.
    '''
    def __init__(self,ithist_x,itcorr_y,ituncorr_y,itearlywindow,itlatewindow,earlysac, latesac):
        self.itx = ithist_x
        self.corry = itcorr_y
        self.uncorry = ituncorr_y
        self.it_ewin = itearlywindow
        self.it_lwin = itlatewindow
        self.es = earlysac
        self.ls = latesac

    #WHAT WE NEED:
    #  - A function that will do the integral in the early and late windows 
    #    for the uncorrelated and correlated terms
    #  - The actual calculation of the delayed sacrifice
    #  - EVENTUALLY, we need to figure out the uncertainties here...


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
