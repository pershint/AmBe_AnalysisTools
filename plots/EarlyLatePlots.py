import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as spc
import pandas
import time


class EarlyLateSacComparer(object):
    '''This class takes in two sac_percut pandas dataframes and two metadata
    dictionaries that are returned after running the AmBeDCSacrifice method
    AnalyzeData on SNO+ ntules.'''

    def __init__(self,dataframe1, dataframe2, data1meta, data2meta):
        self.d1 = dataframe1
        self.d2 = dataframe2
        self.m1 = data1meta
        self.m2 = data2meta
        self.xlabel_dict = {'energy': 'Energy (MeV)', 'udotr': r'$U . R$',
                'posr3': r'$(R/R_{AV})^{3}$','posr':'Radius (mm)'}

    def _weighted_stdev(self,vals,valavg, uncs):
        weights = 1/(uncs**2)
        return np.sqrt((len(weights)*np.sum(weights*(vals-valavg)**2))/((len(weights)-1)*np.sum(weights)))
    
    def _flatline(self,x,b):
        return b
    
    def PlotDataSacTotalDifference(self,fitdiff=True,title=None,evtype='delayed',savedir="."):
        '''Plots the weighted difference of the sacrifices of ofdataframe1 and dataframe2.
        the error is taken as the two errors in quadrature'''
        if self.d1 is None or self.d2 is None:
            print("You must first give the class dataframes to parse!")
            return
        sns.set_style("whitegrid")
        xkcd_colors = ['blue','black',
                'yellowish orange', 'warm pink', 'light eggplant', 'clay', 'red', 'leaf',
                'aqua blue','vomit', 'black','twilight']
        sns.set_palette(sns.xkcd_palette(xkcd_colors))#,len(self.sac_percut)))
        fig = plt.figure()
        fig.set_size_inches(18.5,11.5)
        ax = fig.add_subplot(1,1,1)

        #Load the sacrifices, and the total number of events of each type
        d1_sac = self.d1[evtype]['data']["total"].fractional_sacrifice
        d1_sac_unc = self.d1[evtype]['data']['total'].fs_uncertainty
        d2_sac = self.d2[evtype]['data']["total"].fractional_sacrifice 
        d2_sac_unc = self.d2[evtype]['data']['total'].fs_uncertainty
        d1_acc = np.ones(len(d1_sac)) - d1_sac
        d2_acc = np.ones(len(d2_sac)) - d2_sac
        num_d1 = self.m1["numdataevts_precuts_%s"%(evtype)]
        num_d2 = self.m2["numdataevts_precuts_%s"%(evtype)]
        num_tot = num_d1 + num_d2
        #Do a weighted subtraction of the two
        subterm = (float(num_d2)/float(num_tot))*d2_sac
        subterm_unc = subterm*np.sqrt((1/float(num_d2)) + (1/float(num_tot)) + (d2_sac_unc/d2_sac)**2)
        d1corr_sac = d1_sac - (float(num_d2)/float(num_tot))*d2_sac
        d1corr_sac_unc = np.sqrt(d1_sac_unc**2 + subterm_unc**2)
        d1corr_sac_unc = pandas.Series(d1corr_sac_unc)
        plt.errorbar(x=self.d1[evtype]['data']["total"].vardat, 
                y=d1corr_sac,
                yerr=d1corr_sac_unc,
                linestyle='none', marker='o', label="Delayed sacrifice", markersize=12,
                elinewidth=5, capsize=0)
        if fitdiff is True:
            popt, pcov = spc.curve_fit(self._flatline, self.d1[evtype]['data']['total'].vardat, d1corr_sac,
                    p0=[0.98], sigma=d1corr_sac_unc)
            #one standard deviation
            print("BEST FIT: " + str(popt))
            print("PCOVARIANCE: " + str(pcov))
            stdev = np.sqrt(np.diag(pcov))
            plt.axhline(popt, linewidth=3, alpha=0.8, color='k',label= r'fit: $\mu = %f,unc = %f$' % (float(popt[0]), float(stdev[0])))
            print("WEIGHTED STANDARD DEVIATION FROM FLAT: " + str(self._weighted_stdev(d1corr_sac,
                float(popt[0]),d1corr_sac_unc)))
        legend = plt.legend(loc=3,frameon=1)
        frame = legend.get_frame()
        frame.set_facecolor("white")
        y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
        ax.yaxis.set_major_formatter(y_formatter)
        #plt.yscale("log")
        plt.ylabel("Fractional sacrifice",fontsize=34)
        varindict = False
        for var in self.xlabel_dict:
            if self.m1["variable"] == var:
                plt.xlabel(self.xlabel_dict[var],fontsize=32)
                varindict = True
            else:
                if varindict is False: 
                    plt.xlabel(self.m1["variable"],fontsize=32)
        if title is None:
            plt.title("Delayed data cleaning sacrifice corrected for background estimate \n"+ \
                    "event type %s"%(evtype),fontsize=36)
        else: 
            plt.title(title,fontsize=36)
        plt.tick_params(labelsize=32)
        variable = self.m1["variable"]
        plt.savefig(savedir+"/DelayedCorrSacrifice_%s_%s.pdf"%(evtype,variable))
        plt.show()
        plt.close()
