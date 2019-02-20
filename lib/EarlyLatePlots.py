import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as spc
import pandas
import time


class SacComparer(object):
    '''This class takes in two sac_percut pandas dataframes and two metadata
    dictionaries that are returned after running the AmBeDCSacrifice method
    AnalyzeData on SNO+ ntules. Makes several plots specific to the AmBe analysis.'''

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

    def PlotSacrifices(self,fitdiff=True,title=None,dattype="data",evtype='delayed',savedir=".",label1="DataFrame 1", label2="DataFrame 2"):
        '''Plots the two total sacrifices from each dataframe on the same plot..'''
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
        d1_sac = self.d1[evtype][dattype]["total"].fractional_sacrifice
        d1_sac_unc = self.d1[evtype][dattype]['total'].fs_uncertainty
        d2_sac = self.d2[evtype][dattype]["total"].fractional_sacrifice 
        d2_sac_unc = self.d2[evtype][dattype]['total'].fs_uncertainty
        #plt.errorbar(x=self.d1[evtype][dattype]["total"].vardat, 
        #        y=d1_sac,
        #        yerr=d1_sac_unc,
        #        linestyle='none', marker='o', label=label1, markersize=12,
        #        elinewidth=5, capsize=0)
        plt.errorbar(x=self.d2[evtype][dattype]["total"].vardat, 
                y=d2_sac,
                yerr=d2_sac_unc,
                linestyle='none', marker='o', label=label2, markersize=12,
                elinewidth=5, capsize=0)
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
            plt.title("Comparison of fractional data cleaning sacrifices in AmBe internal data \n"+ \
                    "event type %s"%(evtype),fontsize=36)
        else: 
            plt.title(title,fontsize=36)
        plt.tick_params(labelsize=32)
        variable = self.m1["variable"]
        plt.savefig(savedir+"/SacrificeComparison_%s_%s.pdf"%(evtype,variable))
        plt.show()
        plt.close()

