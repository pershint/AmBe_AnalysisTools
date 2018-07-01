import json
import numpy as np
import sys

def init_data(f_arg,*args):
    print("SCRIPT NAME: %s"%(f_arg))
    print("FIRST ARG: %s"%(str(args[0][1])))
    data_filename = args[0][1]
    with open(data_filename,"r") as f:
        data_file = json.load(f)
    return data_file

class sac_calculator(object):
    def __init__(self, datafile=None):
        self.d = datafile
        self.sacrifice = None
        self.sac_uertainty = None
        self.top = None
        self.bot = None
        self.top_u = None
        self.bot_u = None


    def calculate_sacrifice(self):
        '''Calculate the sacrifice of neutron captures in AmBe events.
        input is a JSON file with data cleaning sacrifices and rates of events
        in N16 central data, AmBe central data, and background runs'''
        d = self.d
        t1 = d['e_PT']['v']*d['e_DT']['v']*d['N_T']['v']
        t2 = d['e_PA']['v']*d['e_DB']['v']*d['T_c']*\
                (d['R_PT']['v']-d['R_PB']['v'])*d['R_DB']['v']*d['dt']
        t3 = d['e_PB']['v']*d['e_DB']['v']*d['N_PBDB']['v']
        b1 = d['e_PA']['v']*d['PE']['v']*d['N_T']['v']
        b2 = d['e_PB']['v']*((1.0-d['PE']['v'])*d['N_T']['v'])
        b3 = d['e_PB']['v']*d['T_c']*(d['R_PT']['v']-d['R_PB']['v'])* \
                d['R_DB']['v']*d['dt']
        b4 = d['e_PB']['v']*d['N_PBDB']['v']
        self.top = t1 - t2 - t3
        print("TOP: %s"%(str(self.top)))
        self.bot = b1 + b2 - b3 - b4
        print("BOT: %s"%(str(self.bot)))
        self.sacrifice = self.top/self.bot
        print("SACRIFICE: %s"%(str(self.top/self.bot)))
    
    def calculate_uncertainty(self):
        '''Calculates the uertainty of the numerator and denominator for
        the delayed event sacrifice, then the sacrifice uertainty itself'''
    
        d = self.d
        #First, re-calculate the vues for the top and bottom terms
        t1 = d['e_PT']['v']*d['e_DT']['v']*d['N_T']['v']
        t2 = d['e_PA']['v']*d['e_DB']['v']*d['T_c']*\
                d['R_DB']['v']*d['dt']*(d['R_PT']['v']-d['R_PB']['v'])
        t3 = d['e_PB']['v']*d['e_DB']['v']*d['N_PBDB']['v']
        b1 = d['e_PA']['v']*d['PE']['v']*d['N_T']['v']
        b2 = d['e_PB']['v']*((1.0-d['PE']['v'])*d['N_T']['v'])
        b3 = d['e_PB']['v']*d['T_c']*(d['R_PT']['v']-d['R_PB']['v'])* \
                d['R_DB']['v']*d['dt']
        b4 = d['e_PB']['v']*d['N_PBDB']['v']
        #Now, get the uertainty for each term
        t1_u = t1 * np.sqrt(
                (d['e_PT']['u']/d['e_PT']['v'])**2 + 
                (d['e_DT']['u']/d['e_DT']['v'])**2 + 
                (d['N_T']['u']/d['N_T']['v'])**2)
        t2_u = t2 * np.sqrt(
                (d['e_PA']['u']/d['e_PA']['v'])**2 +
                (d['e_DB']['u']/d['e_DB']['v'])**2 +
                (d['R_DB']['u']/d['R_DB']['v'])**2 +
                (np.sqrt(d['R_PT']['u']**2 + d['R_PB']['u']**2)/ 
                        (d['R_PT']['v']-d['R_PB']['v']))**2)
        t3_u = t3 * np.sqrt(
                (d['e_PB']['u']/d['e_PB']['v'])**2 +
                (d['e_DB']['u']/d['e_DB']['v'])**2 +
                (d['N_PBDB']['u']/d['N_PBDB']['v'])**2)

        b1_u = b1 * np.sqrt(
                (d['e_PA']['u']/d['e_PA']['v'])**2 +
                (d['PE']['u']/d['PE']['v'])**2 + 
                (d['N_T']['u']/d['N_T']['v'])**2)
        b2_u = b2 * np.sqrt(
                (d['e_PB']['u']/d['e_PB']['v'])**2 +
                (d['PE']['u']/(1.0-d['PE']['v']))**2 + 
                (d['N_T']['u']/d['N_T']['v'])**2)
        b3_u = b3 * np.sqrt(
                (d['e_PB']['u']/d['e_PB']['v'])**2 +
                (np.sqrt(d['R_PT']['u']**2 + d['R_PB']['u']**2)/ 
                        (d['R_PT']['v']-d['R_PB']['v']))**2 +
                (d['R_DB']['u']/d['R_DB']['v'])**2)
        b4_u = b4 * np.sqrt(
                (d['e_PB']['u']/d['e_PB']['v'])**2 +
                (d['N_PBDB']['u']/d['N_PBDB']['v'])**2)

        #Now, calculate the top and bottom uncertainties
        self.top_unc = np.sqrt(t1_u**2 + t2_u**2 + t3_u**2)
        self.bot_unc = np.sqrt(b1_u**2 + b2_u**2 + b3_u**2 + b4_u**2)
        print("TOP TERM UNCERTAINTY: %s"%(str(self.top_unc)))
        print("BOTTOM TERM UNCERTAINTY: %s"%(str(self.bot_unc)))

        #Finally, the sacrifice uncertainty
        self.sac_uncertainty = self.sacrifice * np.sqrt(
                (self.top_unc/self.top)**2 + (self.bot_unc/self.bot)**2)
        print("SAC UNCERTAINTY: %s"%(str(self.sac_uncertainty)))


if __name__=="__main__":
    DF = init_data(sys.argv[0],sys.argv)
    mycalculator = sac_calculator(datafile=DF)
    mycalculator.calculate_sacrifice()
    mycalculator.calculate_uncertainty()
