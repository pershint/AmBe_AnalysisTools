#Builds IBD prompt-delayed candidates from given datasets; but, once a pair is
#Found, those two events are removed from the running for other candidate pairings
import copy as cp
import math
import time
import numpy as np
from numpy import *

from sys import stdout
import glob

import os
from sys import argv

from decimal import *
import argparse

parser = argparse.ArgumentParser(description='Parser to decide what analysis to do')
parser.add_argument('--debug', dest='DEBUG',action='store_true',
        help='Run in debug mode and get some extra outputs and stuff')
parser.add_argument('--nhitpcut', dest='NPROMPT',action='store',type=int,
        help='Cut all events as a prompt candidate above this nhit')
parser.add_argument('--nhitdcut', dest='NDELAYED', action='store',type=int,
        help='Cut all events as a delayed candidate below this nhit')
parser.add_argument('--datadir', dest='DATADIR', action='store',type=str,
        help='Directory where files to analyze are stored')
parser.add_argument('--timecut', dest='TIMETHRESH', action='store',type=float,
        help='Only consider IBD candidates with an interevent time <= input')
parser.add_argument('--distcut', dest='INTERDIST', action='store',type=float,
        help='Only consider IBD candidates with an interevent_dist <= input.'+\
                'Will include events with no fitValid, but cut those with' +\
                'fitValid and failing the cut')
parser.add_argument('--fitValid', dest='FITVALID', action='store_true',
        help='Only include events that have a valid fit for IBD candidates')
parser.add_argument('--posrcut', dest='RADIUSCUT', action='store',
        help='Only consider IBD candidates with an interevent distance <= input.'+\
                'Will include events with no fitValid, but cut those with' +\
                'fitValid and failing the cut')
parser.set_defaults(DEBUG=False,NPROMPT=16, NDELAYED=6, TIMETHRESH=5.0E6, 
        INTERDIST=None, RADIUSCUT=None, FITVALID=False,
        DATADIR= "/home/onetrueteal/share/AmBe/AmBe_Control")
args = parser.parse_args()

DEBUG = args.DEBUG
NPROMPT=args.NPROMPT
NDELAYED=args.NDELAYED
TIMETHRESH=args.TIMETHRESH
INTERDIST=args.INTERDIST
RADIUSCUT=args.RADIUSCUT
FITVALID=args.FITVALID
DATADIR=args.DATADIR
import ROOT
from ROOT import TChain, TFile, gROOT

if NPROMPT < NDELAYED:
    print("Script logic needs NPROMPT > NDELAYED.  Exiting")
    sys.exit(0)

def innerDist(prev_x, prev_y, prev_z, x, y, z):
    if prev_x < -900000 or prev_y < -900000 or prev_z < -900000 or \
           x < -900000 or y < -900000 or z < -900000:
        return -999999
    else:
        return np.sqrt((prev_x - x)**2 + (prev_y - y)**2 + (prev_z - z)**2)


def deleteEventsFromBuff(i1,i2,Buffer_nhits, Buffer_times, Buffer_entrynums):
    del Buffer_nhits[i1]
    del Buffer_times[i1]
    del Buffer_entrynums[i1]
    del Buffer_nhits[i2 -1]
    del Buffer_times[i2 -1]
    del Buffer_entrynums[i2 - 1]
    return Buffer_nhits, Buffer_times, Buffer_entrynums

def deleteBuffEv(Buffer_nhits, Buffer_times, Buffer_entrynums):
    del Buffer_nhits[len(Buffer_nhits)-1]
    del Buffer_times[len(Buffer_times)-1]
    del Buffer_entrynums[len(Buffer_entrynums)-1]
    return Buffer_nhits, Buffer_times, Buffer_entrynums
#---------- TUNABLE PARAMETERS ------------#
#######DIRECTORY OF FILES#############
basepath = os.path.dirname(__file__)
FilesInBunch = 1
#######NAME OF OUTPUT FILE##########
fileN = 'IBD_Candidates_AmBeBkg_allfits'  #.ntuple.root is appended later

#######NAME OF TREE WITH DATA########
dattree = "output"
#------------- END TUNABLE PARAMETERS --------------#



def loadNewEvent(evnum, tchain, buff_nhits, buff_entrynums, buff_times):
    '''
    Using the fed-in TChain, fill buff_entrynums and get that entry's time
    using the uTNsecs and uTSecs ntuple variables'''
    haveValidev = False
    if evnum % 10000 == 0:
        print("EVNUM: " + str(evnum))
        print("BUFF_ENTRYNUMS: " + str(buff_entrynums))
        print("BUFF_TIMES: " + str(buff_times))
    while not haveValidev:
        tchain.GetEntry(evnum)
        if evnum > tchain.GetEntries():
            print("TCHAIN EXHAUSTED.  BREAKING AND CLOSING FILE.")
            break
        if FITVALID is True:
            if tchain.fitValid==0:
                evnum+=1
                continue
        if RADIUSCUT is not None:
            if tchain.fitValid==1 and (float(tchain.posr) > float(RADIUSCUT)):
                    evnum+=1
                    continue
        if tchain.nhits < min(NPROMPT,NDELAYED):
            evnum+=1
            continue
        haveValidev = True
    if haveValidev is True:
        buff_times.append(tchain.uTSecs*1.0E9 + (tchain.uTNSecs))
        buff_nhits.append(tchain.nhits)
        buff_entrynums.append(evnum)
    else:
        print("WHY THE F**K ARE WE HERE")
    #Assuming Bkg_rates is a numpy array
    return evnum, buff_nhits, buff_entrynums, buff_times

def DeleteOutOfWindow(buff_nhits,buff_entrynums, buff_times):
    allinwindow = False
    while not allinwindow:
        if len(buff_times) < 1:
            print("WE GOT A NO LENGTH BUFF IN ALLWINDOW LOOP...")
        if buff_times[len(buff_times)-1] - buff_times[0] > TIMETHRESH:
            del buff_times[0]
            del buff_nhits[0]
            del buff_entrynums[0]
        else:
            allinwindow=True
    return buff_nhits, buff_entrynums, buff_times


if __name__ == '__main__':
    '''Set up variables for output root tree'''
    uTNSecs_p        = np.zeros(1,dtype=float64)
    uTSecs_p      = np.zeros(1,dtype=float64)
    posx_p       = np.zeros(1,dtype=float64)
    posy_p        = np.zeros(1,dtype=float64)
    posz_p      = np.zeros(1,dtype=float64)
    dirx_p      = np.zeros(1,dtype=float64)
    diry_p      = np.zeros(1,dtype=float64)
    dirz_p      = np.zeros(1,dtype=float64)
    energy_p    = np.zeros(1,dtype=float64)
    nhits_p     = np.zeros(1,dtype=int)
    nhitsCleaned_p  = np.zeros(1,dtype=int)
    dcApplied_p = np.zeros(1,dtype=long)
    dcFlagged_p = np.zeros(1,dtype=long)
    runID_p     = np.zeros(1,dtype=int)
    waterFit_p  = np.zeros(1,dtype=int)
    fitValid_p  = np.zeros(1,dtype=int)
    posr_p      = np.zeros(1,dtype=float64)
    beta14_p    = np.zeros(1,dtype=float64)
    itr_p       = np.zeros(1,dtype=float64)

    uTNSecs_d        = np.zeros(1,dtype=float64)
    uTSecs_d      = np.zeros(1,dtype=float64)
    posx_d       = np.zeros(1,dtype=float64)
    posy_d        = np.zeros(1,dtype=float64)
    posz_d      = np.zeros(1,dtype=float64)
    dirx_d      = np.zeros(1,dtype=float64)
    diry_d      = np.zeros(1,dtype=float64)
    dirz_d      = np.zeros(1,dtype=float64)
    energy_d    = np.zeros(1,dtype=float64)
    nhits_d     = np.zeros(1,dtype=int)
    nhitsCleaned_d  = np.zeros(1,dtype=int)
    dcApplied_d = np.zeros(1,dtype=long)
    dcFlagged_d = np.zeros(1,dtype=long)
    runID_d     = np.zeros(1,dtype=int)
    waterFit_d  = np.zeros(1,dtype=int)
    fitValid_d  = np.zeros(1,dtype=int)
    posr_d      = np.zeros(1,dtype=float64)
    beta14_d    = np.zeros(1,dtype=float64)
    itr_d       = np.zeros(1,dtype=float64)

    interevent_dist = np.zeros(1,dtype=float64)
    interevent_time = np.zeros(1,dtype=float64)

    event_number = np.zeros(1,dtype=int)

    #Variables for filling ProcSummary
    nhit_p_cut = np.zeros(1,dtype=int)
    nhit_d_cut = np.zeros(1,dtype=int)
    allevnum = np.zeros(1,dtype=int)
    allpairnum = np.zeros(1,dtype=int)
    if INTERDIST is not None:
        interdistcut = np.zeros(1,dtype=float64)
    if RADIUSCUT is not None:
        posrcut = np.zeros(1,dtype=float64)
    #FIXME: want to save the list of files that were analyzed

    AmBe_files = glob.glob(DATADIR+"*.ntuple.root")
    print("AMBE FILES BEING USED: " + str(AmBe_files))
    #Do the IBD searches in batches of 3
    AmBe_files_split = []
    AmBe_Bunch = []
    for j,f in enumerate(AmBe_files):
        AmBe_Bunch.append(f)
        if j/FilesInBunch == float(j)/float(FilesInBunch):
            AmBe_files_split.append(AmBe_Bunch)
            AmBe_Bunch=[]
    if len(AmBe_Bunch)>0:
        AmBe_files_split.append(AmBe_Bunch)

    for j,fileBunch in enumerate(AmBe_files_split):
        print("MAKING IBD CANDIDATES FOR FILEBUNCH %s" % (str(fileBunch)))
        AmBeChain = ROOT.TChain(dattree)
        for f in fileBunch:
            AmBeChain.Add(f)
        '''Open a root file with name of dataType'''
    
        f_root = ROOT.TFile("%s_%i.ntuple.root"%(fileN,j),"recreate")
        
        m_root = ROOT.TTree("ProcSummary","Cuts applied and some meta information")
        m_root.Branch('nhit_p_cut', nhit_p_cut, 'nhit_p_cut/I')
        m_root.Branch('nhit_d_cut', nhit_d_cut, 'nhit_d_cut/I')
        m_root.Branch('allevnum', allevnum, 'allevnum/I')
        m_root.Branch('allpairnum', allpairnum, 'allpairnum/I')
        if INTERDIST is not None:
            m_root.Branch('interdistcut', interdistcut, 'interdistcut/D')
        if RADIUSCUT is not None:
            m_root.Branch('posrcut', posrcut, 'posrcut/D')
        nhit_p_cut[0] = NPROMPT
        nhit_d_cut[0] = NDELAYED
        if INTERDIST is not None:
            interdistcut[0] = INTERDIST
        if RADIUSCUT is not None:
            posrcut[0] = RADIUSCUT
    
        t_root = ROOT.TTree("CombinedOutput","Combined Prompt & Delayeds in AmBe Data")
        t_root.Branch('uTNSecs_p',       uTNSecs_p,    'uTNSecs_p/D')
        t_root.Branch('uTSecs_p',       uTSecs_p,    'uTSecs_p/D')
        t_root.Branch('posz_p',      posz_p,   'posz_p/D')
        t_root.Branch('posy_p',      posy_p,   'posy_p/D')
        t_root.Branch('posx_p',      posx_p,   'posx_p/D')
        t_root.Branch('posr_p',      posr_p,   'posr_p/D')
        t_root.Branch('dirz_p',      dirz_p,   'dirz_p/D')
        t_root.Branch('diry_p',      diry_p,   'diry_p/D')
        t_root.Branch('dirx_p',      dirx_p,   'dirx_p/D')
        t_root.Branch("energy_p",    energy_p, 'energy_p/D')
        t_root.Branch("nhits_p",     nhits_p,  'nhits_p/I')
        t_root.Branch("nhitsCleaned_p", nhitsCleaned_p, 'nhitsCleaned_p/I')
        t_root.Branch("dcApplied_p", dcApplied_p,   'dcApplied_p/L')
        t_root.Branch("dcFlagged_p", dcFlagged_p,   'dcFlagged_p/L')
        t_root.Branch("runID_p", runID_p,   'runID_p/I')
        t_root.Branch("waterFit_p", waterFit_p, 'waterFit_p/I')
        t_root.Branch("fitValid_p", fitValid_p, "fitValid_p/I")
        t_root.Branch("beta14_p", beta14_p, "beta14_p/D")
        t_root.Branch("itr_p", itr_p, "itr_p/D")
    
        t_root.Branch('uTNSecs_d',       uTNSecs_d,    'uTNSecs_d/D')
        t_root.Branch('uTSecs_d',       uTSecs_d,    'uTSecs_d/D')
        t_root.Branch('posz_d',      posz_d,   'posz_d/D')
        t_root.Branch('posy_d',      posy_d,   'posy_d/D')
        t_root.Branch('posx_d',      posx_d,   'posx_d/D')
        t_root.Branch('posr_d',      posr_d,   'posr_d/D')
        t_root.Branch('dirz_d',      dirz_d,   'dirz_d/D')
        t_root.Branch('diry_d',      diry_d,   'diry_d/D')
        t_root.Branch('dirx_d',      dirx_d,   'dirx_d/D')
        t_root.Branch("energy_d",    energy_d, 'energy_d/D')
        t_root.Branch("nhits_d",     nhits_d,  'nhits_d/I')
        t_root.Branch("nhitsCleaned_d", nhitsCleaned_d, 'nhitsCleaned_d/I')
        t_root.Branch("dcApplied_d", dcApplied_d,   'dcApplied_d/L')
        t_root.Branch("dcFlagged_d", dcFlagged_d,   'dcFlagged_d/L')
        t_root.Branch("runID_d", runID_d,   'runID_d/I')
        t_root.Branch("waterFit_d", waterFit_d, 'waterFit_d/I')
        t_root.Branch("fitValid_d", fitValid_d, "fitValid_d/I")
        t_root.Branch("beta14_d", beta14_d, "beta14_d/D")
        t_root.Branch("itr_d", itr_d, "itr_d/D")
    
        t_root.Branch('interevent_time', interevent_time,  'interevent_time/D')
        t_root.Branch('interevent_dist', interevent_dist,  'interevent_dist/D')
        t_root.Branch('pair_number', event_number, 'pair_number/I')
    
        #As we go through the AmBe files, store the entrynums and times.  If
        #We get two entries w/in TIMETHRESH, store them as a single IBD pair.
        Buffer_entrynums = []
        Buffer_times = []
        Buffer_nhits = []
        #The following continues as long as the selected file still has entrys left
        pairnum = 0
        indexnum = 0
        num_noorder = 0
        while (pairnum < 100000000):
            if DEBUG is True:
                print("CURRENT EV IN BUFFER: " + str(len(Buffer_times)))
            if float(pairnum) / 20000.0 == int(pairnum / 20000.0):
                print("ENTRYNUM: " + str(pairnum))
    
            #load a new event into our buffer
            indexnum, Buffer_nhits, Buffer_entrynums, Buffer_times = loadNewEvent(indexnum, AmBeChain,\
                    Buffer_nhits, Buffer_entrynums, Buffer_times)
            if indexnum >= AmBeChain.GetEntries():
                break
            indexnum+=1
    
            #If no entries or one entry, 
            if len(Buffer_entrynums) <= 1:
                continue
    
            #Remove events at start of buffer outside the time width range
            Buffer_nhits, Buffer_entrynums, Buffer_times = \
                    DeleteOutOfWindow(Buffer_nhits, Buffer_entrynums, Buffer_times)
            delayedindex = len(Buffer_entrynums) - 1
            #checkForOld(Buffer_entrynums, Buffer_times)
            invaliddel = False
            outoforder = False
            for i in xrange(delayedindex):
                if Buffer_times[delayedindex] - Buffer_times[i] < 0:
                    num_noorder+=1
                    outoforder = True
                    break
                if Buffer_times[delayedindex] - Buffer_times[i] < TIMETHRESH:
                    #Found a pair
                    AmBeChain.GetEntry(Buffer_entrynums[i])
                    if AmBeChain.nhits < NPROMPT:
                        continue
                    nhits_p[0]     = AmBeChain.nhits
                    fitValid_p[0]  = AmBeChain.fitValid
                    uTNSecs_p[0]        = AmBeChain.uTNSecs
                    uTSecs_p[0]      = AmBeChain.uTSecs
                    posx_p[0]       = AmBeChain.posx
                    posy_p[0]        = AmBeChain.posy
                    posz_p[0]      = AmBeChain.posz
                    posr_p[0]        = AmBeChain.posr
                    dirx_p[0]      = AmBeChain.dirx
                    diry_p[0]      = AmBeChain.diry 
                    dirz_p[0]      = AmBeChain.dirz 
                    energy_p[0]    = AmBeChain.energy 
                    nhitsCleaned_p[0]  = AmBeChain.nhitsCleaned 
                    dcApplied_p[0] = AmBeChain.dcApplied 
                    dcFlagged_p[0] = AmBeChain.dcFlagged 
                    runID_p[0]     = AmBeChain.runID 
                    waterFit_p[0]  = AmBeChain.waterFit
                    beta14_p[0]  = AmBeChain.beta14
                    itr_p[0]  = AmBeChain.itr
    
                    AmBeChain.GetEntry(Buffer_entrynums[delayedindex])
    
                    #First, fill only variables to check delayed and interev cuts
                    if AmBeChain.nhits < NDELAYED:
                        invaliddel = True
                        break
                    nhits_d[0]     = AmBeChain.nhits 
                    waterFit_d[0]  = AmBeChain.waterFit
                    posx_d[0]       = AmBeChain.posx
                    posy_d[0]        = AmBeChain.posy
                    posz_d[0]      = AmBeChain.posz
                    fitValid_d[0]  = AmBeChain.fitValid
                    uTNSecs_d[0]        = AmBeChain.uTNSecs
                    uTSecs_d[0]      = AmBeChain.uTSecs
                    posr_d[0]        = AmBeChain.posr
                    dirx_d[0]      = AmBeChain.dirx
                    diry_d[0]      = AmBeChain.diry 
                    dirz_d[0]      = AmBeChain.dirz 
                    energy_d[0]    = AmBeChain.energy 
                    nhitsCleaned_d[0]  = AmBeChain.nhitsCleaned 
                    dcApplied_d[0] = AmBeChain.dcApplied 
                    dcFlagged_d[0] = AmBeChain.dcFlagged 
                    runID_d[0]     = AmBeChain.runID 
                    beta14_d[0]  = AmBeChain.beta14
                    itr_d[0]  = AmBeChain.itr
    
    
                    if waterFit_p[0] == 0 or waterFit_d[0] == 0:
                        interevent_dist[0] = 0
                    else:
                        interevent_dist[0] = innerDist(posx_p[0], posy_p[0],
                            posz_p[0],posx_d[0],posy_d[0],posz_d[0])
                    if INTERDIST is not None:
                        if fitValid_p[0] == 1 and fitValid_d[0] == 1 and \
                                float(interevent_dist[0]) > float(INTERDIST):
                            continue
                    interevent_time[0] = Buffer_times[delayedindex] - Buffer_times[i]
 
                    #nsp = cp.deepcopy(uTNSecs_p[0])
                    #sp = cp.deepcopy(uTSecs_p[0])
                    #nsd = cp.deepcopy(uTNSecs_d[0])
                    #sd = cp.deepcopy(uTSecs_d[0])
                    #if sp == sd:
                    #    interevent_time[0] = nsd - nsp
                    #elif (sd - sp) > 0:
                    #    interevent_time[0] = (sd-sp)*1.0E9 + nsd - nsp
                    #else:
                    #    continue
                    event_number[0] = pairnum
                    t_root.Fill()
                    found_pair = True
                    pairnum+=1
                    deleteEventsFromBuff(i,delayedindex,Buffer_nhits,
                            Buffer_times,Buffer_entrynums)
                    break
                    #now, delete the most current event if it's not above
                    #prompt threshold
            if outoforder is True:
                print("OUT OF ORDER EVENT: CLEARING BUFFER...")
                Buffer_times = []
                Buffer_entrynums = []
                Buffer_nhits = []
                continue
            if len(Buffer_nhits)>0 and int(Buffer_nhits[len(Buffer_nhits)-1]) < NPROMPT:
                Buffer_nhits, Buffer_times, Buffer_entrynums = \
                        deleteBuffEv(Buffer_nhits, Buffer_times, Buffer_entrynums)
            if invaliddel is True:
                Buffer_nhits, Buffer_times, Buffer_entrynums = \
                        deleteBuffEv(Buffer_nhits, Buffer_times, Buffer_entrynums)
        allevnum[0] = indexnum
        allpairnum[0] = pairnum
        m_root.Fill()
        f_root.cd()
        t_root.Write()
        m_root.Write()
        f_root.Close()
        print("NUMBER OUT OF ORDER: " + str(num_noorder))

