import multiprocessing as mp
import numpy as np
import subprocess
import concurrent.futures
import os
from integrations import long_integration, short_integration, kozai_check


# make function needed for multiprocessing
def do_integration(intN):

        """ when running this code, change the parameters of long_integration"""
        # first do long integration 
        # arguments are i, minA, maxA, minE, maxE, maxI, totalTime (exponent), Nparticles, max_exit_distance
        # returns minA, maxA, minE, maxE, maxI, Nparticles, totalTime, filename

        amin, amax, emin, emax, imax, maxD, nParticles, totTime, filname = long_integration(intN, 38.81, 39.95, 0.0, 0.6, 35.0, 9, 5000, 85)

        #now do short integration
        # arguments:i, minA, maxA, minE, maxE, maxI, shortTime, fileName, snapshotSlice (0 = time 0, -2 = second to last, -1 = last)     
        # return sim right now
        # we want to do 3 short integrations starting with time 0, 10%totalTime, and totalTime

        #amin = 38.81
        #amax = 39.81
        #emin = 0.0
        #emax = 0.6
        #imax = 35
        #maxD = 85
        #filname = 'Jul162020.22.27_part500_A_38.810-39.810_Q_15.524-63.696_I_0-0.611_E_0.000-0.600_even_q_{}'.format(intN)
        shortSim0 = short_integration(intN, amin, amax, emin, emax, imax, maxD, 1e5, filname, 0)
        #shortSim3 = short_integration(intN, amin, amax, emin, emax, imax, maxD, 1e5, filname, -3) #pay attention to these indices
        shortSim1 = short_integration(intN, amin, amax, emin, emax, imax, maxD, 1e5, filname, -2)
        shortSim2 = short_integration(intN, amin, amax, emin, emax, imax, maxD, 1e5, filname, -1)
        kozaiShortInt1 = kozai_check(intN, amin, amax, emin, emax, imax, maxD, 1e7, filname, 0)
        kozaiShortInt2 = kozai_check(intN, amin, amax, emin, emax, imax, maxD, 1e7, filname, -1) 
        kozaiShortInt3 = kozai_check(intN, amin, amax, emin, emax, imax, maxD, 1e7, filname, -2) 

        return filname


#multiprocessing 
#when we run this code, change this number. It determines how many parallel processes are happening
lenInt = 1

args = np.arange(lenInt)

with mp.Pool(lenInt) as pool:
    results = pool.map(do_integration,args)

