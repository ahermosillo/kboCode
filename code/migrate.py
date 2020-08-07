import multiprocessing as mp
import numpy as np
import subprocess
import concurrent.futures
import os
from integrations import long_integration, short_integration



# create files that we are going to append data to
# that is done here



# run integrations
# this works but they run once the other one is done, so using multiprocess instead
# p0 = subprocess.run(['python', 'do_integration.py'])
# p1 = subprocess.run(['python', 'do_integration.py'])
# p2 = subprocess.run(['python', 'do_integration.py'])
# p3 = subprocess.run(['python', 'do_integration.py'])
# p4 = subprocess.run(['python', 'do_integration.py'])
# p5 = subprocess.run(['python', 'do_integration.py'])

#try multiprocessing

# make function needed for multiprocessing
def do_integration(intN):
        # first do long integration 
        # arguments are i, minA, maxA, minE, maxE, maxI, totalTime, intervalTime, Nparticles
        # returns minA, maxA, minE, maxE, maxI, Nparticles, totalTime, filename
        amin, amax, emin, emax, imax, nParticles, totTime, filname = long_integration(intN, 38.81, 39.81, 0.0, 0.6, 35.0, 6, 20)

        #now do short integration
        # arguments:i, xminA, maxA, minE, maxE, maxI, nParticles, shortTime, fileName, snapshotSlice (0 = time 0, -2 = second to last, -1 = last)
        # return sim right now
        # we want to do 3 short integrations starting with time 0, 10%totalTime, and totalTime
        shortSim0 = short_integration(intN, amin, amax, emin, emax, imax, 1e5, filname, 0)
        shortSim1 = short_integration(intN, amin, amax, emin, emax, imax, 1e5, filname, -2)
        shortSim2 = short_integration(intN, amin, amax, emin, emax, imax, 1e5, filname, -1)

        return filname


#multiprocessing 
lenInt = 1

args = np.arange(lenInt)

with mp.Pool(lenInt) as pool:
    results = pool.map(do_integration,args)

