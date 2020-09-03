import rebound
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import csv
import pdb
import glob
import time 
#import tqdm    
import os as osp       
import shutil
import itertools
from rebound import hash as h
from ctypes import c_uint32



def moveAndCreateDir(sPath, dstDir):
    if not osp.path.isdir(dstDir):
        osp.makedirs(dstDir)
        shutil.move(sPath,dstDir)
        #print("created directory {} and moved it to {}".format(sPath, dstDir))
    else:
    	# print("directory already exists")
    	pass

def prob(i):
    sigma_i = 12.0*np.pi/180.  # average 15 degrees
    prob_i = np.sin(i)*np.exp(-i**2/(2.0*sigma_i**2))
    return(prob_i)

def get_i():
    keep = False
    
    while keep == False:
        num = np.random.uniform()*np.pi  # pick an inclination
        
        if np.random.uniform() < prob(num):  # decide whether to keep the inclination
            keep = True

    return(num)

def long_integration(i, minA, maxA, minE, maxE, maxI, totalExpTime, Nparticles, maxDistance):
    

    """runs an nbody integrationfor the giant planets + Nparticles test particles. 
    It goes for totalTime year
    it returns minA, maxA, minE, maxE and creates a bin file with the simulation archive
    """
    
    start_time = time.time()
    tm = time.gmtime()
    dat= time.strftime("%b%d%Y.%H.%M", tm)        

    base = 30*(3./2.)**(2./3.) #~39AU

    minQ = (minA*(1.-maxE))            #Perihelion distance
    maxQ = (maxA*(1.+maxE))            #Apehelion distance
    maxI =  maxI*(np.pi/180)           #maximum inclination argument 

    sim = rebound.Simulation()
    sim.units = ('yr', 'AU', 'Msun')
    print("starting integration {}".format(i))

    totalTime = 10**totalExpTime

    print("set units as yr, AU, Msun")
    print("Integration is running for {} years for giant planets plus {} test particles".format(totalTime, Nparticles))

    #this is a file with m, r, x, y, z, vx, vy, vz for sun and giant planets
    # they are barycenter values from NASA Horizons
    with open('planetParamsCartesianAU_yr.txt', 'r') as file:
        data = file.read()

    sim.add_particles_ascii(data)

    print("---added sun, jupiter, saturn, uranus, neptune---")


    sun=sim.particles[0]
    jupiter=sim.particles[1]
    saturn=sim.particles[2]
    uranus=sim.particles[3]
    neptune=sim.particles[4]

    sun.hash = 1
    jupiter.hash = 2
    saturn.hash = 3
    uranus.hash = 4
    neptune.hash = 5

    np.random.seed(i)

    """
    #even e range
    for p in range(Nparticles):
        sim.add(a=np.random.uniform(minA, maxA), e=np.random.uniform(minE,maxE), inc=np.random.uniform(0,maxI), omega=np.random.rand()*(2.*np.pi),Omega=np.random.rand()*(2.*np.pi),M=np.random.rand()*(2.*np.pi))
    """

    #even q range
    for p in range(Nparticles):
        sem = np.random.uniform(minA,maxA)
        ecc = 1-(np.random.uniform(15.5, sem))/sem 
        sim.add(a=sem, e=ecc, inc= get_i(), omega=np.random.uniform(0,2*np.pi),Omega=np.random.uniform(0, 2*np.pi),M=np.random.uniform(0,2*np.pi), hash = p+6) 

    sim.dt = 0.2
    sim.integrator = 'whfast'
    #sim.ri_whfast.safe_mode = 0
    #sim.ri_whfast.corrector = 11 
    sim.move_to_com()            # move particles to center of momentum frame
    sim.N_active = 5             # number of active massive particles, includes Sun
    sim.exit_max_distance = maxDistance
    #sim.status()

    #the following code should be set up on first use to locate and store your simulations
    sourcePath = '/data/galadriel/ahr/kboTest/code/' 
    destinPath = '/data/galadriel/ahr/kboTest/Long_Integrations/'

    #def heartbeat(sim):         # verification that sim is running
    #   print(sim.contents.t)
    #sim.heartbeat = heartbeat

    filename = '{}_part{}_time{}_A_{:.3f}-{:.3f}_Q_{:.3f}-{:.3f}_I_0-{:.3f}_E_{:.3f}-{:.3f}_even_q_{}'.format(dat,Nparticles,totalTime,minA,maxA,minQ,maxQ,maxI,minE,maxE,i)


    #--------Automated Simulation Archive--------
    #sim.automateSimulationArchive('{}.bin'.format(filename), interval = intervalTime, deletefile = True)
    #sim.integrate(totalTime, exact_finish_time = 0)

    # ------- Manual Snapshots --------------------------

    realtimes = np.logspace(1, totalExpTime, totalExpTime)
    realtimes = np.insert(realtimes, 0, 0)

    EscParticles = []
    for i,t in enumerate(realtimes):
        while sim.t < t:
            try:
                #print('Integrating for {} years'.format(t))
                #sim.ri_whfast.recalculate_coordinates_this_timestep = 1
                sim.integrate(t, exact_finish_time=0)
                #print(sim.t)
              
            except rebound.Escape as error:
                #print('A particle escaped at {} years'.format(sim.t))
                for j in range(sim.N):
                    p = sim.particles[j]
                    dist = p.x*p.x + p.y*p.y + p.z*p.z
                    if dist > sim.exit_max_distance**2:
                        index = p.hash
                sim.remove(hash=index)
                #sim.ri_whfast.recalculate_coordinates_this_timestep = 1
                EscParticles.append(index)
                #sim.integrator_synchronize()
        print("done: {}".format(t))
        print("time passed: {}".format(time.time() - start_time))
        sim.simulationarchive_snapshot('{}.bin'.format(filename))
    #print("these particles escaping {}".format(EscParticles))

    print("long integration is done")


    sourceFile = '{}{}.bin'.format(sourcePath,filename)
    destDir    = '{}{}'.format(destinPath,filename)
    moveAndCreateDir(sourceFile,destDir)

    data = np.transpose([minA, maxE, minE, maxE, maxI, Nparticles, totalTime, filename])
   # np.savetxt("{}_{}.txt".format(dat, i), data, fmt = '%s')
    
    print("long Integration took {}".format(time.time() - start_time))

    return minA, maxA, minE, maxE, maxI, maxDistance, Nparticles, totalTime, filename



def short_integration(integrationN, minA, maxA, minE, maxE, maxI, maxDistance, shortTime, fileName, indexSimulation):

    start_simul_time = time.time()

    #the following code should be set up on first use to locate and store your simulations

    destinPath = '/data/galadriel/ahr/kboTest/Long_Integrations/'

    longInt    = '{}{}/{}'.format(destinPath,fileName,fileName)

    sa = rebound.SimulationArchive("{}.bin".format(longInt)) #names archive object 
    sim = sa[indexSimulation] ## see comment above for this 
    ST = sim.t             #(the snapshot time we're using as initial conditions)
    # doing this just for naming purposes

    sTemp    = 'TemporaryDirectory_time_{}'.format(np.round(ST))
    iRes     = 'In_Resonance'
    nRes     = 'Not_In_Resonance'
    sInt     = 'Short_Integrations'

    mainDir    = '{}{}/'.format(destinPath,fileName)
    dirName    = '{}{}/{}'.format(destinPath,fileName,sInt)
    subDirTemp = '{}{}/{}/{}'.format(destinPath,fileName,sInt,sTemp)
    irDir      = '{}{}/{}/{}/{}'.format(destinPath,fileName,sInt,sTemp,iRes)
    nrDir      = '{}{}/{}/{}/{}'.format(destinPath,fileName,sInt,sTemp,nRes)

    try: 
        osp.mkdir(destinPath)
        #print('Directory',destinPath,'created.')
    except FileExistsError: 
        pass
        #print('Directory',destinPath,'already exists.')    

    try:
        osp.mkdir(dirName)
        #print ('Directory',dirName,'created.')
    except FileExistsError:
        pass
        #print("Directory",dirName,"already exists.")

    try:
        osp.mkdir(subDirTemp)
        #print ('Directory',subDirTemp,'created.')
    except FileExistsError:
        pass
        #print("Directory",subDirTemp,"already exists.")

    try:
        osp.mkdir(irDir)
        #print ('Directory',irDir,'created.')
    except FileExistsError:
        pass
        #print("Directory",irDir,"already exists.")

    try:
        osp.mkdir(nrDir)
        #print ('Directory',nrDir,'created.')
    except FileExistsError:
        pass
        #print("Directory",nrDir,"already exists.")


    IT = shortTime          #this is the short integration run time.
    ET = ST + IT
    Nout = 1000
    Ntotal = sim.N -1 #number of objects without the sun
    npart = sim.N - sim.N_active # number of test particles
    sim.exit_max_distance = maxDistance
    #sim.ri_whfast.safe_mode = 0
    #sim.ri_whfast.corrector = 11 

    # setting hashes to the particles just for now since doing 1e9 integration
    ### TAKE THIS OFF LATER ############################################3

    #for partN in range(sim.N):
     #   sim.particles[partN].hash = partN + 1

    ### TAKE THIS OFF LATER. JUST HAVE IT FOR THE 1E9 LONG INT THAT RAN BEFORE REMOVAL OF PARTICLES ###########

    print("starting short integration {} with start time {}".format(integrationN, ST))
    
    ecc = np.zeros(npart)
    ax = np.zeros(npart)
    inc = np.zeros(npart)
    
    for i in range(npart):
        ecc[i] = sim.particles[i+5].e
        ax[i] = sim.particles[i+5].a
        inc[i] = sim.particles[i+5].inc
        
        
    arange = np.linspace(37,50)
    fig, (p1, p2) = plt.subplots(1,2, figsize = (11,5))
    p1.plot(ax, ecc, '.')
    p2.plot(ax, inc, '.')
    p1.plot(arange, 1-20/arange , '--', label = "r = 20AU")
    p1.plot(arange, 1-25/arange , '--', label = "r = 25AU")
    p1.plot(arange, 1-30/arange, '--', label = "r = 30AU")
    p1.plot(arange, 1-35/arange, '--', label = "r = 35AU")
    p1.vlines(minA,0, 1.1)
    p1.vlines(maxA, 0,1.1)
    p1.hlines(0.6, 37, 50)
    p1.set_xlim(37, 50)
    p1.set_ylim(0,1.1)
    p1.set_title("a v e; t = {}".format(ST))
    p1.set_xlabel("a [AU]")
    p1.set_ylabel("e ")
    p2.vlines(minA, 0, np.pi/2)
    p2.vlines(maxA, 0, np.pi/2)
    p2.set_xlim(37, 50)
    p2.set_ylim(0,np.pi/2)
    p2.set_title("a v i; t = {}".format(ST))
    p2.set_xlabel("a [AU]")
    p2.set_ylabel("i")
    p1.legend()
    plt.savefig('{}/ave_avi_plot'.format(subDirTemp))

    # pointer to get the hashes of particles still alive

    hashesParticles = np.zeros(sim.N,dtype="uint32")
    sim.serialize_particle_data(hash=hashesParticles)
    #print(hashesParticles)

    intTimes = np.linspace(ST, ET, Nout)

    ### ------------ where we are going to store values ----------------- #####
    # these are 2D matrices that store the values for each particle, for each timestep
    # made them have 9999 so it's noticeable where particles escaped (values didn't store)

    l = np.zeros((Ntotal,Nout))
    p = np.zeros((Ntotal,Nout))
    a = np.ones((Ntotal,Nout))*9999
    e = np.ones((Ntotal,Nout))*9999
    lasc_node = np.ones((Ntotal,Nout))*9999
    arg_peri = np.ones((Ntotal,Nout))*9999
    t_anom = np.ones((Ntotal,Nout))*9999
    incl = np.ones((Ntotal,Nout))*9999
    phi = np.zeros((Ntotal,Nout))
    M_anom = np.ones((Ntotal,Nout))*9999
    xvals = np.ones((Ntotal,Nout))*9999
    yvals = np.ones((Ntotal,Nout))*9999
    timeArray = np.ones(Nout)*9999

    mln_arr = []
    mlp_arr = []
    pj_arr = []
    a_arr = []
    e_arr = []
    inc_arr = []
    lasc_node_arr = []
    arg_peri_arr = []
    t_anom_arr = []
    M_anom_arr = []
    x_arr = []
    y_arr = []

    print(sim.exit_max_distance)


    for i,times in enumerate(intTimes):
        while sim.t < times:
            try:
                #sim.ri_whfast.recalculate_coordinates_this_timestep = 1
                sim.integrate(times, exact_finish_time = 0)
            except rebound.Escape as error:
                #print(error)
                for part in range(sim.N):
                    psim = sim.particles[part]
                    dist = psim.x*psim.x + psim.y*psim.y + psim.z*psim.z
                    if dist > sim.exit_max_distance**2:
                        index = psim.hash
    #                     print("particle with hash {} escaped".format(index))
                sim.remove(hash=index)
                #sim.ri_whfast.recalculate_coordinates_this_timestep = 1
        #sim.integrator_synchronize()
        for j, j_hash in enumerate(hashesParticles[1:]):
    #         print("j: {} ; hash: {}".format(j, j_hash))
            try: 
                ps = sim.particles[h(c_uint32(j_hash))]
                l[j][i] = ps.calculate_orbit().l
                p[j][i] = ps.calculate_orbit().pomega
                a[j][i] = ps.calculate_orbit().a
                e[j][i] = ps.calculate_orbit().e
                incl[j][i] = ps.calculate_orbit().inc
                lasc_node[j][i] = ps.calculate_orbit().Omega 
                arg_peri[j][i] = ps.calculate_orbit().omega
                t_anom[j][i] = ps.calculate_orbit().f
                M_anom[j][i] = ps.calculate_orbit().M
                xvals[j][i] = ps.x
                yvals[j][i] = ps.y
            except rebound.ParticleNotFound as error: 
                # since particles escaping as we store/integrate
                pass
    #             print("idk {}".format(error))

            #renaming values
            mlp = l[j][i]
            pj = p[j][i]
            mln = l[3][i]
            sem = a[j][i]
            ecc = e[j][i]
            inc = incl[j][i]
            lan = lasc_node[j][i]
            ap = arg_peri[j][i]
            ta = t_anom[j][i]
            ma = M_anom[j][i]
            x = xvals[j][i]
            y = yvals[j][i]


            #appending to cleaned up arrays
            mln_arr.append(mln)
            mlp_arr.append(mlp)
            pj_arr.append(pj)
            a_arr.append(sem)
            e_arr.append(ecc)
            inc_arr.append(inc)
            lasc_node_arr.append(lan)
            arg_peri_arr.append(ap)
            t_anom_arr.append(ta)
            M_anom_arr.append(ma)
            x_arr.append(x)
            y_arr.append(y)

            phi_temp = 3.*mlp - 2.*mln - pj   
            phi[j][i] = phi_temp%(2*np.pi)

    print("done: after short int {} particles left".format(sim.N))


    resonant_particles = []
    count = 0        
    for i in range(Ntotal):
        phi_i = phi[i]
        if (all(phi[i] < 355*np.pi/180) and all(phi[i] > 5*np.pi/180)):
            #print(phi_i)
            resonant_particles.append(i)
            count +=1

            #plt.figure(figsize=(15,10))
            #plt.title('Resonant angle libration', fontsize = 24)
            #plt.xlabel('Time(years)', fontsize = 18)
            #plt.ylabel('Resonant argument (degrees)', fontsize = 18)
            #plt.scatter(intTimes,phi_i, marker = '.',s = 10)
            #plt.ylim(0, 2*np.pi)
            #plt.savefig('{}/Particle {} Phi vs Time Plot.png'.format(irDir,i))  
            #plt.clf()
        with open("{}/Particles_in_resonance_{}.txt".format(subDirTemp, np.round(ST)), "w+") as my_file:                               
            for i in resonant_particles:
                 my_file.write(str(i)+"\n")



    nonresonant_particles = []
    count_n = 0        
    for j in range(Ntotal):
        phi_n = phi[j]
        if (any(phi[j] > 355*np.pi/180) and any(phi[j] < 5*np.pi/180)):
            #print(phi_n)
            nonresonant_particles.append(j)
            count_n +=1

            #plt.figure(figsize=(15,10))
            #plt.title('Resonant angle circulation', fontsize = 24)
            #plt.xlabel('Time(years)', fontsize = 18)
            #plt.ylabel('Resonant argument (degrees)', fontsize = 18)
            #plt.scatter(intTimes,phi_n, marker = '.',s = 10)
            #plt.ylim(0,2*np.pi)
            #plt.savefig('{}/Particle {} Phi vs Time Plot.png'.format(nrDir,i))  
            #plt.clf()
        with open("{}/Particles_not_in_resonance_{}.txt".format(subDirTemp, np.round(ST)), "w+") as my_file:                               
            for i in nonresonant_particles:
                my_file.write(str(i)+"\n")



    deltaTheta = np.zeros((Ntotal, Nout))
    deltaThetaArr = []

    rotate_diff = []
    for i in range(Nout):
        for j in range(Ntotal):
            Neptune_xsim = xvals[3][i]
            Neptune_ysim = yvals[3][i]

            xsim = xvals[j][i]
            ysim = yvals[j][i]
            theta_sim = np.arctan2(ysim,xsim)


            #xsur = 29.06807239827766   
            #ysur = -7.125912195043246
            xsur = 26.85758046958696
            ysur = -13.32890006819031
            # for julian date ? 
            theta_sur = np.arctan2(ysur,xsur)

            #difference of both angles
            new_theta = theta_sim - theta_sur
            deltaTheta[j][i] = new_theta
            deltaThetaArr.append(new_theta)
            #new theta is added to the lasc of all the giant planets and the test particles

            rotate_diff.append(new_theta)

    #Applying rotation to array of lasc values by adding rotate_diff values to corresponding values in rotated_longitude array
    rotated_longitude = np.zeros((Ntotal, Nout))

    for i in range(len(lasc_node[0])):
        for j in range(len(lasc_node)):
            rotated_longitude[j][i] = lasc_node[j][i] - deltaTheta[3][i]





    #data array
    data_arr = []

    for particle in range(len(l)):
        data_arr.append([])   
        for integration in range(len(l[0])):
            #print("Integration: " + str(integration) + "Num: " + str(num))
            data_arr[particle].append(hashesParticles[particle])
            data_arr[particle].append(a[particle][integration]) 
            data_arr[particle].append(e[particle][integration])
            data_arr[particle].append(incl[particle][integration])
            data_arr[particle].append(lasc_node[particle][integration])
            data_arr[particle].append(arg_peri[particle][integration])
            data_arr[particle].append(t_anom[particle][integration])
            data_arr[particle].append(phi[particle][integration])
            data_arr[particle].append(deltaTheta[particle][integration])
            data_arr[particle].append(xvals[particle][integration])
            data_arr[particle].append(yvals[particle][integration])
            if particle in resonant_particles:
                data_arr[particle].append(True)
            else:
                data_arr[particle].append(False)   

    data_arr = np.array(data_arr).reshape(len(l)*len(l[0]), 12) # (numParticles*numTimesteps, numOutputs)     

    with open('{}/{}_data_array.csv'.format(subDirTemp,fileName), mode = 'w') as file:
        datawriter = csv.writer(file, delimiter = ',')
        datawriter.writerow(['pnumber', 'a', 'e', 'i', 'Omega', 'w', 'f', 'phi', 'dTheta', 'x', 'y', 'resonance'])
        for d in data_arr:
            datawriter.writerow(d)



    A = []
    E = []
    I = []
    N = []
    P = []
    M = []
    F = []


    X = []
    Y = []
    

    for semimajor in a[resonant_particles]:
        A.append(semimajor)
    for eccentricity in (e[resonant_particles]):
        E.append(eccentricity)
    for inclination in (incl[resonant_particles]):
        I.append(inclination)
    for longRot in rotated_longitude[resonant_particles]:
        N.append(longRot) 
    for argument in (arg_peri[resonant_particles]):
        P.append(argument)    
    for mean in (M_anom[resonant_particles]):
        M.append(mean)
    for trueAn in (t_anom[resonant_particles]):
        F.append(trueAn)

    for val1 in xvals[resonant_particles]:
        X.append(val1)
    for val2 in yvals[resonant_particles]:
        Y.append(val2)


    A_iter = itertools.chain.from_iterable(A)
    E_iter = itertools.chain.from_iterable(E)
    I_iter = itertools.chain.from_iterable(I)
    N_iter = itertools.chain.from_iterable(N)
    P_iter = itertools.chain.from_iterable(P)
    M_iter = itertools.chain.from_iterable(M)
    F_iter = itertools.chain.from_iterable(F)

    X_iter = itertools.chain.from_iterable(X)
    Y_iter = itertools.chain.from_iterable(Y)


    A1 = list(A_iter)
    E1 = list(E_iter)
    I1 = list(I_iter)
    N1 = list(N_iter)
    P1 = list(P_iter)
    M1 = list(M_iter)
    F1 = list(F_iter)

    X1 = list(X_iter)
    Y1 = list(Y_iter)

    data=np.transpose([A1,E1,I1,N1,P1,M1])
    np.savetxt('{}/Survey_data_{}'.format(subDirTemp, integrationN),data,fmt='%s')

    cartesian_data=np.transpose([A1, E1, I1, N1, P1, F1, X1,Y1])
    np.savetxt('{}/data_for_checking_resonance_{}'.format(subDirTemp, integrationN),cartesian_data,fmt='%s')


    #Writing a,e, and i values to read in Plotting script

    res_plot_data = []
    nonres_plot_data = []    

    ar = []
    er = []
    ir = []
    for i in a[resonant_particles]:
        ar.append(i[0])
    for j in e[resonant_particles]:
        er.append(j[0])
    for k in incl[resonant_particles]:
        ir.append(k[0])  
    res_plot_data.append(ar) 
    res_plot_data.append(er) 
    res_plot_data.append(ir)  
    np.savetxt('{}/res_plot_data.txt'.format(subDirTemp), res_plot_data, fmt = '%s')

    anr = []
    enr = []
    inr = []
    for i in a[nonresonant_particles]:
        anr.append(i[0])
    for j in e[nonresonant_particles]:
        enr.append(j[0])
    for k in incl[nonresonant_particles]:
        inr.append(k[0])
    nonres_plot_data.append(anr) 
    nonres_plot_data.append(enr) 
    nonres_plot_data.append(inr)  
    np.savetxt('{}/nonres_plot_data.txt'.format(subDirTemp), nonres_plot_data, fmt = '%s')

    print('DONE!')
    print("short integration {} took {}".format(integrationN, time.time() - start_simul_time))

    return sim



def kozai_check(integrationN, minA, maxA, minE, maxE, maxI, maxDistance, shortTime, fileName, indexSimulation):
    start_simul_time = time.time()

    #the following code should be set up on first use to locate and store your simulations

    destinPath = '/data/galadriel/ahr/kboTest/Long_Integrations/'
    
    longInt    = '{}{}/{}'.format(destinPath,fileName,fileName)
    print("{}.bin".format(longInt))
    sa = rebound.SimulationArchive("{}.bin".format(longInt)) #names archive object 
    sim = sa[indexSimulation] ## see comment above for this 
    ST = np.round(sim.t)             #(the snapshot time we're using as initial conditions)
    print("simulation time: {}".format(ST))
    print("particles N: {}".format(sim.N))
    # doing this just for naming purposes

    sTemp    = 'TemporaryDirectory_time_{}'.format(np.round(ST))
    iRes     = 'In_Resonance'
    nRes     = 'Not_In_Resonance'
    sInt     = 'Short_Integrations'
    kozFiles = 'Kozai_Resonance'
    iKoz     = 'In_Kozai_Resonance'
    nKoz     = 'Not_In_Kozai_Resonance'
    
    
    mainDir    = '{}{}/'.format(destinPath,fileName)
    dirName    = '{}{}/{}'.format(destinPath,fileName,sInt)
    subDirTemp = '{}{}/{}/{}'.format(destinPath,fileName,sInt,sTemp)
    irDir      = '{}{}/{}/{}/{}'.format(destinPath,fileName,sInt,sTemp,iRes)
    nrDir      = '{}{}/{}/{}/{}'.format(destinPath,fileName,sInt,sTemp,nRes)
    kozDir     = '{}{}/{}/{}/{}/{}'.format(destinPath,fileName,sInt,sTemp,iRes,kozFiles)
    irKoz      = '{}{}/{}/{}/{}/{}/{}'.format(destinPath,fileName,sInt,sTemp,iRes,kozFiles,iKoz)
    nrKoz      = '{}{}/{}/{}/{}/{}/{}'.format(destinPath,fileName,sInt,sTemp,iRes,kozFiles,nKoz)

    try: 
        osp.mkdir(destinPath)
        print('Directory',destinPath,'created.')
    except FileExistsError: 
        print('Directory',destinPath,'already exists.')    

    try:
        osp.mkdir(dirName)
        print ('Directory',dirName,'created.')
    except FileExistsError:
        print("Directory",dirName,"already exists.")

    try:
        osp.mkdir(subDirTemp)
        print ('Directory',subDirTemp,'created.')
    except FileExistsError:
        print("Directory",subDirTemp,"already exists.")

    try:
        osp.mkdir(irDir)
        print ('Directory',irDir,'created.')
    except FileExistsError:
        print("Directory",irDir,"already exists.")

    try:
        osp.mkdir(nrDir)
        print ('Directory',nrDir,'created.')
    except FileExistsError:
        print("Directory",nrDir,"already exists.")
        
        
    try:
        osp.mkdir(kozDir)
        print ('Directory',kozDir,'created.')
    except FileExistsError:
        print("Directory",kozDir,"already exists.")    
    try:
        osp.mkdir(irKoz)
        print ('Directory',irKoz,'created.')
    except FileExistsError:
        print("Directory",irKoz,"already exists.")

    try:
        osp.mkdir(nrKoz)
        print ('Directory',nrKoz,'created.')
    except FileExistsError:
        print("Directory",nrKoz,"already exists.")

    IT = shortTime          #this is the short integration run time.
    ET = ST + IT
    Ntotal = sim.N -1 #number of objects without the sun
    npart = sim.N - sim.N_active # number of test particles
    Nout = 1000
    sim.exit_max_distance = maxDistance
    t = np.linspace(ST,ET,Nout)

    hashesParticles = np.zeros(sim.N,dtype="uint32")
    sim.serialize_particle_data(hash=hashesParticles)
    #print(hashesParticles)

    print("now running short integration from {} with {} particles + giant planets + sun".format(sim.t, npart))

    intTimes = np.linspace(ST, ET, Nout)

    ### ------------ where we are going to store values ----------------- #####
    # these are 2D matrices that store the values for each particle, for each timestep
    # made them have 9999 so it's noticeable where particles escaped (values didn't store)

    l = np.zeros((Ntotal,Nout))
    p = np.ones((Ntotal,Nout))*9999
    a = np.ones((Ntotal,Nout))*9999
    e = np.ones((Ntotal,Nout))*9999
    lasc_node = np.ones((Ntotal,Nout))*9999
    arg_peri = np.ones((Ntotal,Nout))*9999
    t_anom = np.ones((Ntotal,Nout))*9999
    incl = np.ones((Ntotal,Nout))*9999
    phi = np.zeros((Ntotal,Nout))
    M_anom = np.ones((Ntotal,Nout))*9999
    xvals = np.ones((Ntotal,Nout))*9999
    yvals = np.ones((Ntotal,Nout))*9999
    timeArray = np.ones(Nout)*9999

    mln_arr = []
    mlp_arr = []
    pj_arr = []
    a_arr = []
    e_arr = []
    inc_arr = []
    lasc_node_arr = []
    arg_peri_arr = []
    t_anom_arr = []
    M_anom_arr = []
    x_arr = []
    y_arr = []

    print(sim.exit_max_distance)

    for i,times in enumerate(intTimes):
        while sim.t < times:
            try:
                sim.integrate(times, exact_finish_time = 0)
            except rebound.Escape as error:
                for part in range(sim.N):
                    psim = sim.particles[part]
                    dist = psim.x*psim.x + psim.y*psim.y + psim.z*psim.z
                    if dist > sim.exit_max_distance**2:
                        index = psim.hash
    #                     print("particle with hash {} escaped".format(index))
                sim.remove(hash=index)
        for j, j_hash in enumerate(hashesParticles[1:]):
    #         print("j: {} ; hash: {}".format(j, j_hash))
            try: 
                ps = sim.particles[h(c_uint32(j_hash))]
                l[j][i] = ps.calculate_orbit().l
                p[j][i] = ps.calculate_orbit().pomega
                a[j][i] = ps.calculate_orbit().a
                e[j][i] = ps.calculate_orbit().e
                incl[j][i] = ps.calculate_orbit().inc
                lasc_node[j][i] = ps.calculate_orbit().Omega 
                arg_peri[j][i] = ps.calculate_orbit().omega
                t_anom[j][i] = ps.calculate_orbit().f
                M_anom[j][i] = ps.calculate_orbit().M
                xvals[j][i] = ps.x
                yvals[j][i] = ps.y
            except rebound.ParticleNotFound as error: 
                # since particles escaping as we store/integrate
                pass
    #             print("idk {}".format(error))

            #renaming values
            mlp = l[j][i]
            pj = p[j][i]
            mln = l[3][i]
            sem = a[j][i]
            ecc = e[j][i]
            inc = incl[j][i]
            lan = lasc_node[j][i]
            ap = arg_peri[j][i]
            ta = t_anom[j][i]
            ma = M_anom[j][i]
            x = xvals[j][i]
            y = yvals[j][i]


            #appending to cleaned up arrays
            mln_arr.append(mln)
            mlp_arr.append(mlp)
            pj_arr.append(pj)
            a_arr.append(sem)
            e_arr.append(ecc)
            inc_arr.append(inc)
            lasc_node_arr.append(lan)
            arg_peri_arr.append(ap)
            t_anom_arr.append(ta)
            M_anom_arr.append(ma)
            x_arr.append(x)
            y_arr.append(y)

            phi_temp = 3.*mlp - 2.*mln - pj   
            phi[j][i] = phi_temp%(2*np.pi)

    #print("done: after short int {} particles left".format(sim.N))



    #change paths to wherever you want the images saved to
    #Testing for 3:2 MMR and Kozai 

    resonant_particles = []
    nonresonant_particles = []   
    kozai_particles = []
    nonkozai_particles = []
           
    for i in range(sim.N-1):
        phi_i = phi[i]*180/np.pi
        kozai = arg_peri[i]*180/np.pi
        if (all(phi[i] < 355/180*np.pi) and all(phi[i] > 5/180*np.pi)):            
            #resonant_particles.append(i)
            if ((all(kozai < 175) and all(kozai > 5)) or (all(kozai < 355) and all(kozai > 185))):
                #print('KOZAI PARTICLE')
                kozai_particles.append(i)
            else:
                nonkozai_particles.append(i)
        else:
            #nonresonant_particles.append(i)
            pass
            
    with open("{}/Particles_in_Kozai_{}.txt".format(kozDir, ST), "w+") as my_file:                               
        for i in kozai_particles:
            my_file.write(str(i)+"\n")
    with open("{}/Particles_not_in_Kozai_{}.txt".format(kozDir, ST), "w+") as my_file:                               
        for i in nonkozai_particles:
            my_file.write(str(i)+"\n")
    
    
    #e and i values for ei plot of kozai particles
    koz_plot_data = []
    nonkoz_plot_data = []    

    ek = []
    ik = []
    for i in e[kozai_particles]:
        ek.append(i[0])
    for j in incl[kozai_particles]:
        ik.append(j[0])
    koz_plot_data.append(ek) 
    koz_plot_data.append(ik)  
    np.savetxt('{}/koz_plot_data.txt'.format(kozDir), koz_plot_data, fmt = '%s')
    
    
    enk = []
    ink = []
    
    for j in e[nonkozai_particles]:
        enk.append(j[0])
    for k in incl[nonkozai_particles]:
        ink.append(k[0])
    nonkoz_plot_data.append(enk) 
    nonkoz_plot_data.append(ink) 
    np.savetxt('{}/nonkoz_plot_data.txt'.format(kozDir), nonkoz_plot_data, fmt = '%s')

    
    #omega values for Kozai libration plots
    ok = []
    for l in arg_peri[kozai_particles]:
        ok.append(l)
    np.savetxt('{}/omega_koz_data.txt'.format(kozDir), ok, fmt = '%s')
    onk = []
    for l in arg_peri[nonkozai_particles]:
        onk.append(l)
    np.savetxt('{}/omega_nonkoz_data.txt'.format(kozDir), onk, fmt = '%s')
    
    
    with open("{}/time_values{}.txt".format(kozDir, ST), "w+") as my_file:                               
        my_file.write(str(t))
        
        
        
        
    print('DONE!')    
    
    
    
    return sim
