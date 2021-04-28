from   numpy.random import Generator, PCG64     # numpy randon number generator routines
import numpy             as np
import matplotlib.pyplot as plt
import random           
from itertools import accumulate

#------------------------ This is the function for Next Reaction  -----------------------#

def nextRact( x0, jps, kap, T, r):
    """
    This is the funciton of Next Reaction algorithm. 
    It will simulate Xt, the number of compounds in the reactions. 
    Arguments
    x0: vector of size of number of different types of compounds. Intial status.
    jps: set of vectors of same size of x0. Jumps of the reactions. 
    kap: vector of size of number of reactions. Stores the reaction rate. 
    T: (type = float) final time of the reaction 
    return
    x: the list that contains all steps for the reactions
    t: a list of jump time 
    """

    t = [0]
    x = [x0]
    xt = x0
    tau = np.array([0, 0, 0, 0, 0, 0, 0 ,0])


    while (t[-1] < 10):
        print(t[-1])
        G = invG( t[-1], T, tau, xt, kap, r)
        m = G.argmin()
        xt = np.add(xt, jps[m])
        x.append(xt.tolist())
        T[m] = T[m] - np.log(np.random.uniform(0,1))
        tau = tau + ( G[m] - t[-1] ) * reactRate( xt, kap, r)
        dt = G[m]
        t.append(dt)

    return x, t
#------------------------- Helper Function Computing Reaction Rate ---------------------#

def reactRate( xt, kap, r):
    rate = np.zeros(r)
    rate[0] = kap[0] * xt[0]
    rate[1] = kap[1] * xt[1]
    rate[2] = kap[2] * xt[1]
    rate[3] = kap[3] * xt[2]
    rate[4] = kap[4] * xt[2] * (xt[2] - 1)
    rate[5] = kap[5] * xt[3]
    rate[6] = kap[6] * xt[0] * xt[3]
    rate[7] = kap[7] * xt[4]
    return np.array(rate)


def invG( t, T, tau, xt, kap, r):
     return(t + np.divide((T - tau), reactRate(xt, kap, r)))

################################### Main Program #########################################

#x0 = [1, 10, 50, 10, 0]
x0 = [20, 200, 1000, 200, 0]

#kap = [200, 10, 25, 1, 0.01, 1, 0, 0]
kap = [200, 10, 25, 1, 0.01, 1, 2, 0.1]

jps = [[0,1,0,0,0],[0,0,1,0,0],[0,-1,0,0,0],[0,0,-1,0,0],[0,0,-2,1,0],[0,0,0,-1,0],
        [-1,0,0,-1,1],[1,0,0,1,-1]]
r = 8
T = - np.log(np.random.uniform(0,1,r))

x, t = nextRact( x0, jps, kap, T, r)

x = np.array(x)
numG = x[:,0]
numM = x[:,1]
numP = x[:,2]
numD = x[:,3]
numB = x[:,4]
xt = x[-1]

################################# Meaningless Stuff #######################################

print('Final State:', xt)

fig, ax = plt.subplots()                                
#ax.plot(time, numG, '-', label = 'numG')    
ax.plot(t, numM, '-', label = 'numM') 
ax.plot(t, numP, '-', label = 'numP')   
ax.plot(t, numD, '-', label = 'numD')
#ax.plot(time, numB, '-', label = 'numB')            
ax.legend()                                             
ax.set_ylabel('num of compounds')                                     
ax.set_xlabel('Time')
title = " Trajectory of Compounds with k =", kap
ax.set_title(title)                                  
plt.show() 
