from   numpy.random import Generator, PCG64     # numpy randon number generator routines
import numpy             as np
import matplotlib.pyplot as plt
import random           
from itertools import accumulate

#------------------------ This is the function for Gillespie's  ------------------------#
def gillespie( x0, jps, kap, T, r):
    
    """ This is the funciton of Gillespie' algorithm. 
    It will simulate Xt, the number of compounds in the reactions. 
    Arguments
    x0: vector of size of number of different types of compounds. Intial status.
    jps: set of vectors of same size of x0. Jumps of the reactions. 
    kap: vector of size of number of reactions. Stores the reaction rate. 
    T: (type = float) final time of the reaction 
    return
    x: the list that contains all steps for the reactions
    t: a list of jump time """

    xt = x0
    x = [x0]
    t = [0]

    # simulate the process while within final time T
    while (t[-1] < T):
        # generate an exp. r.v. R 
        rate = reactRate( xt, kap, r)
        #print(rate)
        s = np.sum(rate)
        R = np.random.exponential(1/s)
        t.append(t[-1] + R)


        # generate the next reaction 
        rand = np.random.uniform(0,1)
        reactProb = [x/s for x in rate]
        accumProb = list(accumulate(map(float,reactProb)))
        lStar = compare(accumProb, rand)
        xt = np.add(xt, jps[lStar])
        x.append(xt.tolist())
        
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
    return rate

#------------------------- Helper Function Finding the Next Jump -----------------------#
def compare(lst, num):
    for i in range(len(lst) - 1):
        if num <= lst[i]: return i 
    return 0


################################### Main Program #########################################

# initialize the parameters
#x0 = [1, 10, 50, 10, 0]
x0 = [20, 200, 1000, 200, 0]
jps = [[0,1,0,0,0],[0,0,1,0,0],[0,-1,0,0,0],[0,0,-1,0,0],[0,0,-2,1,0],[0,0,0,-1,0],
        [-1,0,0,-1,1],[1,0,0,1,-1]]
kap = [200, 10, 25, 1, 0.01, 1, 0, 0]
#kap = [200, 10, 25, 1, 0.01, 1, 2, 0.1]

T = 10
r = 8

# run the algorithm 
x, t = gillespie(x0, jps, kap, T, r)

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
title = "G Trajectory of Compounds with k =", kap
ax.set_title(title)                                  
plt.show() 