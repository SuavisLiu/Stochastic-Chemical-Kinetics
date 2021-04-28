from   numpy.random import Generator, PCG64     # numpy randon number generator routines
import numpy             as np
import matplotlib.pyplot as plt
import random           
from array import *
from itertools import accumulate

#------------------------ This is the function for Tau Leaping  ------------------------#
def tauleap( x0, kap, jps, T, h, r):
    """This is the funciton of tau leaping algorithm. 
    It will simulate Xt, the number of compounds in the reactions. 
    Arguments
    x0: vector of size of number of different types of compounds. Intial status.
    jps: set of vectors of same size of x0. Jumps of the reactions. 
    kap: vector of size of number of reactions. Stores the reaction rate. 
    T: (type = float) final time of the reaction 
    r: number of possible reactions
    h: step size
    return
    x: the list that contains all steps for the reactions
    t: a list of jump time """

    # initial status
    xt = x0
    x = [x0]
    t = [0]
    ks = np.array(jps).transpose()
    while (t[-1] < T):
        # generate tau 
        tau = np.array(reactRate( xt, kap, r))
        #tau = (t[-1] - h * (len(t)-1)) * np.array(xt)
        tau[tau < 0] = 0
        #print(tau)
        # geenrate Y
        Y = np.random.poisson(tau * h)
    
        #print(Y)

        if (t[-1] + h > T):
            h = T - t[-1]

        # the change over the current interval
        dx = np.matmul(ks, Y)  # return an array 
        #print(dx)
        xt = np.add(xt, dx)
        x.append(xt.tolist())
        t.append(t[-1] + h)


    return x, t

#------------------------- Helper Function Computing Reaction Rate ---------------------#
def reactRate( xt, kap, r):

    """ This function returns a r-vector for the reaction rate """
    xt = np.array(xt)
    xt[xt < 0] = 0
    xt.tolist()
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

################################### Main Program #########################################

# initialize the parameters
#x0 = [1, 10, 50, 10, 0]
x0 = [20, 200, 1000, 200, 0]
jps = [[0,1,0,0,0],[0,0,1,0,0],[0,-1,0,0,0],[0,0,-1,0,0],[0,0,-2,1,0],[0,0,0,-1,0],
        [-1,0,0,-1,1],[1,0,0,1,-1]]
kap = [200, 10, 25, 1, 0.01, 1, 2, 0.1]
#kap = [200, 10, 25, 1, 0.01, 1, 0, 0]
T = 10
r = 8
h = 0.01

# run the algorithm 
x, t = tauleap(x0, kap, jps, T, h, r)


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
title = "Tau Trajectory of Compounds with k =", kap
ax.set_title(title)                                  
plt.show() 