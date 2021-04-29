from   numpy.random import Generator, PCG64     # numpy randon number generator routines
import numpy             as np
import matplotlib.pyplot as plt
import random           
from itertools import accumulate

import gillespie
import tau_leap


#------------------------- Helper Function Finding the index  -----------------------#
def compare(lst, num):
    for i in range(len(lst) - 1):
        if num <= lst[i]: return i 
    return 0

################################### Main Program #########################################
# initialize the parameters
x0 = [1, 10, 50, 10, 0]
jps = [[0,1,0,0,0],[0,0,1,0,0],[0,-1,0,0,0],[0,0,-1,0,0],[0,0,-2,1,0],[0,0,0,-1,0],
        [-1,0,0,-1,1],[1,0,0,1,-1]]
#kap = [200, 10, 25, 1, 0.01, 1, 2, 0.1]
kap = [200, 10, 25, 1, 0.01, 1, 0, 0]
T = 10
r = 8
error = []
hlist = []
h = 0.0001


for i in range (3):

    h = h/2
    maxE = 0

    x_gillespie, t_gillespie = gillespie.gillespie(x0, jps, kap, T, r)
    x_tau, t_tau = tau_leap.tauleap( x0, kap, jps, T, h, r)

    t = 0
    while (t < T):

        idx1 = compare(t_gillespie, t)
        idx2 = compare(t_tau, t)

        x = np.array(x_gillespie[idx1])
        y = np.array(x_tau[idx2])

        z = sum(np.absolute(x - y))
        if maxE < z:
            maxE  = z
        
        t = t + 0.25

    error = error + [maxE]
    hlist = hlist + [h]

p = np.log(error[0]/error[1]) / np.log(hlist[0]/hlist[1])

print(error)
print(hlist)
print(p)

hlist = np.log(hlist)
error = np.log(error)

fig, ax = plt.subplots()                                
ax.plot(hlist , error , '*-')     
ax.legend()                                             
ax.set_ylabel('log error')                                     
ax.set_xlabel('log h')
title = "Error vs. h in loglog Plot"
ax.set_title(title)                                  
plt.show()