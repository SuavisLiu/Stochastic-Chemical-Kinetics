# Stochastic Chemical Kinetics

## Introduction 
This project contains three algorithm modelling the chemical reactions.

## Next Reaction

** nextReaction.py contains the following main function. **

Inspired by Bill

This is the funciton of Next Reaction algorithm. 
It will simulate Xt, the number of compounds in the reactions. 

** Arguments **
x0: vector of size of number of different types of compounds. Intial status.
jps: set of vectors of same size of x0. Jumps of the reactions. 
kap: vector of size of number of reactions. Stores the reaction rate. 
T: (type = float) final time of the reaction 

** Returns **
x: the list that contains all steps for the reactions
t: a list of jump time 
    
    
## Gillespie Algorithm

** gillespie.py contains the following main function. **

This is the funciton of Gillespie' algorithm. 
It will simulate Xt, the number of compounds in the reactions. 
    
    
** Arguments **
x0: vector of size of number of different types of compounds. Intial status.
jps: set of vectors of same size of x0. Jumps of the reactions. 
kap: vector of size of number of reactions. Stores the reaction rate. 
T: (type = float) final time of the reaction 
r: number of possible reactions
n: numebr of different kinds of compounds
    
    
** Returns **
x: the list that contains all steps for the reactions
t: a list of jump time 
    



## Tau Leaping algorithm 

** tau_leap.py contains the following main function. **

   This is the funciton of tau leaping algorithm. 
   It will simulate Xt, the number of compounds in the reactions. 
    
    
   * Arguments
    x0: vector of size of number of different types of compounds. Intial status.
    jps: set of vectors of same size of x0. Jumps of the reactions. 
    kap: vector of size of number of reactions. Stores the reaction rate. 
    T: (type = float) final time of the reaction 
    r: number of possible reactions
    h: step size
    
    
   * return
    x: the list that contains all steps for the reactions
    t: a list of jump time

## Order of Accurarcy 

   orderAcc.py contains the major code for computing the weak order of accurarcy to tau leaping.
   
   Here, we compare the approximate solution --- tau leaping, with the true solution --- gillespie. 
   
   
