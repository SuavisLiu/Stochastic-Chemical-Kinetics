# Stochastic Chemical Kinetics

## Intro
This project contains three algorithm modelling the chemical reactions. 

## Next Reaction
    
## Gillespie Algorithm
    This is the funciton of Gillespie' algorithm. 
    It will simulate Xt, the number of compounds in the reactions. 
    Arguments
    x0: vector of size of number of different types of compounds. Intial status.
    jps: set of vectors of same size of x0. Jumps of the reactions. 
    kap: vector of size of number of reactions. Stores the reaction rate. 
    T: (type = float) final time of the reaction 
    r: number of possible reactions
    n: numebr of different kinds of compounds
    return: xT --- the final status up to time T.


## Tau Leaping algorithm 
    This is the funciton of tau leaping algorithm. 
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
    t: a list of jump time
