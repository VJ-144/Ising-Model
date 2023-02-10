"""
This file contains all the functions required to run and visualise an implmentation of the Ising Model
using both Glauber and Kawasaki dynamics. These functions are used to run the simulation in 
the run.ising.simulation.py file.

The functions require the existence of the /Data/Glauber and /Data/Kawasaki directories
to store calculated data and run error analysis. The directory names are spelling and case sensitive.

"""


import os
import sys
import math
import time
import random
import itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def spin_NN(idx, spin, N):
    """
    Calculates the spin value of a single point in the spin matrix configuration given a single index coordinate.
    Also returns nearest neighbor spin values of the indexed coordinate.

    Parameters:
    idx (tuple) - integer indices for spin value of interest (i, j)
    spin (2D array) - spin configuration matrix
    N (int) - length of one axis of the spin matrix configuration

    Returns:
    spin_initial (int) - spin value of given index coordinates
    spin_left (int) - spin vaue of the nearest neighbor to the left of the index coordinate
    spin_right (int) - spin vaue of the nearest neighbor to the right of the index coordinate
    spin_top (int) - spin vaue of the nearest neighbor above the index coordinate
    spin_bottom (int) - spin vaue of the nearest neighbor below the index coordinate
    """

    # setting spin matrix i and j indices for spin value of interest. Extracting spin value from spin configuation.
    itrial, jtrial = idx
    spin_initial = spin[itrial,jtrial]

    # calculating nearest neighbor spin values
    spin_left=spin[np.mod(itrial-1, N),jtrial]
    spin_right=spin[np.mod(itrial+1, N),jtrial]
    spin_top=spin[itrial,np.mod(jtrial+1, N)]
    spin_bottom=spin[itrial,np.mod(jtrial-1, N)]

    # returning spin value for indexed point and nearest neighbor spin values
    return spin_initial, spin_left, spin_right, spin_top, spin_bottom


def GenerateRandom_Idx(N):
    """
    Generates random number which act as the indices to sample the spin configuration for a complete sweep
    of the Glauber or Kawasaki method

    Parameters:
    N (int) - length of one axis of the spin matrix configuration

    Returns:
    i (1D array) - array of random indices to be used in sampling 1 axes (x-axes) of the spin configuration
    j (1D array) - array of random indices to be used in sampling 1 axes (y-axes) of the spin configuration
    """

    # array of random i and j indices to be used for sampling the spin matrix per sweep 
    i = np.random.choice(N, N*N)
    j = np.random.choice(N, N*N)

    # returns tuple of indice i and j arrays to sample the spin configuration 
    return (i, j)


def Glauber(idx, spin, kT, N):
    """
    Spin configuration update for simulating Glauber Dynamics. Returns the change in energy given
    a sampled spin and the updated spin configuration.

    Parameters:
    idx (tuple) - integer indices for spin value of interest (i, j)
    spin (2D array) - spin configuration matrix
    kT (float) - initialised temperature of the simulation
    N (int) - length of one axis of the spin matrix configuration

    Returns:
    deltaE (float) - change in energy for a sampled spin given Glauber Dynamics
    spin (2D array) - updated spin configuration matrix
    """

    # Extracting spin value of the sampled index and it's nearest neighbor from spin configuation.
    itrial, jtrial = idx
    spin_initial, spin_left, spin_right, spin_top, spin_bottom = spin_NN(idx, spin, N)
    
    # calculating the change in energy given the sampled spin value is flipped
    deltaE = 2 * J * spin_initial * (spin_left + spin_right + spin_bottom + spin_top)

    # setting probability criteria for spin flipping
    Prob_threshold = np.random.random()

    # flipping sampled spin value if energy is always less than zero or satisfies the proability criteria
    if (deltaE < 0 or np.exp(-(deltaE)/kT) > Prob_threshold):
        spin[itrial,jtrial] *= -1

        # returns change in energy due to flipped spin and new spin matrix configuration
        return deltaE, spin
    else:  
        # return no change in energy and old spin matrix configuration if update conditions unsatisfied
        return 0, spin 

    
def Kawasaki(idx1, idx2, spin, kT, N):
    """
    Spin configuration update for simulating Kawasaki Dynamics. Returns the change in energy given
    a sampled spin and the updated spin configuration.

    Parameters:
    idx1 (tuple) - integer indices for 1st spin value of interest to be swapper (i1, j1)
    idx2 (tuple) - integer indices for 2nd spin value of interest to be swapper (i2, j2)
    spin (2D array) - spin configuration matrix
    kT (float) - initialised temperature of the simulation
    N (int) - length of one axis of the spin matrix configuration

    Returns:
    deltaE (float) - change in energy for a sampled set of spins given Kawasaki Dynamics
    spin (2D array) - updated spin configuration matrix
    """
    
    # Extracting spin value of the 1st sampled index and it's nearest neighbor from spin configuation.
    spin_initial1, spin_left1, spin_right1, spin_top1, spin_bottom1 = spin_NN(idx1, spin, N)
    itrial1, jtrial1 = idx1

    # Extracting spin value of the 2nd sampled index and it's nearest neighbor from spin configuation.
    spin_initial2, spin_left2, spin_right2, spin_top2, spin_bottom2 = spin_NN(idx2, spin, N)
    itrial2, jtrial2 = idx2

    # calculating total energy change when switching the positions of the sampled spins
    deltaE1 = 2 * J * spin_initial1 * (spin_left1 + spin_right1 + spin_bottom1 + spin_top1)
    deltaE2 = 2 * J * spin_initial2 * (spin_left2 + spin_right2 + spin_bottom2 + spin_top2)
    deltaE_Total = deltaE1 + deltaE2

    # total energy correction for nearest neighbor overlap
    if ( [itrial1,jtrial1] == [np.mod(itrial2-1, N), jtrial2] ): deltaE_Total += 4 * J
    if ( [itrial1,jtrial1] == [np.mod(itrial2+1, N), jtrial2] ): deltaE_Total += 4 * J
    if ( [itrial1,jtrial1] == [itrial2, np.mod(jtrial2+1, N)] ): deltaE_Total += 4 * J
    if ( [itrial1,jtrial1] == [itrial2, np.mod(jtrial2-1, N)] ): deltaE_Total += 4 * J

    # setting probability criteria for spin flipping
    Prob_threshold = np.random.random()

    # switching sampled spin positions if energy is always less than zero or satisfies the proability criteria
    if (deltaE_Total < 0 or np.exp(-(deltaE_Total)/kT) > Prob_threshold): 
        spin[itrial1,jtrial1] *= -1
        spin[itrial2,jtrial2] *= -1

        # returns change in energy due to flipped spin and new spin matrix configuration
        return deltaE_Total, spin
    else:
        # return no change in energy and old spin matrix configuration if update conditions unsatisfied
        return 0, spin


def GetEnergy(spin, N):
    """
    Calculates the total energy of a given spin configuration.

    Parameters:
    spin (2D array) - spin configuration matrix
    N (int) - length of one axis of the spin matrix configuration

    Returns:
    energy - total energy of the spin configuration
    """

    # energy counter to track total contributions of spin energies
    energy = 0

    # looping over 2D spin configuration
    for ij in itertools.product(range(N), repeat=2):

        # setting spin index to sample
        idx = ij

        # Extracting spin value of the sampled index and it's nearest neighbor from spin configuation.
        initial_spin, spin_left, spin_right, spin_top, spin_bottom = spin_NN(idx, spin, N)

        # calculating spin contribution from a single spin value
        energy += - J * initial_spin * (spin_left + spin_right + spin_top + spin_bottom) 

    # returning total energy of the spin configuration and accounting for overcounting
    return energy/2


def initialise_simulation():
    """
    Initialise command line arguments given the file is run from command line. 
    Also checks if the command line arguments are valid.

    Parameters:
    None

    Returns:
    N (int) - length of one axis of the spin matrix configuration
    kT (float) - initialised temperature of the simulation
    model (string) - Dynamics model to be used when updating spin matrix
    BatchRun (bool True/False) - Runs a batch simulation at varied temperatures (1.0-3.0)K in increments of 0.1K if True
    """

    # checks command line arguments are correct otherwise stops simulation
    if (len(sys.argv) != 5):
        print("Error Input Usage: python ising.animation.py N kT model BatchRun")
        sys.exit()

    # reads input arguments from terminal
    N=int(sys.argv[1]) 
    kT=float(sys.argv[2])
    model=str(sys.argv[3])
    BatchRun=str(sys.argv[4])

    if (model != 'Kawasaki' and  model != 'Glauber'):
        print('Model Parameter Error: choose from model arguments\n1 -- Glauber\n2 -- Kawasaki' )
        sys.exit()
    elif (BatchRun != 'False' and  BatchRun != 'True'):
        print('BatchRun Parameter Error: choose to simulate a range/batch of temperatures:\n Parameter Options:\n1 -- True\n2 -- False\nOptions are case sensitive' )
        sys.exit()        

    # returns Dynamics model, simulation temperature and 0-axes length of spin configuration
    return model, kT, N, BatchRun



def update_SpinConfig(model, kT, spin, N):
    """
    Updates spin configuration based upon input dynamics model and simulation temperature. 
    Saves energy and magnetism data to directory /{currentPath}/Data/{model} every 10 sweeps when sweeps>100

    Parameters:
    model (string) - Dynamics model to be used when updating spin (matrix) configuration
    kT (float) - initialised temperature of the simulation
    spin (2D array) - spin configuration matrix
    N (int) - length of one axis of the spin matrix configuration

    Returns:
    spin (2D array) - updated spin configuration matrix after 10100 sweeps
    """

    # assign non-changing J variable as global
    global J
    J=1

    # setting up animantion figure
    fig = plt.figure()
    im=plt.imshow(spin, animated=True)

    # number of sweeps for simulation
    nstep=500

    # path to directory for energy and magnetism data to be stored every 10 sweeps when sweeps > 100
    outFilePath = os.getcwd() + f'/Data/{model}/{N}N_Temp{kT}_{model}Model.dat'
    data=open( outFilePath,'w')

    # calculates energy of the inital spin matrix configuration
    energy = GetEnergy(spin, N)

    # sweeps counter
    sweeps = 0

    # starts loop to apply dynamics model and sample spins
    for n in range(nstep):
        
        # generates new random indices to sample the location of 1 spin value with every sweep.
        # generates additional random indices to sample the location of a second spin for Kawasaki dynamics
        rand_x1, rand_y1 = GenerateRandom_Idx(N)
        if (model=='Kawasaki'): rand_x2, rand_y2 = GenerateRandom_Idx(N)

        # loops over entire spin configuration number of elements
        for i in range(N*N):

            # applies gluber dynamics update rules if model specified
            if (model=='Glauber'):

                # select indice values to sample spin configuration from randomly generated numbers
                idx = ( rand_x1[i], rand_y1[i] )

                # calculate energy change and new spin configuration using glauber dynamics
                deltaE, spin = Glauber(idx, spin, kT, N)

            # applies kawasaki dynamics update rules if model specified
            elif (model=='Kawasaki'):

                # select 2 indices to sample spin location and skips interation if they are the same (as no energy change)

                # selects 2 indice values to sample spin configuration from randomly generated numbers
                idx1 = ( rand_x1[i], rand_y1[i] )
                idx2 = ( rand_x2[i], rand_y2[i] )

                # skips iteration if both sampled spin values are the same, i.e. deltaE=0
                if ( spin[rand_x1[i], rand_y1[i]] == spin[rand_x2[i], rand_y2[i]] ): continue 

                # calculate energy change and new spin configuration using kawasaki dynamics
                deltaE, spin = Kawasaki(idx1, idx2, spin, kT, N)

            # calculates total energy based off incremental changes in energy
            energy += deltaE
            
        # plots animated update of spin configuration and record measurements every 10 sweeps and then sweeps>100
        if(n%10==0 and n>100): 

            # prints current number of sweep to terminal
            sweeps +=10
            print(f'sweeps={sweeps}', end='\r')

            # calculates new spin configuration magnetism
            magnetism = np.abs(np.sum(spin))

            # writes new spin configuration energy and magnetism data to file 
            data.write('{0:5.5e} {1:5.5e}\n'.format(energy, magnetism))

            # animates spin configuration 
            plt.cla()
            im=plt.imshow(spin, animated=True)
            plt.draw()
            plt.pause(0.0001)    

    # closes data file
    data.close()

    # returns updated spin matrix after 10100 sweeps
    return spin



# def main():
#     """
#     Runs the isling simulation for 10100 sweeps for a given temperature
#     and model from the command line and visualises
#     """

#     # records the start of simulation time
#     start_time = time.time()

#     # reads arguments from the command line, i.e. model, temperature and spin matrix size
#     model, kT, N, BatchRun = initialise_simulation()

#     # setting 2D spin matrix configuration all up for glauber model and spin up/down for kawasaki model
#     if (model == 'Glauber'): spin = np.random.choice([1], size=(N, N))
#     else: spin = np.random.choice([1, -1], size=(N, N))

#     # updates given spin configuration for 10100 sweeps
#     new_spin = update_SpinConfig(model, kT, spin, N)

#     # prints simulation time to terminal
#     print("--- %s seconds ---" % (time.time() - start_time))


# main()