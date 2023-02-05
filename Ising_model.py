import itertools

import os
import sys
import math
import time
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from tqdm import tqdm

import cProfile

start_time = time.time()
random.seed(444)

# calculates random spin and nearest neighbours
def spin_NN(idx, spin, N):

    # set spin indices and value from previously generated random numbers
    itrial, jtrial = idx
    spin_initial=spin[itrial,jtrial]

    # nearest neigbour spins
    spin_left=spin[np.mod(itrial-1, N),jtrial]
    spin_right=spin[np.mod(itrial+1, N),jtrial]
    spin_top=spin[itrial,np.mod(jtrial+1, N)]
    spin_bottom=spin[itrial,np.mod(jtrial-1, N)]

    return spin_initial, spin_left, spin_right, spin_top, spin_bottom


# generates random number for all spins to randomly sample
def GenerateRandom_Idx(N):
    # position of spin to be switched
    i = np.random.choice(N, N*N)
    j = np.random.choice(N, N*N)
    return (i, j)


# updates system using Glauber Dynamics
def Glauber(idx, spin, kT, N):

    # sets spin positions and nearest neaboughr positions and indices
    spin_initial, spin_left, spin_right, spin_top, spin_bottom = spin_NN(idx, spin, N)
    itrial, jtrial = idx
    
    # energy change if spins where flipped
    deltaE = 2 * J * spin_initial * (spin_left + spin_right + spin_bottom + spin_top)

    # flip spin if the probability criteria is met
    Prob_threshold = np.random.random()

    # spin is always flipped if energy change is negative
    if (deltaE < 0 or np.exp(-(deltaE)/kT) > Prob_threshold):
        spin[itrial,jtrial] *= -1
        return deltaE, spin
    # return 0 energy if spin not flipped
    else:  
        return 0, spin 

    


# updates system using Kawasaki Dynamics
def Kawasaki(lx, ly, spin):
    
    # selects 2 random spin locations
    spin_positions1, spin_idx1 = RandSpin(lx, ly)
    itrial1, jtrial1 = spin_idx1
    spin_initial1, spin_left1, spin_right1, spin_top1, spin_bottom1 = spin_positions1

    spin_positions2, spin_idx2 = RandSpin(lx, ly)
    itrial2, jtrial2 = spin_idx2
    spin_initial2, spin_left2, spin_right2, spin_top2, spin_bottom2 = spin_positions2


    # ensures we do not randomly sample the exact same spin 
    while ( (spin_initial1) == (spin_initial2) ): 
        spin_positions2, spin_idx2 = RandSpin(lx, ly)
        itrial2, jtrial2 = spin_idx2
        spin_initial2, spin_left2, spin_right2, spin_top2, spin_bottom2 = spin_positions2

    # calculating flipped energies of seperate spins
    deltaE1 = 2 * spin_initial1 * (spin_left1 + spin_right1 + spin_bottom1 + spin_top1)
    deltaE2 = 2 * spin_initial2 * (spin_left2 + spin_right2 + spin_bottom2 + spin_top2)

    # calculating total energy
    deltaE_Total = deltaE1 + deltaE2

    # correction for nearest neubhours
    if ( [itrial1,jtrial1] == [np.mod(itrial2-1, lx), jtrial2] ): deltaE_Total += 4
    if ( [itrial1,jtrial1] == [np.mod(itrial2+1, lx), jtrial2] ): deltaE_Total += 4
    if ( [itrial1,jtrial1] == [itrial2, np.mod(jtrial2+1, lx)] ): deltaE_Total += 4
    if ( [itrial1,jtrial1] == [itrial2, np.mod(jtrial2-1, lx)] ): deltaE_Total += 4

    # flip spin if the probability criteria is met
    Prob_threshold = np.random.random()

    # spin is always flipped if energy change is negative
    if (deltaE_Total < 0 or np.exp(-(deltaE_Total)/kT) > Prob_threshold): 
        spin[itrial1,jtrial1]=spin_initial2
        spin[itrial2,jtrial2]=spin_initial1
    else:
        return 0, spin

    return deltaE_Total, spin


def GetEnergy(spin, N):

    energy = 0
    for ij in itertools.product(range(N), repeat=2):

        #indices for look
        i,j = ij

        # finding nearest neibours
        initial_spin = spin[i,j]
        spin_left=spin[np.mod(i-1, N),j]
        spin_right=spin[np.mod(i+1, N),j]
        spin_top=spin[i,np.mod(j+1, N)]
        spin_bottom=spin[i,np.mod(j-1, N)]

        energy += - J * initial_spin * (spin_left + spin_right + spin_top + spin_bottom) 
    return energy/2


def initialise_simulation():

    # reads inputs from terminal
    N=int(sys.argv[1]) 
    kT=float(sys.argv[2])
    model=str(sys.argv[3])

    # checks command line arguments are correct
    if (len(sys.argv) != 4):
        print("Usage python ising.animation.py model N T_inital")
        sys.exit()
    elif (model != 'Kawasaki' and  model != 'Glauber'):
        print('Error: choose from model arguments\n1 -- Glauber\n2 -- Kawasaki' )
        sys.exit()    

    return model, kT, N



def update_SpinConfig(model, kT, spin, N):

    # assign non changing variables as global
    global J
    J=1

    # setting up animantion figure
    # fig = plt.figure()
    # im=plt.imshow(spin, animated=True)

    # get initial energy of spin configuration
    energy = GetEnergy(spin, N)

    # number of sweeps for simulation
    nstep=10100

    # list for average magnetism and energy data every 10 sweeps when sweeps > 100
    total_mag = []
    total_energy = []

    outFilePath = os.getcwd() + f'/Data/{model}/{N}N_Temp{kT}_{model}Model.dat'
    data=open( outFilePath,'w')
    # data.write("{0} {1}\n".format('Energy', 'Magnetism'))

    # sweeps counter
    sweeps = 0

    for n in range(nstep ):
        
        # generates new random indices to sample the spin config with every sweep
        rand_x1, rand_y1 = GenerateRandom_Idx(N)
        if (model=='Kawasaki'): rand_x2, rand_y2 = GenerateRandom_Idx(N)

        # loops over entire spin config number of elements
        for i in range(N*N):

            # uses different algirithms to calculate deltaE dependent on input
            if (model=='Glauber'):

                # select indice to sample spin location
                idx = ( rand_x1[i], rand_y1[i] )

                # calculate energy change and new spin config using glauber
                deltaE, spin = Glauber(idx, spin, kT, N)

            # uses Kawasaki model to calculate deltaE if specified
            elif (model=='Kawasaki'):

                # select 2 indices to sample spin location
                idx1 = ( rand_x1[i], rand_y1[i] )
                idx2 = ( rand_x2[i], rand_y2[i] )

                # calculate energy change and new spin config using kawasaki
                deltaE, spin = Kawasaki(idx1, idx2, spin)

            # change total energy
            energy += deltaE

        # plot animated update of spin configuration and record measurements every 10 sweeps
        if(n%10==0 and n>100): 

            # prints current sweep to terminal
            sweeps +=10
            print(f'sweeps={sweeps}', end='\r')

            # calculating new spin config matrix energy and magnetism
            magnetism = np.sum(spin)
            
            # recording calculated magnetism and energy to list
            total_mag.append(magnetism)
            total_energy.append(energy)


            # taking the average of all recorded data
            # mean_energy = np.mean(total_energy)
            # mean_mag = np.mean(total_mag)



            # writing data to file
            data.write('{0:.0e} {1:5.5e}\n'.format(energy, magnetism))

            # closing data file


            # animates spin configuration 
            # plt.cla()
            # im=plt.imshow(spin, animated=True)
            # plt.draw()
            # plt.pause(0.0001)    

    data.close()
    # converting magnetism and energy list to numpy arrays
    # total_mag = np.asarray(total_mag)
    # total_energy = np.asarray(total_energy)

    # calculating heat capacity and energy per spin
    # h_capac = (1/(N*N*kT*kT)) * np.var(total_energy)
    # suscept = (1/(N*N*kT)) * np.var(total_mag)

    # # taking the average of all recorded data
    # mean_energy = np.mean(total_energy)
    # mean_mag = np.mean(total_mag)


    # opening file to print all measurement data


    return spin



def main():

    start_time = time.time()

    model, kT, N = initialise_simulation()


    # set 2D matrix for spins
    # setting up/down spins depending on model
    if (model == 'Glauber'): spin = np.random.choice([1], size=(N, N))
    else: spin = np.random.choice([1, -1], size=(N, N))

    new_spin = update_SpinConfig(model, kT, spin, N)

    print("--- %s seconds ---" % (time.time() - start_time))

# main()
#random comment