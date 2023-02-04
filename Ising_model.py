import matplotlib
matplotlib.use('TKAgg')

import os
import sys
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

random.seed(10)

# calculates random spin and nearest neighbours
def RandSpin(lx, ly):

    #select spin indices randomly
    itrial=np.random.randint(0,lx)
    jtrial=np.random.randint(0,ly)

    # initial and flipped spins
    spin_initial=spin[itrial,jtrial]

    # nearest neigbour spins
    spin_left=spin[np.mod(itrial-1, lx),jtrial]
    spin_right=spin[np.mod(itrial+1, lx),jtrial]
    spin_top=spin[itrial,np.mod(jtrial+1, lx)]
    spin_bottom=spin[itrial,np.mod(jtrial-1, lx)]

    return (spin_initial, spin_left, spin_right, spin_top, spin_bottom), (itrial, jtrial)


# updates system using Glauber Dynamics
def Glauber(lx, ly, spin):

    # selects random spin elements and indices
    spin_positions, spin_idx = RandSpin(lx, ly)
    itrial, jtrial = spin_idx

    # sets spin positions and nearest neaboughr positions 
    spin_initial, spin_left, spin_right, spin_top, spin_bottom = spin_positions
    
    # energy change if spins where flipped
    deltaE = 2* spin_initial * (spin_left + spin_right + spin_bottom + spin_top)

    # flip spin if the probability criteria is met
    Prob_threshold = random.random()


    # spin is always flipped if energy change is negative
    if (deltaE < 0 or np.exp(-(deltaE)/kT) > Prob_threshold): spin[itrial,jtrial]=-spin[itrial,jtrial]
    # return 0 energy if spin not flipped
    else: return 0, spin 

    return deltaE, spin


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
    Prob_threshold = random.random()

    # spin is always flipped if energy change is negative
    if (deltaE_Total < 0 or np.exp(-(deltaE_Total)/kT) > Prob_threshold): 
        spin[itrial1,jtrial1]=spin_initial2
        spin[itrial2,jtrial2]=spin_initial1
    else:
        return 0, spin

    return deltaE_Total, spin


def GetEnergy(N, spin):

    energy = 0

    for i in range(lx):
        for j in range(ly):

            initial_spin = spin[i,j]
            spin_left=spin[np.mod(i-1, lx),j]
            spin_right=spin[np.mod(i+1, lx),j]
            spin_top=spin[i,np.mod(j+1, lx)]
            spin_bottom=spin[i,np.mod(j-1, lx)]

            energy += - initial_spin * (spin_left + spin_right + spin_top + spin_bottom) 

    return energy

# calculates spin config susceptability
def getSuscpt(N, magnetism, kT):

    mean_mag = (1/N) * magnetism
    mean_mag_sq = (1/N) * magnetism**2
    suscept = (1/N*kT) * ( mean_mag_sq - mean_mag**2 )

    return suscept

# calculates heat capacity data
def getHeatCap(N, energy, kT):

    mean_energy = (1/N) * energy
    mean_energy_sq = (1/N) * energy**2
    h_capac = (1/N*kT**2) * ( mean_energy_sq - mean_energy**2 )

    return h_capac


def plot(spin, obs, n, data):

    # records data once 100 sweeps has been exceeded
    if (True):           
        suscept, h_capac, energy, magnetism = obs 
        print(obs)
        data.write('{0:.0f} {1:5.5e} {2:5.5e} {3:5.5e} {4:5.5e}\n'.format(n, suscept, h_capac, energy, magnetism))

    plt.cla()
    im=plt.imshow(spin, animated=True)
    plt.draw()
    plt.pause(0.0001)

# END OF FUNCTIONS

J=1.0
nstep=10000

#input

if(len(sys.argv) != 4):
    print ("Usage python ising.animation.py N T model")
    sys.exit()

lx=int(sys.argv[1]) 
ly=lx 
kT=float(sys.argv[2])
model=str(sys.argv[3]) 

N = lx

# exits code if incorrect algorithm inputs are set
if (model != 'Kawasaki' and  model != 'Glauber'):
    print('Error: choose from model arguments\n1 -- Glauber\n2 -- Kawasaki' )
    exit()

# set 2D matrix for spins
# setting up/down spins randomly
spin = np.random.randint(2, size=(N, N))
spin[spin ==0] = -1

# get initial energy of spin configuration
energy = GetEnergy(N, spin)

# setting up animantion figure
fig = plt.figure()
im=plt.imshow(spin, animated=True)

# opening file to print measurement data
outFilePath = os.getcwd() + f'\\Data\\{model}\\{lx}N_Temp{kT}_{model}Model.dat'
data=open( outFilePath,'w' )

#update loop here - for Glauber dynamic
for n in range(nstep):
    for i in range(lx):
        for j in range(ly):

            # uses different algirithms to calculate deltaE dependent on input
            if (model=='Glauber'): deltaE, spin = Glauber(lx, ly, spin)
            if (model=='Kawasaki'): deltaE, spin = Kawasaki(lx, ly, spin)

            # adding energy to inital energy config to determine changed values
            energy += deltaE
            magnetism = np.sum(spin)

            # calculating heat capacity
            h_capac = getHeatCap(N, energy, kT)

            # calculates susceptibility
            suscept= getSuscpt(N, magnetism, kT)

    # plot animated update of spin configuration and record measurements every 10 sweeps
    if(n%10==0): 
        obs = (suscept, h_capac, energy, magnetism)
        plot(spin, obs, n, data)
    
data.close()