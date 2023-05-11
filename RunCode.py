"""
This file is used to run the ising model simulation for a single temperature and
for batches of temperatures which reuse spin configurations from previous calculations.
This file require the existence of the /Data/Glauber and /Data/Kawasaki directories in the same directory
to store calculated data and run error analysis. The directory names are spelling and case sensitive.

This file must be run through terminal with the input arguments as follows:

'python run.ising.simulation.py N kT, model, BatchRun'

A full explanation of the parameters can be found in the README.txt file
"""

# importing ising model functions required
import Functions as func
import numpy as np
import time
import sys


def main():
    """
    Runs the Ising Model simulation for a range of temperatures and continuously 
    uses/feeds back the previously generated spin matrix for new temperatures   
    """

    # spin matrix configuration size (NxN) and model type to simulate dynamics
    # N = 50
    # model = 'Glauber'
    # # model = 'Kawasaki'


    # reads arguments from the command line, i.e. model, temperature and spin matrix size
    model, kT, N, BatchRun = im.initialise_simulation()

    # setting 2D spin matrix configuration all up for glauber model and spin up/down for kawasaki model
    if (model == 'Glauber'): spin = np.random.choice([1], size=(N, N))
    else: spin = np.random.choice([1, -1], size=(N, N))

    # runs the simulation for a range of temperature values from 1->3 in increments of 0.1
    if (BatchRun=='True'):

        # range of temperatures to be simulated
        kT_range = np.arange(1, 3.1, 0.1)        

        # starting simulation loop over temperatures
        for i in range(len(kT_range)):

            # updates spin matrix at temperature T for 10100 sweeps
            T = np.round(kT_range[i], 2)
            new_spinConfig = im.update_SpinConfig(model, T, spin, N)

            # feeds newly evolved spin matrix configuration back into simulation with new temperature
            spin = new_spinConfig

            # prints notification to the terminal upon a succesful 10100 sweeps
            print(f'Succesful Simulation @ Parameters: N={N} T={np.round(kT_range[i], 2)} Model={model}')

    # runs the simulation for a given temperature value
    elif (BatchRun=='False'):

        # records the start of simulation time
        start_time = time.time()

        # updates given spin configuration for 10100 sweeps
        new_spin = im.update_SpinConfig(model, kT, spin, N)

        # prints simulation time to terminal
        print("--- %s seconds ---" % (time.time() - start_time))


main()