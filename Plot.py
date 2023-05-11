"""
This file is used to critically analyse the output data from the ising simulation and 
calculate the specific heat and susceptability measurements per kT (temperature point).
The data plotted includes average energy, magnetism, specific heat capacity per spin and susceptability.

This file require the existence of the /Data/Glauber and /Data/Kawasaki directories in the same directory
the simulation is extract the saved data and run the bootstap error analyses. The directory names are 
spelling and case sensitive.
"""

import string
import sys
import os
import numpy as np
import matplotlib.pyplot as plt


def BootstrapError(energy, mag, N, kT):
    """
    Calculates the errors for the specific heat capacity and susceptability using the 
    bootstrap method

    Parameters:
    energy (1D array) - All simulated energy data for a specified temperature kT per 10 sweeps when sweeps>100
    mag (1D array) - All simulated magnetism data for a specified temperature kT per 10 sweeps when sweeps>100
    N (int) - length of one axis of the spin matrix configuration
    kT (float) - initialised temperature of the simulation

    Returns:
    h_capa_err (float) - Error on the heat capacity run at the specified temperature kT
    suscept_err (float) - Error on the susceptability run at the specified temperature kT
    """
    
    # number of groups the data will be split into for sampling
    Nsplit = 10
    data_size = len(energy)

    # list to store sampled heat capacities and suseptabilites
    h_cap_samples = []
    suscept_samples = []

    # looping over number of groups
    for i in range(Nsplit):

        # generating random indices to sample subsets
        sample_idx = np.random.randint(0, data_size, size=(int(data_size/Nsplit)))

        # calculating energy and magnetism subsets
        energy_subset = energy[sample_idx]
        mag_subset = mag[sample_idx]

        # calculating heat capacity and susceptability from subset data
        h_capac_subset = (1/(N*N*kT*kT)) * np.var(energy_subset)
        suscept_subset = (1/(N*N*kT)) * np.var(mag_subset)

        # adding calculated subset data to list
        h_cap_samples.append(h_capac_subset)
        suscept_samples.append(suscept_subset)

    # taking standard deviation of sampled heat capacity and susceptability calculated from subsets
    h_capa_err = np.std(h_cap_samples)
    suscept_err = np.std(suscept_samples)

    # returning errors for heat capacity and susceptability
    return h_capa_err, suscept_err

def main():
    
    # checks command line arguments are correct otherwise stops simulation
    if (len(sys.argv) != 2):
        print("Error Input Usage: python PlotObservables.py model")
        sys.exit()

    # reads terminal argument as Glauber or Kawasaki
    model = str(sys.argv[1])      

    if (model != 'Kawasaki' and  model != 'Glauber'):
        print('Model Parameter Error: choose from model arguments\n1 -- Glauber\n2 -- Kawasaki' )
        sys.exit()

     
    # reads data from directory with raw data of 50x50 spin configuration previously simulated
    pathToFile = os.getcwd() + f'/Raw_Submission_Results/Good_{model}/'
    directory = os.fsencode(pathToFile)

    # data lists for energy, magnetism, heat capacity and temperature
    tot_energy = []
    tot_mag = []
    tot_hCap = []
    tot_sus = []
    tot_kT = []

    # errors list for energy, magnetism, heat capacity and temperature
    tot_hCap_err = []
    tot_suscept_err = []
    tot_energy_err = []
    tot_mag_err = []
    
    # loops over data files in the Raw_Submission_Results directory to read results (1 loop == 1 kT simulation)
    for file in os.listdir(directory):

        # finds the path to data files
        filename = os.fsdecode(file)
        directory = os.fsdecode(directory)
        path = os.path.join(directory, filename)

        # calculates temperature and matrix size from data file name
        kT = float(filename[8:11])
        N = int(filename[0:2])

        # reads in stored energy and magnetism data from file
        rawData = np.loadtxt(path)
        energy = rawData[:,0]
        mag = rawData[:,1]

        # calculates average energy and magnetism
        aver_En = np.mean(energy)
        aver_Mag = np.mean(mag)
        
        # calculates specific heat capacity and susceptability
        h_capac = (1/(N*N*kT*kT)) * np.var(energy)
        suscept = (1/(N*N*kT)) * np.var(mag)

        # calculates error on heat capacity and suscepability using bootstrap method
        hCap_err, suscept_err = BootstrapError(energy, mag, N, kT)

        # stores averaged and calculated data in list
        tot_energy.append(aver_En)
        tot_mag.append(aver_Mag)
        tot_hCap.append(h_capac)
        tot_sus.append(suscept)
        tot_kT.append(kT)

        # stores averaged and calculated errors in list
        tot_hCap_err.append(hCap_err)
        tot_suscept_err.append(suscept_err)
        tot_energy_err.append(np.std(energy))
        tot_mag_err.append(np.std(mag))


    np.savetxt(f'{model}Processed_DataPoints.csv', np.array([tot_kT, tot_energy, tot_mag, tot_hCap, tot_sus]).T, delimiter='\t', fmt="%s", header="kT Energy    Magnetism   Heat_Cap   Susceptability ")

    # setting canvas for plotting
    fig, ax = plt.subplots(2, 2, figsize=(7, 5))

    # setting figure title
    fig.suptitle(f'{model} Model {N}x{N}N at Varied Temperature', fontsize=16)
    fig.subplots_adjust(top=0.8, hspace=0.55, wspace=0.4)

    # plotting susceptibility
    ax[0,0].set_title('Susceptibility', pad=16)
    ax[0,0].errorbar(tot_kT, tot_sus, marker='o', markersize = 4, linestyle='--', yerr=tot_suscept_err, color='black', capsize=3)
    ax[0,0].set_xlabel('kT [K]')
    ax[0,0].set_ylabel('$\chi(M)$ [-]')


    # plotting heat capacity
    ax[0,1].set_title('Heat Capacity', pad=16)
    ax[0,1].errorbar(tot_kT, tot_hCap, marker='o', markersize = 4, linestyle='--', yerr=tot_hCap_err, color='black', capsize=3)
    ax[0,1].set_ylabel('$C/N$ [-]')
    ax[0,1].set_xlabel('kT [K]')


    # plotting magnetism
    ax[1,0].set_title('Magnetism', pad=4)
    ax[1,0].errorbar(tot_kT, tot_mag, marker='o', markersize = 4, linestyle='--', yerr=tot_mag_err, color='black', capsize=3)
    ax[1,0].set_ylabel('$\mid M \mid$ [-]')
    ax[1,0].set_xlabel('kT [K]')


    # plotting energy
    ax[1,1].set_title('Energy', pad=4)
    ax[1,1].errorbar(tot_kT, tot_energy,marker='o', markersize = 4, linestyle='--', yerr=tot_energy_err, color='black', capsize=3)
    ax[1,1].set_ylabel('Average E [-]')
    ax[1,1].set_xlabel('kT [K]')

    # displays and saves plot
    plt.savefig(f'{model}Plots.png')
    plt.show()
    


# runs script
if __name__ == "__main__":
    main()