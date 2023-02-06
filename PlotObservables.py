import os
import numpy as np
import matplotlib.pyplot as plt


def main():

    # model = 'Kawasaki'
    model = 'Glauber'

    pathToFile = os.getcwd() + f'/Data/{model}/'
    directory = os.fsencode(pathToFile)

    tot_energy = []
    tot_mag = []
    tot_hCap = []
    tot_sus = []

    tot_kT = []
        
    for file in os.listdir(directory):

        filename = os.fsdecode(file)
        directory = os.fsdecode(directory)
        path = os.path.join(directory, filename)

        kT = float(filename[8:11])
        N = int(filename[0:2])

        rawData = np.loadtxt(path)
        energy = rawData[:,0]
        mag = rawData[:,1]

        aver_En = np.mean(energy)
        aver_Mag = np.mean(mag)
        
        h_capac = (1/(N*N*kT*kT)) * np.var(energy)
        suscept = (1/(N*N*kT)) * np.var(mag)


        tot_energy.append(aver_En)
        tot_mag.append(aver_Mag)
        tot_hCap.append(h_capac)
        tot_sus.append(suscept)
        tot_kT.append(kT)


    # setting canvas for plotting
    fig, ax = plt.subplots(2, 2, figsize=(7, 5))

    # setting figure title
    fig.suptitle(f'{model} Model at Varied T', fontsize=16)
    fig.subplots_adjust(top=0.8, hspace=0.4)

    # plotting susceptibility + title
    ax[0,0].set_title('Susceptibility', pad=16)
    ax[0,0].plot(tot_kT, tot_sus, marker='x', linestyle='--')
    # ax[0,0].ticklabel_format(axis='x', style='sci', scilimits=(4,4))

    # plotting heat capacity + title
    ax[0,1].set_title('Heat Capacity', pad=16)
    ax[0,1].plot(tot_kT, tot_hCap, marker='x', linestyle='--')
    # ax[0,1].ticklabel_format(axis='both', style='sci', scilimits=(4,4))

    # plotting magnetism + title
    ax[1,0].set_title('Magnetism', pad=4)
    ax[1,0].plot(tot_kT, tot_mag, marker='x', linestyle='--')
    # ax[1,0].ticklabel_format(axis='both', style='sci', scilimits=(4,4))

    # plotting energy + title
    ax[1,1].set_title('Energy', pad=4)
    ax[1,1].plot(tot_kT, tot_energy, marker='x', linestyle='--')
    # ax[1,1].ticklabel_format(axis='both', style='sci', scilimits=(4,4))

    # displays plot
    plt.show()


# runs script
if __name__ == "__main__":
    main()