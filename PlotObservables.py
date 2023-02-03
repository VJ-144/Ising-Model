import os
import numpy as np
import matplotlib.pyplot as plt


def main():

    # file with Susceptibility and Heat Capacity Data
    file = '50N_Temp1.0_KawasakiModel.dat'

    # selects the correct directiory the data is in based off file name
    if (len(file) > 28):
        model = 'Kawasaki'
    else:
        model = 'Glauber'    

    # finds data from correct directories
    # pathToFile = f'C:\\Users\\Vijay\\OneDrive\\Documents\\Univeristy Work\\Year 5\\MVP\\checkpoint1\\Data\\{model}\\'
    pathToFile = os.getcwd() + f'\\Data\\{model}\\'
    Data = np.loadtxt(pathToFile + file)

    # Susceptibility and Heat Capacity Data
    suscept = Data[:,0]
    heat_cap = Data[:,1]
    energy = Data[:,2]
    magnetism = Data[:,3]

    # setting canvas for plotting
    fig, ax = plt.subplots(2, 2, figsize=(7, 5))

    # setting figure title
    fig.suptitle(f'{file}', fontsize=16)
    fig.subplots_adjust(top=0.8, hspace=0.4)

    # plotting susceptibility + title
    ax[0,0].set_title('Susceptibility', pad=16)
    ax[0,0].plot(suscept)
    ax[0,0].ticklabel_format(axis='x', style='sci', scilimits=(4,4))

    # plotting heat capacity + title
    ax[0,1].set_title('Heat Capacity', pad=16)
    ax[0,1].plot(heat_cap)
    ax[0,1].ticklabel_format(axis='both', style='sci', scilimits=(4,4))

    # plotting magnetism + title
    ax[1,0].set_title('Magnetism', pad=4)
    ax[1,0].plot(magnetism)
    ax[1,0].ticklabel_format(axis='both', style='sci', scilimits=(4,4))

    # plotting energy + title
    ax[1,1].set_title('Energy', pad=4)
    ax[1,1].plot(energy)
    ax[1,1].ticklabel_format(axis='both', style='sci', scilimits=(4,4))

    # displays plot
    plt.show()


# runs script
if __name__ == "__main__":
    main()