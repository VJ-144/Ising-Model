import os
import numpy as np


def main():

    model = 'Glauber'
    # model = Kawasaki

    N = 50
    T_range = np.arange(1,3.1, 0.1)
    
    print('Batch Simulation Starting \n')

    for i in range(len(T_range)):
        cmd = f'python Ising_model.py {N} {np.round(T_range[i], 2)} {model}'
        os.system(cmd)

        # prints completed simulation parameters to screen
        print(f'Succesful Simulation @ Parameters: N={N} T={np.round(T_range[i], 2)} Model={model}')

    print('\nBatch Simulations Complete')

    return 0


main()