import os
import numpy as np


def main():

    N = 50
    model = 'Glauber'
    # model = Kawasaki
    T_range = np.arange(1,3.1, 0.1)
    
    print('Batch Simulation Starting \n')

    for i in range(len(T_range)):
        cmd = f'python checkpoint1.py {N} {np.round(T_range[i], 2)} {model}'
        os.system(cmd)

        # print to completed simulations to screen
        print(f'Successful Simulation: N={N} T={np.round(T_range[i], 2)} Model={model}')

    print('\nBatch Simulations Complete')

    return 0


main()