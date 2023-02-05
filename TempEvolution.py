import Ising_model as im
import numpy as np
from tqdm import tqdm

def main():

    N = 50
    model = 'Glauber'
    # model = 'Kawasaki'

    if (model == 'Glauber'): spin = np.random.choice([1], size=(N, N))
    else: spin = np.random.choice([1, -1], size=(N, N))

    kT_range = np.arange(1,3.1, 0.1)

    for i in range(len(kT_range)):

        T = np.round(kT_range[i], 2)

        new_spinConfig = im.update_SpinConfig(model, T, spin, N)
        spin = new_spinConfig

        print(f'Succesful Simulation @ Parameters: N={N} T={np.round(kT_range[i], 2)} Model={model}')




main()