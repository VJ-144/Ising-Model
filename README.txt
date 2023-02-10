INSTRUCTIONS

To run this ising simulation we must require conditions that the correct directories are created
for output data to be stored. The simulation is run for 10100 sweeps and the first 100 are discarded
to allow the system to reach equilibrium.

The simulation conditions are as follows:

1 -- There are directories /Data/Glauber/ and /Data/Kawasaki/ in the same directories as the run.isling.simation.py, 
     PlotObservables.py and IsingModel_Functions.py files.

     This is nessesary as data is stored to the files in these directories every 10 sweeps when sweeps>100.
     The names of these directories are case sensitive.

FURTHER CODEBASE EXPLANATION 

PYTHON FILES 

----    The run.isling.simation.py runs the isling simulation. The following parameters must be specified in terminal
        when running the code in this order.

        'python run.isling.simulation.py N kT Model BatchRun'

        Parameters:
            N - length of axis for square spin matrix configuration
            kT - desired temperature to run simulation
            Model - must specify either 'Glauber' or 'Kawasaki'
            BatchRun - must specify either 'false' or 'true'


            The BatchRun parameter runs the code for varied temperatures from 1.0-3.0K in increments of 0.1K
            while feeding the previously calculated spin matrix back into the simulation

            example: 'python run.isling.simulation.py 50 1.5 Glauber false'
                     This would run the ising simulation for a single temperature 1.5 at matrix size 50x50


----    The PlotObservables.py file takes the stored data files to calculate heat capacity and susceptabilty data   
        as well as plot the graphs at varied temperatures. This file uses previously run and saved data of the ising
        simulation for both Glauber and Kawasaki dynamics at 50x50N to create plots.
        The code is run with the arguments:

        'python PlotObservables.py Model'

        Parameters:
            Model - must specify either 'Glauber' or 'Kawasaki'


----    The IsingModel_Functions.py contains all the functions for the ising simulation.



DATA AND GRAPHS FOR MARKING


----    The raw magnetism and energy data calculated with every simulation is stored in the {currentPath}/Data/Glauber 
        and {currentPath}/Data/Kawasaki directories and is overwritten every simulation. This contains all sampled 
        energy and magnetism data for every 10 sweeps when sweeps>100.

        The file naming format is '{int}N_Temp{float}_{ModelName}Model.dat'


----    The total quantative points for the plotted graphs magnetism, specific heat, energy and susceptability against kT
        can be located in files KawasakiProcessed_DataPoints.csv and GlauberProcessed_DataPoints.csv


----    The Raw_Submission_Results directory contain previously sampled energy and magnetism for a batch simulation at 50x50N
        matrix size and temperature ranging from 1.0->3.0K in increments of 0.1K. Plots are run from this simulated data.

        The file naming format is '{int}N_Temp{value}_{ModelName}Model.dat'


----    The plotted graphs for grading can be located within this directory and are called 'GlauberPlots.png'
        and 'KawasakiPlots.png'

