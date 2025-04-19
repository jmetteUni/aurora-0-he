# Aurora Hydrothermal Venting Application for the Regional Ocean Modelling System
## Description
This application is used for the simulation of the Aurora hydrothermal vent system in the Arctic Ocean. It is intended to be used with a modified version of the Regional Ocean Modelling System (ROMS) and OpenMPI.

The application uses a 256x256 horizontal grid with 30 vertical layers. Initial and boundary conditions are inferred from the ERA5 Interim reanalysis product. Surface fluxes are set to zero. Input of heat flux and a passive tracer, to simulate the hydtrothermal vent is prescirbed with forcing files. Tides are simulated based on tidal constituents from the Arctic Ocean Tidal Inverse Model. 

The simulation is run with a timestep of 4 sec for one month. 
## How to use
- Install ROMS from this fork [https://github.com/jmetteUni/roms](https://github.com/jmetteUni/roms) or modify the [official version](https://github.com/myroms/roms) accordingly. Documentation on how to install ROMS can be found in the associated README files.
- Install OpenMPI.
- Clone this directory outside of the directory where ROMS is installed.  
- Build the application by executing the build script with
'''
./build_roms.sh
'''
Optionally use
'''
./build_roms.sh -j N
'''
where N is the number of processes used for buidling.
- Run the application with 
'''
mpirun -np N romsG aurora-1.in 
'''
where N is again the number of processe, "romsG" is the ROMS executable and "aurora-1.in" is the input script.
 
