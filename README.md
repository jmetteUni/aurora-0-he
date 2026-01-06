# Aurora Hydrothermal Venting Application<br> for the Regional Ocean Modelling System (ROMS)

## Description
This application is used for the simulation of the Aurora hydrothermal vent system in the Arctic Ocean. It is intended to be used with a modified version of the Regional Ocean Modelling System (ROMS) (documented [here](https://github.com/jmetteUni/roms/tree/bottom-tracer) and [OpenMPI](https://www.open-mpi.org/).

The application uses a 256x256 horizontal grid with 30 vertical layers. Initial and boundary conditions are inferred from the ERA5 Interim reanalysis product. Surface fluxes are set to zero. Input of heat flux and a passive tracer, to simulate the hydrothermal vent is prescribed with forcing files. Tides are simulated based on tidal constituents from the Arctic Ocean Tidal Inverse Model.

The simulation is run with a timestep of 4 sec for one month.

For general information on the ROMS model check the official [wiki](https://www.myroms.org/wiki/Documentation_Portal) and the [forum](https://www.myroms.org/forum/), but note that the wiki is outdated in some places. Most of the files have also some documentation as comments inside.

# Preparing the input files


# How to run the model
This is  short tutorial on how to run this application with the ROMS model. For more in-depth documentation see the official ROMS documents. The application in it's structure follows roughly the test case [UPWELLING](https://www.myroms.org/wiki/UPWELLING_CASE). The test cases can be obtained [here](https://github.com/myroms/roms_test). It is recommended that you test your model with the UPWELLING test case first.

## Prerequisites ([wiki](https://www.myroms.org/wiki/Getting_Started))
Software you will need to have installed:<br>
- git<br>
- (an implementation of MPI for running on multiple cores, recommended)<br>
- fortran90 or fortran95 compiler<br>
- version 3.81 or higher<br>
- cpp for C-processing<br>
- Perl<br>
- NetCDF

## Setting up the model and directory structure
There are many different ways how to structure your source code and your application. This way here is not neccesarily the best, but the way it worked for me.

1. Clone the modified ROMS [source code](https://github.com/jmetteUni/roms/tree/bottom-tracer) in one directory.
2. Clone this application into a different directory with a characteristic application name.

Here, the different subdirectories and files which make up the application are briefly explained:<br>
### _Functionals_ ([wiki](https://www.myroms.org/wiki/Functionals))
contains header (.h) files. These are copied from the source code under `ROMS/Functionals/` and are therefore mostly similar. But some are modified specifically for this setup such as _ana_srflux.h_ for example, which specifies the analytical (not prescribed from forcing files) surface forcing. In this file surface forcings are set to 0 for this application to represent a simplified ice cover. The directory also contains header files which are not used by the application (because there function is replaced by forcing files or they just not neccesary) but for simplicity just all files from the source code are copied and they don't hurt.

### _matlab_
This contains matlab scripts which are used to create the input forcing files. The forcing files are NetCDF files which prescribe certain conditions for the model instead of conditions from analytical functions (the header files). They are used for example to simulate realistic boundaries (That is why the UPWELLING test case used no NetCDF files) and therefore can get quite large in terms of disk space. This application uses four different scripts and one _.mat_-file:<br>
- _grid-M256L256.mat_ is a _MATLAB_ file representing the grid of the model domain<br>
- _get_bottom_forcing_jm.m_ to create the forcing file for the input of the passive bottom tracer<br>
- _get_initial_condition_from_GLORYS_jm.m_ to interpolate ERA5 data to the ROMS grid for the initial conditions<br>
- _get_boundary_condition_GLORYS_jm.m_ to interpolate ERA5 data to the ROMS grid for the boundary conditions<br>
- _otps2frc_wrapper.m_ to create the forcing file for the tidal motion from the OTPS tidal model<br>

### _Data_ (not present in the repository)
This subdirectory should contain all the input and forcing files, so the grid file, etc. It is not part of the repository due to size limitations. For this application you should have the input files mentioned above.

### _varinfo.yaml_ ([wiki](https://www.myroms.org/wiki/varinfo.yaml))
This config file contains information about all the different variable metadata for the input and output. [Todo: Check if varinfo was modified]

### _build_roms.sh_ ([wiki](https://www.myroms.org/wiki/build_roms))
This script is used for building the ROMS executable. After running you should find binary called _romsG_ (or _romsS_ when running in serial without MPI) in your application directory. The script depends on _aurora-0.h_ and the header files in `/Functionals/` so if either of these are changed (or you updated the source code) you should re-build the executable.

### _aurora_0.h_
This file sets which c-preprocessor flags are activated or not and which functionals are used.

### _aurora_0.in_ ([wiki](https://www.myroms.org/wiki/roms.in))
This file defines the inputs for a model run.

## Adjusting the setup
## Running



