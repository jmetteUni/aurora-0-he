# Aurora Hydrothermal Venting Application<br> for the Regional Ocean Modelling System (ROMS)

## Description
This application is used for the simulation of the Aurora hydrothermal vent system in the Arctic Ocean. It is intended to be used with a modified version of the Regional Ocean Modelling System (ROMS) (documented [here](https://github.com/jmetteUni/roms/tree/bottom-tracer) and [OpenMPI](https://www.open-mpi.org/).

The application uses a 256x256 horizontal grid with 30 vertical layers. Initial and boundary conditions are inferred from the ERA5 Interim reanalysis product. Surface fluxes are set to zero. Input of heat flux and a passive tracer, to simulate the hydrothermal vent is prescribed with forcing files. Tides are simulated based on tidal constituents from the Arctic Ocean Tidal Inverse Model. The simulation is run with a timestep of 4 sec for one month.

It is heavily inspired and uses a similar approach as the setup in [Xu and German, 2023](http://dx.doi.org/10.3389/fmars.2023.1213470). G. Xu also helped extensively with this setup here. The simulation results were published in the master thesis [Plume Dispersal in the Arctic Ocean](https://doi.org/10.26092/elib/4441).

For general information on the ROMS model check the official [wiki](https://www.myroms.org/wiki/Documentation_Portal) (I linked some important articles directly below) and the [forum](https://www.myroms.org/forum/), but note that the wiki is outdated in some places. Most of the files have also some documentation as comments inside.

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
- NetCDF<br>
- the [ROMS MATLAB tools](https://github.com/myroms/roms_matlab/tree/main) to prepare the input files.

## Setting up the model and directory structure
There are many different ways how to structure your source code and your application. This way here is not neccesarily the best, but the way it worked for me.

1. Clone the modified ROMS [source code](https://github.com/jmetteUni/roms/tree/bottom-tracer) in one directory.
2. Clone this application into a different directory with a characteristic application name.

In the following part, the different subdirectories and files which make up the application are briefly explained:<br>
### _Functionals_ ([wiki](https://www.myroms.org/wiki/Functionals))
This subdirectory contains header (.h) files. These are copied from the source code under `ROMS/Functionals/` and are therefore mostly similar. But some are modified specifically for this setup such as _ana_srflux.h_ for example, which specifies the analytical (not prescribed from forcing files) surface forcing. In this file surface forcings are set to 0 for this application to represent a simplified ice cover. The directory also contains header files which are not used by the application (because there function is replaced by forcing files or they just not neccesary) but for simplicity just all files from the source code are copied and they don't hurt. In case you update the source code, it is necessary to check if you also need to update header files in your application.

### _matlab_
This subdirectory contains matlab scripts which are used to create the input forcing files. The forcing files are NetCDF files which prescribe certain conditions for the model instead of conditions from analytical functions (the header files). They are used for example to simulate realistic boundaries (That is why the UPWELLING test case used no NetCDF files) and therefore can get quite large in terms of disk space. This application uses four different scripts and one _.mat_-file. One additional script is the example of how the ROMS MATLAB tools need to be modified for preparing the input files.<br>
- _grid-M256L256.mat_ is a _MATLAB_ file representing the grid of the model domain.<br>
- _get_bottom_forcing_jm.m_ to create the forcing file for the input of the passive bottom tracer.<br>
- _get_initial_condition_from_GLORYS_jm.m_ to interpolate ERA5 data to the ROMS grid for the initial conditions.<br>
- _get_boundary_condition_GLORYS_jm.m_ to interpolate ERA5 data to the ROMS grid for the boundary conditions.<br>
- _otps2frc_wrapper.m_ to create the forcing file for the tidal motion from the OTPS tidal model.<br>
- _roms_metadata.m_ which should replace the same file in the _MATLAB_ tools in `/roms_matlab/netcdf/_roms_metadata.m`.
### _Data_ (not present in the repository)
This subdirectory should contain all the input and forcing files, so the grid file, etc. It is not part of the repository due to size limitations. For this application you should have the input files mentioned above.

### _varinfo.yaml_ ([wiki](https://www.myroms.org/wiki/varinfo.yaml))
This config file contains information about all the different variable metadata for the input and output. This is in principal the same file as in the ROMS source code in `/roms/ROMS/Externals/varinfo.yaml`. Although when the source code gets updated it is recommended to also copy the new version of the file to the application as it could change.

### _build_roms.sh_ ([wiki](https://www.myroms.org/wiki/build_roms))
This script is used for building the ROMS executable. After running you should find binary called _romsG_ (or _romsS_ when running in serial without MPI) in your application directory. The script depends on _aurora-0.h_ and the header files in `/Functionals/` so if either of these are changed (or if you updated the source code) you should re-build the executable. Important sections here are:

    export     MY_PROJECT_DIR=${PWD}

This sets the project directory as the directory, where the file is executed. If this build script lives in the base directory of your application you should be fine.

    export     MY_ROMS_SRC=${MY_ROOT_DIR}/roms

This sets the directory of the model source code.

    export           USE_MPI=on            # distributed-memory parallelism
    export        USE_MPIF90=on            # compile with mpif90 script
    #export         which_MPI=intel         # compile with mpiifort library
    export         which_MPI=mpich         # compile with MPICH library
    #export         which_MPI=mpich2        # compile with MPICH2 library
    #export         which_MPI=mvapich2      # compile with MVAPICH2 library
    #export         which_MPI=openmpi       # compile with OpenMPI library

    #export        USE_OpenMP=on            # shared-memory parallelism

    #export              FORT=ifort
    export              FORT=gfortran
    #export              FORT=pgi

This sets the compiler options. `USE_MPI` controlles, if the model is running with MPI and below you can set the used MPI application. `FOR` sets your fortran compiler. For all this options you can just uncomment the one suiting your installation.

### _aurora_0.h_
This file sets which c-preprocessor flags are activated or not and which functionals are used. They are either activated (`#define`) or deactivated (`#undef`), similar effect as if you would delete it). Most of these flags are also used as in the setup G. Xu. Some particular important ones are explained now.

    #undef ANA_GRID
    #undef ANA_INITIAL
    #undef ANA_BTFLUX`

This block (de-)activates the use of analytical grid, initial conditions and bottom flux forcing. If deactivated they need to be supplied by NetCDF-files in the `.in`-file, if not they are constructed from the analytical functions in the header-files. This options can be especially useful for testing your setup and your input files.

    #define T_PASSIVE
    #undef ANA_PASSIVE

This two lines activate the use of the custom passive tracer and deactivates the use of analytical passive tracer conditions respectively. As above they then have to be provided as a input file.

    #define TIDE_GENERATING_FORCES
    #if defined TIDE_GENERATING_FORCES
    #define SSH_TIDES
    #define UV_TIDES

    #define ADD_FSOBC
    #define ADD_M2OBC
    #undef RAMP_TIDES
    #endif

This block turns tidal forcing on or off (again, similar principle as above, except there are now analytical tidal forcings). Note the if-condition, which activates the other flags.

### _aurora_0.in_ ([wiki](https://www.myroms.org/wiki/roms.in))
This file defines the inputs and all other settings for a model run. In this section the most important adjustments are highlighted and explained. Otherwise this follows also mostly the example file in the UPWELLING test case. This file is extensively documented with inline comments and every available setting option is explained in the section at the end of the file.

    TITLE = Aurora hydrothermal vent setup, with helium tracer
    MyAppCPP = AURORA0  !!! Check if this makes is required

This names and describes your application

    VARNAME = ./varinfo.yaml

Here you have to specify the path to your _varinfo.yaml_ file (see above).

    Lm == 254            ! Number of I-direction INTERIOR RHO-points
    Mm == 254            ! Number of J-direction INTERIOR RHO-points
     N == 30             ! Number of vertical levels

This sets the grid dimensions of your application. They have to correspond to your grid input file where `Lm = dimension in your grid file - 2` and the same for Mm. The dimension of your grid file should ideally scale properly with the amount of MPI processes your application will run with.

From line `112`and onwards you can set the advection algorithms and the lateral boundary condition types. They are set as in the setup by Xu and German, 2023. Some of the lateral boundary conditions require a NetCDF input file for the boundary condiitons while others (like Clo = Closed or Per = Periodic) don't. If you want to test your setup with only analytical boundary conditions you have to adjust these. See for comparison the UPWELLING test case and the [wiki](https://www.myroms.org/wiki/Boundary_Conditions).

    NTIMES == 669600			!15sec×60min×24h×31d
    DT == 4.0d0     			!4sec
    NDTFAST == 30

Here the first line set the number of baroclinic timesteps and the second line the delta of the time step, so combined the total time the model will run. `NDTFAST` sets the number of barotropic time steps betweean each baroclinic time step.

    NRST == 43200
    ...
    NAVG == 900
    ...
    NDIA == 1800

These settings define how often (in number of timesteps) data is written to the output file, namely the RESTART file, the AVERAGE file and the DIAGNOSTICS file. Usually the AVG file is the one you would use to analyse the results.



## Adjusting the setup
## Running
## Output


