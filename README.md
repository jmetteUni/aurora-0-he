# Aurora Hydrothermal Venting Application<br> for the Regional Ocean Modelling System (ROMS)

# Table of contents
- [Description](#Description)
- [How to run the model](#How-to-run-the-model)
    - [Prerequisites](#Prerequisites)
    - [Preparing the input files](#Preparing-the-input-files)
    - [Setting up the model and directory structure](#Setting-up-the-model-and-directory-structure)
    - [Running](#Running)
    - [Output](#Output)
- [References](#References)

# Description
This application is used for the simulation of the Aurora hydrothermal vent system in the Arctic Ocean. It is intended to be used with a modified version of the Regional Ocean Modelling System (ROMS) (documented [here](https://github.com/jmetteUni/roms/tree/bottom-tracer) and [OpenMPI](https://www.open-mpi.org/).

The application uses a 256x256 horizontal grid with 30 vertical layers. Initial and boundary conditions are inferred from the ERA5 Interim reanalysis product. Surface fluxes are set to zero. Input of heat flux and a passive tracer, to simulate the hydrothermal vent is prescribed with forcing files. Tides are simulated based on tidal constituents from the Arctic Ocean Tidal Inverse Model. The simulation is run with a timestep of 4 sec for one month.

It is heavily inspired and uses a similar approach as the setup in _Xu and German, 2023_[^1]. G. Xu also helped extensively with this setup here. The simulation results were published in the master thesis _Plume Dispersal in the Arctic Ocean_[^2].

For general information on the ROMS model check the official [wiki](https://www.myroms.org/wiki/Documentation_Portal) (I linked some important articles directly below) and the [forum](https://www.myroms.org/forum/), but note that the wiki is outdated in some places. Most of the files have also some documentation as comments inside.

# How to run the model
This is a tutorial on how to use and run this application with the ROMS model. For more in-depth documentation see the official ROMS documents. The application in it's structure follows roughly the test case [UPWELLING](https://www.myroms.org/wiki/UPWELLING_CASE). The test cases can be obtained [here](https://github.com/myroms/roms_test). It is recommended that you test your model with the UPWELLING test case first.

## Prerequisites ([wiki](https://www.myroms.org/wiki/Getting_Started))
Software you will need to have installed:<br>
- git<br>
- (an implementation of MPI for running on multiple cores, recommended)<br>
- fortran90 or fortran95 compiler<br>
- gnu make version 3.81 or higher<br>
- cpp for C-processing<br>
- Perl<br>
- NetCDF<br>
- the ROMS MATLAB tools [^6] to prepare the input files.

## Preparing the input files

The input files are prepared with MATLAB. Most of them rely on functions from the ROMS MATLAB collection [^6], so it should be added to the MATLAB paths.

### Grid
There are various ways the a grid can be created. In this case it was done using the MATLAB toolbox _GridBuilder_[^5]. It is interactive and GUI based. The grid can be saved either as an editable `.mat` file or as NetCDF for the use in the ROMS setup. A detailed manual is also available on the website of [^5]. The `.mat` file is provided in the `matlab/` directory. The parameters used for the grid creation are documented and explained in [^1]. The grid uses mutlibeam bathymetry data of the area as a base which is imported in GridBuilder. The bathymetry data is based on Warnke et al., 2024 [^4].

### Bottom forcing

The bottom forcing representing the vent is constructed in MATLAB using the `get_bottom_forcing_jm.m` script. It uses the grid NetCDF file as an input and outputs two NetCDF files, one for the heat flux and one for the passive tracer flux. The paths, the vent position, the start and end time and the amount of energy flux and tracer material (called "dye_amount") has to be specified at the top of the script. Additionally it needs the grid parameters used during the grid creation (see above). The energy flux and the dye amount is the constant input provided to the model domain per timestep.

By default, the flux is only taking place in the grid cell, where the specified position is located. But the script also has options for including several vents at different positions. There is also code present for adding the volume flux, but this is not fully functional and was not used in this setup.

In principal the time period of the forcing file can be longer (here: one year) then the model run lasts (here: one month).

### Initial conditions

The initial conditions prescribe the state inside the model domain at the start of the run. They are interpolated from the GLORYS reanalysis product [^3] onto the roms grid, so it also needs the grid as an input. The paths, the start time, the reference time and the grid parameters have to be provided at the top of the script. It is easiest to set the start and reference time equal but for other use cases it can be different. For more information on the reference time see the [wiki](https://www.myroms.org/wiki/time_ref).

### Boundary condiitons

The boundary conditions prescribe the state at the margins of the domain during the run. How they are effecting the domain is specified by the boundary condition options in the `aurora.in` file (see below and in the [wiki](https://www.myroms.org/wiki/Boundary_Conditions)). They are interpolated from the GLORYS reanalysis product [^3] to the outer boundaries of the ROMS grid, so it also needs the grid file as an input. The paths, the start and end time and the grid parameters have to be provided at the top of the script.

The variabel `nemo_time` sets the amount of timestemps in the output forcing file.
> [!CAUTION]<br>
> In the given script `nemo_time = 14` also equals the number of timesteps in the reanalysis product. This way interpolation of the reanalysis data is avoided, but due to the model run lasting for only one month but the 14 time steps are given for one year of reanalysis data, this results in only two boundary states during the model run. It is recommended, that the script should be modified to use reanalysis data with an appropriate time resoultion corresponding to the planned model run.

In principal the time period of the forcing file can be longer (here: one year) then the model run lasts (here: one month).

### Tidal forcing

This provides tidal forcing for the domain during the model run. It is enacted as changes in sea surface elevation. The forcing file is created with a wrapper script which needs the grid file as an input and depends on the TMD MATLAB toolbox and the OSU tidal prediction software. The process is explained in detail [here](https://www.myroms.org/wiki/Tidal_Forcing) in the ROMS wiki.

## Setting up the model and directory structure
There are many different ways how to structure your source code and your application. This way here is not neccesarily the best, but the way it worked for me.

1. Clone the modified ROMS [source code](https://github.com/jmetteUni/roms/tree/bottom-tracer) in one directory.
2. Clone this application into a different directory with a characteristic application name.

In the following part, the different subdirectories and files which make up the application are briefly explained:<br>
### _Functionals_ ([wiki](https://www.myroms.org/wiki/Functionals))
This subdirectory contains header (.h) files. These are copied from the source code under `ROMS/Functionals/` and are therefore mostly similar. But some are modified specifically for this setup such as `ana_srflux.h` for example, which specifies the analytical (not prescribed from forcing files) surface forcing. In this file surface forcings are set to 0 for this application to represent a simplified ice cover. The directory also contains header files which are not used by the application (because there function is replaced by forcing files or they're just not neccesary) but for simplicity just all files from the source code are copied and they don't hurt. In case you update the source code, it is necessary to check if you also need to update header files in your application.

### _matlab_
This subdirectory contains matlab scripts which are used to create the input forcing files. The forcing files are NetCDF files which prescribe certain conditions for the model instead of conditions from analytical functions (the header files). They are used for example to simulate realistic boundaries (That is why the UPWELLING test case used no NetCDF files) and therefore can get quite large in terms of disk space. This application uses four different scripts and also provides a one `.mat`-file for the grid. One additional script is the example of how the ROMS MATLAB tools need to be modified for preparing the input files.<br>
- `grid-M256L256.mat` is a MATLAB file representing the grid of the model domain.<br>
- `get_bottom_forcing_jm.m` to create the forcing file for the input of the passive bottom tracer.<br>
- `get_initial_condition_from_GLORYS_jm.m` to interpolate ERA5 data to the ROMS grid for the initial conditions.<br>
- `get_boundary_condition_GLORYS_jm.m` to interpolate ERA5 data to the ROMS grid for the boundary conditions.<br>
- `otps2frc_wrapper.m` to create the forcing file for the tidal motion from the OTPS tidal model.<br>
- `roms_metadata.m` which should replace the same file in the MATLAB tools in `/roms_matlab/netcdf/roms_metadata.m`.

### _Data_ (not present in the repository)
This subdirectory should contain all the input and forcing files, so the grid file, etc. It is not part of the repository due to size limitations. For this application you should have the following input files:<br>
- grid-M256L256.nc<br>
  The roms grid file with 256 × 256 cells, constructed with the MATLAB toolbox GirdBuilder: https://austides.com/downloads/. <br>
- grid-M256L256_bhflux_180MW_20220726T00.nc<br>
  The forcing file providing the bottom heat flux representing the vent heat input with 130MW. <br>
- grid-bpflux_20220726T00.nc<br>
  The forcing file prociding the passive tracer input representing helium istope concentration ratio delta3He. <br>
- grid-M256L256_bry_nemo_20220726T00_20230701T00.nc<br>
  Forcing file for boundary values around the model domain, interpolated from GLORYS reanalysis data. <br>
- grid-M256L256_ini_nemo_20220726T00_64layer.nc<br>
  Initial conditions inside the model domain, interpolated from GLORYS reanalysis data. <br>
- grid-M256L256_tides<br>
  Forcing file for the tidal constituents, with data from the OSUS tidal prediction model. <br>

### _varinfo.yaml_ ([wiki](https://www.myroms.org/wiki/varinfo.yaml))
This config file contains information about all the different variable metadata for the input and output. This is in principal the same file as in the ROMS source code in `/roms/ROMS/Externals/varinfo.yaml`. It doesn't matter if you use this one or the one from the source code, as long as they are equal and the correct path is specified in the `aurora_0.in` file. Although when the source code gets updated it is recommended to also copy the new version of the file to the application as it could change.

### _build_roms.sh_ ([wiki](https://www.myroms.org/wiki/build_roms))
This script is used for building the ROMS executable. After running you should find binary called `romsG` (this may be named slightly different in different setups) in your application directory. The script depends on `aurora-0.h` and the header files in `/Functionals/` so if either of these are changed (or if you updated the source code) you should re-build the executable. Important sections here are:

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

This sets the compiler options. `USE_MPI` controlles, if the model is running with MPI and below you can set the used MPI application. `FORT` sets your fortran compiler. For all this options you can just uncomment the one suiting your installation.

### _aurora_0.h_
This file sets, which c-preprocessor flags are activated or not and which functionals are used. They are either activated (`#define`) or deactivated (`#undef`), similar effect as if you would delete it). Most of these flags are also used as in the setup by G. Xu. Some particular important ones are explained now.

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
This file defines the inputs and all other settings for a model run. In this section the most important adjustments are highlighted and explained. Otherwise this follows also mostly the example file in the UPWELLING test case and the setup by G. Xu. This file is extensively documented with inline comments and every available setting option is explained in an extra section at the end of the file.

    TITLE = Aurora hydrothermal vent setup, with helium tracer
    MyAppCPP = AURORA0

This names and describes your application

    VARNAME = ./varinfo.yaml

Here you have to specify the path to your _varinfo.yaml_ file (see above).

    Lm == 254            ! Number of I-direction INTERIOR RHO-points
    Mm == 254            ! Number of J-direction INTERIOR RHO-points
     N == 30             ! Number of vertical levels

This sets the grid dimensions of your application. They have to correspond to your grid input file where `Lm = dimension in your grid file - 2` and the same for `Mm`. The dimension of your grid file should ideally scale properly with the tiling variables (see below) and the amount of MPI processes your application will run with.

      NtileI == 8                           ! I-direction partition
      NtileJ == 4                           ! J-direction partition

The tiling variables, which partition your grid into subsets for parallel computation with MPI. It is required that `NtileI * NtileJ = N ` where N is the number of MPI processes used when running the setup.

    NPT =  1             ! Number of inactive passive tracers

This activates the passive tracer representing helium input from the vent.

From line `112` and onwards you can set the advection algorithms and the lateral boundary condition types. They are set as in the setup by Xu and German, 2023. Some of the lateral boundary conditions require a NetCDF input file for the boundary condiitons while others (like Clo = Closed or Per = Periodic) don't. If you want to test your setup with only analytical boundary conditions you have to adjust these. See for comparison the UPWELLING test case and the [wiki](https://www.myroms.org/wiki/Boundary_Conditions).

    NTIMES == 669600            !15sec×60min×24h×31d
    DT == 4.0d0                 !4sec
    NDTFAST == 30

Here the first line sets the number of baroclinic timesteps and the second line the delta of the time step, so combined the total time the model will run. `NDTFAST` sets the number of barotropic time steps between each baroclinic time step.

    NRST == 43200
    ...
    NAVG == 900
    ...
    NDIA == 1800

These settings define how often (in number of timesteps) data is written to the output file, namely the RESTART file, the AVERAGE file and the DIAGNOSTICS file. Usually the AVG file is the one you would use to analyse the results.

    Vtransform == 2             ! transformation equation
    Vstretching == 4            ! stretching function

    THETA_S == 0.0d0            ! surface stretching parameter
    THETA_B == 3.0d0            ! bottom  stretching parameter
    TCLINE == 20.0d0            ! critical depth (m)

This are properties of your grid. They have to match the ones you used when generating your grid file and the other input files.

    TIME_REF =  20220726.0d0    ! yyyymmdd.dd

This sets a reference time for your application. If a date is used here all other dates in the input files (for example tidal forcing) have to correspond to this reference date.

     GRDNAME == ./Data/grid-M256L256.nc
     ININAME == ./Data/grid-M256L256_ini_nemo_20220726T00_64layer.nc
     ! ININAME == roms_ini.nc
     ...
     BRYNAME == ./Data/grid-M256L256_bry_nemo_20220726T00_20230701T00.nc
     ...
     TIDENAME == ./Data/grid-M256L256_tides.nc

This sets the NetCDF files for the grid, the initial condiitons, the boundary condtions and the tidal forcing files. If the application should be run without them (e.g. for testing) set the parameter as in the option which is commented out.

     NFFILES == 2                                                       ! number of unique forcing files
     FRCNAME == ./Data/grid-M256L256_bhflux_130MW_20220726T00.nc \      ! forcing file 1, grid 1
                ./Data/grid-M256L256_bpflux_20220726T00.nc

This sets the forcing files for the bottom flux. It is seperated into the heatflux (bhflux) and the passive tracer flux (bpflux), so `NFFILES == 2`.

## Running
After installing the prerequisites and setting everything up as described above, you start by building the excutable. Navigate into the application directory, where the build script `build_roms.sh` lives, and run it with

    ./build_roms.sh

Optionally you can speed up the build using

    ./build_roms.sh -j N

where `N` is the number of processes.

If the build was succesfull you will find a ROMS executable in the application directory. In this case it will be named _romsG_ but for example if you are building withouth MPI (or for other cases) it can be named slightly different. Then you can run the model with

    mpirun -np N romsG aurora-0.in

using MPI, where `N`is again the number of MPI processes, `romsG` the executable and `aurora-0.in` the input script. `N` has to be the product of the tiling variables in the `aurora-0.in` file (see above). The header file `aurora-0.h` is not called directly but was already used in the build process.

## Output

The model produces several output files. The most relevant one is the `roms_avg.nc` which holds averaged values of all relevant physical quantities. For an inspection, if it contains the espected data I would recommend to have ncview installed on the computing host, so you can check quickly if everything worked. For more detailed diagnostics the `roms_dia.nc` file contains more data and more model parameters. The `roms_rst.nc` file contains a full state of the model, and can be used to restart a model run from this point.

While running, the model prints a lot of useful information to the console (see [wiki](https://www.myroms.org/wiki/Standard_Output)), especially for debugging, if a run does not work as espected. It can be useful to save this output to a text file or something similar. It is also helpful to run the model with something like the terminal manager [screen](https://www.gnu.org/software/screen/manual/screen.html), to not accidently stopping your model run.

# References

[^1]: Xu G and German CR (2023) Dispersion of deep-sea hydrothermal plumes at the Endeavour Segment of the Juan de Fuca Ridge: a multiscale numerical study. Front. Mar. Sci. 10:1213470. doi: https://doi.org/10.3389/fmars.2023.1213470
[^2]: Mette, Jonathan. “Plume Dispersal in the Arctic Ocean - the Aurora Site at Gakkel Ridge,” May 16, 2025. https://media.suub.uni-bremen.de/handle/elib/22668.
[^3]: Global Ocean Physics Reanalysis, E.U. Copernicus Marine Service Information (CMEMS). Marine Data Store (MDS). https://doi.org/10.48670/moi-00021
[^4]: Warnke, Fynn; Unland, Ellen; Höppner, Laura; Dorschel, Boris; Dreutter, Simon; Schlindwein, Vera (2024): Multibeam bathymetry raw data (Atlas Hydrosweep DS 3 echo sounder entire dataset) of RV POLARSTERN during cruise PS137 \[dataset\]. PANGAEA, https://doi.org/10.1594/PANGAEA.963721
[^5]: James, C. : GridBuilder v1.3.3. https://austides.com/downloads/. Version: March  2023
[^6]: Arango, H. G.: ROMS Matlab Processing Scripts. https://github.com/myroms/roms_matlab. Version: Oct. 2023
