# HOMME_NGGPS_benchmark

## Configuring a 13km AVEC benchmark for the High Order Methods Modeling Environment (HOMME)

  The AVEC benchmark was defined and used as part of the Next Generation Global Prediction System (NGGPS) 
dynamica core testing plan. A 13km and a 3km resoution tests case were defined.  This document describes 
how to configure the High Order Methods Modeling Environment (HOMME) for the 13km benchmark configuration.

  The following files are provided to simplify the execution of this benchmark.  

```
  README.md                A readme file to describe the benchmark
  nggps_ne32.nl            A low resolution 100km version of the benchmark to simplify testing.
  nggps_ne256.nl           A full resolution 13km version of the benchmark.
  homme1_3_26.tar.gz       A tar file containing the version of the HOMME source code.
  baroclinic_inst_mod.F90  A patched version of a source file necessary to initialize the problem.  
```


 Untar the file homme1_3_26.tar.gz which will create a directory structure the root of which will be reffered to as $HOMMEDIR.  You will need to replace the existing file src/baroclinic_inst_mod.F90 with the patched version provided inorder to initialize the problem correctly.  

 A recent version of the netcdf library needs to be installed on the platform in which you intend to install
HOMME.  You will also need a cmake version that is at least 2.8.10.2 or newer. General directions for configuring and building HOMME using cmake can be located at the following URL

    https://wiki.ucar.edu/display/homme/The+HOMME+CMake+build+and+testing+system

  Note that HOMME can be built with a wide variety of options.  We suggest the following cmake configuration 
that  explicitly build the benchmark configuration. In a separate directory different from the source code, create 
a cmake configuration script with the following contents.  Note that we refer to this directory as $BUILDDIR

```#------------------------------------------------------------
#!/bin/csh 
rm -rf CMakeFIles CMakeCache.txt

# load The necessary module files
# module load cmake
# module load netcdf

cmake -C $HOMMEDIR/cmake/machineFiles/$MACH.cmake -DQSIZE_D=10 -DPREQX_PLEV=128 -DPREQX_NP=4 \
      -DBUILD_HOMME_SWDGX=FALSE \
      -DBUILD_HOMME_SWEQX=FALSE \
      -DBUILD_HOMME_PRIMDGX=FALSE \
      -DPREQX_USE_ENERGY=FALSE preqx \
      $HOMMEDIR
#------------------------------------------------------------
```

   Note that $MACH is one of the supported cmake machine files found in the directory $HOMMEDIR/cmake/makeFILes/ and 
$HOMMEDIR is the installation path where the contents of homme1_3_26.tar.gz was placed.  Once the cmake has been configured is is possible to build the executible using the following command:

```
   make -j 6 preqx   
```


   Once the executible has been built we recommend creating the following directory $BUILDDIR/tests/NGGPS. Copy the provide namelist files nggps_ne32.nl and nggps_ne256.nl into the $BUILDDIR/tests/NGGPS directory, and create any necessary runscripts in this directory as well.  You will also need to create a "movies" and "vcoord" subdirectory in the $BUILDDIR/tests/NGGPS directory.  In the vcoord directory you will need to copy the vertical coordinate files which can be created in the following manner:
    
```
cd $HOMMEDIR/test/vcoord
make sigeqx
execute sigeqx and enter 128, which is the number of vertical levels, when prompted.  
```

Copy the two ascii files that are generated into the $BUILDDIR/tests/NGGPS/vcoord directory.  

You should now be ready to submit and run the NGGPS benchmark from the $BUILDDIR/tests/NGGPS directory. The namelist 
files for the benchmarks are configured to run for 16 hours of model time, and will generate a HommeTime_stats file at the completion of a run. The total time to complete the timestep piece of HOMME is indicated by the prim_run timer in the HommeTime_stats file.  For example:

```
  > grep prim_run HommeTime_stats 
 prim_run             1024     1024 4.915200e+05   4.824771e+06  4711.804 (   264      0)  4711.623 (   770      0)
```

   The number 4711.804 is the maximum amount of time in seconds that the run took on the slowest MPI rank.  When reporting the execution time for the NGGPS benchmark please divide this number by 8 to indcate the time it takes to execute 2 hours.  So a value of 588.97 seconds would be reported. When reporting execution time please include both the number of cores, nodes and elapsed time as well as the date and platform on which it was collected.  Some sample timings are provided below.
  

## Sample Results:
-----------------

System: Edison 
Date:   June 2015
Version: 

| Cores   | Nodes   | Elapsed time (seconds) (2h, average) |
| ------- |---------| -------------------------------------|
|6144     | 256     | 66.54 |
|12288    | 512     | 33.38 |
|24576    | 1024    | 16.34 |
|49152    | 2048    |  7.79 |
|98304    | 4096    |  3.28 |
|131072   | 5462    |  2.48 |


System:  Yellowstone 

Date:    November 2016

Version tag: dungeon12 

| Cores   |   Nodes  | Elapsed time (seconds (2h, average) |
|---------|----------|-------------------------------------|
| 1024    |  64      | 588.97 |

 
## A full description of the dynamical core and configuration is provided below for reference.

```
Model: HOMME
Numerical Method: Spectral element 
precision: real*8
Resolution: 13km (cubed sphere, each panel with 256*256 elements, polynomials of degree 3 in each element)
Grid points (horizontal physics columns): 3,538,946
degrees of freedom along the equator: 3072
Average grid spacing along equator: 12.9km
vertical layers: 128
Timestep: 40s (RK5 dynamics, RK3-SSP tracer transport)
	  120s (vertical remaping)

Simulation: J&W baroclinic instability
10 tracers, initiated in a 0/1 checkboard pattern per AVEC report 
hyperviscosity nu = 1e12
"limiter8" - quasi-monotone limiter
```
