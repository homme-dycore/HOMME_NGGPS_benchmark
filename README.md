# HOMME_NGGPS_benchmark

## Configuring a 13km AVEC benchmark for the High Order Methods Modeling Environment (HOMME)

  The AVEC benchmark was defined and used as part of the Next Generation Global Prediction System (NGGPS) 
dynamical core testing plan. A 13km and a 3km resoution tests case were defined.  This document describes 
how to configure the High Order Methods Modeling Environment (HOMME) for the 13km benchmark configuration.



### Simple configuration (*for supported systems - Cori (NERSC), Cheyenne (NCAR) & Yellowstone (NCAR)*)

```
 1. Download the code from Github (https://github.com/homme-dycore/HOMME_NGGPS_benchmark)
 2. Edit the 'setup.sh' script in the main directory:
    a) Set the system name, all lowercase (eg, 'cori', 'cheyenne' or 'yellowstone')
    b) Select which cases you wish to run (eg, 'small, medium')
 3. Run the script by providing it a directory in which to build and run (eg, './setup.sh ~/HOMME_benchmark')

  When the jobs finish, you should receive results files for each one in the directory provided above.
```

### Custom configuration 

  To run on an unsupported system requires one additional step, between steps 2 and 3 in the simple 
configuration.  A user needs to create a '<system name>.config' file in the ./systems directory of
the benchmark with the relevant information on configuration size, build, environment and run options.
Copying an existing file that is similar to your architecture and making changes is probably the 
easiest approach.  The three test configurations correspond to a total of 384 elements (small),
6144 elements (medium) and 393,216 elements (large), and the nodes & ranks in the CONFIG line should
be 

```
 1. Download the code from Github (https://github.com/homme-dycore/HOMME_NGGPS_benchmark)
 2. Edit the 'setup.sh' script in the main directory:
    a) Set the system name, all lowercase (eg, 'cori', 'cheyenne' or 'yellowstone')
    b) Select which cases you wish to run (eg, 'small, medium')
 3. Create a '<system name>.config' file in the 'systems' directory
    a) Existing tests are all provided in 1 element/core - to set the number
       of nodes, just divide the total number of ranks by your ranks per node.
 3. Run the script by providing it a directory in which to build and run (eg, './setup.sh ~/HOMME_benchmark')

  When the jobs finish, you should receive results files for each one in the directory provided above.
```

### Results

  Upon successful completion of a job, the benchmark will create a results file in the directory
in which it was created.  An example file for a small case run on Cheyenne, verifying against Yellowstone,
is shown below:

```
Verification using 't' to a tolerance of 1.0e-11:
Original:   [PASS]      t Relative Differences     Min:  2.300000e-14     Max:  1.213000e-12     Sum:  0.000000e+00     
Optimized:  [PASS]      t Relative Differences     Min:  9.000000e-15     Max:  1.068000e-12     Sum:  0.000000e+00     

Improvement: 1.5357x   Optimized:    1.288   Original:    1.978 
```


### Verification 

  The verification is based off a perturbation test using a machine epsilon
of ~2.22e-16 in a run on the 'Yellowstone' (NCAR) system.  This was run using 
the medium-sized case, and resulted in relative differences in temperature from 
the unperturbed cases of < 1e-12, as seen below:

```
Original Code:
    t Relative Differences     Min:  9.000000e-15     Max:  3.060000e-13    Sum:  0.000000e+00

Optimized Code:
    t Relative Differences     Min:  1.400000e-14     Max:  1.060000e-13    Sum:  0.000000e+00
```

  The automated testing uses a slightly higher value, 1e-11, to verify
results.  This passes on Cheyenne, Yellowstone and Cori using the Intel 17
compilers.

### Advanced use

  The 'setup.sh' script creates build scripts in the build directories, and run scripts with preset values in
the case directories.  The build scripts can be re-run if code changes are introduced and the scripts will 
still properly link to the executables, and the run scripts can be modified to run on different numbers of
ranks/nodes to establish scaling numbers.  


### Requirements

 A recent version of the netcdf library needs to be installed on the platform in which you intend to install
HOMME.  You will also need a cmake version that is at least 2.8.10.2 or newer.  The scripts make use of a bash
shell and GNU 'bc'.

General directions for configuring and building HOMME using cmake can be located at the following URL

    https://wiki.ucar.edu/display/homme/The+HOMME+CMake+build+and+testing+system

### Sample Results
-----------------

The sample results show timings for various versions of the optimized code.  The current version is 
dungeon28.  The previous results are in reverse chronological order.

System: Cori
Date:  January 2017
Version tag: dungeon16

| HwCores   |   Nodes  | Elapsed time (seconds (2h, average) |
|---------|----------|-------------------------------------|
| 4096    |  64      | 158.92 |
| 6144    |  96      | 103.53 |
| 8192    |  128     |  77.25 |
| 16384   |  256     |  41.13 |
| 24576   |  384     |  28.06 |
| 32768   |  512     |  23.64 |
| 49152   |  768     |  14.60 |
| 65536   |  1024    |  11.33 |
| 98304   |  1536    |   7.87 |
| 196608  |  3072    |   4.57 |
| 393216  |  6144    |   2.87 |


System: Cheyenne
Date: January 2017
Version tag: dungeon15

| HwCores   |   Nodes  | Elapsed time (seconds (2h, average) |
|-----------|----------|-------------------------------------|
| 6144      |  171     | 96.72 |
| 8192      |  228     | 83.32 |
| 12288     |  342     | 55.90 |
| 16384     |  456     | 41.37 |
| 24576     |  683     | 26.89 |
| 32768     |  911     | 19.22 |
| 49152     |  1366    | 11.19 |
| 65536     |  1821    |  7.75 |
| 98304     |  2731    |  5.42 |
 

System:  Yellowstone
Date:    November 2016
Version tag: dungeon12

| HwCores   |   Nodes  | Elapsed time (seconds (2h, average) |
|---------|----------|-------------------------------------|
| 1024    |  64      | 588.97 |


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
