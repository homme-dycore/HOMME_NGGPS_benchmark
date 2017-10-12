#!/bin/tcsh -f
#
###############################################################
# Explicit leapfrog default benchmark
###############################################################
#
# Spectral Element -- swtc5
# NE=30, dt=90, nu=1.5e15, limiter=0, filter_freq=0, NP=4
#
###############################################################

# start_dir   = reg_test directory (used to copy namelist and
#               link to namelist)
set start_dir = $PWD/..

# out_dir = where output will be saved
if ( $?HOMME_OUT ) then
	set out_dir = $HOMME_OUT/sweqx/swtc5
	echo Found \$HOMME_OUT, output will be stored in $out_dir
else
	set out_dir = ~/scratch1/reg-test/sweqx/swtc5
	echo Did not find \$HOMME_OUT, output will be stored in
	echo $out_dir instead
endif

echo

# bld_dir = $HOMME/build/sweqx
if ( ! $?HOMME_ROOT ) then
        set HOMME_ROOT = $PWD/../../..
	echo \$HOMME_ROOT is undefined, using $HOMME_ROOT
else
	if ( ! -d $HOMME_ROOT/build/sweqx ) then
		echo $HOMME_ROOT does not contain build system for HOMME,
		set HOMME_ROOT = ~/codes/homme
		echo using $HOMME_ROOT instead
		echo
	endif
endif

if ( -d $HOMME_ROOT/build/sweqx ) then
	set bld_dir = $HOMME_ROOT/build/sweqx
else
	echo $HOMME_ROOT does not contain build system for HOMME! 
	exit
endif

echo Building from $HOMME_ROOT

set NCPU = 16
if ( ${?PBS_NODEFILE} ) then
   set NCPU = `wc $PBS_NODEFILE | awk '{print $1}' - `
endif
if ( ${?SLURM_NNODES} ) then
   set NCPU = ${SLURM_JOB_NUM_NODES}
   set cores = 8
   set MPI_RUN = "mpiexec --npernode $cores numa_wrapper --ppn $cores"
else
	set MPI_RUN = "mpirun -np $NCPU"
endif
echo NCPU = $NCPU

# Configure and Make
# MNL note: don't need to completely clean / rebuild if config.h
#           and dependencies are already set correctly...
cd $bld_dir

# CFG contains flags to use w/ ./configure
set CFG = "--enable-blas --enable-lapack --with-netcdf=$NETCDF_PATH --with-pnetcdf=$PNETCDF_PATH NP=4 PLEV=1"
if ( -e config.status ) then
	echo config.status exists, running it!
	./config.status
	set NP = `sed -n 's/#define NP \([0-9]\)/\1/p' config.log`
	set PLEV = `sed -n 's/#define PLEV \([0-9]\)/\1/p' config.log`
	set PIO = `sed -n 's/#define PIO \([0-9]\)/\1/p' config.log`
	echo NP = $NP
	echo PLEV = $PLEV
	echo PIO = $PIO
	if ( ( $NP == 4 ) && ( $PLEV == 1 ) && ( $PIO != 1) ) then
		echo Already configured with NP = 4, PLEV = 1, and PIO_INTERP...
		echo Skipping configure stage
		rm -f sweqx
		make -j 4
	else
		echo Need to reconfigure / rebuild to change NP and PLEV
		make distclean
		./configure $CFG 
		make -j4 depends ; make clean
		make -j 4
	endif
else
	echo No config.status exists, need to configure / build from start
	./configure $CFG 
	make -j4 depends ; make clean
	make -j 4
endif

# Create directories for output
mkdir -p $out_dir/movies
cp $start_dir/namelists/swtc5.nl $out_dir/input.nl
cd $out_dir
rm -f swtc5ref.ncl
ln -s $start_dir/ncl/swtc5ref.ncl ./

# Run
date
$MPI_RUN $bld_dir/sweqx < input.nl | tee sweqx.out
date
ncl swtc5ref.ncl | tee -a sweqx.out
