#!/bin/tcsh
# Execute the DCMIP 1-2 test on a laptop or desktop computer

set EXE       = preqx_L64
set TEST_NAME = DCMIP1-2_NE8
set NCPU      = 8

#_______________________________________________________________________
# set path variables

set TEST_DIR  = `pwd`
cd ../../..
set HOMME     = `pwd`
make -j6 ${EXE}
if($status) exit

set  EXE_PATH = ${HOMME}/test_execs/${EXE}/${EXE}
echo HOMME    = ${HOMME}
echo EXE_PATH = ${EXE_PATH}

#_______________________________________________________________________
# create directories for simulation output

set JOBID   = `date "+%y%m%d_%H%M%S"`
set RUN_DIR = ./RUN_${TEST_NAME}_$JOBID
echo "RUN_DIR= $RUN_DIR"

cd $TEST_DIR
mkdir -p $RUN_DIR/movies

cp VCOORD/*.ascii            $RUN_DIR
cp SCRIPTS/*.ncl             $RUN_DIR
cp NAMELISTS/${TEST_NAME}.nl $RUN_DIR

#_______________________________________________________________________
# launch the executable

cd $RUN_DIR
date
echo "mpirun -n $NCPU $EXE < namelists/${TEST_NAME}.nl"
mpirun -np $NCPU $EXE_PATH < ${TEST_NAME}.nl
if($status) exit

# print timing info
cat HommeTime_stats | grep walltotal
echo "DCMIP1-2 `cat HommeTime_stats | grep prim_run`"
echo

# print error norms
ncl dcmip1-2_error_norms.ncl
echo

# plot some results
ncl dcmip1-2_lat_height.ncl
echo

exit




