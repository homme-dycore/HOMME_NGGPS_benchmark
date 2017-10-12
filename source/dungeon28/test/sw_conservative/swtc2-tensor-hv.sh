#!/bin/bash -f

set -e

#MESH FILE LOCATION IS SET IN *.NL FILE

tstep=30 
nu=1.2e-6 #3.5e-8 for hv_scaling=3.2, 1.2e-6 for hv_scaling=4.0
hvscaling=4

test_case=swtc2

#CHANGE PATHS TO YOUR LOCALS
wdir=$HOME/runhomme/sweqx
input=$HOME/homme7/test/sw_conservative

NCPU=46

name=${test_case}-tensor

cd ${input}/../../build/sweqx
./configure PLEV=1 NP=4 --with-netcdf=$NETCDF_PATH  --enable-blas --enable-lapack
	make  -j 4 depends
	make clean

rm -f sweqx
make -j 4

cd $wdir


let sfreq=24*3600
sfreq=`echo "$sfreq / $tstep" | bc`

	sed s/tstep.\*/"tstep = $tstep"/  $input/swtc2-tensor-hv.nl |\
	sed s/hypervis_scaling.\*/"hypervis_scaling = $hvscaling"/  |\
	sed s/nu=.\*/"nu= $nu"/  |\
	sed s/nu_s=.\*/"nu_s= $nu"/  |\
	sed s/statefreq.\*/"statefreq = $sfreq"/  \
	> swtc2.nl

	time mpirun -np $NCPU  ${input}/../../build/sweqx/sweqx < swtc2.nl | tee  sweq.out

	mv -f sweq.mass $name.mass
	mv -f sweq.out $name.out
	mv -f movies/swtc21.nc movies/$name.nc


exit
