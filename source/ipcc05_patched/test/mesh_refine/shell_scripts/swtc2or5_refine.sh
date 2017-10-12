#!/bin/tcsh -f
#XPBS -l nodes=100:ppn=4
#PBS -l nodes=200:ppn=2
#PBS -l walltime=1:00:00
#PBS -N swtc1
#PBS -j oe
#PBS -A FY081407
#PBX -a 1758
#XXX -W depend=afterany:jobid


###############################################################
#script for running swtc2 or swtc5 with mesh refinement
#be careful to change wdir, input, builddir, AND the whole path to grid
###############################################################


set wdir = ~/runhomme/sweqx
set input = ~/homme7/test/mesh_refine/shell_scripts
set builddir = ~/homme7/build/sweqx
set NCPU = 46
echo NCPU = $NCPU

cd builddir

#with pio
#./configure PLEV=6 NP=4 --with-netcdf=$NETCDF_PATH  --with-pnetcdf=$PNETCDF_PATH --enable-pio  --enable-blas --enable-lapack

#with pio or interp and MESH
#./configure PLEV=1 NP=4 --with-netcdf=$NETCDF_PATH  --with-pnetcdf=$PNETCDF_PATH --enable-blas --enable-lapack --enable-mesh-refine #--enable-pio
#make  -j 4 depends
#make clean

rm -f sweqx
make -j 4

if ( $status ) exit

mkdir $wdir
cd $wdir
mkdir movies

# defaults:

# output units: 0,1,2 = timesteps, days, hours
set OUTUNITS = 2 #freq in hours
set OUTFREQ = 10 #10 hours

set test_case = swtc2
set ndays = 12

set NE = 0
set grid="\/home\/onguba\/homme7\/test\/mesh_refine\/grids\/mountain_10_x2"
set tstep = 150
set nu = 0.01 #default hv scaling is 3.2 in the namelist: nu=0.01 for hv_power=3.2 and nu=4.5e-7 for hv_power=4 were tested
set nu_s = 0
set limiter = 0

set name = ${test_case}

set sfreq=6
@ sfreq *= 3600
set sfreq=`echo "$sfreq / $tstep" | bc`


sed s/ne=.\*/"ne = $NE"/  $input/../namelists/swtc_2and5.nl |\
sed s/tstep.\*/"tstep = $tstep"/  |\
sed s/ndays.\*/"ndays = $ndays"/  |\
sed s/test_case.\*/"test_case = ${test_case}"/  |\
sed s/limiter_option.\*/"limiter_option = $limiter"/  |\
sed s/nu=.\*/"nu= $nu"/  |\
sed s/nu_s=.\*/"nu_s= $nu_s"/  |\
sed s/statefreq.\*/"statefreq = $sfreq"/  |\
sed s/mesh_file.\*/"mesh_file = \'${grid}.g\'"/ \
> swtc.nl

#

date
rm -f movies/swtc?.nl
mpirun -np $NCPU ~/homme7/build/sweqx/sweqx < swtc.nl | tee  sweq.out
date

mv -f sweq.mass $name.mass
mv -f sweq.out $name.out
mv -f movies/${name}1.nc movies/$name.nc
#  ncl $input/geo.ncl


