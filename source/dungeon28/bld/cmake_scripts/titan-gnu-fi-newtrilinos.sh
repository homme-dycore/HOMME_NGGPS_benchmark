#!/bin/csh

##### Before running:
#module unload PrgEnv-gnu 

#module load cmake/2.8.11.2
#module load PrgEnv-gnu
#module load gcc/4.9.0
#module load cray-shmem
#module load cray-mpich
#module load cray-netcdf-hdf5parallel/4.3.3.1
#module load cray-parallel-netcdf/1.5.0
#module load python
#module load p-netcdf/1.5.0
#module load cray-hdf5-parallel/1.8.14

# export CRAYPE_LINK_TYPE='dynamic' 

# do either
#setenv Trilinos_DIR ${CRAY_TRILINOS_PREFIX_DIR}
# or set Trilinos_DIR to your own trilinos install
setenv HOMME_ROOT /ccs/home/$USER/trunk
setenv HOMME_PROJID cli106
setenv Trilinos_DIR /lustre/atlas/world-shared/cli900/cesm/software/Trilinos/Trilinos-11.12.1_gptl/titan-gnu-ci-nophal_4.9.0/install

# also, set environment variable HOMME_ROOT to point to 
# your HOMME repo, or change the argument below in your 
# copy of this script.

cmake \
  -D CMAKE_BUILD_TYPE=RELEASE \
  -D CMAKE_C_COMPILER=cc \
  -D CMAKE_CXX_COMPILER=CC \
  -D CMAKE_Fortran_COMPILER=ftn \
  -D HOMME_PROJID=cli106 \
  -D BUILD_HOMME_SWIM=ON \
  -D BUILD_HOMME_PRIM=ON \
  -D BUILD_HOMME_SWEQX=OFF \
  -D BUILD_HOMME_PREQX=OFF \
  -D BUILD_HOMME_SWDGX=OFF \
  -D BUILD_HOMME_PRIMDGX=OFF \
  -D HOMME_BASELINE_DIR=/lustre/atlas/scratch/4ue/cli106/trunk/reg_test \
  -D ENABLE_HORIZ_OPENMP=ON \
  -D ENABLE_OPENMP=ON \
  -D REFSOLN=ON \
  -D SWIM_PLEV=1 \
  -D SWIM_NP=4 \
  -D SWIM_NC=4 \
  -D PRIM_PLEV=26 \
  -D PRIM_NP=4 \
  -D PRIM_NC=4 \
  -D CMAKE_Fortran_FLAGS="-ffree-line-length-none" \
  -D PIO_FILESYSTEM_HINTS="lustre" \
  -D CMAKE_INSTALL_PREFIX=/$MEMBERWORK/cli106/test \
  -D NETCDF_DIR=${NETCDF_DIR} \
  -D PNETCDF_DIR=${PNETCDF_DIR} \
  $HOMME_ROOT
