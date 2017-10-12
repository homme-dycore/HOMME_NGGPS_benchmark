#!/bin/csh

##### Before running:
# 
# module load cmake subversion 
# module swap PrgEnv-pgi PrgEnv-gnu
# module load netcdf pnetcdf hdf5 
# export XTPE_LINK_TYPE='dynamic' 
# NOTE: trilinos link below should be accessable to all users

# do either
# module load cray-trilinos
#setenv Trilinos_DIR ${CRAY_TRILINOS_PREFIX_DIR}
# or set Trilinos_DIR to your own trilinos install
setenv HOMME_ROOT /global/homes/e/$USER/trunk
setenv HOMME_PROJID 1
setenv Trilinos_DIR /project/projectdirs/ccsm1/Trilinos/trilinos-10.12.2/hopper-gnu/install 

# also, set environment variable HOMME_ROOT to point to 
# your HOMME repo, or change the argument below in your 
# copy of this script.

cmake \
  -D CMAKE_BUILD_TYPE=RELEASE \
  -D CMAKE_C_COMPILER=cc \
  -D CMAKE_CXX_COMPILER=CC \
  -D CMAKE_Fortran_COMPILER=ftn \
  -D FI_ONLY=ON \
  -D REFSOLN=OFF \
  -D SWIM_PLEV=1 \
  -D SWIM_NP=16 \
  -D SWIM_NC=16 \
  -D PRIM_PLEV=26 \
  -D PRIM_NP=4 \
  -D PRIM_NC=4 \
  -D CMAKE_Fortran_FLAGS="-ffree-line-length-none" \
  -D PIO_FILESYSTEM_HINTS="lustre" \
  -D CMAKE_INSTALL_PREFIX=/$GSCRATCH/$USER/test \
  -D NETCDF_DIR=${NETCDF_DIR} \
  -D PNETCDF_DIR=${PNETCDF_DIR} \
  -D HDF5_DIR=${HDF5_DIR} \
  $HOMME_ROOT

# --prefix="$HOME/seacism"


