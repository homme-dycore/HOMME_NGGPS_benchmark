#!/bin/csh

##### Before running:
# 
# module unload netcdf-hdf5parallel/4.3.0 
# module swap PrgEnv-pgi PrgEnv-gnu
# module load cmake/2.8.6 subversion python netcdf-hdf5parallel/4.3.0
# module load p-netcdf boost/1.49.0
# export XTPE_LINK_TYPE='dynamic' 
# NOTE: trilinos link below has piro, but is in a directory for cli054 users

# do either
# module load cray-trilinos
#setenv Trilinos_DIR ${CRAY_TRILINOS_PREFIX_DIR}
# or set Trilinos_DIR to your own trilinos install
setenv HOMME_ROOT /ccs/home/$USER/trunk
#setenv HOMME_PROJID cli064
setenv Trilinos_DIR /ccs/proj/cli062/worley/Trilinos/titan-gnu-albany-ci-nophal/install

# also, set environment variable HOMME_ROOT to point to 
# your HOMME repo, or change the argument below in your 
# copy of this script.

cmake \
  -D CMAKE_BUILD_TYPE=RELEASE \
  -D CMAKE_C_COMPILER=cc \
  -D CMAKE_CXX_COMPILER=CC \
  -D CMAKE_Fortran_COMPILER=ftn \
  -D HOMME_PROJID=cli064 \
  -D BUILD_HOMME_SWIM=ON \
  -D BUILD_HOMME_PRIM=ON \
  -D BUILD_HOMME_SWEQX=OFF \
  -D BUILD_HOMME_PREQX=OFF \
  -D BUILD_HOMME_SWDGX=OFF \
  -D BUILD_HOMME_PRIMDGX=OFF \
  -D REFSOLN=ON \
  -D SWIM_PLEV=1 \
  -D SWIM_NP=4 \
  -D SWIM_NC=4 \
  -D PRIM_PLEV=26 \
  -D PRIM_NP=4 \
  -D PRIM_NC=4 \
  -D CMAKE_Fortran_FLAGS="-ffree-line-length-none" \
  -D PIO_FILESYSTEM_HINTS="lustre" \
  -D CMAKE_INSTALL_PREFIX=/tmp/work/$USER/test \
  -D NETCDF_DIR=${NETCDF_DIR} \
  -D PNETCDF_DIR=${PNETCDF_DIR} \
  -D HDF5_DIR=${HDF5_DIR} \
  -D CURL_INCLUDE_DIR=/usr/bin/curl \
  $HOMME_ROOT
