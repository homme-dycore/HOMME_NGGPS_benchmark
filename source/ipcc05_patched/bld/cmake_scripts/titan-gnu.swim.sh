#!/bin/csh

##### Before running:
# module unload netcdf-hdf5parallel/4.2.0 
# module swap PrgEnv-pgi PrgEnv-gnu
# module swap gcc/4.8.1 gcc/4.7.2
# module load cmake/2.8.6 subversion python netcdf-hdf5parallel/4.2.0
# module load p-netcdf boost/1.49.0
# export XTPE_LINK_TYPE='dynamic' 
# NOTE: trilinos link below has piro, but is in a directory for cli054 users

# do either
# module load cray-trilinos
# setenv Trilinos_DIR ${CRAY_TRILINOS_PREFIX_DIR}
# or set Trilinos_DIR to your own trilinos install

# also, edit environment variable HOMME_ROOT below to point to 
# your HOMME repo, or change the argument below in your 
# copy of this script. Also edit HOMME_PROJID below.

# EDIT the following three variables according to directions above!
setenv HOMME_ROOT /your_homme_root_dir
setenv HOMME_PROJID your_projID
setenv Trilinos_DIR /your_trilinos_dir

cmake \
  -D CMAKE_BUILD_TYPE=RELEASE \
  -D CMAKE_C_COMPILER=cc \
  -D CMAKE_CXX_COMPILER=CC \
  -D CMAKE_Fortran_COMPILER=ftn \
  -D BUILD_HOMME_SWIM=ON \
  -D CMAKE_Fortran_FLAGS="-ffree-line-length-none" \
  -D PIO_FILESYSTEM_HINTS="lustre" \
  -D CMAKE_INSTALL_PREFIX=/tmp/work/$USER/test \
  -D NETCDF_DIR=${NETCDF_DIR} \
  -D PNETCDF_DIR=${PNETCDF_DIR} \
  -D HDF5_DIR=${HDF5_DIR} \
  $HOMME_ROOT
