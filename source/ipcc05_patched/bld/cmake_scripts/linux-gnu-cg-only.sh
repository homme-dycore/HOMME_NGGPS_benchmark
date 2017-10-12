#!/bin/bash

if [ "$HOMME_ROOT" == "" ] ; then 
  echo "Must specify location of HOMME source environment variable HOMME_ROOT" 
  return -1
fi 

DATE=`date +%Y-%m-%d`

cmake \
  -D CMAKE_BUILD_TYPE=RELEASE \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_Fortran_COMPILER=mpif90 \
  -D CMAKE_Fortran_FLAGS="-ffree-form -ffree-line-length-none" \
  -D NUM_PLEV=1 \
  -D SWIM_NP=8 \
  -D SWIM_NC=8 \
  -D CG_ONLY=ON \
  -D REFSOLN=ON \
  -D CMAKE_INSTALL_PREFIX=/home/$USER/homme-trunk-$DATE \
  -D NETCDF_DIR=${NETCDF_DIR} \
  -D PNETCDF_DIR=${PNETCDF_DIR} \
  -D HDF5_DIR=${HDF5_DIR} \
  -D HOMME_FIND_BLASLAPACK=ON \
  $HOMME_ROOT

