#!/bin/bash

if [ "$HOMME_ROOT" == "" ] ; then 
  echo "Must specify location of HOMME source environment variable HOMME_ROOT" 
  return -1
fi 

DATE=`date +%Y-%m-%d`
export HOMME_ROOT=/home/4ue/trunk

cmake \
  -D CMAKE_BUILD_TYPE=DEBUG \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D CMAKE_Fortran_COMPILER=mpif90 \
  -D CMAKE_INSTALL_PREFIX=/home/$USER/homme-trunk-$DATE \
  -D CMAKE_Fortran_FLAGS="-ffree-form -ffree-line-length-none" \
  -D FI_ONLY=ON \
  -D REFSOLN=ON \
  -D SWIM_PLEV=1 \
  -D SWIM_NP=4 \
  -D SWIM_NC=4 \
  -D PRIM_PLEV=26 \
  -D PRIM_NP=4 \
  -D PRIM_NC=4 \
  -D NETCDF_DIR=${NETCDF_DIR} \
  -D PNETCDF_DIR=${PNETCDF_DIR} \
  -D HDF5_DIR=${HDF5_DIR} \
  -D HOMME_FIND_BLASLAPACK=ON \
  -D BUILD_HOMME_PRIMDGX=FALSE \
  $HOMME_ROOT

