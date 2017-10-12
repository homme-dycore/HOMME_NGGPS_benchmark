#!/bin/bash

if [ "$HOMME_ROOT" == "" ] ; then 
  echo "Must specify location of HOMME source environment variable HOMME_ROOT" 
  return -1
fi 

export  Pnetcdf_DIR=/sw/xk6/p-netcdf/1.3.1/cle4.1_pgi12.8.0

cmake \
-DCMAKE_BUILD_TYPE=RELEASE \
  -DCMAKE_C_COMPILER=cc \
  -DCMAKE_CXX_COMPILER=CC \
  -DCMAKE_Fortran_COMPILER=ftn \
  -DCMAKE_Fortran_FLAGS="-ta=nvidia -Mcuda=cc20,ptxinfo" \
  -DPIO_FILESYSTEM_HINTS="lustre" \
  -DCMAKE_EXE_LINKER_FLAGS="-Wl,--defsym,main=MAIN_" \
  -DCMAKE_CXX_FLAGS_RELEASE:STRING="-fast -O3 -DNDEBUG" \
  -DCMAKE_C_FLAGS_RELEASE:STRING="-fast -O3 -DNDEBUG" \
  -DCMAKE_Fortran_FLAGS_RELEASE:STRING="-fast -O3" \
  -DHOMME_ARCH=Linux \
  -DNUM_PLEV=1 \
  -DBUILD_HOMME_SWEQX=ON \
  -DBUILD_HOMME_SWIM=ON \
  -DCMAKE_INSTALL_PREFIX=/tmp/work/$USER/homme-pgi-gpu \
  $HOMME_ROOT



#  -DCMAKE_Fortran_FLAGS="-ta=nvidia -Mcuda=cuda4.0,cc20,ptxinfo" \
