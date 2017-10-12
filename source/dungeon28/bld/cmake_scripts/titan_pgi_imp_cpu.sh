
#!/bin/ksh


source /opt/modules/default/init/ksh

export CRAY_CPU_TARGET=istanbul
export CRAYPE_LINK_TYPE='dynamic'

source ../env_mach_specific.openacc.new
module list 

#  set environment
export HOMME_ROOT=/ccs/home/4ue/trunk
export HOMME_PROJID=cli106ms
export Trilinos_DIR=/lustre/atlas1/cli900/world-shared/cesm/software/Trilinos/Trilinos-11.12.1_gptl/titan-pgi-acme-ci-nophal_15.10.lustre/install
#export Trilinos_DIR=/ccs/home/2ya/trilinos_install

cmake \
  -C $HOMME_ROOT/cmake/machineFiles/titan.cmake  \
  -D CMAKE_Fortran_COMPILER=ftn \
  -D CMAKE_C_COMPILER=cc  \
  -D CMAKE_CXX_COMPILER=CC \
  -D OPT_FLAGS="-O2 -Minline -Mipa=fast"         \
  -D LDFLAGS= "-Wl,-z muldefs --allow-multiple-definition " \
  -D DEBUG_FLAGS= " " \
  -D HOMME_PROJID=cli106ms  \
  -D BUILD_HOMME_SWIM=ON  \
  -D BUILD_HOMME_PRIM=ON \
  -D BUILD_HOMME_SWEQX=OFF \
  -D BUILD_HOMME_PREQX=OFF \
  -D BUILD_HOMME_SWDGX=OFF \
  -D BUILD_HOMME_PRIMDGX=OFF \
  -D PIO_FILESYSTEM_HINTS="lustre" \
  -D REFSOLN=ON \
  -D SWEQX_NP=4 \
  -D PREQX_NP=4 \
  -D CMAKE_Fortran_FLAGS=" " \
  -D CURL_INCLUDE_DIR=/usr/bin/curl \
  -D NETCDF_DIR=${NETCDF_DIR} \
  -D PNETCDF_DIR=${PNETCDF_DIR} \
  -D HDF5_DIR=${HDF5_DIR} \
  -D ENABLE_OPENMP=FALSE \
  $HOMME_ROOT

make -j16


