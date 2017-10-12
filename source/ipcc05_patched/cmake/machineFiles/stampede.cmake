# CMake initial cache file for stampede
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{TACC_NETCDF_DIR} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{TACC_HDF5_DIR} CACHE FILEPATH "")
SET (SZIP_DIR $ENV{TACC_HDF5_DIR} CACHE FILEPATH "")
SET (PNETCDF_DIR $ENV{TACC_PNETCDF_DIR} CACHE FILEPATH "")
SET (USE_MPIEXEC ibrun CACHE FILEPATH "")
