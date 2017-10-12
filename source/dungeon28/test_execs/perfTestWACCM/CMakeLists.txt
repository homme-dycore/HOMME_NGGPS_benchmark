SET_SOURCE_FILES_PROPERTIES(${HOMME_SOURCE_DIR}/src/share/viscosity_mod.F90 PROPERTIES COMPILE_FLAGS "-O3 -check bounds")
# Set the qsize_d in dimensions_mod.F90
SET(QSIZE_D 135)
# Set the variables for this test executable
#                                  NP  NC PLEV USE_PIO  WITH_ENERGY
createTestExec(perfTestWACCM preqx 4   4   70   FALSE        TRUE)
