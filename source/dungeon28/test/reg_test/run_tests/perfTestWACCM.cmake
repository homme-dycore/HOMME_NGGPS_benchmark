###############################################################
# RK + PIO_INTERP 
###############################################################
#
# Spectral Element -- 9 days of ASP baroclinic test
# (Jablonowski and Williamson test + 134 tracers)
# NE=15, dt=150, nu=1e16, filter_freq=0, NV=4, PLEV=70
# (explicit RK with subcycling)
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME perfTestWACCM)
# The type of run (preqx,sweqx,swdgx,etc.)
SET(TEST_TYPE preqx)
# The specifically compiled executable that this test uses
SET(EXEC_NAME perfTestWACCM)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}.nl)
SET(VCOORD_FILES ${HOMME_ROOT}/test/vcoord/*70*)
SET(REFSOLN_FILES ${HOMME_ROOT}/test/reg_test/ref_sol/T340ref.nc)

SET(EXTRA_NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}-ne30.nl
                         ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}-ne60.nl
                         ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}-ne120.nl)

SET(NC_OUTPUT_FILES 
  camBaroMoist-asp_baroclinic1.nc 
  camBaroMoist-asp_baroclinic2.nc)


