###############################################################
#
# Discontinuous Galerkin -- swtc5
# NE=6, dt=30, nu=7e5, limiter=0, filter_freq=1, NP=6
#
###############################################################

# The name of this test (should be the basename of this file)
SET(TEST_NAME swtc5-dg)
# The type of run (preqx,sweqx,swdgx,etc.)
SET(TEST_TYPE swdgx)
# The specifically compiled executable that this test uses
SET(EXEC_NAME swtc-dgA)

SET(NUM_CPUS 16)

SET(NAMELIST_FILES ${HOMME_ROOT}/test/reg_test/namelists/${TEST_NAME}.nl)
SET(NCL_FILES ${HOMME_ROOT}/test/reg_test/ncl/swtc5ref.ncl)

SET(NC_OUTPUT_FILES swtc51.nc)

