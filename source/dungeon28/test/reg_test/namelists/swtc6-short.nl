&prof_inparm
profile_detail_limit = 2
/

&ctl_nl
partmethod    = 0
topology      = "cube"
test_case     = "swtc6"
ne            = 2
ndays         = 1
statefreq     = 1440
restartfreq   = -100
restartfile   = "./R0001"
runtype       = 0
tstep         = 6.
smooth        = 0.00
integration   = "explicit"
nu = 0e16
nu_s = 0e16
hypervis_order = 2
hypervis_subcycle = 1
/
&solver_nl
precon_method = "block_jacobi"
maxits        = 100
tol           = 1.e-12
debug_level   = 0
/
&filter_nl
transfer_type = "bv"
filter_type   = "taylor"
filter_freq   = 0
filter_mu     = 0.04D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
kcut_fm       = 2
wght_fm       = 0.10D0
/
&analysis_nl
output_timeunits = 2
output_frequency = 6
output_varnames1 = 'u', 'v', 'ps'
interp_type      = 0          ! native high order
output_type      = 'netcdf'
io_stride        = 8
/
