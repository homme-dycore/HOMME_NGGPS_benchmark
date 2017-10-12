&ctl_nl
NThreads      = @TPR@
nthreads_accel      = @TRA@
partmethod    = 4
topology      = "cube"
test_case     = "asp_baroclinic"
rotate_grid = 0
ne            = @NELEM@
qsize         = 25
tracer_advection_formulation=1
tstep_type    = 5
compute_mean_flux=1  
ndays         = 1
statefreq     = 60
accumfreq     = -1
accumstart    = 200
accumstop     = 1200
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
tstep         = 360
rsplit=5
qsplit=1
psurf_vis=0  
integration   = "explicit"
smooth = 0
nu            = 2e16
nu_s          = -1        ! use same value as nu
nu_q          = 2e16    
nu_p = 0
nu_div = -1
npdg=0
limiter_option = 4
energy_fixer = -1
hypervis_order = 2
u_perturb = 1
vert_remap_q_alg = 1
disable_diagnostics = .true.
/
&solver_nl
precon_method = "identity"
maxits        = 500
tol           = 1.e-9
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 0
filter_mu     = 0.04D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "../vcoord/camm-26.fbin.littleendian"
vfile_int     = "../vcoord/cami-26.fbin.littleendian"
/

&prof_inparm
profile_outpe_num = 100
profile_single_file	= .true.
/


&analysis_nl
! to compare with EUL ref solution:
! interp_nlat = 512
! interp_nlon = 1024
interp_gridtype = 2
 
output_timeunits = 1,1
output_frequency = 1,1
output_start_time = 1,1
output_end_time = 30000,30000
output_varnames1 = ps,zeta,T,u,v
output_varnames2 = Q
io_stride=8
output_type = 'netcdf' 
output_dir = './'
/
