&ctl_nl
horz_num_threads = 1
vert_num_threads = 1
tracer_num_threads = 1
partmethod    = 4
topology      = "cube"
test_case     = "jw_baroclinic"
u_perturb     = 1
rotate_grid   = 0
ne            = 256
qsize         = 10
!ndays          = 1
!  nmax = 2160  24h  assuming tstep=40s    
!  nmax = 1440  16h    
!  nmax = 1080  12h
!  nmax = 720   8h
!  nmax = 180   2h
nmax = 1440
statefreq     = 90
restartfreq   = 43200
restartfile   = "./R0001"
runtype       = 0
tstep         = 40
rsplit        = 3
qsplit        = 1
tstep_type    = 5
energy_fixer  = -1
integration   = "explicit"
smooth        = 0
nu            = 1e12
nu_div        = 1e12
nu_p          = 1e12
nu_q          = 1e12
nu_s          = -1
nu_top        = 1.0e5
se_ftype      = 0
limiter_option = 8
vert_remap_q_alg = 1
hypervis_order = 2
hypervis_subcycle=3
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
vfile_mid     = "./sabm-128.ascii"
vfile_int     = "./sabi-128.ascii"
/

&prof_inparm
profile_outpe_num = 100
profile_single_file		= .true.
/

&analysis_nl
! to compare with EUL ref solution:
 interp_nlat = 32
 interp_nlon = 64
 interp_gridtype=2
 
 output_timeunits=1,1
 output_frequency=0,0
 output_start_time=0,0
 output_end_time=30000,30000
 output_varnames1='ps','zeta','dp3d'
 output_varnames2='Q','Q2','Q3','Q4','Q5'
 io_stride=8
 output_type = 'netcdf' 
/

