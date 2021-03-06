&prof_inparm
profile_detail_limit		= 2
/

&ctl_nl
partmethod    = 4
topology      = "cube"
test_case     = "baroclinic"
ne            = 9
qsize         = 0
ndays         = 12
statefreq     = 1440
restartfreq   = -100
restartfile   = "./R0001"
runtype       = 0
tstep         = 60
integration   = 'explicit'
tracer_advection_formulation = 0
smooth        = 0.05
nu            = 7e5
nu_s          = 7e5
hypervis_order = 1
hypervis_subcycle = 1
/
&solver_nl
precon_method = "identity"
maxits        = 50
tol           = 1.e-7
/
&filter_nl
filter_type   = "taylor"
transfer_type = "bv"
filter_freq   = 0
filter_mu     = 0.04
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&vert_nl
vform         = "ccm"
vfile_mid     = "vcoord/sabm-20.fbin"
vfile_int     = "vcoord/sabi-20.fbin"
/
&analysis_nl
 output_timeunits=1,
 output_frequency=12,
 output_varnames1='T','zeta'
 interp_type = 0          ! native high order
! interp_type = 1         ! bilinear
 output_type='netcdf'
 io_stride = 8
/

