!=======================================================!
! 	1 day = 1 * 24 * 3600 = 86400 sec		!
! 	nmax  = ndays * 86400 / tstep 			!
! 	12 days at 30.0 stepsize: nmax= 34560		!
!=======================================================!
&ctl_nl
NThreads      = 1
partmethod    = 4
topology      = "cube"
test_case     = "swtc5"
ne            = 6
ndays         = 15
statefreq     = 1000
tasknum       = 0
restartfreq   = -1
restartfile   = "./restart/R000000050"
runtype       = 0
tstep         = 30.0
integration   = "explicit"
smooth        = 0.05
nu            = 7.0e5
/
&solver_nl
precon_method = "block_jacobi"
maxits        = 100
tol           = 1.e-12
/
&filter_nl
transfer_type = "bv"
filter_type   = "taylor"
filter_freq   = 1
filter_mu     = 0.05D0
p_bv          = 12.0D0
s_bv          = .666666666666666666D0
wght_fm       = 0.10D0
kcut_fm       = 2
/
&analysis_nl
!=======================================================!
!  currently up to 5 streams are allowed		!
!  output_stream_count=1				!
!							!
!  timunits: 0= steps, 1=days, 2=hours			!
!  output_timeunits=1,2 				!
!  output_start_time=0,1176				!			
!  output_end_time=-1,-1				!
!  output_frequency=1,1 				!
!  output_dir ="./movies/"				!
!							!
!  allowed variables: 'ps   ','geop ','u    ','v    ',	!
!                     'latp ','lonp ','latv ','lonv ',	!
!                     'elem ','Time ' 			!
!							!
!  output_varnames1-5					!
!=======================================================!
output_start_time= 0
output_end_time  =-1
output_frequency = 1
output_timeunits = 1
output_varnames1 = 'geop','u','v'
/
&dg_nl
riemanntype= 0
alphatype= 0
alpha_dg = 0.0D0
/









