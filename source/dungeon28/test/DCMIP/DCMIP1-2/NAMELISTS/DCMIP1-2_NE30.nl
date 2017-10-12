&ctl_nl
  NThreads          = 1
  vert_num_threads  = 1
  partmethod        = 4
  topology          = "cube"
  test_case         = "DCMIP1-2"
  ne                = 30
  qsize             = 4
  ndays             = 1
  statefreq         = 48
  restartfreq       = -100
  restartfile       = "./R0001"
  runtype           = 0
  tstep             = 300
  tstep_type        = 1
  qsplit            = 1
  rsplit            = 3
  integration       = 'explicit'
  smooth            = 0
  nu                = 1E15
  nu_p              = 0
  nu_s              = 0
  nu_q              = 0
  limiter_option    = 8
  hypervis_order    = 2
  hypervis_subcycle = 1
  prescribed_wind   = 1
  energy_fixer      = -1
  vert_remap_q_alg  = 1
/
&filter_nl
  filter_type       = "taylor"
  transfer_type     = "bv"
  filter_freq       = 0
  filter_mu         = 0.04
  p_bv              = 12.0D0
  s_bv              = .80
  wght_fm           = 0.10D0
  kcut_fm           = 2
/
&vert_nl
  vform             = "ccm"
  vfile_mid         = "./12k_top-64m.ascii"
  vfile_int         = "./12k_top-64i.ascii"
/
&analysis_nl
 output_dir         = "./movies/"     ! destination dir for netcdf file
 output_timeunits   = 2,              ! 1=days, 2=hours, 0=timesteps
 output_frequency   = 12,             ! interval between outputs
 output_varnames1   = 'Q2','geo'      ! variables to write to file
 output_type        ='netcdf'         ! netcdf or pnetcdf
 num_io_procs       = 16
/
&prof_inparm
  profile_outpe_num   = 512
  profile_single_file	= .true.
/