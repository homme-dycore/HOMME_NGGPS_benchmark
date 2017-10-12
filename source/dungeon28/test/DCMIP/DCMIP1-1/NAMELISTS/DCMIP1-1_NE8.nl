#DCMIP Test 1-1, NE=8 (3.75 degree) configuration

&ctl_nl
  NThreads          = 1
  vert_num_threads  = 1
  partmethod        = 4
  topology          = "cube"
  test_case         = "DCMIP1-1"
  ne                = 8
  qsize             = 4
  ndays             = 12
  rotate_grid       = 0
  statefreq         = 36
  restartfreq       = 43200
  restartfile       = "./R0001"
  runtype           = 0
  tstep             = 1200
  tstep_type        = 1
  qsplit            = 1
  rsplit            = 3
  integration       = "explicit"
  smooth            = 0.00        ! disabled
  nu                = 6E16
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
  filter_mu         = 0.04D0
  filter_freq_advection = 0
  filter_mu_advection   = 0.00
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
  output_dir       = "./movies/"                    ! destination dir for netcdf file
  output_timeunits = 1                              ! 1=days, 2=hours, 0=timesteps
  output_frequency = 1                              ! interval between outputs
  output_varnames1 = 'Q','Q3','geo'                 ! tracers Q2 and Q4 are suppressed for this test
  output_type      ='netcdf'                        ! netcdf or pnetcdf
  num_io_procs     = 16
/
&prof_inparm
  profile_outpe_num   = 512
  profile_single_file	= .true.
/
