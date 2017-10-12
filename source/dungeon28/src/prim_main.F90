#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

program prim_main
  ! -----------------------------------------------
#ifdef _PRIM
  use prim_driver_mod, only : prim_init1, prim_init2, prim_run, prim_finalize,&
      leapfrog_bootstrap, prim_run_subcycle, prim_commtest
  use hybvcoord_mod, only : hvcoord_t, hvcoord_init
#else
! use these for dg3d

#endif

  ! -----------------------------------------------
  use parallel_mod, only : parallel_t, initmpi, syncmp, haltmp, abortmp, iam
  ! -----------------------------------------------
  use hybrid_mod, only : hybrid_t, PrintHybrid, init_loop_ranges, get_loop_ranges, config_thread_region
  ! -----------------------------------------------
  use thread_mod, only : max_num_threads, vert_num_threads, tracer_num_threads, horz_num_threads, &
                         omp_get_thread_num, omp_set_num_threads, omp_get_nested, &
                         omp_get_num_threads, omp_get_max_threads, initomp
  ! ----------------------------------------------- 
  use time_mod, only : tstep, nendstep, timelevel_t, TimeLevel_init, time_at, dt_phys
  ! -----------------------------------------------
  use dimensions_mod, only : nelemd, qsize, ntrac
  ! -----------------------------------------------
  use control_mod, only : restartfreq, vfile_mid, vfile_int, runtype, integration, statefreq, tstep_type, rsplit, qsplit, commtest
  ! -----------------------------------------------
  use element_mod, only : element_t
  !-----------------------------------------------
  use chemistry, only : compute_chemistry_FQ
  use fvm_control_volume_mod, only : fvm_struct
  use fvm_control_volume_mod, only : n0_fvm
  use fvm_mod, only : fill_halo_fvm
  use spelt_mod, only : spelt_struct
  ! -----------------------------------------------
  use common_io_mod, only:  output_dir
  ! -----------------------------------------------
  use umjs_baroclinic_wave_mod, only : prim_printstate_par_terminator
#ifdef _REFSOLN
  use prim_state_mod, only : prim_printstate_par
  ! -----------------------------------------------
#endif


#ifdef PIO_INTERP
  use interp_movie_mod, only : interp_movie_output, interp_movie_finish, interp_movie_init
  use interpolate_driver_mod, only : interpolate_driver
#else
  use prim_movie_mod, only : prim_movie_output, prim_movie_finish,prim_movie_init
#endif
  use common_movie_mod, only : nextoutputstep
  use perf_mod, only : t_initf, t_prf, t_finalizef, t_startf, t_stopf ! EXTERNAL
  use perf_utils, only : t_detail_minimal, t_detail_low, t_detail_medium, t_detail_high, t_detail_max ! EXTERNAL
  !-----------------
  use restart_io_mod , only : restartheader_t, writerestart
  !-----------------
  use hybrid_mod, only : config_thread_region

  implicit none
  type (element_t), pointer  :: elem(:)
#if defined(_SPELT)
    type (spelt_struct), pointer   :: fvm(:)
#else
     type (fvm_struct), pointer   :: fvm(:)    
#endif

  type (hybrid_t)       :: hybrid ! parallel structure for shared memory/distributed memory
  type (parallel_t)                    :: par  ! parallel structure for distributed memory programming
  type (RestartHeader_t)                  :: RestartHeader
  type (TimeLevel_t)    :: tl             ! Main time level struct
  type (hvcoord_t)     :: hvcoord         ! hybrid vertical coordinate struct

  real*8 timeit, et, st, tmp, cly_mass_init(2)
  integer ithr
  integer :: nets, nete
  integer ierr
  integer nstep
  integer nstep_end
  integer maxthreads

  character (len=20)                          :: numproc_char
  character (len=20)                          :: numtrac_char
  
  logical :: dir_e ! boolean existence of directory where output netcdf goes
  
  ! =====================================================
  ! Begin executable code set distributed memory world...
  ! =====================================================
  par=initmpi()
  call initomp()

  cly_mass_init = 0 !only used in terminator test

  ! =====================================
  ! Set number of threads...
  ! =====================================

  call t_initf(LogPrint=par%masterproc, &
       Mpicom=par%comm, MasterTask=par%masterproc, &
       maxthreads=max_num_threads)
  call t_startf('Total',t_detail_minimal)
  call t_startf('prim_init1', t_detail_high)
  call prim_init1(elem,  fvm, par,tl)
  call t_stopf('prim_init1',t_detail_high)


  ! =====================================
  ! Begin threaded region so each thread can print info
  ! =====================================
  ! =====================================
  ! Begin threaded region so each thread can print info
  ! =====================================
  call init_loop_ranges(nelemd)

!  print *,'_main: ID: ', iam, hybrid%ithr,' sum(rspheremp): ',sum(elem(1)%rspheremp)

  ! ================================================
  ! Initialize thread decomposition
  ! Note: The OMP Critical is required for threading since the Fortran 
  !   standard prohibits multiple I/O operations on the same unit.
  ! ================================================
!  if (par%rank<100) then 
!     write(6,9) par%rank,ithr,nets,nete 
!  endif
!9 format("process: ",i2,1x,"thread: ",i2,1x,"element limits: ",i5," - ",i5)
! back to single threaded
  hybrid = config_thread_region(par,'serial')

  ! ==================================
  ! Initialize the vertical coordinate  (cam initializes hvcoord externally)
  ! ==================================
  hvcoord = hvcoord_init(vfile_mid, vfile_int, .true., hybrid%masterthread, ierr)
  if (ierr /= 0) then
     call haltmp("error in hvcoord_init")
  end if

#ifdef PIO_INTERP
  if(runtype<0) then
     ! Interpolate a netcdf file from one grid to another
     call interpolate_driver(elem, hybrid)
     call haltmp('interpolation complete')
  end if
#endif

  if(par%masterproc) print *,"Primitive Equation Initialization..."
  call t_startf('prim_init2', t_detail_high)
  !$OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete)
  hybrid = config_thread_region(par,'horizontal')
  call get_loop_ranges(hybrid,ibeg=nets,iend=nete)

  call prim_init2(elem, fvm, hybrid, nets, nete, tl, hvcoord)

  !$OMP END PARALLEL
  call t_stopf('prim_init2', t_detail_high)

  hybrid = config_thread_region(par,'serial')
  
  ! Here we get sure the directory specified
  ! in the input namelist file in the 
  ! variable 'output_dir' does exist.
  ! this avoids a abort deep within the PIO 
  ! library (SIGABRT:signal 6) which in most
  ! architectures produces a core dump.
  if (hybrid%masterthread) then 
     open(unit=447,file=trim(output_dir) // "/output_dir_test",iostat=ierr)
     if ( ierr==0 ) then
        print *,'Directory ',trim(output_dir), ' does exist: initialing IO'
        close(447)
     else
        print *,'Error creating file in directory ',trim(output_dir)
        call haltmp("Please be sure the directory exist or specify 'output_dir' in the namelist.")
     end if
  endif
#if 0
  this ALWAYS fails on lustre filesystems.  replaced with the check above
  inquire( file=output_dir, exist=dir_e )
  if ( dir_e ) then
     if(hybrid%masterthread) print *,'Directory ',output_dir, ' does exist: initialing IO'
  else
     if(hybrid%masterthread) print *,'Directory ',output_dir, ' does not exist: stopping'
     call haltmp("Please get sure the directory exist or specify one via output_dir in the namelist file.")
  end if
#endif

#ifdef PIO_INTERP
  ! initialize history files.  filename constructed with restart time
  ! so we have to do this after ReadRestart in prim_init2 above
  call interp_movie_init( elem, hybrid, 1, nelemd, hvcoord, tl )
#else
  call prim_movie_init( elem, fvm, par, hvcoord, tl )
#endif


  ! output initial state for NEW runs (not restarts or branch runs)
  if (runtype == 0 ) then
#ifdef PIO_INTERP
     call interp_movie_output(elem, tl, hybrid, 0d0, 1, nelemd,fvm=fvm, hvcoord=hvcoord)
#else
     call prim_movie_output(elem, tl, hvcoord, hybrid, 1,nelemd, fvm)
#endif
  endif


  ! advance_si not yet upgraded to be self-starting.  use leapfrog bootstrap procedure:
  if(integration == 'semi_imp') then
     if (runtype /= 1 ) then
        if(hybrid%masterthread) print *,"Leapfrog bootstrap initialization..."
        call leapfrog_bootstrap(elem, hybrid,1,nelemd,tstep,tl,hvcoord)
     endif
  endif
  
  if(par%masterproc) then 
    if (dt_phys>0) then
      print *,"chemistry time-step is ",dt_phys

      if (tstep*rsplit.NE.dt_phys) then
         write(*,*) "dt_phys must be remap time-step",dt_phys
         dt_phys = tstep*rsplit
         write(*,*) "changing dt_phys to            ",dt_phys
!         stop
      end if
    else
      print *,"no idealized physics/chemistry",dt_phys
    end if
  end if

  if(par%masterproc) print *,"Entering main timestepping loop"
  call t_startf('prim_main_loop', t_detail_minimal)
  do while(tl%nstep < nEndStep)

     !$OMP PARALLEL NUM_THREADS(horz_num_threads), DEFAULT(SHARED), PRIVATE(hybrid,nets,nete)
     hybrid = config_thread_region(par,'horizontal')
     call get_loop_ranges(hybrid,ibeg=nets,iend=nete)
     if(tl%nstep==0) then
       if (par%rank<100) then 
           write(6,9) par%rank,hybrid%ithr,nets,nete
       endif
     endif
 9 format("process: ",i2,1x,"thread: ",i2,1x,"element limits: ",i5," - ",i5)

     nstep = nextoutputstep(tl)

     !
     ! intial condition terminator chemistry diagnostics
     !
     if (dt_phys>0.and.tl%nstep==0) call prim_printstate_par_terminator(elem, tl,hybrid,hvcoord,nets,nete, fvm,cly_mass_init,.true.) !terminator

!JMD     print *,'prim_main: nets:nete ',nets,nete
     do while(tl%nstep<nstep)
        call t_startf('prim_run', t_detail_low)
        !
        ! physics forcing from HOMME
        !
!JMD         call PrintHybrid(hybrid,'prim_main: before compute_chemistry_FQ')
        if (dt_phys>0) then
           tmp = time_at(tl%nstep)/dt_phys
           if (tmp==NINT(tmp)) then !only call physics every dt_phys seconds
              call compute_chemistry_FQ(elem,hybrid,hvcoord,nets,nete, fvm, tl%n0,n0_fvm)
           end if
        end if
!JMD         call PrintHybrid(hybrid,'prim_main: after compute_chemistry_FQ')

        nstep_end = tl%nstep + qsplit
        if (rsplit>0) then
           nstep_end = tl%nstep + qsplit*rsplit
        endif

        if (tstep_type>0) then  ! forward in time subcycled methods
           if(commtest) then  
               call prim_commtest(hybrid,nets,nete,tl)
           else
               call prim_run_subcycle(elem, fvm, hybrid,nets,nete, tstep, tl, hvcoord,1)
           endif
        else  ! leapfrog
           call prim_run(elem, hybrid,nets,nete, tstep, tl, hvcoord, "leapfrog")
        endif
        !
        ! compute diagnostics for terminator chemistry
        !
        if (MODULO(nstep_end,statefreq)==0 .or. nstep_end==tl%nstep0) then
           if (dt_phys>0) call prim_printstate_par_terminator(elem, tl,hybrid,hvcoord,nets,nete,  fvm,cly_mass_init,.false.) !terminator
        endif

        call t_stopf('prim_run', t_detail_low)
     end do
     !$OMP END PARALLEL 
     hybrid = config_thread_region(par,'serial')
#ifdef PIO_INTERP
     if (ntrac>0) call fill_halo_fvm(elem,fvm,hybrid,nets,nete,n0_fvm)
     call interp_movie_output(elem, tl, hybrid, 0d0, 1, nelemd,fvm=fvm, hvcoord=hvcoord)
#else
     call prim_movie_output(elem, tl, hvcoord, hybrid, 1,nelemd, fvm)
#endif
     

#ifdef _REFSOLN
     call prim_printstate_par(elem, tl,hybrid,hvcoord,nets,nete, par)
#endif 
     
     ! ============================================================
     ! Write restart files if required 
     ! ============================================================
     if((restartfreq > 0) .and. (MODULO(tl%nstep,restartfreq) ==0)) then 
        call WriteRestart(elem, hybrid%ithr,1,nelemd,tl)
     endif
  end do
  call t_stopf('prim_main_loop', t_detail_minimal)

  if(par%masterproc) print *,"Finished main timestepping loop",tl%nstep
  call prim_finalize(hybrid)
  if(par%masterproc) print *,"closing history files"
#ifdef PIO_INTERP
  call interp_movie_finish
#else
  call prim_movie_finish
#endif


  call t_stopf('Total', t_detail_minimal)
  if(par%masterproc) print *,"writing timing data"
!   write(numproc_char,*) par%nprocs
!   write(numtrac_char,*) ntrac
!   call system('mkdir -p '//'time/'//trim(adjustl(numproc_char))//'-'//trim(adjustl(numtrac_char))) 
!   call t_prf('time/HommeFVMTime-'//trim(adjustl(numproc_char))//'-'//trim(adjustl(numtrac_char)),par%comm)
  call t_prf('HommeTime', par%comm)
  if(par%masterproc) print *,"calling t_finalizef"
  call t_finalizef()
  call haltmp("exiting program...")
end program prim_main

