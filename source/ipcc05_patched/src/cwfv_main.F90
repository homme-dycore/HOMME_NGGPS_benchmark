!----------------------------------------------------------------------------!
! Main File for the CWFV project in HOMME                                   !
! Author: Christoph Erath                                                    !
! Date: 13.January 2012                                                      !
!----------------------------------------------------------------------------!

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


program cwfv_main
  ! -----------------------------------------------
  use cwfv_mod, only: fvm_init1, fvm_init2, cwfv_struct
  use cwfv_bench_mod, only: cwfv_run_bench
  ! -----------------------------------------------
  ! -----------------------------------------------
  use dimensions_mod, only : nelemd, ntrac
  ! -----------------------------------------------
  use domain_mod, only : domain1d_t, decompose
  ! -----------------------------------------------
  use edge_mod, only : EdgeBuffer_t
  ! -----------------------------------------------
  use element_mod, only : element_t
  ! -----------------------------------------------
  use fvm_init_mod, only : fvm_init
  ! -----------------------------------------------
  use time_mod, only : timelevel_t, timelevel_init
  ! -----------------------------------------------
  use parallel_mod, only : parallel_t, initmp, syncmp, haltmp
  ! -----------------------------------------------
  use reduction_mod, only : ReductionBuffer_ordered_1d_t
  ! -----------------------------------------------
  use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads
  use hybrid_mod, only : hybrid_create, hybrid_t 
  ! ------EXTERNAL----------------
  use perf_mod, only : t_initf, t_prf, t_finalizef, t_startf, t_stopf ! _EXTERNAL
  ! -----------------------------------------------

  implicit none
  type (element_t), pointer :: elem(:)
  type (cwfv_struct), pointer                :: fvm(:)
  type (TimeLevel_t)                          :: tl              ! time level struct
  
  type (EdgeBuffer_t)  :: edgecell,edgepoints                 ! 1 component edge buffer (1, 3d scalar field)
  
  type (ReductionBuffer_ordered_1d_t)  :: red    ! reduction buffer for cg
  type (parallel_t)    :: par              ! parallel structure for distributed memory programming
  type (hybrid_t) :: hybrid
  type (domain1d_t), allocatable :: dom_mt(:)
  

  character (len=20)                          :: numproc_char
  character (len=20)                          :: numtrac_char

  integer                                     :: nets,nete
  integer                                     :: ithr
  integer                                     :: ierr


  ! =====================================================
  ! Begin executable code set distributed memory world...
  ! =====================================================

  par=initmp()
  ! initialisation for parallel
  call t_initf('input.nl',LogPrint=par%masterproc, &
       Mpicom=par%comm, MasterTask=par%masterproc)
  call t_startf('Total')
  call fvm_init(elem,red,par)
  ! =====================================================
  ! Allocate state variables
  ! =====================================================

  if(par%masterproc) print *,"allocating state variables..."
  !JMD allocate(state(nelemd))
  if(par%masterproc) print *,"number of elements=",SIZE(elem)

  ! =====================================
  ! Set number of threads...
  ! =====================================

  if(par%masterproc) print *,"Main:NThreads=",NThreads


  ! =====================================
  !  Sync-up to make sure timing is clean
  ! =====================================

  call syncmp(par)
  call fvm_init1(par)


  ! initialize time management
  call timeLevel_init(tl)



  call omp_set_num_threads(NThreads)
  allocate(dom_mt(0:NThreads-1))
  do ithr=0,NThreads-1
     dom_mt(ithr)=decompose(1,nelemd,NThreads,ithr)
  end do

  ! =====================================
  ! Begin threaded region...
  ! =====================================
#if (defined HORIZ_OPENMP)
  !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ithr,nets,nete)
#endif
  ithr=omp_get_thread_num()
  hybrid = hybrid_create(par,ithr,NThreads)  
  nets=dom_mt(ithr)%start
  nete=dom_mt(ithr)%end
  print *,"nets: ",nets,"nete: ",nete, "ithr: ",ithr
  !
  ! ================================================
  ! Initialize thread decomposition
  ! ================================================
  !
  write(6,9) par%rank,ithr,nets,nete 
9 format("process: ",i2,1x,"thread: ",i2,4x,"element limits: ",i4,i4)




  ! ================================================
  ! Initialize ghost cell code (should be done in threaded region, so not in fvm_init1)
  ! ================================================
    

  !-Create new CWFV Mesh
  allocate(fvm(nelemd))
  call fvm_init2(elem,fvm,hybrid,nets,nete,tl)




  call cwfv_run_bench(elem,fvm,hybrid,nets,nete,tl)

#if (defined HORIZ_OPENMP)
  !$OMP END PARALLEL
#endif
  ! ================================================
  ! End distributed memory region
  ! ================================================
  call t_stopf('Total')
  if(par%masterproc) print *,"writing timing data"
  call t_prf('HommeTimeCWFV', par%comm)
  if(par%masterproc) print *,"calling t_finalizef"
  call t_finalizef()
  call haltmp("exiting program...")
  deallocate(elem)
end program cwfv_main

