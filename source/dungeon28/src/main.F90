#ifdef HAVE_CONFIG_H
#include "config.h"
#endif


program main
  use init_mod, only : init
  ! -----------------------------------------------
  use sweq_mod, only : sweq, sweq_rk
  ! -----------------------------------------------
  !      use kinds
  ! -----------------------------------------------
  use parallel_mod, only : parallel_t, initmpi, syncmp, haltmp
  ! -----------------------------------------------
  use thread_mod, only : initomp, max_num_threads, omp_get_thread_num, omp_set_num_threads
  use hybrid_mod, only : init_loop_ranges
  ! -----------------------------------------------
  !      use time_mod
  ! -----------------------------------------------
  use dimensions_mod, only : nelemd
  ! -----------------------------------------------
  ! use domain_mod, only : domain1d_t, decompose
  ! -----------------------------------------------
  use element_mod, only : element_t
  ! -----------------------------------------------
  use fvm_control_volume_mod, only : fvm_struct
  ! -----------------------------------------------
  !      use state_mod
  ! -----------------------------------------------
  use edgetype_mod, only : EdgeBuffer_t
  ! -----------------------------------------------
  use reduction_mod, only : ReductionBuffer_ordered_1d_t
  ! -----------------------------------------------
  use perf_mod, only : t_initf, t_prf, t_finalizef, t_startf, t_stopf ! EXTERNAL
  use perf_utils, only : t_detail_minimal, t_detail_low, t_detail_medium, t_detail_high, t_detail_max ! EXTERNAL
  ! -----------------------------------------------
  use control_mod, only : integration

  implicit none
  type (element_t), pointer :: elem(:)
  type (fvm_struct), pointer  :: fvm(:)
  
  type (EdgeBuffer_t)  :: edge1            ! 1 component edge buffer (1, 3d scalar field)
  type (EdgeBuffer_t)  :: edge2            ! 2 component edge buffer (1, 3d vector field)
  type (EdgeBuffer_t)  :: edge3            ! 3 component edge buffer (1, 3d vector + 1 3d scalar field)
  type (ReductionBuffer_ordered_1d_t)  :: red    ! reduction buffer for cg
  type (parallel_t)    :: par              ! parallel structure for distributed memory programming

  integer nets,nete
  integer ithr
  integer ierr
  integer :: maxthreads

  ! =====================================================
  ! Begin executable code set distributed memory world...
  ! =====================================================

  par=initmpi()
  call initomp()

  call t_initf(LogPrint=par%masterproc, &
       Mpicom=par%comm, MasterTask=par%masterproc, &
       maxthreads=max_num_threads)
  call t_startf('Total',t_detail_minimal)
  
  call init(elem,edge1,edge2,edge3,red,par,fvm)
  ! =====================================================
  ! Allocate state variables
  ! =====================================================

  if(par%masterproc) print *,"allocating state variables..."
  !JMD allocate(state(nelemd))

  ! =====================================
  ! Set number of threads...
  ! =====================================
  call init_loop_ranges(nelemd)
  
  nets = 1
  nete = nelemd 
    

  if(par%masterproc) print *,"Main:max_num_threads=",max_num_threads

  ! =====================================
  !  Sync-up to make sure timing is clean
  ! =====================================

  call syncmp(par)

  ! =====================================
  ! Begin threaded region...
  ! =====================================
#if (defined HORIZ_OPENMP)
  !$OMP PARALLEL DEFAULT(SHARED), PRIVATE(ithr,nets,nete)
#endif

  !
  ! ================================================
  ! Initialize thread decomposition
  ! ================================================
  !
  write(6,9) par%rank,ithr,nets,nete 
9 format("process: ",i2,1x,"thread: ",i2,1x,"element limits: ",i4," - ",i4)

  if(integration == "runge_kutta")then
     call sweq_rk(elem,edge1,edge2,edge3,red,par,ithr,nets,nete)
  else
     call sweq(elem,fvm,edge1,edge2,edge3,red,par,ithr,nets,nete)
  endif

#if (defined HORIZ_OPENMP)
  !$OMP END PARALLEL
#endif
  ! ================================================
  ! End distributed memory region
  ! ================================================
  call t_stopf('Total', t_detail_minimal)
  call t_prf('HommeSWTime',par%comm)
  call t_finalizef()
  call haltmp("exiting program...")
  deallocate(elem)
end program main
