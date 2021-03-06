#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

program dg_main
!=======================================================================================================!
        use init_mod, only : init
! -----------------------------------------------
use dg_sweq_mod, only : sweq_dg
! -----------------------------------------------
!      use kinds
! -----------------------------------------------
      use parallel_mod, only : parallel_t, initmp, syncmp, haltmp
! -----------------------------------------------
      use thread_mod, only : nthreads, omp_get_thread_num, omp_set_num_threads
! -----------------------------------------------
!      use time_mod
! -----------------------------------------------
      use dimensions_mod, only : nelemd, np, ne, nlev
! -----------------------------------------------
      use domain_mod, only : domain1d_t, decompose
! -----------------------------------------------
      use element_mod, only : element_t
! -----------------------------------------------
!      use state_mod
! -----------------------------------------------
      use edge_mod, only : EdgeBuffer_t
! -----------------------------------------------
      use reduction_mod, only : ReductionBuffer_ordered_1d_t
! -----------------------------------------------
!=======================================================================================================!
      implicit none
      type (element_t), pointer :: elem(:)
      type (EdgeBuffer_t)  :: edge1            ! 1 component edge buffer (1, 3d scalar field)
      type (EdgeBuffer_t)  :: edge2            ! 2 component edge buffer (1, 3d vector field)
      type (EdgeBuffer_t)  :: edge3            ! 3 component edge buffer (1, 3d vector + 1 3d scalar field)
      type (ReductionBuffer_ordered_1d_t)  :: red    ! reduction buffer for cg
      type (parallel_t)    :: par              ! parallel structure for distributed memory programming
      type (domain1d_t), allocatable :: dom_mt(:)

      integer nets,nete
      integer ithr
      integer ierr
!=======================================================================================================!

!=======================================================================================================!
!	Begin executable code set distributed memory world...						!
!=======================================================================================================!      
      par=initmp()

      call init(elem, edge1,edge2,edge3,red,par)
  
!      open (unit=80,file='cs_geometry.dat')
!      write(80,*) np,nv,ne,nlev
!      close(80)
!=======================================================================================================!
!	Allocate state variables									!
!=======================================================================================================!
      if(par%masterproc) print *,'allocating state variables...'
 
!=======================================================================================================!
!	Set number of threads...									!
!=======================================================================================================!
      if(par%masterproc) print *,'Main:NThreads=',NThreads

      call omp_set_num_threads(NThreads)
     
      allocate(dom_mt(0:NThreads-1))
      do ithr=0,NThreads-1
         dom_mt(ithr)=decompose(1,nelemd,NThreads,ithr)
      end do
!=======================================================================================================!
!	Sync-up to make sure timing is clean								!
!=======================================================================================================!
      call syncmp(par)

!=======================================================================================================!
!	Begin threaded region...									!
!=======================================================================================================!

      ithr=omp_get_thread_num()
      nets=dom_mt(ithr)%start
      nete=dom_mt(ithr)%end
!=======================================================================================================!
!	Initialize thread decomposition									!
!=======================================================================================================!
      write(6,9) par%rank,ithr,nets,nete 
  9   format('process: ',i5,1x,'thread: ',i2,1x,'element limits: ',i4," - ",i4)

      call sweq_dg(elem,edge1,edge2,edge3,red,par,ithr,nets,nete)

      if(par%masterproc) then 
         print *,'Finishing DG 3D Computation...'
      endif
!=======================================================================================================!
!	Sync-up to make sure timing is clean								!
!=======================================================================================================!
      call syncmp(par)
      deallocate(elem)

!=======================================================================================================!
!	End distributed memory region									!
!=======================================================================================================!

      call haltmp('exiting program...')
      
!=======================================================================================================!
end program dg_main
