!-----------------------------------------------------------------------------------!
! MODULE FVM_MOD----------------------------------------------------------CE-for FVM!
! fvm_MOD File for the fvm project in HOMME                                         !
! Author: Christoph Erath                                                           !
! Date: 25.January 2011                                                             !
! MAIN module to run fvm on HOMME                                                   !
! 14.November 2011: reorganisation done                                             !
!-----------------------------------------------------------------------------------!
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module fvm_bench_mod
  use kinds, only : real_kind, int_kind, longdouble_kind,iulog
  use edge_mod, only : freeghostbuffertr, ghostVpack, ghostVunpack, &
                       edgeVpack, edgeVunpack, freeedgebuffer 
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, np, ne, nc, nhe, nlev, ntrac
  use time_mod, only : timelevel_t
  use element_mod, only : element_t, timelevels
  use hybrid_mod, only : hybrid_t
  use perf_mod, only : t_startf, t_stopf, t_barrierf ! EXTERNAL
  use perf_utils, only : t_detail_low, t_detail_medium, t_detail_high, t_detail_max ! EXTERNAL


contains


subroutine cslam_run_bench(elem,fvm,red,hybrid,nets,nete,tl)
  ! ---------------------------------------------------------------------------------  
  use fvm_bsp_mod, only: fvm_bsp, get_boomerang_velocities_gll, get_solidbody_velocities_gll
  ! ---------------------------------------------------------------------------------  
  use fvm_control_volume_mod, only: fvm_struct
  ! ---------------------------------------------------------------------------------
  use fvm_line_integrals_flux_mod, only: cslam_runflux
  use fvm_mod, only: fill_halo_fvm, cellghostbuf, edgeveloc
  use fvm_line_integrals_mod, only : cslam_runairdensity, fvm_mcgregor,fvm_mcgregordss, &
         fvm_rkdss
  use fvm_consistent_se_cslam, only : run_consistent_se_cslam
  ! ---------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------
  use derivative_mod, only : derivative_t, derivative_stag_t, derivinit, deriv_print
  ! ---------------------------------------------------------------------------------
  use reduction_mod, only : ReductionBuffer_ordered_1d_t
  ! ---------------------------------------------------------------------------------
  use time_mod, only : tstep, nmax, time_at, timelevel_update, timelevel_init
  ! ---------------------------------------------------------------------------------
  use coordinate_systems_mod, only : spherical_polar_t,spherical_to_cart, &
                                     cart2cubedspherexy
  ! ---------------------------------------------------------------------------------
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest
  ! ---------------------------------------------------------------------------------
  use parallel_mod, only: global_shared_buf, global_shared_sum
  ! ---------------------------------------------------------------------------------
  use global_norms_mod, only: wrap_repro_sum
  use physical_constants    , only : dd_pi
  ! ---------------------------------------------------------------------------------
  use reduction_mod, only : parallelmax, parallelmin
  ! ---------------------------------------------------------------------------------  
  use physical_constants, only : g
  ! -----------------------------------------------
  use control_mod, only : TRACERTRANSPORT_LAGRANGIAN_FVM, TRACERTRANSPORT_FLUXFORM_FVM
  use control_mod, only : TRACERTRANSPORT_CONSISTENT_SE_FVM
  use control_mod, only : tracer_transport_type, qsplit
  use parallel_mod, only: abortmp
  use fvm_control_volume_mod     , only: n0_fvm, np1_fvm
  use common_io_mod, only: output_frequency
#ifdef PIO_INTERP
     use interp_movie_mod, only : interp_movie_init, interp_movie_output, interp_movie_finish
#else
     use shal_movie_mod, only : shal_movie_init, shal_movie_output, shal_movie_finish
#endif

  implicit none
  type (element_t), intent(inout)                :: elem(:)
  type (fvm_struct), intent(inout)             :: fvm(:)
  type (ReductionBuffer_ordered_1d_t),intent(in)    :: red   ! reduction buffer         (shared)
  type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
  integer, intent(in)                         :: nets  ! starting thread element number (private)
  integer, intent(in)                         :: nete  ! ending thread element number   (private)

  real (kind=real_kind)                       :: massstart, mass, maxc, maxcstart,minc, mincstart, tmp, tmpref  
  real (kind=real_kind)                       :: tmp1(nets:nete), tmp2(nets:nete), tmpt
  real (kind=real_kind)                       :: l1,l2, lmax 
  
  integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h
  type (TimeLevel_t)                          :: tl              ! time level struct
  type (derivative_t)                         :: deriv           ! derivative struct
  type (spherical_polar_t)                    :: tmpsphereincart   
 
  character (len=99)                          :: filename
  
  real (kind=real_kind)   , dimension(10*(nc+2*nhe)*(nc+2*nhe),6)  :: weights_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_eul_index_all
  integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2)  :: weights_lgr_index_all
  integer (kind=int_kind)                                          :: jall
    
  real (kind=real_kind), dimension(5,1-nhe:nc+nhe,1-nhe:nc+nhe)      :: recons  
  real (kind=real_kind), dimension(nelemd,1-nhe:nc+nhe,1-nhe:nc+nhe) :: area    
  real (kind=real_kind)                                              :: xtmp
  real (kind=longdouble_kind)                                        :: fvm_nodes(nc+1)
  
  real (kind=real_kind), dimension(np,np,2)    :: vstar, vhat
  real (kind=real_kind)                        :: maxcflx, maxcfly  
  !
  ! fvm diagnostics
  !
  real (kind=real_kind) :: csum(ntrac+1),relative_mass_change(ntrac+1),cmin(ntrac),cmax(ntrac)
  real (kind=real_kind) :: psc_mass, psc_min, psc_max,dp_fvm_mass, dp_fvm_min, dp_fvm_max 

  
  integer                                     :: ierr
  
  integer  chooselev   !for test reason the output

 !-----------------------------------------------------------------------------------!  
 chooselev=1
 
  if(hybrid%masterthread) then 
    print *,"!-----------------------------------------------------------------------!"
    print *,"!  Test CASE for fvm, Christoph Erath                                 !" 
    print *,"!-----------------------------------------------------------------------!" 
  endif
  qsplit=1
!  tracer_transport_type = TRACERTRANSPORT_LAGRANGIAN_FVM


  ! Initialize derivative structure
  ! fvm nodes are equally spaced in alpha/beta
  ! HOMME with equ-angular gnomonic projection maps alpha/beta space
  ! to the reference element via simple scale + translation
  ! thus, fvm nodes in reference element [-1,1] are a tensor product of
  ! array 'fvm_nodes(:)' computed below:
  xtmp=nc 
  do i=1,nc+1
    fvm_nodes(i)= 2*(i-1)/xtmp - 1
  end do
  call derivinit(deriv,fvm_corners=fvm_nodes)
!-----------------------------------------------------------------------------------! 
!  call TimeLevel_Qdp(tl, qsplit, n0_fvm, np1_fvm)

  do ie=nets,nete
     fvm(ie)%dp_ref=1.0D0
     !
     ! Initialize fields
     !
     call fvm_bsp(fvm(ie),tl)
     do j=1,nc
        do i=1,nc
           fvm(ie)%psc(i,j) = fvm(ie)%dp_fvm(i,j,1,n0_fvm)
        end do
     end do
     !
     ! reset the new unknown
     !
     fvm(ie)%c     (:,:,:,:,np1_fvm)=-10000.0D0 !dbg
     fvm(ie)%dp_fvm(:,:,:  ,np1_fvm)=-10000.0D0 !dbg
  end do

  
  !first exchange of the initial values
  call fill_halo_fvm(elem,fvm,hybrid,nets,nete,n0_fvm)
!-----------------------------------------------------------------------------------!     

!Initialize Output via geopotential (should be changed, separate output for fvm
!write first time step to IO 
#ifdef PIO_INTERP
  call interp_movie_init(elem,hybrid,nets,nete,tl=tl)    
  call interp_movie_output(elem,tl, hybrid, 0D0, nets, nete,fvm)
#else
    call shal_movie_init(elem,hybrid,fvm)
    call shal_movie_output(elem,tl, hybrid, 0D0, nets, nete,deriv,fvm)
#endif 
!
!-----------------------------------------------------------------------------------!
!
  if(hybrid%masterthread) then 
    print *
    print *,"Arrival grid created , interpolation points calculated, initialization done. " 
    print *
  endif
  tmp=0
  
  call t_barrierf('fvm time loop', hybrid%par%comm)
  call t_startf('fvm', t_detail_high)


  !
  ! initialize for fvm diagnostics
  !
  if (ntrac>0) then
     do itr=1,ntrac
        do ie=nets,nete
           tmp1(ie) = MINVAL(fvm(ie)%c(1:nc,1:nc,:,itr,n0_fvm)) 
        enddo
        cmin(itr) = ParallelMin(tmp1,hybrid)
        do ie=nets,nete
           tmp1(ie) = MAXVAL(fvm(ie)%c(1:nc,1:nc,:,itr,n0_fvm))
        enddo
        cmax(itr) = ParallelMax(tmp1,hybrid)
        !
        ! compute total tracer mass
        !
        global_shared_buf(:,1) = 0.0D0
        do k=1,nlev
           do ie=nets,nete
              global_shared_buf(ie,1) = global_shared_buf(ie,1)+&
                   SUM(fvm(ie)%c(1:nc,1:nc,k,itr,n0_fvm)*&
                   fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)*&
                   fvm(ie)%area_sphere(1:nc,1:nc))
           end do
        enddo
        call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
        csum(itr) = global_shared_sum(1)
     enddo
     !
     ! psC diagnostics (PSC = surface pressure implied by fvm)
     !
     do ie=nets,nete
        tmp1(ie) = MINVAL(fvm(ie)%psc(1:nc,1:nc))
     enddo
     psc_min = ParallelMin(tmp1,hybrid)
     do ie=nets,nete
        tmp1(ie) = MAXVAL(fvm(ie)%psc(1:nc,1:nc))
     enddo
     psc_max = ParallelMax(tmp1,hybrid)
     global_shared_buf(:,1) = 0.0D0
     do ie=nets,nete
        global_shared_buf(ie,1) = SUM(fvm(ie)%psc(1:nc,1:nc)*fvm(ie)%area_sphere(1:nc,1:nc))
     enddo
     call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
     psc_mass = global_shared_sum(1)     
     !
     ! dp_fvm
     !
     do ie=nets,nete
        tmp1(ie) = MINVAL(fvm(ie)%dp_fvm(1:nc,1:nc,:,n0_fvm))
     enddo
     dp_fvm_min = ParallelMin(tmp1,hybrid)
     do ie=nets,nete
        tmp1(ie) = MAXVAL(fvm(ie)%dp_fvm(1:nc,1:nc,:,n0_fvm))
     enddo
     dp_fvm_max = ParallelMax(tmp1,hybrid)     
     global_shared_buf(:,1) = 0.0D0
     do k=1,nlev
        do ie=nets,nete
           global_shared_buf(ie,1) = global_shared_buf(ie,1)+&
                SUM(fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)*fvm(ie)%area_sphere(1:nc,1:nc))
        end do
     enddo
     call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
     dp_fvm_mass = global_shared_sum(1)
     
     !
     ! save mass of initial condition
     !
     do ie=nets,nete
        fvm(ie)%mass(1:ntrac) = csum(1:ntrac)
        fvm(ie)%mass(ntrac+1)  = psc_mass        
     end do

     if(hybrid%masterthread) then
        !
        ! fvm diagnostics
        !
        write(iulog,'(A37)') "-------------------------------------"
        write(iulog,'(A37)') "fvm diagnostics for initial condition"
        write(iulog,'(A37)') "-------------------------------------"
        do itr=1,ntrac
           write(iulog,'(A29,I1,4(E23.15))')&
                "#c, min, max, ave = ",itr,cmin(itr), cmax(itr), csum(itr)/(4.0D0*DD_PI)
        enddo
        write(iulog,'(A37,3(E23.15))')&
             "   min(dp_), max(dp_), ave(dp_) =  ",dp_fvm_min, dp_fvm_max, dp_fvm_mass/(4.0D0*DD_PI)
        write(iulog,'(A37,3(E23.15))')&
             "   min(psC), max(psC), ave(psC) =  ",psc_min, psc_max, psC_mass/(4.0D0*DD_PI)
        write(iulog,'(A36)') "                                   "
        
     end if
  end if

    do ie=nets,nete
       do j=1,nc
          do i=1,nc
             fvm(ie)%cstart(1:nc,1:nc,1:ntrac)=fvm(ie)%c(1:nc,1:nc,1,1:ntrac,n0_fvm)
          enddo
       enddo
    enddo


  
  !BEGIN TIME LOOP, start at 0, calculate then next step
  DO WHILE(tl%nstep< nmax)
     ! start old mcgregor----------------------
     !     do ie=nets,nete
     !       do k=1,nlev
     !         vstar = get_boomerang_velocities_gll(elem(ie), time_at(tl%nstep+1))
     !         vhat= (get_boomerang_velocities_gll(elem(ie), time_at(tl%nstep)) + vstar) / 2.0D0
     !         ! calculate high order approximation
     !         call fvm_mcgregor(elem(ie), deriv, tstep, vhat, vstar, 3)
     !      
     !         ! apply DSS to make vstar C0
     !         elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%spheremp(:,:)*vstar(:,:,1) 
     !         elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%spheremp(:,:)*vstar(:,:,2) 
     !       enddo
     !       call edgeVpack(edgeveloc,elem(ie)%derived%vstar(:,:,1,:),nlev,0,ie)
     !       call edgeVpack(edgeveloc,elem(ie)%derived%vstar(:,:,2,:),nlev,nlev,ie)
     !     enddo 
     !     call bndry_exchangeV(hybrid,edgeveloc)
     !     do ie=nets,nete
     !        call edgeVunpack(edgeveloc,elem(ie)%derived%vstar(:,:,1,:),nlev,0,ie)
     !        call edgeVunpack(edgeveloc,elem(ie)%derived%vstar(:,:,2,:),nlev,nlev,ie)
     !        do k=1, nlev  
     !          elem(ie)%derived%vstar(:,:,1,k)=elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
     !          elem(ie)%derived%vstar(:,:,2,k)=elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
     !        end do
     !     end do
     !end old mcgegor-----------------------
     ! ! start mcgregordss
     do ie=nets,nete
        do k=1,nlev
           elem(ie)%derived%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie), time_at(tl%nstep+1))
           fvm(ie)%vn0(:,:,:,k)=get_boomerang_velocities_gll(elem(ie),time_at(tl%nstep))
           !         elem(ie)%derived%vstar(:,:,:,k)=get_solidbody_velocities_gll(elem(ie), time_at(tl%nstep+1))
           !         fvm(ie)%vn0(:,:,:,k)=get_solidbody_velocities_gll(elem(ie),time_at(tl%nstep))
        end do
     end do
     !     call fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, tstep, 3)
     call fvm_rkdss(elem,fvm,nets,nete, hybrid, deriv, tstep, 3)
     
     !     tmpt=(time_at(tl%nstep+1)-time_at(tl%nstep))/3
     !     do ie=nets,nete
     ! !       do i=1,4
     !         do k=1,nlev
     !           fvm(ie)%vstar(:,:,:,k)=get_boomerang_velocities_gll(elem(ie),time_at(tl%nstep+1)-tmpt*(i-1))
     ! !         elem(ie)%derived%vstar(:,:,:,k)=get_solidbody_velocities_gll(elem(ie), time_at(tl%nstep+1))
     ! !         fvm(ie)%vn0(:,:,:,k)=get_solidbody_velocities_gll(elem(ie),time_at(tl%nstep))
     !         end do
     ! !       end do
     !     end do    
     
     
     ! ! end mcgregordss   
        
     if (tracer_transport_type == TRACERTRANSPORT_CONSISTENT_SE_FVM) then
        if(mod(tl%nstep,1)==0.and.hybrid%masterthread) write(iulog,*) "running consistent se-cslam"
        call run_consistent_se_cslam(elem,fvm,hybrid,deriv,tstep,tl,nets,nete,0.0D0)        
     else if (tracer_transport_type == TRACERTRANSPORT_FLUXFORM_FVM) then
        call cslam_runflux      (elem,fvm,hybrid,deriv,tstep,tl,nets,nete,0.0D0,.FALSE.)
        if(mod(tl%nstep,1)==0.and.hybrid%masterthread) write(iulog,*) "running ff-cslam"
     else if (tracer_transport_type == TRACERTRANSPORT_LAGRANGIAN_FVM) then
        call cslam_runairdensity(elem,fvm,hybrid,deriv,tstep,tl,nets,nete,0.0D0) !run regular CSLAM
        if(mod(tl%nstep,1)==0.and.hybrid%masterthread) write(iulog,*) "running cslam"
     else
        call abortmp('Bad tracer_transport_type in fvm_bench')
     end if

     call TimeLevel_update(tl,"forward")
!    call TimeLevel_Qdp(tl, qsplit, n0_fvm, np1_fvm)
     write(iulog,*) "tl%nstep",tl%nstep
     if (tl%nstep==0.or.mod(tl%nstep,1)==0) then  
        !
        ! fvm diagnostics
        !
        if (ntrac>0) then
           do itr=1,ntrac
              do ie=nets,nete
                 tmp1(ie) = MINVAL(fvm(ie)%c(1:nc,1:nc,:,itr,n0_fvm)) 
              enddo
              cmin(itr) = ParallelMin(tmp1,hybrid)
              do ie=nets,nete
                 tmp1(ie) = MAXVAL(fvm(ie)%c(1:nc,1:nc,:,itr,n0_fvm))
              enddo
              cmax(itr) = ParallelMax(tmp1,hybrid)
              !
              ! compute total tracer mass
              !
              global_shared_buf(:,1) = 0.0D0
              do k=1,nlev
                 do ie=nets,nete
                    global_shared_buf(ie,1) = global_shared_buf(ie,1)+&
                         SUM(fvm(ie)%c(1:nc,1:nc,k,itr,n0_fvm)*&
                         fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)*&
                         fvm(ie)%area_sphere(1:nc,1:nc))
                 end do
              enddo
              call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
              csum(itr) = global_shared_sum(1)
           enddo
           do itr=1,ntrac              
              ie=nets
              if (ABS(fvm(ie)%mass(itr))<1.0E-14) then
                 relative_mass_change(itr) = csum(itr) - fvm(ie)%mass(itr)
              else
                 relative_mass_change(itr) = (csum(itr) - fvm(ie)%mass(itr))/fvm(ie)%mass(itr)
              end if
           end do
           !
           ! psC diagnostics (PSC = surface pressure implied by fvm)
           !
           do ie=nets,nete
              tmp1(ie) = MINVAL(fvm(ie)%psc(1:nc,1:nc))
           enddo
           psc_min = ParallelMin(tmp1,hybrid)
           do ie=nets,nete
              tmp1(ie) = MAXVAL(fvm(ie)%psc(1:nc,1:nc))
           enddo
           psc_max = ParallelMax(tmp1,hybrid)
           do ie=nets,nete
              global_shared_buf(ie,1) = SUM(fvm(ie)%psc(1:nc,1:nc)*fvm(ie)%area_sphere(1:nc,1:nc))
           enddo
           call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
           psc_mass = global_shared_sum(1)
           !
           ! dp_fvm
           !
           do ie=nets,nete
              tmp1(ie) = MINVAL(fvm(ie)%dp_fvm(1:nc,1:nc,:,n0_fvm))
           enddo
           dp_fvm_min = ParallelMin(tmp1,hybrid)
           do ie=nets,nete
              tmp1(ie) = MAXVAL(fvm(ie)%dp_fvm(1:nc,1:nc,:,n0_fvm))
           enddo
           dp_fvm_max = ParallelMax(tmp1,hybrid)
           
           global_shared_buf(:,1) = 0.0D0
           do k=1,nlev
              do ie=nets,nete
                 global_shared_buf(ie,1) = global_shared_buf(ie,1)+&
                      SUM(fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)*fvm(ie)%area_sphere(1:nc,1:nc))
              end do
           enddo
           call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
           dp_fvm_mass = global_shared_sum(1)
           
           !
           !
           !
           if (tl%nstep==0) then
              do ie=nets,nete
                 fvm(ie)%mass(:)       = 0.0D0
                 fvm(ie)%mass(1:ntrac) = csum(1:ntrac)
              end do
           end if
           
           !    end if
           !
           !
           maxcflx = parallelmax(fvm(:)%maxcfl(1,chooselev),hybrid)
           maxcfly = parallelmax(fvm(:)%maxcfl(2,chooselev),hybrid)

           if(hybrid%masterthread) then
              !
              ! fvm diagnostics
              !
              write(iulog,'(A36)') "-----------------------------------"
              write(iulog,'(A36)') "fvm diagnostics                    "
              write(iulog,'(A36)') "-----------------------------------"

              write(iulog,*) 'time=', time_at(tl%nstep), 'timeatmax',Time_at(nmax)
              write(iulog,*) 'STEP',tl%nstep,'MAXSTEP',nmax,'t0', n0_fvm, 't1', np1_fvm
              write(iulog,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly 
              print *

              do itr=1,ntrac
                 write(iulog,'(A29,I1,4(E23.15))')&
                      "#c, min, max, ave, change = ",itr,cmin(itr), cmax(itr), csum(itr)/(4.0D0*DD_PI),relative_mass_change(itr)
              enddo
              write(iulog,'(A37,3(E23.15))')&
                   "   min(dp_), max(dp_), ave(dp_) =  ",dp_fvm_min, dp_fvm_max, dp_fvm_mass/(4.0D0*DD_PI)
                   
                   write(iulog,'(A37,4(E23.15))')&
                   "   min(psC), max(psC), ave(psC) =  ",psc_min, psc_max, psC_mass/(4.0D0*DD_PI),&
                   (psc_mass-fvm(nets)%mass(ntrac+1))/fvm(nets)%mass(ntrac+1)
              write(iulog,'(A36)') "                                   "
                 
           end if

        endif
     endif
    !-----------------------------------------------------------------------------------!  
    
#ifdef PIO_INTERP
     call interp_movie_output(elem, tl, hybrid, 0D0, nets, nete,fvm)
#else     
     call shal_movie_output(elem, tl, hybrid, 0D0, nets, nete,deriv,fvm)
#endif
     !-----------------------------------------------------------------------------------!  
  END DO
!------------END TIME LOOP-------------END TIME LOOP--------------END TIME LOOP-----!
!-----------------------------------------------------------------------------------! 
  call t_stopf('fvm', t_detail_high)


  call freeghostbuffertr(cellghostbuf)
  call freeedgebuffer(edgeveloc)
#ifdef PIO_INTERP
    call interp_movie_finish
#else
    call shal_movie_finish
#endif
!-----------------------------------------------------------------------------------!  

    !SUMMARY
    if(hybrid%masterthread) then 
       print *
       print *,"!-----------------------------------------------------------------------!"
       print *,"!  Test CASE for FVM, Christoph Erath                                   !" 
       print *,"!-----------------------------------------------------------------------!"
       print *, "Summary" 
       write(iulog,*) 'number of elements', 6*ne*ne*nc*nc
       
       print *
       
       write(iulog,*) "CFL: maxcflx=", maxcflx, "maxcfly=", maxcfly 
       write(iulog,*) "l1 = ", l1, "l2 = ", l2, "lmax = ", lmax
       write(iulog,*) "ne*nc = ", ne*nc, "timestep = ", tstep
    endif

! Error analysis/ complicated, but for a first try o.k.
    do itr=1,ntrac
       do ie=nets,nete
          tmp=0.0D0
          tmpref=0.0D0
          global_shared_buf(ie,:)=0.0D0
          do j=1,nc
             do i=1,nc
                global_shared_buf(ie,1)=global_shared_buf(ie,1)+ &
                     fvm(ie)%area_sphere(i,j)*abs(fvm(ie)%c(i,j,chooselev,itr,n0_fvm)-fvm(ie)%cstart(i,j,itr))
                global_shared_buf(ie,2)=global_shared_buf(ie,2)+fvm(ie)%area_sphere(i,j)*abs(fvm(ie)%cstart(i,j,itr))

                global_shared_buf(ie,3)=global_shared_buf(ie,3)+ &
                     fvm(ie)%area_sphere(i,j)*(fvm(ie)%c(i,j,chooselev,itr,n0_fvm)-fvm(ie)%cstart(i,j,itr))* &
                     (fvm(ie)%c(i,j,chooselev,itr,n0_fvm)-fvm(ie)%cstart(i,j,itr))
                global_shared_buf(ie,4)=global_shared_buf(ie,4)+fvm(ie)%area_sphere(i,j)*&
                     (fvm(ie)%cstart(i,j,itr))*(fvm(ie)%cstart(i,j,itr))
          
                tmp=max(tmp,abs(fvm(ie)%c(i,j,chooselev,itr,n0_fvm)-fvm(ie)%cstart(i,j,itr)))
                tmpref=max(tmpref,abs(fvm(ie)%cstart(i,j,itr)))
             end do
          end do
          tmp1(ie)=tmp
          tmp2(ie)=tmpref
       end do
       call wrap_repro_sum(nvars=4, comm=hybrid%par%comm)
       l1=global_shared_sum(1)/global_shared_sum(2)
       l2=sqrt(global_shared_sum(3)/global_shared_sum(4))
   
       lmax = parallelmax(tmp1,hybrid)/parallelmax(tmp2,hybrid)


       !SUMMARY
       if(hybrid%masterthread) then 
          write(iulog,*) "Error norms for tracer ",itr
          write(iulog,*) "l1 = ", l1, "l2 = ", l2, "lmax = ", lmax
          write(iulog,*) "  "
       endif
    end do
  0817 format("*****ELEMENT ",I6,2x,I6,2x,I1)
end subroutine cslam_run_bench


end module fvm_bench_mod
