!-----------------------------------------------------------------------------------!
!MODULE FVM_MOD-----------------------------------------------------------CE-for FVM!
! FVM_MOD File for the fvm project in HOMME                                         !
! Author: Christoph Erath                                                           !
! Date: 25.January 2011                                                             !
! MAIN module to run fvm on HOMME                                                   !
! 14.November 2011: reorganisation done                                             !
! 7.Februar 2012: cslam_run and cslam_runair                                        !
!-----------------------------------------------------------------------------------!

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module fvm_mod    
  use kinds, only : real_kind
  use edge_mod, only : initghostbufferTR, freeghostbuffertr, &
       ghostVpack, ghostVunpack,  initEdgebuffer
  use edgetype_mod, only : ghostbuffertr_t, edgebuffer_t
  use bndry_mod, only: ghost_exchangeV                     
  use dimensions_mod, only: nlev, ne, nc, nhe, nlev, ntrac, np, ntrac_d,ns, nhr, nhc
  use time_mod, only : timelevel_t
  use element_mod, only : element_t, timelevels
  use fvm_control_volume_mod, only: fvm_struct
  use hybrid_mod, only : hybrid_t
  use perf_mod, only : t_startf, t_stopf ! EXTERNAL
  use perf_utils, only : t_detail_low, t_detail_medium, t_detail_high, t_detail_max ! EXTERNAL
  
  implicit none
  private
  save
  
  type (ghostBuffertr_t)                      :: cellghostbuf
  type (EdgeBuffer_t)                         :: edgeveloc

  public :: cellghostbuf, edgeveloc, fvm_init1,fvm_init2, fvm_init_flux, fill_halo_fvm
contains
  subroutine fill_halo_fvm(elem,fvm,hybrid,nets,nete,tnp0)
    use dimensions_mod, only: irecons_tracer, irecons_air
    implicit none
    type (element_t),intent(inout)            :: elem(:)                 
    type (fvm_struct),intent(inout)           :: fvm(:)  
    type (hybrid_t),intent(in)                :: hybrid                  
    
    integer,intent(in)                        :: nets,nete,tnp0       
    integer                                   :: ie                 
    !
    ! note "call initghostbufferTR(cellghostbuf,nlev,ntrac+1,nhc,nc)" in fvm_init1.
    ! should initghostbuffer be called here?
    !
    call t_startf('FVM pack', t_detail_high)
    do ie=nets,nete
       call ghostVpack(cellghostbuf, fvm(ie)%dp_fvm(:,:,:,tnp0),nhc,nc,nlev,1,    0,elem(ie)%desc)
       call ghostVpack(cellghostbuf, fvm(ie)%c(:,:,:,:,tnp0)   ,nhc,nc,nlev,ntrac,1,elem(ie)%desc)
    end do
    call t_stopf('FVM pack', t_detail_high)
    call t_startf('FVM Communication', t_detail_medium)
    call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc,ntrac+1)
    call t_stopf('FVM Communication', t_detail_medium)
    !-----------------------------------------------------------------------------------!                        
    call t_startf('FVM Unpack', t_detail_high)
    do ie=nets,nete
       call ghostVunpack(cellghostbuf, fvm(ie)%dp_fvm(:,:,:,tnp0), nhc, nc,nlev,1    ,0,elem(ie)%desc)
       call ghostVunpack(cellghostbuf, fvm(ie)%c(:,:,:,:,tnp0),    nhc, nc,nlev,ntrac,1,elem(ie)%desc)
    enddo
    call t_stopf('FVM Unpack', t_detail_high)
!    call freeghostbuffertr(cellghostbuf)    
  end subroutine fill_halo_fvm


  
  ! initialize global buffers shared by all threads
  subroutine fvm_init1(par,elem)
    use parallel_mod, only : parallel_t, haltmp
    use control_mod, only : tracer_transport_type, tracer_grid_type, rsplit
    use control_mod, only : TRACERTRANSPORT_LAGRANGIAN_FVM, TRACERTRANSPORT_FLUXFORM_FVM, TRACER_GRIDTYPE_FVM
    use control_mod, only : TRACERTRANSPORT_CONSISTENT_SE_FVM
    use fvm_control_volume_mod, only: n0_fvm, np1_fvm, fvm_supercycling
    use dimensions_mod, only: qsize, irecons_air, irecons_tracer
    type (parallel_t) :: par
    type (element_t)  :: elem(:)

    !
    ! initialize fvm time-levels
    !
    n0_fvm  = 1
    np1_fvm = 2
    !
    if (par%masterproc) then 
       print *, "                                       "
       print *, "|-------------------------------------|"
       print *, "| Tracer transport scheme information |"
       print *, "|-------------------------------------|"
       print *, "                                       "
    end if
    if (tracer_transport_type == TRACERTRANSPORT_LAGRANGIAN_FVM) then
       if (par%masterproc) then 
          print *, "Running Lagrangian CSLAM, Lauritzen et al., (2010), J. Comput. Phys."
          print *, "Possibly with `Enforcement Of Consistency (EOC)', Erath et al., (2013), Mon. Wea. Rev."
          print *, "(EOC is hardcoded in fvm_lineintegrals_mod with logical EOC"
          print *, "CSLAM = Conservative Semi-LAgrangian Multi-tracer scheme"
          print *, "  "
       end if
    else if (tracer_transport_type == TRACERTRANSPORT_FLUXFORM_FVM) then
       if (par%masterproc) then 
          print *, "Running Flux-form CSLAM, Harris et al. (2011), J. Comput. Phys."
          print *, "CSLAM = Conservative Semi-LAgrangian Multi-tracer scheme"
          print *, "Lauritzen et al., (2010), J. Comput. Phys."
          print *, "  "
       end if
    else if (tracer_transport_type == TRACERTRANSPORT_CONSISTENT_SE_FVM) then
       if (par%masterproc) then 
          print *, "Running consistent SE-CSLAM, Lauritzen et al. (2015, in prep)."
          print *, "Air flux prescribed by SE (Taylor, Overfelt, Ullrich)"
          print *, "SE = Spectral Element"
          print *, "CSLAM = Conservative Semi-LAgrangian Multi-tracer scheme"
          print *, "Lauritzen et al., (2010), J. Comput. Phys."
          print *, "  "
       end if
    else
       call haltmp("going into fvm_init1 with inconsistent tracer_transport_type")
    end if

    if (par%masterproc) print *, "fvm resolution is nc*nc in each element: nc = ",nc

    if (tracer_grid_type.ne.TRACER_GRIDTYPE_FVM) then
       if (par%masterproc) then 
         print *, "ERROR: tracer_grid_type is not TRACER_GRIDTYPE_FVM"
         print *, "tracer_grid_type = ",tracer_grid_type
       end if
       call haltmp("going into fvm_init1 with inconsistent tracer_grid_type")
    end if

    if (nc<3) then
       if (par%masterproc) then 
          print *, "NUMBER OF CELLS ERROR for fvm: Number of cells parameter"
          print *, "parameter nc at least 3 (nc>=3), nc*nc cells per element. This is"
          print *, "needed for the cubic reconstruction, which is only implemented yet! STOP"
       endif
       call haltmp("stopping")
    end if

    if (par%masterproc) then
       print *, "  "
       if (ns==1) then
          print *, "ns==1: using no interpolation for mapping cell averages values across edges"
          print *, "Note: this is not a recommended setting - large errors at panel edges!"
       else if (ns==2) then
          print *, "ns==2: using linear interpolation for mapping cell averages values across edges"
          print *, "Note that ns=4 is default CSLAM setting used in Lauritzen et al. (2010)"
          print *, "so this option is slightly less accurate (but the stencil is smaller near panel edges!)"

       else if (ns==3) then
          print *, "ns==3: using quadratic interpolation for mapping cell averages values across edges"
          print *, "Note that ns=4 is default CSLAM setting used in Lauritzen et al. (2010)"
          print *, "so this option is slightly less accurate (but the stencil is smaller near panel edges!)"
       else if (ns==4) then
          print *, "ns==4: using cubic interpolation for mapping cell averages values across edges"
          print *, "This is default CSLAM setting used in Lauritzen et al. (2010)"
       else 
          print *, "Not a tested value for ns but it should work! You choose ns = ",ns
       end if

!       if (ns.NE.3) then
!         write(*,*) "In fvm_reconstruction_mod function matmul_w has been hard-coded for ns=3 for performance"
!         write(*,*) "Revert to general code - outcommented above"
!         call haltmp("stopping")
!       end if
    end if

    if (MOD(ns,2)==0.and.nhr+(nhe-1)+ns/2>nc+nc) then
       print *, "to run this combination of ns and nhr you need to increase nc to ",nhr+ns/2+nhe-1
       print *, "You choose (ns,nhr,nc,nhe)=",ns,nhr,nc,nhe
       call haltmp("stopping")
    end if
    if (MOD(ns,2)==1.and.nhr+(ns-1)/2+(nhe-1)>nc+nc) then
       print *, "to run this combination of ns and nhr you need to increase nc to ",nhr+(ns-1)/2+nhe-1
       print *, "You choose (ns,nhr,nc,nhe)=",ns,nhr,nc,nhe
       call haltmp("stopping")
    end if

    if (nc==3.and.ns.ne.3) then
       if (par%masterproc) then
          print *, "Recommended setting for nc=3 is ns=3 (linear interpolation in halo)"
          print *, "You choose ns=",ns
          print *, "Goto dimensions_mod to change value of ns"
          print *, "or outcomment call haltmop below (i.e. you know what you are doing!)"
       endif
       call haltmp("stopping")
    end if

    if (nc==4.and.ns.ne.4) then
       if (par%masterproc) then
          print *, "Recommended setting for nc=4 is ns=4 (cubic interpolation in halo)"
          print *, "You choose ns=",ns
          print *, "Goto dimensions_mod to change value of ns"
          print *, "or outcomment call haltmop below (i.e. you know what you are doing!)"
       endif
       call haltmp("stopping")
    end if

    if (nhe .ne. 1) then
       if (par%masterproc) then
          print *, "PARAMETER ERROR for fvm: Number of halo zone for the extended"
          print *,"element nhe has to be 1, only this is available now! STOP!"
       endif
       call haltmp("stopping")
    end if

    if (irecons_air>irecons_tracer) then
       write(*,*) "configuration irecons_air>irecons_tracer", irecons_air,irecons_tracer
       write(*,*) "not supported - ABORT"
       call haltmp("stopping")
    end if
    !


    if (ntrac>ntrac_d) then
       if (par%masterproc) print *,'ntrac,ntrac_d=',ntrac,ntrac_d
       call haltmp("PARAMETER ERROR for fvm: ntrac > ntrac_d")
    endif

    if (qsize>0) then
       if (par%masterproc) then
          print *, 'FYI: running both SE and fvm tracers!'
       end if
    end if

    if (qsize>0.and.mod(rsplit,fvm_supercycling).ne.0) then
       if (par%masterproc) then
          print *,'cannot supercycle fvm tracers with respect to se tracers'
          print *,'with this choice of rsplit =',rsplit
          print *,'rsplit must be a multiple of fvm_supercycling=',fvm_supercycling
          call haltmp("PARAMETER ERROR for fvm: mod(rsplit,fvm_supercycling<>0")
       end if
    endif


    if (par%masterproc) then 
       print *, "                                            "
       print *, "Done Tracer transport scheme information    "
       print *, "                                            "
    end if
    
    call initghostbufferTR(cellghostbuf,nlev,ntrac+1,nhc,nc)
    call initEdgebuffer(par,edgeveloc,elem,2*nlev)
  end subroutine fvm_init1
  
  
  
  ! initialization that can be done in threaded regions
  subroutine fvm_init2(elem,fvm,hybrid,nets,nete,tl)
    use fvm_control_volume_mod, only: fvm_mesh_ari
    use fvm_analytic_mod, only: computexytosphere_moments
    use bndry_mod, only: compute_ghost_corner_orientation
    use bndry_mod, only : ghost_exchangevfull
    
    !  use edge_mod, only : ghostvpackfull, ghostvunpackfull
    
    type (timelevel_t) :: tl
    type (fvm_struct) :: fvm(:)
    type (element_t) :: elem(:)
    type (hybrid_t)                             :: hybrid
    integer :: ie,nets,nete
    
    call compute_ghost_corner_orientation(hybrid,elem,nets,nete)
    ! run some tests:
    !    call test_ghost(hybrid,elem,nets,nete)
    
    do ie=nets,nete
       call fvm_mesh_ari(elem(ie),fvm(ie),tl)
       call computexytosphere_moments(fvm(ie),elem(ie)%desc)
    enddo    
  end subroutine fvm_init2

  
  subroutine fvm_init_flux(elem,fvm,hybrid,nets,nete,irecons)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
!    use dimensions_mod, only: irecons_tracer, irecons_air
    implicit none
    type (element_t) ,intent(inout)  :: elem(:)
    type (fvm_struct),intent(inout)  :: fvm(:) 
    type (hybrid_t)  ,intent(in)     :: hybrid                      
    integer          ,intent(in)     :: nets,nete,irecons
    !
    type (ghostBuffertr_t)  :: cellghostbuf_tmp
    integer                 :: ie, ixy, ivertex, i, j,istart,itot,ishft,imin,imax
    integer, dimension(1:4) :: ipermute = (/1,2,3,4/)
    integer, dimension(2,4) :: unit_vec
    integer                 :: rot90_matrix(2,2), iside, count, m, n
    real (kind=real_kind)   :: cartx,carty    

    imin=1-nhe
    imax=nc+nhe
    !
    ! fill halo start
    !
    itot=9
    call initghostbufferTR(cellghostbuf_tmp,1,itot,nhe,nc)
    do ie=nets,nete
       istart = 0
       do ixy=1,2
          do ivertex=1,4
             !
             ! phl_opt: memory usage is not optimal here
             !
             call ghostVpack(cellghostbuf_tmp, fvm(ie)%flux_vertex_cart(:,:,ixy,ivertex) ,nhe,nc,1,1,istart,elem(ie)%desc)
             istart = istart+1
          end do
       end do
       call ghostVpack(cellghostbuf_tmp, fvm(ie)%flux_orient(1,:,:) ,nhe,nc,1,1,istart,elem(ie)%desc)
    end do
    call ghost_exchangeV(hybrid,cellghostbuf_tmp,nhe,nc,itot)
    do ie=nets,nete
       istart = 0
       do ixy=1,2
          do ivertex=1,4
             call ghostVunpack(cellghostbuf_tmp, fvm(ie)%flux_vertex_cart(:,:,ixy,ivertex) ,nhe,nc,1,1,istart,elem(ie)%desc)
             istart = istart+1
          end do
       end do
       call ghostVunpack(cellghostbuf_tmp, fvm(ie)%flux_orient(1,:,:) ,nhe,nc,1,1,istart,elem(ie)%desc)
    enddo
    call freeghostbuffertr(cellghostbuf_tmp)
    
    !
    ! indicator for non-existing cells 
    ! set flux_vertex_cart to corner value in non-existent cells
    !
    do ie=nets,nete
       if (fvm(ie)%cubeboundary==nwest) then
          fvm(ie)%flux_orient     (:  ,1-nhe:0     ,nc+1 :nc+nhe) = -1
          fvm(ie)%flux_vertex_cart(1-nhe:0     ,nc+1 :nc+nhe,1,:) = fvm(ie)%flux_vertex_cart(1,nc,1,4)
          fvm(ie)%flux_vertex_cart(1-nhe:0     ,nc+1 :nc+nhe,2,:) = fvm(ie)%flux_vertex_cart(1,nc,2,4)
       else if (fvm(ie)%cubeboundary==swest) then
          fvm(ie)%flux_orient (:,1-nhe:0     ,1-nhe:0     ) = -1
          fvm(ie)%flux_vertex_cart(1-nhe:0     ,1-nhe:0     ,1,:) = fvm(ie)%flux_vertex_cart(1,1,1,1)
          fvm(ie)%flux_vertex_cart(1-nhe:0     ,1-nhe:0     ,2,:) = fvm(ie)%flux_vertex_cart(1,1,2,1)
       else if (fvm(ie)%cubeboundary==neast) then
          fvm(ie)%flux_orient (:,nc+1 :nc+nhe,nc+1 :nc+nhe) = -1
          fvm(ie)%flux_vertex_cart(nc+1 :nc+nhe,nc+1 :nc+nhe,1,:) = fvm(ie)%flux_vertex_cart(nc,nc,1,3)
          fvm(ie)%flux_vertex_cart(nc+1 :nc+nhe,nc+1 :nc+nhe,2,:) = fvm(ie)%flux_vertex_cart(nc,nc,2,3)
       else if (fvm(ie)%cubeboundary==seast) then
          fvm(ie)%flux_orient (:,nc+1 :nc+nhe,1-nhe:0     ) = -1
          fvm(ie)%flux_vertex_cart(nc+1 :nc+nhe,1-nhe:0     ,1,:) = fvm(ie)%flux_vertex_cart(nc,1,1,2)
          fvm(ie)%flux_vertex_cart(nc+1 :nc+nhe,1-nhe:0     ,2,:) = fvm(ie)%flux_vertex_cart(nc,1,2,2)
       end if
    end do
    
    !
    ! set vectors for perpendicular flux vector
    !
    rot90_matrix(1,1) = 0; rot90_matrix(2,1) =  1 !counter-clockwise rotation matrix
    rot90_matrix(1,2) =-1; rot90_matrix(2,2) =  0 !counter-clockwise rotation matrix 
    
    iside = 1
    unit_vec(1,iside) = 0 !x-component of displacement vector for side 1
    unit_vec(2,iside) = 1 !y-component of displacement vector for side 1
    
    do iside=2,4
       unit_vec(:,iside) = MATMUL(rot90_matrix(:,:),unit_vec(:,iside-1))
    end do
    
    !
    ! fill halo done
    !
    !-------------------------------

    do ie=nets,nete
       !       fvm(ie)%weight_displ(:,:,:,:,:) = 1.0E9!dbg
       do j=imin,imax
          do i=imin,imax
             ishft = NINT(fvm(ie)%flux_orient(2,i,j))
             do ixy=1,2
                !
                ! rotate coordinates if needed through permutation
                !
                fvm(ie)%flux_vertex_cart(i,j,ixy,1:4) = cshift(fvm(ie)%flux_vertex_cart(i,j,ixy,1:4),shift=ishft)
                fvm(ie)%flux_vec        (ixy,i,j,1:4) = cshift(unit_vec                (ixy,1:4    ),shift=ishft)
                !
                ! set flux vector to zero in non-existent cells (corner halo) 
                !
                fvm(ie)%flux_vec        (ixy,i,j,1:4) = fvm(ie)%ifct(i,j)*fvm(ie)%flux_vec(ixy,i,j,1:4)
             end do
             
          end do
       end do
              !
       ! pre-compute weights for reconstruction at cell vertices
       !
       !  ! Evaluate constant order terms
       !  value = fcube(a,b) + &
       !  ! Evaluate linear order terms
       !          recons(a,b,1) * (cartx - centroid(a,b,1)) + &
       !          recons(a,b,2) * (carty - centroid(a,b,2)) + &
       !  ! Evaluate second order terms
       !          recons(a,b,3) * (centroid(a,b,1)**2 - centroid(a,b,3)) + &
       !          recons(a,b,4) * (centroid(a,b,2)**2 - centroid(a,b,4)) + &
       !          recons(a,b,5) * (centroid(a,b,1) * centroid(a,b,2) - centroid(a,b,5)) + &
       !
       !          recons(a,b,3) * (cartx - centroid(a,b,1))**2 + &
       !          recons(a,b,4) * (carty - centroid(a,b,2))**2 + &
       !          recons(a,b,5) * (cartx - centroid(a,b,1)) * (carty - centroid(a,b,2))
    end do
    select case (irecons)
    case(3)
       do ie=nets,nete
          do j= 1-nhe,nc+nhe
             do i=1-nhe,nc+nhe  
                count = 1
                do n = j, j+1
                   do m = i, i+1
                      cartx = fvm(ie)%flux_vertex_cart(i,j,1,count); carty = fvm(ie)%flux_vertex_cart(i,j,2,count);
                      
                      fvm(ie)%vertex_recons_weights(count,1,i,j) = cartx - fvm(ie)%spherecentroid(i,j,1)
                      fvm(ie)%vertex_recons_weights(count,2,i,j) = carty - fvm(ie)%spherecentroid(i,j,2)

                      count=count+1
                   end do
                enddo
             end do
          end do
       end do
    case(6)
       do ie=nets,nete
          do j= 1-nhe,nc+nhe
             do i=1-nhe,nc+nhe  
                do count=1,4
                   cartx = fvm(ie)%flux_vertex_cart(i,j,1,count); carty = fvm(ie)%flux_vertex_cart(i,j,2,count);
                   
                   fvm(ie)%vertex_recons_weights(count,1,i,j) = cartx - fvm(ie)%spherecentroid(i,j,1)
                   fvm(ie)%vertex_recons_weights(count,2,i,j) = carty - fvm(ie)%spherecentroid(i,j,2)
                   
                   fvm(ie)%vertex_recons_weights(count,3,i,j) = (fvm(ie)%spherecentroid(i,j,1)**2 - &
                        fvm(ie)%spherecentroid(i,j,3))   + &
                        (cartx - fvm(ie)%spherecentroid(i,j,1))**2 
                   fvm(ie)%vertex_recons_weights(count,4,i,j) = (fvm(ie)%spherecentroid(i,j,2)**2 - &
                        fvm(ie)%spherecentroid(i,j,4)) + &
                        (carty - fvm(ie)%spherecentroid(i,j,2))**2 

                   fvm(ie)%vertex_recons_weights(count,5,i,j) = (cartx - fvm(ie)%spherecentroid(i,j,1))*     &
                        (carty - fvm(ie)%spherecentroid(i,j,2))+     &
                        (fvm(ie)%spherecentroid(i,j,1) *             &
                        fvm(ie)%spherecentroid(i,j,2) - &
                        fvm(ie)%spherecentroid(i,j,5))
                end do
             end do
          end do
       end do

       do ie=nets,nete
          do j= 1-nhe,nc+nhe
             do i=1-nhe,nc+nhe  
                   fvm(ie)%recons_metrics(i,j,1) = fvm(ie)%spherecentroid(i,j,1)**2 -fvm(ie)%spherecentroid(i,j,3)
                   fvm(ie)%recons_metrics(i,j,2) = fvm(ie)%spherecentroid(i,j,2)**2 -fvm(ie)%spherecentroid(i,j,4)
                   fvm(ie)%recons_metrics(i,j,3) = fvm(ie)%spherecentroid(i,j,1)*fvm(ie)%spherecentroid(i,j,2)-&
                                                   fvm(ie)%spherecentroid(i,j,5)

                   fvm(ie)%recons_metrics_integral(i,j,1) = &
                        2.0D0*fvm(ie)%spherecentroid(i,j,1)**2 -fvm(ie)%spherecentroid(i,j,3)
                   fvm(ie)%recons_metrics_integral(i,j,2) = &
                        2.0D0*fvm(ie)%spherecentroid(i,j,2)**2 -fvm(ie)%spherecentroid(i,j,4)
                   fvm(ie)%recons_metrics_integral(i,j,3) = &
                        2.0D0*fvm(ie)%spherecentroid(i,j,1)*fvm(ie)%spherecentroid(i,j,2)-&
                                                   fvm(ie)%spherecentroid(i,j,5)
             end do
          end do
       end do


    case default
       write(*,*) "irecons out of range",irecons
       stop
    end select


  end subroutine fvm_init_flux
  

end module fvm_mod
