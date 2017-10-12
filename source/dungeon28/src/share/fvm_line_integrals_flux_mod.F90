#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE FVM_LINE_INTEGRALS_MOD_FLUX-------------------------------------------------!
!                                                                                   !
!-----------------------------------------------------------------------------------!
module fvm_line_integrals_flux_mod

  use kinds, only               : int_kind, real_kind
  use dimensions_mod, only      : nc, nhe, ngpc
  use fvm_line_integrals_mod, only: ldbg
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, ne, nc, nhe, nlev, ntrac, np, ntrac_d,ns, nhr, nhc
  use perf_mod, only : t_startf, t_stopf ! _EXTERNAL
  use perf_utils, only: t_detail_low, t_detail_medium, t_detail_high, t_detail_max  ! EXTERNAL

  implicit none
  private
  real (kind=real_kind),parameter, public   :: bignum = 1.0D20
  real (kind=real_kind),parameter, public   :: tiny   = 1.0D-12
  real (kind=real_kind),parameter           :: fuzzy_width = 10.0*tiny
  logical                                   :: lexact_horizontal_line_integrals=.FALSE.
  integer, private, parameter :: num_weights_flux = 4*(nc+2*nhe)*(nc+nhe)    

  public :: compute_weights_fluxform,cslam_runflux,ff_cslam_remap
contains
  !
  !**************************************************************************************
  !
  !
  ! ff-cslam subroutines
  !
  !
  !**************************************************************************************
  !  
  subroutine cslam_runflux(elem,fvm,hybrid,deriv,dt_fvm,tl,nets,nete,p_top,loverwrite_se_flux)
    ! ---------------------------------------------------------------------------------
    use fvm_control_volume_mod     , only: n0_fvm, np1_fvm
    ! ---------------------------------------------------------------------------------
    use fvm_reconstruction_mod, only: reconstruction_gradient
    ! ---------------------------------------------------------------------------------
    use derivative_mod, only : derivative_t
    ! -----------------------------------------------
    use edge_mod, only :  ghostVpack2d_level, ghostVunpack2d_level,initghostbufferTR,freeghostbuffertr
    use edgetype_mod, only : ghostBuffertr_t
    use control_mod, only : qsplit
    use time_mod   , only : TimeLevel_Qdp,timelevel_t
    use element_mod, only : element_t, timelevels
    use fvm_control_volume_mod, only: fvm_struct
    use hybrid_mod, only : hybrid_t
    use fvm_line_integrals_mod, only : fvm_mesh_dep
    use fvm_mod, only : fill_halo_fvm
    
    implicit none
    type (element_t), intent(inout)                :: elem(:)
    type (fvm_struct), intent(inout)             :: fvm(:)
    type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
    integer, intent(in)                         :: nets  ! starting thread element number (private)
    integer, intent(in)                         :: nete  ! ending thread element number   (private)
    real (kind=real_kind), intent(in)           :: dt_fvm
    real (kind=real_kind), intent(in)           :: p_top
    logical, intent(in) :: loverwrite_se_flux

    integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h, ntmp
    type (TimeLevel_t)                          :: tl              ! time level struct
    type (derivative_t)                         :: deriv           ! derivative struct
    
    real (kind=real_kind)   , dimension(num_weights_flux,6,2)  :: weights_all_flux
    integer (kind=int_kind),  dimension(num_weights_flux,2,2)  :: weights_eul_index_all_flux
    integer (kind=int_kind),  dimension(num_weights_flux,2,2)  :: weights_lgr_index_all_flux
    integer (kind=int_kind), dimension(2)                            :: jall

    integer (kind=int_kind)                            :: jall_max !dbg
    
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,6)      :: recons
    
    real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer0 
    
    real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer_air0   
    real (kind=real_kind), dimension(1:nc+1,1:nc+1,2)                  :: flux_air
    real (kind=real_kind), dimension(1:nc+1,1:nc+1,2)                  :: flux_tracer
    real (kind=real_kind) :: q0,q1,rho0,rho1,area,diff_dbg 

    real (kind=real_kind)                                              :: flux_area(1:nc+1,1:nc+1,2)
    
    type (ghostBuffertr_t)                      :: buflatlon
    
!xx    call initghostbufferTR(buflatlon,nlev,2,2,nc+1)    ! use the tracer entry 2 for lat lon
    
    call t_startf('ff-cslam scheme', t_detail_high) 
    
!    call TimeLevel_Qdp(tl, qsplit, n0_fvm, np1_fvm)    

    do ie=nets, nete
       do k=1,nlev
          call fvm_mesh_dep(elem(ie),deriv,fvm(ie),dt_fvm,tl,k)
       end do
       !       fvm(ie)%dsphere(:,:,:)%r=1.0D0  !!! RADIUS IS ASSUMED TO BE 1.0DO !!!!       
    end do
    !
    ! fill halo for dp_fvm and c
    !
    call fill_halo_fvm(elem,fvm,hybrid,nets,nete,n0_fvm)

    jall_max=0
    do ie=nets, nete
       do k=1,nlev
          !
          ! zero velocity - for debugging
          !
          !          if (fvm(ie)%cubeboundary > 0) then
          !             do j=1,nc+1                                                                             
          !                do i=1,nc+1                                                                              
          !                   if (i==1.or.j==1.or.i==nc+1.or.j==nc+1) then
          !                      fvm%dsphere(i,j,k)=fvm%asphere(i,j) !zero velocity
          !                   end if
          !                end do
          !             end do
          !          end if
          
          call compute_weights_fluxform(fvm(ie),6,weights_all_flux,weights_eul_index_all_flux, &
               weights_lgr_index_all_flux,k,jall)
          
          if (jall(1)>num_weights_flux) then
             write(*,*) "jall(1)>num_weights_flux"
             stop
          endif
          if (jall(2)>num_weights_flux) then
             write(*,*) "jall(2)>num_weights_flux"
             stop
          endif
          
          tracer_air0=fvm(ie)%dp_fvm(:,:,k,n0_fvm)       
          call reconstruction_gradient(tracer_air0, fvm(ie),recons,6,.false.)
!          recons=0.0D0
          !
          ! do remapping for x (j=1) and y (j=2) fluxes
          !
          call ff_cslam_remap(tracer_air0,flux_air,weights_all_flux,recons, &
               fvm(ie)%spherecentroid, weights_eul_index_all_flux,&
               weights_lgr_index_all_flux, jall)  
          !
          ! compute flux-areas (to be used for tracer advection)
          !
          call ff_cslam_flux_area(flux_area,weights_all_flux, &
               weights_eul_index_all_flux, weights_lgr_index_all_flux, jall)
          !
          ! forecast equation for air density
          !
          do j=1,nc
             do i=1,nc
!                fvm(ie)%dp_fvm(i,j,k,tl%np1)=tracer_air0(i,j)-(&  
                fvm(ie)%dp_fvm(i,j,k,np1_fvm)=tracer_air0(i,j)-(&  
                     flux_air(i+1,j,1)-flux_air(i,j,1)+flux_air(i,j+1,2)-flux_air(i,j,2)&
                     )/fvm(ie)%area_sphere(i,j)
             end do
          end do
          !
          ! Interface with J.Overfelt SE fluxes
          !
          if (loverwrite_se_flux) then
             elem(ie)%sub_elem_mass_flux(1:nc,1:nc,1,k) =   flux_air(1:nc  ,1:nc  ,2)
             elem(ie)%sub_elem_mass_flux(1:nc,1:nc,2,k) =  -flux_air(2:nc+1,1:nc  ,1)
             elem(ie)%sub_elem_mass_flux(1:nc,1:nc,3,k) =  -flux_air(1:nc  ,2:nc+1,2)
             elem(ie)%sub_elem_mass_flux(1:nc,1:nc,4,k) =   flux_air(1:nc  ,1:nc  ,1)
          end if



          !
          !loop through all tracers
          !
          do itr=1,ntrac
!             tracer0=fvm(ie)%c(:,:,k,itr,tl%n0)
             tracer0=fvm(ie)%c(:,:,k,itr,n0_fvm)
             call reconstruction_gradient(tracer0, fvm(ie),recons,6,.true.)
!             recons=0.0
             !
             ! do remapping for x-y fluxes
             !
             call ff_cslam_remap_q(tracer0,flux_tracer,weights_all_flux, &
                  recons, &
                  fvm(ie)%spherecentroid, weights_eul_index_all_flux,&
                  weights_lgr_index_all_flux, jall, flux_area)  
             !
             ! forecast equation for tracer (kg/kg)
             !
             do j=1,nc
                do i=1,nc
                   q0   = fvm(ie)%c(i,j,k,itr,n0_fvm)
                   rho0 = fvm(ie)%dp_fvm(i,j,k,n0_fvm)
                   rho1 = fvm(ie)%dp_fvm(i,j,k,np1_fvm)
                   area = fvm(ie)%area_sphere(i,j)
                   !
                   q1=q0*rho0-(&  
                        !
                        flux_air(i+1,j  ,1)*flux_tracer(i+1,j  ,1)-&
                        flux_air(i  ,j  ,1)*flux_tracer(i  ,j  ,1)+&
                        flux_air(i  ,j+1,2)*flux_tracer(i  ,j+1,2)-&
                        flux_air(i  ,j  ,2)*flux_tracer(i  ,j  ,2)&
                        )/area
                   fvm(ie)%c(i,j,k,itr,np1_fvm)=q1/rho1
                end do
             end do
          enddo  !End Tracer
       end do  !End Level
       !
       ! Compute surface pressure implied by fvm
       !
       do j=1,nc
          do i=1,nc
             fvm(ie)%psc(i,j) = sum(fvm(ie)%dp_fvm(i,j,:,np1_fvm)) +  p_top
          end do
       end do
    end do
    !
    ! advance fvm time-levels
    !
    if (.not.loverwrite_se_flux) then
       ntmp     = np1_fvm
       np1_fvm  = n0_fvm
       n0_fvm   = ntmp
    end if


    call t_stopf('ff-cslam scheme', t_detail_high)
  end subroutine cslam_runflux

  subroutine ff_cslam_remap(tracer0,flux_air,weights_all_flux,recons,centroid, &
       weights_eul_index_all_flux, weights_lgr_index_all_flux, jall)
    
    real (kind=real_kind)                                  ,   intent(in):: tracer0(1-nhc:nc+nhc,1-nhc:nc+nhc)
    real (kind=real_kind)  , dimension(1:nc+1,1:nc+1,2)    , intent(inout):: flux_air
    integer (kind=int_kind), dimension(2)                  , intent(in)  :: jall  
    real (kind=real_kind)                                  , intent(in)  :: recons(1-nhe:nc+nhe,1-nhe:nc+nhe,6)
    real (kind=real_kind)  ,                                 intent(in)  :: centroid(1-nhe:nc+nhe,1-nhe:nc+nhe,5)
    real (kind=real_kind)  , dimension(num_weights_flux,6,2), intent(in)  :: weights_all_flux
    integer (kind=int_kind), dimension(num_weights_flux,2,2), intent(in)  :: weights_eul_index_all_flux
    integer (kind=int_kind), dimension(num_weights_flux,2,2), intent(in)  :: weights_lgr_index_all_flux
    !
    ! should make  dimension(num_weights_flux,2,2) to dimension(2,num_weights_flux,2)
    !
    integer :: h, jx, jy, jdx, jdy, k
    
    flux_air = 0.0D0
    do k=1,2
       do h=1,jall(k)
          jx  = weights_lgr_index_all_flux(h,1,k)
          jy  = weights_lgr_index_all_flux(h,2,k)
          jdx = weights_eul_index_all_flux(h,1,k)
          jdy = weights_eul_index_all_flux(h,2,k)
          
          flux_air(jx,jy,k) = flux_air(jx,jy,k)+weights_all_flux(h,1,k)*(&
               ! all constant terms 
               tracer0(jdx,jdy) - recons(jdx,jdy,2)*centroid(jdx,jdy,1) &
               - recons(jdx,jdy,3)*centroid(jdx,jdy,2) &
               + recons(jdx,jdy,4)*(2.0D0*centroid(jdx,jdy,1)**2 -centroid(jdx,jdy,3)) &
               + recons(jdx,jdy,5)*(2.0D0*centroid(jdx,jdy,2)**2 -centroid(jdx,jdy,4)) &
               + recons(jdx,jdy,6)*(2.0D0*centroid(jdx,jdy,1)*centroid(jdx,jdy,2)-centroid(jdx,jdy,5))) + &
               ! linear terms
               weights_all_flux(h,2,k)*&
               (recons(jdx,jdy,2) - recons(jdx,jdy,4)*2.0D0*centroid(jdx,jdy,1) &
               - recons(jdx,jdy,6)*centroid(jdx,jdy,2)) + &
               weights_all_flux(h,3,k)*&
               (recons(jdx,jdy,3) - recons(jdx,jdy,5)*2.0D0*centroid(jdx,jdy,2) &
               - recons(jdx,jdy,6)*centroid(jdx,jdy,1)) + &
               ! quadratic terms
               weights_all_flux(h,4,k)*recons(jdx,jdy,4)+&
               weights_all_flux(h,5,k)*recons(jdx,jdy,5)+&
               weights_all_flux(h,6,k)*recons(jdx,jdy,6)
       end do
    end do
  end subroutine ff_cslam_remap
  !
  ! compute flux areas
  !
  subroutine ff_cslam_flux_area(flux_area,weights_all_flux, &
       weights_eul_index_all_flux, weights_lgr_index_all_flux, jall)
    
    real (kind=real_kind)                                   , intent(inout) :: flux_area(1:nc+1,1:nc+1,2)
    integer (kind=int_kind), dimension(2)                   , intent(in)    :: jall  
    real (kind=real_kind)  , dimension(num_weights_flux,6,2), intent(in)    :: weights_all_flux
    integer (kind=int_kind), dimension(num_weights_flux,2,2), intent(in)    :: weights_eul_index_all_flux
    integer (kind=int_kind), dimension(num_weights_flux,2,2), intent(in)    :: weights_lgr_index_all_flux
 
    integer                                     :: h, jx, jy, jdx, jdy, k
    
    flux_area = 0.0D0
    do k=1,2
       do h=1,jall(k)
          jx  = weights_lgr_index_all_flux(h,1,k)
          jy  = weights_lgr_index_all_flux(h,2,k)
          jdx = weights_eul_index_all_flux(h,1,k)
          jdy = weights_eul_index_all_flux(h,2,k)
          
          flux_area(jx,jy,k) =  flux_area(jx,jy,k)+weights_all_flux(h,1,k)
       end do
    end do
  end subroutine ff_cslam_flux_area

  subroutine ff_cslam_remap_q(tracer0,flux,weights_all_flux,recons,centroid, &
       weights_eul_index_all_flux, weights_lgr_index_all_flux, jall,flux_area)
    
    real (kind=real_kind), intent(in)           :: tracer0(1-nhc:nc+nhc,1-nhc:nc+nhc)
    real (kind=real_kind), intent(inout)        :: flux(1:nc+1,1:nc+1,2)
    integer (kind=int_kind), dimension(2), intent(in)         :: jall  
    real (kind=real_kind), intent(in)           :: recons(1-nhe:nc+nhe,1-nhe:nc+nhe,6)
    real (kind=real_kind), intent(in)           :: centroid(1-nhe:nc+nhe,1-nhe:nc+nhe,5)
    real (kind=real_kind)  , dimension(num_weights_flux,6,2), intent(in)  :: weights_all_flux
    integer (kind=int_kind), dimension(num_weights_flux,2,2), intent(in)  :: weights_eul_index_all_flux
    integer (kind=int_kind), dimension(num_weights_flux,2,2), intent(in)  :: weights_lgr_index_all_flux
    real (kind=real_kind)                                   , intent(in)  :: flux_area(1:nc+1,1:nc+1,2)

    integer                                     :: h, jx, jy, jdx, jdy, k
    
    flux   = 0.0D0
    do k=1,2
       do h=1,jall(k)
          jx  = weights_lgr_index_all_flux(h,1,k)
          jy  = weights_lgr_index_all_flux(h,2,k)
          jdx = weights_eul_index_all_flux(h,1,k)
          jdy = weights_eul_index_all_flux(h,2,k)
          
          flux(jx,jy,k) = flux(jx,jy,k)+weights_all_flux(h,1,k)*(&
               ! all constant terms 
               tracer0(jdx,jdy) - recons(jdx,jdy,2)*centroid(jdx,jdy,1) &
               - recons(jdx,jdy,3)*centroid(jdx,jdy,2) &
               + recons(jdx,jdy,4)*(2.0D0*centroid(jdx,jdy,1)**2 -centroid(jdx,jdy,3)) &
               + recons(jdx,jdy,5)*(2.0D0*centroid(jdx,jdy,2)**2 -centroid(jdx,jdy,4)) &
               + recons(jdx,jdy,6)*(2.0D0*centroid(jdx,jdy,1)*centroid(jdx,jdy,2)-centroid(jdx,jdy,5))) + &
               ! linear terms
               weights_all_flux(h,2,k)*&
               (recons(jdx,jdy,2) - recons(jdx,jdy,4)*2.0D0*centroid(jdx,jdy,1) &
               - recons(jdx,jdy,6)*centroid(jdx,jdy,2)) + &
               weights_all_flux(h,3,k)*&
               (recons(jdx,jdy,3) - recons(jdx,jdy,5)*2.0D0*centroid(jdx,jdy,2) &
               - recons(jdx,jdy,6)*centroid(jdx,jdy,1)) + &
               ! quadratic terms
               weights_all_flux(h,4,k)*recons(jdx,jdy,4)+&
               weights_all_flux(h,5,k)*recons(jdx,jdy,5)+&
               weights_all_flux(h,6,k)*recons(jdx,jdy,6)
       end do
       do jy=1,nc+1
          do jx=1,nc+1
             if (abs(flux_area(jx,jy,k))<1.0E-12) then
!                flux(jx,jy,k) = 0.0D0
                flux(jx,jy,k) = tracer0(jx,jy)
             else
                flux(jx,jy,k) = flux(jx,jy,k)/flux_area(jx,jy,k)
             end if
          end do
       end do
    end do
  end subroutine ff_cslam_remap_q

  
  subroutine fluxform_cslam_area(area,fvm_area,weights,weights_eul_index_all, &
       weights_lgr_index_all, jall)
    implicit none
    real (kind=real_kind), intent(inout) :: area(1:nc,1:nc)    
    real (kind=real_kind), intent(in)    :: fvm_area(1:nc,1:nc)    
    real (kind=real_kind), intent(in)    :: weights(10*(nc+2*nhe)*(nc+2*nhe),6,2)
    integer (kind=int_kind), intent(in)  :: weights_eul_index_all(10*(nc+2*nhe)*(nc+2*nhe),2,2)
    integer (kind=int_kind), intent(in)  :: weights_lgr_index_all(10*(nc+2*nhe)*(nc+2*nhe),2,2)
    integer (kind=int_kind), dimension(2), intent(in)        :: jall  
    
    real (kind=real_kind)                :: flux_area(1:nc+1,1:nc+1,2)
    integer                              :: h, jx, jy, jdx, jdy,j
    
    flux_area = 0.0D0
    do j=1,2
       do h=1,jall(j)
          jx  = weights_lgr_index_all(h,1,j)
          jy  = weights_lgr_index_all(h,2,j)
          jdx = weights_eul_index_all(h,1,j)
          jdy = weights_eul_index_all(h,2,j)
          
          flux_area(jx,jy,j) = flux_area(jx,jy,j)+weights(h,1,j)
       end do
    end do
    
    area = 0.0D0
    do jy=1,nc
       do jx=1,nc
          area(jx,jy)=fvm_area(jx,jy)-&
               (flux_area(jx+1,jy  ,1)-flux_area(jx,jy,1)+&
               flux_area(jx  ,jy+1,2)-flux_area(jx,jy,2))
       end do
    end do
  end subroutine fluxform_cslam_area
  


  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE COMPUTE_WEIGHTS_FLUXFORM------------------------------------------------!
  !                                                                                   !
  ! THIS CODE IS BASED ON COMPUTE_WEIGHTS CODED BY CHRISTOPH ERATH AND MODIFIED FOR   !
  ! FLUX-FORM CSLAM BY PETER HJORT LAURITZEN                                          !
  !                                                                                   !
  ! DESCRIPTION:                                                                      !
  !-----------------------------------------------------------------------------------!
  subroutine compute_weights_fluxform(fvm,nreconstruction,weights_all,weights_eul_index_all, &
       weights_lgr_index_all,klev,jall)  
    use fvm_control_volume_mod, only:  fvm_struct                                         
    use coordinate_systems_mod,  only :  cartesian2D_t, &
         cart2cubedspherexy, spherical_to_cart
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    
    use fvm_line_integrals_mod, only: gauss_points

    use dimensions_mod, only: ip_dbg
    
    implicit none
    type (fvm_struct), intent(inout)                                :: fvm
    integer (kind=int_kind), intent(in)                            :: nreconstruction
    ! arrays for collecting cell data
    !
    !
    ! maximum number of overlaps for CN<1 is 4 per flux edge
    !
    ! nhe is CEILING(CN) which is 1 here
    !
    !
    ! dimension(4*(nc+2*nhe)*(nc+nhe),nreconstruction,2)
    !
    real (kind=real_kind),dimension(:,:,:), intent(out) :: weights_all
    !
    ! dimension(4*(nc+2*nhe)*(nc+nhe),2,2)
    !
    integer (kind=int_kind), dimension(:,:,:), intent(out) :: weights_eul_index_all, weights_lgr_index_all
    integer (kind=int_kind), intent(in)                       :: klev
    integer (kind=int_kind), dimension(:), intent(out)        :: jall !dimension(2)
    
    ! local workspace
    ! max number of line segments is:
    ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
    !
    
    integer (kind=int_kind)                     :: jx,jy
    integer                                     :: jx_min, jx_max, jy_min, jy_max
    integer                                     :: jx_min1, jx_max1, jy_min1, jy_max1
    integer                                     :: jx_min2, jx_max2, jy_min2, jy_max2
    logical                                     :: swap1, swap2
    
    integer (kind=int_kind)                     :: i, j
    
    type (cartesian2D_t)                        :: dcart(-1:nc+3,-1:nc+3)       ! Cartesian coordinates 

    type (cartesian2D_t)                        :: acart(-1:nc+3,-1:nc+3)       ! Cartesian coordinates 
    
    real (kind=real_kind), dimension(0:5)       :: xflux_x,xflux_y,yflux_x,yflux_y
    integer (kind=int_kind)                     :: inttmp
    real (kind=real_kind)                       :: tmp
    ! for Gaussian quadrature
    real (kind=real_kind), dimension(ngpc)      :: gsweights, gspts
    ! weight-variables for individual cells
    integer (kind=int_kind) :: jmax_segments_cell
    real (kind=real_kind)   , dimension(nhe*50,nreconstruction,2)   :: weights_flux_cell
    integer (kind=int_kind),  dimension(nhe*50,2,2)                 :: weights_eul_index_cell
    integer (kind=int_kind), dimension(2)                           :: jcollect_cell
    logical :: lxflux,lyflux
    !
    !
    !    x--x--x--x--x
    !    |  |  |  |  |
    !    x--x--x--x--x
    !    |  |  |  |  |
    !    x--x--x--x--x
    !    |  |  |  |  |
    !    x--x--x--x--x
    !    |  |  |  |  |
    !    x--x--x--x--x
    !
    !
    !            =========
    !            |       |
    !            |       |
    !            |       |
    !            |       |
    !            |       |
    !    =================================
    !    |       |       |       |       |
    !    |       |       |       |       |
    !    |       |       |       |       |
    !    |       |       |       |       |
    !    |       |       |       |       |
    !    =================================
    !            |       |
    !            |       |
    !            |       |
    !            |       |
    !            |       |
    !            =========
    !
    jx_min=fvm%jx_min; jx_max=fvm%jx_max; 
    jy_min=fvm%jy_min; jy_max=fvm%jy_max;
    !
    ! fvm%cubeboundary=0 means interior element
    !
    if (fvm%cubeboundary > 0) then
       !
       ! element is on panel side
       !
       jx_min1=fvm%jx_min1; jx_max1=fvm%jx_max1; 
       jy_min1=fvm%jy_min1; jy_max1=fvm%jy_max1;
       swap1=fvm%swap1
       if (fvm%cubeboundary > 4) then
          !
          ! element is on a panel corner
          !
          jx_min2=fvm%jx_min2; jx_max2=fvm%jx_max2;
          jy_min2=fvm%jy_min2; jy_max2=fvm%jy_max2;
          swap2=fvm%swap2
       endif
    endif
    
    jmax_segments_cell = nhe*50
    
    call gauss_points(ngpc,gsweights,gspts)
    tmp =0.0D0
    jall = 1
    ! 
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=1,nc+1
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%faceno,dcart(jx,jy))  
          acart(jx,jy)%x = fvm%acartx(jx)
          acart(jx,jy)%y = fvm%acarty(jy)
       end do
    end do
    ip_dbg = fvm%faceno !dbg
    do jy=1, nc+1
       do jx=1, nc+1 
          
          !
          !    o-y-o-y-o-y-o-y-o-  jy=nc+1         
          !    |   |   |   |   |
          !    x   x   x   x   x
          !    |   |   |   |   |
          !    o-y-o-y-o-y-o-y-o-  jy=nc
          !    |   |   |   |   |
          !    x   x   x   x   x
          !    |   |   |   |   |                    x=xflux
          !    o-y-o-y-o-y-o-y-o-  jy=3             y=yflux
          !    |   |   |   |   |
          !    x   x   x   x   x
          !    |   |   |   |   |
          !    o-y-o-y-o-y-o-y-o-  jy=2
          !    |   |   |   |   |
          !    x   x   x   x   x
          !    |   |   |   |   |
          !    o-y-o-y-o-y-o-y-o-  jy=1
          !
          !
          ! lxflux and lyflux is to deal with cells where not both x and y flux are computed!
          !
          lxflux=jy<nc+1
          lyflux=jx<nc+1
          if (lxflux.or.lyflux) then
             if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
                  acart,dcart)     
             if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
                  acart,dcart)     

             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,&
                  yflux_y,jx,jy,&
                  nreconstruction,fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,&
                  jmax_segments_cell) 
             
             
             do j=1,2
                if (jcollect_cell(j)>0) then
                   weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                        weights_flux_cell(1:jcollect_cell(j),:,j)
                   weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                        weights_eul_index_cell(1:jcollect_cell(j),:,j)
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                   jall(j) = jall(j)+jcollect_cell(j)          
                endif
             end do
          end if
       end do
    end do 

    !WEST SIDE
    if (fvm%cubeboundary == west) then
       !
       !
       ! This Figure shows the element to the East 
       ! of the element in question 
       ! (with no index swapping)
       !       
       !                    o-y-o---o---o---o-  jy=nc+1         
       !                    |   |   |   |   |
       !                    x   |   |   |   |
       !                    |   |   |   |   |
       !                    o-y-o---o---o---o-  jy=nc
       !                    |   |   |   |   |
       !                    x   |   |   |   |
       !                    |   |   |   |   |                    x=xflux
       !                    o-y-o---o---o---o-  jy=3             y=yflux
       !                    |   |   |   |   |
       !                    x   |   |   |   |
       !                    |   |   |   |   |
       !                    o-y-o---o---o---o-  jy=2
       !                    |   |   |   |   |
       !                    x   |   |   |   |
       !                    |   |   |   |   |
       !                    o-y-o---o---o---o-  jy=1
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       !
       jx=1
       do jy=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(west),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
               fvm%nbrsface(west),dcart(jx+1,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(west),acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx+1,jy)),&
               fvm%nbrsface(west),acart(jx+1,jy))                  
          ip_dbg = fvm%nbrsface(west) !dbg
       end do
       
       do jy=1,nc+1
          lxflux=jy<nc+1
          lyflux=.true.
          if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
               acart,dcart)     
          if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
               acart,dcart)     

          if(swap1) then  !flip orientation
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min1,jx_min1,&
                  nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else  
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)

          end if
          if (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)+nhe-1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=jy_max1-jy_min1-weights_eul_index_cell(i,2,j)-nhe-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy

                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
    endif


    !EAST SIDE
    if (fvm%cubeboundary == east) then
       !
       !
       ! This Figure show the element to the East 
       ! of the element in question 
       ! (with no index swapping)
       !       
       !                    o---o---o---o-y-o-  jy=nc+1         
       !                    |   |   |   |   |
       !                    |   |   |   |   x
       !                    |   |   |   |   |
       !                    o---o---o---o-y-o-  jy=nc
       !                    |   |   |   |   |
       ! x=xflux            |   |   |   |   x               
       ! y=yflux            |   |   |   |   |                    
       !                    o---o---o---o-y-o-  jy=3        
       !                    |   |   |   |   |
       !                    |   |   |   |   x
       !                    |   |   |   |   |
       !                    o---o---o---o-y-o-  jy=2
       !                    |   |   |   |   |
       !                    |   |   |   |   x
       !                    |   |   |   |   |
       !                    o---o---o---o-y-o-  jy=1
       !        
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face 
       do jx=nc,nc+1
          do jy=1,nc+1
             call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
                  fvm%nbrsface(east),dcart(jx,jy))                  
!             call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
!                  fvm%nbrsface(east),dcart(jx+1,jy))                  
             call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
                  fvm%nbrsface(east),acart(jx,jy))                  
!             call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx+1,jy)),&
!                  fvm%nbrsface(east),acart(jx+1,jy))                  
             ip_dbg = fvm%nbrsface(east) !dbg
          end do
       end do
       do jx=nc,nc+1
          do jy=1,nc+1
             lxflux=jy<nc+1.AND.jx==nc+1
             lyflux=(jx==nc)

             if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
                  acart,dcart)     
             if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
                  acart,dcart)     

             if(swap1) then !flip orientation
                call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min1,&
                     jx_min1,nreconstruction,&
                     fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                     tmp,ngpc,gsweights,gspts,&
                     weights_flux_cell,weights_eul_index_cell,jcollect_cell,&
                     jmax_segments_cell)


                do j=1,2
                   do i=1,jcollect_cell(j)
                      inttmp=weights_eul_index_cell(i,1,j)
                      weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                      weights_eul_index_cell(i,2,j)=inttmp
                   end do
                end do
             else
                call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                     fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                     tmp,ngpc,gsweights,gspts,&
                     weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)



             end if
             !I have to correct the number
             if (fvm%faceno==5) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,2,j)=jy_max1-jy_min1-weights_eul_index_cell(i,2,j)-nhe-nhe+1
                   end do
                end do
             end if
             if (fvm%faceno==6) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)-nhe+1
                   end do
                end do
             end if
             do j=1,2
                if (jcollect_cell(j)>0) then
                   weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                   weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                        weights_eul_index_cell(1:jcollect_cell(j),:,j)
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                   jall(j) = jall(j)+jcollect_cell(j)          
                endif
             end do
          end do
       end do
    endif


    !NORTH SIDE 
    if (fvm%cubeboundary == north) then
       !
       ! no "swapping case"
       !
       !
       !                    o-y-o-y-o-y-o-y-o-  jy=nc+1         
       !                    |   |   |   |   |
       !                    x   x   x   x   x
       !                    |   |   |   |   |
       !                    o---o---o---o---o-  jy=nc
       !                    |   |   |   |   |
       ! x=xflux            |   |   |   |   |               
       ! y=yflux            |   |   |   |   |                    
       !                    o---o---o---o---o-  jy=3        
       !                    |   |   |   |   |
       !                    |   |   |   |   |
       !                    |   |   |   |   |
       !                    o---o---o---o---o-  jy=2
       !                    |   |   |   |   |
       !                    |   |   |   |   |
       !                    |   |   |   |   |
       !                    o---o---o---o---o-  jy=1
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       !
       do jy=nc,nc+1
          do jx=1,nc+1
             call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
                  fvm%nbrsface(north),dcart(jx,jy))                  
!             call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
!                  fvm%nbrsface(north),dcart(jx,jy+1))                  
             call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
                  fvm%nbrsface(north),acart(jx,jy))                  
!             call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy+1)),&
!                  fvm%nbrsface(north),acart(jx,jy+1))                  
             ip_dbg = fvm%nbrsface(north) !dbg
          end do
       end do
       
       do jy=nc,nc+1
          do jx=1,nc+1
             
             lxflux=(jy==nc)
             lyflux=(jy==nc+1).AND.jx<nc+1

             if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
                  acart,dcart)     
             if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
                  acart,dcart)     

             if(swap1) then !flip orientation
                call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min1,&
                     jx_min1,nreconstruction,&
                     fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                     tmp,ngpc,gsweights,gspts,&
                     weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)




                do j=1,2
                   do i=1,jcollect_cell(j)
                      inttmp=weights_eul_index_cell(i,1,j)
                      weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                      weights_eul_index_cell(i,2,j)=inttmp
                   end do
                end do
             else
                call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                     fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                     tmp,ngpc,gsweights,gspts,&
                     weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  




             end if
             if (fvm%faceno==2) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                   end do
                end do
             end if
             if (fvm%faceno==3) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j))-nhe-nhe+1
                      weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                   end do
                end do
             end if
             if (fvm%faceno==4) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                   end do
                end do
             end if
             if (fvm%faceno==6) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                      weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                   end do
                end do
             end if
             do j=1,2
                if (jcollect_cell(j)>0) then
                   weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                   weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                        weights_eul_index_cell(1:jcollect_cell(j),:,j)
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                   jall(j) = jall(j)+jcollect_cell(j)          
                endif
             end do
          end do
       end do
    end if
    !SOUTH SIDE

    if (fvm%cubeboundary == south) then
       !
       ! no "swapping case"
       !
       !
       !                    o---o---o---o---o-  jy=nc+1         
       !                    |   |   |   |   |
       !                    |   |   |   |   |
       !                    |   |   |   |   |
       !                    o---o---o---o---o-  jy=nc
       !                    |   |   |   |   |
       ! x=xflux            |   |   |   |   |               
       ! y=yflux            |   |   |   |   |                    
       !                    o---o---o---o---o-  jy=3        
       !                    |   |   |   |   |
       !                    |   |   |   |   |
       !                    |   |   |   |   |
       !                    o---o---o---o---o-  jy=2
       !                    |   |   |   |   |
       !                    x   x   x   x   x
       !                    |   |   |   |   |
       !                    o-y-o-y-o-y-o-y-o-  jy=1
       !
       !
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       !
       !begin dbg
       jy=1
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%faceno,dcart(jx,jy))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
               fvm%faceno,dcart(jx,jy+1))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%faceno,acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy+1)),&
               fvm%faceno,acart(jx,jy+1))                  
       end do

       jy=1
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(south),dcart(jx,jy))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
               fvm%nbrsface(south),dcart(jx,jy+1))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(south),acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy+1)),&
               fvm%nbrsface(south),acart(jx,jy+1))         
          ip_dbg = fvm%nbrsface(south) !dbg         
       end do
       
       do jx=1,nc+1
          lxflux=.TRUE.
          lyflux=(jx<nc+1)
          
          if (lxflux) then
             call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,acart,dcart)
          end if
          if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
               acart,dcart)     
             
          if(swap1) then !flip orientation
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)

             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  


          end if
          if  (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                end do
             end do
          end if
          if  (fvm%faceno==3) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          if  (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          if  (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             end if
          end do
       end do
       
    endif
    !SOUTHWEST Corner
    if (fvm%cubeboundary == swest) then

       !
       ! no "swapping case"
       !
       !                    o-y-o---o---o---o-  jy=nc+1         
       !                    |   |   |   |   |
       !                    x   |   |   |   |
       !                    |   |   |   |   |
       !                    o-y-o---o---o---o-  jy=nc
       !                    |   |   |   |   |
       ! x=xflux            x   |   |   |   |               
       ! y=yflux            |   |   |   |   |                    
       !                    o-y-o---o---o---o-  jy=3        
       !                    |   |   |   |   |
       !                    x   |   |   |   |
       !                    |   |   |   |   |
       !                    o-y-o---o---o---o-  jy=2
       !                    |   |   |   |   |
       !                    x   x   x   x   x
       !                    |   |   |   |   |
       !                    o-y-o-y-o-y-o-y-o-  jy=1
       !
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face

       !
       ! start with South side
       !
       jy=1
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(south),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
               fvm%nbrsface(south),dcart(jx,jy+1))                  

          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(south),acart(jx,jy))
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy+1)),&
               fvm%nbrsface(south),acart(jx,jy+1))
          ip_dbg = fvm%nbrsface(south) !dbg
       end do


       do jx=1,nc+1

          lxflux=.TRUE.
          lyflux=(jx<nc+1)

          if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
               acart,dcart)     
          if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
               acart,dcart)     

          if(swap1) then !flip orientation
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)


             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)


          end if
          if ((fvm%faceno==3) .OR. (fvm%faceno==5)) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)+1
                   weights_eul_index_cell(i,2,j)=(weights_eul_index_cell(i,2,j)+nhe-1)
                end do
             end do
          end if
          if (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)+1
                end do
             end do
          end if
          if (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
       
       !
       ! West side 
       !       
       ! calculate xy Cartesian on the cube of departure points on the corresponding face
       jx=1
       do jy=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(west),dcart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
               fvm%nbrsface(west),dcart(jx+1,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(west),acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx+1,jy)),&
               fvm%nbrsface(west),acart(jx+1,jy))                  
          ip_dbg = fvm%nbrsface(west) !dbg
       end do

       do jy=1,nc+1
          lxflux=jy<nc+1
          lyflux=.true.

          if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
               acart,dcart)     
          if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
               acart,dcart)     

          if(swap2) then !flip orientation
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)


             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)


          end if
          if (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)+nhe-1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2,j)+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell (1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
    endif
    
    ! SOUTHEAST Corner
    if (fvm%cubeboundary == seast) then
       !
       ! no "swapping case"
       !
       !                    o---o---o---o-y-o-  jy=nc+1         
       !                    |   |   |   |   |
       !                    |   |   |   |   x
       !                    |   |   |   |   |
       !                    o---o---o---o-y-o-  jy=nc
       !                    |   |   |   |   |
       ! x=xflux            |   |   |   |   x               
       ! y=yflux            |   |   |   |   |                    
       !                    o---o---o---o-y-o-  jy=3        
       !                    |   |   |   |   |
       !                    |   |   |   |   x
       !                    |   |   |   |   |
       !                    o---o---o---o-y-o-  jy=2
       !                    |   |   |   |   |
       !                    x   x   x   x   x
       !                    |   |   |   |   |
       !                    o-y-o-y-o-y-o-y-o-  jy=1
       !
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face 
       !
       !
       ! SOUTH side
       !
       jy=1
       do jx=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(south),dcart(jx,jy))  
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
               fvm%nbrsface(south),dcart(jx,jy+1))  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(south),acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy+1)),&
               fvm%nbrsface(south),acart(jx,jy+1))                  
          ip_dbg = fvm%nbrsface(south) !dbg
       end do
       do jx=1,nc+1

          lxflux=.TRUE.
          lyflux=(jx<nc+1)

          if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
               acart,dcart)     
          if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
               acart,dcart)     

          if(swap1) then !flip orientation
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)



             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if


          !I have to correct the number   
          if  (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                end do                   
             end do
          end if
          if  (fvm%faceno==3) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          if  (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          if  (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j)+nhe)-nhe+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)+nhe-1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
       
       ! calculate xy Cartesian on the cube of departure points on the corresponding face  
       !
       ! EAST
       !
       jx=nc
       do jy=1,nc+1
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(east),dcart(jx,jy))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
               fvm%nbrsface(east),dcart(jx+1,jy))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(east),acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx+1,jy)),&
               fvm%nbrsface(east),acart(jx+1,jy))                  
          ip_dbg = fvm%nbrsface(east) !dbg
       end do

       do jx=nc,nc+1
          do jy=1,nc+1

             lxflux=jy<nc+1.AND.jx==nc+1
             lyflux=(jx==nc)

             if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
                  acart,dcart)     
             if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
                  acart,dcart)     
             
             if(swap2) then !flip orientation
                call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min2,jx_min2,nreconstruction,&
                     fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                     tmp,ngpc,gsweights,gspts,&
                     weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
                do j=1,2
                   do i=1,jcollect_cell(j)
                      inttmp=weights_eul_index_cell(i,1,j)
                      weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                      weights_eul_index_cell(i,2,j)=inttmp
                   end do
                end do
             else
                call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                     fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                     tmp,ngpc,gsweights,gspts,&
                     weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             end if


             if (fvm%faceno==5) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,2,j)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2,j)+1
                   end do
                end do
             end if
             if (fvm%faceno==6) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)-nhe+1
                   end do
                end do
             end if
             do j=1,2
                if (jcollect_cell(j)>0) then
                   weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                   weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                        weights_eul_index_cell(1:jcollect_cell(j),:,j)
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                   jall(j) = jall(j)+jcollect_cell(j)          
                endif
             end do
          end do
       end do
    endif
    
    !NORTHEAST Corner
    if (fvm%cubeboundary == neast) then
       !
       !                    o-y-o-y-o-y-o-y-o-  jy=nc+1         
       !                    |   |   |   |   |
       !                    x   x   x   x   x
       !                    |   |   |   |   |
       !                    o---o---o---o-y-o-  jy=nc
       !                    |   |   |   |   |
       ! x=xflux            |   |   |   |   x               
       ! y=yflux            |   |   |   |   |                    
       !                    o---o---o---o-y-o-  jy=3        
       !                    |   |   |   |   |
       !                    |   |   |   |   x
       !                    |   |   |   |   |
       !                    o---o---o---o-y-o-  jy=2
       !                    |   |   |   |   |
       !                    |   |   |   |   x
       !                    |   |   |   |   |
       !                    o---o---o---o-y-o-  jy=1
       !
       ! calculate xy Cartesian on the cube of departure points on the corresponding face  
       !
       ! NORTH part
       !
       jy=nc
       do jx=1,nc+1
             call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
                  fvm%nbrsface(north),dcart(jx,jy))                
             call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
                  fvm%nbrsface(north),dcart(jx,jy+1))                
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(north),acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy+1)),&
               fvm%nbrsface(north),acart(jx,jy+1))                  
          ip_dbg = fvm%nbrsface(north) !dbg
       end do
       
       do jy=nc,nc+1
          do jx=1,nc+1
             
             lxflux=(jy==nc)
             lyflux=(jy==nc+1).AND.jx<nc+1

             if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
                  acart,dcart)     
             if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
                  acart,dcart)     

             if(swap1) then !flip orientation
                call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min1,jx_min1,nreconstruction,&
                     fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                     tmp,ngpc,gsweights,gspts,&
                     weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
                do j=1,2
                   do i=1,jcollect_cell(j)
                      inttmp=weights_eul_index_cell(i,1,j)
                      weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                      weights_eul_index_cell(i,2,j)=inttmp
                   end do
                end do
             else
                call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                     fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                     tmp,ngpc,gsweights,gspts,&
                     weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             end if


             if (fvm%faceno==2) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,2,1)=(weights_eul_index_cell(i,2,j))-nhe+1
                   end do
                end do
             end if
             if (fvm%faceno==3) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                      weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                   end do
                end do
             end if
             if (fvm%faceno==6) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                      weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                   end do
                end do
             end if
             if (fvm%faceno==4) then
                do j=1,2
                   do i=1,jcollect_cell(j)
                      weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)-nhe-nhe+1
                   end do
                end do
             end if
             do j=1,2
                if (jcollect_cell(j)>0) then
                   weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                   weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                        weights_eul_index_cell(1:jcollect_cell(j),:,j)
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                   weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                   jall(j) = jall(j)+jcollect_cell(j)          
                endif
          end do
       end do
    end do
       
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    !
    ! EAST part
    !
    jx=nc
    do jy=1,nc+1
       call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
            fvm%nbrsface(east),dcart(jx,jy))                  
       call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
            fvm%nbrsface(east),dcart(jx+1,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(east),acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx+1,jy)),&
               fvm%nbrsface(east),acart(jx+1,jy))                  
          ip_dbg = fvm%nbrsface(east) !dbg
    end do

    do jx=nc,nc+1
       do jy=1,nc+1
          lxflux=jy<nc+1.AND.jx==nc+1
          lyflux=(jx==nc)   

          if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
               acart,dcart)     
          if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
               acart,dcart)     

          if(swap2) then !flip orientation
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min2,jx_min2,nreconstruction,&
                  fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                  fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if


          if (fvm%faceno==5) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=jy_max2-jy_min2-weights_eul_index_cell(i,2,j)-nhe-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
    end do
 endif
    
 !NORTH WEST CORNER 
 if (fvm%cubeboundary == nwest) then
    !
    !
    !                    o-y-o-y-o-y-o-y-o-  jy=nc+1         
    !                    |   |   |   |   |
    !                    x   x   x   x   x
    !                    |   |   |   |   |
    !                    o-y-o---o---o---o-  jy=nc
    !                    |   |   |   |   |
    ! x=xflux            x   |   |   |   |               
    ! y=yflux            |   |   |   |   |                    
    !                    o-y-o---o---o---o-  jy=3        
    !                    |   |   |   |   |
    !                    x   |   |   |   |
    !                    |   |   |   |   |
    !                    o-y-o---o---o---o-  jy=2
    !                    |   |   |   |   |
    !                    x   |   |   |   |
    !                    |   |   |   |   |
    !                    o-y-o---o---o---o-  jy=1
    !
    !
    ! calculate xy Cartesian on the cube of departure points on the corresponding face
    !
    ! NORTH side
    !
    jy=nc
    do jx=1,nc+1
       call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
            fvm%nbrsface(north),dcart(jx,jy))                   
       call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy+1,klev)),&
            fvm%nbrsface(north),dcart(jx,jy+1))                   
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(north),acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy+1)),&
               fvm%nbrsface(north),acart(jx,jy+1))                  
          ip_dbg = fvm%nbrsface(north) !dbg
    end do
    
    do jy=nc,nc+1
       do jx=1,nc+1
          lxflux=(jy==nc)
          lyflux=(jy==nc+1).AND.jx<nc+1
          
          if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
               acart,dcart)     
          if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
               acart,dcart)     

          if(swap1) then !flip orientation
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min1,jx_min1,nreconstruction,&
                  fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
             do j=1,2
                do i=1,jcollect_cell(j)
                   inttmp=weights_eul_index_cell(i,1,j)
                   weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                   weights_eul_index_cell(i,2,j)=inttmp
                end do
             end do
          else
             call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
                  fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                  tmp,ngpc,gsweights,gspts,&
                  weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if


          if (fvm%faceno==2) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==3) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1,j))+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          if (fvm%faceno==4) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)+1
                end do
             end do
          end if
          if (fvm%faceno==6) then
             do j=1,2
                do i=1,jcollect_cell(j)
                   weights_eul_index_cell(i,1,j)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1,j)+1
                   weights_eul_index_cell(i,2,j)=weights_eul_index_cell(i,2,j)-nhe+1
                end do
             end do
          end if
          do j=1,2
             if (jcollect_cell(j)>0) then
                weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
                weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                     weights_eul_index_cell(1:jcollect_cell(j),:,j)
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
                weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
                jall(j) = jall(j)+jcollect_cell(j)          
             endif
          end do
       end do
    end do
    
     ! calculate xy Cartesian on the cube of departure points on the corresponding face
    !
    ! WEST
    !
    jx=1
    do jy=1,nc+1
       call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
            fvm%nbrsface(west),dcart(jx,jy))                  
       call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx+1,jy,klev)),&
            fvm%nbrsface(west),dcart(jx+1,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx,jy)),&
               fvm%nbrsface(west),acart(jx,jy))                  
          call cart2cubedspherexy(spherical_to_cart(fvm%asphere(jx+1,jy)),&
               fvm%nbrsface(west),acart(jx+1,jy))                  
          ip_dbg = fvm%nbrsface(west) !dbg
    end do


    do jy=1,nc+1
       lxflux=jy<nc+1
       lyflux=.true.
       
       if (lxflux) call getdep_cellboundariesxyvec_xflux(xflux_x,xflux_y,jx,jy,&
            acart,dcart)     
       if (lyflux) call getdep_cellboundariesxyvec_yflux(yflux_x,yflux_y,jx,jy,&
            acart,dcart)     
       
       if(swap2) then !flip orientation
          call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jy_min2,jx_min2,nreconstruction,&
               fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
               tmp,ngpc,gsweights,gspts,&
               weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          do j=1,2
             do i=1,jcollect_cell(j)
                inttmp=weights_eul_index_cell(i,1,j)
                weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,2,j)
                weights_eul_index_cell(i,2,j)=inttmp
             end do
          end do
       else
          call compute_weights_flux_cell(lxflux,lyflux,xflux_x,xflux_y,yflux_x,yflux_y,jx,jy,nreconstruction,&
               fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
               tmp,ngpc,gsweights,gspts,&
               weights_flux_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
       end if
       
       
       
       if (fvm%faceno==5) then
          do j=1,2
             do i=1,jcollect_cell(j)
                weights_eul_index_cell(i,1,j)=weights_eul_index_cell(i,1,j)+nhe-1
             end do
          end do
       end if
       if (fvm%faceno==6) then
          do j=1,2
             do i=1,jcollect_cell(j)
                weights_eul_index_cell(i,2,j)=jy_max2-jy_min2-weights_eul_index_cell(i,2,j)-nhe-nhe+1
             end do
          end do
       end if
       do j=1,2
          if (jcollect_cell(j)>0) then
             weights_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = weights_flux_cell(1:jcollect_cell(j),:,j)
             weights_eul_index_all(jall(j):jall(j)+jcollect_cell(j)-1,:,j) = &
                  weights_eul_index_cell(1:jcollect_cell(j),:,j)
             weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,1,j) = jx
             weights_lgr_index_all(jall(j):jall(j)+jcollect_cell(j)-1,2,j) = jy
             jall(j) = jall(j)+jcollect_cell(j)                 
          endif
       end do
    end do
 endif
 !end if!dbg
 jall=jall-1

    
  end subroutine compute_weights_fluxform


  subroutine compute_weights_flux_cell(ldo_xflux,ldo_yflux,xflux_x,xflux_y,yflux_x,yflux_y,&
       jx,jy,nreconstruction,xgno,ygno,&
       jx_min, jx_max, jy_min, jy_max,tmp,&
       ngauss,gauss_weights,abscissae,weights,weights_eul_index,jcollect,jmax_segments)

    use fvm_line_integrals_mod, only : compute_weights_cell
    implicit none
    logical, intent(in) :: ldo_xflux, ldo_yflux
    integer (kind=int_kind)                 , intent(in):: nreconstruction, jx,jy,ngauss,jmax_segments
    real (kind=real_kind)   ,  dimension(0:), intent(inout):: xflux_x,xflux_y !dimension(0:5)
    real (kind=real_kind)   ,  dimension(0:), intent(inout):: yflux_x,yflux_y !dimension(0:5)
    !
    integer (kind=int_kind), intent(in)               :: jx_min, jy_min, jx_max, jy_max
    real (kind=real_kind), dimension(-nhe:), intent(inout) :: xgno, ygno !dimension(-nhe:nc+2+nhe)
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(:), intent(in) :: gauss_weights, abscissae !dimension(ngauss)
    !
    ! boundaries of domain
    !
    real (kind=real_kind):: tmp
    !
    ! Number of Eulerian sub-cell integrals for the cell in question
    !
    integer (kind=int_kind), dimension(:), intent(out)       :: jcollect !dimension(2)
    !
    ! local workspace
    !
    !
    ! max number of line segments is:
    !
    ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
    !
    real (kind=real_kind)  , dimension(:,:,:), intent(out) :: weights !dimension(jmax_segments,nreconstruction,2)
    integer (kind=int_kind), dimension(:,:,:), intent(out) :: weights_eul_index !dimension(jmax_segments,2,2)
    !
    ! local workspace
    !
    integer :: nvertex,jcollect1,jcollect2
    real (kind=real_kind)   ,  dimension(0:5) :: xcell_flux,ycell_flux
    real (kind=real_kind)   ,  dimension(0:5) :: xcell_flux2,ycell_flux2
    real (kind=real_kind)                     :: weight_sign,weight_sign2
    logical :: lzero_flux
    real (kind=real_kind)   ,  &
         dimension(jmax_segments,nreconstruction,2) :: weights2
    integer (kind=int_kind),  &
         dimension(jmax_segments,2)      :: weights_eul_index2

    lzero_flux = .false.
    !
    ! figure out how to loop over flux-edges (in calling routine)
    !
    !
    !               gno(jx ,jy+1)-----------------gno(jx+1,jy+1)
    !              /|                             |  
    !             / |                             |  
    !            /  |                             |  
    ! cell(jx,jy+1) |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |   gno(jx ,jy  )-----------------gno(jx+1,jy  )
    !           |  /                             /
    !           | /                             /
    !           |/                             /
    !       cell(jx,jy)-------------------cell(jx+1,jy)
    !
    !    
    jcollect=0

    if (ldo_xflux) then
       !
       ! constuct xflux-cell
       !
       xcell_flux = xflux_x; ycell_flux = xflux_y;

       call make_flux_area(jx,jy,xcell_flux,ycell_flux,xcell_flux2,ycell_flux2,&
            lzero_flux,nvertex,weight_sign,weight_sign2)
       
       if (lzero_flux) then
          jcollect1=0
       else
!          if (ldbg) write(*,*) "nvertex is ",nvertex
!          if (ldbg) write(*,*) "flux area after make_flux_area is:"
!          do j=1,nvertex
!             if (ldbg) write(*,*) xcell_flux(j),ycell_flux(j)
!          end do
!          if (ldbg) write(*,*) "going into compute weights cell xflux",jx,jy
          call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
               xcell_flux(1:nvertex),ycell_flux(1:nvertex),jx,jy,nreconstruction,xgno,ygno,&
               jx_min, jx_max, jy_min, jy_max,tmp,&
               ngauss,gauss_weights,abscissae,weights(:,:,1),weights_eul_index(:,:,1),jcollect1,jmax_segments)
          weights(1:jcollect1,:,1) = weight_sign*weights(1:jcollect1,:,1)
!          if (ldbg) write(*,*) "jcollect1 is ",jcollect1
!          if (ldbg) write(*,*) "weights are (from inside)",weights(1:jcollect1,1,1)
!          if (ldbg) write(*,*) "ixweights are (from inside)",weights_eul_index(1:jcollect1,1,1)
!          if (ldbg) write(*,*) "iyweights are (from inside)",weights_eul_index(1:jcollect1,2,1)
!
          if (weight_sign2>-2) then
!             if (ldbg) write(*,*) "hour glass second triangle xflux"
             !
             ! hour-glass flow situation
             !
             call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
                  xcell_flux2(1:nvertex),ycell_flux2(1:nvertex),jx,jy,nreconstruction,xgno,ygno,&
                  jx_min, jx_max, jy_min, jy_max,tmp,&
                  ngauss,gauss_weights,abscissae,weights2(:,:,1),weights_eul_index2,jcollect2,jmax_segments)
             weights(jcollect1+1:jcollect1+jcollect2,:,1) = weight_sign2*weights2(1:jcollect2,:,1)
             weights_eul_index(jcollect1+1:jcollect1+jcollect2,:,1) = weights_eul_index2(1:jcollect2,:)
             jcollect1 = jcollect1+jcollect2
          end if
       end if
       jcollect(1) = jcollect1
    end if

    if (ldo_yflux) then
!       if (ldbg) write(*,*) "doing yflux"
       !
       ! constuct yflux-cell
       !
       xcell_flux = yflux_x; ycell_flux = yflux_y;

!       if (ldbg) write(*,*) "yflux before after make_flux_cell is: jx,jy",jx,jy

       call make_flux_area(jx,jy,xcell_flux,ycell_flux,xcell_flux2,ycell_flux2,&
            lzero_flux,nvertex,weight_sign,weight_sign2)
       
       if (lzero_flux) then
          jcollect(2)=0
!          if (ldbg) write(*,*) "zero yflux"
       else
!          if (ldbg) write(*,*) "nvertex",nvertex
!          if (ldbg) write(*,*) "yflux area after make_flux_cell is: jx,jy",jx,jy
!          do j=1,nvertex
!             if (ldbg) write(*,*) xcell_flux(j),ycell_flux(j)
!          end do
          call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
               xcell_flux(1:nvertex),ycell_flux(1:nvertex),jx,jy,nreconstruction,xgno,ygno,&
               jx_min, jx_max, jy_min, jy_max,tmp,&
               ngauss,gauss_weights,abscissae,weights(:,:,2),weights_eul_index(:,:,2),&
               jcollect1,jmax_segments)
          weights(1:jcollect1,:,2) = weight_sign*weights(1:jcollect1,:,2)
!          if (ldbg) write(*,*) "jcollect1 is ",jcollect1
!          if (ldbg) write(*,*) "weights are (from inside)",weights(1:jcollect1,1,2)
!          if (ldbg) write(*,*) "ixweights are (from inside)",weights_eul_index(1:jcollect1,1,2)
!          if (ldbg) write(*,*) "iyweights are (from inside)",weights_eul_index(1:jcollect1,2,2)
!          if (ldbg) write(*,*) "compute weights cell"
          if (weight_sign2>-2) then
!             if (ldbg) write(*,*) "hour glass y-flux"
             !
             ! hour-glass flow situation
             !
             call compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
                  xcell_flux2(1:nvertex),ycell_flux2(1:nvertex),jx,jy,nreconstruction,xgno,ygno,&
                  jx_min, jx_max, jy_min, jy_max,tmp,&
                  ngauss,gauss_weights,abscissae,weights2(:,:,2),weights_eul_index2,jcollect2,jmax_segments)


!             if (ldbg) write(*,*) "jcollect2 is ",jcollect2
!             if (ldbg) write(*,*) "weights2 are (from inside)",weights2(1:jcollect2,1,2)
!             if (ldbg) write(*,*) "ixweights2 are (from inside)",weights_eul_index2(1:jcollect2,1)
!             if (ldbg) write(*,*) "iyweights2 are (from inside)",weights_eul_index2(1:jcollect2,2)

             weights(jcollect1+1:jcollect1+jcollect2,:,2) = weight_sign2*weights2(1:jcollect2,:,2)
             weights_eul_index(jcollect1+1:jcollect1+jcollect2,:,2) = weights_eul_index2(1:jcollect2,:)
             jcollect1 = jcollect1+jcollect2
          end if
          jcollect(2) = jcollect1
       end if
    end if
  end subroutine compute_weights_flux_cell

  subroutine make_flux_area(jx,jy,xcell_flux,ycell_flux,xcell_flux2,ycell_flux2,&
       lzero_flux,nvertex,weight_sign,weight_sign2)
    implicit none
    integer (kind=int_kind)              , intent(in   ):: jx,jy
    real (kind=real_kind), dimension(0:), intent(inout):: xcell_flux   ,ycell_flux !dimension(0:5)
    real (kind=real_kind), dimension(0:), intent(  out):: xcell_flux2  ,ycell_flux2!dimension(0:5)
    logical                , intent(out) :: lzero_flux
    integer (kind=int_kind), intent(out):: nvertex
    real (kind=real_kind), intent(out) :: weight_sign,weight_sign2
    !
    ! local workspace
    !
    integer (kind=int_kind) :: intersect
    real (kind=real_kind), dimension(0:5) :: xcell_flux_tmp, ycell_flux_tmp
    real (kind=real_kind)                 :: xcross, ycross
    integer (kind=int_kind) :: isLeft1, isLeft2
    integer (kind=int_kind) :: j!dbg

    !
    !
    ! Graphical illustration of flux_cell arrays
    !
    !
    !              (3)-----------------------------
    !              /|                             |  
    !             / |                             |  
    !            /  |                             |  
    !          (2)  |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !    x-flux |   |                             |  
    !           |   |                             |  
    !           |   |                             |  
    !           |  (4)----------------------------|
    !           |  /
    !           | /
    !           |/  
    !          (1)
    !
    !
    !               -------------------------------
    !               |                             |  
    !               |                             |  
    !               |                             |  
    !               |                             |  
    !               |                             |  
    !               |                             |  
    !    y-flux     |                             |  
    !               |                             |  
    !               |                             |  
    !              (3)---------------------------(4)
    !              /                             /
    !             /                             /
    !            /                             /
    !          (2)----------------------------(1)
    !
    weight_sign  = -999999.99
    weight_sign2 = -999999.99

    !
    ! check if point 1 is to the left (isLeft1==1) or right (isLeft1==-1) of flux cell
    ! side (line from point 3 to 4)
    !
    isLeft1 = isLeft(xcell_flux(4),ycell_flux(4),xcell_flux(3),ycell_flux(3),&
         xcell_flux(1),ycell_flux(1))
    !
    ! do the same for point 2
    !
    isLeft2 = isLeft(xcell_flux(4),ycell_flux(4),xcell_flux(3),ycell_flux(3),&
         xcell_flux(2),ycell_flux(2))

!    if (ldbg) write(*,*) "isLeft1,isLeft2",isLeft1,isLeft2


    lzero_flux = (isLeft1==1000.and.isLeft2==1000)

    if (lzero_flux) then
       !
       ! flux area is empty (point 1 and 2 on line that connects point 3 and 4)
       !
       xcell_flux=-999999.99; ycell_flux=-999999.99; nvertex=-1
    else
       if (isLeft1==1000) then
!          if (ldbg) write(*,*) "trajec1 zero case"
          !
          ! flux-area a triangle; point 1 lies on flux-side (point 3 to 4)
          !
          ! 3 and 4 define y-flux side (similar for x-flux; rotate picture -90 degrees
          !
          !
          !                       2
          !                      / \
          ! 2\                  /   \
          ! | \         or     /     \
          ! |  \              /       \
          ! 3===1==4         3======4--1
          !
          ! In this example isLeft2=-1
          !
          !
          nvertex = 3          
          if (isLeft2>0) then
             !
             ! point 2 is to the left/under flux side (3 to 4)
             !
             weight_sign = 1.0
          else
             call reverse(xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),nvertex)
             weight_sign = -1.0
          end if
       else if (isLeft2==1000) then
          !
          ! same as previous if-statement but for point 2
          !
 !         if (ldbg) write(*,*) "trajec2 zero case"
          !
          ! flux-area a triangle
          !
          nvertex = 3
          xcell_flux(3) = xcell_flux(4); ycell_flux(3) = ycell_flux(4);
          if (isLeft1>0) then
             weight_sign = 1.0
          else
             call reverse(xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),nvertex)
             weight_sign = -1.0
          end if
       else if (ABS(isLeft1+isLeft2)==2) then          
          !
          ! Both point 1 and 2 are to the left or right of flux side (point 3 to 4)
          ! flux-area is a simply connected non-convex polygon
          !
!          if (ldbg) write(*,*) "flux-area is a quadrilateral"
          nvertex = 4
          if (isLeft1==1) then
             weight_sign=1.0
          else
             call reverse(xcell_flux,ycell_flux,nvertex)
             weight_sign  = -1.0
          end if
       else if (isLeft1+isLeft2==0) then
!          if (ldbg) write(*,*) "complex case"
          xcell_flux_tmp = xcell_flux; ycell_flux_tmp = ycell_flux 
          !
          ! hour-glass flow situation or non-simple (possibly slef-intersecting) polygon
          !
          ! There are 6 possible cases
          ! ++++++++++++++++++++++++++
          !
          ! There are 3 isLeft==1 cases:
          ! 
          !
          !                           2
          !                          /  \
          !                         /    \
          !    2                   /      \
          !   / \                 /        \            2
          !  /   \               /          \          \ \
          ! 3=====x=====4       3==========4 x          x 3==========4
          !        \   /                    \ \          \          /
          !         \ /                      \1           \        /
          !          1                                     \      /
          !                                                 \    /
          !                                                  \  /
          !                                                    1
          !
          ! intersect == 0       intersect == 1        intersect == -1
          ! clockwise            counter clockwise     clockwise
          ! (1 to x to 4)        (1 to x to 4)         (1 to x to 4)
          ! weight_sign==1       weight_sign==-1       weight_sign==1
          !
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          !
          ! There are 3 isLeft==-1 cases:
          ! 
          !
          !                                                    1
          !                                                   /  \
          !                                                  /    \
          !          1                                      /      \ 
          !         / \                      1             /        \
          !        /   \                    / /           /          \
          ! 3=====x=====4       3==========4  x          x 3==========4
          !  \   /               \           /           //
          !   \ /                 \         /            2
          !    2                   \       /              
          !                         \     /               
          !                          \   /                
          !                           \ /                 
          !                            2
          !                           
          !       
          ! intersect == 0       intersect == 1        intersect == -1
          ! counter clockwise    clockwise             counter clockwise
          ! (1 to x to 4)        (1 to x to 4)         (1 to x to 4)
          ! weight_sign ==-1     weight_sign ==1       weight_sign ==-1
          ! weight_sign2== 1     weight_sign2==1       weight_sign2==-1
          !

          call line_intersect(xcell_flux(1:4),ycell_flux(1:4),xcross,ycross,&
               intersect)
!          if (ldbg) write(*,*) "intersect",intersect
          nvertex = 3
          xcell_flux(1) = xcell_flux_tmp(1); ycell_flux(1) = ycell_flux_tmp(1);
          xcell_flux(2) = xcross            ; ycell_flux(2) = ycross;
          xcell_flux(3) = xcell_flux_tmp(4); ycell_flux(3) = ycell_flux_tmp(4);
          
          weight_sign = 1.0
          if ((isLeft1>0.and.intersect==1).or.((isLeft1<0.and.intersect.NE.1))) then
             call reverse(xcell_flux(0:nvertex+1),ycell_flux(0:nvertex+1),nvertex)
             weight_sign = -1.0
          end if
                             
!          do j=1,nvertex
!             if (ldbg) write(*,*) xcell_flux(j),ycell_flux(j)
!          end do
          
          
          xcell_flux2(1) = xcross            ; ycell_flux2(1) = ycross;
          xcell_flux2(2) = xcell_flux_tmp(2) ; ycell_flux2(2) = ycell_flux_tmp(2);
          xcell_flux2(3) = xcell_flux_tmp(3) ; ycell_flux2(3) = ycell_flux_tmp(3);
          
          if (intersect==0) then
             weight_sign2 = -weight_sign
          else
             weight_sign2 =  weight_sign
          end if
          
          if ((isLeft1>0.and.intersect>-1).or.((isLeft1<0.and.intersect==-1))) &
               call reverse(xcell_flux2(0:nvertex+1),ycell_flux2(0:nvertex+1),nvertex)
          
       end if
    end if
    
  end subroutine make_flux_area

  subroutine orient(x,y,nvertex,weight_sign)
    use fvm_line_integrals_mod, only: area
    implicit none
    integer (kind=int_kind)                      , intent(in):: nvertex
    real (kind=real_kind), dimension(0:), intent(inout):: x,y !dimension(0:nvertex+1)
    real (kind=real_kind)                        , intent(out  ):: weight_sign
    !
    real (kind=real_kind), dimension(0:nvertex+1) :: xtmp,ytmp 

    if (area(x(1:nvertex),y(1:nvertex),nvertex)<0) then
       xtmp(1:nvertex)=x(1:nvertex); ytmp(1:nvertex)=y(1:nvertex)

       x(1:nvertex) = xtmp(nvertex:1:-1); y(1:nvertex) = ytmp(nvertex:1:-1);
       x(0        ) = x   (nvertex     ); y(0        ) = y   (nvertex     );
       x(nvertex+1) = x   (1           ); y(nvertex+1) = y   (1           );
       weight_sign  = -1.0
    else
       xtmp(1:nvertex)=x(1:nvertex); ytmp(1:nvertex)=y(1:nvertex)
       weight_sign =  1.0
    end if
  end subroutine orient


  subroutine reverse(x,y,nvertex)
    implicit none
    integer (kind=int_kind)           , intent(in):: nvertex
    real (kind=real_kind), dimension(0:), intent(inout):: x,y !dimension(0:nvertex+1)
    !
    real (kind=real_kind), dimension(0:nvertex+1) :: xtmp,ytmp

    xtmp(1:nvertex)=x(1:nvertex); ytmp(1:nvertex)=y(1:nvertex)

    x(1:nvertex) = xtmp(nvertex:1:-1); y(1:nvertex) = ytmp(nvertex:1:-1);
    x(0        ) = x   (nvertex     ); y(0        ) = y   (nvertex     );
    x(nvertex+1) = x   (1           ); y(nvertex+1) = y   (1           );
  end subroutine reverse


subroutine getdep_cellboundariesxyvec_xflux(xcell,ycell,jx,jy,acart,dcart)
use coordinate_systems_mod, only : cartesian2D_t
  implicit none
  real (kind=real_kind), dimension(0:), intent(out)       :: xcell,ycell !dimension(0:5)
  integer (kind=int_kind), intent(in)                      :: jx, jy
  type (cartesian2D_t), intent(in)                         :: dcart(-1:,-1:),acart(-1:,-1:) !(-1:nc+3,-1:nc+3)

  xcell(1) = dcart(jx,jy  )%x ; ycell(1) = dcart(jx,jy  )%y
  xcell(2) = dcart(jx,jy+1)%x ; ycell(2) = dcart(jx,jy+1)%y
  xcell(3) = acart(jx,jy+1)%x ; ycell(3) = acart(jx,jy+1)%y
  xcell(4) = acart(jx,jy  )%x ; ycell(4) = acart(jx,jy  )%y

  !
  ! truncate
  !
  if (ABS(xcell(1)-xcell(4))<1.0E-9) xcell(1) = xcell(4)
  if (ABS(xcell(2)-xcell(3))<1.0E-9) xcell(2) = xcell(3)
  if (ABS(ycell(1)-ycell(4))<1.0E-9) ycell(1) = ycell(4)
  if (ABS(ycell(2)-ycell(3))<1.0E-9) ycell(2) = ycell(3)
  !
  !
  ! are these lines necessary?
  !
  xcell(5) = xcell(1)         ; ycell(5) = ycell(1)          
  xcell(0) = xcell(4)         ; ycell(0) = ycell(4)
end subroutine getdep_cellboundariesxyvec_xflux

subroutine getdep_cellboundariesxyvec_yflux(xcell,ycell,jx,jy,acart,dcart)
use coordinate_systems_mod, only : cartesian2D_t
  implicit none
  real (kind=real_kind), dimension(0:), intent(out)       :: xcell,ycell !dimension(0:5)
  integer (kind=int_kind), intent(in)                      :: jx, jy
  type (cartesian2D_t), intent(in)                         :: dcart(-1:,-1:), acart(-1:,-1:)!(-1:nc+3,-1:nc+3)

  xcell(1) = dcart(jx+1,jy)%x ; ycell(1) = dcart(jx+1,jy)%y
  xcell(2) = dcart(jx  ,jy)%x ; ycell(2) = dcart(jx  ,jy)%y
  xcell(3) = acart(jx  ,jy)%x ; ycell(3) = acart(jx  ,jy)%y
  xcell(4) = acart(jx+1,jy)%x ; ycell(4) = acart(jx+1,jy)%y

  !
  ! truncate
  !
  if (ABS(xcell(1)-xcell(4))<1.0E-9) xcell(1) = xcell(4)
  if (ABS(xcell(2)-xcell(3))<1.0E-9) xcell(2) = xcell(3)
  if (ABS(ycell(1)-ycell(4))<1.0E-9) ycell(1) = ycell(4)
  if (ABS(ycell(2)-ycell(3))<1.0E-9) ycell(2) = ycell(3)
  !
  ! are these lines necessary?
  !
  xcell(5) = xcell(1)         ; ycell(5) = ycell(1)          
  xcell(0) = xcell(4)         ; ycell(0) = ycell(4)
end subroutine getdep_cellboundariesxyvec_yflux


integer(kind=int_kind) function isLeft(ax,ay,bx,by,cx,cy)
  implicit none
  real (kind=real_kind), INTENT(IN) :: ax,ay,bx,by,cx,cy
  
  if (((bx - ax)*(cy - ay) - (by - ay)*(cx - ax)) > tiny) then
     isLeft = 1
  else if (((bx - ax)*(cy - ay) - (by - ay)*(cx - ax)) < -tiny) then
     isLeft = -1
  else
     isLeft = 1000
  end if
end function isLeft



!
! this algorithm is from Paul Bourke:
!
! http://paulbourke.net/geometry/pointlineplane/pdb.c
!
! returns point of intersection between the straight lines that connect:
!
! Point 1 (x(1),y(1)) and point 2 (x(2),y(2)) 
!
!                     and
!
! point 3 (x(3),y(3)) and point 4 (x(4),y(4))
!
! Intersect is +/-1 if intersection is to the right/left of line segment
! that goes from point 3 to 4. If intersect is 0 the intersection is on
! the line segment.
!
subroutine line_intersect(x,y,xintersect,yintersect,intersect)
  implicit none
  real (kind=real_kind), DIMENSION(:), INTENT(IN ) :: x,y !dimension(4)
  real (kind=real_kind), INTENT(OUT) :: xintersect,yintersect
  integer (kind=int_kind), intent(out) :: intersect

  real (kind=real_kind) :: mua,mub
  real (kind=real_kind) :: denom,numera,numerb

   denom  = (y(4)-y(3)) * (x(2)-x(1)) - (x(4)-x(3)) * (y(2)-y(1));
   numera = (x(4)-x(3)) * (y(1)-y(3)) - (y(4)-y(3)) * (x(1)-x(3));
   numerb = (x(2)-x(1)) * (y(1)-y(3)) - (y(2)-y(1)) * (x(1)-x(3));

   intersect = -100000
   ! Are the line coincident?
   if (ABS(numera) < tiny .AND. ABS(numerb) < tiny .AND. ABS(denom) < tiny) then
      write(*,*) "line coincident - should not happen!"
      stop
   else if (ABS(denom) < tiny) then
      ! Are the line parallel?
      write(*,*) "lines parallel  - should not happen!"
      stop
   else
      mua = numera / denom
      mub = numerb / denom
      xintersect = x(1) + mua * (x(2) - x(1));
      yintersect = y(1) + mua * (y(2) - y(1));

      ! Is the intersection along the segments
      if (mub<0) then
         intersect=-1
      else if (mub>1.0) then
         intersect= 1
      else
         intersect = 0
      end if
   end if
 end subroutine line_intersect

subroutine debugging_print_cells(jx,jy,acart,dcart)
  use coordinate_systems_mod, only : cartesian2D_t
  implicit none
  integer (kind=int_kind), intent(in) :: jx,jy
  type (cartesian2D_t), intent(in) :: dcart(-1:nc+3,-1:nc+3)
  type (cartesian2D_t), intent(in) :: acart(-1:nc+3,-1:nc+3)

  write(*,*) "jx,jy",jx,jy
  write(*,*) "------------------------------"
  write(*,*) "dep cell"
  write(*,*) dcart(jx,jy  )%x,dcart(jx,jy)%y
  write(*,*) dcart(jx,jy+1)%x,dcart(jx,jy+1)%y
  write(*,*) dcart(jx+1,jy+1)%x,dcart(jx+1,jy+1)%y
  write(*,*) dcart(jx+1,jy)%x,dcart(jx+1,jy)%y
  write(*,*) dcart(jx,jy)%x,dcart(jx,jy)%y
  write(*,*) "-----------------"
  write(*,*) "jx,jy",jx,jy
  write(*,*) "------------------------------"
  write(*,*) "arr cell"
  write(*,*) acart(jx,jy  )%x,acart(jx,jy)%y
  write(*,*) acart(jx,jy+1)%x,acart(jx,jy+1)%y
  write(*,*) acart(jx+1,jy+1)%x,acart(jx+1,jy+1)%y
  write(*,*) acart(jx+1,jy)%x,acart(jx+1,jy)%y
  write(*,*) acart(jx,jy)%x,acart(jx,jy)%y
  write(*,*) "-----------------"


end subroutine debugging_print_cells
end module fvm_line_integrals_flux_mod
