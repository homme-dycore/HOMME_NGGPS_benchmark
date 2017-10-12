#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!MODULE FVM_LINE_INTEGRALS_MOD--------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! Given four vertices of a simply connected cell (x(i),y(i)), i=1,4 ordered in a    !
! counter clockwise manner, compute coefficients for line integral.                 !               
! This module contains everything  to do that and is base on the weights computation!
! from Peter Lauritzens code, here adapted for HOMME                                ! 
!                                                                                   !
!-----------------------------------------------------------------------------------!
module fvm_line_integrals_mod

  use kinds, only               : int_kind, real_kind
  use dimensions_mod, only      : nc, nhe, ngpc
  use dimensions_mod, only: nelem, nelemd, nelemdmax, nlev, ne, nc, nhe, nlev, ntrac, np, ntrac_d,ns, nhr, nhc
  use dimensions_mod, only: ie_dbg, pr_dbg
  use time_mod, only : timelevel_t
  use element_mod, only : element_t, timelevels
  use fvm_control_volume_mod, only: fvm_struct
  use hybrid_mod, only : hybrid_t
  use edge_mod, only : initghostbufferTR, freeghostbuffertr, &
       ghostVpack, ghostVunpack,  initEdgebuffer
  use edgetype_mod, only : ghostbuffertr_t, edgebuffer_t
  use perf_mod, only : t_startf, t_stopf ! EXTERNAL
  use perf_utils, only: t_detail_low, t_detail_medium, t_detail_high, t_detail_max  ! EXTERNAL

  implicit none
  private
  type (EdgeBuffer_t)                         :: edgeveloc

  real (kind=real_kind),parameter, public   :: bignum = 1.0D20
  real (kind=real_kind),parameter, public   :: tiny   = 1.0D-12
  real (kind=real_kind),parameter           :: fuzzy_width = 10.0*tiny
  ! turn on/off EOC (Enforcement of Consistency) -> Erath et al. MWR, 2013
  logical                                   :: EOC=.TRUE.
  integer, private, parameter :: num_weights      = 10*(nc+2*nhe)*(nc+2*nhe)
  
  logical, public :: ldbg=.false.!dbg xxx

  ! namelist variables for testing
  integer, public, parameter            :: IDEAL_TEST_OFF = 0
  integer, public, parameter            :: IDEAL_TEST_ANALYTICAL_DEPARTURE = 1
  integer, public, parameter            :: IDEAL_TEST_ANALYTICAL_WINDS = 2
  integer, public                       :: fvm_ideal_test = IDEAL_TEST_OFF
  integer, public, parameter            :: IDEAL_TEST_BOOMERANG = 1
  integer, public, parameter            :: IDEAL_TEST_SOLIDBODY = 2
  integer, public                       :: fvm_test_type = IDEAL_TEST_BOOMERANG








  public :: compute_weights, compute_weights_cell, gauss_points, getdep_cellboundariesxyvec
  public :: compute_slope,y_cross_eul_lon,x_cross_eul_lat,area, truncate_vertex
  public :: cslam_runairdensity
  public ::  fvm_mcgregor, fvm_mcgregordss, fvm_rkdss
  public :: fvm_mesh_dep ! for fvm_debugging

contains

  !
  !**************************************************************************************
  !
  !
  ! Lagrangian cslam subroutines
  !
  !
  !**************************************************************************************
  !


  ! use this subroutine for benchmark tests, couple airdensity with tracer concentration
  subroutine cslam_runairdensity(elem,fvm,hybrid,deriv,dt_fvm,tl,nets,nete,p_top)
    ! ---------------------------------------------------------------------------------
    use fvm_control_volume_mod, only: n0_fvm, np1_fvm
    ! ---------------------------------------------------------------------------------  
    use fvm_reconstruction_mod, only: reconstruction_gradient
    ! ---------------------------------------------------------------------------------
    use derivative_mod, only : derivative_t
    ! ---------------------------------------------------------------------------------
    use edge_mod, only :  ghostVpack2d_level, ghostVunpack2d_level,initghostbufferTR,freeghostbuffertr
    use edgetype_mod, only : ghostBuffertr_t
    use fvm_mod, only : fill_halo_fvm
    use bndry_mod, only : ghost_exchangevfull
    use bndry_mod, only: ghost_exchangeV                     
    implicit none
    type (element_t), intent(inout)                :: elem(:)
    type (fvm_struct), intent(inout)             :: fvm(:)
    type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
    integer, intent(in)                         :: nets  ! starting thread element number (private)
    integer, intent(in)                         :: nete  ! ending thread element number   (private)
    real (kind=real_kind), intent(in)           :: dt_fvm
    real (kind=real_kind), intent(in)           :: p_top
    
    integer                                     :: i,j,k,ie,itr, jx, jy, jdx, jdy, h, ntmp
    type (TimeLevel_t)                          :: tl              ! time level struct
    type (derivative_t)                         :: deriv           ! derivative struct
    
    real (kind=real_kind)   , dimension(num_weights,6)  :: weights_all
    integer (kind=int_kind),  dimension(num_weights,2)  :: weights_eul_index_all
    integer (kind=int_kind),  dimension(num_weights,2)  :: weights_lgr_index_all
    integer (kind=int_kind)                                          :: jall
    
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,6)      :: recons
    
    real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer0 
    
    real (kind=real_kind), dimension(1-nhc:nc+nhc,1-nhc:nc+nhc)        :: tracer_air0   
    real (kind=real_kind), dimension(1:nc,1:nc)                        :: tracer1, tracer_air1 
    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,6)      :: recons_air   

    type (ghostBuffertr_t)                      :: buflatlon

    call initghostbufferTR(buflatlon,nlev,2,2,nc+1)    ! use the tracer entry 2 for lat lon
    call t_startf('cslam scheme', t_detail_high) 

    do ie=nets, nete
       do k=1,nlev
          call fvm_mesh_dep(elem(ie),deriv,fvm(ie),dt_fvm,tl,k)
       end do
    end do

    !
    ! is this boundary exchange necessary? only of .EOC.=.TRUE. in fvm_lineintegrals_mod
    !    
    do ie=nets,nete
       ! if changing to nhe>1 the 3rd argument has to be changed to nhe+1
       call ghostVpack2d_level(buflatlon,fvm(ie)%dsphere(:,:,:)%lat,1,2, nc+1,nlev,elem(ie)%desc) !kptr = 1 for lat
       call ghostVpack2d_level(buflatlon,fvm(ie)%dsphere(:,:,:)%lon,2,2, nc+1,nlev,elem(ie)%desc) !kptr =2 for lon
    end do
    !-----------------------------------------------------------------------------------! 
    call ghost_exchangeV(hybrid,buflatlon,2,nc+1,2)
    !-----------------------------------------------------------------------------------!  
    do ie=nets,nete
       call ghostVunpack2d_level(buflatlon,fvm(ie)%dsphere(:,:,:)%lat,1,2, nc+1,nlev,elem(ie)%desc)
       call ghostVunpack2d_level(buflatlon,fvm(ie)%dsphere(:,:,:)%lon,2,2, nc+1,nlev,elem(ie)%desc)
       fvm(ie)%dsphere(:,:,:)%r=1.0D0  !!! RADIUS IS ASSUMED TO BE 1.0DO !!!!       
    end do

    !
    ! fill halo for dp_fvm and c --- call fill_halo_fvm instead of code below!
    !
    call fill_halo_fvm(elem,fvm,hybrid,nets,nete,n0_fvm)

    do ie=nets, nete
       do k=1,nlev
          !-Departure fvm Meshes, initialization done                                                               
          call compute_weights(fvm(ie),6,weights_all,weights_eul_index_all, &
               weights_lgr_index_all,k,jall)     
          tracer_air0=fvm(ie)%dp_fvm(:,:,k,n0_fvm)     
          call reconstruction_gradient(tracer_air0, fvm(ie),recons_air,6,.false.)
          tracer_air1=0.0D0   
          call cslam_remap(tracer_air0,tracer_air1,weights_all, recons_air, &
               fvm(ie)%spherecentroid, weights_eul_index_all, weights_lgr_index_all, jall,0)
          ! finish scheme
          do j=1,nc
             do i=1,nc
                tracer_air1(i,j)=tracer_air1(i,j)/fvm(ie)%area_sphere(i,j)
                fvm(ie)%dp_fvm(i,j,k,np1_fvm)=tracer_air1(i,j)
             end do
          end do
          !loop through all tracers
          do itr=1,ntrac
             tracer0=fvm(ie)%c(:,:,k,itr,n0_fvm)
             call reconstruction_gradient(tracer0, fvm(ie),recons,6,.true.)

             tracer1=0.0D0   
             call cslam_remap_air(tracer0,tracer1,tracer_air0,weights_all, recons,recons_air,&
                  fvm(ie)%spherecentroid,weights_eul_index_all, weights_lgr_index_all, jall)  

             ! finish scheme
             do j=1,nc
                do i=1,nc
                   fvm(ie)%c(i,j,k,itr,np1_fvm)= &
                        (tracer1(i,j)/fvm(ie)%area_sphere(i,j))/tracer_air1(i,j)
                   !             fvm(ie)%c(i,j,k,itr,np1_fvm)=tracer1(i,j)/fvm(ie)%area_sphere(i,j)
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

       !note write tl%np1 in buffer                                                                 

!       call ghostVpack(cellghostbuf, fvm(ie)%dp_fvm(:,:,:,tl%np1),nhc,nc,nlev,1,    0,   elem(ie)%desc)
!       call ghostVpack(cellghostbuf, fvm(ie)%c(:,:,:,:,tl%np1),   nhc,nc,nlev,ntrac,1,elem(ie)%desc)
    end do

    !
    ! advance fvm time-levels
    !
    ntmp     = np1_fvm
    np1_fvm  = n0_fvm
    n0_fvm   = ntmp



    call t_stopf('cslam scheme', t_detail_high)
!    call t_startf('FVM Communication')
!    call ghost_exchangeV(hybrid,cellghostbuf,nhc,nc,ntrac+1)
!    call t_stopf('FVM Communication')
    !-----------------------------------------------------------------------------------!                         
!    call t_startf('FVM Unpack')
!    do ie=nets,nete
!       call ghostVunpack(cellghostbuf, fvm(ie)%dp_fvm(:,:,:,tl%np1), nhc, nc,nlev,1,    0,   elem(ie)%desc)
!       call ghostVunpack(cellghostbuf, fvm(ie)%c(:,:,:,:,tl%np1),    nhc, nc,nlev,ntrac,1,elem(ie)%desc)
!    enddo
!    call t_stopf('FVM Unpack')
    call freeghostbuffertr(buflatlon)

    
  end subroutine cslam_runairdensity


  subroutine cslam_remap(tracer0,tracer1,weights_all,recons,centroid, &
       weights_eul_index_all, weights_lgr_index_all, jall, iflux)
    
    real (kind=real_kind),   intent(in)        :: tracer0(1-nhc:nc+nhc,1-nhc:nc+nhc)
    integer (kind=int_kind), intent(in)        :: iflux !to accomodate flux call
    real (kind=real_kind),   intent(inout)     :: tracer1(1:nc+iflux,1:nc+iflux)
    integer (kind=int_kind), intent(in)        :: jall  
    real (kind=real_kind),   intent(in)        :: recons(1-nhe:nc+nhe,1-nhe:nc+nhe,6)
    real (kind=real_kind),   intent(in)        :: centroid(1-nhe:nc+nhe,1-nhe:nc+nhe,5)
    real (kind=real_kind)   , dimension(num_weights,6), intent(in)  :: weights_all
    integer (kind=int_kind),  dimension(num_weights,2), intent(in)  :: weights_eul_index_all
    integer (kind=int_kind),  dimension(num_weights,2), intent(in)  :: weights_lgr_index_all

 
    integer                                     :: h, jx, jy, jdx, jdy
    
    tracer1 = 0.0D0
    do h=1,jall
       jx  = weights_lgr_index_all(h,1)
       jy  = weights_lgr_index_all(h,2)
       jdx = weights_eul_index_all(h,1)
       jdy = weights_eul_index_all(h,2)
       
       tracer1(jx,jy) = tracer1(jx,jy)+weights_all(h,1)*(&
            ! all constant terms 
            tracer0(jdx,jdy) - recons(jdx,jdy,2)*centroid(jdx,jdy,1) &
            - recons(jdx,jdy,3)*centroid(jdx,jdy,2) &
            + recons(jdx,jdy,4)*(2.0D0*centroid(jdx,jdy,1)**2 -centroid(jdx,jdy,3)) &
            + recons(jdx,jdy,5)*(2.0D0*centroid(jdx,jdy,2)**2 -centroid(jdx,jdy,4)) &
            + recons(jdx,jdy,6)*(2.0D0*centroid(jdx,jdy,1)*centroid(jdx,jdy,2)-centroid(jdx,jdy,5))) + &
            ! linear terms
            weights_all(h,2)*&
            (recons(jdx,jdy,2) - recons(jdx,jdy,4)*2.0D0*centroid(jdx,jdy,1) &
            - recons(jdx,jdy,6)*centroid(jdx,jdy,2)) + &
            weights_all(h,3)*&
            (recons(jdx,jdy,3) - recons(jdx,jdy,5)*2.0D0*centroid(jdx,jdy,2) &
            - recons(jdx,jdy,6)*centroid(jdx,jdy,1)) + &
            ! quadratic terms
            weights_all(h,4)*recons(jdx,jdy,4)+&
            weights_all(h,5)*recons(jdx,jdy,5)+&
            weights_all(h,6)*recons(jdx,jdy,6)
    end do
  end subroutine cslam_remap


  
  ! do remapping with air (i.e. conserve mass of air density * concentration),
  ! see Nair et.al 2010 in JCP: A class of deformational flow test cases for linear transport
  ! schemes on the sphere, Appendix B
  subroutine cslam_remap_air(tracer0,tracer1,tracer_air, weights_all,recons, recons_air, centroid, &
       weights_eul_index_all, weights_lgr_index_all, jall)
    
    real (kind=real_kind), intent(in)           :: tracer0(1-nhc:nc+nhc,1-nhc:nc+nhc)
    real (kind=real_kind), intent(inout)        :: tracer1(1:nc,1:nc)
    real (kind=real_kind), intent(in)           :: tracer_air(1-nhc:nc+nhc,1-nhc:nc+nhc)
    real (kind=real_kind), intent(in)           :: recons(1-nhe:nc+nhe,1-nhe:nc+nhe,6)
    real (kind=real_kind), intent(in)           :: recons_air(1-nhe:nc+nhe,1-nhe:nc+nhe,6)
 
    integer (kind=int_kind), intent(in)         :: jall  
    real (kind=real_kind), intent(in)           :: weights_all(num_weights,6)
    real (kind=real_kind), intent(in)           :: centroid(1-nhe:nc+nhe,1-nhe:nc+nhe,5)
    integer (kind=int_kind), intent(in)         :: weights_eul_index_all(num_weights,2)
    integer (kind=int_kind), intent(in)         :: weights_lgr_index_all(num_weights,2)
 
    integer                                     :: h, jx, jy, jdx, jdy    
 
    do h=1,jall
       jx  = weights_lgr_index_all(h,1)
       jy  = weights_lgr_index_all(h,2)
       jdx = weights_eul_index_all(h,1)
       jdy = weights_eul_index_all(h,2)                     
       tracer1(jx,jy) = tracer1(jx,jy)+&
            ! air density times tracer reconstruction
            tracer_air(jdx,jdy)*(weights_all(h,1)*(&      ! 1 is for air
            ! all constant terms 
            tracer0(jdx,jdy) - recons(jdx,jdy,2)*centroid(jdx,jdy,1) &
            - recons(jdx,jdy,3)*centroid(jdx,jdy,2) &
            + recons(jdx,jdy,4)*(2.0D0*centroid(jdx,jdy,1)**2 -centroid(jdx,jdy,3)) &
            + recons(jdx,jdy,5)*(2.0D0*centroid(jdx,jdy,2)**2 -centroid(jdx,jdy,4)) &
            + recons(jdx,jdy,6)*(2.0D0*centroid(jdx,jdy,1)*centroid(jdx,jdy,2)-centroid(jdx,jdy,5))) + &
            ! linear terms
            weights_all(h,2)* &
            (recons(jdx,jdy,2)- recons(jdx,jdy,4)*2.0D0*centroid(jdx,jdy,1) &
            - recons(jdx,jdy,6)*centroid(jdx,jdy,2)) + &
            weights_all(h,3)* &
            (recons(jdx,jdy,3) - recons(jdx,jdy,5)*2.0D0*centroid(jdx,jdy,2) &
            - recons(jdx,jdy,6)*centroid(jdx,jdy,1)) + &
            ! quadratic terms
            weights_all(h,4)*recons(jdx,jdy,4) + &
            weights_all(h,5)*recons(jdx,jdy,5) + &
            weights_all(h,6)*recons(jdx,jdy,6)) + &
            
            !tracer times air reconstruction
            tracer0(jdx,jdy)*(weights_all(h,1)*(&      
            ! all constant terms 
            !       tracer_air &  this term cancels it out
            - recons_air(jdx,jdy,2)*centroid(jdx,jdy,1) - recons_air(jdx,jdy,3)*centroid(jdx,jdy,2) &
            + recons_air(jdx,jdy,4)*(2.0D0*centroid(jdx,jdy,1)**2 -centroid(jdx,jdy,3)) &
            + recons_air(jdx,jdy,5)*(2.0D0*centroid(jdx,jdy,2)**2 -centroid(jdx,jdy,4)) &
            + recons_air(jdx,jdy,6)*(2.0D0*centroid(jdx,jdy,1)*centroid(jdx,jdy,2)-centroid(jdx,jdy,5))) + &
            ! linear terms
            weights_all(h,2)* &
            (recons_air(jdx,jdy,2) - recons_air(jdx,jdy,4)*2.0D0*centroid(jdx,jdy,1) &
            - recons_air(jdx,jdy,6)*centroid(jdx,jdy,2)) + &
            weights_all(h,3)* &
            (recons_air(jdx,jdy,3) - recons_air(jdx,jdy,5)*2.0D0*centroid(jdx,jdy,2) &
            - recons_air(jdx,jdy,6)*centroid(jdx,jdy,1)) + &
            ! quadratic terms
            weights_all(h,4)*recons_air(jdx,jdy,4)+&
            weights_all(h,5)*recons_air(jdx,jdy,5)+&
            weights_all(h,6)*recons_air(jdx,jdy,6))
    end do
  end subroutine cslam_remap_air

  
  subroutine cslam_get_area(fvm,area,k)
    use fvm_control_volume_mod, only:  fvm_struct
    type (fvm_struct), intent(inout)                               :: fvm  
    real (kind=real_kind), intent(out)        :: area(1:nc,1:nc)
    integer, intent(in)                                            :: k
    integer                                                        :: h, jx, jy, jdx, jdy
    real (kind=real_kind)   , dimension(10*(nc+2*nhe)*(nc+2*nhe),6):: weights_all
    integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2):: weights_eul_index_all
    integer (kind=int_kind),  dimension(10*(nc+2*nhe)*(nc+2*nhe),2):: weights_lgr_index_all
    integer (kind=int_kind)                                        :: jall
    
    call compute_weights(fvm,6,weights_all,weights_eul_index_all, &
         weights_lgr_index_all,k,jall)
    
    area = 0.0D0
    do h=1,jall
       jx  = weights_lgr_index_all(h,1)
       jy  = weights_lgr_index_all(h,2)
       jdx = weights_eul_index_all(h,1)
       jdy = weights_eul_index_all(h,2)
       
       area(jx,jy) = area(jx,jy)+weights_all(h,1)
    end do
  end subroutine cslam_get_area
  




! ----------------------------------------------------------------------------------!
!SUBROUTINE COMPUTE_WEIGHTS-----------------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
! DESCRIPTION: computes the coefficients for line integral, set up everything to use!
!              the existing subroutines for the extended (incl. halo zone) element  !
!                                                                                   !
! CALLS: compute_weights_cell, getdep_cellboundaries                                !
! INPUT:  fvm  ... structure, see fvm_control_volume_mod.F90                    ! 
!         nreconstruction ... choose the reconstruction: 1-linear, 3-quadratic      !
!                             6-cubic                                               !
! OUTPUT: weights_all ... value of the weights, number depends of the different     !
!                         intersection/overlap areas to the Eulerian cell           !
!         weights_eul_index_all...correspoding Eulerian cell index of the overlap   !
!         weights_lgr_index_all...which arrival cell is conntected with the         !
!                                 departure cell                                    !
!         jall...number of intersections of an element                              !
!-----------------------------------------------------------------------------------
subroutine compute_weights(fvm,nreconstruction,weights_all,weights_eul_index_all, &
                                           weights_lgr_index_all,klev,jall)  
  use fvm_control_volume_mod, only:  fvm_struct                                         
  use coordinate_systems_mod,  only :  cartesian2D_t, spherical_polar_t, &
                                       cart2cubedspherexy, spherical_to_cart
  use physical_constants, only : DD_PI
  use control_mod, only : north, south, east, west, neast, nwest, seast, swest

  implicit none
  type (fvm_struct), intent(inout)                                :: fvm
  integer (kind=int_kind), intent(in)                            :: nreconstruction
  ! arrays for collecting cell data
  !
  ! dimension(10*(nc+2*nhe)*(nc+2*nhe),2)
  !
  real (kind=real_kind)  , dimension(:,:), intent(out) :: weights_all
  integer (kind=int_kind), dimension(:,:), intent(out) :: weights_eul_index_all, weights_lgr_index_all
  integer (kind=int_kind), intent(in)                       :: klev
  integer (kind=int_kind), intent(out)                      :: jall

  ! local workspace
  ! max number of line segments is:
  ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
  !

  integer (kind=int_kind)                     :: jx,jy
  integer                                     :: jx_min, jx_max, jy_min, jy_max
  integer                                     :: jx_min1, jx_max1, jy_min1, jy_max1
  integer                                     :: jx_min2, jx_max2, jy_min2, jy_max2
  logical                                     :: swap1, swap2
  
  integer (kind=int_kind)                     :: i, jtmp, k
  
  type (cartesian2D_t)                        :: dcart(-1:nc+3,-1:nc+3)       ! Cartesian coordinates 
  
  real (kind=real_kind), dimension(0:5)       :: xcell,ycell
  integer (kind=int_kind)                     :: inttmp
  real (kind=real_kind)                       :: tmp
  logical                                     :: swap
  ! for Gaussian quadrature
  real (kind=real_kind), dimension(ngpc)      :: gsweights, gspts
  ! weight-variables for individual cells
  integer (kind=int_kind) :: jmax_segments_cell
  real (kind=real_kind)   , dimension(nhe*50,nreconstruction)   :: weights_cell
  integer (kind=int_kind),  dimension(nhe*50,2)                 :: weights_eul_index_cell
  integer (kind=int_kind)                                       :: jcollect_cell
  
  integer (kind=int_kind)                    :: jallactual, jallactual_eul, ja
  real (kind=real_kind)                      :: da_cslam(1-nhe:nc+nhe,1-nhe:nc+nhe), centroid_cslam(1-nhe:nc+nhe,1-nhe:nc+nhe,5)
  real (kind=real_kind)                      :: area

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
  do jy=-1,nc+3
    do jx=-1,nc+3  
      call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
           fvm%faceno,dcart(jx,jy))  
    end do
  end do
  
  do jy=1, nc
     do jx=1, nc            
        !
        ! define departure cell
        !
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     

        call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 

        if (jcollect_cell>0) then
           weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
           weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                weights_eul_index_cell(1:jcollect_cell,:)
           weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
           weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
           jall = jall+jcollect_cell          
        endif
     end do
  end do
  jallactual=jall
  
!!!!  WEIGHTS CORRECTION FOR THE interior element/cells, i.e.
  ! xphl need to compute all overlaps that span the Eulerian halo cells
  ! xphl 
  ! xphl
  ! 
  ! xphl Erath C., P.H. Lauritzen, and H.M Tufo. 2013: On mass-conservation in high-order high-resolution 
  ! xphl rigorous remapping schemes on the sphere. Mon. Wea. Rev.
  !
  if (EOC) then
     do jy=jy_min-1, 0
        do jx=jx_min-1, jx_max  
           if ((fvm%cubeboundary == swest) .and. (jx<1) .and. (jy<1)) then
              !
              ! xphl no cells in "south-west" halo
              !
              cycle
           endif
           if ((fvm%cubeboundary == seast) .and. (jx>nc) .and. (jy<1)) then
              cycle
           endif
           call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)  
           call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
                fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
                tmp,ngpc,gsweights,gspts,&
                weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
           if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
           endif
        end do
     end do
     
     do jy=nc+1, jy_max
        do jx=jx_min-1, jx_max      
           if ((fvm%cubeboundary == nwest) .and. (jx<1) .and. (jy>nc)) then
              cycle
           endif
           if ((fvm%cubeboundary == neast) .and. (jx>nc) .and. (jy>nc)) then
              cycle
           endif
           call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)  
           call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
                fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
                tmp,ngpc,gsweights,gspts,&
                weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
           if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
           endif
        end do
     end do
     do jx=jx_min-1, 0
        do jy=1, nc             
           call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)  
           call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
                fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
                tmp,ngpc,gsweights,gspts,&
                weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
           if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
           endif
        end do
     end do
     do jx=nc+1, jx_max
        do jy=1, nc             
           call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)  
           call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
                fvm%acartx,fvm%acarty,jx_min, jx_max, jy_min, jy_max, &
                tmp,ngpc,gsweights,gspts,&
                weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell) 
           if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
           endif
        end do
     end do
     !
     ! xphl done weight correction
     !
     da_cslam=0.0D0
     centroid_cslam=0.0D0
     do ja=1,jall-1
        jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
        da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
        centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
     end do
     
     jall=jallactual
     jallactual_eul=jall
  endif !end EOC
    
  !WEST SIDE
  if (fvm%cubeboundary == west) then
     ! calculate xy Cartesian on the cube of departure points on the corresponding face
     do jx=-1,2      
        do jy=-1,nc+3
           call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
                fvm%nbrsface(west),dcart(jx,jy))                  
        end do
     end do
     jx=1
     do jy=1,nc
        !       call getdep_cellboundaries(xcell,ycell,jx,jy,fvm%nbrsface(west),fvm%dsphere) 
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
        
        if(swap1) then  !flip orientation
           call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                tmp,ngpc,gsweights,gspts,&
                weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
           do i=1,jcollect_cell
              !
              ! xphl swap i and j indices for ovelap areas on panel to the west - what panel?
              !
              inttmp=weights_eul_index_cell(i,1)
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
              weights_eul_index_cell(i,2)=inttmp
           end do
        else  
           call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
                fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                tmp,ngpc,gsweights,gspts,&
                weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        end if
        !I have to correct the number - xphl????
        if (fvm%faceno==5) then
           do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
           end do
        end if
        if (fvm%faceno==6) then
           do i=1,jcollect_cell
              weights_eul_index_cell(i,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2)-nhe-nhe+1
           end do
        end if
        if (jcollect_cell>0) then
           weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
           weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                weights_eul_index_cell(1:jcollect_cell,:)
           weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
           weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
           jall = jall+jcollect_cell          
        endif
     end do
     
     ! for Eulerian Correction  (Erath et al., MWR, 2013)
     if(EOC) then
        jallactual=jall
        !
        do jx=-1,1      
          do jy=-1,nc+2
            if ((jx==1) .and. (jy>0) .and. (jy<nc+1)) then 
              cycle
            endif
            !       call getdep_cellboundaries(xcell,ycell,jx,jy,fvm%nbrsface(west),fvm%dsphere) 
            call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          
            if(swap1) then  !flip orientation
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                   fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
              do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1)
                weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
                weights_eul_index_cell(i,2)=inttmp
              end do
            else  
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min1,jy_min1,nreconstruction,&
                   fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            end if
            !I have to correct the number
            if (fvm%faceno==5) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
              end do
            end if
            if (fvm%faceno==6) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2)-nhe-nhe+1
              end do
            end if
            if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
            endif
          end do
        end do
      
        do ja=jallactual_eul,jall-1
          jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
          da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
          centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
        end do
        ! do not save weights only used for EOC    
        jall=jallactual
      endif
    endif
    !EAST SIDE
    if (fvm%cubeboundary == east) then
      ! calculate xy Cartesian on the cube of departure points on the corresponding face 
      do jx=nc,nc+3 
        do jy=-1,nc+3
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(east),dcart(jx,jy))                  
        end do
      end do
      jx=nc
      do jy=1,nc
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
        if(swap1) then !flip orientation
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
               fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          do i=1,jcollect_cell
            inttmp=weights_eul_index_cell(i,1)
            weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
            weights_eul_index_cell(i,2)=inttmp
          end do
        else
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
               fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        end if
        !I have to correct the number
        if (fvm%faceno==5) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2)-nhe-nhe+1
          end do
        end if
        if (fvm%faceno==6) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
          end do
        end if
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
               weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
      end do
      
      if (EOC) then
        jallactual=jall
        ! for Eulerian Correction  
        do jx=nc,nc+2      
          do jy=-1,nc+2
            if ((jx==nc) .and. (jy>0) .and. (jy<nc+1)) then 
              cycle
            endif
            call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
            if(swap1) then !flip orientation
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                   fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
              do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1)
                weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
                weights_eul_index_cell(i,2)=inttmp
              end do
            else
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min1,jy_min1,nreconstruction,&
                   fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            end if
            !I have to correct the number - xphl?????
            if (fvm%faceno==5) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,2)=jy_max1-jy_min1-weights_eul_index_cell(i,2)-nhe-nhe+1
              end do
            end if
            if (fvm%faceno==6) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
              end do
            end if
            if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
            endif
          end do
        end do
      
        do ja=jallactual_eul,jall-1
          jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
          da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
          centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
        end do
      
        jall=jallactual
      endif
    endif
    
    !NORTH SIDE 
    if (fvm%cubeboundary == north) then
      ! calculate xy Cartesian on the cube of departure points on the corresponding face  
      do jy=nc,nc+3
        do jx=-1,nc+3
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(north),dcart(jx,jy))                  
        end do
      end do
      jy=nc
      do jx=1,nc
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
        if(swap1) then !flip orientation
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
               fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          do i=1,jcollect_cell
            inttmp=weights_eul_index_cell(i,1)
            weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
            weights_eul_index_cell(i,2)=inttmp
          end do
        else
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
               fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
        end if
        !I have to correct the number - xphl?
        if (fvm%faceno==2) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
          end do
        end if
        if (fvm%faceno==3) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1))-nhe-nhe+1
            weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
          end do
        end if
        if (fvm%faceno==4) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
          end do
        end if
        if (fvm%faceno==6) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
            weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
          end do
        end if
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
               weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
      end do
      
      if(EOC) then
        jallactual=jall
        ! for Eulerian Correction      
        do jy=nc,nc+2      
          do jx=-1,nc+2
            if ((jy==nc) .and. (jx>0) .and. (jx<nc+1)) then 
              cycle
            endif
            call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
            if(swap1) then !flip orientation
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                   fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
              do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1)
                weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
                weights_eul_index_cell(i,2)=inttmp
              end do
            else
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min1,jy_min1,nreconstruction,&
                   fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
            end if
            !I have to correct the number
            if (fvm%faceno==2) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
              end do
            end if
            if (fvm%faceno==3) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1))-nhe-nhe+1
                weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
              end do
            end if
            if (fvm%faceno==4) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
              end do
            end if
            if (fvm%faceno==6) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
                weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
              end do
            end if
            if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
            endif
          end do
        end do
      
        do ja=jallactual_eul,jall-1
          jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
          da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
          centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
        end do
  !     
        jall=jallactual    
      end if
     endif
    !SOUTH SIDE
    if (fvm%cubeboundary == south) then
      ! calculate xy Cartesian on the cube of departure points on the corresponding face  
      do jy=-1,2
        do jx=-1,nc+3
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(south),dcart(jx,jy))                   
        end do
      end do
      jy=1
      do jx=1,nc
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap1) then !flip orientation
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
               fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          do i=1,jcollect_cell
            inttmp=weights_eul_index_cell(i,1)
            weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
            weights_eul_index_cell(i,2)=inttmp
          end do
        else
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
               fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
        end if
        !I have to correct the number - xphl????
        if  (fvm%faceno==2) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          end do
        end if
        if  (fvm%faceno==3) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
            weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
          end do
        end if
        if  (fvm%faceno==4) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
          end do
        end if
        if  (fvm%faceno==5) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
            weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
          end do
        end if
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
               weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
      end do
 
      if(EOC) then     
        jallactual=jall
        !for Eulerian correction (Erath et al., 2013)  
        do jy=-1,1      
          do jx=-1,nc+2
            if ((jy==1) .and. (jx>0) .and. (jx<nc+1)) then 
              cycle
            endif
            call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
            if(swap1) then !flip orientation
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                   fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
              do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1)
                weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
                weights_eul_index_cell(i,2)=inttmp
              end do
            else
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min1,jy_min1,nreconstruction,&
                   fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)  
            end if
            !I have to correct the number - xphl???
            if  (fvm%faceno==2) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
              end do
            end if
            if  (fvm%faceno==3) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
                weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
              end do
            end if
            if  (fvm%faceno==4) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
              end do
            end if
            if  (fvm%faceno==5) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
                weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
              end do
            end if
            if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
            endif
          end do
        end do
      
        do ja=jallactual_eul,jall-1
          jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
          da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
          centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
        end do
      
        jall=jallactual    
      end if
    endif
    !SOUTHWEST Corner
    if (fvm%cubeboundary == swest) then
      ! calculate xy Cartesian on the cube of departure points on the corresponding face  
      do jy=-1,2
        do jx=-1,nc+3
          if ((jy<1) .and. (jx<1)) then   ! in the southwest corner are no values!!!
            cycle
          end if
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(south),dcart(jx,jy))                  
        end do
      end do
      jy=1
      do jx=1,nc
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap1) then !flip orientation
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
               fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          do i=1,jcollect_cell
            inttmp=weights_eul_index_cell(i,1)
            weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
            weights_eul_index_cell(i,2)=inttmp
          end do
        else
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
               fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        end if
        !I have to correct the number - xphl????
        if ((fvm%faceno==3) .OR. (fvm%faceno==5)) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
            weights_eul_index_cell(i,2)=(weights_eul_index_cell(i,2)+nhe-1)
          end do
        end if
        if (fvm%faceno==2) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
          end do
        end if
        if (fvm%faceno==4) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
          end do
        end if
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
               weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
      end do
      
      if(EOC) then
        jallactual=jall
        !for Eulerian correction   
        do jy=-1,1      
          do jx=-1,nc+2
            if (((jy<1) .and. (jx<1)) .or. ((jy==1) .and. (jx>0) .and. (jx<nc+1))) then 
              cycle
            endif
            call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
            if(swap1) then !flip orientation
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                   fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
              do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1)
                weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
                weights_eul_index_cell(i,2)=inttmp
              end do
            else
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min1, jy_min1,nreconstruction,&
                   fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            end if
            !I have to correct the number - xphl???
            if ((fvm%faceno==3) .OR. (fvm%faceno==5)) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
                weights_eul_index_cell(i,2)=(weights_eul_index_cell(i,2)+nhe-1)
              end do
            end if
            if (fvm%faceno==2) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
              end do
            end if
            if (fvm%faceno==4) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
              end do
            end if
            if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
            endif
          end do
        end do
      
        do ja=jallactual_eul,jall-1
          jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
          da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
          centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
        end do
      
      
        jall=jallactual    
        jallactual_eul=jall  
      endif
      
      ! calculate xy Cartesian on the cube of departure points on the corresponding face  
      do jx=-1,2
        do jy=-1,nc+3
          if ((jy<1) .and. (jx<1)) then
            cycle
          end if
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
               fvm%nbrsface(west),dcart(jx,jy))                  
        end do
      end do
      jx=1
      do jy=1,nc
        call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
        if(swap2) then !flip orientation
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
               fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          do i=1,jcollect_cell
            inttmp=weights_eul_index_cell(i,1)
            weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
            weights_eul_index_cell(i,2)=inttmp
          end do
        else
          call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
               fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
               tmp,ngpc,gsweights,gspts,&
               weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        end if
        !I have to correct the number - xphl????
        if (fvm%faceno==5) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
          end do
        end if
        if (fvm%faceno==6) then
          do i=1,jcollect_cell
            weights_eul_index_cell(i,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2)+1
          end do
        end if
        if (jcollect_cell>0) then
          weights_all(jall:jall+jcollect_cell-1,:) = weights_cell (1:jcollect_cell,:)
          weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
               weights_eul_index_cell(1:jcollect_cell,:)
          weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
          weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
          jall = jall+jcollect_cell          
        endif
      end do
      
      if(EOC) then
        jallactual=jall
        !for Eulerian correction   
      
        do jx=-1,1      
          do jy=-1,nc+2
            if (((jx<1) .and. (jy<1)) .or. ((jx==1) .and. (jy>0) .and. (jy<nc+1))) then 
              cycle
            endif
            call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
            if(swap2) then !flip orientation
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                   fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
              do i=1,jcollect_cell
                inttmp=weights_eul_index_cell(i,1)
                weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
                weights_eul_index_cell(i,2)=inttmp
              end do
            else
              call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min2, jy_min2,nreconstruction,&
                   fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                   tmp,ngpc,gsweights,gspts,&
                   weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            end if
            !I have to correct the number - xphl???
            if (fvm%faceno==5) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
              end do
            end if
            if (fvm%faceno==6) then
              do i=1,jcollect_cell
                weights_eul_index_cell(i,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2)+1
              end do
            end if
            if (jcollect_cell>0) then
              weights_all(jall:jall+jcollect_cell-1,:) = weights_cell (1:jcollect_cell,:)
              weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                   weights_eul_index_cell(1:jcollect_cell,:)
              weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
              weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
              jall = jall+jcollect_cell          
            endif
          end do
        end do
      
        do ja=jallactual_eul,jall-1
          jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
          da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
          centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
        end do
        !     
        jall=jallactual      
      
      endif
    endif
    
    ! SOUTHEAST Corner
    if (fvm%cubeboundary == seast) then
      ! calculate xy Cartesian on the cube of departure points on the corresponding face  
      do jy=-1,2
        do jx=-1,nc+3
          if ((jy<1) .and. (jx>nc+1)) then   ! in the southwest corner are no values!!!
            cycle
          end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(south),dcart(jx,jy))  
      end do
    end do
    jy=1
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
             fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number   
      if  (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
        end do
      end if
      if  (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if  (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
             weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
    
    if(EOC) then
      jallactual=jall
      ! Eulerian weight correction    
      do jy=-1,1
        do jx=-1,nc+2
          if (((jy<1) .and. (jx>nc)) .or. ((jy==1) .and. (jx>0) .and. (jx<nc+1))) then 
            cycle
          endif
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap1) then !flip orientation
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                 fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            do i=1,jcollect_cell
              inttmp=weights_eul_index_cell(i,1)
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
              weights_eul_index_cell(i,2)=inttmp
            end do
          else
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min1, jy_min1,nreconstruction,&
                 fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number - xphl?
          if  (fvm%faceno==2) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
            end do
          end if
          if  (fvm%faceno==3) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
              weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
            end do
          end if
          if  (fvm%faceno==4) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
            end do
          end if
          if  (fvm%faceno==5) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1)+nhe)-nhe+1
              weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)+nhe-1
            end do
          end if
          if (jcollect_cell>0) then
            weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
            weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                 weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
            weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
            jall = jall+jcollect_cell          
          endif
        end do
      end do
    
      do ja=jallactual_eul,jall-1
        jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
        da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
        centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
      end do
    
      jall=jallactual    
      jallactual_eul=jall
    endif
    
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=nc,nc+3
      do jy=-1,nc+3
        if ((jy<1) .and. (jx>nc+1)) then
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(east),dcart(jx,jy))                   
      end do
    end do
    jx=nc
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
      if(swap2) then !flip orientation
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
             fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number - xphl?
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2)+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
             weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
    
    if(EOC) then
      jallactual=jall
      ! Eulerian weight correction        
      do jx=nc,nc+2
        do jy=-1,nc+2
          if (((jx>nc) .and. (jy<1)) .or. ((jx==nc) .and. (jy>0) .and. (jy<nc+1))) then 
            cycle
          endif
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)     
          if(swap2) then !flip orientation
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                 fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            do i=1,jcollect_cell
              inttmp=weights_eul_index_cell(i,1)
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
              weights_eul_index_cell(i,2)=inttmp
            end do
          else
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min2, jy_min2,nreconstruction,&
                 fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number - xphl ???
          if (fvm%faceno==5) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,2)=(jy_max2-jy_min2)-weights_eul_index_cell(i,2)+1
            end do
          end if
          if (fvm%faceno==6) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
            end do
          end if
          if (jcollect_cell>0) then
            weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
            weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                 weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
            weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
            jall = jall+jcollect_cell          
          endif
        end do
      end do
    
      do ja=jallactual_eul,jall-1
        jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
        da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
        centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
      end do
      !     
      jall=jallactual
    
    endif
  endif
  
  !NORTHEAST Corner
  if (fvm%cubeboundary == neast) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=nc,nc+3
      do jx=-1,nc+3
        if ((jy>nc+1) .and. (jx>nc+1)) then   ! in the southwest corner are no values!!!
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(north),dcart(jx,jy))                
      end do
    end do
    jy=nc
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
             fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number - xphl????
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=(weights_eul_index_cell(i,2))-nhe+1
        end do
      end if
      if (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
             weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
    
    if(EOC) then
      jallactual=jall
      !Eulerian correction
      do jy=nc,nc+2      
        do jx=-1,nc+2
          if (((jy>nc) .and. (jx>nc)) .or. ((jy==nc) .and. (jx>0) .and. (jx<nc+1))) then 
            cycle
          endif
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap1) then !flip orientation
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                 fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            do i=1,jcollect_cell
              inttmp=weights_eul_index_cell(i,1)
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
              weights_eul_index_cell(i,2)=inttmp
            end do
          else
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min1, jy_min1,nreconstruction,&
                 fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number - xphl????
          if (fvm%faceno==2) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,2)=(weights_eul_index_cell(i,2))-nhe+1
            end do
          end if
          if (fvm%faceno==3) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
              weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
            end do
          end if
          if (fvm%faceno==6) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
              weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
            end do
          end if
          if (fvm%faceno==4) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)-nhe-nhe+1
            end do
          end if
          if (jcollect_cell>0) then
            weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
            weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                 weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
            weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
            jall = jall+jcollect_cell          
          endif
        end do
      end do
    
      do ja=jallactual_eul,jall-1
        jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
        da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
        centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
      end do
    
    
      jall=jallactual    
      jallactual_eul=jall  
    endif
    
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=nc,nc+3
      do jy=-1,nc+3
        if ((jy>nc+1) .and. (jx>nc+1)) then
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(east),dcart(jx,jy))                  
      end do
    end do
    jx=nc
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
             fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number - xphl????
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2)-nhe-nhe+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
        end do
      end if
      
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
             weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
    
    if(EOC) then
      jallactual=jall
      !Eulerian correction
      do jx=nc, nc+2
        do jy=-1,nc+2
          if (((jx>nc) .and. (jy>nc)) .or. ((jx==nc) .and. (jy>0) .and. (jy<nc+1))) then 
            cycle
          endif
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap2) then !flip orientation
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                 fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            do i=1,jcollect_cell
              inttmp=weights_eul_index_cell(i,1)
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
              weights_eul_index_cell(i,2)=inttmp
            end do
          else
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min2, jy_min2,nreconstruction,&
                 fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number - xphl????
          if (fvm%faceno==5) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2)-nhe-nhe+1
            end do
          end if
          if (fvm%faceno==6) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)-nhe+1
            end do
          end if
        
          if (jcollect_cell>0) then
            weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
            weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                 weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
            weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
            jall = jall+jcollect_cell          
          endif
        end do
      end do
    
      do ja=jallactual_eul,jall-1
        jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
        da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
        centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
      end do
      !     
      jall=jallactual    
    
    endif
  endif
  
  !NORTH WEST CORNER 
  if (fvm%cubeboundary == nwest) then
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jy=nc, nc+3
      do jx=-1,nc+3
        if ((jy>nc+1) .and. (jx<1)) then   ! in the southwest corner are no values!!!
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(north),dcart(jx,jy))                   
      end do
    end do
    jy=nc
    do jx=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap1) then !flip orientation
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
             fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number - xphl???
      if (fvm%faceno==2) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==3) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1))+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (fvm%faceno==4) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
          weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
             weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
    
    if(EOC) then
      jallactual=jall
      !Eulerian correction
      do jy=nc,nc+2
        do jx=-1,nc+2
          if (((jy>nc) .and. (jx<1)) .or. ((jy==nc) .and. (jx>0) .and. (jx<nc+1))) then 
            cycle
          endif
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap1) then !flip orientation
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min1,jx_min1,nreconstruction,&
                 fvm%acarty1,fvm%acartx1,jy_min1, jy_max1, jx_min1, jx_max1, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            do i=1,jcollect_cell
              inttmp=weights_eul_index_cell(i,1)
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
              weights_eul_index_cell(i,2)=inttmp
            end do
          else
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min1, jy_min1,nreconstruction,&
                 fvm%acartx1,fvm%acarty1,jx_min1, jx_max1, jy_min1, jy_max1, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number - xphl????
          if (fvm%faceno==2) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
            end do
          end if
          if (fvm%faceno==3) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-(weights_eul_index_cell(i,1))+1
              weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
            end do
          end if
          if (fvm%faceno==4) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
            end do
          end if
          if (fvm%faceno==6) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=(jx_max1-jx_min1)-weights_eul_index_cell(i,1)+1
              weights_eul_index_cell(i,2)=weights_eul_index_cell(i,2)-nhe+1
            end do
          end if
          if (jcollect_cell>0) then
            weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
            weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                 weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
            weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
            jall = jall+jcollect_cell          
          endif
        end do
      end do
      do ja=jallactual_eul,jall-1
        jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
        da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
        centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
      end do
      jall=jallactual    
      jallactual_eul=jall  
    endif
    
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    do jx=-1,2
      do jy=-1,nc+3
        if ((jy>nc+1) .and. (jx<1)) then
          cycle
        end if
        call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(jx,jy,klev)),&
             fvm%nbrsface(west),dcart(jx,jy))                  
      end do
    end do
    jx=1
    do jy=1,nc
      call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
      if(swap2) then !flip orientation
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
             fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
        do i=1,jcollect_cell
          inttmp=weights_eul_index_cell(i,1)
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
          weights_eul_index_cell(i,2)=inttmp
        end do
      else
        call compute_weights_cell(4,.not.EOC,xcell,ycell,jx,jy,nreconstruction,&
             fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
             tmp,ngpc,gsweights,gspts,&
             weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
      end if
      !I have to correct the number - xphl????
      if (fvm%faceno==5) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
        end do
      end if
      if (fvm%faceno==6) then
        do i=1,jcollect_cell
          weights_eul_index_cell(i,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2)-nhe-nhe+1
        end do
      end if
      if (jcollect_cell>0) then
        weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
        weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
             weights_eul_index_cell(1:jcollect_cell,:)
        weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
        weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
        jall = jall+jcollect_cell          
      endif
    end do
    
    if(EOC) then
      jallactual=jall
      !Eulerian correction
      do jx=-1,1
        do jy=-1,nc+2
          if (((jx<1) .and. (jy>nc)) .or. ((jx==1) .and. (jy>0) .and. (jy<nc+1))) then 
            cycle
          endif
          call getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)      
          if(swap2) then !flip orientation
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jy_min2,jx_min2,nreconstruction,&
                 fvm%acarty2,fvm%acartx2,jy_min2, jy_max2, jx_min2, jx_max2, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
            do i=1,jcollect_cell
              inttmp=weights_eul_index_cell(i,1)
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,2)
              weights_eul_index_cell(i,2)=inttmp
            end do
          else
            call compute_weights_cell(4,.not.EOC,xcell,ycell,jx_min2, jy_min2,nreconstruction,&
                 fvm%acartx2,fvm%acarty2,jx_min2, jx_max2, jy_min2, jy_max2, &
                 tmp,ngpc,gsweights,gspts,&
                 weights_cell,weights_eul_index_cell,jcollect_cell,jmax_segments_cell)
          end if
          !I have to correct the number - xphl????
          if (fvm%faceno==5) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,1)=weights_eul_index_cell(i,1)+nhe-1
            end do
          end if
          if (fvm%faceno==6) then
            do i=1,jcollect_cell
              weights_eul_index_cell(i,2)=jy_max2-jy_min2-weights_eul_index_cell(i,2)-nhe-nhe+1
            end do
          end if
          if (jcollect_cell>0) then
            weights_all(jall:jall+jcollect_cell-1,:) = weights_cell(1:jcollect_cell,:)
            weights_eul_index_all(jall:jall+jcollect_cell-1,:) = &
                 weights_eul_index_cell(1:jcollect_cell,:)
            weights_lgr_index_all(jall:jall+jcollect_cell-1,1) = jx
            weights_lgr_index_all(jall:jall+jcollect_cell-1,2) = jy
            jall = jall+jcollect_cell          
          endif
        end do
      end do
    
      do ja=jallactual_eul,jall-1
        jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
        da_cslam(jx,jy) = da_cslam(jx,jy)+weights_all(ja,1)
        centroid_cslam(jx,jy,:) = centroid_cslam(jx,jy,:)+weights_all(ja,2:6)
      end do
      !     
      jall=jallactual       
    endif
  endif
  
  jall=jall-1

  if(EOC) then
     !
     ! here is the correction (Erath et al., 2013), uncomment it if you want to run the scheme 
     ! without weight correction
     !  
     do ja=1,jall
        jx = weights_eul_index_all(ja,1); jy = weights_eul_index_all(ja,2);
        area=fvm%area_sphere(jx,jy)
        weights_all(ja,1) = weights_all(ja,1)*abs(area/da_cslam(jx,jy))
        do k=2,6
           if (ABS(centroid_cslam(jx,jy,k-1))>tiny) &
                weights_all(ja,k) = weights_all(ja,k)*&
                abs(area*fvm%spherecentroid(jx,jy,k-1)/centroid_cslam(jx,jy,k-1))
        end do
     end do
  endif
 
end subroutine compute_weights  


!END SUBROUTINE COMPUTE_WEIGHTS-------------------------------------------CE-for FVM!

! ----------------------------------------------------------------------------------!
!SUBROUTINE GETDEP_CELLBOUNDARIESXYVEC------------------------------------CE-for FVM!
! AUTHOR: CHRISTOPH ERATH, 20.March 2012                                            !
! DESCRIPTION: needed to apply the search subroutines                               !
!                                                                                   !
! CALLS: cart2cubedspherexy                                                         !
! INPUT:  jx,jy   ... index of the cell                                             !
!         dcart ... Cartesian coordinates of the cell on the corresponding face     !
! OUTPUT: xcell, ycell ... x and y Cartesian coordinates on the cube of the cell    !
!                          dsphere                                                  !
!-----------------------------------------------------------------------------------!
subroutine getdep_cellboundariesxyvec(xcell,ycell,jx,jy,dcart)
use coordinate_systems_mod, only : cartesian2D_t
  implicit none
  !
  ! dimension(0:5)
  !
  real (kind=real_kind), dimension(0:), intent(inout)     :: xcell,ycell
  integer (kind=int_kind), intent(in)                      :: jx, jy
!  type (cartesian2D_t), intent(in)                         :: dcart(-1:nc+3,-1:nc+3)
  type (cartesian2D_t), intent(in)                         :: dcart(-1:,-1:)

  xcell(1) = dcart(jx  ,jy  )%x 
  ycell(1) = dcart(jx  ,jy  )%y
  xcell(2) = dcart(jx  ,jy+1)%x 
  ycell(2) = dcart(jx  ,jy+1)%y
  xcell(3) = dcart(jx+1,jy+1)%x 
  ycell(3) = dcart(jx+1,jy+1)%y
  xcell(4) = dcart(jx+1,jy  )%x 
  ycell(4) = dcart(jx+1,jy  )%y

  xcell(5) = xcell(1)        
  ycell(5) = ycell(1)          
  xcell(0) = xcell(4)        
  ycell(0) = ycell(4)
end subroutine getdep_cellboundariesxyvec

!END SUBROUTINE GETDEP_CELLBOUNDARIESXYVEC--------------------------------CE-for FVM!
    
! ----------------------------------------------------------------------------------!
!NEXT SUBROUTINES ARE TAKEN FROM------------------------------------------CE-for FVM!
! PETER LAURITZENs CODE adapted to use in HOMME                                    !
! ----------------------------------------------------------------------------------!

  subroutine compute_weights_cell(nvertex,lexact_horizontal_line_integrals,&
       xcell_in,ycell_in,jx,jy,nreconstruction,xgno,ygno,&
       jx_min, jx_max, jy_min, jy_max,tmp,&
       ngauss,gauss_weights,abscissae,weights,weights_eul_index,jcollect,jmax_segments)

    implicit none
    integer (kind=int_kind), intent(in) :: nvertex
    logical, intent(in) :: lexact_horizontal_line_integrals
    integer (kind=int_kind)                  , intent(in):: nreconstruction, jx,jy,ngauss,jmax_segments
    !
    ! dimension(nvertex)
    !
    real (kind=real_kind)   ,  dimension(:), intent(in):: xcell_in,ycell_in
    !
    integer (kind=int_kind), intent(in)               :: jx_min, jy_min, jx_max, jy_max
    !
    ! dimension(-nhe:nc+2+nhe)
    !
    real (kind=real_kind), dimension(-nhe:), intent(in) :: xgno, ygno
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
    integer (kind=int_kind), intent(out)       :: jcollect
    !
    ! local workspace
    !
    !
    ! max number of line segments is:
    !
    ! (number of longitudes)*(max average number of crossings per line segment = 3)*ncube*2
    !
    real (kind=real_kind)   ,  &
         dimension(jmax_segments,nreconstruction), intent(out) :: weights
    integer (kind=int_kind),  &
         dimension(jmax_segments,2), intent(out)      :: weights_eul_index
    
    integer (kind=int_kind) :: jsegment,i
    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer (kind=int_kind)  :: jcross_lat, iter
    !
    ! max. crossings per side is 2*nhe
    !
    real (kind=real_kind), &
         dimension(8*nhe,2) :: r_cross_lat
    integer (kind=int_kind), &
         dimension(8*nhe,2) :: cross_lat_eul_index
    real (kind=real_kind)   ,  dimension(nvertex) :: xcell,ycell

    xcell = xcell_in(1:nvertex)
    ycell = ycell_in(1:nvertex)

    jsegment   = 0
    weights    = 0.0D0
    jcross_lat = 0
          
!    if (jx==4.and.jy==3) then
!       ldbg=.true.
!    else
!       ldbg=.true.
!    end if
    if (ldbg) write(*,*) "going into side_integral"
!    if (ldbg) write(*,*) "cell is:"
!    do i=1,nvertex
!       if (ldbg) write(*,*) xcell(i),ycell(i)
!    end do
!    if (ldbg) write(*,*) "are line-integrals being computed exactly?",lexact_horizontal_line_integrals
    call side_integral(lexact_horizontal_line_integrals,xcell,ycell,nvertex,jsegment,jmax_segments,&
         weights,weights_eul_index,nreconstruction,jx,jy,xgno,ygno,jx_min, jx_max, jy_min, jy_max,&
         ngauss,gauss_weights,abscissae,&
         jcross_lat,r_cross_lat,cross_lat_eul_index)
    !
    !**********************
    ! 
    ! Do inner integrals
    !
    !**********************
    !    
    if (ldbg) write(*,*) "going into compute_inner_line_integrals_lat"
    call compute_inner_line_integrals_lat(lexact_horizontal_line_integrals,&
         r_cross_lat,cross_lat_eul_index,&
         jcross_lat,jsegment,jmax_segments,xgno,jx_min, jx_max, jy_min, jy_max,&
         weights,weights_eul_index,&
         nreconstruction,ngauss,gauss_weights,abscissae)
    
    
    IF (ABS((jcross_lat/2)-DBLE(jcross_lat)/2.0)>tiny) then
      WRITE(*,*) "number of latitude crossings are not even: ABORT",jcross_lat,jx,jy
      STOP
    END IF
    
    !
    ! collect line-segment that reside in the same Eulerian cell
    !
    if (jsegment>0) then
      call collect(weights,weights_eul_index,nreconstruction,jcollect,jsegment,jmax_segments)
      !
      ! DBG
      !
!      tmp=0.0
!      !       WRITE(*,*) "max area",maxval(weights(1:jcollect,1))
!      !       WRITE(*,*) "min area",minval(weights(1:jcollect,1))
!      !       stop
!      do i=1,jcollect     
!        IF (weights(i,1)<-tiny) THEN
!           write(*,*) "weights from compute_weights_cell negative:",weights(i,1),i
!           write(*,*) "index jx,jy ",jx,jy
!        END IF
!        !       end do
!        
!        tmp=tmp+weights(i,1)
!      enddo
!
!      IF (abs(tmp)>0.04) THEN
!        WRITE(*,*) "sum of weights seems too large",tmp
!!dbg        stop
!      END IF
!      IF (tmp<-1.0E-9) THEN
!        WRITE(*,*) "sum of weights is negative - negative area?",tmp,jx,jy
! !       stop
!      END IF
    else
      jcollect = 0
    end if
  end subroutine compute_weights_cell


  !
  !****************************************************************************
  !
  ! organize data and store it
  !
  !****************************************************************************
  !
  subroutine collect(weights,weights_eul_index,nreconstruction,jcollect,jsegment,jmax_segments)
    implicit none
    integer (kind=int_kind),                                  INTENT(IN   ) :: jsegment,jmax_segments
    integer (kind=int_kind)                                 , intent(in)    :: nreconstruction
    !
    real (kind=real_kind)  , dimension(:,:), intent(inout) :: weights !dimension(jmax_segments,nreconstruction)
    integer (kind=int_kind), dimension(:,:), intent(inout) :: weights_eul_index !dimension(jmax_segments,2)
    integer (kind=int_kind),                                  INTENT(OUT  ) :: jcollect
    !
    ! local workspace
    !
    integer (kind=int_kind) :: imin, imax, jmin, jmax, i,j,k,h
    logical                 :: ltmp

    real (kind=real_kind)  , dimension(jmax_segments,nreconstruction) :: weights_out
    integer (kind=int_kind), dimension(jmax_segments,2     ) :: weights_eul_index_out

    weights_out           = 0.0D0
    weights_eul_index_out = -100

    imin = MINVAL(weights_eul_index(1:jsegment,1))
    imax = MAXVAL(weights_eul_index(1:jsegment,1))
    jmin = MINVAL(weights_eul_index(1:jsegment,2))
    jmax = MAXVAL(weights_eul_index(1:jsegment,2))

    ltmp = .FALSE.

    jcollect = 1

    do j=jmin,jmax
       do i=imin,imax
          do k=1,jsegment
             if (weights_eul_index(k,1)==i.AND.weights_eul_index(k,2)==j) then
                weights_out(jcollect,1:nreconstruction) = &
                weights_out(jcollect,1:nreconstruction) + weights(k,1:nreconstruction)
                if (ldbg) write(*,*) "qqqq eul index ",i,j,weights(k,1)
                ltmp = .TRUE.
                h = k
             endif
          enddo
          if (ltmp) then
             weights_eul_index_out(jcollect,:) = weights_eul_index(h,:)
             jcollect = jcollect+1
          endif
          ltmp = .FALSE.
       enddo
    enddo
    jcollect = jcollect-1
    weights           = weights_out
    weights_eul_index = weights_eul_index_out
  end subroutine collect

  !    
  ! compute crossings with Eulerian latitudes and longitudes
  !
  !

  !
  !*****************************************************************************************
  !
  ! 
  !
  !*****************************************************************************************
  !
  subroutine compute_inner_line_integrals_lat(lexact_horizontal_line_integrals,r_cross_lat,&
       cross_lat_eul_index,&
       jcross_lat,jsegment,jmax_segments,xgno,jx_min,jx_max,jy_min, jy_max,weights,weights_eul_index,&
       nreconstruction,ngauss,gauss_weights,abscissae)    
    implicit none
    logical, intent(in) :: lexact_horizontal_line_integrals
    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer (kind=int_kind),         intent(in):: jcross_lat, jmax_segments,nreconstruction,ngauss
    integer (kind=int_kind),         intent(inout):: jsegment
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae
    !
    ! max. crossings per side is 2*nhe
    !
    
    real (kind=real_kind)  , dimension(:,:), intent(in):: r_cross_lat ! dimension(8*nhe,2)
    integer (kind=int_kind), dimension(:,:), intent(in):: cross_lat_eul_index ! ! dimension(8*nhe,2)
    integer (kind=int_kind)                , intent(in):: jx_min, jx_max, jy_min, jy_max

    real (kind=real_kind), dimension(-nhe:), intent(in)  :: xgno !dimension(-nhe:nc+2+nhe)
    !
    ! dimension(jmax_segments,nreconstruction)
    !
    real (kind=real_kind), dimension(:,:), intent(inout) :: weights
    !
    ! dimension(jmax_segments,2)
    !
    integer (kind=int_kind), dimension(:,:), intent(inout)               :: weights_eul_index

    real (kind=real_kind)   , dimension(nreconstruction)         :: weights_tmp
    integer (kind=int_kind) :: imin, imax, jmin, jmax, i,j,k, isgn, h, eul_jx, eul_jy
    integer (kind=int_kind) :: idx_start_y,idx_end_y
    logical                 :: ltmp,lcontinue
    real (kind=real_kind), dimension(2)  :: rstart,rend,rend_tmp
    real (kind=real_kind), dimension(2)  :: xseg, yseg
    5   FORMAT(10e14.6)
    lcontinue = .TRUE.
    if (jcross_lat>0) then
       do i=MINVAL(cross_lat_eul_index(1:jcross_lat,2)),MAXVAL(cross_lat_eul_index(1:jcross_lat,2))
          !
          ! find "first" crossing with Eulerian cell i
          !
          do k=1,jcross_lat
             if (cross_lat_eul_index(k,2)==i) exit
          enddo
          do j=k+1,jcross_lat
             !
             ! find "second" crossing with Eulerian cell i
             !
             if (cross_lat_eul_index(j,2)==i) then
                if (r_cross_lat(k,1)<r_cross_lat(j,1)) then
                   rstart = r_cross_lat(k,1:2)
                   rend   = r_cross_lat(j,1:2)
                   imin   = cross_lat_eul_index(k,1)
                   imax   = cross_lat_eul_index(j,1)
                else
                   rstart = r_cross_lat(j,1:2)
                   rend   = r_cross_lat(k,1:2)
                   imin   = cross_lat_eul_index(j,1)
                   imax   = cross_lat_eul_index(k,1)
                endif
!                write(*,*) "from inner",rstart,rend
!                write(*,*) "from h,i  ",h,i
                do h=imin,imax
                   if (h==imax) then
                      rend_tmp = rend
                   else
                      rend_tmp(1) = xgno(h+1)
                      rend_tmp(2) = r_cross_lat(k,2)
                   endif
                   xseg(1) = rstart(1)
                   xseg(2) = rend_tmp(1)
                   yseg(1) = rstart(2)
                   yseg(2) = rend_tmp(2)
                   call get_weights_exact(lexact_horizontal_line_integrals, weights_tmp,xseg,yseg,&
                        nreconstruction,ngauss,gauss_weights,abscissae)

!phl                   if (i.LE.nc.AND.i.GE.1.AND.h.LE.nc.AND.h.GE.1) then
                   if (i.LE.jy_max-1.AND.i.GE.jy_min.AND.h.LE.jx_max-1.AND.h.GE.jx_min) then
                      jsegment=jsegment+1
                      weights_eul_index(jsegment,1) = h 
                      weights_eul_index(jsegment,2) = i
                      weights(jsegment,1:nreconstruction) = -weights_tmp
                      if (ldbg) write(*,*) "from/to ",xseg(1),yseg(1),xseg(2),yseg(2)
                      if (ldbg) write(*,*) "from inner",weights_tmp
                   endif
                   !
                   ! subtract the same weights on the west side of the line
                   !
                   if (i.LE.jy_max.AND.i.GE.jy_min+1.AND.h.LE.jx_max-1.AND.h.GE.jx_min) then
!phl                   if (i.GE.2.AND.i.LE.nc+1.AND.h.LE.nc.AND.h.GE.1) then
                      jsegment = jsegment+1
                      weights_eul_index(jsegment,1) = h 
                      weights_eul_index(jsegment,2) = i-1
                      weights(jsegment,1:nreconstruction) = weights_tmp
                   endif
                   !
                   ! prepare for next iteration
                   !
                   rstart = rend_tmp
                enddo
             endif
          enddo
       enddo
    endif
  end subroutine compute_inner_line_integrals_lat

  !
  ! line integral from (a1_in,a2_in) to (b1_in,b2_in)
  ! If line is coniciding with an Eulerian longitude or latitude the routine
  ! needs to know where an adjacent side is located to determine which
  ! reconstruction must be used. therefore (c1,c2) is passed to the routine
  !
  !   

  subroutine side_integral(lexact_horizontal_line_integrals,&
       x_in,y_in,nvertex,jsegment,jmax_segments,&
       weights,weights_eul_index,nreconstruction,jx,jy,xgno,ygno,jx_min,jx_max,jy_min,jy_max,&
       ngauss,gauss_weights,abscissae,&!)!phl add jx_min etc.
       jcross_lat,r_cross_lat,cross_lat_eul_index)
    implicit none


    logical, intent(in) :: lexact_horizontal_line_integrals
    integer (kind=int_kind),            intent(in)    :: nreconstruction,jx,jy,jmax_segments,ngauss
    integer (kind=int_kind), intent(in)               :: nvertex
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(:), intent(in) :: gauss_weights, abscissae !dimension(ngauss)
    real (kind=real_kind), dimension(:)        , intent(in)    :: x_in,y_in !dimension(1:nvertex)

    integer (kind=int_kind), intent(in)               :: jx_min, jy_min, jx_max, jy_max
    real (kind=real_kind), dimension(-nhe:), intent(in) :: xgno, ygno !dimension(-nhe:nc+2+nhe)
    integer (kind=int_kind),            intent(inout) :: jsegment
!    integer (kind=int_kind),dimension(0:2),intent(in)    :: jx_eul_in, jy_eul_in
    real (kind=real_kind)   , dimension(:,:), intent(out) :: weights !dimension(jmax_segments,nreconstruction)
    integer (kind=int_kind),  &
         dimension(jmax_segments,2), intent(out) :: weights_eul_index

    !
    ! variables for registering crossings with Eulerian latitudes and longitudes
    !
    integer (kind=int_kind),         intent(inout):: jcross_lat
    !
    ! max. crossings per side is 2*nhe
    !
    real (kind=real_kind), &
         dimension(8*nhe,2), intent(inout):: r_cross_lat
    integer (kind=int_kind), &
         dimension(8*nhe,2), intent(inout):: cross_lat_eul_index
    !
    ! local variables
    !
    real (kind=real_kind) :: dist_lon,dist_lat, tmp_a1, tmp_a2, tmp_x(1), tmp_b2, a1, a2, b2
    real (kind=real_kind) :: dist
    real (kind=real_kind), dimension(2) :: xseg,yseg 
    real (kind=real_kind), dimension(0:3) :: x,y
    real (kind=real_kind)               :: lon_x,lat_y,lon_y,lat_x
    real (kind=real_kind)               :: xeul,yeul,xcross,ycross,slope
    integer (kind=int_kind) ::    jx_eul_tmp,jy_eul_tmp
    integer (kind=int_kind)            :: xsgn1,ysgn1,xsgn2,ysgn2
    integer (kind=int_kind) :: ifrom_left, iter,previous_jy_eul_cross
    logical :: lcontinue, lregister_cross, lsame_cell_x, lsame_cell_y

    integer (kind=int_kind) :: jx_eul, jy_eul, side_count,jdbg
    real (kind=real_kind), dimension(0:nvertex+2)  :: xcell,ycell
    real (kind=real_kind), dimension(0:nvertex+2)  :: xcell_tmp,ycell_tmp
    

5   FORMAT(10e14.6)
    !
    !***********************************************
    !
    ! find jx_eul and jy_eul for (x(1),y(1))
    !
    !***********************************************
    !
    jx_eul = jx; jy_eul = jy    
    xcell(1:nvertex)=x_in; ycell(1:nvertex)=y_in
    DO iter=1,nvertex
      CALL truncate_vertex(xcell(iter),jx_eul,xgno)
      CALL truncate_vertex(ycell(iter),jy_eul,ygno)
    END DO
    xcell(0) = xcell(nvertex); xcell(nvertex+1)=xcell(1); xcell(nvertex+2)=xcell(2);
    ycell(0) = ycell(nvertex); ycell(nvertex+1)=ycell(1); ycell(nvertex+2)=ycell(2);


    IF ((&
         MAXVAL(xcell).LE.xgno(jx_min).OR.MINVAL(xcell).GE.xgno(jx_max).OR.&
         MAXVAL(ycell).LE.ygno(jy_min).OR.MINVAL(ycell).GE.ygno(jy_max))&
         .OR.area(xcell(1:nvertex),ycell(1:nvertex),nvertex).EQ.0.0) THEN
         !
         ! the area call is technically only needed for flux-form CSLAM
         !
!       write(*,*) "area ",area(xcell(1:nvertex),ycell(1:nvertex),nvertex)
!       write(*,*) "xcell(1:nvertex)",xcell(1:nvertex)
!       write(*,*) "ycell(1:nvertex)",ycell(1:nvertex)
      !
      ! entire cell off panel
      !
    ELSE             
      jx_eul = MIN(MAX(jx,jx_min),jx_max)
      jy_eul = MIN(MAX(jy,jy_min),jy_max)
      CALL which_eul_cell(xcell(1:3),jx_eul,xgno)
      CALL which_eul_cell(ycell(1:3),jy_eul,ygno)
      
      side_count = 1
      DO WHILE (side_count<nvertex+1)
        jdbg = 0
        iter = 0
        lcontinue = .TRUE.
        x(0:3) = xcell(side_count-1:side_count+2); y(0:3) = ycell(side_count-1:side_count+2); 
        DO while (lcontinue)
          iter = iter+1
          IF (iter>10) THEN
            WRITE(*,*) "search not converging",iter
            STOP
          END IF
          lsame_cell_x = (x(2).GE.xgno(jx_eul).AND.x(2).LE.xgno(jx_eul+1))
          lsame_cell_y = (y(2).GE.ygno(jy_eul).AND.y(2).LE.ygno(jy_eul+1))
          IF (lsame_cell_x.AND.lsame_cell_y) THEN
            !
            !****************************
            !
            ! same cell integral
            !
            !****************************
            !
            xseg(1) = x(1); yseg(1) = y(1); xseg(2) = x(2); yseg(2) = y(2)
            jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
            lcontinue = .FALSE.
            !
            ! prepare for next side if (x(2),y(2)) is on a grid line
            !
            IF (x(2).EQ.xgno(jx_eul+1).AND.x(3)>xgno(jx_eul+1)) THEN
              !
              ! cross longitude jx_eul+1
              !
              jx_eul=jx_eul+1
            ELSE IF (x(2).EQ.xgno(jx_eul  ).AND.x(3)<xgno(jx_eul)) THEN
              !
              ! cross longitude jx_eul
              !
              jx_eul=jx_eul-1
            END IF
            IF (y(2).EQ.ygno(jy_eul+1).AND.y(3)>ygno(jy_eul+1)) THEN
              !
              ! register crossing with latitude: line-segments point Northward
              !
              jcross_lat = jcross_lat + 1
              jy_eul     = jy_eul     + 1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              cross_lat_eul_index(jcross_lat,2) = jy_eul
              r_cross_lat(jcross_lat,1) = x(2)
              r_cross_lat(jcross_lat,2) = y(2)
!              write(*,*) "A register crossing with latitude",x(2),y(2),jx_eul,jy_eul
            ELSE IF (y(2).EQ.ygno(jy_eul  ).AND.y(3)<ygno(jy_eul)) THEN
              !
              ! register crossing with latitude: line-segments point Southward
              !
              jcross_lat = jcross_lat+1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              cross_lat_eul_index(jcross_lat,2) = jy_eul
              r_cross_lat(jcross_lat,1) = x(2)
              r_cross_lat(jcross_lat,2) = y(2)
!              write(*,*) "B register crossing with latitude",x(2),y(2),jx_eul,jy_eul
              !
              jy_eul=jy_eul-1
            END IF
            lcontinue=.FALSE.
          ELSE
            !
            !****************************
            !
            ! not same cell integral
            !
            !****************************
            !
            IF (lsame_cell_x) THEN
              ysgn1 = (1+INT(SIGN(1.0D0,y(2)-y(1))))/2 !"1" if y(2)>y(1) else "0"
              ysgn2 = INT(SIGN(1.0D0,y(2)-y(1)))       !"1" if y(2)>y(1) else "-1"
              !
              !*******************************************************************************
              !
              ! there is at least one crossing with latitudes but no crossing with longitudes
              !
              !*******************************************************************************
              !
              yeul   = ygno(jy_eul+ysgn1)
              IF (x(1).EQ.x(2)) THEN
                !
                ! line segment is parallel to longitude (infinite slope)
                !
                xcross = x(1)
              ELSE
                slope  = (y(2)-y(1))/(x(2)-x(1))
                xcross = x_cross_eul_lat(x(1),y(1),yeul,slope)
                !
                ! constrain crossing to be "physically" possible
                !
                xcross = MIN(MAX(xcross,xgno(jx_eul)),xgno(jx_eul+1))
                !
                ! debugging
                !
                IF (xcross.GT.xgno(jx_eul+1).OR.xcross.LT.xgno(jx_eul)) THEN
                  WRITE(*,*) "xcross is out of range",jx,jy
                  WRITE(*,*) "xcross-xgno(jx_eul+1), xcross-xgno(jx_eul))",&
                       xcross-xgno(jx_eul+1), xcross-ygno(jx_eul)
                  STOP
                END IF
              END IF
              xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xcross; yseg(2) = yeul
              jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
              !
              ! prepare for next iteration
              !
              x(0) = x(1); y(0) = y(1); x(1) = xcross; y(1) = yeul; jy_eul = jy_eul+ysgn2
              !
              ! register crossing with latitude
              !
              jcross_lat = jcross_lat+1
              cross_lat_eul_index(jcross_lat,1) = jx_eul
              if (ysgn2>0) then                
                cross_lat_eul_index(jcross_lat,2) = jy_eul
              else
                cross_lat_eul_index(jcross_lat,2) = jy_eul+1
              end if
              r_cross_lat(jcross_lat,1) = xcross
              r_cross_lat(jcross_lat,2) = yeul
!              write(*,*) "C register crossing with latitude",xcross,yeul,jx_eul,cross_lat_eul_index(jcross_lat,2)
            ELSE IF (lsame_cell_y) THEN
              !
              !*******************************************************************************
              !
              ! there is at least one crossing with longitudes but no crossing with latitudes
              !
              !*******************************************************************************
              !
              xsgn1 = (1+INT(SIGN(1.0D0,x(2)-x(1))))/2 !"1" if x(2)>x(1) else "0"
              xsgn2 = INT(SIGN(1.0D0,x(2)-x(1))) !"1" if x(2)>x(1) else "-1"
              xeul   = xgno(jx_eul+xsgn1)
              IF (ABS(x(2)-x(1))<fuzzy_width) THEN
                ! fuzzy crossing
                ycross = 0.5*(y(2)-y(1))
              ELSE
                slope  = (y(2)-y(1))/(x(2)-x(1))
                ycross = y_cross_eul_lon(x(1),y(1),xeul,slope)
              END IF
              !
              ! constrain crossing to be "physically" possible
              !
              ycross = MIN(MAX(ycross,ygno(jy_eul)),ygno(jy_eul+1))
              !
              ! debugging
              !
              IF (ycross.GT.ygno(jy_eul+1).OR.ycross.LT.ygno(jy_eul)) THEN
                WRITE(*,*) "ycross is out of range"
                WRITE(*,*) "jx,jy,jx_eul,jy_eul",jx,jy,jx_eul,jy_eul
                WRITE(*,*) "ycross-ygno(jy_eul+1), ycross-ygno(jy_eul))",&
                     ycross-ygno(jy_eul+1), ycross-ygno(jy_eul)
                STOP
              END IF
              xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xeul; yseg(2) = ycross
              jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
              !
              ! prepare for next iteration
              !
              x(0) = x(1); y(0) = y(1); x(1) = xeul; y(1) = ycross; jx_eul = jx_eul+xsgn2
            ELSE
              !
              !*******************************************************************************
              !
              ! there are crossings with longitude(s) and latitude(s)
              !
              !*******************************************************************************
              ! 
              xsgn1 = (1+INT(SIGN(1.0D0,x(2)-x(1))))/2 !"1" if x(2)>x(1) else "0"
              xsgn2 = (INT(SIGN(1.0D0,x(2)-x(1)))) !"1" if x(2)>x(1) else "0"
              xeul   = xgno(jx_eul+xsgn1) 
              ysgn1 = (1+INT(SIGN(1.0D0,y(2)-y(1))))/2 !"1" if y(2)>y(1) else "0"
              ysgn2 = INT(SIGN(1.0D0,y(2)-y(1)))       !"1" if y(2)>y(1) else "-1"
              yeul   = ygno(jy_eul+ysgn1)
              
              slope  = (y(2)-y(1))/(x(2)-x(1))
              IF (ABS(x(2)-x(1))<fuzzy_width) THEN
                ycross = 0.5*(y(2)-y(1))
              ELSE
                ycross = y_cross_eul_lon(x(1),y(1),xeul,slope)
              END IF
              xcross = x_cross_eul_lat(x(1),y(1),yeul,slope)

              
              IF ((xsgn2>0.AND.xcross.LE.xeul).OR.(xsgn2<0.AND.xcross.GE.xeul)) THEN
                !
                ! cross latitude
                !
                xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xcross; yseg(2) = yeul
                jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
                !
                ! prepare for next iteration
                !
                x(0) = x(1); y(0) = y(1); x(1) = xcross; y(1) = yeul; jy_eul = jy_eul+ysgn2
                !
                ! register crossing with latitude
                !
                jcross_lat = jcross_lat+1
                cross_lat_eul_index(jcross_lat,1) = jx_eul
                if (ysgn2>0) then                
                  cross_lat_eul_index(jcross_lat,2) = jy_eul
                else
                  cross_lat_eul_index(jcross_lat,2) = jy_eul+1
                end if
                r_cross_lat(jcross_lat,1) = xcross
                r_cross_lat(jcross_lat,2) = yeul
!              write(*,*) "D register crossing with latitude",xcross,yeul,jx_eul,cross_lat_eul_index(jcross_lat,2)
              ELSE
                !
                ! cross longitude
                !
                xseg(1) = x(1); yseg(1) = y(1); xseg(2) = xeul; yseg(2) = ycross
                jx_eul_tmp = jx_eul; jy_eul_tmp = jy_eul; 
                !
                ! prepare for next iteration
                !
                x(0) = x(1); y(0) = y(1); x(1) = xeul; y(1) = ycross; jx_eul = jx_eul+xsgn2
              END IF
              
            END IF
          END IF
          !
          ! register line-segment (don't register line-segment if outside of panel)
          !
          if (jx_eul_tmp>=jx_min.AND.jy_eul_tmp>=jy_min.AND.&
               jx_eul_tmp<=jx_max-1.AND.jy_eul_tmp<=jy_max-1) then
            jsegment=jsegment+1
            weights_eul_index(jsegment,1) = jx_eul_tmp
            weights_eul_index(jsegment,2) = jy_eul_tmp

!            call get_weights_exact(lexact_horizontal_line_integrals, weights_tmp,xseg,yseg,&
!                 nreconstruction,ngauss,gauss_weights,abscissae)

             
            !
            ! debugging phl
            !
!            write(*,*) "start test"
!            xseg(1) = 4.3660942908512038D-002; yseg(1) = 0.83909966789439261D0
!            xseg(2) = 8.7488663525923979D-002; yseg(2) = 0.83909963117727981D0
!!## weight is    2.8010300820513669E-002
!
!            call get_weights_exact(.false.,&
!                 weights(jsegment,1:nreconstruction),&
!                 xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
!            write(*,*) "weights is",weights(jsegment,1)
!            jsegment = jsegment+1
!!#
!            xseg(1)= 8.7488663525923979D-002; yseg(1) =   0.83909963117727981D0
!            xseg(2) = 4.3660942908512038D-002; yseg(2) =  0.83909963117727981D0
!
!
!
!            call get_weights_exact(.false.,&
!                 weights(jsegment,1:nreconstruction),&
!                 xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
!            write(*,*) "weights is",weights(jsegment,1)
!            write(*,*) "sum=",weights(jsegment-1,1)+weights(jsegment,1)
!            write(*,*) "end test"
!            stop
!
            call get_weights_exact(lexact_horizontal_line_integrals.AND.ABS(yseg(2)-yseg(1))<tiny,&
                 weights(jsegment,:),&
                 xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
            if (ldbg) then
               write(*,*) "from inside side-integral"
               write(*,*) "line-integral from/to ",xseg(1),yseg(1),xseg(2),yseg(2)
               write(*,*) "weight is ",weights(jsegment,1)
            end if



!old            call get_weights_gauss(weights(jsegment,1:nreconstruction),&
!old                 xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
          ELSE
            !
            ! segment outside of panel
            !
          END IF
          
        END DO
        side_count = side_count+1
      END DO
    END IF
  end subroutine side_integral
 

  real (kind=real_kind) function compute_slope(x,y)
    implicit none
    real (kind=real_kind), dimension(:), intent(in) :: x,y !dimension(2)
    if (fuzzy(ABS(x(2)-x(1)),fuzzy_width)>0) THEN
      compute_slope = (y(2)-y(1))/(x(2)-x(1))
    else
      compute_slope = bignum
    end if
  end function compute_slope

  real (kind=real_kind) function y_cross_eul_lon(x,y,xeul,slope)
    implicit none
    real (kind=real_kind), intent(in) :: x,y
    real (kind=real_kind)              , intent(in) :: xeul,slope
    !    
    ! line: y=a*x+b
    !
    real (kind=real_kind) :: a,b
    
    b = y-slope*x 
    y_cross_eul_lon = slope*xeul+b
  end function y_cross_eul_lon

  real (kind=real_kind) function x_cross_eul_lat(x,y,yeul,slope)
    implicit none
    real (kind=real_kind), intent(in) :: x,y
    real (kind=real_kind)              , intent(in) :: yeul,slope

    if (fuzzy(ABS(slope),fuzzy_width)>0) THEN
        x_cross_eul_lat = x+(yeul-y)/slope
    ELSE
      x_cross_eul_lat = bignum
    END IF
  end function x_cross_eul_lat

  subroutine get_weights_exact(lexact_horizontal_line_integrals,weights,xseg,yseg,nreconstruction,&
       ngauss,gauss_weights,abscissae)
    use fvm_analytic_mod, only: I_00, I_10, I_01, I_20, I_02, I_11
    use dimensions_mod, only: ip_dbg, io_dbg, ie_dbg !dbg
    implicit none
    logical, intent(in) :: lexact_horizontal_line_integrals
    integer (kind=int_kind), intent(in) :: nreconstruction, ngauss
    real (kind=real_kind), intent(out) :: weights(:)
    real (kind=real_kind), dimension(:), intent(in) :: gauss_weights, abscissae !dimension(ngauss)
    
    
    real (kind=real_kind), dimension(:), intent(in) :: xseg,yseg !dimension(2)
    !
    ! compute weights
    !
    real (kind=real_kind) :: tmp,slope,b,integral,dx2,xc
    integer (kind=int_kind) :: i

    !
    ! dbg start
    !
!    if (io_dbg.and.ie_dbg==9) then       
    if (io_dbg) then       
       write(ip_dbg+40,*) xseg(1),yseg(1),xseg(2)-xseg(1),yseg(2)-yseg(1)
       write(ip_dbg+40,*) "   "
    end if

    !
    ! dbg end
    !



    if(lexact_horizontal_line_integrals) then
      weights(1) = ((I_00(xseg(2),yseg(2))-I_00(xseg(1),yseg(1))))
      if (ABS(weights(1))>1.0) THEN
        WRITE(*,*) "1 exact weights(jsegment)",weights(1),xseg,yseg
        stop
      end if
      if (nreconstruction>1) then
         weights(2) = ((I_10(xseg(2),yseg(2))-I_10(xseg(1),yseg(1))))
         weights(3) = ((I_01(xseg(2),yseg(2))-I_01(xseg(1),yseg(1))))
      endif
      if (nreconstruction>3) then
         weights(4) = ((I_20(xseg(2),yseg(2))-I_20(xseg(1),yseg(1))))
         weights(5) = ((I_02(xseg(2),yseg(2))-I_02(xseg(1),yseg(1))))
         weights(6) = ((I_11(xseg(2),yseg(2))-I_11(xseg(1),yseg(1))))
      endif
    else
      call get_weights_gauss(weights,xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
    endif
  end subroutine get_weights_exact



  subroutine get_weights_gauss(weights,xseg,yseg,nreconstruction,ngauss,gauss_weights,abscissae)
    use fvm_analytic_mod, only: I_00, I_10, I_01, I_20, I_02, I_11
    
    implicit none
    integer (kind=int_kind), intent(in) :: nreconstruction,ngauss
    real (kind=real_kind), intent(out) :: weights(:)
    real (kind=real_kind), dimension(2     ), intent(in) :: xseg,yseg
    real (kind=real_kind) :: slope
    !
    ! compute weights
    !
    !
    ! for Gaussian quadrature
    !
    real (kind=real_kind), dimension(ngauss), intent(in) :: gauss_weights, abscissae

    ! if line-segment parallel to x or y use exact formulaes else use qudrature
    !
    real (kind=real_kind) :: tmp,b,integral,dx2,xc,x,y
    integer (kind=int_kind) :: i

!    if (fuzzy(abs(xseg(1) -xseg(2)),fuzzy_width)==0)then
    if (xseg(1).EQ.xseg(2))then
      weights = 0.0D0
    else
      slope    = (yseg(2)-yseg(1))/(xseg(2)-xseg(1))
      b        = yseg(1)-slope*xseg(1)
      dx2      = 0.5D0*(xseg(2)-xseg(1))
      xc       = 0.5D0*(xseg(1)+xseg(2))
      integral = 0.0D0
      do i=1,ngauss
        x        = xc+abscissae(i)*dx2
        y        = slope*x+b
        integral = integral+gauss_weights(i)*F_00(x,y)
      enddo
      weights(1) = integral*dx2  
      if (nreconstruction>1) then
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_10(x,y)
        enddo
        weights(2) = integral*dx2  
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_01(x,y)
        enddo
        weights(3) = integral*dx2  
      endif
      if (nreconstruction>3) then
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_20(x,y)
        enddo
        weights(4) = integral*dx2  
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_02(x,y)
        enddo
        weights(5) = integral*dx2  
        integral = 0.0D0
        do i=1,ngauss
          x        = xc+abscissae(i)*dx2
          y        = slope*x+b
          integral = integral+gauss_weights(i)*F_11(x,y)
        enddo
        weights(6) = integral*dx2  
      endif
    end if
  end subroutine get_weights_gauss

  real (kind=real_kind) function F_00(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp
    !
    x = x_in
    y = y_in
    F_00 =y/((1.0D0+x*x)*SQRT(1.0D0+x*x+y*y))
  end function F_00

  real (kind=real_kind) function F_10(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp

    x = x_in
    y = y_in

    F_10 =x*y/((1.0D0+x*x)*SQRT(1.0D0+x*x+y*y))
  end function F_10

  real (kind=real_kind) function F_01(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y

    x = x_in
    y = y_in

    F_01 =-1.0D0/(SQRT(1.0D0+x*x+y*y))
  end function F_01

  real (kind=real_kind) function F_20(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp

    x = x_in
    y = y_in

    F_20 =x*x*y/((1.0D0+x*x)*SQRT(1.0D0+x*x+y*y))
  end function F_20

  real (kind=real_kind) function F_02(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,alpha, tmp

    x = x_in
    y = y_in

    alpha = ATAN(x)
!     F_02 =-y/SQRT(1.0D0+x*x+y*y)+ASINH(y*COS(alpha))
    tmp=y*COS(alpha)
    F_02 =-y/SQRT(1.0D0+x*x+y*y)+log(tmp+sqrt(tmp*tmp+1))
    
    !
    ! cos(alpha) = 1/sqrt(1+x*x)
    !
  end function F_02

  real (kind=real_kind) function F_11(x_in,y_in)
    implicit none
    real (kind=real_kind), intent(in) :: x_in,y_in
    real (kind=real_kind) :: x,y,tmp

    x = x_in
    y = y_in

    F_11 =-x/(SQRT(1.0D0+x*x+y*y))
  end function F_11

  subroutine which_eul_cell(x,j_eul,gno)
    implicit none
    integer (kind=int_kind)                               , intent(inout) :: j_eul
    real (kind=real_kind), dimension(:)                    , intent(in)    :: x !dimension(3)
    real (kind=real_kind), dimension(-nhe:), intent(in)    :: gno ! dimension(-nhe:nc+2+nhe)
!    real (kind=real_kind), intent(in)    :: eps
    
    real (kind=real_kind) :: d1,d2,d3,d1p1
    logical                 :: lcontinue
    integer :: iter
    
    lcontinue = .TRUE.
    iter = 0 

    DO WHILE (lcontinue)
      iter = iter+1 
      IF (x(1).GE.gno(j_eul).AND.x(1).LT.gno(j_eul+1)) THEN
        lcontinue = .FALSE.
        !
        ! special case when x(1) is on top of grid line
        !
        IF (x(1).EQ.gno(j_eul)) THEN
          !
          ! x(1) is on top of gno(J_eul)
          !
          IF (x(2).GT.gno(j_eul)) THEN
            j_eul = j_eul
          ELSE IF (x(2).LT.gno(j_eul)) THEN
            j_eul = j_eul-1
          ELSE
            !
            ! x(2) is on gno(j_eul) grid line; need x(3) to determine Eulerian cell 
            !
            IF (x(3).GT.gno(j_eul)) THEN
              !
              ! x(3) to the right
              !
              j_eul = j_eul
            ELSE IF (x(3).LT.gno(j_eul)) THEN
              !
              ! x(3) to the left
              !
              j_eul = j_eul-1
            ELSE
              WRITE(*,*) "inconsistent cell: x(1)=x(2)=x(3)",x(1),x(2),x(3)
              STOP
            END IF
          END IF
        END IF
      ELSE
        ! 
        ! searching - prepare for next iteration
        !
        IF (x(1).GE.gno(j_eul+1)) THEN
          j_eul = j_eul + 1
        ELSE
          !
          ! x(1).LT.gno(j_eul)
          !
          j_eul = j_eul - 1
        END IF
      END IF
      IF (iter>1000.OR.j_eul<-nhe.OR.j_eul>nc+2+nhe) THEN
        WRITE(*,*) "search is which_eul_cell not converging!", iter, j_eul,nhe,nc+2+nhe
        WRITE(*,*) "gno", gno(nc+2+nhe), gno(-nhe)
        write(*,*) gno
        STOP
      END IF
    END DO
  END subroutine which_eul_cell


  subroutine truncate_vertex(x,j_eul,gno)
    implicit none
    integer (kind=int_kind)                               , intent(inout) :: j_eul
    real (kind=real_kind)                    , intent(inout)    :: x
    real (kind=real_kind), dimension(-nhe:), intent(in)    :: gno !dimension(-nhe:nc+2+nhe)
!    real (kind=real_kind), intent(in)    :: eps
    
    logical                 :: lcontinue
    integer :: iter, xsgn
    real (kind=real_kind) :: dist,dist_new,tmp
    
    lcontinue = .TRUE.
    iter = 0 
    dist = bignum
    xsgn     = INT(SIGN(1.0D00,x-gno(j_eul)))
    
    DO WHILE (lcontinue)
      if ((j_eul<-nhe) .or. (j_eul>nc+2+nhe)) then
        write(*,*) 'somthing is wrong', j_eul, -nhe,nc+2+nhe, iter
        stop
      endif
      iter     = iter+1 
      tmp      = x-gno(j_eul)
      dist_new = ABS(tmp)
      IF (dist_new>dist) THEN
        lcontinue = .FALSE.
!      ELSE IF (ABS(tmp)<1.0E-11) THEN
      ELSE IF (ABS(tmp)<1.0E-9) THEN
!      ELSE IF (ABS(tmp)<tiny) THEN
!      ELSE IF (ABS(tmp)<1.0E-4) THEN
        x = gno(j_eul)
        lcontinue = .FALSE.
      ELSE
        j_eul = j_eul+xsgn
        dist = dist_new
      END IF
      IF (iter>100) THEN
        WRITE(*,*) "truncate vertex not converging"
        STOP
      END IF
    END DO
  END subroutine truncate_vertex




!********************************************************************************
!
! Gauss-Legendre quadrature
!
! Tabulated values
!
!********************************************************************************
subroutine gauss_points(n,weights,points)
  implicit none
  integer (kind=int_kind)           , intent(in ) :: n
  real (kind=real_kind), dimension(:), intent(out) :: weights, points !dimension(n)
  
  select case (n)
!    CASE(1)
!       abscissae(1) = 0.0D0
!       weights(1)   = 2.0D0
  case(2)
     points(1)    = -sqrt(1.0D0/3.0D0)
     points(2)    =  sqrt(1.0D0/3.0D0)
     weights(1)   =  1.0D0
     weights(2)   =  1.0D0
  case(3)
     points(1)    = -0.774596669241483377035853079956D0
     points(2)    =  0.0D0
     points(3)    =  0.774596669241483377035853079956D0
     weights(1)   =  0.555555555555555555555555555556D0
     weights(2)   =  0.888888888888888888888888888889D0
     weights(3)   =  0.555555555555555555555555555556D0
  case(4)
     points(1)    = -0.861136311594052575223946488893D0
     points(2)    = -0.339981043584856264802665659103D0
     points(3)    =  0.339981043584856264802665659103D0
     points(4)    =  0.861136311594052575223946488893D0
     weights(1)   =  0.347854845137453857373063949222D0
     weights(2)   =  0.652145154862546142626936050778D0 
     weights(3)   =  0.652145154862546142626936050778D0 
     weights(4)   =  0.347854845137453857373063949222D0      
  case(5)
     points(1)    = -(1.0D0/3.0D0)*sqrt(5.0D0+2.0D0*sqrt(10.0D0/7.0D0))
     points(2)    = -(1.0D0/3.0D0)*sqrt(5.0D0-2.0D0*sqrt(10.0D0/7.0D0))
     points(3)    =  0.0D0
     points(4)    =  (1.0D0/3.0D0)*sqrt(5.0D0-2.0D0*sqrt(10.0D0/7.0D0))
     points(5)    =  (1.0D0/3.0D0)*sqrt(5.0D0+2.0D0*sqrt(10.0D0/7.0D0))
     weights(1)   = (322.0D0-13.0D0*sqrt(70.0D0))/900.0D0
     weights(2)   = (322.0D0+13.0D0*sqrt(70.0D0))/900.0D0
     weights(3)   = 128.0D0/225.0D0
     weights(4)   = (322.0D0+13.0D0*sqrt(70.0D0))/900.0D0
     weights(5)   = (322.0D0-13.0D0*sqrt(70.0D0))/900.0D0
  case default
     write(*,*) 'n out of range in glwp of module gll. n=',n
     write(*,*) '0<n<5'
     stop
  end select

end subroutine gauss_points

!------------------------------------------------------------------------------
! FUNCTION SIGNUM_FUZZY
!
! Description:
!   Gives the sign of the given real number, returning zero if x is within 
!     a small amount from zero.
!------------------------------------------------------------------------------
  function signum_fuzzy(x)
    implicit none

    real (kind=real_kind) :: signum_fuzzy
    real (kind=real_kind) :: x

    IF (x > fuzzy_width) THEN
      signum_fuzzy = 1.0D0
    ELSEIF (x < fuzzy_width) THEN
      signum_fuzzy = -1.0D0
    ELSE
      signum_fuzzy = 0.0D0
    ENDIF
  end function

  function fuzzy(x,epsilon)
    implicit none

    integer (kind=int_kind) :: fuzzy
    real (kind=real_kind), intent(in) :: epsilon
    real (kind=real_kind) :: x

    IF (ABS(x)<epsilon) THEN
      fuzzy = 0
    ELSE IF (x >epsilon) THEN
      fuzzy = 1
    ELSE !IF (x < fuzzy_width) THEN
      fuzzy = -1
    ENDIF
  end function


  real (kind=real_kind) function area(xseg,yseg,nvertex)
    implicit none
    integer (kind=int_kind)                    , intent(in):: nvertex
    real (kind=real_kind)  , dimension(nvertex), intent(in):: xseg,yseg
    !
    integer (kind=int_kind):: i
    area = 0.0
    do i=1,nvertex-1
       !
       ! factor 0.5 missing for area computation, however,
       ! we are only interested in the sign of the "area"
       !
       area = area - (yseg(i+1)-yseg(i))*(xseg(i+1)+xseg(i))
!       area = area - (yseg(i+1)-yseg(i))*(xseg(i+1)+xseg(i))
    end do
    area = area - (yseg(1)-yseg(nvertex))*(xseg(1)+xseg(nvertex))
    if (abs(area)< tiny) area = 0.0
    if (ldbg) write(*,*) "area is ",area
  end function area

  !
  ! below is departure point code
  !
  !
  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE FVM_MESH_DEP--------------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 30.October 2011                                          !
  ! DESCRIPTION: Calculates the departure mesh                                        !
  !                                                                                   !
  ! CALLS: solidbody or boomerang, cart2cubedspherexy, analytical_function            !
  ! INTPUT/OUTPUT:                                                                    !
  !        fvm   ...  structure                                                       !
  ! INPUT: nstep   ... actual step                                                    !
  !-----------------------------------------------------------------------------------!
  subroutine fvm_mesh_dep(elem, deriv, fvm, dt, tl,klev)
    use coordinate_systems_mod, only : cartesian2D_t
    use element_mod, only : element_t
    use fvm_control_volume_mod, only: fvm_struct, fvm_supercycling
    use time_mod, only : timelevel_t, time_at
    use parallel_mod, only : haltmp
    use control_mod, only : test_cfldep
    
    use derivative_mod, only : derivative_t
    use fvm_bsp_mod, only: boomerang, solidbody
    
    implicit none
    type (derivative_t)  , intent(in) :: deriv
    type (fvm_struct), intent(inout)   :: fvm
    type (timelevel_t),intent(in)        :: tl
    integer,intent(in)                   :: klev
    
    type (element_t), intent(inout)      :: elem
    real (kind=real_kind)                :: time,dt, dx, dy, cflx, cfly, maxcflx, maxcfly  
    integer                              :: i,j
    type (cartesian2D_t)                 :: dcart(nc+1,nc+1)
    
    
!phl    ! for the benchmark test, use more accurate departure point creation
    if (fvm_ideal_test == IDEAL_TEST_OFF) then
  ! for given velocities in the element structure
      call fvm_dep_from_gll(elem, deriv, fvm%asphere,fvm%dsphere,dt*fvm_supercycling,tl,klev)    
    else
!phl    !CE: define new mesh for fvm fvm on an equal angular grid
!phl    ! go from alpha,beta -> cartesian xy on cube -> lat lon on the sphere
      do j = 1, nc+1
        do i = 1, nc+1               
          if (fvm_test_type == IDEAL_TEST_BOOMERANG) then
             !
             ! broken
             ! 
            call boomerang(fvm%asphere(i,j), fvm%dsphere(i,j,klev),tl%nstep)
          else if (fvm_test_type == IDEAL_TEST_SOLIDBODY) then
            call solidbody(fvm%asphere(i,j), fvm%dsphere(i,j,klev))
          else
            call haltmp("Unknown FVM ideal test type")
          end if
        end do
      end do
    end if
    
    if (test_cfldep) then
       call check_departurecell(fvm,klev) 
    endif
    
  end subroutine fvm_mesh_dep
  
  ! subroutine fvm_dep_from_gll_iter(elem, deriv, fvm, dt, tl, klev)
  !   use coordinate_systems_mod, only : cartesian2D_t
  !   use element_mod, only : element_t
  !   use fvm_control_volume_mod, only: fvm_struct
  !   use time_mod, only : timelevel_t, time_at
  !   use parallel_mod, only : haltmp
  !   use control_mod, only : test_cfldep
  ! 
  !   use derivative_mod, only : derivative_t
  ! 
  !   implicit none
  !   type (element_t), intent(in)          :: elem
  !   type (derivative_t)  , intent(in) :: deriv
  !   type (spherical_polar_t),intent(in)   :: asphere(nc+1,nc+1)
  !   type (spherical_polar_t),intent(out)  :: dsphere(nc+1,nc+1)
  !   type (timelevel_t),intent(in)        :: tl
  !   real (kind=real_kind),intent(in)      :: dt
  !   integer,intent(in)                   :: klev
  !   
  !   real (kind=real_kind)                 :: uxyz_gll(np,np,3),uxyz(nc+1,nc+1,3)
  !   real (kind=real_kind)                 :: un,ue,ur,clon,clat,slon,slat
  !   type(cartesian3D_t)                   :: acart
  !   integer :: i,j,k
  !   
  !   real (kind=real_kind)  :: vstar(np,np,2)
  !   
  !    ! convert velocity from lat/lon to cartesian 3D
  !    vstar = elem%derived%vstar(:,:,:,klev)
  !     
  !    do i=1,np
  !      do j=1,np
  !         clon = cos(elem%spherep(i,j)%lon)
  !         slon = sin(elem%spherep(i,j)%lon)
  !         clat = cos(elem%spherep(i,j)%lat)
  !         slat = sin(elem%spherep(i,j)%lat)
  ! 
  !         ur = 0
  !         ue = vstar(i,j,1) 
  !         un = vstar(i,j,2)
  ! 
  !         uxyz_gll(i,j,1)= clon*clat*ur - clon*slat*un - slon*ue
  !         uxyz_gll(i,j,2)= slon*clon*ur - slon*slat*un + clon*ue
  !         uxyz_gll(i,j,3)= slat          *ur + clat          *un
  !     
  !      enddo
  !    enddo
  !    ! interpolate velocity to fvm nodes
  !    do i=1,3
  !       uxyz(:,:,i)=interpolate_gll2fvm_corners(uxyz_gll(:,:,i),deriv)
  !    end do 
  !    ! compute departure point 
  !    ! crude, 1st order accurate approximation.  to be improved
  !    do i=1,nc+1
  !       do j=1,nc+1
  !          ! note: asphere()%r=1, so we need to convert velocity to radians/sec:
  !          acart = change_coordinates(asphere(i,j))  
  !          acart%x = acart%x - dt*uxyz(i,j,1)/rearth
  !          acart%y = acart%y - dt*uxyz(i,j,2)/rearth
  !          acart%z = acart%z - dt*uxyz(i,j,3)/rearth
  !          dsphere(i,j) = change_coordinates(acart)
  !          dsphere(i,j)%r = asphere(i,j)%r
  !       enddo
  !    enddo
  !   
  ! end subroutine fvm_dep_from_gll_iter
  
  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE FVM_DEP_FROM_GLL----------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR 14. December 2011                            !
  ! DESCRIPTION: calculates the deparute grid for fvm coming from the gll points      !
  !                                                                                   !
  ! CALLS: 
  ! INPUT: 
  !        
  ! OUTPUT: 
  !-----------------------------------------------------------------------------------!
  subroutine fvm_dep_from_gll(elem, deriv, asphere,dsphere,dt,tl,klev)
    use physical_constants, only : DD_PI, rearth
    use coordinate_systems_mod, only : spherical_polar_t, cartesian3D_t, change_coordinates
    use time_mod, only : timelevel_t
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t, interpolate_gll2fvm_corners
    
    implicit none
    type (element_t), intent(in)          :: elem
    type (derivative_t)  , intent(in) :: deriv
    type (spherical_polar_t),intent(in)   :: asphere(nc+1,nc+1)
    type (spherical_polar_t),intent(out)  :: dsphere(-1:nc+3,-1:nc+3,nlev)
    type (timelevel_t),intent(in)        :: tl
    real (kind=real_kind),intent(in)      :: dt
    integer,intent(in)                   :: klev
    
    real (kind=real_kind)                 :: uxyz_gll(np,np,3),uxyz(nc+1,nc+1,3)
    real (kind=real_kind)                 :: un,ue,ur,clon,clat,slon,slat
    type(cartesian3D_t)                   :: acart
    integer :: i,j,k
    
    real (kind=real_kind)  :: vstar(np,np,2)
    
    ! convert velocity from lat/lon to cartesian 3D
    vstar = elem%derived%vstar(:,:,:,klev)
    
    do i=1,np
       do j=1,np
          clon = cos(elem%spherep(i,j)%lon)
          slon = sin(elem%spherep(i,j)%lon)
          clat = cos(elem%spherep(i,j)%lat)
          slat = sin(elem%spherep(i,j)%lat)
          
          ur = 0
          ue = vstar(i,j,1) 
          un = vstar(i,j,2)
          
          uxyz_gll(i,j,1)= clon*clat*ur - clon*slat*un - slon*ue
          uxyz_gll(i,j,2)= slon*clon*ur - slon*slat*un + clon*ue
          uxyz_gll(i,j,3)= slat          *ur + clat          *un
          
       enddo
    enddo
    ! interpolate velocity to fvm nodes
    do i=1,3
       uxyz(:,:,i)=interpolate_gll2fvm_corners(uxyz_gll(:,:,i),deriv)
    end do
    ! compute departure point 
    ! crude, 1st order accurate approximation.  to be improved
    do i=1,nc+1
       do j=1,nc+1
          ! note: asphere()%r=1, so we need to convert velocity to radians/sec:
          acart = change_coordinates(asphere(i,j))  
          acart%x = acart%x - dt*uxyz(i,j,1)/rearth
          acart%y = acart%y - dt*uxyz(i,j,2)/rearth
          acart%z = acart%z - dt*uxyz(i,j,3)/rearth
          dsphere(i,j,klev) = change_coordinates(acart)
          dsphere(i,j,klev)%r = asphere(i,j)%r
       enddo
    enddo
    
  end subroutine fvm_dep_from_gll
  !END SUBROUTINE FVM_DEP_FROM_GLL------------------------------------------CE-for FVM!
  
  !SUBROUTINE CHECK_DEPARTURECELL-------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
  ! DESCRIPTION: Check if the departure cell is degenerated or if the                 !
  !              CFL number > 1.0D0                                                   !
  !         IF THE CHECK FAILS, THE PROGRAM WILL STOP                                 !
  !                                                                                   !
  ! INPUT:  fvm  ... fvm structur                                                 ! 
  !         klev   ... Level (vertical)                                               !
  !-----------------------------------------------------------------------------------!
  subroutine check_departurecell(fvm,klev)
    use coordinate_systems_mod, only: cartesian2D_t, cart2cubedspherexy, spherical_to_cart
    use fvm_control_volume_mod, only:  fvm_struct                                         
    implicit none
    
    type (fvm_struct), intent(inout)   :: fvm
    integer,intent(in)                   :: klev
    
    type (cartesian2D_t)                 :: dcart(nc+1,nc+1)
    real (kind=real_kind)                :: cflx, cfly, maxcflx, maxcfly  
    integer                              :: i,j
    logical                              :: crossline
    
    ! calculate xy Cartesian on the cube of departure points on the corresponding face  
    ! the calculation of dcart is also done in fvm_line_integrals_mod.F90, one could 
    ! save the points here...
    ! BUT: this subroutine should only be used for test reasons
    do j=1,nc+1
       do i=1,nc+1               
          call cart2cubedspherexy(spherical_to_cart(fvm%dsphere(i,j,klev)),&
               fvm%faceno,dcart(i,j))              
       end do
    end do
    cflx=0.0D0
    cfly=0.0D0
    fvm%maxcfl(1,klev)=cflx
    fvm%maxcfl(2,klev)=cfly
    crossline=.FALSE.
    do j=1,nc
       do i=1,nc  
          ! equidistant mesh in alpha/beta coordinates
          cflx=abs((atan(dcart(i,j)%x)-atan(fvm%acartx(i)))/fvm%dalpha)
          cfly=abs((atan(dcart(i,j)%y)-atan(fvm%acarty(j)))/fvm%dbeta)
          if(cflx>fvm%maxcfl(1,klev)) fvm%maxcfl(1,klev)=cflx
          if(cfly>fvm%maxcfl(2,klev)) fvm%maxcfl(2,klev)=cfly
          
          !one could stop immediately here, if crossline=.TRUE., but want to calculate 
          ! all CFL first
          call check_lines_cross(dcart(i,j),dcart(i+1,j),dcart(i,j+1),dcart(i+1,j+1),crossline)  
          call check_lines_cross(dcart(i,j),dcart(i,j+1),dcart(i+1,j),dcart(i+1,j+1),crossline)
       end do
    end do
    ! nodes on the north and east boundary
    do i=1,nc+1  
       cflx=abs((atan(dcart(i,nc+1)%x)-atan(fvm%acartx(i)))/fvm%dalpha)
       cfly=abs((atan(dcart(nc+1,j)%y)-atan(fvm%acarty(j)))/fvm%dbeta)
       if(cflx>fvm%maxcfl(1,klev)) fvm%maxcfl(1,klev)=cflx
       if(cfly>fvm%maxcfl(2,klev)) fvm%maxcfl(2,klev)=cfly
    end do
    
    if (fvm%maxcfl(1,klev) > 1.0D0 .OR. fvm%maxcfl(2,klev) > 1.0D0) then
       write(*,*) "Error in fvm_mod.F90: CFL number too big"
       write(*,*) "CFL has to be < 1.0D0"
       write(*,*) "Choose a smaller time step!"
       write(*,*) "max CFL in this element: maxcflx", fvm%maxcfl(1,klev), "maxcfly", fvm%maxcfl(2,klev) 
       write(*,*) "klev=",klev
       STOP "Exit program!"
    endif
    
    if (crossline) then
       write(*,*) "FATAL Error in fvm_mod.F90: departure cell is degenerated!"
       write(*,*) "Choose a smaller time step! max CFL in this element: maxcflx", fvm%maxcfl(1,klev), "maxcfly", fvm%maxcfl(2,klev)
       STOP 'Exit program!'
    endif
  end subroutine check_departurecell
  
  
  ! see, e.g., http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
  subroutine check_lines_cross(p1,p2,q1,q2,crossline)
    use coordinate_systems_mod, only: cartesian2D_t
    implicit none
    type(cartesian2D_t), intent (in)  :: p1,p2,q1,q2
    logical            , intent(inout)  :: crossline
    !
    ! local workspace
    !
    REAL (kind=real_kind)    :: dn,tp,tq
    
    dn = (q2%y-q1%y)*(p2%x-p1%x)-(q2%x-q1%x)*(p2%y-p1%y)
    
    if (abs(dn)>1.0D-12) then
       tp = ((q2%x-q1%x)*(p1%y-q1%y)-(q2%y-q1%y)*(p1%x-q1%x))/dn
       tq = ((p2%x-p1%x)*(p1%y-q1%y)-(p2%y-p1%y)*(p1%x-q1%x))/dn
       ! implement the next two lines faster!!!!!!
       if (tp>=0.0D0 .and. tp<=1.0D00 .and. &
            tq>=0.0D0 .and. tq<=1.0D00) then
          crossline=.TRUE.
       endif
    endif
  end subroutine check_lines_cross
  
  !SUBROUTINE CHECK_DEPARTURECELLEST----------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 17.October 2011                                          !
  ! DESCRIPTION: Check if the departure cell is degenerated or if the                 !
  !              CFL number > 1.0D0                                                   !
  !         IF THE CHECK FAILS, THE PROGRAM WILL STOP                                 !
  !                                                                                   !
  ! INPUT:  acartx, acarty  ... arrival grid coordinates in the element               ! 
  !         dcart           ... departure grid coordinates in the element             !                                         !
  !                                                                                   !
  !OLD AND DIFFERENT CHECK FOR THE SHAPE OF THE DEPARTURE CELL------------------------!
  subroutine check_departurecellest(acartx,acarty,dcart)
    use coordinate_systems_mod, only: cartesian2D_t
    implicit none
    
    real (kind=real_kind), intent(in)    :: acartx(-nhe:nc+2+nhe), acarty(-nhe:nc+2+nhe)
    type (cartesian2D_t), intent(in)     :: dcart(nc+1,nc+1)
    real (kind=real_kind)                :: dx, dy, cflx, cfly, maxcflx, maxcfly  
    integer                              :: i,j
    logical                              :: crossline
    
    cflx=0.0D0
    cfly=0.0D0
    maxcflx=cflx
    maxcfly=cfly
    do j=1,nc
       do i=1,nc  
          !estimate (rough) CFL number
          dx=acartx(i)-acartx(i+1)
          dy=acarty(j)-acarty(j+1)
          
          cflx=abs((dcart(i,j)%x-acartx(i))/dx)
          cfly=abs((dcart(i,j)%y-acarty(j))/dy)
          
          if(cflx>maxcflx) maxcflx=cflx
          if(cfly>maxcfly) maxcfly=cfly
          call check_lines_cross(dcart(i,j),dcart(i+1,j),dcart(i,j+1),dcart(i+1,j+1),crossline)  
          call check_lines_cross(dcart(i,j),dcart(i,j+1),dcart(i+1,j),dcart(i+1,j+1),crossline)
       end do
    end do
    
    do i=1,nc  
       !estimate (rough) CFL number
       dx=acartx(i)-acartx(i+1)
       dy=acarty(j)-acarty(j+1)
       
       cflx=abs((dcart(i,nc+1)%x-acartx(i))/dx)
       cfly=abs((dcart(nc+1,j)%y-acarty(j))/dy)
       
       if(cflx>maxcflx) maxcflx=cflx
       if(cfly>maxcfly) maxcfly=cfly
    end do
    
    if (maxcflx > 1.0D0 .OR. maxcfly > 1.0D0) then
       write(*,*) "FATAL Error in fvm_line_integrals.F90: CFL number too big"
       write(*,*) "CFL has to be < 1.0D0"
       write(*,*) "Choose a smaller time step!"
       write(*,*) "max CFL in this element: maxcflx", maxcflx, "maxcfly", maxcfly 
       STOP "Exit program!"
    endif
    
    if (crossline) then
       write(*,*) "FATAL Error in fvm_line_integrals.F90: departure cell is degenerated!"
       write(*,*) "Choose a smaller time step!"
       write(*,*) "max CFL in this element: maxcflx", maxcflx, "maxcfly", maxcfly 
       STOP "Exit program!"
    endif
  end subroutine check_departurecellest
  
  
  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE FVM_MCGREGOR--------------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 09. January 2012                                         !
  ! DESCRIPTION: ! using McGregor AMS 1993 scheme: Economical Determination of        !
  !                Departure Points for Semi-Lagrangian Models                        !
  !                                                                                   !
  ! CALLS: 
  ! INPUT: 
  !        
  ! OUTPUT: 
  !-----------------------------------------------------------------------------------!
  subroutine fvm_mcgregor(elem, deriv, dt_fvm, vhat, vstar,order)
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t, gradient_sphere_routine
    use derivative_mod, only : ugradv_sphere_routine, vorticity_sphere_routine 
    implicit none
    
    type (element_t), intent(in)                                :: elem
    type (derivative_t), intent(in)                             :: deriv      ! derivative struct
    real (kind=real_kind), intent(in)                           :: dt_fvm
    real (kind=real_kind), dimension(np,np,2), intent(inout)    :: vstar
    real (kind=real_kind), dimension(np,np,2), intent(in)       :: vhat
    
    integer, intent(in)                                         :: order
    
    integer                                            :: i
    real (kind=real_kind), dimension(np,np,2)          :: ugradv
    real (kind=real_kind)                              :: timetaylor
    
    ugradv=vstar
    timetaylor=1
    do i=1,order
       !     tmp = 0.5D0*(vgradv(:,:,1)**2 + vgradv(:,:,2)**2) 
       ! 
       !     gradvstar = gradient_sphere(tmp,deriv,elem%Dinv)    ! scalar -> latlon vector
       !     tmp = vorticity_sphere(vgradv,deriv,elem)                 ! latlon vector -> scalar 
       ! 
       !     ! v \nabla v expressed through gradient of mean scalar and relative velocity
       !     ! see (17) in Williamson et.al. 1992, JCP 102 (211-224)
       !     vgradv(:,:,1)= gradvstar(:,:,1) - vstarold(:,:,2)*tmp(:,:)
       !     vgradv(:,:,2)= gradvstar(:,:,2) + vstarold(:,:,1)*tmp(:,:)
       !     
       !     timetaylor=-timetaylor*dt_fvm/(i+1)
       !     vstar=vstar + timetaylor*vgradv
       
       call ugradv_sphere_routine(vhat,ugradv,deriv,elem,ugradv)
       timetaylor=-timetaylor*dt_fvm/(i+1)
       
       vstar=vstar + timetaylor*ugradv  
    end do
  end subroutine fvm_mcgregor
  !END SUBROUTINE FVM_MCGREGOR----------------------------------------------CE-for FVM!
  
  
  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE FVM_MCGREGORDSS-----------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, 26. May 2012                                             !
  ! DESCRIPTION: ! using McGregor AMS 1993 scheme: Economical Determination of        !
  !                Departure Points for Semi-Lagrangian Models                        !
  !                McGegror version with DSS every ugradv                             !
  ! CALLS: 
  ! INPUT: 
  !        
  ! OUTPUT: 
  !-----------------------------------------------------------------------------------!
  subroutine fvm_mcgregordss(elem,fvm,nets,nete, hybrid, deriv, dt_fvm, ordertaylor)
    use derivative_mod, only : derivative_t, ugradv_sphere_routine 
    use edge_mod, only : edgevpack, edgevunpack
    use bndry_mod, only : bndry_exchangev
    
    implicit none
    
    type (element_t), intent(inout)                :: elem(:)
    type (fvm_struct), intent(in)                :: fvm(:)
    
    integer, intent(in)                         :: nets  ! starting thread element number (private)
    integer, intent(in)                         :: nete  ! ending thread element number   (private)
    type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
    
    type (derivative_t), intent(in)                             :: deriv      ! derivative struct
    real (kind=real_kind), intent(in)                           :: dt_fvm
    integer, intent(in)                                         :: ordertaylor
    
    real (kind=real_kind), dimension(nets:nete,np,np,2,nlev)    :: ugradv
    real (kind=real_kind), dimension(nets:nete,np,np,2,nlev)    :: vhat
    integer                                                     :: ie, k, order
    real (kind=real_kind), dimension(np,np,2)                   :: ugradvtmp
    real (kind=real_kind)                                       :: timetaylor
    
    !------------------------------------------------------------------------------------
    timetaylor=1  
    do  order=1,ordertaylor
       timetaylor=-timetaylor*dt_fvm/(order+1)  
       do ie=nets,nete
          if (order==1)then
             ugradv(ie,:,:,:,:)=elem(ie)%derived%vstar(:,:,:,:) 
             vhat(ie,:,:,:,:)=(fvm(ie)%vn0(:,:,:,:) + ugradv(ie,:,:,:,:))/2 
          endif
          do k=1,nlev
             call ugradv_sphere_routine(vhat(ie,:,:,:,k),ugradv(ie,:,:,:,k),deriv,elem(ie),ugradvtmp)
             ugradv(ie,:,:,1,k) = elem(ie)%spheremp(:,:)*ugradvtmp(:,:,1) 
             ugradv(ie,:,:,2,k) = elem(ie)%spheremp(:,:)*ugradvtmp(:,:,2) 
          enddo
          call edgeVpack(edgeveloc,ugradv(ie,:,:,1,:),nlev,0,ie)
          call edgeVpack(edgeveloc,ugradv(ie,:,:,2,:),nlev,nlev,ie)
       enddo
       call t_startf('bndry_exchangeV.edgeveloc', t_detail_medium)
       call bndry_exchangeV(hybrid,edgeveloc)
       call t_stopf('bndry_exchangeV.edgeveloc', t_detail_medium)
       do ie=nets,nete
          call edgeVunpack(edgeveloc,ugradv(ie,:,:,1,:),nlev,0,ie)
          call edgeVunpack(edgeveloc,ugradv(ie,:,:,2,:),nlev,nlev,ie)
          do k=1, nlev  
             ugradv(ie,:,:,1,k)=ugradv(ie,:,:,1,k)*elem(ie)%rspheremp(:,:)
             ugradv(ie,:,:,2,k)=ugradv(ie,:,:,2,k)*elem(ie)%rspheremp(:,:)
             elem(ie)%derived%vstar(:,:,:,k)=elem(ie)%derived%vstar(:,:,:,k) + timetaylor*ugradv(ie,:,:,:,k)
          end do
       end do
    end do
    
  end subroutine fvm_mcgregordss
  !END SUBROUTINE FVM_MCGREGORDSS-------------------------------------------CE-for FVM!
  
  ! ----------------------------------------------------------------------------------!
  !SUBROUTINE FVM_RKDSS-----------------------------------------------CE-for FVM!
  ! AUTHOR: CHRISTOPH ERATH, MARK TAYLOR, 06. December 2012                                             !
  ! DESCRIPTION: ! create a runge kutta taylor serios mixture to calculate the departure grid                            !
  ! CALLS: 
  ! INPUT: 
  !        
  ! OUTPUT: 
  !-----------------------------------------------------------------------------------!
  subroutine fvm_rkdss(elem,fvm,nets,nete, hybrid, deriv, dt_fvm, ordertaylor)
    use derivative_mod, only : derivative_t, ugradv_sphere_routine 
    use edge_mod, only : edgevpack, edgevunpack
    use bndry_mod, only : bndry_exchangev
    
    implicit none
    
    type (element_t), intent(inout)                :: elem(:)
    type (fvm_struct), intent(in)                :: fvm(:)
    
    integer, intent(in)                         :: nets  ! starting thread element number (private)
    integer, intent(in)                         :: nete  ! ending thread element number   (private)
    type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)
    
    type (derivative_t), intent(in)                             :: deriv      ! derivative struct
    real (kind=real_kind), intent(in)                           :: dt_fvm
    integer, intent(in)                                         :: ordertaylor
    
    integer                                                     :: ie, k, order
    real (kind=real_kind), dimension(np,np,2)                   :: ugradvtmp
    real (kind=real_kind)                                       :: timetaylor
    !
    ! RK-SSP 2 stage 2nd order:
    !     x*(t+1) = x(t) + U(x(t),t)                          
    !     x(t+1) = x(t) +  1/2 ( U(x*(t+1),t+1) + U(x(t),t) )       
    ! apply taylor series:
    !  U(x*(t+1),t+1) = U(x(t),t+1) + (x*(t+1)-x(t)) gradU(x(t),t+1)
    !
    !  (x(t+1)-x(t))/dt =  1/2( U(x(t),t+1)+U(x(t),t)) + dt 1/2 U(x(t),t) gradU(x(t),t+1)  
    !
    ! suppose dt = -dt_fvm (we go backward)
    !  (x(t-dt_fvm)-x(t))/-dt_fvm =  1/2( U(x(t),t-dt_fvm)+U(x(t),t)) - dt_fvm 1/2 U(x(t),t) gradU(x(t),t-dt_fvm)  
    !
    !  x(t-dt_fvm) = x(t)) -dt_fvm * [ 1/2( U(x(t),t-dt_fvm)+U(x(t),t)) - dt_fvm 1/2 U(x(t),t) gradU(x(t),t-dt_fvm) ]  
    !
    !    !------------------------------------------------------------------------------------
    do ie=nets,nete
       ! vn0 = U(x,t)
       ! vstar = U(x,t+1)
       do k=1,nlev
          !ugradvtmp(:,:,:)=ugradv_sphere(fvm(ie)%vn0(:,:,:,k),elem(ie)%derived%vstar(:,:,:,k),deriv,elem(ie))
          call ugradv_sphere_routine (elem(ie)%derived%vstar(:,:,:,k),fvm(ie)%vn0(:,:,:,k),deriv,&
                                                                        elem(ie),ugradvtmp(:,:,:))
          
          elem(ie)%derived%vstar(:,:,:,k) = &
               (elem(ie)%derived%vstar(:,:,:,k) + fvm(ie)%vn0(:,:,:,k))/2   - dt_fvm*ugradvtmp(:,:,:)/2
          
          
          
          elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%spheremp(:,:)
          elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%spheremp(:,:)
       enddo
       call edgeVpack(edgeveloc,elem(ie)%derived%vstar,2*nlev,0,ie)
    enddo
    call t_startf('bndry_exchangeV.edgeveloc', t_detail_medium)
    call bndry_exchangeV(hybrid,edgeveloc)
    call t_stopf('bndry_exchangeV.edgeveloc', t_detail_medium)
    do ie=nets,nete
       call edgeVunpack(edgeveloc,elem(ie)%derived%vstar,2*nlev,0,ie)
       do k=1, nlev  
          elem(ie)%derived%vstar(:,:,1,k) = elem(ie)%derived%vstar(:,:,1,k)*elem(ie)%rspheremp(:,:)
          elem(ie)%derived%vstar(:,:,2,k) = elem(ie)%derived%vstar(:,:,2,k)*elem(ie)%rspheremp(:,:)
       end do
    end do
    
  end subroutine fvm_rkdss
  !END SUBROUTINE FVM_rkdss-------------------------------------------CE-for FVM!

end module fvm_line_integrals_mod
