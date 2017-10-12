#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module fvm_consistent_se_cslam
  use kinds, only : real_kind, int_kind
  use edge_mod, only : initghostbufferTR, freeghostbuffertr, &
       ghostVpack, ghostVunpack
  use edgetype_mod, only: ghostbuffertr_t
  use bndry_mod, only: ghost_exchangeV                     
  use dimensions_mod, only: nelem, nelemd, nelemdmax, ne, nc, nhe, nlev, ntrac, np, nhr, nhc, ngpc
  use dimensions_mod, only: irecons_tracer, irecons_air
  use dimensions_mod, only: ie_dbg,pr_dbg
  use time_mod, only : timelevel_t
  use element_mod, only : element_t
  use fvm_control_volume_mod, only: fvm_struct
  use hybrid_mod, only : hybrid_t
  use perf_mod, only : t_startf, t_stopf ! EXTERNAL
  use perf_utils, only: t_detail_low, t_detail_medium, t_detail_high, t_detail_max  ! EXTERNAL

  implicit none
  private
  save
  
  type (ghostBuffertr_t)                      :: cellghostbuf

  real (kind=real_kind), dimension(ngpc), private :: gsweights, gspts
  real (kind=real_kind),parameter       , private :: eps=2.0D-14

  integer, parameter :: ie_limit=730, ilev_limit=30 !dbgxxx
!  integer :: i_dbg,j_dbg,k_dbg,iside_dbg !dbgxxx
!  logical :: ldbg                        !dbgxxx
  public :: run_consistent_se_cslam
contains
  
  
  !
  !**************************************************************************************
  !
  !
  !
  !
  !**************************************************************************************
  !  
  subroutine run_consistent_se_cslam(elem,fvm,hybrid,deriv,dt_fvm,tl,nets,nete,p_top)
    ! ---------------------------------------------------------------------------------
    use fvm_control_volume_mod     , only: n0_fvm, np1_fvm
    use fvm_line_integrals_flux_mod, only: ff_cslam_remap, cslam_runflux !dbg
    use fvm_line_integrals_mod, only: IDEAL_TEST_OFF, fvm_ideal_test!dbg
    use fvm_mod, only: fill_halo_fvm
    ! ---------------------------------------------------------------------------------  
    use fvm_reconstruction_mod, only: reconstruction
    use fvm_line_integrals_mod, only: gauss_points
    ! ---------------------------------------------------------------------------------
    ! -----------------------------------------------
    use derivative_mod, only : derivative_t
    use dimensions_mod, only : lprint !xxx debugging
    
    implicit none
    type (element_t), intent(inout)             :: elem(:)
    type (fvm_struct), intent(inout)            :: fvm(:)
    type (hybrid_t), intent(in)                 :: hybrid   ! distributed parallel structure (shared)   
    type (ghostBuffertr_t)                      :: fluxghostbuf    ! buffer for se_flux halo exchange
    type (TimeLevel_t)                          :: tl              ! time level struct
    type (derivative_t)                         :: deriv           ! derivative struct - only used for dbg
    
    integer, intent(in)                         :: nets  ! starting thread element number (private)
    integer, intent(in)                         :: nete  ! ending thread element number   (private)
    real (kind=real_kind), intent(in)           :: dt_fvm
    real (kind=real_kind), intent(in)           :: p_top
    
    !high-order air density reconstruction
    real (kind=real_kind) :: cair    (1-nhe:nc+nhe,1-nhe:nc+nhe,irecons_air,nlev,nets:nete) 
    real (kind=real_kind) :: ctracer(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons_tracer,ntrac)    
    real (kind=real_kind) :: inv_area(nc,nc)
    integer               :: i,j,k,ie,itr, ntmp
    

    !
    ! fill halo for dp_fvm and c
    !
    do ie=nets,nete       
       do k=1,nlev
          elem(ie)%sub_elem_mass_flux(:,:,:,k) = dt_fvm*elem(ie)%sub_elem_mass_flux(:,:,:,k)*fvm(ie)%dp_ref_inverse(k)
          fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)   =         fvm(ie)%dp_fvm (1:nc,1:nc,k,n0_fvm)*fvm(ie)%dp_ref_inverse(k)
       end do
    end do
!    call set_flow_cases(fvm,elem,nets,nete)!dbgxxx
!    call plot_element_numbers(fvm,elem,nets,nete) !dbgxx
    !
    ! phl-opt: does not scale well with number of tracers (Wallmax 0.095 for 126 tracers and 0.022 for 40 tracers)
    !
    call fill_halo_fvm(elem,fvm,hybrid,nets,nete,n0_fvm)   
    
    if (fvm_ideal_test /= IDEAL_TEST_OFF) then
      !
      !       ! overwrite elem(ie)%sub_elem_mass_flux with ff-cslam fluxes
      !       !
      call cslam_runflux(elem,fvm,hybrid,deriv,dt_fvm,tl,nets,nete,p_top,.TRUE.)
      !!       call unit_test_for_flow_cases(elem,fvm,nets,nete)!dbg
      !       write(*,*) "IDEAL_TEST"
      !!       stop
    end if
    
    call gauss_points(ngpc,gsweights,gspts) !set gauss points/weights
    gspts = 0.5D0*(gspts+1.0D0) !shift location so in [0:1] instead of [-1:1]
    
    call t_startf('fvm prep 2', t_detail_high)    
    call initghostbufferTR(fluxghostbuf,4,nlev,nhe,nc)!
    do ie=nets,nete
      fvm(ie)%se_flux    (1:nc,1:nc,:,:) = elem(ie)%sub_elem_mass_flux(:,:,:,:) 
      !
      ! phl-opt: better first guess displacement
      !
      do k=1,nlev
        if (irecons_air>1) then
          call reconstruction(fvm(ie)%dp_fvm(:,:,k,n0_fvm), fvm(ie),cair(:,:,:,k,ie),irecons_air,.false.)
        else
           !
           ! default branch
           !
          cair(1-nhe:nc+nhe,1-nhe:nc+nhe,1,k,ie) = fvm(ie)%dp_fvm(1-nhe:nc+nhe,1-nhe:nc+nhe,k,n0_fvm)
        end if
      end do
      ie_dbg=ie
      call compute_displacements_for_swept_areas (fvm(ie),cair(:,:,:,:,ie),irecons_air)!new
      call ghostVpack  (fluxghostbuf, fvm(ie)%se_flux(:,:,:,:),nhe,nc,4,nlev,0,elem(ie)%desc)
    end do
    call ghost_exchangeV(hybrid,fluxghostbuf,nhe,nc,nlev)
    do ie=nets,nete
      call ghostVunpack  (fluxghostbuf, fvm(ie)%se_flux(:,:,:,:),nhe,nc,4,nlev,0,elem(ie)%desc)
      call ghost_flux_unpack(fvm(ie))
    enddo
    call freeghostbuffertr(fluxghostbuf)
    call t_stopf('fvm prep 2', t_detail_high)
    !write(*,*) "done computing displacements"
    do ie=nets,nete
       inv_area=1.0D0/fvm(ie)%area_sphere(1:nc,1:nc)  
       do k=1,nlev
          call t_startf('fvm tracers reconstruct', t_detail_high)
          call reconstruction(fvm(ie)%c(1-nhc:nc+nhc,1-nhc:nc+nhc,k,1:ntrac,n0_fvm),&
                  fvm(ie),ctracer(:,:,:,:),irecons_tracer,.true.)
          call t_stopf('fvm tracers reconstruct', t_detail_high)
          call t_startf('fvm swept_flux', t_detail_high)
          !
          ! tracers and air updated in swept_flux (MASTER SUBROUTINE)
          !
          !          k_dbg=k;ie_dbg=ie
          call swept_flux(elem(ie),fvm(ie),k,ie,ctracer,inv_area) !new version of air_flux_remap
          call t_stopf('fvm swept_flux', t_detail_high)          
       end do
       !
       ! surface pressure implied by fvm
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
    ntmp     = np1_fvm
    np1_fvm  = n0_fvm
    n0_fvm   = ntmp
  end subroutine run_consistent_se_cslam
  
  
  subroutine swept_flux(elem,fvm,ilev,ie,ctracer,inv_area)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    use fvm_control_volume_mod     , only: n0_fvm, np1_fvm
    implicit none
    type (element_t) , intent(in)   :: elem
    type (fvm_struct), intent(inout):: fvm
    integer          , intent(in) :: ilev
    integer          , intent(in) :: ie !dbg
    real (kind=real_kind), intent(inout) :: ctracer(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons_tracer,ntrac)
    real (kind=real_kind), intent(in) :: inv_area(nc,nc)

    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides) :: base_vtx
    integer                  , dimension(2,num_area, imin:imax,imin:imax,num_sides) :: idx
    real (kind=real_kind)    , dimension(imin:imax,imin:imax,num_sides)             :: mass_flux_se
    real (kind=real_kind)    , dimension(irecons_tracer,num_area) :: weights


    real (kind=real_kind)                     :: gamma
    integer :: i,j,iside,iarea,ix,iw

    integer, parameter :: num_seg_max=5
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)               :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8) :: x_start, dgam_vec

    REAL(KIND=real_kind) :: flux,flux_tracer(ntrac)

    REAL(KIND=real_kind), dimension(num_area) :: dp_area

    REAL(KIND=real_kind), dimension(1:nc,1:nc) :: invtmp         

    logical :: tl1,tl2,tr1,tr2

    integer, dimension(4), parameter :: imin_side = (/1   ,0   ,1   ,1   /)
    integer, dimension(4), parameter :: imax_side = (/nc  ,nc  ,nc  ,nc+1/)
    integer, dimension(4), parameter :: jmin_side = (/1   ,1   ,0   ,1   /)
    integer, dimension(4), parameter :: jmax_side = (/nc+1,nc  ,nc  ,nc  /)

    integer :: iseg, iseg_tmp,flowcase,ii,jj,itr

    call define_swept_areas(fvm,ilev,displ,base_vec,base_vtx,idx)

    mass_flux_se = 0.0D0
    mass_flux_se(1:nc,1:nc,1:4)  = -elem%sub_elem_mass_flux(1:nc,1:nc,1:4,ilev)
    mass_flux_se(0   ,1:nc,2  )  =  elem%sub_elem_mass_flux(1   ,1:nc,4  ,ilev)
    mass_flux_se(nc+1,1:nc,4  )  =  elem%sub_elem_mass_flux(nc  ,1:nc,2  ,ilev)
    mass_flux_se(1:nc,0   ,3  )  =  elem%sub_elem_mass_flux(1:nc,1   ,1  ,ilev)
    mass_flux_se(1:nc,nc+1,1  )  =  elem%sub_elem_mass_flux(1:nc,nc  ,3  ,ilev)
    !
    ! prepare for air/tracer update
    ! 
    fvm%dp_fvm(1:nc,1:nc,ilev,np1_fvm) = fvm%dp_fvm(1:nc,1:nc,ilev,n0_fvm)*fvm%area_sphere(1:nc,1:nc)
    do itr=1,ntrac
       fvm%c(1:nc,1:nc,ilev,itr,np1_fvm) = fvm%c(1:nc,1:nc,ilev,itr,n0_fvm)*fvm%dp_fvm(1:nc,1:nc,ilev,np1_fvm)
       do iw=1,irecons_tracer
          ctracer(1-nhe:nc+nhe,1-nhe:nc+nhe,iw,itr)=ctracer(1-nhe:nc+nhe,1-nhe:nc+nhe,iw,itr)*&
                                                fvm%dp_fvm(1-nhe:nc+nhe,1-nhe:nc+nhe,ilev,n0_fvm)
       end do
    end do

    do iside=1,4
       do j=jmin_side(iside),jmax_side(iside)
          do i=imin_side(iside),imax_side(iside)
             !             if (fvm%se_flux(i,j,iside,ilev)>eps) then
             if (mass_flux_se(i,j,iside)>eps) then
                !
                !        ||             ||
                !  tl1   ||             || tr1   
                !        ||             ||        
                !  =============================
                !        ||             ||      
                !  tl2   ||             || tr2
                !        ||             ||      
                !      
                tl1 = displ(3,i,j,iside)<0.0D0.and.displ(6,i,j,iside).ge.0.0D0 !departure point in tl1 quadrant
                tl2 = displ(6,i,j,iside)<0.0D0.and.displ(7,i,j,iside)   >0.0D0 !departure point in tl2 quadrant
                tr1 = displ(2,i,j,iside)<0.0D0.and.displ(4,i,j,iside).ge.0.0D0 !departure point in tr1 quadrant
                tr2 = displ(4,i,j,iside)<0.0D0.and.displ(5,i,j,iside)   >0.0D0 !departure point in tr2 quadrant

                !
                ! pathological cases
                !
                !        |  ||           ||                      ||           ||
                !        |  ||-----------||                      ||-----------||
                !        |  ||           ||                      ||           ||
                !  ================================     =================================
                !           ||           ||                   |  ||           ||
                !  ---------||           ||             ------|--||           ||
                !           ||           ||                   |  ||           ||
                !
!                tl1=tl1.or.tl2
!                tr1=tr1.or.tr2
!                tl1=displ(3,i,j,iside)<0.0D0.and..not.(tl1.and.tl2)
!                tr1=displ(2,i,j,iside)<0.0D0.and..not.(tr1.and.tr2)
                
                num_seg=-1; num_seg_static=-1 !initialization
                if (.not.tl1.and..not.tl2.and..not.tr1.and..not.tr2) then
                   flowcase=0
                   !
                   !        ||             ||                 ||             ||                ||             ||         
                   !        ||  *       *  ||                 ||  *----------*                 |*----------*  ||         
                   !        || /         \ ||                 || /           ||                ||           \ ||       
                   !        ||/           \||                 ||/            ||                ||            \||
                   !  =============================     =============================     ============================= 
                   !        ||             ||                 ||             ||                ||             ||
                   !
                   !
                   call define_area3_center (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                        num_seg_static,x_start, dgam_vec,fvm%se_flux(i,j,iside,ilev))

                   gamma=fvm%se_flux(i,j,iside,ilev)
                else
                   if (tl1.and.tr1) then
                      flowcase=1
                      !
                      !
                      !  tl1   ||             || tr1             ||             ||                ||             ||
                      !     *--||-------------||--*           *--||-------------||                ||-------------||--*
                      !      \ ||             || /             \ ||             ||\              /||             || /   
                      !       \||             ||/               \||             || \            / ||             ||/
                      !  =============================     =========================*===     ==*========================== 
                      !        ||             ||                 ||             ||                ||             ||
                      !
                      call define_area2           (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static,&
                           num_seg, num_seg_static,x_start, dgam_vec)
                      call define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static,&
                           num_seg, num_seg_static,x_start, dgam_vec)
                      call define_area4           (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static,&
                           num_seg, num_seg_static,x_start, dgam_vec)
                      gamma=1.0D0
                   else if (tl1.and..not.tr1.and..not.tr2) then
                      flowcase=2
                      !
                      !        ||             ||                 ||             ||                ||             ||
                      !     *--||----------*  ||                /||----------*  ||             *--||-------------*
                      !      \ ||           \ ||               / ||           \ ||              \ ||             ||       
                      !       \||            \||              /  ||            \||               \||             ||
                      !  =============================     ==*==========================     ============================= 
                      !        ||             ||                 ||             ||                ||             ||
                      !
                      call define_area2     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,&
                           x_start, dgam_vec)
                      call define_area3_left(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,&
                           x_start, dgam_vec)
                      gamma=1.0D0
                   else if (tr1.and..not.tl1.and..not.tl2) then !displ(3).ge.0.0D0) then
                      flowcase=3
                      !
                      !        ||  *----------||--*              ||  *----------||\                *-------------||--*
                      !        || /           || /               || /           || \              ||             || /      
                      !        ||/            ||/                ||/            ||  \             ||             ||/
                      !  =============================     ==========================*==     ============================= 
                      !        ||             ||                 ||             ||                ||             ||
                      !        ||             ||                 ||             ||                ||             ||
                      !        ||             ||                 ||             ||                ||             ||
                      !
                      call define_area3_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
                           num_seg_static, x_start, dgam_vec)
                      call define_area4      (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
                           num_seg_static, x_start, dgam_vec)
                      gamma=1.0D0
                   else if (tl2.and..not.tr1.and..not.tr2) then !displ(2).ge.0.0D0) then
                      flowcase=4
                      !
                      !        ||----------*  ||                 ||-------------*        
                      !       /||           \ ||                /||             ||       
                      !      / ||            \||               / ||             ||       
                      !  ===/=========================     ===/========================= 
                      !     | /||             ||              | /||             ||       
                      !     |/ ||             ||              |/ ||             ||       
                      !     *  ||             ||              *  ||             ||       
                      !
                      call define_area1_area2(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      call define_area3_left (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,&
                           x_start, dgam_vec)
                      gamma = 1.0D0
                   else if (tr2.and..not.tl1.and..not.tl2) then !displ(3).ge.0.0D0) then
                      flowcase=5
                      !                case(5)
                      !
                      !
                      !        ||  *-----2----||       
                      !        || /1         3||\      
                      !        ||/      4     || \     
                      !  ============================= 
                      !        ||             ||\ |    
                      !        ||             || \|    
                      !        ||             ||  *    
                      !
                      call define_area3_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      call define_area4_area5(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      gamma=1.0D0
                   else if (tl2.and.tr1.and..not.tr2) then
                      flowcase=6
                      !                case(6)
                      !
                      !
                      !        ||-------------||--*    
                      !       /||             || /     
                      !      / ||             ||/      
                      !  ===/========================= 
                      !     | /||             ||       
                      !     |/ ||             ||       
                      !     *  ||             ||       
                      !
                      !
                      call define_area1_area2     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      call define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      call define_area4           (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)

                      gamma=1.0D0
                   else if (tr2.and.tl1.and..not.tl2) then
                      flowcase=7
                      !                case(7)
                      !
                      !
                      !     *--||-------------||       
                      !      \ ||             ||\      
                      !       \||             || \     
                      !  ============================= 
                      !        ||             ||\ |    
                      !        ||             || \|    
                      !        ||             ||  *    
                      !
                      !
                      call define_area2           (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      call define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      call define_area4_area5     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      gamma =  1.0D0
                   else if (tl2.and.tr2) then
                      flowcase=8
                      !                case(8)
                      !
                      !
                      !        ||-------------||       
                      !       /||             ||\      
                      !      / ||             || \     
                      !  ============================= 
                      !     | /||             ||\ |    
                      !     |/ ||             || \|    
                      !     *  ||             ||  *    
                      !
                      !
                      !
                      !
                      !
                      call define_area1_area2     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      call define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      call define_area4_area5     (i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg,&
                           num_seg_static,x_start, dgam_vec)
                      gamma =  1.0D0
                   else
                      !                case default
                      write(*,*) "ERROR - unknown flow case",flowcase !dbgxxx
                      write(*,*) "ie",ie
                      stop
                   end if
                end if
                !
                ! iterate to get flux area
                !
                call t_startf('fvm swept_area: get_gamma', t_detail_high)
                do iarea=1,num_area
                   dp_area(iarea) = fvm%dp_fvm(idx(1,iarea,i,j,iside),idx(2,iarea,i,j,iside),ilev,n0_fvm)
                end do
                !                i_dbg=i;j_dbg=j;iside_dbg=iside!dbg_xxx
                call get_flux_segments_area_iterate(x,x_static,dx_static,dx,x_start,dgam_vec,num_seg,num_seg_static,&
                     num_seg_max,num_area,dp_area,flowcase,gamma,mass_flux_se(i,j,iside),-1.0D0,2.0D0)
                call t_stopf('fvm swept_area: get_gamma', t_detail_high)
                !
                ! pack segments for high-order weights computation
                !
                do iarea=1,num_area
                   do iseg=1,num_seg_static(iarea)
                      iseg_tmp=num_seg(iarea)+iseg
                      x (:,iseg_tmp,iarea)  = x_static (:,iseg,iarea)
                      dx(:,iseg_tmp,iarea)  = dx_static(:,iseg,iarea)
                   end do
                   num_seg(iarea)=num_seg(iarea)+MAX(0,num_seg_static(iarea))
                end do
!                call write_cells(fvm,x(:,1:num_seg_max,1:num_area),dx(:,1:num_seg_max,1:num_area),& !dbg
!                     num_area,num_seg_max,num_seg,ie,i,j,iside,ilev,flowcase)                       !dbg
                !
                ! compute higher-order weights
                !
                call t_startf('fvm swept_area: get_high_order_w', t_detail_high)
                call get_high_order_weights_over_areas(x,dx,num_seg,num_seg_max,num_area,weights)
                call t_stopf('fvm swept_area: get_high_order_w', t_detail_high)
                !
                !**************************************************
                !
                ! remap air and tracers
                !
                !**************************************************
                ! 
                call t_startf('fvm swept_area: remap', t_detail_high)
                flux=0.0D0; flux_tracer=0.0D0
                do iarea=1,num_area
                   if (num_seg(iarea)>0) then
                      ii=idx(1,iarea,i,j,iside); jj=idx(2,iarea,i,j,iside)
                      flux=flux+weights(1,iarea)*fvm%dp_fvm(ii,jj,ilev,n0_fvm)
                      do itr=1,ntrac
                         do iw=1,irecons_tracer
                            flux_tracer(itr) = flux_tracer(itr)+weights(iw,iarea)*ctracer(ii,jj,iw,itr)
                         end do
                      end do
                   end if
                end do
                fvm%dp_fvm(i  ,j  ,ilev        ,np1_fvm) = fvm%dp_fvm(i  ,j  ,ilev        ,np1_fvm)-flux
                fvm%     c(i  ,j  ,ilev,1:ntrac,np1_fvm) = fvm%     c(i  ,j  ,ilev,1:ntrac,np1_fvm)-flux_tracer(1:ntrac)
                !
                ! update flux in nearest neighbor cells
                !
                if (iside==1) then
                   fvm%dp_fvm(i,j-1,ilev        ,np1_fvm) = fvm%dp_fvm(i,j-1,ilev        ,np1_fvm)+flux
                   fvm%     c(i,j-1,ilev,1:ntrac,np1_fvm) = fvm%     c(i,j-1,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
                end if
                if (iside==2) then
                   fvm%dp_fvm(i+1,j,ilev        ,np1_fvm) = fvm%dp_fvm(i+1,j,ilev        ,np1_fvm)+flux
                   fvm%     c(i+1,j,ilev,1:ntrac,np1_fvm) = fvm%     c(i+1,j,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
                end if
                if (iside==3) then
                   fvm%dp_fvm(i,j+1,ilev        ,np1_fvm) = fvm%dp_fvm(i,j+1,ilev        ,np1_fvm)+flux
                   fvm%     c(i,j+1,ilev,1:ntrac,np1_fvm) = fvm%     c(i,j+1,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
                end if
                if (iside==4) then
                   fvm%dp_fvm(i-1,j,ilev        ,np1_fvm) = fvm%dp_fvm(i-1,j,ilev        ,np1_fvm)+flux
                   fvm%     c(i-1,j,ilev,1:ntrac,np1_fvm) = fvm%     c(i-1,j,ilev,1:ntrac,np1_fvm)+flux_tracer(1:ntrac)
                end if
                call t_stopf('fvm swept_area: remap', t_detail_high)
             end if
          end do
       end do
    end do
    !
    ! convert to mixing ratio
    !
    invtmp=1.0D0/fvm%dp_fvm(1:nc,1:nc,ilev,np1_fvm)         
    do itr=1,ntrac
       fvm%c(1:nc,1:nc,ilev,itr,np1_fvm) = fvm%c(1:nc,1:nc,ilev,itr,np1_fvm)*invtmp
    end do    
    !
    ! scale back pressure
    !
    fvm%dp_fvm(1:nc,1:nc,ilev,np1_fvm) = fvm%dp_fvm(1:nc,1:nc,ilev,np1_fvm)*fvm%dp_ref(ilev)*inv_area(:,:)
  end subroutine swept_flux

  subroutine ghost_flux_unpack(fvm)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    implicit none
    type (fvm_struct), intent(inout) :: fvm
    
    integer :: i,j,k,ishft
    real (kind=real_kind) :: tmp(nc,nlev), tmp2(nlev)
    
    !
    ! rotate coordinates if needed
    !       
    if (fvm%cubeboundary.NE.0) then
       do k=1,nlev
          do j=1-nhe,nc+nhe
             do i=1-nhe,nc+nhe
                ishft = NINT(fvm%flux_orient(2,i,j))
                fvm%se_flux(i,j,1:4,k) = cshift(fvm%se_flux(i,j,1:4,k),shift=ishft)
             end do
          end do
       end do
       !
       ! non-existent cells in physical space - necessary?
       !
       if (fvm%cubeboundary==nwest) then
          fvm%se_flux(1-nhe:0,nc+1 :nc+nhe,:,:) = 0.0D0
       else if (fvm%cubeboundary==swest) then
          fvm%se_flux(1-nhe:0,1-nhe:0     ,:,:) = 0.0D0
       else if (fvm%cubeboundary==neast) then
          fvm%se_flux(nc+1 :nc+nhe,nc+1 :nc+nhe,:,:) = 0.0D0
       else if (fvm%cubeboundary==seast) then
          fvm%se_flux(nc+1 :nc+nhe,1-nhe:0,:,:) = 0.0D0
       end if
    end if
  end subroutine ghost_flux_unpack

  subroutine compute_displacements_for_swept_areas(fvm,cair,irecons)
    use fvm_control_volume_mod     , only: n0_fvm
    implicit none
    type (fvm_struct), intent(inout)     :: fvm
    integer, intent(in) :: irecons
    real (kind=real_kind)                :: cair(1-nhe:nc+nhe,1-nhe:nc+nhe,irecons,nlev) !high-order air density reconstruction
    !
    !   flux iside 1                     flux iside 3                    flux iside 2       flux iside 4
    !
    !   |          |                     |  ---1--> |                    |    --2-->|       |--1-->    |
    !  -4----------3-   /\              -4----------3-                  -4----------3-     -4----------3-   ||
    !   |          |   /||\              |\\\\\\\\\\|    ||              |   |\\\\\\|       |\\\\\\|   |
    !   |  --2-->  |    || dv(1)         |\\\\\\\\\\|    ||              |   |\\\\\\|       |\\\\\\|   |
    !   |----------|    ||               |----------|    || dv(3)        |   |\\\\\\|       |\\\\\\|   |
    !   |\\\\\\\\\\|    ||               | <--2---  |   \||/             |   |\\\\\\|       |\\\\\\|   |
    !   |\\\\\\\\\\|    ||               |          |    \/              |   |\\\\\\|       |\\\\\\|   |
    !  -1----------2-                   -1----------2-                  -1----------2-     -1----------2-
    !   |  <--1--  |                     |          |                    |    <--1--|       |<--2--
    !
    !                                                                     /                          \
    !   line-integral                                                    <==========         =========>
    !   from vertex 2                                                     \  dv(2)              dv(4)/
    !   to 1
    ! 
    !   Note vertical
    !   lines have 
    !   zero line-
    !   integral!
    !
    integer               :: i,j,k,iside,ix
    integer, parameter :: num_area=1, num_seg_max=2
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area,4,nc,nc) :: x_static, dx_static
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area,4,nc,nc) :: x, dx
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area)         :: x_tmp, dx_tmp
    integer             , dimension(              num_area,4      ) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8,                   4,nc,nc) :: x_start, dgam_vec
    REAL(KIND=real_kind), dimension(num_area) :: dp_area
    integer, dimension(4) :: flowcase
    REAL(KIND=real_kind)  :: gamma, flux_se
    
    num_seg_static(1,1) =  1; num_seg(1,1) = 1; flowcase(1) = -1
    num_seg_static(1,2) =  0; num_seg(1,2) = 2; flowcase(2) = -2
    num_seg_static(1,3) =  1; num_seg(1,3) = 1; flowcase(3) = -1
    num_seg_static(1,4) =  0; num_seg(1,4) = 2; flowcase(4) = -4

    do j=1,nc
       do i=1,nc
          do ix=1,2
             iside=1;
             x_static (ix,1,1,iside,i,j) = fvm%flux_vertex_cart(i,j,ix,2)
             dx_static(ix,1,1,iside,i,j) = fvm%flux_vertex_cart(i,j,ix,1)-fvm%flux_vertex_cart(i,j,ix,2)
             x_start  (ix,1,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,1)
             x_start  (ix,2,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,2)
             dgam_vec (ix,1,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,4)-fvm%flux_vertex_cart(i,j,ix,1)
             !
             ! compute first guess
             !
             gamma                       = 0.5D0
             x        (ix,1,1,iside,i,j) = x_start(ix,1,iside,i,j)+gamma*dgam_vec(ix,1,iside,i,j)
             dx       (ix,1,1,iside,i,j) = -dx_static(ix,1,1,iside,i,j)
             !
             ! side 2
             !
             iside=2;
             x_start  (ix,1,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,2)
             x_start  (ix,2,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,3)
             dgam_vec (ix,1,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,1)-fvm%flux_vertex_cart(i,j,ix,2)
             x        (ix,1,1,iside,i,j) = x_start(ix,1,iside,i,j)
             !
             ! compute first guess - gamma=1
             !
             gamma                       = 0.5D0
             dx       (ix,1,1,iside,i,j) =  gamma*dgam_vec (ix,1,  iside,i,j)
             x        (ix,2,1,iside,i,j) =  x_start(ix,2,iside,i,j)+gamma*dgam_vec(ix,1,iside,i,j)
             dx       (ix,2,1,iside,i,j) = -gamma*dgam_vec (ix,1,  iside,i,j)
             !
             ! side 3
             !
             iside=3;
             x_static (ix,1,1,iside,i,j) = fvm%flux_vertex_cart(i,j,ix,4)
             dx_static(ix,1,1,iside,i,j) = fvm%flux_vertex_cart(i,j,ix,3)-fvm%flux_vertex_cart(i,j,ix,4)
             x_start  (ix,1,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,3)
!             x_start  (ix,2,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,4)
             dgam_vec (ix,1,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,2)-fvm%flux_vertex_cart(i,j,ix,3)
             !
             ! compute first guess - gamma=1
             !
             gamma                       = 0.5D0
             x        (ix,1,1,iside,i,j) = x_start(ix,1,iside,i,j)+gamma*dgam_vec(ix,1,iside,i,j)
             dx       (ix,1,1,iside,i,j) = -dx_static(ix,1,1,iside,i,j)
             !
             ! side 4
             !
             iside=4;
             x_start  (ix,1,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,1)
             x_start  (ix,2,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,4)
             dgam_vec (ix,1,  iside,i,j) = fvm%flux_vertex_cart(i,j,ix,2)-fvm%flux_vertex_cart(i,j,ix,1)
             x        (ix,2,1,iside,i,j) = x_start(ix,2,iside,i,j)
             !
             ! compute first guess - gamma=1
             !
             gamma                       = 0.5D0
             dx       (ix,2,1,iside,i,j) =  gamma*dgam_vec (ix,1,  iside,i,j)
             x        (ix,1,1,iside,i,j) =  x_start(ix,1,iside,i,j)+gamma*dgam_vec(ix,1,iside,i,j)
             dx       (ix,1,1,iside,i,j) = -gamma*dgam_vec (ix,1,  iside,i,j)
          end do
       end do
    end do

    do k=1,nlev
       do j=1,nc
          do i=1,nc
             dp_area = cair(i,j,1,k)
             do iside=1,4
                flux_se = -fvm%se_flux(i,j,iside,k)
                if (flux_se>eps) then
                   gamma=0.5D0
!                   iside_dbg=iside!dbgxxx
!                   i_dbg=i;j_dbg=j;k_dbg=k;
                   !
                   ! this copying is necessary since get_flux_segments_area_iterate change x and dx
                   !
                   x_tmp (:,1:num_seg(1,iside),:)=x (:,1:num_seg(1,iside),:,iside,i,j)
                   dx_tmp(:,1:num_seg(1,iside),:)=dx(:,1:num_seg(1,iside),:,iside,i,j)

                   call get_flux_segments_area_iterate(&
                        x_tmp(:,:,:),x_static(:,:,:,iside,i,j),dx_static(:,:,:,iside,i,j),dx_tmp(:,:,:),&
                        x_start(:,:,iside,i,j),dgam_vec(:,:,iside,i,j),num_seg(:,iside),num_seg_static(:,iside),&
                        num_seg_max,num_area,dp_area,flowcase(iside),gamma,flux_se,0.0D0,1.0D0)
                   fvm%se_flux(i,j,iside,k) = ABS(SUM(gamma*dgam_vec(:,1,iside,i,j)))
                else
                   fvm%se_flux(i,j,iside,k) = 0.0D0
                end if
             enddo
          end do
       end do
    end do
  end subroutine compute_displacements_for_swept_areas



  subroutine get_flux_segments_area_iterate(x,x_static,dx_static,dx,x_start,dgam_vec,num_seg,num_seg_static,&
       num_seg_max,num_area,c,flow_case,gamma,flux,gamma_min,gamma_max)
    implicit none
    integer                                                , intent(in)    :: num_area, num_seg_max
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(in)    :: x_static, dx_static
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx
    integer             , dimension(num_area              ), intent(in) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8)                   , intent(in) :: x_start, dgam_vec
    REAL(KIND=real_kind)                                   , intent(inout) :: gamma
    REAL(KIND=real_kind)                                   , intent(in) :: flux,gamma_min,gamma_max
    integer                                                , intent(in) :: flow_case

    real (kind=real_kind), dimension(num_area)             , intent(in) :: c

    real (kind=real_kind)                                :: flux_static
    real (kind=real_kind)                                :: weight_area, xtmp(2), xtmp2(2)
    real (kind=real_kind)                                :: gamma1, gamma2, gamma3, dgamma, f1, f2
    real (kind=real_kind), dimension(  ngpc  ) :: xq,yq    
    real (kind=real_kind), dimension(  ngpc,1) :: F !linear

    integer :: iseg,iarea,iter,ipt
    integer, parameter :: iter_max=10
    !
    ! compute static line-integrals (not necessary to recompute them for every iteration)
    !
    flux_static = 0.0D0  
    do iarea=1,num_area
       weight_area=0.0D0
       do iseg=1,num_seg_static(iarea)

!rck vector directive needed here
!DIR$ SIMD
          do ipt=1,ngpc                
             xq(ipt) = x_static(1,iseg,iarea)+dx_static(1,iseg,iarea)*gspts(ipt)! create quadrature point locations
             yq(ipt) = x_static(2,iseg,iarea)+dx_static(2,iseg,iarea)*gspts(ipt)                
             F(ipt,1) = yq(ipt)/(SQRT(1.0D0+xq(ipt)*xq(ipt) + yq(ipt)*yq(ipt))*(1.0D0+xq(ipt)*xq(ipt)))! potential ! potential             
          enddo
          weight_area = weight_area+sum(gsweights(:)*F(:,1))*dx_static(1,iseg,iarea) !integral
       end do
       flux_static = flux_static+weight_area*0.5D0*c(iarea)      !add to swept flux
    end do
    !
    ! initilization
    !
    gamma1=0.0D0; f1=-flux   ! zero flux guess 1
    !
    ! compute flux integrals of first guess passed to subroutine
    !
    gamma2=gamma
    f2 = 0.0D0  
    do iarea=1,num_area
       weight_area=0.0D0
       do iseg=1,num_seg(iarea)
!rck vector directive needed here
!DIR$ SIMD
          do ipt=1,ngpc                
             xq(ipt) = x(1,iseg,iarea)+dx(1,iseg,iarea)*gspts(ipt)! create quadrature point locations
             yq(ipt) = x(2,iseg,iarea)+dx(2,iseg,iarea)*gspts(ipt)                
             F(ipt,1) = yq(ipt)/(SQRT(1.0D0+xq(ipt)*xq(ipt) + yq(ipt)*yq(ipt))*(1.0D0+xq(ipt)*xq(ipt)))! potential ! potential             
          enddo
          weight_area = weight_area+sum(gsweights(:)*F(:,1))*dx(1,iseg,iarea)! integral
       end do
       f2 = f2+weight_area*0.5D0*c(iarea)
    end do
    f2 = f2+flux_static-flux !integral error
    dgamma=(gamma2-gamma1)*f2/(f2-f1);
    gamma3 = gamma2-dgamma;                    ! Newton "guess" for gamma
    gamma1 = gamma2; f1 = f2; gamma2 = gamma3; ! prepare for iteration
    do iter=1,iter_max
       !
       ! update vertex location: flow_case dependent to avoid many zero operations
       !
       select case(flow_case)
       case(-4)
          iarea=1
          dx       (:,2,1) =  gamma3*dgam_vec (:,1)
          x        (:,1,1) =  x_start(:,1)+gamma3*dgam_vec(:,1)
          dx       (:,1,1) = -gamma3*dgam_vec (:,1)

       case(-2)
          iarea=1
          dx       (:,1,iarea) =  gamma3*dgam_vec (:,1)
          x        (:,2,iarea) =  x_start(:,2)+gamma3*dgam_vec(:,1)
          dx       (:,2,iarea) = -gamma3*dgam_vec (:,1)
       case(-1)
          !
          ! to compute first-guess perpendicular displacements for iside=1
          !
          iarea=1
          x        (:,1,iarea) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx       (:,1,iarea) = -dx_static(:,1,iarea)
          x        (:,2,iarea) = x_start(:,2)+gamma3*dgam_vec(:,1)
          dx       (:,2,iarea) = x_start(:,2)-x(:,2,iarea)
       case(0)
          iarea=3
          xtmp = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx       (:,1,iarea) = xtmp(:  )-x(:,1,iarea)           !dynamic - line 2
          x        (:,2,iarea) = xtmp(:  )                        !dynamic - line 3
          dx       (:,2,iarea) = x_static(:,2,iarea)-x(:,2,iarea) !dynamic - line 3
       case(1)
          iarea=2
          xtmp(:        ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)        !dynamic - line 2
          x   (:,2,iarea) = xtmp(:)                     !dynamic  - line 3
          dx  (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 3

          iarea            = 3
          xtmp (:  )       = x_start(:,4)+gamma3*dgam_vec(:,4)
          xtmp2(:  )       = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic 
          x    (:,2,iarea) = xtmp (:)         !dynamic 
          dx   (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic 
          x    (:,3,iarea) = xtmp2(:)              !dynamic
          dx   (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic

          iarea         = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic - line 2
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(2)
          iarea=2
          xtmp(:        ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)        !dynamic - line 2
          x   (:,2,iarea) = xtmp(:)                     !dynamic  - line 3
          dx  (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 3
          
          iarea=3
          xtmp(:        ) = x_start(:,4)+gamma3*dgam_vec(:,4)!   
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)        !dynamic - line 1
          x   (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx  (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(3)
          iarea         = 3
          xtmp    (:  ) = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic - line 2                   
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,2,iarea)-xtmp(:) !dynamic - line 2

          iarea         = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic - line 2
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(4)
          iarea           = 1
          xtmp(:        ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x (:,2,iarea) = xtmp(:)                      !dynamic
          dx(:,2,iarea) = x_static(:,1,iarea)-xtmp(:)  !dynamic

          iarea         = 2
          xtmp    (:  ) = x_start(:,2)+gamma3*dgam_vec(:,2)
          xtmp2   (:  ) = x_start(:,3)+gamma3*dgam_vec(:,3)

          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic

          x (:,2,iarea) = xtmp (:)          !dynamic
          dx(:,2,iarea) = xtmp2(:)-xtmp(:)  !dynamic

          x (:,3,iarea) = xtmp2(:)                !dynamic
          dx(:,3,iarea) = x(:,1,iarea)-xtmp2(:)   !dynamic

          iarea            = 3
          xtmp (:        ) = x_start(:,4)+gamma3*dgam_vec(:,4)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic - line 1                   
          x    (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx   (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(5)
          iarea                = 3
          xtmp    (:  )        = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic - line 2
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,2,iarea)-xtmp(:) !dynamic - line 2

          iarea         = 4 
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          xtmp2   (:  ) = x_start(:,7)+gamma3*dgam_vec(:,7)

          dx(:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
          x (:,2,iarea) = xtmp(:)          !dynamic -line 2
          dx       (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic - line 2
          x        (:,3,iarea) = xtmp2(:)               !dynamic  -line 1
          dx       (:,3,iarea) = x(:,1,iarea)-xtmp2(:)  !dynamic - line 1

          iarea             = 5
          xtmp  (:  )       = x_start(:,8)+gamma3*dgam_vec(:,8)

          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
          x        (:,2,iarea) = xtmp(:)                     !dynamic -line 2
          dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(6)
          iarea = 1
          xtmp(:  ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x (:,2,iarea) = xtmp(:)                      !dynamic
          dx(:,2,iarea) = x_static(:,1,iarea)-xtmp(:)  !dynamic

          iarea         = 2
          xtmp    (:  ) = x_start(:,2)+gamma3*dgam_vec(:,2)                   
          xtmp2   (:  ) = x_start(:,3)+gamma3*dgam_vec(:,3)
                   
          dx(:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic
          x (:,2,iarea) = xtmp (:)          !dynamic
          dx(:,2,iarea) = xtmp2(:)-xtmp(:)  !dynamic
          x (:,3,iarea) = xtmp2(:)                !dynamic
          dx(:,3,iarea) = x(:,1,iarea)-xtmp2(:)   !dynamic

          iarea            = 3
          xtmp (:  )       = x_start(:,4)+gamma3*dgam_vec(:,4)
          xtmp2(:  )       = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic 
          x    (:,2,iarea) = xtmp (:)         !dynamic 
          dx   (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic 
          x    (:,3,iarea) = xtmp2(:)              !dynamic
          dx   (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic

          iarea         = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic - line 2
          x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
          dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(7)
          iarea=2
          xtmp(:        ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)        !dynamic - line 2
          x   (:,2,iarea) = xtmp(:)                     !dynamic  - line 3
          dx  (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 3

          iarea            = 3
          xtmp (:  )       = x_start(:,4)+gamma3*dgam_vec(:,4)
          xtmp2(:  )       = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic 
          x    (:,2,iarea) = xtmp (:)         !dynamic 
          dx   (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic 
          x    (:,3,iarea) = xtmp2(:)              !dynamic
          dx   (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic

          iarea      = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          xtmp2   (:  ) = x_start(:,7)+gamma3*dgam_vec(:,7)

          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic
          x        (:,2,iarea) = xtmp(:)              !dynamic
          dx       (:,2,iarea) = xtmp2(:)-xtmp(:)     !dynamic
          x        (:,3,iarea) = xtmp2(:)               !dynamic
          dx       (:,3,iarea) = x(:,1,iarea)-xtmp2(:)  !dynamic

          iarea      = 5
          xtmp (:  ) = x_start(:,8)+gamma3*dgam_vec(:,8)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
          x    (:,2,iarea) = xtmp(:)                     !dynamic -line 2
          dx   (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case(8)
          iarea = 1
          xtmp(:  ) = x_start(:,1)+gamma3*dgam_vec(:,1)
          dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
          x (:,2,iarea) = xtmp(:)                      !dynamic
          dx(:,2,iarea) = x_static(:,1,iarea)-xtmp(:)  !dynamic

          iarea         = 2
          xtmp    (:  ) = x_start(:,2)+gamma3*dgam_vec(:,2)                   
          xtmp2   (:  ) = x_start(:,3)+gamma3*dgam_vec(:,3)
                   
          dx(:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic
          x (:,2,iarea) = xtmp (:)          !dynamic
          dx(:,2,iarea) = xtmp2(:)-xtmp(:)  !dynamic
          x (:,3,iarea) = xtmp2(:)                !dynamic
          dx(:,3,iarea) = x(:,1,iarea)-xtmp2(:)   !dynamic

          iarea            = 3
          xtmp (:  )       = x_start(:,4)+gamma3*dgam_vec(:,4)
          xtmp2(:  )       = x_start(:,5)+gamma3*dgam_vec(:,5)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic 
          x    (:,2,iarea) = xtmp (:)         !dynamic 
          dx   (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic 
          x    (:,3,iarea) = xtmp2(:)              !dynamic
          dx   (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic

          iarea      = 4
          xtmp    (:  ) = x_start(:,6)+gamma3*dgam_vec(:,6)
          xtmp2   (:  ) = x_start(:,7)+gamma3*dgam_vec(:,7)

          dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic
          x        (:,2,iarea) = xtmp(:)              !dynamic
          dx       (:,2,iarea) = xtmp2(:)-xtmp(:)     !dynamic
          x        (:,3,iarea) = xtmp2(:)               !dynamic
          dx       (:,3,iarea) = x(:,1,iarea)-xtmp2(:)  !dynamic

          iarea      = 5
          xtmp (:  ) = x_start(:,8)+gamma3*dgam_vec(:,8)
          dx   (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
          x    (:,2,iarea) = xtmp(:)                     !dynamic -line 2
          dx   (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
       case default
          write(*,*) "flow case not defined in get_flux_segments_area_iterate",flow_case
!          stop
       end select
       !
       ! compute flux integral
       !
       f2 = 0.0D0  
       do iarea=1,num_area
          weight_area=0.0D0
          do iseg=1,num_seg(iarea)
             xq(:) = x(1,iseg,iarea)+dx(1,iseg,iarea)*gspts(:)! create quadrature point locations
             yq(:) = x(2,iseg,iarea)+dx(2,iseg,iarea)*gspts(:)
             F(:,1) = yq/(SQRT(1.0D0+xq*xq + yq*yq)*(1.0D0+xq*xq))! potential
             weight_area = weight_area+sum(gsweights(:)*F(:,1))*dx(1,iseg,iarea)! integral
          end do
          f2 = f2+weight_area*0.5D0*c(iarea)
       end do
       f2 = f2+flux_static-flux !integral error
       if (ABS(f2)<eps) then
          gamma=gamma3
          !
          ! converged so exit
          !
          exit
       else
          !
          ! Newton increment
          !          
          dgamma=(gamma2-gamma1)*f2/(f2-f1);
          if (ABS(dgamma)>eps) then
             gamma3 = gamma2-dgamma;
          else
             !
             ! dgamma set to minimum displacement to avoid f2-f1=0
             !
             gamma3=gamma2-SIGN(1.0D0,dgamma)*eps
          end if
          !
          ! prepare for next iteration
          !
          gamma1 = gamma2; f1 = f2; gamma2 = gamma3;
       endif
       
       
    end do
    !
    ! DEBUGGING
    !
!    if (gamma<-1.0D-10.or.gamma>2.2D0) then
!       write(*,*) "WARNING: gamma perhaps out of bounds",gamma
!       write(*,*) "flux",flux
!       write(*,*) "i_dbg,j_dbg,k_dbg,iside_dbg,ie_dbg",i_dbg,j_dbg,k_dbg,iside_dbg,ie_dbg
!       write(*,*) "c",c
!    end if
!    if (iter>9) then
!       write(*,*) "WARNING: iter number high",iter,gamma3,iter,f2,flux,flow_case,f2-f1
!       write(*,*) "flux",flux
!       write(*,*) "i_dbg,j_dbg,k_dbg,iside_dbg,ie_dbg",i_dbg,j_dbg,k_dbg,iside_dbg,ie_dbg
!       write(*,*) "c",c
!    end if
!   if (ABS(f2)>eps) then
!       write(*,*) "WARNING: iteration not converged:",gamma3,iter,flux,f1,f2,f1-f2
!       write(*,*) "flux",flux
!       write(*,*) "i_dbg,j_dbg,k_dbg,iside_dbg,ie_dbg",i_dbg,j_dbg,k_dbg,iside_dbg,ie_dbg
!       write(*,*) "c",c
!    end if
!    if (gamma>3.0D0) then
!       write(*,*) "warning gamma large",gamma,iter,flux,f1,f2
!       write(*,*) "i_dbg,j_dbg,k_dbg,iside_dbg,ie_dbg",i_dbg,j_dbg,k_dbg,iside_dbg,ie_dbg
!       write(*,*) "c",c
!    end if
  end subroutine get_flux_segments_area_iterate


  subroutine get_high_order_weights_over_areas(x,dx,num_seg,num_seg_max,num_area,weights)
    implicit none
    integer                                                 , intent(in)    :: num_area, num_seg_max
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area ), intent(inout) :: x, dx
    integer             , dimension(num_area               ), intent(in)    :: num_seg
    REAL(KIND=real_kind), dimension(irecons_tracer,num_area), intent(out)   :: weights

    real (kind=real_kind), dimension(ngpc,num_seg_max               ) :: xq,yq        !quadrature points along line segments
    real (kind=real_kind), dimension(ngpc,num_seg_max,irecons_tracer) :: F            !potentials
    real (kind=real_kind), dimension(                 irecons_tracer) :: weights_area
    real (kind=real_kind), dimension(ngpc,num_seg_max) :: xq2, yrh, rho, tmp !intermediate variables for optimization
    REAL(KIND=real_kind) , dimension(ngpc,num_seg_max) :: xq2ir, xq2i, rhoi  !intermediate variables for optimization

    integer :: iseg,iarea,i,j,k

    weights(1:irecons_tracer,1:num_area) = 0.0D0 !may not be necessary dbgxxx
    do iarea=1,num_area
       do iseg=1,num_seg(iarea)
          xq(:,iseg) = x(1,iseg,iarea)+dx(1,iseg,iarea)*gspts(:)
          yq(:,iseg) = x(2,iseg,iarea)+dx(2,iseg,iarea)*gspts(:)
       end do
       !
       ! potentials (equation's 23-28 in CSLAM paper; Lauritzen et al., 2010):
       !
       ! (Rory Kelly optimization)
       !
       do j=1,num_seg(iarea)
!DIR$ SIMD
          do i=1,ngpc
             xq2(i,j)   =  xq(i,j)*xq(i,j)
             xq2i(i,j)  =  1.0D0/(1.0D0+xq2(i,j))   
             xq2ir(i,j) =  SQRT(xq2i(i,j))
             rho(i,j)   =  SQRT(1.0D0+xq2(i,j)+yq(i,j)*yq(i,j))
             rhoi(i,j)  =  1.0D0/rho(i,j)
             yrh(i,j)   =  yq(i,j)*rhoi(i,j)
             tmp(i,j)   =  yq(i,j)*xq2ir(i,j)
             F(i,j,1)   =  yrh(i,j)*xq2i(i,j)                 !F_00 !F_00
             F(i,j,2)   =  xq(i,j)*yrh(i,j)*xq2i(i,j)         !F_10 !F_10
             F(i,j,3)   = -1.0D0*rhoi(i,j)                    !F_01 !F_01
             F(i,j,4)   =  xq2(i,j)*yrh(i,j)*xq2i(i,j)        !F_20 !F_20
             F(i,j,6)   = -xq(i,j)*rhoi(i,j)                  !F_11 !F_11
          enddo
          !
          ! take F(i,j,5) out of loop above since it prevents vectorization
          !
          do i=1,ngpc
             F(i,j,5)   = -yq(i,j)*rhoi(i,j)+log(tmp(i,j)+rho(i,j)*xq2ir(i,j))  !F_02 !F_02
          end do
       enddo
       weights_area = 0.0D0
       do k=1,irecons_tracer
          do iseg=1,num_seg(iarea)             
             weights_area(k) = weights_area(k) + sum(gsweights(:)*F(:,iseg,k))*0.5D0*dx(1,iseg,iarea)
          end do
       end do
       weights(1:irecons_tracer,iarea) = weights_area(1:irecons_tracer)
    end do
  end subroutine get_high_order_weights_over_areas

  subroutine define_swept_areas(fvm,ilev,displ,base_vec,base_vtx,idx)
    use control_mod, only : north, south, east, west, neast, nwest, seast, swest
    implicit none
    type (fvm_struct), intent(inout) :: fvm
    integer          , intent(in)    :: ilev


    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(out) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(out) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(out) :: base_vtx
    integer                  , dimension(2,num_area, imin:imax,imin:imax,num_sides), intent(out) :: idx

    real (kind=real_kind) :: flux_sum     (0:nc+1,0:nc+1,2)
    integer               :: degenerate   (1:nc+1,1:nc+1  )
    integer               :: circular_flow(1:nc+1,1:nc+1  )
    integer               :: illcond      (1:nc+1,1:nc+1)  
    integer               :: ib,i,j,sgn, iside, iarea

    !
    ! set where reconstruction function is as a function of area and side
    !
    integer, dimension(num_area*4), parameter :: idx_shift_tmp = (/-1,-1, 0, 1, 1,&  !iside=1
                                                                    1, 0, 0, 0, 1,&  !iside=2
                                                                    1, 1, 0,-1,-1,&  !iside=3
                                                                   -1, 0, 0, 0,-1/)  !iside=4

    integer, dimension(num_area*4), parameter :: idy_shift_tmp = (/-1, 0, 0, 0,-1,&  !iside=1
                                                                   -1,-1, 0, 1, 1,&  !iside=2
                                                                    1, 0, 0, 0, 1,&  !iside=3
                                                                    1, 1, 0,-1,-1/)  !iside=4

    integer, dimension(num_area,4), parameter :: idx_shift = RESHAPE(idx_shift_tmp,(/num_area,4/))
    integer, dimension(num_area,4), parameter :: idy_shift = RESHAPE(idy_shift_tmp,(/num_area,4/))

    integer, dimension(4), parameter :: iside_m1 = (/4,1,2,3/)
    integer, dimension(4), parameter :: iside_p1 = (/2,3,4,1/)
    integer, dimension(4), parameter :: iside_p2 = (/3,4,1,2/)
    integer, dimension(4), parameter :: iside_p3 = (/4,1,2,3/)

    integer, dimension(4), parameter :: imin_side = (/1   ,0   ,1   ,1   /)
    integer, dimension(4), parameter :: imax_side = (/nc  ,nc  ,nc  ,nc+1/)
    integer, dimension(4), parameter :: jmin_side = (/1   ,1   ,0   ,1   /)
    integer, dimension(4), parameter :: jmax_side = (/nc+1,nc  ,nc  ,nc  /)



    integer :: iur,jur,ilr,jlr,iul,jul,ill,jll

    ib = fvm%cubeboundary 
    flux_sum(0:nc+1,1:nc+1,1) = fvm%se_flux(0:nc+1,0:nc  ,3,ilev)-fvm%se_flux(0:nc+1,1:nc+1,1,ilev)
    flux_sum(1:nc+1,0:nc+1,2) = fvm%se_flux(0:nc  ,0:nc+1,2,ilev)-fvm%se_flux(1:nc+1,0:nc+1,4,ilev)

    !        
    ! Degenerate case ("two departure points")
    !
    !           ||  |                        || no change in this situation ||  no change in this situation
    !           ||  |                        ||                             ||
    !           ||--------                   ||----------                   ||----------
    !           ||  |                        ||                             ||  
    ! =======================      =======================         =====================
    !       |   ||                       |   ||                             ||
    !  -----|---||                 ------|---||                    ---------||
    !       |   ||                       |   ||                             ||
    !       |   ||                       |   ||                             ||
    !
    !
    where (flux_sum(0:nc,1:nc+1,1)*flux_sum(1:nc+1,1:nc+1,1)<0.0D0.and.flux_sum(1:nc+1,0:nc,2)*flux_sum(1:nc+1,1:nc+1,2)<0.0D0)
       degenerate(:,:) = 0
    elsewhere
       degenerate(:,:) = 1
    end where

    if (ib>0) then
       if (ib==swest) degenerate(1   ,1   ) = 1 
       if (ib==nwest) degenerate(1   ,nc+1) = 1 
       if (ib==neast) degenerate(nc+1,nc+1) = 1 
       if (ib==seast) degenerate(nc+1,1   ) = 1 
    end if
    
    do j=1,nc+1
       do i=1,nc+1
          do sgn=-1,1,2
             if (&
                  sgn*flux_sum(i-1,j,1)<0.0D0.and.sgn*flux_sum(i,j-1,2)>0.0D0.and.&
                  sgn*flux_sum(i  ,j,1)>0.0D0.and.sgn*flux_sum(i,j  ,2)<0.0D0) then
                circular_flow(i,j) = 0
             else
                circular_flow(i,j) = 1
             end if
          end do
       end do
    end do
    !
    ! wrap around corners
    !
    if (ib==nwest) then      
       flux_sum(0,nc+1,1) = fvm%se_flux(0,nc,3,ilev)-fvm%se_flux(1,nc+1,4,ilev)
       flux_sum(1,nc+1,2) = fvm%se_flux(0,nc,3,ilev)-fvm%se_flux(1,nc+1,4,ilev) 

       i=1;j=nc+1;
       circular_flow(i,j) = 1
       do sgn=-1,1,2
          if (&
               sgn*flux_sum(i,j-1,2)>0.0D0.and.&
               sgn*flux_sum(i  ,j,1)>0.0D0.and.sgn*flux_sum(i,j  ,2)<0.0D0) then
             circular_flow(i,j) = 0
          end if
       end do
    else if (ib==swest) then      
       flux_sum(0,1,1) = fvm%se_flux(1,0,4,ilev)-fvm%se_flux(0,1,1,ilev)
       flux_sum(1,0,2) = fvm%se_flux(0,1,1,ilev)-fvm%se_flux(1,0,4,ilev)
       i=1;j=1;
       circular_flow(i,j) = 1
       do sgn=-1,1,2
          if (&
               sgn*flux_sum(i-1,j,1)<0.0D0.and.&
               sgn*flux_sum(i  ,j,1)>0.0D0.and.sgn*flux_sum(i,j  ,2)<0.0D0) then
             circular_flow(i,j) = 0
          end if
       end do
    else if (ib==neast) then      
       flux_sum(nc+1,nc+1,1) = fvm%se_flux(nc+1,nc,3,ilev)-fvm%se_flux(nc,nc+1,2,ilev)
       flux_sum(nc+1,nc+1,2) = fvm%se_flux(nc,nc+1,2,ilev)-fvm%se_flux(nc+1,nc,3,ilev) 
       i=nc+1;j=nc+1;
       circular_flow(i,j) = 1
       do sgn=-1,1,2
          if (&
               sgn*flux_sum(i-1,j,1)<0.0D0.and.sgn*flux_sum(i,j-1,2)>0.0D0.and.&
               sgn*flux_sum(i,j  ,2)<0.0D0) then
             circular_flow(i,j) = 0
          end if
       end do
    else if (ib==seast) then      
       flux_sum(nc+1,1   ,1) = fvm%se_flux(nc,0,2,ilev)-fvm%se_flux(nc+1,1,1,ilev)
       flux_sum(nc+1,0   ,2) = fvm%se_flux(nc,0,2,ilev)-fvm%se_flux(nc+1,1,1,ilev) 
       i=nc+1;j=1;
       circular_flow(i,j) = 1
       do sgn=-1,1,2
          if (&
               sgn*flux_sum(i-1,j,1)<0.0D0.and.sgn*flux_sum(i,j-1,2)>0.0D0.and.&
               sgn*flux_sum(i,j  ,2)<0.0D0) then
             circular_flow(i,j) = 0
          end if
       end do
    end if
    illcond = circular_flow*degenerate
    !
    !
    !
    !
    do iside=1,4
       do j=jmin_side(iside),jmax_side(iside)
          do i=imin_side(iside),imax_side(iside)
             if (fvm%se_flux(i,j,iside,ilev)>eps) then
                iur = i+idx_shift(4,iside); jur = j+idy_shift(4,iside) !(i,j) index of upper right quadrant
                ilr = i+idx_shift(5,iside); jlr = j+idy_shift(5,iside) !(i,j) index of lower left  quadrant
                iul = i+idx_shift(2,iside); jul = j+idy_shift(2,iside) !(i,j) index of upper right quadrant
                ill = i+idx_shift(1,iside); jll = j+idy_shift(1,iside) !(i,j) index of lower left  quadrant

                !iside=1
                if (iside==1) then
                displ(0,i,j,iside) = -flux_sum   (i  ,j  ,1)*illcond(i,j)     !center left
                displ(1,i,j,iside) = -flux_sum   (i  ,j  ,1)*illcond(i+1,j)   !center right
                displ(2,i,j,iside) =  flux_sum   (i+1,j  ,2)*illcond(i+1,j)   !c2
                displ(3,i,j,iside) = -flux_sum   (i  ,j  ,2)*illcond(i  ,j)   !c3
                displ(4,i,j,iside) = -flux_sum   (i+1,j  ,1)*illcond(i+1,j)   !r1
                displ(5,i,j,iside) = -flux_sum   (i+1,j-1,2)*illcond(i+1,j)   !r2
                displ(6,i,j,iside) = -flux_sum   (i-1,j  ,1)*illcond(i  ,j)   !l1
                displ(7,i,j,iside) =  flux_sum   (i  ,j-1,2)*illcond(i  ,j)   !l2
                end if
                if (iside==2) then
                !iside=2
                displ(0,i,j,iside) =  flux_sum   (i+1,j  ,2)*illcond(i+1,j  )     !center left
                displ(1,i,j,iside) =  flux_sum   (i+1,j  ,2)*illcond(i+1,j+1)   !center right
                displ(2,i,j,iside) =  flux_sum   (i  ,j+1,1)*illcond(i+1,j+1)   !c2
                displ(3,i,j,iside) = -flux_sum   (i  ,j  ,1)*illcond(i+1,j  )   !c3
                displ(4,i,j,iside) =  flux_sum   (i+1,j+1,2)*illcond(i+1,j+1)   !r1 
                displ(5,i,j,iside) = -flux_sum   (i+1,j+1,1)*illcond(i+1,j+1)   !r2
                displ(6,i,j,iside) =  flux_sum   (i+1,j-1,2)*illcond(i+1,j)   !l1
                displ(7,i,j,iside) =  flux_sum   (i+1,j  ,1)*illcond(i+1,j)   !l2
                end if
                !iside=3
                if (iside==3) then
                displ(0,i,j,iside) =  flux_sum   (i  ,j+1,1)*illcond(i+1,j+1)     !center left
                displ(1,i,j,iside) =  flux_sum   (i  ,j+1,1)*illcond(i  ,j+1)   !center right
                displ(2,i,j,iside) = -flux_sum   (i  ,j  ,2)*illcond(i  ,j+1)   !c2
                displ(3,i,j,iside) =  flux_sum   (i+1,j  ,2)*illcond(i+1,j+1)   !c3
                displ(4,i,j,iside) =  flux_sum   (i-1,j+1,1)*illcond(i  ,j+1)   !r1
                displ(5,i,j,iside) =  flux_sum   (i  ,j+1,2)*illcond(i  ,j+1)   !r2
                displ(6,i,j,iside) =  flux_sum   (i+1,j+1,1)*illcond(i+1,j+1)   !l1
                displ(7,i,j,iside) = -flux_sum   (i+1,j+1,2)*illcond(i+1,j+1)   !l2
                end if
                if (iside==4) then
                !iside=4
                displ(0,i,j,iside) = -flux_sum   (i  ,j  ,2)*illcond(i  ,j+1)     !center left
                displ(1,i,j,iside) = -flux_sum   (i  ,j  ,2)*illcond(i  ,j  )   !center right
                displ(2,i,j,iside) = -flux_sum   (i  ,j  ,1)*illcond(i  ,j  )   !c2
                displ(3,i,j,iside) =  flux_sum   (i  ,j+1,1)*illcond(i  ,j+1)   !c3
                displ(4,i,j,iside) = -flux_sum   (i  ,j-1,2)*illcond(i  ,j  )   !r1
                displ(5,i,j,iside) =  flux_sum   (i-1,j  ,1)*illcond(i  ,j  )   !r2
                displ(6,i,j,iside) = -flux_sum   (i  ,j+1,2)*illcond(i  ,j+1)   !l1
                displ(7,i,j,iside) = -flux_sum   (i-1,j+1,1)*illcond(i  ,j+1)   !l2
                end if

!                if (ilev<5.and.i==2.and.j==2) &
!                     write(*,*) "displ",i,j,iside,ilev,displ(:,i,j,iside)


                base_vtx(:,1,i,j,iside) = fvm%flux_vertex_cart(i  ,j  ,:,iside          )       !vertex center left
                base_vtx(:,2,i,j,iside) = fvm%flux_vertex_cart(i  ,j  ,:,iside_p1(iside))       !vertex center right
                base_vtx(:,3,i,j,iside) = fvm%flux_vertex_cart(iur,jur,:,iside          )       !vertex upper right
                base_vtx(:,4,i,j,iside) = fvm%flux_vertex_cart(ilr,jlr,:,iside_p3(iside))       !vertex lower right
                base_vtx(:,5,i,j,iside) = fvm%flux_vertex_cart(iul,jul,:,iside_p1(iside))       !vertex upper left
                base_vtx(:,6,i,j,iside) = fvm%flux_vertex_cart(ill,jll,:,iside_p2(iside))       !vertex lower left

                base_vec(:, 1,i,j,iside) = fvm%flux_vec    (:,i  ,j  ,iside          )      !vector center
                base_vec(:, 2,i,j,iside) = fvm%flux_vec    (:,i  ,j  ,iside_p1(iside))      !vector center right
                base_vec(:, 3,i,j,iside) = fvm%flux_vec    (:,i  ,j  ,iside_p3(iside))      !vector center left
                base_vec(:, 4,i,j,iside) = fvm%flux_vec    (:,iur,jur,iside          )      !vector upper right 1
                base_vec(:, 5,i,j,iside) = fvm%flux_vec    (:,iur,jur,iside_p3(iside))      !vector upper right 2
                base_vec(:, 6,i,j,iside) = fvm%flux_vec    (:,ilr,jlr,iside_p3(iside))      !vector lower right 1
                base_vec(:, 7,i,j,iside) = fvm%flux_vec    (:,ilr,jlr,iside_p2(iside))      !vector lower right 2
                base_vec(:, 8,i,j,iside) = fvm%flux_vec    (:,iul,jul,iside          )      !vector upper left 1
                base_vec(:, 9,i,j,iside) = fvm%flux_vec    (:,iul,jul,iside_p1(iside))      !vector upper left 2
                base_vec(:,10,i,j,iside) = fvm%flux_vec    (:,ill,jll,iside_p1(iside))      !vector lower left 1
                base_vec(:,11,i,j,iside) = fvm%flux_vec    (:,ill,jll,iside_p2(iside))      !vector lower left 2

                do iarea=1,5
                   idx(1,iarea,i,j,iside) = i+idx_shift(iarea,iside)
                   idx(2,iarea,i,j,iside) = j+idy_shift(iarea,iside)
                end do
             else
                displ(:,i,j,iside) = 9D99!dbgxxx
             end if
          end do
       end do
    end do
    !
    ! wrap around corners here
    !

  end subroutine define_swept_areas


  !
  ! Notation conventions used in define_area subroutines
  !
  !
  !
  !   ^    ||--->   ^   <---||    ^
  !  /|\   || 3    /|\    2 ||   /|\
  !   | 6  ||     1 |       ||    | 4
  !   |    ||       |       ||    | 
  ! =================================
  !        ||               ||    
  !        ||               ||    
  !      7 ||               || 5   
  !    <---||               ||--->    
  !

  subroutine define_area1_area2(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,&
       x_start, dgam_vec)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8)                   , intent(inout):: x_start, dgam_vec


    real (kind=real_kind)    , dimension(2,3) :: xdep !departure points
    real (kind=real_kind)                     :: gamma
    integer :: iarea


    REAL(KIND=real_kind) :: dcross, tmp(2), static_integral, xtmp(2),xtmp2(2)
    !
    !
    !        ||-----        ||
    !       /||             ||
    !      / ||             ||
    !  ===X========================= 
    !     | /||             ||
    !     |/ ||             ||
    !     *  ||             ||
    !
    !
    ! crossing X
    if (SUM(ABS(base_vec(:,9,i,j,iside))).NE.0) then
       gamma = displ(0,i,j,iside)*displ(7,i,j,iside)/(displ(0,i,j,iside)-displ(6,i,j,iside))
       gamma = MAX(MIN(gamma,displ(7,i,j,iside),-displ(3,i,j,iside)),0.0D0)
    else
       !
       ! corner case
       !
       gamma=displ(0,i,j,iside)
    end if

    
    xdep    (:,1) = base_vtx(:, 6,i,j,iside)+displ(7,i,j,iside)*base_vec(:,10,i,j,iside)-displ(6,i,j,iside)*base_vec(:,11,i,j,iside)
    x_start (:,1) = base_vtx(:, 6,i,j,iside)
    dgam_vec(:,1) = base_vec(:,10,i,j,iside)*gamma
    
    xdep(:,2) = base_vtx(:,2,i,j,iside)+displ(1,i,j,iside)*base_vec(:, 1,i,j,iside)+displ(2,i,j,iside)*base_vec(:, 2,i,j,iside)
    
    iarea                  = 1
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 1
    
    x_static (:,1,iarea) = base_vtx(:,6,i,j,iside)       !static
    dx_static(:,1,iarea) = xdep(:,1)-x_static(:,1,iarea) !static
    
    xtmp(:        ) = x_start(:,1)+dgam_vec(:,1)
    x   (:,1,iarea) = xdep(:,1)                  !static
    dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic
    
    x (:,2,iarea) = xtmp(:)                      !dynamic
    dx(:,2,iarea) = x_static(:,1,iarea)-xtmp(:)  !dynamic
    !
    !
    !
    iarea                  = 2
    num_seg       (iarea)  = 3
    
    x_start (:,2) = base_vtx(:,5,i,j,iside)
    dgam_vec(:,2) = base_vec(:,9,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,2)+dgam_vec(:,2)
    
    x_start (:,3) = base_vtx(:,5,i,j,iside)
    dgam_vec(:,3) = base_vec(:,8,i,j,iside)*displ(0,i,j,iside)
    xtmp2   (:  ) = x_start(:,3)+dgam_vec(:,3)
    
    x   (:,1,iarea) = base_vtx(:,5,i,j,iside) !static
    dx  (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic
    
    x (:,2,iarea) = xtmp (:)          !dynamic
    dx(:,2,iarea) = xtmp2(:)-xtmp(:)  !dynamic
    
    x (:,3,iarea) = xtmp2(:)                !dynamic
    dx(:,3,iarea) = x(:,1,iarea)-xtmp2(:)   !dynamic
  end subroutine define_area1_area2


  subroutine define_area2(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,x_start, dgam_vec)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8)                   , intent(inout):: x_start, dgam_vec


    real (kind=real_kind)    , dimension(2,3) :: xdep !departure points
    real (kind=real_kind)                     :: gamma
    integer :: iarea


    REAL(KIND=real_kind) :: dcross, tmp(2), static_integral, xtmp(2),xtmp2(2)
    ! *: xdep(:,1)
    ! x: xtmp
    !
    !      2 ||             ||      
    !     *--x              ||      
    !     1\3||1            ||      
    !       \||             ||      
    !  =============================
    !        ||             ||      
    !
    !
    ! compute departure points (xdep(1) is left; xdep(3) is right and xdep(2) is midway
    !
    xdep(:,1) = base_vtx(:,5,i,j,iside)+&
         MAX(0.0D0,displ(6,i,j,iside))*base_vec(:,8,i,j,iside)-displ(3,i,j,iside)*base_vec(:,9,i,j,iside)
    x_start (:,1) = base_vtx(:,5,i,j,iside)
    gamma         = displ(0,i,j,iside)
    dgam_vec(:,1) = base_vec(:,8,i,j,iside)*gamma

    iarea                  = 2
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 1
                  
    x_static (:,1,iarea) = base_vtx(:,5,i,j,iside)       !static  - line 1
    dx_static(:,1,iarea) = xdep(:,1)-x_static(:,1,iarea) !static  - line 1
    
    xtmp     (:        ) = x_start(:,1)+dgam_vec(:,1)
    x        (:,1,iarea) = xdep(:,1)                  !static  - line 2
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic - line 2
    
    x        (:,2,iarea) = xtmp(:)                     !dynamic  - line 3
    dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 3
  end subroutine define_area2


  subroutine define_area3_left(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, &
       num_seg, num_seg_static,x_start, dgam_vec)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8)                   , intent(inout):: x_start, dgam_vec


    real (kind=real_kind)    , dimension(2,3) :: xdep !departure points
    real (kind=real_kind)                     :: gamma
    integer :: iarea


    REAL(KIND=real_kind) :: dcross, tmp(2), static_integral, xtmp(2),xtmp2(2)

    ! iarea = 3
    !-------------------------------------------------------------------------------------------
    !
    !          xtmp         xdep(2)
    !           |x-----2------*   ||
    !           ||             \  ||
    !           |1              3 ||
    !           ||               \||
    !        ===========4==============
    !
    !
    xdep(:,2) = base_vtx(:,2,i,j,iside)+displ(1,i,j,iside)*base_vec(:,1,i,j,iside)&
         +MAX(0.0D0,displ(2,i,j,iside))*base_vec(:,2,i,j,iside)
    x_start (:,4) = base_vtx(:,1,i,j,iside)
    gamma         = displ(0,i,j,iside)
    dgam_vec(:,4) = base_vec(:,1,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,4)+dgam_vec(:,4)

    iarea                  = 3
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 2
                  
    x_static (:,1,iarea) = xdep(:,2)                         !static  - line 3
    dx_static(:,1,iarea) = base_vtx(:,2,i,j,iside)-xdep(:,2) !static  - line 3

    x_static (:,2,iarea) = base_vtx(:,2,i,j,iside)                         !static  - line 4
    dx_static(:,2,iarea) = base_vtx(:,1,i,j,iside)-base_vtx(:,2,i,j,iside) !static  - line 4
    
    x        (:,1,iarea) = base_vtx(:,1,i,j,iside)    !static  - line 1
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic - line 1
                   
    x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
    dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
  end subroutine define_area3_left

  subroutine define_area3_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
       num_seg_static,x_start, dgam_vec)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8)                   , intent(inout):: x_start, dgam_vec


    real (kind=real_kind)    , dimension(2,3) :: xdep !departure points
    real (kind=real_kind)                     :: gamma
    integer :: iarea

    REAL(KIND=real_kind) :: dcross, tmp(2), static_integral, xtmp(2),xtmp2(2)
    !
    !
    !        ||  *-----2----||\       
    !        || /1         3|| \      
    !        ||/      4     ||       
    !  ============================= 
    !        ||             ||
    !        ||             ||
    !        ||             ||
    !
    xdep(:,1) = base_vtx(:,1,i,j,iside)+displ(0,i,j,iside)*base_vec(:,1,i,j,iside)&
         +MAX(0.0D0,displ(3,i,j,iside))*base_vec(:,3,i,j,iside)
    x_start (:,5) = base_vtx(:,2,i,j,iside)
    gamma         = displ(1,i,j,iside)
    dgam_vec(:,5) = base_vec(:,1,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,5)+dgam_vec(:,5)
    
    iarea                  = 3
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 2
    
    x_static (:,1,iarea) = base_vtx(:,1,i,j,iside)           !static  - line 1
    dx_static(:,1,iarea) = xdep(:,1)-base_vtx(:,1,i,j,iside) !static  - line 1
    
    x_static (:,2,iarea) = base_vtx(:,2,i,j,iside)                         !static  - line 4
    dx_static(:,2,iarea) = base_vtx(:,1,i,j,iside)-base_vtx(:,2,i,j,iside) !static  - line 4
    
    x        (:,1,iarea) = xdep(:,1)            !static  - line 2
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea) !dynamic - line 2
    
    x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
    dx       (:,2,iarea) = x_static(:,2,iarea)-xtmp(:) !dynamic - line 2   
  end subroutine define_area3_right


  subroutine define_area3_left_right(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
       num_seg_static,x_start, dgam_vec)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8)                   , intent(inout):: x_start, dgam_vec
    
    
    real (kind=real_kind)    , dimension(2,3) :: xdep !departure points
    real (kind=real_kind)                     :: gamma
    integer :: iarea
    
    REAL(KIND=real_kind) :: dcross, tmp(2), static_integral, xtmp(2),xtmp2(2)
    !
    !        ||-------------||       
    !       /||             ||\      
    !        ||             ||      
    !  ============================= 
    !        ||             ||
    !        ||             ||    
    !        ||             ||   
    !
    x_start (:,4) = base_vtx(:,1,i,j,iside)
    x_start (:,5) = base_vtx(:,2,i,j,iside)
    gamma         = displ(0,i,j,iside)
    dgam_vec(:,4) = base_vec(:,1,i,j,iside)*gamma
    dgam_vec(:,5) = base_vec(:,1,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,4)+dgam_vec(:,4)
    xtmp2   (:  ) = x_start(:,5)+dgam_vec(:,5)
    
    iarea                  = 3
    num_seg       (iarea)  = 3
    num_seg_static(iarea)  = 1
    
    x_static (:,1,iarea) = base_vtx(:,2,i,j,iside)                         !static
    dx_static(:,1,iarea) = base_vtx(:,1,i,j,iside)-base_vtx(:,2,i,j,iside) !static
    
    x        (:,1,iarea) = base_vtx(:,1,i,j,iside)    !static  
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)       !dynamic 
    
    x        (:,2,iarea) = xtmp (:)         !dynamic 
    dx       (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic 
    
    x        (:,3,iarea) = xtmp2(:)              !dynamic
    dx       (:,3,iarea) = x_start(:,5)-xtmp2(:) !dynamic
  end subroutine define_area3_left_right

  subroutine define_area4_area5(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
       num_seg_static,x_start, dgam_vec)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8)                   , intent(inout):: x_start, dgam_vec
    
    
    real (kind=real_kind)    , dimension(2,3) :: xdep !departure points
    real (kind=real_kind)                     :: gamma
    integer :: iarea
    
    REAL(KIND=real_kind) :: dcross, tmp(2), static_integral, xtmp(2),xtmp2(2)
    !
    !        ||     --------||       
    !        ||             ||\      
    !        ||             || \     
    !  ============================= 
    !        ||             ||\ |    
    !        ||             || \|    
    !        ||             ||  *    
    !
    !
    ! iarea  = 4
    !
    iarea                  = 4
    num_seg       (iarea)  = 3

    if (SUM(ABS(base_vec(:,5,i,j,iside))).NE.0) then        
       gamma = displ(1,i,j,iside)*displ(5,i,j,iside)/(displ(1,i,j,iside)-displ(4,i,j,iside))
       gamma = MAX(MIN(gamma,displ(5,i,j,iside),-displ(2,i,j,iside)),0.0D0)
    else
       !
       ! corner case
       !
       gamma = displ(1,i,j,iside)
    end if
    x_start (:,6) = base_vtx(:,3,i,j,iside)
    dgam_vec(:,6) = base_vec(:,4,i,j,iside)*displ(1,i,j,iside)
    xtmp    (:  ) = x_start(:,6)+dgam_vec(:,6)
    x_start (:,7) = base_vtx(:,3,i,j,iside)
    dgam_vec(:,7) = base_vec(:,5,i,j,iside)*gamma
    xtmp2   (:  ) = x_start(:,7)+dgam_vec(:,7)
    
    x        (:,1,iarea) = base_vtx(:,3,i,j,iside)!static   -line 1
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
    
    x        (:,2,iarea) = xtmp(:)          !dynamic -line 2
    dx       (:,2,iarea) = xtmp2(:)-xtmp(:) !dynamic - line 2
    
    x        (:,3,iarea) = xtmp2(:)               !static   -line 1
    dx       (:,3,iarea) = x(:,1,iarea)-xtmp2(:)  !dynamic - line 1
    !
    !iarea = 5
    !
    xdep(:,1) = base_vtx(:,4,i,j,iside)+displ(5,i,j,iside)*base_vec(:,6,i,j,iside)&
         -displ(4,i,j,iside)*base_vec(:,7,i,j,iside)

!    if (SUM(ABS(base_vec(:,5,i,j,iside))).NE.0) then        
!       gamma = displ(1,i,j,iside)*displ(5,i,j,iside)/(displ(1,i,j,iside)-displ(4,i,j,iside))
!       gamma = MAX(MIN(gamma,displ(5,i,j,iside),-displ(2,i,j,iside)),0.0D0)
!    else
!       !
!       ! corner case
!       !
!       gamma = displ(1,i,j,iside)
!    end if
    
    x_start (:,8) = base_vtx(:,4,i,j,iside)
    dgam_vec(:,8) = base_vec(:,6,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,8)+dgam_vec(:,8)
    
    iarea                  = 5
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 1
    
    x        (:,1,iarea) = base_vtx(:,4,i,j,iside)!static   -line 1
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)   !dynamic - line 1
    
    x_static (:,1,iarea) = xdep(:,1)                        !static - line 1
    dx_static(:,1,iarea) = x(:,1,iarea)-x_static(:,1,iarea) !static - line 1
    
    x        (:,2,iarea) = xtmp(:)                     !dynamic -line 2
    dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
  end subroutine define_area4_area5


  subroutine define_area4(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, &
       num_seg_static,x_start, dgam_vec)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8)                   , intent(inout):: x_start, dgam_vec
    
    
    real (kind=real_kind)    , dimension(2,3) :: xdep !departure points
    real (kind=real_kind)                     :: gamma
    integer :: iarea
    
    REAL(KIND=real_kind) :: dcross, tmp(2), static_integral, xtmp(2),xtmp2(2)                   

    iarea                  = 4
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 1
    
    xdep(:,1) = base_vtx(:,3,i,j,iside)+MAX(0.0D0,displ(4,i,j,iside))*base_vec(:,4,i,j,iside)&
         -displ(2,i,j,iside)*base_vec(:,5,i,j,iside)
    x_start (:,6) = base_vtx(:,3,i,j,iside)
    gamma         = displ(1,i,j,iside)
    dgam_vec(:,6) = base_vec(:,4,i,j,iside)*gamma
    xtmp    (:  ) = x_start(:,6)+dgam_vec(:,6)
    
    x_static (:,1,iarea) = xdep(:,1)                         !static 
    dx_static(:,1,iarea) = base_vtx(:,3,i,j,iside)-xdep(:,1) !static
    
    x        (:,1,iarea) = base_vtx(:,3,i,j,iside) !static  - line 2
    dx       (:,1,iarea) = xtmp(:)-x(:,1,iarea)    !dynamic - line 2
    
    x        (:,2,iarea) = xtmp(:)                     !dynamic  -line 2
    dx       (:,2,iarea) = x_static(:,1,iarea)-xtmp(:) !dynamic - line 2
  end subroutine define_area4

  subroutine define_area3_center(i,j,iside,displ,base_vec,base_vtx,x, dx, x_static, dx_static, num_seg, num_seg_static,&
       x_start, dgam_vec,se_flux_center)
    implicit none
    integer, intent(in) :: i,j,iside
    integer, parameter :: num_area=5, num_sides=4, imin= 0, imax=nc+1
    real (kind=real_kind)    , dimension(0:7       , imin:imax,imin:imax,num_sides), intent(inout) :: displ
    integer (kind=real_kind) , dimension(1:2,11    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vec
    real (kind=real_kind)    , dimension(1:2, 6    , imin:imax,imin:imax,num_sides), intent(inout) :: base_vtx
    integer, parameter :: num_seg_max=5
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(inout) :: x, dx, x_static, dx_static
    integer             , dimension(num_area)              , intent(inout) :: num_seg, num_seg_static
    REAL(KIND=real_kind), dimension(2,8)                   , intent(inout):: x_start, dgam_vec
    REAL(KIND=real_kind) , intent(in   ):: se_flux_center
    
    
    real (kind=real_kind)    , dimension(2,3) :: xdep !departure points
    real (kind=real_kind)                     :: gamma
    integer :: iarea
    
    REAL(KIND=real_kind) :: dcross, tmp(2), static_integral, xtmp(2),xtmp2(2)                   

    !                 xdep(2)
    !                 ______X______
    !        ||      /             \      ||        
    !        ||  *--/               \--*  ||        
    !        || /xdep(1)         xdep(3)\ ||        
    !        ||/                         \||        
    !  ========================================  
    !        ||                           ||        
    !
    !
    ! compute departure points (xdep(1) is left; xdep(3) is right and xdep(2) is midway
    !
    
    xdep(:,1) = base_vtx(:,1,i,j,iside)+&
         displ(0,i,j,iside)*base_vec(:,1,i,j,iside)+displ(3,i,j,iside)*base_vec(:,3,i,j,iside)
    xdep(:,3) = base_vtx(:,2,i,j,iside)+&
         displ(1,i,j,iside)*base_vec(:,1,i,j,iside)+displ(2,i,j,iside)*base_vec(:,2,i,j,iside)
    xdep(:,2) = 0.5D0*(xdep(:,1)+xdep(:,3))
 
!
! my new code
!   
!    gamma= se_flux_center
!    x_start(:,1) = base_vec(:,3,i,j,iside)*ABS((xdep(:,2)-base_vtx(:,1,i,j,iside)))+&
!         base_vtx(:,1,i,j,iside) !xdep(2) - midway between departure points projected to side 1 
!    
!    dgam_vec(:,1) = base_vec(:,1,i,j,iside)*(SUM(ABS(base_vtx(:,2,i,j,iside)-base_vtx(:,1,i,j,iside))))

    !
    ! same as old code
    !
    !start
    gamma= se_flux_center
    x_start(:,1) = xdep(:,2)
    dgam_vec(1,1) = -(xdep(2,3)-xdep(2,1))
    dgam_vec(2,1) =  (xdep(1,3)-xdep(1,1))
    !end



    xdep(:,2)     = x_start(:,1)+gamma*dgam_vec(:,1)
    iarea                  = 3
    num_seg       (iarea)  = 2
    num_seg_static(iarea)  = 3
    
    !                 ______X______
    !        ||    2 /             \ 3    ||        
    !        ||  *--/               \--*  ||        
    !        || /                       \ ||        
    !        ||/ 1          5           4\||        
    !  ========================================  
    !        ||                           ||        
    !
    x_static (:,1,iarea) = base_vtx(:,1,i,j,iside)       !static  - line 1
    dx_static(:,1,iarea) = xdep(:,1)-x_static(:,1,iarea) !static  - line 1
    
    x        (:,1,iarea) = xdep(:,1)                     !static  - line 2
    dx       (:,1,iarea) = xdep(:,2)-x(:,1,iarea)        !dynamic - line 2
    
    x        (:,2,iarea) = xdep(:,2)                     !dynamic - line 3
    dx       (:,2,iarea) = xdep(:,3)-x(:,2,iarea)        !dynamic - line 3
    
    x_static (:,2,iarea) = xdep(:,3)                                  !static  - line 4
    dx_static(:,2,iarea) = base_vtx(:,2,i,j,iside)-x_static(:,2,iarea)!static  - line 4
    
    x_static (:,3,iarea) = base_vtx(:,2,i,j,iside)                         !static - line 5
    dx_static(:,3,iarea) = base_vtx(:,1,i,j,iside)-base_vtx(:,2,i,j,iside) !static - line 5

  end subroutine define_area3_center
  !
  ! UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST
  ! UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST
  ! UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST-UNIT-TEST
  !
  ! EVERYTHING BELOW IS DEBUGGING SUBROUTINES
  !
  !
  subroutine write_cells(fvm,x,dx,num_area,num_seg_max,num_seg,ie,i_in,j_in,iside,ilev,flowcase)
    implicit none
    type (fvm_struct), intent(in) :: fvm
    integer              , intent(in) :: ie, i_in, j_in, iside, ilev, num_area, num_seg_max, flowcase
    REAL(KIND=real_kind), dimension(2,num_seg_max,num_area), intent(in) :: x, dx
    integer             , dimension(num_area              ), intent(in) :: num_seg


    character(len=1024)               :: grid_dat, grid_fname,gnuplot, output_fname
    character(len=1024)               :: displ_dat,displ_fname, swept_flux_dat, swept_flux_fname
    character(len=1024)               :: swept_flux_old_dat, swept_flux_old_fname

    character(len=100)               :: flux_side_string,ie_string,i_string,j_string,xyflux_string,lev_string
    character(len=100)               :: tmp_string,tmp2_string,displ_string
    character(len=5)                 :: which_flux_string, folder_string

    integer :: ii, jj, ss, iarea, iseg, i,j,i_center,j_center
    integer, dimension(4) :: side=(/1,2,3,4/), side_sgn=(/1,1,1,1/)
    logical :: lcell
    !
    ! unique identifier for flux side
    !
    return
    i=i_in; j=j_in
    lcell=.true.
    if (iside==1) then
       i=i_in; j=j_in;  write(which_flux_string,"(A5)") "yflux" 
!       if (i==i_plot.and.j==j_plot) lcell=.true.
    else if (iside==2) then
       i=i_in+1; j=j_in;  write(which_flux_string,"(A5)") "xflux" 
!       if (i-1==i_plot.and.j==j_plot) lcell=.true.
    else if (iside==3) then
       i=i_in; j=j_in+1;  write(which_flux_string,"(A5)") "yflux" 
!       if (j-1==j_plot.and.i==i_plot) lcell=.true.
    else if (iside==4) then
       i=i_in; j=j_in;  write(which_flux_string,"(A5)") "xflux" 
!       if (i==i_plot.and.j==j_plot) lcell=.true.
    else
       write(*,*) "iside out of range",iside
       stop
    end if
    if (ilev<ilev_limit.and.ie<ie_limit.and.lcell) then
    write(tmp_string ,"(I10)") ie  ;write(ie_string ,"(A2)") "ie";
    write(tmp2_string,"(I10)") ilev;write(lev_string,"(A4)") "-lev"
    write(i_string,"(A2,I1,A1)") "-i",i,"-"; write(j_string,"(A1,I1,A1)") "j",j,"-"; 
    write(flux_side_string,*) &
         trim(adjustl(ie_string )),trim(adjustl(tmp_string )),&
         trim(adjustl(lev_string)),trim(adjustl(tmp2_string)),&
         trim(adjustl(i_string)),trim(adjustl(j_string)),&
         trim(adjustl(which_flux_string))

    write(swept_flux_dat,*) trim(adjustl(flux_side_string)),".dat"
    write(folder_string,"(A5)") "data/"; write(swept_flux_fname,*) folder_string,trim(adjustl(flux_side_string)),".dat"
    write (gnuplot ,*) folder_string,trim(adjustl(flux_side_string)),".gp"
    write (grid_dat,*) folder_string,trim(adjustl(tmp_string)),"-grid.dat"
    write (grid_fname,*) trim(adjustl(tmp_string)),"-grid.dat"

    write(displ_string,*) &
         trim(adjustl(ie_string )),trim(adjustl(tmp_string )),&
         trim(adjustl(lev_string)),trim(adjustl(tmp2_string)),&
         "-displ.dat"
    write (displ_dat,*) folder_string,trim(adjustl(displ_string))
    write (displ_fname,*) trim(adjustl(displ_string))
    write (output_fname ,*) trim(adjustl(flux_side_string)),".eps"
    write(swept_flux_old_dat,*) trim(adjustl(flux_side_string)),"-old.dat"
    write(swept_flux_old_fname,*) folder_string,trim(adjustl(flux_side_string)),"-old.dat"
    ! 
    ! plot grid and displacements for entire elements
    !
    OPEN (UNIT=2, STATUS='REPLACE', FILE=trim(adjustl(grid_dat )))
    OPEN (UNIT=4, STATUS='REPLACE', FILE=trim(adjustl(displ_dat)))
    i_center=2
    j_center=2
    do jj=-1,1
       do ii=-1,1
          do ss=0,3
             side=cshift(side,ss)
             write(2,*) fvm%flux_vertex_cart(i_center+ii,j_center+jj,:,side(1))
             write(2,*) fvm%flux_vertex_cart(i_center+ii,j_center+jj,:,side(2))
             write(2,*) " "
             
             write(4,*) fvm%flux_vertex_cart(i_center+ii,j_center+jj,:,side(1))&
                  +fvm%se_flux(i_center+ii,j_center+jj,side(1),ilev)*&
                                  fvm%flux_vec(:,i_center+ii,j_center+jj,side(1))
             write(4,*) fvm%flux_vertex_cart(i_center+ii,j_center+jj,:,side(2))&
                  +fvm%se_flux(i_center+ii,j_center+jj,side(1),ilev)*&
                   fvm%flux_vec(:,i_center+ii,j_center+jj,side(1))
             write(4,*) " "
          end do
       end do
    end do
    CLOSE(2)
    CLOSE(4)
    !       OPEN (UNIT=3, STATUS='REPLACE', FILE=side_dat)
    !       side=(/1,2,3,4/)
    !       side=cshift(side,iside-1)
    !       write(3,*) fvm%flux_vertex_cart(i,j,:,side(1))
    !       write(3,*) fvm%flux_vertex_cart(i,j,:,side(2))
    !       CLOSE(3)
    !
    ! plot swep area
    !
    OPEN (UNIT=5, STATUS='REPLACE', FILE=trim(adjustl(swept_flux_fname)))
    do iarea=1,num_area
       do iseg=1,num_seg(iarea)
             write(5,*) x(:,iseg,iarea)
             write(5,*) x(:,iseg,iarea)+dx(:,iseg,iarea)
             write(5,*) " "
             !          write(5,*) x(2,iseg,iarea), x(2,iseg,iarea)+dx(2,iseg,iarea)
             !          write(5,*) x(1,iseg,iarea), 
             !          write(5,*) x(2,iseg,iarea), x(2,iseg,iarea)+dx(2,iseg,iarea)
          end do
       end do
       close(5)
       
       OPEN (UNIT=99, STATUS='REPLACE', FILE=trim(adjustl(gnuplot)))
       write(99,*) "set terminal eps color solid"
       write(99,"(A24,I2,A6,I3,A1)") "set title ""Flowcase: ",flowcase,", ie=",ie,""""
       write(99,*) "set output """,TRIM(adjustl(output_fname)),""""
       write(99,*) "plot """,TRIM(adjustl(grid_fname        )),""" w l lw 3 lt 2 notitle,"""          &
                            ,TRIM(adjustl(swept_flux_dat    )),""" w lp lw 4 lt 6 ps 0.75 notitle,""" &
                            ,TRIM(adjustl(swept_flux_old_dat)),""" w lp lw 2 lt 3 ps 0.25 notitle,""" &
                            ,TRIM(adjustl(displ_fname   )),""" w l lw 1 lt -1 notitle"

!            TRIM(swept_flux_dat_old),""" w lp lt 4 ps 0.1 notitle"
       CLOSE(99)
    end if
  end subroutine write_cells


  subroutine plot_element_numbers(fvm,elem,nets,nete)
    use control_mod           , only: north, south, east, west, neast, nwest, seast, swest
    use fvm_control_volume_mod, only: n0_fvm
    implicit none
    type (fvm_struct), intent(inout) :: fvm(:)
    type (element_t) , intent(inout) :: elem(:)
    integer :: ie
    integer, intent(in) :: nets,nete
    real (kind=real_kind), dimension(2):: xtmp

    OPEN (UNIT=99, STATUS='REPLACE', FILE="data/elements.gp")
    write(99,*) "set terminal eps color solid"
    write(99,*) "set output ""elements.eps"""

    OPEN (UNIT=4, STATUS='REPLACE', FILE="data/elements.dat")
    do ie=nets,nete
       write(4,*) transform(fvm(ie)%flux_vertex_cart(1 , 1,:,1),NINT(fvm(ie)%flux_orient(1,1,1)))
       write(4,*) transform(fvm(ie)%flux_vertex_cart(nc, 1,:,2),NINT(fvm(ie)%flux_orient(1,1,1)))
       write(4,*) transform(fvm(ie)%flux_vertex_cart(nc,nc,:,3),NINT(fvm(ie)%flux_orient(1,1,1)))
       write(4,*) transform(fvm(ie)%flux_vertex_cart(1 ,nc,:,4),NINT(fvm(ie)%flux_orient(1,1,1)))
       write(4,*) transform(fvm(ie)%flux_vertex_cart(1 , 1,:,1),NINT(fvm(ie)%flux_orient(1,1,1)))
       write(4,*) "   "
       xtmp = 0.5D0*(transform(fvm(ie)%flux_vertex_cart(1 , 1,:,1),NINT(fvm(ie)%flux_orient(1,1,1)))+&
                     transform(fvm(ie)%flux_vertex_cart(nc,nc,:,3),NINT(fvm(ie)%flux_orient(1,1,1))))
       write(99,"(A11,I4,A5,F10.2,A1,F10.2,A25)") "set label """,ie,""" at ",xtmp(1),",",xtmp(2)," center font ""Times, 4"""
    end do
    write(99,*) "plot ""elements.dat"" w l"
    CLOSE(99)             


  end subroutine plot_element_numbers

  function transform(x,ipanel) result(point)
    implicit none
    real (kind=real_kind), intent(in)  :: x(2)
    real (kind=real_kind), dimension(2):: point
    integer, intent(in) :: ipanel
    select case (ipanel)
       case(1)
          point = x
       case(2)
          point(1) = x(1)+2.0D0
          point(2) = x(2)
       case(3)
          point(1) = x(1)+4.0D0
          point(2) = x(2)
       case(4)
          point(1) = x(1)-2.0D0
          point(2) = x(2)
       case(5)
          point(1) = x(1)
          point(2) = x(2)-2.0D0
       case(6)
          point(1) = x(1)
          point(2) = x(2)+2.0D0
       case default
          stop
       end select
  end function transform

  subroutine set_flow_cases(fvm,elem,nets,nete)
    use control_mod           , only: north, south, east, west, neast, nwest, seast, swest
    use fvm_control_volume_mod, only: n0_fvm
    implicit none
    type (fvm_struct), intent(inout) :: fvm(:)
    type (element_t) , intent(inout) :: elem(:)

    integer, intent(in)                         :: nets  ! starting thread element number (private)
    integer, intent(in)                         :: nete  ! ending thread element number   (private)
    integer :: ib,i,j, d1,d2,d3,d4,d5,d6,d7, ie_case, ie

    integer, parameter :: num_area=6
    integer, dimension(num_area*4), parameter :: idx_shift_tmp = (/-1,-1, 0, 1, 1,0,&  !iside=1
                                                                    1, 0, 0, 0, 1,1,&  !iside=2
                                                                    1, 1, 0,-1,-1,0,&  !iside=3
                                                                   -1, 0, 0, 0,-1,-1/)  !iside=4

    integer, dimension(num_area*4), parameter :: idy_shift_tmp = (/-1, 0, 0, 0,-1,-1,&  !iside=1
                                                                   -1,-1, 0, 1, 1,0,&  !iside=2
                                                                    1, 0, 0, 0, 1,1,&  !iside=3
                                                                    1, 1, 0,-1,-1,0/)  !iside=4

    integer, dimension(num_area,4), parameter :: idx_shift = RESHAPE(idx_shift_tmp,(/num_area,4/))
    integer, dimension(num_area,4), parameter :: idy_shift = RESHAPE(idy_shift_tmp,(/num_area,4/))

    integer, dimension(4), parameter :: iside_m1 = (/4,1,2,3/)
    integer, dimension(4), parameter :: iside_p1 = (/2,3,4,1/)
    integer, dimension(4), parameter :: iside_p2 = (/3,4,1,2/)
    integer, dimension(4), parameter :: iside_p3 = (/4,1,2,3/)

    integer :: iur,jur,ilr,jlr,iul,jul,ill,jll,ilc,jlc,iside

    integer   :: ilev,ie1,ie2

    real (kind=real_kind) :: scale, fraction
    !    ib = fvm%cubeboundary 
    ie_case = 1
    i=2; j=2;
    do ie=nets,nete
       elem(ie)%sub_elem_mass_flux(:,:,:,1:nlev) = 0.0D0
    enddo
    fraction = -0.2D0
    if (.false.) then
    do ilev=1,4
       ie=1
       iside=ilev
       do d2=-1,1
          do d3=-1,1
             do d4=-1,1
                do d5=-1,1
                   do d6=-1,1
                      do d7=-1,1


                         iur = i+idx_shift(4,iside); jur = j+idy_shift(4,iside) !(i,j) index of upper right quadrant
                         ilr = i+idx_shift(5,iside); jlr = j+idy_shift(5,iside) !(i,j) index of lower left  quadrant
                         iul = i+idx_shift(2,iside); jul = j+idy_shift(2,iside) !(i,j) index of upper right quadrant
                         ill = i+idx_shift(1,iside); jll = j+idy_shift(1,iside) !(i,j) index of lower left  quadrant
                         ilc = i+idx_shift(6,iside); jlc = j+idy_shift(6,iside) !(i,j) index of lower center  quadrant
                         
                         scale = fvm(ie)%dp_fvm(i,j,ilev,n0_fvm)*fvm(ie)%area_sphere(i,j)*fraction
                         !
                         ! in all cases positive 1,1,1 flux
                         !
                         elem(ie)%sub_elem_mass_flux(i,j    ,iside,ilev)          = scale
                         elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p2(iside),ilev)=-scale
                         if (d2== 1) then
                            elem(ie)%sub_elem_mass_flux(i  ,  j,iside_p1(iside),ilev)= scale
                            ! set flux on other side!
                            elem(ie)%sub_elem_mass_flux(iur,jur,iside_p3(iside),ilev)=-scale
                         else
                            elem(ie)%sub_elem_mass_flux(i  ,  j,iside_p1(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(iur,jur,iside_p3(iside),ilev)= 0.0D0
                         end if
                         if (d2==-1) then
                            elem(ie)%sub_elem_mass_flux(i  ,  j,iside_p1(iside),ilev)=-scale
                            elem(ie)%sub_elem_mass_flux(iur,jur,iside_p3(iside),ilev)= scale
                         else
                            elem(ie)%sub_elem_mass_flux(i  ,  j,iside_p1(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(iur,jur,iside_p3(iside),ilev)= 0.0D0
                         end if
                         if (d3== 1) then
                            elem(ie)%sub_elem_mass_flux(i, j   ,iside_p3(iside),ilev)= scale
                            elem(ie)%sub_elem_mass_flux(iul,jul,iside_p1(iside),ilev)=-scale
                         else
                            elem(ie)%sub_elem_mass_flux(i, j   ,iside_p3(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(iul,jul,iside_p1(iside),ilev)= 0.0D0
                         end if
                         if (d3==-1) then
                            elem(ie)%sub_elem_mass_flux(i, j   ,iside_p3(iside),ilev)=-scale
                            elem(ie)%sub_elem_mass_flux(iul,jul,iside_p1(iside),ilev)= scale
                         else
                            elem(ie)%sub_elem_mass_flux(i, j   ,iside_p3(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(iul,jul,iside_p1(iside),ilev)= 0.0D0
                         end if
                         if (d4==1 ) then
                            elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p2(iside),ilev)=-scale
                            elem(ie)%sub_elem_mass_flux(iur,jur,iside          ,ilev)= scale
                         else
                            elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p2(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(iur,jur,iside          ,ilev)= 0.0D0
                         end if
                         if (d4==-1) then
                            elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p2(iside),ilev)= scale
                            elem(ie)%sub_elem_mass_flux(iur,jur,iside          ,ilev)=-scale
                         else
                            elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p2(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(iur,jur,iside          ,ilev)= 0.0D0
                         end if
                         if (d5==1 ) then
                            elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p3(iside),ilev)= scale
                            elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p1(iside),ilev)=-scale
                         else
                            elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p3(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p1(iside),ilev)= 0.0D0
                         end if
                         if (d5==-1) then
                            elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p3(iside),ilev)=-scale
                            elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p1(iside),ilev)= scale
                         else
                            elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p3(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p1(iside),ilev)= 0.0D0
                         end if
                         if (d6==1)  then
                            elem(ie)%sub_elem_mass_flux(iul,jul,iside          ,ilev)= scale
                            elem(ie)%sub_elem_mass_flux(ill,jll,iside_p2(iside),ilev)=-scale                
                         else
                            elem(ie)%sub_elem_mass_flux(iul,jul,iside          ,ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(ill,jll,iside_p2(iside),ilev)= 0.0D0                
                         end if
                         if (d6==-1) then
                            elem(ie)%sub_elem_mass_flux(iul,jul,iside          ,ilev)=-scale
                            elem(ie)%sub_elem_mass_flux(ill,jll,iside_p2(iside),ilev)= scale                
                         else
                            elem(ie)%sub_elem_mass_flux(iul,jul,iside          ,ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(ill,jll,iside_p2(iside),ilev)= 0.0D0                
                         end if
                         if (d7== 1) then
                            elem(ie)%sub_elem_mass_flux(ill,jll,iside_p1(iside),ilev)= scale
                            elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p3(iside),ilev)=-scale
                         else
                            elem(ie)%sub_elem_mass_flux(ill,jll,iside_p1(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p3(iside),ilev)= 0.0D0
                         end if
                         if (d7==-1) then
                            elem(ie)%sub_elem_mass_flux(ill,jll,iside_p1(iside),ilev)=-scale
                            elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p3(iside),ilev)= scale
                         else
                            elem(ie)%sub_elem_mass_flux(ill,jll,iside_p1(iside),ilev)= 0.0D0
                            elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p3(iside),ilev)= 0.0D0
                         end if

!                         if (d2== 1) elem(ie)%sub_elem_mass_flux(i  ,  j,iside_p1(iside),ilev)=scale
!                         if (d2==-1) elem(ie)%sub_elem_mass_flux(iur,jur,iside_p3(iside),ilev)=scale
!                         if (d3== 1) elem(ie)%sub_elem_mass_flux(i, j   ,iside_p3(iside),ilev)=scale
!                         if (d3==-1) elem(ie)%sub_elem_mass_flux(iul,jul,iside_p1(iside),ilev)=scale
!                         if (d4==1 ) elem(ie)%sub_elem_mass_flux(iur,jur,iside          ,ilev)=scale
!                         if (d4==-1) elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p2(iside),ilev)=scale
!                         if (d5==1 ) elem(ie)%sub_elem_mass_flux(ilr,jlr,iside_p3(iside),ilev)=scale
!                         if (d5==-1) elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p1(iside),ilev)=scale
!                         if (d6==1)  elem(ie)%sub_elem_mass_flux(iul,jul,iside          ,ilev)=scale
!                         if (d6==-1) elem(ie)%sub_elem_mass_flux(ill,jll,iside_p2(iside),ilev)=scale                
!                         if (d7== 1) elem(ie)%sub_elem_mass_flux(ill,jll,iside_p1(iside),ilev)=scale
!                         if (d7==-1) elem(ie)%sub_elem_mass_flux(ilc,jlc,iside_p3(iside),ilev)=scale

                         ie=ie+1
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do
 else

 ! edge case
 
 ilev=1;
 i=2; j=1;
 scale = fvm(2)%dp_fvm(i,j,ilev,n0_fvm)*fvm(2)%area_sphere(i,j)*fraction
 !
 !
 !                             i,j,iside
! ie1=2
! ie2=1112
! elem(ie1)%sub_elem_mass_flux(2,1,1,ilev) =  scale
! elem(ie2)%sub_elem_mass_flux(2,3,3,ilev) = -scale
!
! elem(ie2)%sub_elem_mass_flux(1,3,2,ilev) =  scale
! elem(ie2)%sub_elem_mass_flux(2,3,4,ilev) = -scale
!
! elem(ie2)%sub_elem_mass_flux(1,3,3,ilev) =  scale
! elem(ie1)%sub_elem_mass_flux(1,1,1,ilev) = -scale
!
! elem(ie1)%sub_elem_mass_flux(1,1,2,ilev) =  scale
! elem(ie1)%sub_elem_mass_flux(2,1,4,ilev) = -scale

 !
 ! corner case 1
 !
! ie=211
! ie1=900
! ie2=1126
! elem(ie2)%sub_elem_mass_flux(1,1,1,ilev) =  scale
! elem(ie )%sub_elem_mass_flux(1,3,3,ilev) = -scale
!
! elem(ie1)%sub_elem_mass_flux(3,3,2,ilev) =  scale
! elem(ie )%sub_elem_mass_flux(1,3,4,ilev) = -scale
!
! elem(ie1)%sub_elem_mass_flux(3,3,3,ilev) =  scale
! elem(ie2)%sub_elem_mass_flux(1,1,4,ilev) = -scale

 !
 ! corner case 1
 !
! ie=211
! ie1=900
! ie2=1126
! elem(ie2)%sub_elem_mass_flux(1,1,1,ilev) =  scale
! elem(ie )%sub_elem_mass_flux(1,3,3,ilev) = -scale
!
! elem(ie2)%sub_elem_mass_flux(1,1,4,ilev) =  scale
! elem(ie1)%sub_elem_mass_flux(3,3,3,ilev) = -scale
!
! elem(ie1)%sub_elem_mass_flux(3,3,2,ilev) =  scale
! elem(ie )%sub_elem_mass_flux(1,3,4,ilev) = -scale!!
!
! elem(ie1)%sub_elem_mass_flux(1,1,2,ilev) =  scale
! elem(ie1)%sub_elem_mass_flux(2,1,4,ilev) = -scale
!!
!
! elem(1)%sub_elem_mass_flux(2,2,1,ilev) =   scale
! elem(1)%sub_elem_mass_flux(2,1,3,ilev) =  -scale!
!
! elem(3)%sub_elem_mass_flux(2,2,1,ilev) =  -scale
! elem(3)%sub_elem_mass_flux(2,1,3,ilev) =   scale
!
! elem(2)%sub_elem_mass_flux(2,2,2,ilev) =   scale
! elem(2)%sub_elem_mass_flux(3,2,4,ilev) =  -scale
!

 !
 ! ill-conditioned case
 !
 !dp_area= 1.0000009182857403D0
! elem(266)%sub_elem_mass_flux(2,1,1,30) =   -5.075291944369226033412592115346D-7
! elem(251)%sub_elem_mass_flux(2,3,3,30) =    5.075291944369226033412592115346D-07

    end if
  end subroutine set_flow_cases
end module fvm_consistent_se_cslam
