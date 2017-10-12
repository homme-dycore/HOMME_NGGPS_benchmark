#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef _COLLAPSE_AND_ALIGN
#define OMP_COLLAPSE_SIMD  $OMP SIMD COLLAPSE(2)
#define DIR_VECTOR_ALIGNED DIR$ VECTOR ALIGNED
#endif

#define _B4B 1

!#define _DBG_ print *,"File:",__FILE__," at ",__LINE__
!#define _DBG_ !DBG
!
!
module prim_advance_mod
  use edgetype_mod, only : EdgeDescriptor_t, EdgeBuffer_t
  use kinds, only : real_kind, iulog
  use parallel_mod, only : abortmp, parallel_t, iam, boundaryCommMethod
  use control_mod, only : se_prescribed_wind_2d
  use thread_mod , only : horz_num_threads, vert_num_threads, tracer_num_threads, omp_get_thread_num
  use perf_mod, only: t_startf, t_stopf, t_barrierf ! EXTERNAL
  use perf_utils, only: t_detail_low, t_detail_medium, t_detail_high, t_detail_max  ! EXTERNAL

  implicit none
  private
  save
  public :: prim_advance_exp, prim_advance_si, prim_advance_init, preq_robert3,&
       applyCAMforcing, smooth_phis, overwrite_SEdensity

#ifdef CAM
  public ::  calc_tot_energy_dynamics
#endif

#ifdef TRILINOS
  public :: distribute_flux_at_corners
#endif
  
  type (EdgeBuffer_t) :: edge1
  type (EdgeBuffer_t) :: edge2
  type (EdgeBuffer_t),public :: edge3p1
  type (EdgeBuffer_t),public :: edge3p2

  real (kind=real_kind) :: initialized_for_dt   = 0

  real (kind=real_kind), allocatable :: ur_weights(:)

contains

  subroutine prim_advance_init(par, elem,integration)
    use edge_mod, only : initEdgeBuffer
    use element_mod, only : element_t
    use dimensions_mod, only : nlev, nelemd
    use control_mod, only : qsplit,rsplit
    implicit none
    
    type (parallel_t) :: par
    type (element_t), intent(inout), target   :: elem(:)
    character(len=*)    , intent(in) :: integration
    integer :: i
    integer :: ie
!    integer, allocatable :: globalid(:)

!    print *,'prim_advance_init: nelemd:= ',nelemd
!    allocate(globalid(nelemd))x
!    print *,'prim_advance_init: after allocate '
!    print *,'prim_advance_init: before call to initNewEdgeBuffer rsplit: ',rsplit
    if (rsplit==0) then
       call initEdgeBuffer(par,edge3p1,elem,3*nlev+1,bndry_type=boundaryCommMethod, nthreads=horz_num_threads)
       call initEdgeBuffer(par,edge3p2,elem,3*nlev+1,bndry_type=boundaryCommMethod, nthreads=horz_num_threads*vert_num_threads)
    else
       ! need extra buffer space for dp3d
       call initEdgeBuffer(par,edge3p1,elem,4*nlev+1,bndry_type=boundaryCommMethod, nthreads=horz_num_threads)
       call initEdgeBuffer(par,edge3p2,elem,4*nlev+1,bndry_type=boundaryCommMethod, nthreads=horz_num_threads*vert_num_threads)
    endif

    if(integration == 'semi_imp') then
       call initEdgeBuffer(par,edge1,elem,nlev)
       call initEdgeBuffer(par,edge2,elem,2*nlev)
    end if

    ! compute averaging weights for RK+LF (tstep_type=1) timestepping:
    allocate(ur_weights(qsplit))
    ur_weights(:)=0.0d0

    if(mod(qsplit,2).NE.0)then
       ur_weights(1)=1.0d0/qsplit
       do i=3,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    else
       do i=2,qsplit,2
         ur_weights(i)=2.0d0/qsplit
       enddo
    endif

  end subroutine prim_advance_init


  subroutine prim_advance_exp(elem, deriv, hvcoord, hybrid,dt, tl,  nets, nete)

    use bndry_mod,         only: bndry_exchangev
    use control_mod,       only: prescribed_wind, tstep_type, rsplit, qsplit, moisture, test_case
#ifndef CAM
    use dcmip_wrapper_mod, only: set_dcmip_1_1_fields, set_dcmip_1_2_fields
#endif
    use derivative_mod,    only: derivative_t
    use derivative_mod,    only: vorticity_routine
    use dimensions_mod,    only: np, nlev, nlevp, nvar, nc, nelemd, ldry_mass_vertical_coordinates
    use edge_mod,          only: edgevpack, edgevunpack, initEdgeBuffer
    use edgetype_mod,      only: EdgeBuffer_t
    use element_mod,       only: element_t
    use hybvcoord_mod,     only: hvcoord_t
    use hybrid_mod,        only: hybrid_t, PrintHybrid, config_thread_region
    use reduction_mod,     only: reductionbuffer_ordered_1d_t
    use time_mod,          only: TimeLevel_t,  timelevel_qdp, tevolve
    use diffusion_mod,     only:  prim_diffusion

#ifdef TRILINOS
    use prim_derived_type_mod ,only : derived_type, initialize
    use, intrinsic :: iso_c_binding
#endif

#ifndef CAM
    use asp_tests, only : asp_advection_vertical
#else
    use control_mod, only : prescribed_vertwind
#endif

    implicit none

    type (element_t), intent(inout), target   :: elem(:)
    type (derivative_t)  , intent(in) :: deriv
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t)      , intent(in) :: hybrid
    real (kind=real_kind), intent(in) :: dt
    type (TimeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

    ! =================
    ! Local
    ! =================
    type (hybrid_t)       :: hybridnew
    real (kind=real_kind) ::  dt2, time, dt_vis, x, eta_ave_w
    real (kind=real_kind) ::  eta_dot_dpdn(np,np,nlevp)
    real (kind=real_kind) ::  dp(np,np)
    real (kind=real_kind) ::  tempdp3d(np,np)
    real (kind=real_kind) ::  tempmass(nc,nc)
    real (kind=real_kind) ::  tempflux(nc,nc,4)
    real (kind=real_kind) ::  deta
    integer :: ie,nm1,n0,np1,nstep,method,qsplit_stage,k, qn0
    integer :: n,i,j,lx,lenx,ithr
    integer :: region_num_threads

#ifdef TRILINOS
    real (c_double) ,allocatable, dimension(:) :: xstate(:)

!    type (element_t) :: pc_elem(size(elem))
!    type (element_t) :: jac_elem(size(elem))

! state_object is a derived data type passed thru noxinit as a pointer
    type(derived_type) ,target         :: state_object
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    type(derived_type) ,target         :: pre_object
    type(derived_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre
    type(derived_type) ,target         :: jac_object
    type(derived_type) ,pointer         :: jptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_jac

    integer(c_int) :: ierr = 0

  interface

   subroutine noxsolve(vectorSize,vector,v_container,p_container,j_container,ierr) &
!   subroutine noxsolve(vectorSize,vector,v_container) &
     bind(C,name='noxsolve')
    use ,intrinsic :: iso_c_binding
      integer(c_int)                :: vectorSize
      real(c_double)  ,dimension(*) :: vector
      type(c_ptr)                   :: v_container
      type(c_ptr)                   :: p_container  !precon ptr
      type(c_ptr)                   :: j_container  !analytic jacobian ptr
      integer(c_int)                :: ierr         !error flag
    end subroutine noxsolve

  end interface

#endif

    call t_startf('prim_advance_exp', t_detail_low)
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    ! timelevel to use for accessing Qdp() to compute virtual temperature
    qn0 = -1    ! -1 = disabled (assume dry dynamics)
    if ( moisture /= "dry") then
       call TimeLevel_Qdp(tl, qsplit, qn0)  ! compute current Qdp() timelevel
    endif

    ! integration = "explicit"
    !
    !   tstep_type=0  pure leapfrog except for very first timestep   CFL=1
    !                    typically requires qsplit=4 or 5
    !   tstep_type=1  RK2 followed by qsplit-1 leapfrog steps        CFL=close to qsplit
    !                    typically requires qsplit=4 or 5
    !   tstep_type=2  RK2-SSP 3 stage (as used by tracers)           CFL=.58
    !                    optimal in terms of SSP CFL, but not        CFLSSP=2
    !                    optimal in terms of CFL
    !                    typically requires qsplit=3
    !                    but if windspeed > 340m/s, could use this
    !                    with qsplit=1
    !   tstep_type=3  classic RK3                                    CFL=1.73 (sqrt(3))
    !
    !   tstep_type=4  Kinnmark&Gray RK4 4 stage                      CFL=sqrt(8)=2.8
    !                 should we replace by standard RK4 (CFL=sqrt(8))?
    !                 (K&G 1st order method has CFL=3)
    !   tstep_type=5  Kinnmark&Gray RK3 5 stage 3rd order            CFL=3.87  (sqrt(15))
    !                 From Paul Ullrich.  3rd order for nonlinear terms also
    !                 K&G method is only 3rd order for linear
    !                 optimal: for windspeeds ~120m/s,gravity: 340m/2
    !                 run with qsplit=1
    !                 (K&G 2nd order method has CFL=4. tiny CFL improvement not worth 2nd order)
    !
    ! integration = "full_imp"
    !
    !   tstep_type=11  Backward Euler implicit dynamics, first order
    !   tstep_type=12  BDF2 implicit dynamics, second order
    !

    ! default weights for computing mean dynamics fluxes
    if ((tstep_type==11).or.(tstep_type==12)) then
     eta_ave_w = 1d0  ! don't use eta_ave_w averaging for implicit
    else 
     eta_ave_w = 1d0/qsplit
    end if 

    if(tstep_type==0)then
       method=0                ! pure leapfrog
       if (nstep==0) method=1  ! but use RK2 on first step
    else if (tstep_type==1) then
       method=0                           ! LF
       qsplit_stage = mod(nstep,qsplit)
       if (qsplit_stage==0) method=1      ! RK2 on first of qsplit steps
       ! RK2 + LF scheme has tricky weights:
       eta_ave_w=ur_weights(qsplit_stage+1)
    else
       method = tstep_type                ! other RK variants
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    ! fix dynamical variables, skip dynamics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

    if (1==prescribed_wind .and. .not.se_prescribed_wind_2d) then
       time=tl%nstep*dt

       do ie=nets,nete
#ifdef CAM

          ! Use CAM prescribed eta_dot to get vertical flux

          if (prescribed_vertwind==1) then
             eta_dot_dpdn(:,:,1)      = 0.0d0
             eta_dot_dpdn(:,:,nlev+1) = 0.0d0
             do k = 2,nlev
                dp(:,:) =&
                     (hvcoord%hyam(k) - hvcoord%hyam(k-1))*hvcoord%ps0 + &
                     (hvcoord%hybm(k) - hvcoord%hybm(k-1))*elem(ie)%state%ps(:,:,tl%n0)
                deta = hvcoord%etam(k)-hvcoord%etam(k-1)
                eta_dot_dpdn(:,:,k) = dp(:,:)*elem(ie)%derived%etadot_prescribed(:,:,k)/deta
             end do
          endif
#else

        ! Use HOMME ASP tests to set prescribed fields
        if(test_case(1:4)=="asp_") then
          elem(ie)%state%ps(:,:,np1)  = elem(ie)%state%ps(:,:,n0)
          elem(ie)%state%lnps(:,:,np1)  = elem(ie)%state%lnps(:,:,n0)
          elem(ie)%state%T(:,:,:,np1)   = elem(ie)%state%T(:,:,:,n0)
          elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,n0)
          elem(ie)%derived%div = 0
          call asp_advection_vertical(time,hvcoord,eta_dot_dpdn)

        else if(test_case(1:8)=="DCMIP1-1") then
          call set_dcmip_1_1_fields(elem(ie), eta_dot_dpdn, hybrid,hvcoord,np1,time)

        else if(test_case(1:8)=="DCMIP1-2") then
          call set_dcmip_1_2_fields(elem(ie), eta_dot_dpdn, hybrid,hvcoord,np1,time)

        endif

#endif
          ! Adjust vertical fluxes for qsplit and rsplit values

          if (rsplit==0) then ! eulerian case
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
            do k=1,nlev
              do j=1,np
                do i=1,np
                  elem(ie)%derived%eta_dot_dpdn(i,j,k) = elem(ie)%derived%eta_dot_dpdn(i,j,k) + &
                                                                  eta_dot_dpdn(i,j,k) * eta_ave_w 
                end do
              end do
            end do
          else
             ! lagrangian case
             elem(ie)%derived%eta_dot_dpdn = 0

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
             do k=1,nlev
               do j=1,np
!dir$ ivdep
                 do i=1,np
                   elem(ie)%state%dp3d(i,j,k,np1) = elem(ie)%state%dp3d(i,j,k,n0)  &
                                  + dt*(eta_dot_dpdn(i,j,k+1) - eta_dot_dpdn(i,j,k))
                 enddo
               enddo
             enddo
          end if
        end do



        do ie=nets,nete
          ! subcycling code uses a mean flux to advect tracers 
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(dp)
#endif
          do k=1,nlev
             if (rsplit==0) then
                dp(:,:) =&
                     ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                     ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(:,:,tl%n0) 
             else
                dp(:,:) = elem(ie)%state%dp3d(:,:,k,tl%n0)
             end if

             elem(ie)%derived%vn0(:,:,1,k)=elem(ie)%derived%vn0(:,:,1,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,1,k,n0)*dp(:,:)
             elem(ie)%derived%vn0(:,:,2,k)=elem(ie)%derived%vn0(:,:,2,k)+&
                  eta_ave_w*elem(ie)%state%v(:,:,2,k,n0)*dp(:,:)
          enddo
       end do
       call t_stopf('prim_advance_exp', t_detail_low)
!      call t_adj_detailf(-1)
       return
    endif


    ! ==================================
    ! Take timestep
    ! ==================================

    dt_vis = dt
    if (method==0) then
       ! regular LF step
       dt2 = 2*dt
       if (.not.ldry_mass_vertical_coordinates) then
          call compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       else
          call compute_and_apply_rhs_dry(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       end if
       dt_vis = dt2  ! dt to use for time-split dissipation
    else if (method==1) then
       if (.not.ldry_mass_vertical_coordinates) then
          ! RK2
          ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
          call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       else
          call compute_and_apply_rhs_dry(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       end if
    else if (method==2) then
       if (.not.ldry_mass_vertical_coordinates) then
          ! RK2-SSP 3 stage.  matches tracer scheme. optimal SSP CFL, but
          ! not optimal for regular CFL
          ! u1 = u0 + dt/2 RHS(u0)
          call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w/3)
          ! u2 = u1 + dt/2 RHS(u1)
          call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w/3)
          ! u3 = u2 + dt/2 RHS(u2)
          call compute_and_apply_rhs(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w/3)
       else
          call compute_and_apply_rhs_dry(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w/3)
          call compute_and_apply_rhs_dry(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w/3)
          call compute_and_apply_rhs_dry(np1,np1,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w/3)

       end if

       ! unew = u/3 +2*u3/3  = u + 1/3 (RHS(u) + RHS(u1) + RHS(u2))
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,np1)= elem(ie)%state%v(:,:,:,:,n0)/3 &
               + 2*elem(ie)%state%v(:,:,:,:,np1)/3
          elem(ie)%state%T(:,:,:,np1)= elem(ie)%state%T(:,:,:,n0)/3 &
               + 2*elem(ie)%state%T(:,:,:,np1)/3
          if (rsplit==0) then
             elem(ie)%state%ps(:,:,np1)= elem(ie)%state%ps(:,:,n0)/3 &
                  + 2*elem(ie)%state%ps(:,:,np1)/3
          else
             elem(ie)%state%dp3d(:,:,:,np1)= elem(ie)%state%dp3d(:,:,:,n0)/3 &
                  + 2*elem(ie)%state%dp3d(:,:,:,np1)/3
          endif
       enddo
    else if (method==3) then
       if (.not.ldry_mass_vertical_coordinates) then
          ! classic RK3  CFL=sqrt(3)
          ! u1 = u0 + dt/3 RHS(u0)
          call compute_and_apply_rhs(np1,n0,n0,qn0,dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u2 = u0 + dt/2 RHS(u1)
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u3 = u0 + dt RHS(u2)
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       else
          call compute_and_apply_rhs_dry(np1,n0,n0,qn0,dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       end if
    else if (method==4) then
       if (.not.ldry_mass_vertical_coordinates) then
          ! KG 4th order 4 stage:   CFL=sqrt(8)
          ! low storage version of classic RK4
          ! u1 = u0 + dt/4 RHS(u0)
          call compute_and_apply_rhs(np1,n0,n0,qn0,dt/4,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u2 = u0 + dt/3 RHS(u1)
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u3 = u0 + dt/2 RHS(u2)
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u4 = u0 + dt RHS(u3)
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       else
          call compute_and_apply_rhs_dry(np1,n0,n0,qn0,dt/4,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u2 = u0 + dt/3 RHS(u1)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u3 = u0 + dt/2 RHS(u2)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u4 = u0 + dt RHS(u3)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       end if
    else if (method==5) then
#if 0
       if (.not.ldry_mass_vertical_coordinates) then
          ! KG 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
          ! but nonlinearly only 2nd order
          ! u1 = u0 + dt/5 RHS(u0)
          call compute_and_apply_rhs(np1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u2 = u0 + dt/5 RHS(u1)
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt/5,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u3 = u0 + dt/3 RHS(u2)
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u4 = u0 + dt/2 RHS(u3)
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u5 = u0 + dt RHS(u4)
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       else
          call compute_and_apply_rhs_dry(np1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt/5,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt/2,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w)
       end if
#else
       if (.not.ldry_mass_vertical_coordinates) then
          ! Ullrich 3nd order 5 stage:   CFL=sqrt( 4^2 -1) = 3.87
          ! u1 = u0 + dt/5 RHS(u0)  (save u1 in timelevel nm1)
          
          !
          ! phl: rhs: t=t
          !
          call compute_and_apply_rhs(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w/4)
          ! u2 = u0 + dt/5 RHS(u1)
          
          !
          ! phl: rhs: t=t+dt/5
          !
          call compute_and_apply_rhs(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u3 = u0 + dt/3 RHS(u2)
          !
          ! phl: rhs: t=t+2*dt/5
          !
          
          call compute_and_apply_rhs(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          ! u4 = u0 + 2dt/3 RHS(u3)
          
          !
          ! phl: rhs: t=t+2*dt/5+dt/3
          !
          
          call compute_and_apply_rhs(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
       else
          call compute_and_apply_rhs_dry(nm1,n0,n0,qn0,dt/5,elem,hvcoord,hybrid,&
               deriv,nets,nete,eta_ave_w/4)
          call compute_and_apply_rhs_dry(np1,n0,nm1,qn0,dt/5,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
          call compute_and_apply_rhs_dry(np1,n0,np1,qn0,2*dt/3,elem,hvcoord,hybrid,&
               deriv,nets,nete,0d0)
       end if
          ! compute (5*u1/4 - u0/4) in timelevel nm1:
       do ie=nets,nete
          elem(ie)%state%v(:,:,:,:,nm1)= (5*elem(ie)%state%v(:,:,:,:,nm1) &
               - elem(ie)%state%v(:,:,:,:,n0) ) /4
          elem(ie)%state%T(:,:,:,nm1)= (5*elem(ie)%state%T(:,:,:,nm1) &
               - elem(ie)%state%T(:,:,:,n0) )/4
          if (rsplit==0) then
             elem(ie)%state%ps(:,:,nm1)= ( 5*elem(ie)%state%ps(:,:,nm1) &
                  - elem(ie)%state%ps(:,:,n0) )/4
          else
             elem(ie)%state%dp3d(:,:,:,nm1)= (5*elem(ie)%state%dp3d(:,:,:,nm1) &
                  - elem(ie)%state%dp3d(:,:,:,n0) )/4
          endif
       enddo
       ! u5 = (5*u1/4 - u0/4) + 3dt/4 RHS(u4)
       !
       ! phl: rhs: t=t+2*dt/5+dt/3+3*dt/4         -wrong RK times ...
       !
       if (.not.ldry_mass_vertical_coordinates) then
          call compute_and_apply_rhs(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
               deriv,nets,nete,3*eta_ave_w/4)
       else
          call compute_and_apply_rhs_dry(np1,nm1,np1,qn0,3*dt/4,elem,hvcoord,hybrid,&
               deriv,nets,nete,3*eta_ave_w/4)
       end if
       ! final method is the same as:
       ! u5 = u0 +  dt/4 RHS(u0)) + 3dt/4 RHS(u4)
#endif

    else if ((method==11).or.(method==12)) then
       ! Fully implicit JFNK method (vertically langragian not active yet)
       if (rsplit > 0) then
       call abortmp('ERROR: full_imp integration not yet coded for vert lagrangian adv option')
       end if
!      if (hybrid%masterthread) print*, "fully implicit integration is still under development"

#ifdef TRILINOS
      lenx=(np*np*nlev*3 + np*np*1)*(nete-nets+1)  ! 3 3d vars plus 1 2d vars
      allocate(xstate(lenx))
      xstate(:) = 0d0

      call initialize(state_object, method, elem, hvcoord, &
        qn0, hybrid, deriv, dt, tl, nets, nete)

      call initialize(pre_object, method, elem, hvcoord, &
        qn0, hybrid, deriv, dt, tl, nets, nete)

      call initialize(jac_object, method, elem, hvcoord, &
        qn0, hybrid, deriv, dt, tl, nets, nete)

!      pc_elem = elem
!      jac_elem = elem

        fptr => state_object
        c_ptr_to_object =  c_loc(fptr)
        pptr => state_object
        c_ptr_to_pre =  c_loc(pptr)
        jptr => state_object
        c_ptr_to_jac =  c_loc(jptr)

! create flat state vector to pass through NOX
! use previous time step as the first guess for the new one (because with LF time level update n0=np1)

       np1 = n0

       lx = 1
       do ie=nets,nete
         do k=1,nlev
           do j=1,np
             do i=1,np
               xstate(lx) = elem(ie)%state%v(i,j,1,k,n0)
               lx = lx+1
             end do
           end do
         end do
       end do
       do ie=nets,nete
         do k=1,nlev
           do j=1,np
             do i=1,np
               xstate(lx) = elem(ie)%state%v(i,j,2,k,n0)
               lx = lx+1
             end do
           end do
         end do
       end do
       do ie=nets,nete
         do k=1,nlev
           do j=1,np
             do i=1,np
               xstate(lx) = elem(ie)%state%T(i,j,k,n0)
               lx = lx+1
             end do
           end do
         end do
       end do
       do ie=nets,nete
         do j=1,np
           do i=1,np
             xstate(lx) = elem(ie)%state%ps(i,j,n0)
            lx = lx+1
         end do
       end do
     end do
       
! activate these lines to test infrastructure and still solve with explicit code
!       ! RK2
!       ! forward euler to u(dt/2) = u(0) + (dt/2) RHS(0)  (store in u(np1))
!       call compute_and_apply_rhs(np1,n0,n0,qn0,dt/2,elem,hvcoord,hybrid,&
!            deriv,nets,nete,0d0)
!       ! leapfrog:  u(dt) = u(0) + dt RHS(dt/2)     (store in u(np1))
!       call compute_and_apply_rhs(np1,n0,np1,qn0,dt,elem,hvcoord,hybrid,&
!            deriv,nets,nete,.false.,eta_ave_w)

! interface to use nox and loca solver libraries using JFNK, and returns xstate(n+1)
    call noxsolve(size(xstate), xstate, c_ptr_to_object, c_ptr_to_pre, c_ptr_to_jac, ierr)

    if (ierr /= 0) call abortmp('Error in noxsolve: Newton failed to converge')

      call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
      elem = fptr%base

          lx = 1
          do ie=nets,nete
            do k=1,nlev
              do j=1,np
                do i=1,np
                  elem(ie)%state%v(i,j,1,k,np1) = xstate(lx)
                  lx = lx+1
                end do
              end do
            end do
          end do
          do ie=nets,nete
            do k=1,nlev
              do j=1,np
                do i=1,np
                  elem(ie)%state%v(i,j,2,k,np1) = xstate(lx) 
                  lx = lx+1
                end do
              end do
            end do
          end do
          do ie=nets,nete
            do k=1,nlev
              do j=1,np
                do i=1,np
                  elem(ie)%state%T(i,j,k,np1) = xstate(lx)
                  lx = lx+1
                end do
              end do
            end do
          end do
          do ie=nets,nete
            do j=1,np
              do i=1,np
                elem(ie)%state%ps(i,j,np1) = xstate(lx)
                lx = lx+1
              end do
            end do
          end do

#endif

    else
       call abortmp('ERROR: bad choice of tstep_type')
    endif

    ! ==============================================
    ! Time-split Horizontal diffusion: nu.del^2 or nu.del^4
    ! U(*) = U(t+1)  + dt2 * HYPER_DIFF_TERM(t+1)
    ! ==============================================

    call t_startf('advance_hypervis', t_detail_low)

    ! note:time step computes u(t+1)= u(t*) + RHS.
    ! for consistency, dt_vis = t-1 - t*, so this is timestep method dependent
    if (tstep_type==0) then
       ! leapfrog special case
       call advance_hypervis_lf(edge3p2,elem,hvcoord,hybrid,deriv,nm1,n0,np1,nets,nete,dt_vis)
    else if (method<=10) then ! not implicit
       if (rsplit==0) then
          ! forward-in-time, maybe hypervis applied to PS
          call advance_hypervis(edge3p2,elem,hvcoord,hybrid,deriv,np1,nets,nete,dt_vis,eta_ave_w)
       else
          ! forward-in-time, hypervis applied to dp3d
         if ( vert_num_threads > 1 ) then
!$OMP PARALLEL NUM_THREADS(vert_num_threads), DEFAULT(SHARED), PRIVATE(hybridnew)
           hybridnew = config_thread_region(hybrid,'vertical')
           call advance_hypervis_dp(edge3p2,elem,hvcoord,hybridnew,deriv,np1,qn0,nets,nete,dt_vis,eta_ave_w)
!$OMP END PARALLEL
         else
           call advance_hypervis_dp(edge3p1,elem,hvcoord,hybrid,deriv,np1,qn0,nets,nete,dt_vis,eta_ave_w)
         endif
       endif
    endif
    call t_stopf('advance_hypervis', t_detail_low)

    tevolve=tevolve+dt

    call t_stopf('prim_advance_exp', t_detail_low)
!    call t_adj_detailf(-1)
  end subroutine prim_advance_exp


subroutine prim_advance_si(elem, nets, nete, cg, blkjac, red, &
          refstate, hvcoord, deriv, flt, hybrid, tl, dt)
       use bndry_mod, only : bndry_exchangev
       use cg_mod, only : cg_t, cg_create
       use control_mod, only : filter_freq,debug_level, precon_method
       use derivative_mod, only : derivative_t
       use derivative_mod, only: vorticity_routine
       use derivative_mod, only: divergence_routine
       use derivative_mod, only: gradient_routine
       use derivative_mod, only: gradient_wk_routine
       use dimensions_mod, only: np, nlev, nlevp
       use edge_mod, only : edgevpack, edgevunpack, initEdgeBuffer
       use edgetype_mod, only : EdgeBuffer_t
       use element_mod, only : element_t
       use filter_mod, only : filter_t, preq_filter
       use hybvcoord_mod, only : hvcoord_t
       use hybrid_mod, only : hybrid_t
       use prim_si_ref_mod, only : ref_state_t, set_vert_struct_mat
       use reduction_mod, only : reductionbuffer_ordered_1d_t
       use solver_mod, only : pcg_solver, blkjac_t, blkjac_init
       use time_mod, only : TimeLevel_t
       use prim_si_mod, only : preq_vertadv, preq_omegap, preq_pressure
       use diffusion_mod, only :  prim_diffusion
       use physical_constants, only : kappa, rrearth, rgas, cp, rwater_vapor
       use physics_mod, only : virtual_temperature, virtual_specific_heat
       implicit none

       integer, intent(in)               :: nets,nete
       type (element_t), intent(inout), target :: elem(:)
       type (blkjac_t), allocatable      :: blkjac(:)

       type (cg_t)                       :: cg

       type (ReductionBuffer_ordered_1d_t), intent(inout) :: red

       type (ref_state_t), intent(in), target :: refstate
       type (hvcoord_t), intent(in)      :: hvcoord
       type (derivative_t), intent(in)   :: deriv
       type (filter_t), intent(in)       :: flt
       type (hybrid_t), intent(in)       :: hybrid
       type (TimeLevel_t), intent(in)    :: tl
       real(kind=real_kind), intent(in)  :: dt
       real(kind=real_kind)              :: time_adv
#ifndef _CRAYFTN
       ! ==========================
       ! Local variables...
       ! ==========================

       real(kind=real_kind)                           :: ps0
       real(kind=real_kind)                           :: psref

       real(kind=real_kind), dimension(np,np)         :: ps
       real(kind=real_kind), dimension(np,np)         :: rps
       real(kind=real_kind), dimension(np,np,nlev)    :: rpmid
       real(kind=real_kind), dimension(np,np,nlev)    :: omegap
       real(kind=real_kind), dimension(np,np,nlev)    :: rpdel

       real(kind=real_kind) :: pintref(nlevp)
       real(kind=real_kind) :: pdelref(nlev)
       real(kind=real_kind) :: pmidref(nlev)
       real(kind=real_kind) :: rpdelref(nlev)
       real(kind=real_kind) :: rpmidref(nlev)

       real(kind=real_kind) :: pint(np,np,nlevp)
       real(kind=real_kind) :: pdel(np,np,nlev)
       real(kind=real_kind) :: pmid(np,np,nlev)

       real(kind=real_kind), dimension(np,np,nlevp) :: eta_dot_dp_deta
       real(kind=real_kind), dimension(np,np,nlev)  :: vgrad_ps

       real(kind=real_kind), dimension(np,np,nlev)   :: T_vadv
       real(kind=real_kind), dimension(np,np,2,nlev) :: v_vadv

       real(kind=real_kind), dimension(np,np)      :: HT
       real(kind=real_kind), dimension(np,np)      :: HrefT
       real(kind=real_kind), dimension(np,np)      :: HrefTm1

       real(kind=real_kind), dimension(np,np)      :: Gref0
       real(kind=real_kind), dimension(np,np)      :: Grefm1
       real(kind=real_kind), dimension(np,np)      :: E
       real(kind=real_kind), dimension(np,np)      :: Phi
       real(kind=real_kind), dimension(np,np)      :: dGref

       real(kind=real_kind), dimension(np,np,2)    :: vco
       real(kind=real_kind), dimension(np,np,2)    :: gradT
       real(kind=real_kind), dimension(np,np,2)    :: grad_Phi

       real(kind=real_kind), dimension(:,:), pointer  :: Emat
       real(kind=real_kind), dimension(:,:), pointer  :: Emat_inv
       real(kind=real_kind), dimension(:,:), pointer  :: Amat
       real(kind=real_kind), dimension(:,:), pointer  :: Amat_inv
       real(kind=real_kind), dimension(:), pointer    :: Lambda

       real(kind=real_kind), dimension(:), pointer    :: Tref
       real(kind=real_kind), dimension(:), pointer    :: RTref
       real(kind=real_kind), dimension(:), pointer    :: Pvec
       real(kind=real_kind), dimension(:,:), pointer  :: Href
       real(kind=real_kind), dimension(:,:), pointer  :: Tmat

       real(kind=real_kind) :: Vscript(np,np,2,nlev,nets:nete)
       real(kind=real_kind) :: Tscript(np,np,nlev,nets:nete)
       real(kind=real_kind) :: Pscript(np,np,nets:nete)
       real(kind=real_kind) :: Vtemp(np,np,2,nlev,nets:nete)

       real(kind=real_kind), dimension(np,np)      :: HrefTscript
       real(kind=real_kind), dimension(np,np)      :: suml
       real(kind=real_kind), dimension(np,np,2)    :: gVscript
       real(kind=real_kind), dimension(np,np,nlev) :: div_Vscript

       real(kind=real_kind) :: B(np,np,nlev,nets:nete)
       real(kind=real_kind) :: C(np,np,nlev,nets:nete)
       real(kind=real_kind) :: D(np,np,nlev,nets:nete)

       real(kind=real_kind) :: Gamma_ref(np,np,nlev,nets:nete)

       real(kind=real_kind) :: Gref(np,np,nlev,nets:nete)
       real(kind=real_kind) :: grad_dGref(np,np,2,nlev)
       real(kind=real_kind) :: grad_Gref(np,np,2,nlev)

       real(kind=real_kind) :: div(np,np)
       real(kind=real_kind) :: gv(np,np,2)

       real(kind=real_kind) :: dt2
       real(kind=real_kind) :: rpsref
       real(kind=real_kind) :: rdt
       real(kind=real_kind) :: hkk, hkl
       real(kind=real_kind) :: ddiv

       real(kind=real_kind) :: vgradT
       real(kind=real_kind) :: hybfac
       real(kind=real_kind) :: Crkk
       real(kind=real_kind) :: v1,v2
       real(kind=real_kind) :: term

       real(kind=real_kind) :: Vs1,Vs2
       real(kind=real_kind) :: glnps1, glnps2
       real(kind=real_kind) :: gGr1,gGr2

       real (kind=real_kind),allocatable :: solver_wts(:,:)  ! solver weights array for nonstaggered grid

       integer              :: nm1,n0,np1,nfilt
       integer              :: nstep
       integer              :: i,j,k,l,ie,kptr

       call t_startf('prim_advance_si', t_detail_low)

       nm1   = tl%nm1
       n0    = tl%n0
       np1   = tl%np1
       nstep = tl%nstep


       if ( dt /= initialized_for_dt ) then
          if(hybrid%par%masterproc) print *,'Initializing semi-implicit matricies for dt=',dt

          !$OMP MASTER
          call set_vert_struct_mat(dt, refstate, hvcoord, hybrid%masterthread)
          !$OMP END MASTER

          allocate(solver_wts(np*np,nete-nets+1))
          do ie=nets,nete
             kptr=1
             do j=1,np
                do i=1,np

                   ! so this code is BFB  with old code.  should change to simpler formula below
                   solver_wts(kptr,ie-nets+1) = 1d0/nint(1d0/(elem(ie)%mp(i,j)*elem(ie)%rmp(i,j)))
                   !solver_wts(kptr,ie-nets+1) = elem(ie)%mp(i,j)*elem(ie)%rmp(i,j)

                   kptr=kptr+1
                end do
             end do
          end do
          call cg_create(cg, np*np, nlev, nete-nets+1, hybrid, debug_level, solver_wts)
          deallocate(solver_wts)
          if (precon_method == "block_jacobi") then
             if (.not. allocated(blkjac)) then
                allocate(blkjac(nets:nete))
             endif
             call blkjac_init(elem, deriv,refstate%Lambda,nets,nete,blkjac)
          end if
          initialized_for_dt = dt
       endif


       nfilt = tl%nm1     ! time level at which filter is applied (time level n)
       dt2   = 2.0_real_kind*dt
       rdt   = 1.0_real_kind/dt

       ps0      = hvcoord%ps0
       psref    = refstate%psr

       Emat     => refstate%Emat
       Emat_inv => refstate%Emat_inv
       Amat     => refstate%Amat
       Amat_inv => refstate%Amat_inv
       Lambda   => refstate%Lambda

       RTref    => refstate%RTref
       Tref     => refstate%Tref
       Href     => refstate%Href
       Tmat     => refstate%Tmat
       Pvec     => refstate%Pvec

       ! ============================================================
       ! If the time is right, apply a filter to the state variables
       ! ============================================================

       if (nstep > 0 .and. filter_freq > 0 .and. MODULO(nstep,filter_freq) == 0 ) then
          call preq_filter(elem, edge3p1, flt, cg%hybrid, nfilt, nets, nete)
       end if

       ! ================================================
       ! boundary exchange grad_lnps
       ! ================================================

       do ie = nets, nete

          call gradient_routine(elem(ie)%state%lnps(:,:,n0),deriv,elem(ie)%derived%grad_lnps(:,:,:))
          elem(ie)%derived%grad_lnps(:,:,:) = elem(ie)%derived%grad_lnps(:,:,:) * rrearth

          do k=1,nlevp
             pintref(k)  = hvcoord%hyai(k)*ps0 + hvcoord%hybi(k)*psref
          end do

          do k=1,nlev
             pmidref(k)  = hvcoord%hyam(k)*ps0 + hvcoord%hybm(k)*psref
             pdelref(k)  = pintref(k+1) - pintref(k)
             rpmidref(k) = 1.0_real_kind/pmidref(k)
             rpdelref(k) = 1.0_real_kind/pdelref(k)
          end do

          rpsref   = 1.0_real_kind/psref

          ps(:,:) = EXP(elem(ie)%state%lnps(:,:,n0))
          rps(:,:) = 1.0_real_kind/ps(:,:)

          call preq_pressure(ps0,ps,hvcoord%hyai,hvcoord%hybi,hvcoord%hyam,hvcoord%hybm,pint,pmid,pdel)

          rpmid = 1.0_real_kind/pmid
          rpdel = 1.0_real_kind/pdel

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(v1,v2,vco)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! Contravariant velocities

                   vco(i,j,1) = elem(ie)%Dinv(i,j,1,1)*v1 + elem(ie)%Dinv(i,j,1,2)*v2
                   vco(i,j,2) = elem(ie)%Dinv(i,j,2,1)*v1 + elem(ie)%Dinv(i,j,2,2)*v2

                   vgrad_ps(i,j,k) = ps(i,j)*(vco(i,j,1)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        vco(i,j,2)*elem(ie)%derived%grad_lnps(i,j,2))

                end do
             end do
          end do
          !
          ! WARNING: omega_p changed to omega - need to divide by p !!!!!
          !
          call preq_omegap(elem(ie)%derived%div(:,:,:,n0),vgrad_ps,pdel,rpmid, &
               hvcoord%hybm,hvcoord%hybd,elem(ie)%derived%omega)

          Pscript(:,:,ie)        = 0.0_real_kind
          eta_dot_dp_deta(:,:,1) = 0.0_real_kind

          do k=1,nlev
             do j=1,np
                do i=1,np
                   eta_dot_dp_deta(i,j,k+1) = eta_dot_dp_deta(i,j,k) + &
                        vgrad_ps(i,j,k)*hvcoord%hybd(k) + elem(ie)%derived%div(i,j,k,n0)*pdel(i,j,k)
                   ddiv = elem(ie)%derived%div(i,j,k,n0) - 0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1)
                   Pscript(i,j,ie) = Pscript(i,j,ie) + ddiv*pdelref(k)
                end do
             end do
          end do

          do j=1,np
             do i=1,np
                Pscript(i,j,ie) = elem(ie)%state%lnps(i,j,nm1) + &
                     dt2*( rpsref*Pscript(i,j,ie) - rps(i,j)*eta_dot_dp_deta(i,j,nlev+1) )
             end do
          end do

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
          do k=1,nlev-1
             do j=1,np
                do i=1,np
                   eta_dot_dp_deta(i,j,k+1) = hvcoord%hybi(k+1)*eta_dot_dp_deta(i,j,nlev+1) - &
                        eta_dot_dp_deta(i,j,k+1)
                end do
             end do
          end do

          eta_dot_dp_deta(:,:,nlev+1) = 0.0_real_kind

          call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
               eta_dot_dp_deta,rpdel,T_vadv,v_vadv)

          suml(:,:) = 0.0_real_kind

          do k=1,nlev

             call gradient_routine(elem(ie)%state%T(:,:,k,n0),deriv,gradT(:,:,:))
             gradT(:,:,:) = gradT(:,:,:) * rrearth
             Crkk       = 0.5_real_kind

             do j=1,np
                do i=1,np
                   term = Crkk*(elem(ie)%derived%div(i,j,k,n0) - &
                        0.5_real_kind*elem(ie)%derived%div(i,j,k,nm1))*pdelref(k)
                   suml(i,j)  = suml(i,j) + term

                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! Contravariant velocities

                   vco(i,j,1) = elem(ie)%Dinv(i,j,1,1)*v1 + elem(ie)%Dinv(i,j,1,2)*v2
                   vco(i,j,2) = elem(ie)%Dinv(i,j,2,1)*v1 + elem(ie)%Dinv(i,j,2,2)*v2

                   vgradT = vco(i,j,1)*gradT(i,j,1) + vco(i,j,2)*gradT(i,j,2)

                   !
                   ! WARNING: omega_p changed to omega - need to divide by p !!!!!
                   !
                   Tscript(i,j,k,ie) = elem(ie)%state%T(i,j,k,nm1) &
                        + dt2*(- vgradT - T_vadv(i,j,k)           &
                        + kappa*(elem(ie)%state%T(i,j,k,n0)*elem(ie)%derived%omega(i,j,k) &
                        + Tref(k)*rpmidref(k)*suml(i,j)))
                   suml(i,j)  = suml(i,j) + term
                end do
             end do
          end do

          HrefT(:,:)   = 0.0_real_kind
          HrefTm1(:,:) = 0.0_real_kind
          HT(:,:)      = 0.0_real_kind

          do k=nlev,1,-1

             do j=1,np
                do i=1,np
                   hkl = rpmidref(k)*pdelref(k)
                   hkk = hkl*0.5_real_kind
                   Gref0(i,j)   = HrefT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
                   HrefT(i,j)   = HrefT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
                   Grefm1(i,j)  = HrefTm1(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,nm1)
                   HrefTm1(i,j) = HrefTm1(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,nm1)
                   hkl = rpmid(i,j,k)*pdel(i,j,k)
                   hkk = hkl*0.5_real_kind
                   Phi(i,j) = HT(i,j) + Rgas*hkk*elem(ie)%state%T(i,j,k,n0)
                   HT(i,j)  = HT(i,j) + Rgas*hkl*elem(ie)%state%T(i,j,k,n0)
                end do
             end do

             do j=1,np
                do i=1,np
                   v1     = elem(ie)%state%v(i,j,1,k,n0)
                   v2     = elem(ie)%state%v(i,j,2,k,n0)

                   ! covariant velocity

                   vco(i,j,1) = elem(ie)%D(i,j,1,1)*v1 + elem(ie)%D(i,j,2,1)*v2
                   vco(i,j,2) = elem(ie)%D(i,j,1,2)*v1 + elem(ie)%D(i,j,2,2)*v2

                   E(i,j) = 0.5_real_kind*( v1*v1 + v2*v2 )

                   Gref0(i,j)  =  Gref0(i,j)  + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,n0)
                   Grefm1(i,j) =  Grefm1(i,j) + elem(ie)%state%phis(i,j) + RTref(k)*elem(ie)%state%lnps(i,j,nm1)

                   Phi(i,j)    =  Phi(i,j) + E(i,j) + elem(ie)%state%phis(i,j)
                   dGref(i,j)  =  -(Gref0(i,j)  - 0.5_real_kind*Grefm1(i,j))
                end do
             end do

             call vorticity_routine(vco,deriv,elem(ie)%derived%zeta(:,:,k))
             elem(ie)%derived%zeta(:,:,k)=elem(ie)%derived%zeta(:,:,k)*rrearth

             call gradient_routine(Phi,deriv,grad_Phi(:,:,:))
             grad_Phi(:,:,:) = grad_Phi(:,:,:) * rrearth

             call gradient_wk_routine(dGref,deriv,grad_dGref(:,:,:,k))
             grad_dGref(:,:,:,k) = grad_dGref(:,:,:,k)*rrearth

             do j=1,np
                do i=1,np

                   elem(ie)%derived%zeta(i,j,k) = elem(ie)%rmetdet(i,j)*elem(ie)%derived%zeta(i,j,k)
                   hybfac =  hvcoord%hybm(k)*(ps(i,j)*rpmid(i,j,k))

                   glnps1 = elem(ie)%Dinv(i,j,1,1)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        elem(ie)%Dinv(i,j,2,1)*elem(ie)%derived%grad_lnps(i,j,2)
                   glnps2 = elem(ie)%Dinv(i,j,1,2)*elem(ie)%derived%grad_lnps(i,j,1) + &
                        elem(ie)%Dinv(i,j,2,2)*elem(ie)%derived%grad_lnps(i,j,2)

                   v1 = elem(ie)%Dinv(i,j,1,1)*grad_Phi(i,j,1) + elem(ie)%Dinv(i,j,2,1)*grad_Phi(i,j,2)
                   v2 = elem(ie)%Dinv(i,j,1,2)*grad_Phi(i,j,1) + elem(ie)%Dinv(i,j,2,2)*grad_Phi(i,j,2)

                   Vscript(i,j,1,k,ie) = - v_vadv(i,j,1,k) &
                        + elem(ie)%state%v(i,j,2,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
                        - v1 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps1

                   Vscript(i,j,2,k,ie) = - v_vadv(i,j,2,k) &
                        - elem(ie)%state%v(i,j,1,k,n0) * (elem(ie)%fcor(i,j) + elem(ie)%derived%zeta(i,j,k)) &
                        - v2 - Rgas*hybfac*elem(ie)%state%T(i,j,k,n0)*glnps2

                end do
             end do

          end do

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(Vs1,Vs2)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   Vs1 = elem(ie)%Dinv(i,j,1,1)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(i,j,2,1)*grad_dGref(i,j,2,k)
                   Vs2 = elem(ie)%Dinv(i,j,1,2)*grad_dGref(i,j,1,k) + elem(ie)%Dinv(i,j,2,2)*grad_dGref(i,j,2,k)

                   Vscript(i,j,1,k,ie) = elem(ie)%mp(i,j)*Vscript(i,j,1,k,ie) + Vs1
                   Vscript(i,j,2,k,ie) = elem(ie)%mp(i,j)*Vscript(i,j,2,k,ie) + Vs2

                   Vscript(i,j,1,k,ie) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,1,k,nm1) + dt2*Vscript(i,j,1,k,ie)
                   Vscript(i,j,2,k,ie) = elem(ie)%mp(i,j)*elem(ie)%state%v(i,j,2,k,nm1) + dt2*Vscript(i,j,2,k,ie)
                end do
             end do

          end do

          HrefTscript(:,:) = 0.0_real_kind

          do k=nlev,1,-1

             do j=1,np
                do i=1,np
                   hkl = rpmidref(k)*pdelref(k)
                   hkk = hkl*0.5_real_kind
                   B(i,j,k,ie)      = HrefTscript(i,j) + Rgas*hkk*Tscript(i,j,k,ie)
                   B(i,j,k,ie)      = B(i,j,k,ie) +  elem(ie)%state%phis(i,j) + RTref(k)*Pscript(i,j,ie)
                   HrefTscript(i,j) = HrefTscript(i,j) + Rgas*hkl*Tscript(i,j,k,ie)
                end do
             end do

          end do

          kptr=0
          call edgeVpack(edge2, Vscript(1,1,1,1,ie),2*nlev,kptr,ie)

       end do

       call t_startf('bndry_exchangeV.edge2', t_detail_medium)
       call bndry_exchangeV(cg%hybrid,edge2,location='prim_advance_si #1')
       call t_stopf('bndry_exchangeV.edge2', t_detail_medium)

       do ie = nets, nete

          kptr=0
          call edgeVunpack(edge2, Vscript(1,1,1,1,ie), 2*nlev, kptr, ie)
#ifdef DEBUGOMP
!$OMP BARRIER
#endif

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(gVscript)
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   Vscript(i,j,1,k,ie) = elem(ie)%rmp(i,j)*Vscript(i,j,1,k,ie)
                   Vscript(i,j,2,k,ie) = elem(ie)%rmp(i,j)*Vscript(i,j,2,k,ie)
                end do
             end do

             do j=1,np
                do i=1,np

                   ! Contravariant Vscript

                   gVscript(i,j,1) = elem(ie)%Dinv(i,j,1,1)*Vscript(i,j,1,k,ie) + &
                        elem(ie)%Dinv(i,j,1,2)*Vscript(i,j,2,k,ie)
                   gVscript(i,j,2) = elem(ie)%Dinv(i,j,2,1)*Vscript(i,j,1,k,ie) + &
                        elem(ie)%Dinv(i,j,2,2)*Vscript(i,j,2,k,ie)

                   gVscript(i,j,1) = elem(ie)%metdet(i,j)*gVscript(i,j,1)
                   gVscript(i,j,2) = elem(ie)%metdet(i,j)*gVscript(i,j,2)

                end do
             end do

             call divergence_routine(gVscript,deriv,div_Vscript(:,:,k))
             div_Vscript(:,:,k)=div_Vscript(:,:,k)*rrearth

          end do

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   C(i,j,k,ie) = elem(ie)%metdet(i,j)*B(i,j,k,ie)
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      C(i,j,k,ie) = C(i,j,k,ie) - dt*Amat(l,k)*div_Vscript(i,j,l)
                   end do
                end do
             end do

          end do

          ! ===============================================================
          !  Weight C (the RHS of the helmholtz problem) by the mass matrix
          ! ===============================================================

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
          do k=1,nlev
             do j=1,np
                do i=1,np
                   C(i,j,k,ie) = elem(ie)%mp(i,j)*C(i,j,k,ie)
                end do
             end do
          end do

          ! ===================================
          ! Pack C into the edge1 buffer
          ! ===================================

          kptr=0
          call edgeVpack(edge1,C(1,1,1,ie),nlev,kptr,ie)

       end do

       ! ==================================
       ! boundary exchange C
       ! ==================================

       call t_startf('bndry_exchangeV.edge1', t_detail_medium)
       call bndry_exchangeV(cg%hybrid,edge1,location='prim_advance_si #2')
       call t_stopf('bndry_exchangeV.edge1', t_detail_medium)

       do ie=nets,nete

          ! ===================================
          ! Unpack C from the edge1 buffer
          ! ===================================

          kptr=0
          call edgeVunpack(edge1, C(1,1,1,ie), nlev, kptr, ie)

          ! ===============================================
          ! Complete global assembly by normalizing by rmp
          ! ===============================================

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   D(i,j,k,ie) = elem(ie)%rmp(i,j)*C(i,j,k,ie)
                end do
             end do

          end do

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   C(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      C(i,j,k,ie) = C(i,j,k,ie) + Emat_inv(l,k)*D(i,j,l,ie)
                   end do
                end do
             end do

          end do

       end do
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
       ! ==========================================
       ! solve for Gamma_ref, given C as RHS input
       ! ==========================================

       Gamma_ref = pcg_solver(elem, &
            C,          &
            cg,         &
            red,        &
            edge1,      &
            edge2,      &
            Lambda,     &
            deriv,      &
            nets,       &
            nete,       &
            blkjac)

       ! ================================================================
       ! Backsubstitute Gamma_ref into semi-implicit system of equations
       ! to find prognostic variables at time level n+1
       ! ================================================================

       kptr=0
       do ie = nets, nete

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   Gref(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      Gref(i,j,k,ie) = Gref(i,j,k,ie) + Emat(l,k)*Gamma_ref(i,j,l,ie)
                   end do
                end do
             end do

             do j=1,np
                do i=1,np
                   B(i,j,k,ie) = elem(ie)%mp(i,j) * dt * (B(i,j,k,ie) - Gref(i,j,k,ie))
                end do
             end do

          end do

          call edgeVpack(edge1,B(:,:,:,ie),nlev,kptr,ie)

       end do

       call t_startf('bndry_exchangeV.edge1', t_detail_medium)
       call bndry_exchangeV(cg%hybrid,edge1,location='prim_advance_si #3')
       call t_stopf('bndry_exchangeV.edge1', t_detail_medium)

       do ie = nets, nete

          kptr=0
          call edgeVunpack(edge1, B(:,:,:,ie), nlev, kptr, ie)
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   B(i,j,k,ie) = elem(ie)%rmp(i,j)*B(i,j,k,ie)
                end do
             end do

          end do

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   D(i,j,k,ie) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      D(i,j,k,ie) = D(i,j,k,ie) + Emat_inv(l,k)*B(i,j,l,ie)
                   end do
                end do
             end do

          end do

#if 1
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
                   elem(ie)%derived%div(i,j,k,np1) = 0.0_real_kind
                end do
             end do

             do l=1,nlev
                do j=1,np
                   do i=1,np
                      elem(ie)%derived%div(i,j,k,np1) = elem(ie)%derived%div(i,j,k,np1) + Emat(l,k)*D(i,j,l,ie)/Lambda(l)
                   end do
                end do
             end do

          end do
#endif

          do k=1,nlev

             call gradient_wk_routine(Gref(:,:,k,ie),deriv,grad_Gref(:,:,:,k))
             grad_Gref(:,:,:,k)=grad_Gref(:,:,:,k)*rrearth

             do j=1,np
                do i=1,np
                   gGr1 = grad_Gref(i,j,1,k)
                   gGr2 = grad_Gref(i,j,2,k)
                   elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%Dinv(i,j,1,1)*gGr1 + elem(ie)%Dinv(i,j,2,1)*gGr2
                   elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%Dinv(i,j,1,2)*gGr1 + elem(ie)%Dinv(i,j,2,2)*gGr2
                   Vtemp(i,j,1,k,ie) = elem(ie)%state%v(i,j,1,k,np1)
                   Vtemp(i,j,2,k,ie) = elem(ie)%state%v(i,j,2,k,np1)
                end do
             end do

             do j=1,np
                do i=1,np
                   Pscript(i,j,ie) = Pscript(i,j,ie) - dt*Pvec(k)*elem(ie)%derived%div(i,j,k,np1)
                end do
             end do


             do l=1,nlev
                do j=1,np
                   do i=1,np
                      Tscript(i,j,k,ie) = Tscript(i,j,k,ie) - dt*Tmat(l,k)*elem(ie)%derived%div(i,j,l,np1)
                   end do
                end do
             end do

          end do

          do j=1,np
             do i=1,np
                Pscript(i,j,ie) = elem(ie)%mp(i,j)*Pscript(i,j,ie)
             end do
          end do
          do k=1,nlev
             do j=1,np
                do i=1,np
                   Tscript(i,j,k,ie) = elem(ie)%mp(i,j)*Tscript(i,j,k,ie)
                end do
             end do
          end do

          ! ===============================================
          ! Pack v at time level n+1 into the edge3p1 buffer
          ! ===============================================

          kptr=0
!          call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)
          call edgeVpack(edge3p1, Vtemp(:,:,:,:,ie),2*nlev,kptr,ie)

          kptr=2*nlev
          call edgeVpack(edge3p1, Tscript(:,:,:,ie),nlev,kptr,ie)

          kptr=3*nlev
          call edgeVpack(edge3p1, Pscript(:,:,ie),1,kptr,ie)

       end do

       ! ======================================
       ! boundary exchange v at time level n+1
       ! ======================================

       call t_startf('bndry_exchangeV.edge3p1', t_detail_medium)
       call bndry_exchangeV(cg%hybrid,edge3p1,location='prim_advance_si #4')
       call t_stopf('bndry_exchangeV.edge3p1', t_detail_medium)

!KGEN START(prim_advance_si_bug1)
       do ie=nets,nete

          ! ===================================
          ! Unpack v from the edge2 buffer
          ! ===================================

          kptr=0
          call edgeVunpack(edge3p1, Vtemp(:,:,:,:,ie), 2*nlev, kptr, ie)
!JMD          call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)
!JMD          elem(ie)%state%v(:,:,:,:,np1) = Vtemp(:,:,:,:,ie)


          kptr=2*nlev
          call edgeVunpack(edge3p1, Tscript(:,:,:,ie), nlev, kptr, ie)

          kptr=3*nlev
          call edgeVunpack(edge3p1, Pscript(:,:,ie), 1, kptr, ie)

          ! ==========================================================
          ! Complete global assembly by normalizing velocity by rmp
          ! Vscript = Vscript - dt*grad(Gref)
          ! ==========================================================
          if(iam==1) then           
!BUG  There appears to be a bug in the Intel 15.0.1 compiler that generates
!BUG  incorrect code for this loop if the following print * statement is removed.
             print *,'IAM: ',iam, ' prim_advance_si: after SUM(v(np1)) ',sum(elem(ie)%state%v(:,:,:,:,np1)) 
          endif

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads)
#endif
          do k=1,nlev

             do j=1,np
                do i=1,np
!                   elem(ie)%state%v(i,j,1,k,np1) = Vscript(i,j,1,k,ie) + dt*elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,1,k,np1)
!                   elem(ie)%state%v(i,j,2,k,np1) = Vscript(i,j,2,k,ie) + dt*elem(ie)%rmp(i,j)*elem(ie)%state%v(i,j,2,k,np1)
                   elem(ie)%state%v(i,j,1,k,np1) = Vscript(i,j,1,k,ie) + dt*elem(ie)%rmp(i,j)*Vtemp(i,j,1,k,ie)
                   elem(ie)%state%v(i,j,2,k,np1) = Vscript(i,j,2,k,ie) + dt*elem(ie)%rmp(i,j)*Vtemp(i,j,2,k,ie)
                   elem(ie)%state%T(i,j,k,np1)   = elem(ie)%rmp(i,j)*Tscript(i,j,k,ie)
                end do
             end do

          end do

          do j=1,np
             do i=1,np
                elem(ie)%state%lnps(i,j,np1) = elem(ie)%rmp(i,j)*Pscript(i,j,ie)
             end do
          end do

       end do
!KGEN END(prim_advance_si_bug1)
#ifdef DEBUGOMP
!$OMP BARRIER
#endif

       call prim_diffusion(elem, nets,nete,np1,deriv,dt2,cg%hybrid)
       call t_stopf('prim_advance_si', t_detail_low)
#endif
  end subroutine prim_advance_si


  subroutine preq_robert3(nm1,n0,np1,elem,hvcoord,nets,nete)
  use dimensions_mod, only : np, nlev, qsize
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use time_mod, only: smooth
  use control_mod, only : integration

  implicit none
  integer              , intent(in) :: nm1,n0,np1,nets,nete
  type (hvcoord_t), intent(in)      :: hvcoord
  type (element_t)     , intent(inout) :: elem(:)


  integer :: i,j,k,ie,q
  real (kind=real_kind) :: dp
  logical :: filter_ps = .false.
  if (integration == "explicit") filter_ps = .true.

  call t_startf('preq_robert', t_detail_high)
  do ie=nets,nete
     if (filter_ps) then
        elem(ie)%state%ps(:,:,n0) = elem(ie)%state%ps(:,:,n0) + smooth*(elem(ie)%state%ps(:,:,nm1) &
             - 2.0D0*elem(ie)%state%ps(:,:,n0)   + elem(ie)%state%ps(:,:,np1))
        elem(ie)%state%lnps(:,:,n0) = LOG(elem(ie)%state%ps(:,:,n0))
     else
        elem(ie)%state%lnps(:,:,n0) = elem(ie)%state%lnps(:,:,n0) + smooth*(elem(ie)%state%lnps(:,:,nm1) &
             - 2.0D0*elem(ie)%state%lnps(:,:,n0)   + elem(ie)%state%lnps(:,:,np1))
        elem(ie)%state%ps(:,:,n0) = EXP(elem(ie)%state%lnps(:,:,n0))
     endif

     elem(ie)%state%T(:,:,:,n0) = elem(ie)%state%T(:,:,:,n0) + smooth*(elem(ie)%state%T(:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%T(:,:,:,n0)   + elem(ie)%state%T(:,:,:,np1))
     elem(ie)%state%v(:,:,:,:,n0) = elem(ie)%state%v(:,:,:,:,n0) + smooth*(elem(ie)%state%v(:,:,:,:,nm1) &
          - 2.0D0*elem(ie)%state%v(:,:,:,:,n0) + elem(ie)%state%v(:,:,:,:,np1))

  end do
  call t_stopf('preq_robert', t_detail_high)

  end subroutine preq_robert3




  subroutine applyCAMforcing(elem,fvm,hvcoord,np1,np1_qdp,dt_q,nets,nete,nsubstep)
  use dimensions_mod,         only : np, nc, nlev, qsize, ntrac, ldry_mass_vertical_coordinates,nelemd
  use element_mod,            only : element_t
  use hybvcoord_mod,          only : hvcoord_t
  use time_mod,               only : nsplit
  use control_mod,            only : moisture, tracer_grid_type,ftype
  use control_mod,            only : TRACER_GRIDTYPE_GLL, TRACER_GRIDTYPE_FVM
  use physical_constants,     only : Cp
  use fvm_control_volume_mod, only : fvm_struct, n0_fvm, np1_fvm

  implicit none
  type (element_t)     , intent(inout) :: elem(:)
  type(fvm_struct)     , intent(inout) :: fvm(:)
  real (kind=real_kind), intent(in) :: dt_q
  type (hvcoord_t), intent(in)      :: hvcoord
  integer,  intent(in) :: np1,nets,nete,np1_qdp,nsubstep

  ! local
  integer :: i,j,k,ie,q
  integer :: num_threads
  real (kind=real_kind) :: v1,fdq, fdpq, dp,dt_local, q_physics_updated, q_old 
  real (kind=real_kind) :: ftmp(np,np,nlev,qsize,nelemd) !diagnostics

  logical :: wet

  wet = (moisture /= "dry")
  ! ftype=1: apply all forcings as an adjustment       
  ! ftype=0: apply forcings after every vertical remap
  if (ftype==0) then
    dt_local = dt_q
  else if (ftype==1) then
    !
    ! CAM-FV-stype forcing, i.e. equivalent to updating state once in the
    ! beginning of dynamics
    !
    dt_local = nsplit*dt_q
    if (nsubstep.ne.1) then
#ifdef CAM
      ftmp(:,:,:,:,:) = 0.0D0
      call output_qdp_var_dynamics(ftmp(:,:,:,:,:),nets,nete,'PDC')
      call calc_tot_energy_dynamics(elem,nets,nete,np1,np1_qdp,'dBM')
      call calc_tot_energy_dynamics(elem,nets,nete,np1,np1_qdp,'dBD')
      return !forcing already added
#endif
    end if
  end if

  do ie=nets,nete
    elem(ie)%state%T(:,:,:,np1) = elem(ie)%state%T(:,:,:,np1) + &
         dt_local*elem(ie)%derived%FT(:,:,:,1)
    elem(ie)%state%v(:,:,:,:,np1) = elem(ie)%state%v(:,:,:,:,np1) + &
         dt_local*elem(ie)%derived%FM(:,:,:,:,1)

    ! apply forcing to Qdp
    !
    ! for fvm tracer qsize is usually 1 (qv)
    !
    !
    ! tracers
    !
#if (defined _LOOP_OPENMP)
    num_threads = max(vert_num_threads,tracer_num_threads)
    !$omp parallel do num_threads(num_threads), private(v1) 
#endif
    do q=1,qsize
      do k=1,nlev
        do j=1,np
          do i=1,np
            !
            ! FQ holds mass-weighted q-tendency: (qnew-qold)*dp
            !
            v1 = dt_local*elem(ie)%derived%FQ(i,j,k,q,1)
            if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) + v1 < 0 .and. v1<0) then
              if (elem(ie)%state%Qdp(i,j,k,q,np1_qdp) < 0 ) then
                v1=0  ! Q already negative, dont make it more so
              else
                !v1 = -elem(ie)%state%Qdp(i,j,k,q,np1)
                v1 = -elem(ie)%state%Qdp(i,j,k,q,np1_qdp)
              endif
            endif
            elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)+v1
            ftmp(i,j,k,q,ie) = dt_local*elem(ie)%derived%FQ(i,j,k,q,1)-v1 !Only used for diagnostics!
          enddo
        enddo
      enddo
    enddo
  end do
  if (qsize>0.and.ldry_mass_vertical_coordinates) then
    do ie=nets,nete
      do q=1,qsize          
        do k=1,nlev
          do j=1,np
            do i=1,np
              elem(ie)%state%Q(i,j,k,q) = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)/&
                   elem(ie)%state%dp3d(i,j,k,np1)
            end do
          end do
        end do
      end do
    end do
  end if
  if (ntrac>0) then
    do ie=nets,nete
      !
      ! Repeat for the fvm tracers
      !
      do q = 1, ntrac
        do k = 1, nlev
          do j = 1, nc
            do i = 1, nc
              !
              ! note that fc is mass weighted
              !
              v1 = dt_local*fvm(ie)%fc(i,j,k,q)/fvm(ie)%dp_fvm(i,j,k,n0_fvm)
              if (fvm(ie)%c(i,j,k,q,n0_fvm) + v1 < 0 .and. v1<0) then
                if (fvm(ie)%c(i,j,k,q,n0_fvm) < 0 ) then
                  v1 = 0  ! C already negative, dont make it more so
                else
                  v1 = -fvm(ie)%c(i,j,k,q,n0_fvm)
                end if
              end if
              fvm(ie)%c(i,j,k,q,n0_fvm) = fvm(ie)%c(i,j,k,q,n0_fvm)+ v1                 
            end do
          end do
        end do
      end do    
    enddo
  end if


#ifdef CAM
  call output_qdp_var_dynamics(ftmp(:,:,:,:,:),nets,nete,'PDC')
  call calc_tot_energy_dynamics(elem,nets,nete,np1,np1_qdp,'dBM')
  if (qsize>0.and..not.ldry_mass_vertical_coordinates.and.wet) then
    do ie=nets,nete
      do k=1,nlev
        do j=1,np
          do i=1,np
            q_physics_updated = elem(ie)%state%Qdp(i,j,k,1,np1_qdp)/elem(ie)%state%dp3d(i,j,k,np1)
            q_old             = elem(ie)%state%Q(i,j,k,1)
            !
            ! code quivalent to dme_adjust in physics_types.F90
            !
            dp  = elem(ie)%state%dp3d(i,j,k,np1)
            fdq = 1.0D0+q_physics_updated-q_old
            elem(ie)%state%dp3d(i,j,k,np1) = dp*fdq
            do q=1,qsize
              q_physics_updated = elem(ie)%state%Qdp(i,j,k,q,np1_qdp)/dp
              elem(ie)%state%Q(i,j,k,q) = q_physics_updated/fdq
              elem(ie)%state%Qdp(i,j,k,q,np1_qdp) = &
                   elem(ie)%state%Q(i,j,k,q)*elem(ie)%state%dp3d(i,j,k,np1)
            end do
          end do
        end do
      end do
    end do
  end if
#endif


#ifdef CAM
    call calc_tot_energy_dynamics(elem,nets,nete,np1,np1_qdp,'dBD')
#endif

end subroutine applyCAMforcing

  subroutine advance_hypervis(edge3,elem,hvcoord,hybrid,deriv,nt,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t
  use derivative_mod, only : laplace_sphere_wk_routine
  use derivative_mod, only:  vlaplace_sphere_wk_routine
  use edge_mod, only : edgevpack, edgevunpack
  use edgetype_mod, only : EdgeBuffer_t
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk
  use physical_constants, only: Cp
!  use time_mod, only : TimeLevel_t
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic,nt
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
!JMD  call t_barrierf('sync_advance_hypervis', hybrid%par%comm)
  call t_startf('advance_hypervis', t_detail_low)


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(lap_p,lap_v)
#endif
           do k=1,nlev
              call laplace_sphere_wk_routine (elem(ie)%state%T(:,:,k,nt),  deriv,elem(ie),var_coef=.false.,laplace=lap_p)
              call vlaplace_sphere_wk_routine(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.,laplace=lap_v)
              ! advance in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
        enddo

        call t_startf('bndry_exchangeV.edge3', t_detail_medium)
        call bndry_exchangeV(hybrid,edge3,location='prim_advance_mod:1872')
        call t_startf('bndry_exchangeV.edge3', t_detail_medium)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)

           ! apply inverse mass matrix
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nu_p=0:
!   scale T dissipaton by dp  (conserve IE, dissipate T^2)
! nu_p>0
!   dont scale:  T equation IE dissipation matches (to truncation error)
!                IE dissipation from continuity equation
!                (1 deg: to about 0.1 W/m^2)
!
!  print *,'iam: ',iam,' advance_hypervis: before loop: SUM(elem(1)%state%v(:,:,:,:,nt)): ', SUM(elem(1)%state%v(:,:,:,:,nt))
!  print *,'iam: ',iam,' advance_hypervis: before loop: SUM(elem(1)%state%T(:,:,:,nt)): ', SUM(elem(1)%state%T(:,:,:,nt))
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,1,nlev)
        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
#if 0
              elem(ie)%derived%psdiss_ave(:,:)=&
                   elem(ie)%derived%psdiss_ave(:,:)+eta_ave_w*elem(ie)%state%ps(:,:,nt)/hypervis_subcycle
              elem(ie)%derived%psdiss_biharmonic(:,:)=&
                   elem(ie)%derived%psdiss_biharmonic(:,:)+eta_ave_w*pstens(:,:,ie)/hypervis_subcycle
#else
              do k=1,nlev
                 dptemp1(:,:) = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                      ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(:,:,nt)
                 elem(ie)%derived%dpdiss_ave(:,:,k)=elem(ie)%derived%dpdiss_ave(:,:,k)+eta_ave_w*dptemp1(:,:)/hypervis_subcycle

                 dptemp2(:,:) = (hvcoord%hybi(k+1)-hvcoord%hybi(k))*pstens(:,:,ie)
                 elem(ie)%derived%dpdiss_biharmonic(:,:,k)=&
                      elem(ie)%derived%dpdiss_biharmonic(:,:,k)+eta_ave_w*dptemp2(:,:)/hypervis_subcycle
              enddo
#endif
           endif
           nu_scale=1
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) &
!$omp private(lap_p,lap_v,nu_scale_top,dpdn,dpdn0,nu_scale,utens_tmp,vtens_tmp,ptens_tmp)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 call laplace_sphere_wk_routine (elem(ie)%state%T(:,:,k,nt),  deriv,elem(ie),var_coef=.false.,laplace=lap_p)
                 call vlaplace_sphere_wk_routine(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.,laplace=lap_v)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    if (nu_p==0) then
                       ! normalize so as to conserve IE
                       ! scale by 1/rho (normalized to be O(1))
                       ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                       dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(i,j,nt)
                       dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                       nu_scale = dpdn0/dpdn
                    endif

                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                    endif

                    ptens(i,j,k,ie) = ptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo

           pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
           kptr=0
           call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=3*nlev
           call edgeVpack(edge3,pstens(:,:,ie),1,kptr,ie)
        enddo


        call t_startf('bndry_exchangeV.edge3', t_detail_medium)
        call bndry_exchangeV(hybrid,edge3,location='prim_advance_mod:1998')
        call t_startf('bndry_exchangeV.edge3', t_detail_medium)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)


           ! apply inverse mass matrix, accumulate tendencies
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads)
#endif
           do k=1,nlev
              vtens(:,:,1,k,ie)=dt*vtens(:,:,1,k,ie)*elem(ie)%rspheremp(:,:)
              vtens(:,:,2,k,ie)=dt*vtens(:,:,2,k,ie)*elem(ie)%rspheremp(:,:)
              ptens(:,:,k,ie)=dt*ptens(:,:,k,ie)*elem(ie)%rspheremp(:,:)
           enddo

           ! apply hypervis to u -> u+utens:
           ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
           ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
           ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
           !      X = (u dot utens) + .5 utens dot utens
           !  alt:  (u+utens) dot utens
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    ! update v first (gives better results than updating v after heating)
                    elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                         vtens(i,j,:,k,ie)

                    v1=elem(ie)%state%v(i,j,1,k,nt)
                    v2=elem(ie)%state%v(i,j,2,k,nt)
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         +ptens(i,j,k,ie)-heating/cp

                 enddo
              enddo
           enddo

           if (nu_p>0) then
              kptr=3*nlev
              call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, ie)
              pstens(:,:,ie)=dt*pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
              elem(ie)%state%ps(:,:,nt)=elem(ie)%state%ps(:,:,nt) + pstens(:,:,ie)
           endif

        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo
  endif

  call t_stopf('advance_hypervis', t_detail_low)

  end subroutine advance_hypervis

  subroutine advance_hypervis_dp(edge3,elem,hvcoord,hybrid,deriv,nt,qn0,nets,nete,dt2,eta_ave_w)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev, nc, ntrac, max_corner_elem
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis, swest
  use hybrid_mod, only : hybrid_t, get_loop_ranges
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t
  use derivative_mod, only : laplace_sphere_wk_routine
  use derivative_mod, only : vlaplace_sphere_wk_routine
  use derivative_mod, only : subcell_Laplace_fluxes, subcell_dss_fluxes
  use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
  use edgetype_mod, only : EdgeBuffer_t, EdgeDescriptor_t
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk_dp3d
  use physical_constants, only: Cp
  use derivative_mod, only : subcell_Laplace_fluxes
!  use time_mod, only : TimeLevel_t
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
  real (kind=real_kind) :: eta_ave_w  ! weighting for mean flux terms
  real (kind=real_kind) :: dt2
  integer :: nets,nete, nt, qn0

  ! local
  real (kind=real_kind) :: dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic !qn0 is only for diagnostics
  integer :: kbeg, kend, kblk
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ttens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: dptens
  real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)          :: corners
  real (kind=real_kind), dimension(2,2,2)                       :: cflux
  real (kind=real_kind), dimension(nc,nc,4,nlev,nets:nete)      :: dpflux
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dptemp1,dptemp2
  type (EdgeDescriptor_t)                                       :: desc

! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
!       data is incorrect (offset by a few numbers actually)
!       removed for now.
!       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
!       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_t,lap_dp
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ttens_tmp,dptens_tmp

  real (kind=real_kind)                     :: temp      (np,np,nlev)
  real (kind=real_kind)                     :: laplace_fluxes(nc,nc,4)

  real (kind=real_kind) :: rhypervis_subcycle

  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
  call t_startf('advance_hypervis_dp', t_detail_low)

  call get_loop_ranges(hybrid,kbeg=kbeg,kend=kend)

  kblk = kend - kbeg + 1

  dt=dt2/hypervis_subcycle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (hypervis_order == 1) then
     if (nu_p>0) call abortmp( 'ERROR: hypervis_order == 1 not coded for nu_p>0')
     do ic=1,hypervis_subcycle
        do ie=nets,nete

           do k=kbeg,kend
              call laplace_sphere_wk_routine (elem(ie)%state%T(:,:,k,nt),  deriv,elem(ie),var_coef=.false.,laplace=lap_t)
              call vlaplace_sphere_wk_routine(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.,laplace=lap_v)

              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_t(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr = kbeg - 1
           call edgeVpack(edge3,elem(ie)%state%T(:,:,kbeg:kend,nt),kblk,kptr,ie)

           kptr= kbeg - 1 + nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,1,kbeg:kend,nt),kblk,kptr,ie)

           kptr= kbeg - 1 + 2*nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,2,kbeg:kend,nt),kblk,kptr,ie)
        enddo

        call t_startf('bndry_exchangeV.edge3', t_detail_medium)
        call bndry_exchangeV(hybrid,edge3,location='advance_hypervis_dp1')
        call t_stopf('bndry_exchangeV.edge3', t_detail_medium)

        do ie=nets,nete

           kptr = kbeg - 1
           call edgeVunpack(edge3,elem(ie)%state%T(:,:,kbeg:kend,nt),kblk,kptr,ie)

           kptr= kbeg - 1 + nlev
           call edgeVunpack(edge3,elem(ie)%state%v(:,:,1,kbeg:kend,nt),kblk,kptr,ie)

           kptr= kbeg - 1 + 2*nlev
           call edgeVunpack(edge3,elem(ie)%state%v(:,:,2,kbeg:kend,nt),kblk,kptr,ie)

           ! apply inverse mass matrix
           do k=kbeg,kend
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
     enddo  ! subcycle
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! nu_p=0:
!   scale T dissipaton by dp  (conserve IE, dissipate T^2)
! nu_p>0
!   dont scale:  T equation IE dissipation matches (to truncation error)
!                IE dissipation from continuity equation
!                (1 deg: to about 0.1 W/m^2)
!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
#ifdef CAM
      call calc_tot_energy_dynamics(elem,nets,nete,nt,qn0,'dBH')
#endif

        rhypervis_subcycle=1.0d0/real(hypervis_subcycle,kind=real_kind)
        call biharmonic_wk_dp3d(elem,dptens,dpflux,ttens,vtens,deriv,edge3,hybrid,nt,nets,nete,kbeg,kend)

        do ie=nets,nete

           ! comptue mean flux
           if (nu_p>0) then
              do k=kbeg,kend
#ifdef _B4B
                !OMP_COLLAPSE_SIMD
                !DIR_VECTOR_ALIGNED
                do j=1,np
                  do i=1,np
                    elem(ie)%derived%dpdiss_ave(i,j,k)=elem(ie)%derived%dpdiss_ave(i,j,k)+&
                                 eta_ave_w*elem(ie)%state%dp3d(i,j,k,nt)/hypervis_subcycle
                    elem(ie)%derived%dpdiss_biharmonic(i,j,k)=elem(ie)%derived%dpdiss_biharmonic(i,j,k)+&
                                 eta_ave_w*dptens(i,j,k,ie)/hypervis_subcycle
                  enddo
                enddo
#else
                !OMP_COLLAPSE_SIMD
                !DIR_VECTOR_ALIGNED
                do j=1,np
                  do i=1,np
                    elem(ie)%derived%dpdiss_ave(i,j,k)=elem(ie)%derived%dpdiss_ave(i,j,k)+&
                                 rhypervis_subcycle*eta_ave_w*elem(ie)%state%dp3d(i,j,k,nt)
                    elem(ie)%derived%dpdiss_biharmonic(i,j,k)=elem(ie)%derived%dpdiss_biharmonic(i,j,k)+&
                                 rhypervis_subcycle*eta_ave_w*dptens(i,j,k,ie)
                  enddo
                enddo
#endif
              enddo
           endif
           do k=kbeg,kend
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 call laplace_sphere_wk_routine (elem(ie)%state%T(:,:,k,nt),   deriv,elem(ie),var_coef=.false.,laplace=lap_t)
                 call laplace_sphere_wk_routine (elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),var_coef=.false.,laplace=lap_dp)
                 call vlaplace_sphere_wk_routine(elem(ie)%state%v(:,:,:,k,nt), deriv,elem(ie),var_coef=.false.,laplace=lap_v)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              ! biharmonic terms need a negative sign:
              if (nu_top>0 .and. k<=3) then 
                 !OMP_COLLAPSE_SIMD
                 !DIR_VECTOR_ALIGNED
                 do j=1,np
                   do i=1,np
                     ttens(i,j,k,ie)   = (-nu_s*ttens(i,j,k,ie) + nu_scale_top*nu_top*lap_t(i,j) )
                     dptens(i,j,k,ie)  = (-nu_p*dptens(i,j,k,ie) + nu_scale_top*nu_top*lap_dp(i,j) )
                     vtens(i,j,1,k,ie) = (-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                     vtens(i,j,2,k,ie) = (-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                   enddo
                 enddo
              else
                 !OMP_COLLAPSE_SIMD
                 !DIR_VECTOR_ALIGNED
                 do j=1,np
                   do i=1,np
                     ttens(i,j,k,ie)   = -nu_s*ttens(i,j,k,ie)
                     dptens(i,j,k,ie)  = -nu_p*dptens(i,j,k,ie)
                     vtens(i,j,1,k,ie) = -nu*vtens(i,j,1,k,ie)
                     vtens(i,j,2,k,ie) = -nu*vtens(i,j,2,k,ie)
                   enddo
                 enddo
              endif

              if (0<ntrac) then 
#ifdef _B4B
                !OMP_COLLAPSE_SIMD
                !DIR_VECTOR_ALIGNED
                do j=1,np
                  do i=1,np
                    elem(ie)%sub_elem_mass_flux(i,j,:,k) = elem(ie)%sub_elem_mass_flux(i,j,:,k) - &
                                              eta_ave_w*nu_p*dpflux(i,j,:,k,ie)/hypervis_subcycle
                  enddo
                enddo
                if (nu_top>0 .and. k<=3) then
                  laplace_fluxes=subcell_Laplace_fluxes(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),np,nc)
                  elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                                           eta_ave_w*nu_scale_top*nu_top*laplace_fluxes/hypervis_subcycle
                endif
#else
               !OMP_COLLAPSE_SIMD
               !DIR_VECTOR_ALIGNED
               do j=1,np
                 do i=1,np
                   elem(ie)%sub_elem_mass_flux(i,j,:,k) = elem(ie)%sub_elem_mass_flux(i,j,:,k) - &
                                            rhypervis_subcycle*eta_ave_w*nu_p*dpflux(i,j,:,k,ie)
                 enddo
               enddo
               if (nu_top>0 .and. k<=3) then
                  laplace_fluxes=subcell_Laplace_fluxes(elem(ie)%state%dp3d(:,:,k,nt),deriv,elem(ie),np,nc)
                  elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                                           rhypervis_subcycle*eta_ave_w*nu_scale_top*nu_top*laplace_fluxes
               endif
#endif
              endif

              ! NOTE: we will DSS all tendicies, EXCEPT for dp3d, where we DSS the new state
              !OMP_COLLAPSE_SIMD 
              !DIR_VECTOR_ALIGNED
              do j=1,np
                do i=1,np
                  elem(ie)%state%dp3d(i,j,k,nt) = elem(ie)%state%dp3d(i,j,k,nt)*elem(ie)%spheremp(i,j) &
                                                                                 + dt*dptens(i,j,k,ie)
                enddo
              enddo

           enddo

           kptr = kbeg - 1
           call edgeVpack(edge3,ttens(:,:,kbeg:kend,ie),kblk,kptr,ie)

           kptr = kbeg - 1 + nlev
           call edgeVpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)

           kptr = kbeg - 1 + 2*nlev
           call edgeVpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)

           kptr = kbeg - 1 + 3*nlev
           call edgeVpack(edge3,elem(ie)%state%dp3d(:,:,kbeg:kend,nt),kblk,kptr,ie)
        enddo

        call t_startf('bndry_exchangeV.edge3', t_detail_medium)
        call bndry_exchangeV(hybrid,edge3,location='advance_hypervis_dp2')
        call t_stopf('bndry_exchangeV.edge3', t_detail_medium)

        do ie=nets,nete

           kptr = kbeg - 1
           call edgeVunpack(edge3,ttens(:,:,kbeg:kend,ie),kblk,kptr,ie)

           kptr = kbeg - 1 + nlev
           call edgeVunpack(edge3,vtens(:,:,1,kbeg:kend,ie),kblk,kptr,ie)

           kptr = kbeg - 1 + 2*nlev
           call edgeVunpack(edge3,vtens(:,:,2,kbeg:kend,ie),kblk,kptr,ie)

           kptr = kbeg - 1 + 3*nlev
           call edgeVunpack(edge3,elem(ie)%state%dp3d(:,:,kbeg:kend,nt),kblk,kptr,ie)

           if (0<ntrac) then
             do k=kbeg,kend
               temp(:,:,k) = elem(ie)%state%dp3d(:,:,k,nt) / elem(ie)%spheremp  ! STATE before DSS 
               corners(0:np+1,0:np+1,k) = 0.0d0
               corners(1:np  ,1:np  ,k) = elem(ie)%state%dp3d(1:np,1:np,k,nt) ! fill in interior data of STATE*mass
             enddo
           endif

           if (0<ntrac) then
             desc = elem(ie)%desc
             
             kptr = kbeg - 1 + 3*nlev
             call edgeDGVunpack(edge3,corners(:,:,kbeg:kend),kblk,kptr,ie) 
             do k=kbeg,kend
               !OMP_COLLAPSE_SIMD
               !DIR_VECTOR_ALIGNED
               do j=1,np
                 do i=1,np
                   corners(i,j,k) = corners(i,j,k)/dt

                   temp(i,j,k) =  elem(ie)%rspheremp(i,j)*elem(ie)%state%dp3d(i,j,k,nt) - temp(i,j,k)
                   temp(i,j,k) =  temp(i,j,k)/dt
                 enddo
               enddo

               call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)
 
               cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)  
               cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:) 
               cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:) 
               cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:) 

#ifdef _B4B
               elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                 eta_ave_w*subcell_dss_fluxes(temp(:,:,k), np, nc, elem(ie)%metdet,cflux)/hypervis_subcycle
#else
               elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + &
                 rhypervis_subcycle*eta_ave_w*subcell_dss_fluxes(temp(:,:,k), np, nc, elem(ie)%metdet,cflux)
#endif
             end do
           endif

           ! apply inverse mass matrix, accumulate tendencies
           do k=kbeg,kend
             !OMP_COLLAPSE_SIMD 
             !DIR_VECTOR_ALIGNED
             do j=1,np
               do i=1,np
                 vtens(i,j,1,k,ie)=dt*vtens(i,j,1,k,ie)*elem(ie)%rspheremp(i,j)
                 vtens(i,j,2,k,ie)=dt*vtens(i,j,2,k,ie)*elem(ie)%rspheremp(i,j)
                 ttens(i,j,k,ie)=dt*ttens(i,j,k,ie)*elem(ie)%rspheremp(i,j)

                 elem(ie)%state%dp3d(i,j,k,nt)=elem(ie)%state%dp3d(i,j,k,nt)*elem(ie)%rspheremp(i,j)
               enddo
             enddo
           enddo

           ! apply hypervis to u -> u+utens:
           ! E0 = dpdn * .5*u dot u + dpdn * T  + dpdn*PHIS
           ! E1 = dpdn * .5*(u+utens) dot (u+utens) + dpdn * (T-X) + dpdn*PHIS
           ! E1-E0:   dpdn (u dot utens) + dpdn .5 utens dot utens   - dpdn X
           !      X = (u dot utens) + .5 utens dot utens
           !  alt:  (u+utens) dot utens
           do k=kbeg,kend
              !OMP_COLLAPSE_SIMD 
              !DIR_VECTOR_ALIGNED
              do j=1,np
                 do i=1,np
                    ! update v first (gives better results than updating v after heating)
                    elem(ie)%state%v(i,j,:,k,nt)=elem(ie)%state%v(i,j,:,k,nt) + &
                         vtens(i,j,:,k,ie)

!phl                    v1=elem(ie)%state%v(i,j,1,k,nt)
!phl                    v2=elem(ie)%state%v(i,j,2,k,nt)
!phl                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         +ttens(i,j,k,ie)!phl -heating/cp
                 enddo
              enddo
           enddo
         end do !phl

#ifdef CAM
      call calc_tot_energy_dynamics(elem,nets,nete,nt,qn0,'dCH')
#endif
         do ie=nets,nete
           do k=kbeg,kend
              !OMP_COLLAPSE_SIMD
              !DIR_VECTOR_ALIGNED
              do j=1,np
                 do i=1,np
                    v1=elem(ie)%state%v(i,j,1,k,nt)
                    v2=elem(ie)%state%v(i,j,2,k,nt)
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt) &
                         -heating/cp
                 enddo
              enddo
            enddo
          enddo

#ifdef CAM
      call calc_tot_energy_dynamics(elem,nets,nete,nt,qn0,'dAH')
#endif
    enddo
  endif

  call t_stopf('advance_hypervis_dp', t_detail_low)

  end subroutine advance_hypervis_dp
  
  subroutine advance_hypervis_lf(edge3,elem,hvcoord,hybrid,deriv,nm1,n0,nt,nets,nete,dt2)
  !
  !  take one timestep of:
  !          u(:,:,:,np) = u(:,:,:,np) +  dt2*nu*laplacian**order ( u )
  !          T(:,:,:,np) = T(:,:,:,np) +  dt2*nu_s*laplacian**order ( T )
  !
  !
  !  For correct scaling, dt2 should be the same 'dt2' used in the leapfrog advace
  !
  !
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : nu, nu_div, nu_s, hypervis_order, hypervis_subcycle, nu_p, nu_top, psurf_vis
  use hybrid_mod, only : hybrid_t
  use hybvcoord_mod, only : hvcoord_t
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t
  use derivative_mod, only : laplace_sphere_wk_routine
  use derivative_mod, only : vlaplace_sphere_wk_routine
  use edge_mod, only : edgevpack, edgevunpack
  use edgetype_mod, only : EdgeBuffer_t
  use bndry_mod, only : bndry_exchangev
  use viscosity_mod, only : biharmonic_wk
  use physical_constants, only: Cp
!  use time_mod, only : TimeLevel_t
  implicit none

  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (EdgeBuffer_t)  , intent(inout) :: edge3
  type (derivative_t)  , intent(in) :: deriv
  type (hvcoord_t), intent(in)      :: hvcoord
!  type (TimeLevel_t)   , intent(in) :: tl

  real (kind=real_kind) :: dt2
  integer :: nets,nete

  ! local
  real (kind=real_kind) :: nu_scale, dpdn,dpdn0, nu_scale_top
  integer :: k,kptr,i,j,ie,ic,n0,nt,nm1
  real (kind=real_kind), dimension(np,np,2,nlev,nets:nete)      :: vtens
  real (kind=real_kind), dimension(np,np,nlev,nets:nete)        :: ptens
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(np,np,nlev) :: p
  real (kind=real_kind), dimension(np,np) :: dXdp


! NOTE: PGI compiler bug: when using spheremp, rspheremp and ps as pointers to elem(ie)% members,
  !       data is incorrect (offset by a few numbers actually)
  !       removed for now.
  !       real (kind=real_kind), dimension(:,:), pointer :: spheremp,rspheremp
  !       real (kind=real_kind), dimension(:,:,:), pointer   :: ps

  real (kind=real_kind), dimension(np,np) :: lap_p
  real (kind=real_kind), dimension(np,np,2) :: lap_v
  real (kind=real_kind) :: v1,v2,dt,heating,utens_tmp,vtens_tmp,ptens_tmp


  if (nu_s == 0 .and. nu == 0 .and. nu_p==0 ) return;
!JMD  call t_barrierf('sync_advance_hypervis_lf', hybrid%par%comm)

!JMD  print *,'advance_hypervis_lf: start of subroutine hypervis_order == ',hypervis_order
  call t_startf('advance_hypervis_lf', t_detail_low)

! for non-leapfrog,nt=n0=nmt
!
!  nm1 = tl%nm1   ! heating term uses U,V at average of nt and nm1 levels
!  n0 = tl%n0     ! timelevel used for ps scaling.  use n0 for leapfrog.
!  nt = tl%np1    ! apply viscosity to this timelevel  (np1)


  dt=dt2/hypervis_subcycle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  regular viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 1) then
     if (nu_p>0) stop 'ERROR: hypervis_order == 1 not coded for nu_p>0'
     do ic=1,hypervis_subcycle
        do ie=nets,nete

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(lap_p,lap_v)
#endif
           do k=1,nlev
              call laplace_sphere_wk_routine (elem(ie)%state%T(:,:,k,nt),  deriv,elem(ie),var_coef=.false.,laplace=lap_p)
              call vlaplace_sphere_wk_routine(elem(ie)%state%v(:,:,:,k,nt),deriv,elem(ie),var_coef=.false.,laplace=lap_v)
              ! advace in time.  (note: DSS commutes with time stepping, so we
              ! can time advance and then DSS.  this has the advantage of
              ! not letting any discontinuties accumulate in p,v via roundoff
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)*elem(ie)%spheremp(i,j)  +  dt*nu_s*lap_p(i,j)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,1)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt)*elem(ie)%spheremp(i,j) + dt*nu*lap_v(i,j,2)
                 enddo
              enddo
           enddo

           kptr=0
           call edgeVpack(edge3, elem(ie)%state%T(:,:,:,nt),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,elem(ie)%state%v(:,:,:,:,nt),2*nlev,kptr,ie)
        enddo

        call t_startf('bndry_exchangeV.edge3', t_detail_medium)
        call bndry_exchangeV(hybrid,edge3,location='prim_advance_mod:2520')
        call t_stopf('bndry_exchangeV.edge3', t_detail_medium)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, elem(ie)%state%T(:,:,:,nt), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, elem(ie)%state%v(:,:,:,:,nt), 2*nlev, kptr, ie)

           ! apply inverse mass matrix
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np
                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,nt)
                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,nt)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,nt)
                 enddo
              enddo
           enddo
        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo  ! subcycle
  endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  hyper viscosity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (hypervis_order == 2) then
     do ic=1,hypervis_subcycle
        call biharmonic_wk(elem,pstens,ptens,vtens,deriv,edge3,hybrid,nt,nets,nete,1,nlev)
        do ie=nets,nete

           nu_scale=1
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) &
!$omp private(lap_p,lap_v,nu_scale_top,dpdn,dpdn0,nu_scale,utens_tmp,vtens_tmp,ptens_tmp)
#endif
           do k=1,nlev
              ! advace in time.
              ! note: DSS commutes with time stepping, so we can time advance and then DSS.
              ! note: weak operators alreayd have mass matrix "included"

              ! add regular diffusion in top 3 layers:
              if (nu_top>0 .and. k<=3) then
                 call laplace_sphere_wk_routine (elem(ie)%state%T(:,:,k,nt),  deriv,elem(ie),var_coef=.false.,laplace=lap_p)
                 call laplace_sphere_wk_routine (elem(ie)%state%T(:,:,k,nt),  deriv,elem(ie),var_coef=.false.,laplace=lap_v)
              endif
              nu_scale_top = 1
              if (k==1) nu_scale_top=4
              if (k==2) nu_scale_top=2

              do j=1,np
                 do i=1,np
                    if (psurf_vis==0) then
                       ! normalize so as to conserve IE  (not needed when using p-surface viscosity)
                       ! scale velosity by 1/rho (normalized to be O(1))
                       ! dp/dn = O(ps0)*O(delta_eta) = O(ps0)/O(nlev)
                       dpdn = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(i,j,n0)  ! nt ?
                       dpdn0 = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                            ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*hvcoord%ps0
                       nu_scale = dpdn0/dpdn
                    endif

                    ! biharmonic terms need a negative sign:
                    if (nu_top>0 .and. k<=3) then
                       utens_tmp=(-nu*vtens(i,j,1,k,ie) + nu_scale_top*nu_top*lap_v(i,j,1))
                       vtens_tmp=(-nu*vtens(i,j,2,k,ie) + nu_scale_top*nu_top*lap_v(i,j,2))
                       ptens_tmp=nu_scale*(-nu_s*ptens(i,j,k,ie) + nu_scale_top*nu_top*lap_p(i,j) )
                    else
                       utens_tmp=-nu*vtens(i,j,1,k,ie)
                       vtens_tmp=-nu*vtens(i,j,2,k,ie)
                       ptens_tmp=-nu_scale*nu_s*ptens(i,j,k,ie)
                    endif

                    ptens(i,j,k,ie) = ptens_tmp
                    vtens(i,j,1,k,ie)=utens_tmp
                    vtens(i,j,2,k,ie)=vtens_tmp
                 enddo
              enddo
           enddo

           pstens(:,:,ie)  =  -nu_p*pstens(:,:,ie)
           kptr=0
           call edgeVpack(edge3, ptens(:,:,:,ie),nlev,kptr,ie)
           kptr=nlev
           call edgeVpack(edge3,vtens(:,:,:,:,ie),2*nlev,kptr,ie)
           kptr=3*nlev
           call edgeVpack(edge3,pstens(:,:,ie),1,kptr,ie)
        enddo


        call t_startf('bndry_exchangeV.edge3', t_detail_medium)
        call bndry_exchangeV(hybrid,edge3,location='prim_advance_mod:2618')
        call t_stopf('bndry_exchangeV.edge3', t_detail_medium)

        do ie=nets,nete

           kptr=0
           call edgeVunpack(edge3, ptens(:,:,:,ie), nlev, kptr, ie)
           kptr=nlev
           call edgeVunpack(edge3, vtens(:,:,:,:,ie), 2*nlev, kptr, ie)
           kptr=3*nlev
           call edgeVunpack(edge3, pstens(:,:,ie), 1, kptr, ie)

           if (psurf_vis == 1 ) then
              ! apply p-surface correction
              do k=1,nlev
                 p(:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps(:,:,nt)
              enddo
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(k,dXdp)
#endif
              do k=1,nlev
                 if (k.eq.1) then
                    ! no correction needed
                 else if (k.eq.nlev) then
                    ! one-sided difference
                    dXdp = (elem(ie)%state%T(:,:,k,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                        (p(:,:,k)-p(:,:,k-1))
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 else
                    dXdp = (elem(ie)%state%T(:,:,k+1,nt) - elem(ie)%state%T(:,:,k-1,nt)) / &
                         (p(:,:,k+1)-p(:,:,k-1))
                    ptens(:,:,k,ie) = ptens(:,:,k,ie) - dXdp(:,:)*hvcoord%hybm(k)*pstens(:,:,ie)
                 endif
              enddo
           endif


           ! apply inverse mass matrix, accumulate tendencies
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(v1,v2,heating)
#endif
           do k=1,nlev
              do j=1,np
                 do i=1,np

                    elem(ie)%state%v(i,j,1,k,nt)=elem(ie)%state%v(i,j,1,k,nt) + &
                         dt*elem(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                    elem(ie)%state%v(i,j,2,k,nt)=elem(ie)%state%v(i,j,2,k,nt) +  &
                         dt*elem(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)

                    ! better E conservation if we use v after adding in vtens:
                    v1=.5*(elem(ie)%state%v(i,j,1,k,nt)+elem(ie)%state%v(i,j,1,k,nm1))
                    v2=.5*(elem(ie)%state%v(i,j,2,k,nt)+elem(ie)%state%v(i,j,2,k,nm1))
                    heating = (vtens(i,j,1,k,ie)*v1  + vtens(i,j,2,k,ie)*v2 )

                    elem(ie)%state%T(i,j,k,nt)=elem(ie)%state%T(i,j,k,nt)     + &
                         dt*elem(ie)%rspheremp(i,j)*(cp*ptens(i,j,k,ie) - heating)/cp

                 enddo
              enddo
           enddo
           elem(ie)%state%ps(:,:,nt)=elem(ie)%state%ps(:,:,nt) + dt*elem(ie)%rspheremp(:,:)*pstens(:,:,ie)
        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo
  endif

  call t_stopf('advance_hypervis_lf', t_detail_low)

  end subroutine advance_hypervis_lf


  !
  ! dry-mass vertical coordinate version of compute_and_apply_rhs
  !
  subroutine compute_and_apply_rhs_dry(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,eta_ave_w)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
  !          qn0=-1 for the dry case
  !
  ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
  !
  ! Combining the RHS and DSS pack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! Combining the dt advance and DSS unpack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! note: for prescribed velocity case, velocity will be computed at
  ! "real_time", which should be the time of timelevel n0.
  !
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, nc, nlev, ntrac, max_corner_elem
  use dimensions_mod, only : qsize_condensate_loading
  use hybrid_mod, only : hybrid_t
  use element_mod, only : element_t,PrintElem
  use derivative_mod, only : derivative_t
  use derivative_mod, only : divergence_sphere_routine
  use derivative_mod, only : vorticity_sphere_routine 
  use derivative_mod, only : gradient_sphere_routine  
  use derivative_mod, only : subcell_div_fluxes, subcell_dss_fluxes
  use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
  use edgetype_mod, only : edgedescriptor_t
  use bndry_mod, only : bndry_exchangev
  use control_mod, only : moisture, qsplit, use_cpstar, rsplit, swest
  use hybvcoord_mod, only : hvcoord_t

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa, Rd_on_Rv, cp_liq, cp_ice
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv, preq_omega_ps, preq_hydrostatic
#if ( defined CAM )
  use control_mod, only: se_met_nudge_u, se_met_nudge_p, se_met_nudge_t, se_met_tevolve
#endif

  use time_mod, only : tevolve

  implicit none
  integer, intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8, intent(in) :: dt2

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi

  real (kind=real_kind), dimension(np,np,nlev)   :: omega_full
  real (kind=real_kind), dimension(np,np,nlev)   :: divdp_dry
  real (kind=real_kind), dimension(np,np,nlev)   :: divdp_full
  real (kind=real_kind), dimension(np,np,nlev+1)   :: eta_dot_dpdn  ! half level vertical velocity on p-grid
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind), dimension(np,np,2,nlev):: vdp_dry       !                            
  real (kind=real_kind), dimension(np,np,2,nlev):: vdp_full !dry+vapor
  real (kind=real_kind), dimension(np,np,nlev)   :: vgrad_p_full  
  real (kind=real_kind), dimension(np,np,2     ):: v         !                            
  real (kind=real_kind), dimension(np,np)      :: vgrad_T    ! v.grad(T)
  real (kind=real_kind), dimension(np,np)      :: Ephi       ! kinetic energy + PHI term
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p_full
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_dp_dry
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p_m_pmet  ! gradient(p - p_met)
  real (kind=real_kind), dimension(np,np,nlev)   :: vort       ! vorticity
  real (kind=real_kind), dimension(np,np,nlev)   :: p_dry      ! pressure dry
  real (kind=real_kind), dimension(np,np,nlev)   :: dp_dry     ! delta pressure dry
  real (kind=real_kind), dimension(np,np,nlev)   :: p_full      ! pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: dp_full     ! delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: p_fullmass      ! pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: dp_fullmass     ! delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: rdp        ! inverse of delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: T_vadv     ! temperature vertical advection
  real (kind=real_kind), dimension(np,np,nlev+1) :: ph               ! half level pressures on p-grid
  real (kind=real_kind), dimension(np,np,2,nlev) :: v_vadv   ! velocity vertical advection
  real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)          :: corners
  real (kind=real_kind), dimension(2,2,2)                         :: cflux
  real (kind=real_kind), dimension(np,np)        :: suml
  real (kind=real_kind), dimension(np,np,2)      :: sump
  real (kind=real_kind) ::  kappa_star(np,np,nlev)
  real (kind=real_kind) ::  waterterm, termp(2)
  real (kind=real_kind) ::  vtens1(np,np,nlev)
  real (kind=real_kind) ::  vtens2(np,np,nlev)
  real (kind=real_kind) ::  ttens(np,np,nlev)
  real (kind=real_kind) ::  stashdp3d (np,np,nlev)
  real (kind=real_kind) ::  tempdp3d  (np,np)
  real (kind=real_kind) ::  tempflux  (nc,nc,4)
  real (kind=real_kind) :: inv_epsilon, ckk, term, cp_full_inv, T_v(np,np,nlev)

  real (kind=real_kind) :: rho_inv_u_eqn(np,np), energy_conv_term_T_eqn(np,np)

  type (EdgeDescriptor_t)                                       :: desc

  real (kind=real_kind) :: density_dry, sum_water, density_inv, R_over_cp_full(np,np,nlev)

  real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
  real (kind=real_kind) ::  glnps1,glnps2,gpterm, rho_full
  integer :: i,j,k,kptr,ie,nq
  real (kind=real_kind) ::  u_m_umet, v_m_vmet, t_m_tmet 


!JMD  call t_barrierf('sync_compute_and_apply_rhs', hybrid%par%comm)
  inv_epsilon = 1/Rd_on_Rv

  call t_startf('compute_and_apply_rhs_dry', t_detail_high)
  do ie=nets,nete
     !ps => elem(ie)%state%ps(:,:,n0)
     phi => elem(ie)%derived%phi(:,:,:)

     ! ==================================================
     ! compute pressure (p) on half levels from ps
     ! using the hybrid coordinates relationship, i.e.
     ! e.g. equation (3.a.92) of the CCM-2 description,
     ! (NCAR/TN-382+STR), June 1993, p. 24.
     ! ==================================================
!     do k=1,nlev
!       do j=1,np
!         do i=1,np
!           !
!           ! use ideal gas law for dry air to calculate dry density
!           !
!           sum_water   = 1+SUM(elem(ie)%state%q(i,j,k,1:qsize_condensate_loading))
!           t_v(i,j,k)  = elem(ie)%state%T(i,j,k,n0)*(1+inv_epsilon*elem(ie)%state%q(i,j,k,1))/sum_water
!         end do
!       end do
!     end do

     ! ============================
     ! compute p and delta p
     ! ============================

     do k=1,nlev
       if (rsplit==0) then
         dp_dry(:,:,k) = (hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*elem(ie)%state%ps(:,:,n0)) &
              - (hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps(:,:,n0))
         p_dry (:,:,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps(:,:,n0)
       else
         ! vertically lagrangian code: we advect dp3d instead of ps
         ! we also need grad(p) at all levels (not just grad(ps))
         !p(k)= hyam(k)*ps0 + hybm(k)*ps
         !    = .5*(hyai(k+1)+hyai(k))*ps0 + .5*(hybi(k+1)+hybi(k))*ps
         !    = .5*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
         !
         ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
         !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
         dp_dry(:,:,k) = elem(ie)%state%dp3d(:,:,k,n0)
         if (k==1) then
           p_dry(:,:,k)=hvcoord%hyai(k)*hvcoord%ps0 + dp_dry(:,:,k)/2
         else
           p_dry(:,:,k)=p_dry(:,:,k-1) + dp_dry(:,:,k-1)/2 + dp_dry(:,:,k)/2
         endif
       end if
         do j=1,np
            do i=1,np
               !                                                                                                                               
               ! use ideal gas law for dry air to calculate dry density                                                                        
               !                                                                                                                               
               do nq=1,qsize_condensate_loading
                  elem(ie)%state%q(i,j,k,nq) =elem(ie)%state%Qdp(i,j,k,nq,qn0)/dp_dry(i,j,k)
               end do
               sum_water   = 1+SUM(elem(ie)%state%q(i,j,k,1:qsize_condensate_loading))
               t_v(i,j,k)  = elem(ie)%state%T(i,j,k,n0)*(1+inv_epsilon*elem(ie)%state%q(i,j,k,1))/sum_water
            end do
         end do

       !
       ! convert to gas pressure (dry + water vapor pressure)
       ! 
       dp_full(:,:,k)=(1+elem(ie)%state%q(:,:,k,1))*dp_dry(:,:,k)
       if (k==1) then
         p_full(:,:,k) = hvcoord%hyai(k)*hvcoord%ps0 + dp_full(:,:,k)/2
       else
         p_full(:,:,k)=p_full(:,:,k-1) + dp_full(:,:,k-1)/2 + dp_full(:,:,k)/2
       endif

       call gradient_sphere_routine(dp_dry(:,:,k),deriv,elem(ie)%Dinv,grad_dp_dry(:,:,:,k))
       call gradient_sphere_routine(p_full(:,:,k),deriv,elem(ie)%Dinv,grad_p_full(:,:,:,k))
       
       rdp(:,:,k) = 1.0D0/dp_dry(:,:,k)
       
       ! ==============================
       ! compute vgrad_lnps - for omega_full
        ! ==============================
       do j=1,np
         do i=1,np
           v1 = elem(ie)%state%v(i,j,1,k,n0)
           v2 = elem(ie)%state%v(i,j,2,k,n0)
           vgrad_p_full(i,j,k) = (v1*grad_p_full(i,j,1,k) + v2*grad_p_full(i,j,2,k))
           vdp_dry(i,j,1,k) = v1*dp_dry(i,j,k)
           vdp_dry(i,j,2,k) = v2*dp_dry(i,j,k)
           vdp_full(i,j,1,k) = v1*dp_full(i,j,k)
           vdp_full(i,j,2,k) = v2*dp_full(i,j,k)
         end do
       end do
       
#if ( defined CAM )
       ! ============================
       ! compute grad(P-P_met)
       ! ============================
       if (se_met_nudge_p.gt.0.D0) then
         call gradient_sphere_routine(elem(ie)%derived%ps_met(:,:)+tevolve*elem(ie)%derived%dpsdt_met(:,:),&
                              deriv,elem(ie)%Dinv,vtemp((:,:,:))
         grad_p_m_pmet(:,:,:,k) = grad_p_full(:,:,:,k) - hvcoord%hybm(k)*vtemp(:,:,:)
       endif
#endif
       
       ! ================================
       ! Accumulate mean Vel_rho flux in vn0
       ! ================================
       elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vdp_dry(:,:,:,k)
       
       
       ! =========================================
       !
       ! Compute relative vorticity and divergence
       !
       ! =========================================
       call divergence_sphere_routine( vdp_dry(:,:,:,k),deriv,elem(ie), divdp_dry(:,:,k))
       call divergence_sphere_routine(vdp_full(:,:,:,k),deriv,elem(ie),divdp_full(:,:,k))
       call vorticity_sphere_routine(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie),vort(:,:,k))
       
     enddo

     ! compute T_v for timelevel n0
     !
     ! ====================================================
     ! Compute Hydrostatic equation, modeld after CCM-3
     ! ====================================================
     !
     ! in the dry mass vertical coordinate case the hydrostatic relation 
     ! is dry
     !
     call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p_full,dp_full)
     
     ! ====================================================
     ! Compute omega_full
     ! ====================================================
     ckk       = 0.5d0
     suml(:,:  ) = 0
     sump(:,:,:) = 0
     do k=1,nlev
        do j=1,np 
           do i=1,np
              waterterm   = (1.0D0+elem(ie)%state%q(i,j,k,1))
              term         = -divdp_full(i,j,k)!*waterterm
              termp(:)     = grad_dp_dry(i,j,:,k)*waterterm
              
              v1 = elem(ie)%state%v(i,j,1,k,n0)
              v2 = elem(ie)%state%v(i,j,2,k,n0)
              
              omega_full(i,j,k) = suml(i,j) + ckk*term+&!v1*(ckk*termp(1)+sump(i,j,1)) + v2*(ckk*termp(2)+sump(i,j,2))
                                  vgrad_p_full(i,j,k)
              
              suml(i,j)    = suml(i,j)   + term
              sump(i,j,:)  = sump(i,j,:) + termp(:)
           end do
        end do
     end do

     ! ==================================================
     ! zero partial sum for accumulating sum
     !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
     ! used by eta_dot_dpdn and lnps tendency
     ! ==================================================
     sdot_sum=0


     ! ==================================================
     ! Compute eta_dot_dpdn
     ! save sdot_sum as this is the -RHS of ps equation
     ! ==================================================
     if (rsplit>0) then
        ! VERTICALLY LAGRANGIAN:   no vertical motion
        eta_dot_dpdn=0
        T_vadv=0
        v_vadv=0
     else
        do k=1,nlev
           ! ==================================================
           ! add this term to PS equation so we exactly conserve dry mass
           ! ==================================================
           sdot_sum(:,:) = sdot_sum(:,:) + divdp_dry(:,:,k)
           eta_dot_dpdn(:,:,k+1) = sdot_sum(:,:)
        end do


        ! ===========================================================
        ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
        ! compute at interfaces:
        !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
        ! for reference: at mid layers we have:
        !    omega = v grad p  - integral_etatop^eta[ divdp ]
        ! ===========================================================
#if (defined _LOOP_OPENMP)
        !$omp parallel do num_threads(vert_num_threads) 
#endif
        do k=1,nlev-1
           eta_dot_dpdn(:,:,k+1) = hvcoord%hybi(k+1)*sdot_sum(:,:) - eta_dot_dpdn(:,:,k+1)
        end do

        eta_dot_dpdn(:,:,1     ) = 0.0D0
        eta_dot_dpdn(:,:,nlev+1) = 0.0D0

        ! ===========================================================
        ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
        ! ==============================================
        call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
             eta_dot_dpdn,rdp,T_vadv,v_vadv)
      endif


     ! ================================
     ! accumulate mean vertical flux:
     ! ================================
#if (defined _LOOP_OPENMP)
     !$omp parallel do num_threads(vert_num_threads)
#endif
     do k=1,nlev 
        elem(ie)%derived%eta_dot_dpdn(:,:,k) = &
             elem(ie)%derived%eta_dot_dpdn(:,:,k) + eta_ave_w*eta_dot_dpdn(:,:,k)
        elem(ie)%derived%omega(:,:,k) = &
             elem(ie)%derived%omega(:,:,k) + eta_ave_w*omega_full(:,:,k)
     enddo
     elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) = &
          elem(ie)%derived%eta_dot_dpdn(:,:,nlev+1) + eta_ave_w*eta_dot_dpdn(:,:,nlev+1)
     !
     ! compute R/Cp term for thermodynamic equation
     !
!     if (qsize_condensate_loading==0) then
!        R_over_cp_full = kappa
!     else if (qsize_condensate_loading==1) then
!        do k=1,nlev
!           do j=1,np
!              do i=1,np
!                 R_over_cp_full(i,j,k)   = kappa*(1+      inv_epsilon*elem(ie)%state%q(i,j,k,1))/&
!                      (1+(8.0D0/7.0D0)*inv_epsilon*elem(ie)%state%q(i,j,k,1))
!              end do
!           end do
!        end do
!     else if (qsize_condensate_loading==3) then
!        !
!        ! This code assumes the following tracer index ordering (CAM)
!        !
!        ! 1  Q         Specific humidity
!        ! 2  CLDLIQ    Grid box averaged cloud liquid amount
!        ! 3  CLDICE    Grid box averaged cloud ice amount
!        !
!        do k=1,nlev
!           do j=1,np
!              do i=1,np
!                 R_over_cp_full(i,j,k)   = kappa*(1+      inv_epsilon*elem(ie)%state%q(i,j,k,1))/&
!                      (1+(8.0D0/7.0D0)*inv_epsilon*elem(ie)%state%q(i,j,k,1)+&
!                      (cp_liq/cp)*elem(ie)%state%q(i,j,k,2)+ &
!                      (cp_ice/cp)*elem(ie)%state%q(i,j,k,3))
!              end do
!           end do
!        end do
!     else
!        write(*,*) "condensate loading option not supported"
!        stop
!     endif
!     R_over_cp_full = kappa!phlxxxxx

     ! ==============================================
     ! Compute phi + kinetic energy term: 10*nv*nv Flops
     ! ==============================================
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) &
!$omp private(v1,v2,E,Ephi,vtemp,vgrad_T,gpterm,glnps1,glnps2,sum_water) &
!$omp private(density_dry,density_inv,u_m_umet,v_m_vmet,t_m_tmet)
#endif
     vertloop: do k=1,nlev
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              E = 0.5D0*( v1*v1 + v2*v2 )
              Ephi(i,j)=E+phi(i,j,k)!+elem(ie)%derived%pecnd(i,j,k)
           end do
        end do
        ! ================================================
        ! compute gradp term (ps/p)*(dp/dps)*T
        ! ================================================
        call gradient_sphere_routine(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv,vtemp(:,:,:))
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
           end do
        end do


        ! vtemp = grad ( E + PHI )
        call gradient_sphere_routine(Ephi(:,:),deriv,elem(ie)%Dinv,vtemp) 

        do j=1,np
           do i=1,np
             !
             ! use ideal gas law for dry air to calculate dry density
             !
             sum_water   = 1+SUM(elem(ie)%state%q(i,j,k,1:qsize_condensate_loading))
             density_dry = p_dry(i,j,k)/(Rgas*elem(ie)%state%T(i,j,k,n0))
!             density_inv = 1/(density_dry*sum_water)
             density_inv = Rgas*T_v(i,j,k)/p_full(i,j,k)

             glnps1  = density_inv*grad_p_full(i,j,1,k)
             glnps2  = density_inv*grad_p_full(i,j,2,k)

             v1     = elem(ie)%state%v(i,j,1,k,n0)
             v2     = elem(ie)%state%v(i,j,2,k,n0)
             
             vtens1(i,j,k) =   - v_vadv(i,j,1,k)                           &
                  + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                  - vtemp(i,j,1) - glnps1

             vtens2(i,j,k) =   - v_vadv(i,j,2,k)                            &
                  - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                  - vtemp(i,j,2) - glnps2
             !
             ! for now cp only includes water vapor
             !
             ttens(i,j,k)  = - T_vadv(i,j,k) - vgrad_T(i,j) + &
                               density_inv*omega_full(i,j,k)/cp
!                               R_over_cp_full(i,j,k)*elem(ie)%state%T(i,j,k,n0)*omega_full(i,j,k)/p_full(i,j,k)

             !
             ! phl: add forcing term to T
             !
             
#if ( defined CAM )
              
              if (se_prescribed_wind_2d) then
                 vtens1(i,j,k) = 0.D0
                 vtens2(i,j,k) = 0.D0
                 ttens(i,j,k) = 0.D0
              else
                 if(se_met_nudge_u.gt.0.D0)then
                    u_m_umet = v1 - &
                         elem(ie)%derived%u_met(i,j,k) - &
                         se_met_tevolve*tevolve*elem(ie)%derived%dudt_met(i,j,k)
                    v_m_vmet = v2 - &
                         elem(ie)%derived%v_met(i,j,k) - &
                         se_met_tevolve*tevolve*elem(ie)%derived%dvdt_met(i,j,k)

                    vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)

                    elem(ie)%derived%Utnd(i+(j-1)*np,k) = elem(ie)%derived%Utnd(i+(j-1)*np,k) &
                         + se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)

                    vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)

                    elem(ie)%derived%Vtnd(i+(j-1)*np,k) = elem(ie)%derived%Vtnd(i+(j-1)*np,k) &
                         + se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)

                 endif

                 if(se_met_nudge_p.gt.0.D0)then
                    vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,1,k)  * elem(ie)%derived%nudge_factor(i,j,k)
                    vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,2,k)  * elem(ie)%derived%nudge_factor(i,j,k)
                 endif

                 if(se_met_nudge_t.gt.0.D0)then
                    t_m_tmet = elem(ie)%state%T(i,j,k,n0) - &
                         elem(ie)%derived%T_met(i,j,k) - &
                         se_met_tevolve*tevolve*elem(ie)%derived%dTdt_met(i,j,k)
                    ttens(i,j,k)  = ttens(i,j,k) - se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                    elem(ie)%derived%Ttnd(i+(j-1)*np,k) = elem(ie)%derived%Ttnd(i+(j-1)*np,k) &
                         + se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                 endif
              endif
#endif

           end do
        end do

     end do vertloop

     ! =========================================================
     ! local element timestep, store in np1.
     ! note that we allow np1=n0 or nm1
     ! apply mass matrix
     ! =========================================================
     if (dt2<0) then
        ! calling program just wanted DSS'd RHS, skip time advance
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(v,tempflux)
#endif
        do k=1,nlev
           elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*vtens1(:,:,k)
           elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*vtens2(:,:,k)
           elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremp(:,:)*ttens(:,:,k)
           if (rsplit>0) &
              elem(ie)%state%dp3d(:,:,k,np1) = -elem(ie)%spheremp(:,:)*&
              (divdp_dry(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k))
           if (0<rsplit.and.0<ntrac.and.eta_ave_w.ne.0.) then
              v(:,:,1) =  elem(ie)%Dinv(:,:,1,1)*vdp_dry(:,:,1,k) + elem(ie)%Dinv(:,:,1,2)*vdp_dry(:,:,2,k)
              v(:,:,2) =  elem(ie)%Dinv(:,:,2,1)*vdp_dry(:,:,1,k) + elem(ie)%Dinv(:,:,2,2)*vdp_dry(:,:,2,k)
              tempflux =  eta_ave_w*subcell_div_fluxes(v, np, nc, elem(ie)%metdet)
              elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) - tempflux
           end if
        enddo
        elem(ie)%state%ps(:,:,np1) = -elem(ie)%spheremp(:,:)*sdot_sum
     else
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(v,tempflux) 
#endif
        do k=1,nlev
           elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,1,k,nm1) + dt2*vtens1(:,:,k) )
           elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%v(:,:,2,k,nm1) + dt2*vtens2(:,:,k) )
           elem(ie)%state%T(:,:,k,np1) = elem(ie)%spheremp(:,:)*(elem(ie)%state%T(:,:,k,nm1) + dt2*ttens(:,:,k))
           if (rsplit>0) &
                elem(ie)%state%dp3d(:,:,k,np1) = &
                  elem(ie)%spheremp(:,:) * (elem(ie)%state%dp3d(:,:,k,nm1) - &
                  dt2 * (divdp_dry(:,:,k) + eta_dot_dpdn(:,:,k+1)-eta_dot_dpdn(:,:,k)))


           if (0<rsplit.and.0<ntrac.and.eta_ave_w.ne.0.) then
              v(:,:,1) =  elem(ie)%Dinv(:,:,1,1)*vdp_dry(:,:,1,k) + elem(ie)%Dinv(:,:,1,2)*vdp_dry(:,:,2,k)
              v(:,:,2) =  elem(ie)%Dinv(:,:,2,1)*vdp_dry(:,:,1,k) + elem(ie)%Dinv(:,:,2,2)*vdp_dry(:,:,2,k)
              tempflux =  eta_ave_w*subcell_div_fluxes(v, np, nc, elem(ie)%metdet)
              elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) - tempflux
           end if
        enddo
        elem(ie)%state%ps(:,:,np1) = elem(ie)%spheremp(:,:)*( elem(ie)%state%ps(:,:,nm1) - dt2*sdot_sum )

     endif


     ! =========================================================
     !
     ! Pack ps(np1), T, and v tendencies into comm buffer
     !
     ! =========================================================
     kptr=0
     call edgeVpack(edge3p1, elem(ie)%state%ps(:,:,np1),1,kptr,ie)

     kptr=1
     call edgeVpack(edge3p1, elem(ie)%state%T(:,:,:,np1),nlev,kptr,ie)

     kptr=nlev+1
     call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)

     if (rsplit>0) then
        kptr=kptr+2*nlev
        call edgeVpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     endif
  end do

  ! =============================================================
    ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================
  call t_startf('bndry_exchangeV.edge3p1', t_detail_medium)
  call bndry_exchangeV(hybrid,edge3p1,location='prim_advance_mod:3270')
  call t_stopf('bndry_exchangeV.edge3p1', t_detail_medium)
  do ie=nets,nete
     ! ===========================================================
     ! Unpack the edges for vgrad_T and v tendencies...
     ! ===========================================================
     kptr=0
     call edgeVunpack(edge3p1, elem(ie)%state%ps(:,:,np1), 1, kptr, ie)

     kptr=1
     call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,:,np1), nlev, kptr, ie)

     kptr=nlev+1
     call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

     if (rsplit>0) then
        if (0<ntrac.and.eta_ave_w.ne.0.) then
          do k=1,nlev
             stashdp3d(:,:,k) = elem(ie)%state%dp3d(:,:,k,np1)/elem(ie)%spheremp(:,:)
          end do
        endif

        corners = 0.0d0
        corners(1:np,1:np,:) = elem(ie)%state%dp3d(:,:,:,np1)
        kptr=kptr+2*nlev
        call edgeVunpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,ie)

        if  (0<ntrac.and.eta_ave_w.ne.0.) then
          desc = elem(ie)%desc
          call edgeDGVunpack(edge3p1, corners, nlev, kptr, ie)
          corners = corners/dt2

          do k=1,nlev
            tempdp3d = elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
            tempdp3d = tempdp3d - stashdp3d(:,:,k)
            tempdp3d = tempdp3d/dt2

            call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)
 
            cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)  
            cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:) 
            cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:) 
            cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:) 

            tempflux =  eta_ave_w*subcell_dss_fluxes(tempdp3d, np, nc, elem(ie)%metdet, cflux)
            elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + tempflux
          end do
        end if   
     endif

     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) 
#endif
     do k=1,nlev
        elem(ie)%state%T(:,:,k,np1)   = elem(ie)%rspheremp(:,:)*elem(ie)%state%T(:,:,k,np1)
        elem(ie)%state%v(:,:,1,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,1,k,np1)
        elem(ie)%state%v(:,:,2,k,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%v(:,:,2,k,np1)
     end do

     if (rsplit>0) then
        ! vertically lagrangian: complete dp3d timestep:
        do k=1,nlev
           elem(ie)%state%dp3d(:,:,k,np1)= elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
        enddo
        ! when debugging: also update ps
        !elem(ie)%state%ps(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps(:,:,np1)
     else
        ! vertically eulerian: complete ps timestep:
        elem(ie)%state%ps(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps(:,:,np1)
     endif

  end do

#ifdef DEBUGOMP
!$OMP BARRIER
#endif
  call t_stopf('compute_and_apply_rhs_dry', t_detail_low)
end subroutine compute_and_apply_rhs_dry

  subroutine compute_and_apply_rhs(np1,nm1,n0,qn0,dt2,elem,hvcoord,hybrid,&
       deriv,nets,nete,eta_ave_w)
  ! ===================================
  ! compute the RHS, accumulate into u(np1) and apply DSS
  !
  !           u(np1) = u(nm1) + dt2*DSS[ RHS(u(n0)) ]
  !
  ! This subroutine is normally called to compute a leapfrog timestep
  ! but by adjusting np1,nm1,n0 and dt2, many other timesteps can be
  ! accomodated.  For example, setting nm1=np1=n0 this routine will
  ! take a forward euler step, overwriting the input with the output.
  !
  !    qn0 = timelevel used to access Qdp() in order to compute virtual Temperature
  !          qn0=-1 for the dry case
  !
  ! if  dt2<0, then the DSS'd RHS is returned in timelevel np1
  !
  ! Combining the RHS and DSS pack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! Combining the dt advance and DSS unpack operation in one routine
  ! allows us to fuse these two loops for more cache reuse
  !
  ! note: for prescribed velocity case, velocity will be computed at
  ! "real_time", which should be the time of timelevel n0.
  !
  !
  ! ===================================
  use kinds, only : real_kind
  use dimensions_mod, only : np, nc, nlev, ntrac, max_corner_elem
  use hybrid_mod, only : hybrid_t, get_loop_ranges
  use element_mod, only : element_t,PrintElem
  use derivative_mod, only : derivative_t
  use derivative_mod, only : divergence_sphere_routine
  use derivative_mod, only : vorticity_sphere_routine 
  use derivative_mod, only : gradient_sphere_routine 
  use derivative_mod, only : subcell_div_fluxes, subcell_dss_fluxes
  use edge_mod, only : edgevpack, edgevunpack, edgeDGVunpack
  use edgetype_mod, only : edgedescriptor_t
  use bndry_mod, only : bndry_exchangev
  use control_mod, only : moisture, qsplit, use_cpstar, rsplit, swest
  use hybvcoord_mod, only : hvcoord_t

  use physical_constants, only : cp, cpwater_vapor, Rgas, kappa
  use physics_mod, only : virtual_specific_heat, virtual_temperature
  use prim_si_mod, only : preq_vertadv, preq_omega_ps, preq_hydrostatic
#if ( defined CAM )
  use control_mod, only: se_met_nudge_u, se_met_nudge_p, se_met_nudge_t, se_met_tevolve
#endif

  use time_mod, only : tevolve

  implicit none
  integer, intent(in) :: np1,nm1,n0,qn0,nets,nete
  real*8, intent(in) :: dt2

  type (hvcoord_t)     , intent(in) :: hvcoord
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind) :: eta_ave_w  ! weighting for eta_dot_dpdn mean flux

  ! local
  real (kind=real_kind), pointer, dimension(:,:)      :: ps         ! surface pressure for current tiime level
  real (kind=real_kind), pointer, dimension(:,:,:)   :: phi

  real (kind=real_kind), dimension(np,np,nlev)   :: omega_p
  real (kind=real_kind), dimension(np,np,nlev)   :: T_v
  real (kind=real_kind), dimension(np,np,nlev)   :: divdp
  real (kind=real_kind), dimension(np,np,nlev+1)   :: eta_dot_dpdn  ! half level vertical velocity on p-grid
  real (kind=real_kind), dimension(np,np)      :: sdot_sum   ! temporary field
  real (kind=real_kind), dimension(np,np,2)    :: vtemp     ! generic gradient storage
  real (kind=real_kind), dimension(np,np,2,nlev):: vdp       !                            
  real (kind=real_kind), dimension(np,np,2     ):: v         !                            
  real (kind=real_kind), dimension(np,np)      :: vgrad_T    ! v.grad(T)
  real (kind=real_kind), dimension(np,np)      :: Ephi       ! kinetic energy + PHI term
  real (kind=real_kind), dimension(np,np,2)      :: grad_ps    ! lat-lon coord version
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p
  real (kind=real_kind), dimension(np,np,2,nlev) :: grad_p_m_pmet  ! gradient(p - p_met)
  real (kind=real_kind), dimension(np,np,nlev)   :: vort       ! vorticity
  real (kind=real_kind), dimension(np,np,nlev)   :: p          ! pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: dp         ! delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: rdp        ! inverse of delta pressure
  real (kind=real_kind), dimension(np,np,nlev)   :: T_vadv     ! temperature vertical advection
  real (kind=real_kind), dimension(np,np,nlev)   :: vgrad_p    ! v.grad(p)
  real (kind=real_kind), dimension(np,np,nlev+1) :: ph               ! half level pressures on p-grid
  real (kind=real_kind), dimension(np,np,2,nlev) :: v_vadv   ! velocity vertical advection
  real (kind=real_kind), dimension(0:np+1,0:np+1,nlev)          :: corners
  real (kind=real_kind), dimension(2,2,2)                         :: cflux
  real (kind=real_kind) ::  kappa_star(np,np,nlev)
  real (kind=real_kind) ::  vtens1(np,np,nlev)
  real (kind=real_kind) ::  vtens2(np,np,nlev)
  real (kind=real_kind) ::  ttens(np,np,nlev)
  real (kind=real_kind) ::  stashdp3d (np,np,nlev)
  real (kind=real_kind) ::  tempdp3d  (np,np)
  real (kind=real_kind) ::  tempflux  (nc,nc,4)
  type (EdgeDescriptor_t)                                       :: desc

  real (kind=real_kind) ::  cp2,cp_ratio,E,de,Qt,v1,v2
  real (kind=real_kind) ::  glnps1,glnps2,gpterm
  integer :: i,j,k,kptr,ie
  real (kind=real_kind) ::  u_m_umet, v_m_vmet, t_m_tmet 

!JMD  call t_barrierf('sync_compute_and_apply_rhs', hybrid%par%comm)

  call t_startf('compute_and_apply_rhs', t_detail_low)

  do ie=nets,nete
     !ps => elem(ie)%state%ps(:,:,n0)
     phi => elem(ie)%derived%phi(:,:,:)

     ! ==================================================
     ! compute pressure (p) on half levels from ps
     ! using the hybrid coordinates relationship, i.e.
     ! e.g. equation (3.a.92) of the CCM-2 description,
     ! (NCAR/TN-382+STR), June 1993, p. 24.
     ! ==================================================
     ! vertically eulerian only needs grad(ps)
     if (rsplit==0) &
          call gradient_sphere_routine(elem(ie)%state%ps(:,:,n0),deriv,elem(ie)%Dinv,grad_ps)

     ! ============================
     ! compute p and delta p
     ! ============================
     do k=1,nlev
        if (rsplit==0) then
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
                    dp(i,j,k) = (hvcoord%hyai(k+1)*hvcoord%ps0 + hvcoord%hybi(k+1)*elem(ie)%state%ps(i,j,n0)) &
                              - (hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps(i,j,n0))
                     p(i,j,k) = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps(i,j,n0)
              grad_p(i,j,1,k) = hvcoord%hybm(k)*grad_ps(i,j,1)
              grad_p(i,j,2,k) = hvcoord%hybm(k)*grad_ps(i,j,2)
            enddo
          enddo
        else
           ! vertically lagrangian code: we advect dp3d instead of ps_v
           ! we also need grad(p) at all levels (not just grad(ps))
           !p(k)= hyam(k)*ps0 + hybm(k)*ps
           !    = .5*(hyai(k+1)+hyai(k))*ps0 + .5*(hybi(k+1)+hybi(k))*ps
           !    = .5*(ph(k+1) + ph(k) )  = ph(k) + dp(k)/2
           !
           ! p(k+1)-p(k) = ph(k+1)-ph(k) + (dp(k+1)-dp(k))/2
           !             = dp(k) + (dp(k+1)-dp(k))/2 = (dp(k+1)+dp(k))/2
           dp(:,:,k) = elem(ie)%state%dp3d(:,:,k,n0)
           if (k==1) then
              p(:,:,k)=hvcoord%hyai(k)*hvcoord%ps0 + dp(:,:,k)/2
           else
              p(:,:,k)=p(:,:,k-1) + dp(:,:,k-1)/2 + dp(:,:,k)/2
           endif
           call gradient_sphere_routine(p(:,:,k),deriv,elem(ie)%Dinv,grad_p(:,:,:,k))
        endif

        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
        do j=1,np
          do i=1,np
            rdp(i,j,k) = 1.0D0/dp(i,j,k)
          enddo
        enddo

        ! ============================
        ! compute vgrad_lnps
        ! ============================
        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
        do j=1,np
           do i=1,np
              v1 = elem(ie)%state%v(i,j,1,k,n0)
              v2 = elem(ie)%state%v(i,j,2,k,n0)
!              vgrad_p(i,j,k) = &
!                   hvcoord%hybm(k)*(v1*grad_ps(i,j,1) + v2*grad_ps(i,j,2))
              vgrad_p(i,j,k) = (v1*grad_p(i,j,1,k) + v2*grad_p(i,j,2,k))
              vdp(i,j,1,k) = v1*dp(i,j,k)
              vdp(i,j,2,k) = v2*dp(i,j,k)
           end do
        end do

#if ( defined CAM )
        ! ============================
        ! compute grad(P-P_met)
        ! ============================
        if (se_met_nudge_p.gt.0.D0) then
           call gradient_sphere_routine(elem(ie)%derived%ps_met(:,:)+tevolve*elem(ie)%derived%dpsdt_met(:,:), &
                                deriv,elem(ie)%Dinv,vtemp(:,:,:))
           grad_p_m_pmet(:,:,:,k) = grad_p(:,:,:,k) - hvcoord%hybm(k)*vtemp(:,:,:)
        endif
#endif

        ! ================================
        ! Accumulate mean Vel_rho flux in vn0
        ! ================================
!       elem(ie)%derived%vn0(:,:,:,k)=elem(ie)%derived%vn0(:,:,:,k)+eta_ave_w*vdp(:,:,:,k)
        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
        do j=1,np
          do i=1,np
            elem(ie)%derived%vn0(i,j,1,k)=elem(ie)%derived%vn0(i,j,1,k)+eta_ave_w*vdp(i,j,1,k)
            elem(ie)%derived%vn0(i,j,2,k)=elem(ie)%derived%vn0(i,j,2,k)+eta_ave_w*vdp(i,j,2,k)
          enddo
        enddo

        ! =========================================
        !
        ! Compute relative vorticity and divergence
        !
        ! =========================================
        call divergence_sphere_routine(vdp(:,:,:,k),deriv,elem(ie),divdp(:,:,k))
        call vorticity_sphere_routine(elem(ie)%state%v(:,:,:,k,n0),deriv,elem(ie),vort(:,:,k))

     enddo

     ! compute T_v for timelevel n0
     !if ( moisture /= "dry") then
     if (qn0 == -1 ) then

#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads)
#endif
        do k=1,nlev
           do j=1,np
              do i=1,np
                 T_v(i,j,k) = elem(ie)%state%T(i,j,k,n0)
                 kappa_star(i,j,k) = kappa
              end do
           end do
        end do
     else
 
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(Qt)
#endif
        do k=1,nlev
           do j=1,np
              do i=1,np
               ! Qt = elem(ie)%state%Q(i,j,k,1)
                 Qt = elem(ie)%state%Qdp(i,j,k,1,qn0)/dp(i,j,k)
!!XXgoldyXX
!Qt=0._real_kind
!!XXgoldyXX
                 T_v(i,j,k) = Virtual_Temperature(elem(ie)%state%T(i,j,k,n0),Qt)
                 if (use_cpstar==1) then
                    kappa_star(i,j,k) =  Rgas/Virtual_Specific_Heat(Qt)
                 else
                    kappa_star(i,j,k) = kappa
                 endif
              end do
           end do
        end do

     end if

     ! ====================================================
     ! Compute Hydrostatic equation, modeld after CCM-3
     ! ====================================================
     !call geopotential_t(p,dp,T_v,Rgas,phi)
     call preq_hydrostatic(phi,elem(ie)%state%phis,T_v,p,dp)

     ! ====================================================
     ! Compute omega_p according to CCM-3
     ! ====================================================
     call preq_omega_ps(omega_p,hvcoord,p,vgrad_p,divdp)

     ! ==================================================
     ! zero partial sum for accumulating sum
     !    (div(v_k) + v_k.grad(lnps))*dsigma_k = div( v dp )
     ! used by eta_dot_dpdn and lnps tendency
     ! ==================================================
     sdot_sum=0

     ! ==================================================
     ! Compute eta_dot_dpdn
     ! save sdot_sum as this is the -RHS of ps_v equation
     ! ==================================================
     if (rsplit>0) then
        ! VERTICALLY LAGRANGIAN:   no vertical motion
        eta_dot_dpdn=0
        do k=1,nlev
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              T_vadv(i,j,k)  =0
              v_vadv(i,j,1,k)=0
              v_vadv(i,j,2,k)=0
            enddo
          enddo
        enddo
     else
        do k=1,nlev
           ! ==================================================
           ! add this term to PS equation so we exactly conserve dry mass
           ! ==================================================
           do j=1,np
             do i=1,np
                   sdot_sum(i,j)     = sdot_sum(i,j) + divdp(i,j,k)
               eta_dot_dpdn(i,j,k+1) = sdot_sum(i,j)
             enddo
           enddo
        end do

        ! ===========================================================
        ! at this point, eta_dot_dpdn contains integral_etatop^eta[ divdp ]
        ! compute at interfaces:
        !    eta_dot_dpdn = -dp/dt - integral_etatop^eta[ divdp ]
        ! for reference: at mid layers we have:
        !    omega = v grad p  - integral_etatop^eta[ divdp ]
        ! ===========================================================
        do k=1,nlev-1
          !OMP_COLLAPSE_SIMD 
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
              eta_dot_dpdn(i,j,k+1) = hvcoord%hybi(k+1)*sdot_sum(i,j) - eta_dot_dpdn(i,j,k+1)
            enddo
          enddo
        end do

        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
        do j=1,np
          do i=1,np
            eta_dot_dpdn(i,j,1     ) = 0.0D0
            eta_dot_dpdn(i,j,nlev+1) = 0.0D0
          enddo
        enddo

        ! ===========================================================
        ! Compute vertical advection of T and v from eq. CCM2 (3.b.1)
        ! ==============================================
        call preq_vertadv(elem(ie)%state%T(:,:,:,n0),elem(ie)%state%v(:,:,:,:,n0), &
             eta_dot_dpdn,rdp,T_vadv,v_vadv)
     endif

     ! ================================
     ! accumulate mean vertical flux:
     ! ================================
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads)
#endif
     do k=1,nlev 
       !OMP_COLLAPSE_SIMD 
       !DIR_VECTOR_ALIGNED
       do j=1,np
         do i=1,np
           elem(ie)%derived%eta_dot_dpdn(i,j,k) = &
                  elem(ie)%derived%eta_dot_dpdn(i,j,k) + eta_ave_w*eta_dot_dpdn(i,j,k)
           elem(ie)%derived%omega(i,j,k) = &
                  elem(ie)%derived%omega(i,j,k) + eta_ave_w*omega_p(i,j,k)
         enddo
       enddo
     enddo

     !OMP_COLLAPSE_SIMD 
     !DIR_VECTOR_ALIGNED
     do j=1,np
       do i=1,np
         elem(ie)%derived%eta_dot_dpdn(i,j,nlev+1) = &
           elem(ie)%derived%eta_dot_dpdn(i,j,nlev+1) + eta_ave_w*eta_dot_dpdn(i,j,nlev+1)
       enddo
     enddo

     ! ============================================== 
     ! Compute phi + kinetic energy term: 10*nv*nv Flops
     ! ==============================================
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) &
!$omp private(v1,v2,E,Ephi,vtemp,vgrad_T,gpterm,glnps1,glnps2,u_m_umet,v_m_vmet,t_m_tmet) 
#endif
     do k=1,nlev
! Don't use COLLAPSE directive here as it causes issues with the Intel17-patched compiler
        do j=1,np
           do i=1,np
             Ephi(i,j)= 0.5D0*( elem(ie)%state%v(i,j,1,k,n0)*elem(ie)%state%v(i,j,1,k,n0) + &
                                elem(ie)%state%v(i,j,2,k,n0)*elem(ie)%state%v(i,j,2,k,n0)) + &
                                phi(i,j,k)+elem(ie)%derived%pecnd(i,j,k)
           end do
        end do

        ! ================================================
        ! compute gradp term (ps/p)*(dp/dps)*T
        ! ================================================
        call gradient_sphere_routine(elem(ie)%state%T(:,:,k,n0),deriv,elem(ie)%Dinv,vtemp(:,:,:))
        !OMP_COLLAPSE_SIMD 
        !DIR_VECTOR_ALIGNED
        do j=1,np
           do i=1,np
              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)
              vgrad_T(i,j) =  v1*vtemp(i,j,1) + v2*vtemp(i,j,2)
           end do
        end do

        ! vtemp = grad ( E + PHI )
        call gradient_sphere_routine(Ephi(:,:),deriv,elem(ie)%Dinv,vtemp(:,:,:))
        do j=1,np
           do i=1,np
!              gpterm = hvcoord%hybm(k)*T_v(i,j,k)/p(i,j,k)
!              glnps1 = Rgas*gpterm*grad_ps(i,j,1)
!              glnps2 = Rgas*gpterm*grad_ps(i,j,2)
              gpterm = T_v(i,j,k)/p(i,j,k)
              glnps1 = Rgas*gpterm*grad_p(i,j,1,k)
              glnps2 = Rgas*gpterm*grad_p(i,j,2,k)

              v1     = elem(ie)%state%v(i,j,1,k,n0)
              v2     = elem(ie)%state%v(i,j,2,k,n0)

              vtens1(i,j,k) =   - v_vadv(i,j,1,k)                           &
                   + v2*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,1) - glnps1
              !
              ! phl: add forcing term to zonal wind u
              !
              vtens2(i,j,k) =   - v_vadv(i,j,2,k)                            &
                   - v1*(elem(ie)%fcor(i,j) + vort(i,j,k))        &
                   - vtemp(i,j,2) - glnps2
              !
              ! phl: add forcing term to meridional wind v
              !
              ttens(i,j,k)  = - T_vadv(i,j,k) - vgrad_T(i,j) + kappa_star(i,j,k)*T_v(i,j,k)*omega_p(i,j,k)
              !
              ! phl: add forcing term to T
              !
#if ( defined CAM )
              
              if (se_prescribed_wind_2d) then
                 vtens1(i,j,k) = 0.D0
                 vtens2(i,j,k) = 0.D0
                 ttens(i,j,k) = 0.D0
              else
                 if(se_met_nudge_u.gt.0.D0)then
                    u_m_umet = v1 - &
                         elem(ie)%derived%u_met(i,j,k) - &
                         se_met_tevolve*tevolve*elem(ie)%derived%dudt_met(i,j,k)
                    v_m_vmet = v2 - &
                         elem(ie)%derived%v_met(i,j,k) - &
                         se_met_tevolve*tevolve*elem(ie)%derived%dvdt_met(i,j,k)

                    vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)

                    elem(ie)%derived%Utnd(i+(j-1)*np,k) = elem(ie)%derived%Utnd(i+(j-1)*np,k) &
                         + se_met_nudge_u*u_m_umet * elem(ie)%derived%nudge_factor(i,j,k)

                    vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)

                    elem(ie)%derived%Vtnd(i+(j-1)*np,k) = elem(ie)%derived%Vtnd(i+(j-1)*np,k) &
                         + se_met_nudge_u*v_m_vmet * elem(ie)%derived%nudge_factor(i,j,k)

                 endif

                 if(se_met_nudge_p.gt.0.D0)then
                    vtens1(i,j,k) =   vtens1(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,1,k)  * elem(ie)%derived%nudge_factor(i,j,k)
                    vtens2(i,j,k) =   vtens2(i,j,k) - se_met_nudge_p*grad_p_m_pmet(i,j,2,k)  * elem(ie)%derived%nudge_factor(i,j,k)
                 endif

                 if(se_met_nudge_t.gt.0.D0)then
                    t_m_tmet = elem(ie)%state%T(i,j,k,n0) - &
                         elem(ie)%derived%T_met(i,j,k) - &
                         se_met_tevolve*tevolve*elem(ie)%derived%dTdt_met(i,j,k)
                    ttens(i,j,k)  = ttens(i,j,k) - se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                    elem(ie)%derived%Ttnd(i+(j-1)*np,k) = elem(ie)%derived%Ttnd(i+(j-1)*np,k) &
                         + se_met_nudge_t*t_m_tmet * elem(ie)%derived%nudge_factor(i,j,k)
                 endif
              endif
#endif
           end do
        end do
     end do 

     ! =========================================================
     ! local element timestep, store in np1.
     ! note that we allow np1=n0 or nm1
     ! apply mass matrix
     ! =========================================================
     if (dt2<0) then
        ! calling program just wanted DSS'd RHS, skip time advance
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(v,tempflux)
#endif
        do k=1,nlev
          !OMP_COLLAPSE_SIMD
          !DIR_VECTOR_ALIGNED
          do j=1,np
            do i=1,np
               elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%spheremp(i,j)*vtens1(i,j,k)
               elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%spheremp(i,j)*vtens2(i,j,k)
                 elem(ie)%state%T(i,j,k,np1) = elem(ie)%spheremp(i,j)*ttens(i,j,k)
             enddo
           enddo
           if (rsplit>0) then 
             !OMP_COLLAPSE_SIMD 
             !DIR_VECTOR_ALIGNED
             do j=1,np
               do i=1,np
                 elem(ie)%state%dp3d(i,j,k,np1) = -elem(ie)%spheremp(i,j)*&
                   (divdp(i,j,k) + eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k))
                enddo
             enddo
           endif
           if (0<rsplit.and.0<ntrac.and.eta_ave_w.ne.0.) then 
             !OMP_COLLAPSE_SIMD 
             !DIR_VECTOR_ALIGNED
             do j=1,np
               do i=1,np
                 v(i,j,1) =  elem(ie)%Dinv(i,j,1,1)*vdp(i,j,1,k) + elem(ie)%Dinv(i,j,1,2)*vdp(i,j,2,k)
                 v(i,j,2) =  elem(ie)%Dinv(i,j,2,1)*vdp(i,j,1,k) + elem(ie)%Dinv(i,j,2,2)*vdp(i,j,2,k)
               enddo
             enddo
             tempflux =  eta_ave_w*subcell_div_fluxes(v, np, nc, elem(ie)%metdet)
             elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) - tempflux
           end if
        enddo
        !OMP_COLLAPSE_SIMD
        !DIR_VECTOR_ALIGNED
        do j=1,np
          do i=1,np
            elem(ie)%state%ps(i,j,np1) = -elem(ie)%spheremp(i,j)*sdot_sum(i,j)
          enddo
        enddo
    else
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads) private(v,tempflux)
#endif
        do k=1,nlev
           !OMP_COLLAPSE_SIMD
           !DIR_VECTOR_ALIGNED
           do j=1,np
            do i=1,np
              elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v(i,j,1,k,nm1) + dt2*vtens1(i,j,k) )
              elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%v(i,j,2,k,nm1) + dt2*vtens2(i,j,k) )
              elem(ie)%state%T(i,j,k,np1) = elem(ie)%spheremp(i,j)*(elem(ie)%state%T(i,j,k,nm1) + dt2*ttens(i,j,k))
            enddo
           enddo

           if (rsplit>0) then
             !OMP_COLLAPSE_SIMD
             !DIR_VECTOR_ALIGNED
             do j=1,np
               do i=1,np
                 elem(ie)%state%dp3d(i,j,k,np1) = &
                   elem(ie)%spheremp(i,j) * (elem(ie)%state%dp3d(i,j,k,nm1) - &
                   dt2 * (divdp(i,j,k) + eta_dot_dpdn(i,j,k+1)-eta_dot_dpdn(i,j,k)))
               enddo
             enddo
           endif


           if (0<rsplit.and.0<ntrac.and.eta_ave_w.ne.0.) then
             !OMP_COLLAPSE_SIMD
             !DIR_VECTOR_ALIGNED
             do j=1,np
               do i=1,np
                 v(i,j,1) =  elem(ie)%Dinv(i,j,1,1)*vdp(i,j,1,k) + elem(ie)%Dinv(i,j,1,2)*vdp(i,j,2,k)
                 v(i,j,2) =  elem(ie)%Dinv(i,j,2,1)*vdp(i,j,1,k) + elem(ie)%Dinv(i,j,2,2)*vdp(i,j,2,k)
               enddo
             enddo
             tempflux =  eta_ave_w*subcell_div_fluxes(v, np, nc, elem(ie)%metdet)
             elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) - tempflux
           end if
        enddo
        !OMP_COLLAPSE_SIMD
        !DIR_VECTOR_ALIGNED
        do j=1,np
          do i=1,np
            elem(ie)%state%ps(i,j,np1) = elem(ie)%spheremp(i,j)*( elem(ie)%state%ps(i,j,nm1) - dt2*sdot_sum(i,j) )
          enddo
        enddo

     endif


     ! =========================================================
     !
     ! Pack ps(np1), T, and v tendencies into comm buffer
     !
     ! =========================================================
     kptr=0
     call edgeVpack(edge3p1, elem(ie)%state%ps(:,:,np1),1,kptr,ie)

     kptr=1
     call edgeVpack(edge3p1, elem(ie)%state%T(:,:,:,np1),nlev,kptr,ie)

     kptr=nlev+1
     call edgeVpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1),2*nlev,kptr,ie)

     if (rsplit>0) then
        kptr=kptr+2*nlev
        call edgeVpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr, ie)
     endif
  end do

  ! =============================================================
    ! Insert communications here: for shared memory, just a single
  ! sync is required
  ! =============================================================
  call t_startf('bndry_exchangeV.edge3p1', t_detail_medium)
  call bndry_exchangeV(hybrid,edge3p1,location='prim_advance_mod:3829')
  call t_stopf('bndry_exchangeV.edge3p1', t_detail_medium)
  do ie=nets,nete
     ! ===========================================================
     ! Unpack the edges for vgrad_T and v tendencies...
     ! ===========================================================
     kptr=0
     call edgeVunpack(edge3p1, elem(ie)%state%ps(:,:,np1), 1, kptr, ie)

     kptr=1
     call edgeVunpack(edge3p1, elem(ie)%state%T(:,:,:,np1), nlev, kptr, ie)

     kptr=nlev+1
     call edgeVunpack(edge3p1, elem(ie)%state%v(:,:,:,:,np1), 2*nlev, kptr, ie)

     if (rsplit>0) then
        if (0<ntrac.and.eta_ave_w.ne.0.) then
          do k=1,nlev
            !OMP_COLLAPSE_SIMD
            !DIR_VECTOR_ALIGNED
            do j=1,np
              do i=1,np
                stashdp3d(i,j,k) = elem(ie)%state%dp3d(i,j,k,np1)/elem(ie)%spheremp(i,j)
              enddo
            enddo
          end do
        endif

        corners = 0.0d0
        corners(1:np,1:np,:) = elem(ie)%state%dp3d(:,:,:,np1)
        kptr=kptr+2*nlev
        call edgeVunpack(edge3p1, elem(ie)%state%dp3d(:,:,:,np1),nlev,kptr,ie)

        if  (0<ntrac.and.eta_ave_w.ne.0.) then
          desc = elem(ie)%desc
          call edgeDGVunpack(edge3p1, corners, nlev, kptr, ie)
          corners = corners/dt2

          do k=1,nlev
            tempdp3d = elem(ie)%rspheremp(:,:)*elem(ie)%state%dp3d(:,:,k,np1)
            tempdp3d = tempdp3d - stashdp3d(:,:,k)
            tempdp3d = tempdp3d/dt2

            call distribute_flux_at_corners(cflux, corners(:,:,k), desc%getmapP)
 
            cflux(1,1,:)   = elem(ie)%rspheremp(1,  1) * cflux(1,1,:)  
            cflux(2,1,:)   = elem(ie)%rspheremp(np, 1) * cflux(2,1,:) 
            cflux(1,2,:)   = elem(ie)%rspheremp(1, np) * cflux(1,2,:) 
            cflux(2,2,:)   = elem(ie)%rspheremp(np,np) * cflux(2,2,:) 

            tempflux =  eta_ave_w*subcell_dss_fluxes(tempdp3d, np, nc, elem(ie)%metdet, cflux)
            elem(ie)%sub_elem_mass_flux(:,:,:,k) = elem(ie)%sub_elem_mass_flux(:,:,:,k) + tempflux
          end do
        end if   
     endif

     ! ====================================================
     ! Scale tendencies by inverse mass matrix
     ! ====================================================
#if (defined _LOOP_OPENMP)
!$omp parallel do num_threads(vert_num_threads)
#endif
     do k=1,nlev
       !OMP_COLLAPSE_SIMD
       !DIR_VECTOR_ALIGNED
       do j=1,np
         do i=1,np
           elem(ie)%state%T(i,j,k,np1)   = elem(ie)%rspheremp(i,j)*elem(ie)%state%T(i,j,k,np1)
           elem(ie)%state%v(i,j,1,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,1,k,np1)
           elem(ie)%state%v(i,j,2,k,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%v(i,j,2,k,np1)
         enddo
       enddo
     end do

     if (rsplit>0) then
       ! vertically lagrangian: complete dp3d timestep:
       do k=1,nlev
         !OMP_COLLAPSE_SIMD
         !DIR_VECTOR_ALIGNED
         do j=1,np
           do i=1,np
             elem(ie)%state%dp3d(i,j,k,np1)= elem(ie)%rspheremp(i,j)*elem(ie)%state%dp3d(i,j,k,np1)
           enddo
         enddo
       enddo
       ! when debugging: also update ps_v
       !elem(ie)%state%ps(:,:,np1) = elem(ie)%rspheremp(:,:)*elem(ie)%state%ps(:,:,np1)
     else
       ! vertically eulerian: complete ps_v timestep:
       !OMP_COLLAPSE_SIMD
       !DIR_VECTOR_ALIGNED
       do j=1,np
         do i=1,np
           elem(ie)%state%ps(i,j,np1) = elem(ie)%rspheremp(i,j)*elem(ie)%state%ps(i,j,np1)
         enddo
       enddo
     endif

  end do

  call t_stopf('compute_and_apply_rhs', t_detail_low)

  end subroutine compute_and_apply_rhs
 

  subroutine distribute_flux_at_corners(cflux, corners, getmapP)
    use kinds,          only : int_kind, real_kind
    use dimensions_mod, only : np, max_corner_elem
    use control_mod,    only : swest
    implicit none

    real   (kind=real_kind), intent(out)  :: cflux(2,2,2)
    real   (kind=real_kind), intent(in)   :: corners(0:np+1,0:np+1)
    integer(kind=int_kind),  intent(in)   :: getmapP(:)

    cflux = 0.0d0
    if (getmapP(swest+0*max_corner_elem) /= -1) then
      cflux(1,1,1) =                (corners(0,1) - corners(1,1))     
      cflux(1,1,1) = cflux(1,1,1) + (corners(0,0) - corners(1,1)) / 2.0d0
      cflux(1,1,1) = cflux(1,1,1) + (corners(0,1) - corners(1,0)) / 2.0d0
 
      cflux(1,1,2) =                (corners(1,0) - corners(1,1))     
      cflux(1,1,2) = cflux(1,1,2) + (corners(0,0) - corners(1,1)) / 2.0d0
      cflux(1,1,2) = cflux(1,1,2) + (corners(1,0) - corners(0,1)) / 2.0d0
    else
      cflux(1,1,1) =                (corners(0,1) - corners(1,1))     
      cflux(1,1,2) =                (corners(1,0) - corners(1,1))     
    endif
 
    if (getmapP(swest+1*max_corner_elem) /= -1) then
      cflux(2,1,1) =                (corners(np+1,1) - corners(np,1))     
      cflux(2,1,1) = cflux(2,1,1) + (corners(np+1,0) - corners(np,1)) / 2.0d0
      cflux(2,1,1) = cflux(2,1,1) + (corners(np+1,1) - corners(np,0)) / 2.0d0
 
      cflux(2,1,2) =                (corners(np  ,0) - corners(np,  1))     
      cflux(2,1,2) = cflux(2,1,2) + (corners(np+1,0) - corners(np,  1)) / 2.0d0
      cflux(2,1,2) = cflux(2,1,2) + (corners(np  ,0) - corners(np+1,1)) / 2.0d0
    else
      cflux(2,1,1) =                (corners(np+1,1) - corners(np,1))     
      cflux(2,1,2) =                (corners(np  ,0) - corners(np,1))     
    endif
 
    if (getmapP(swest+2*max_corner_elem) /= -1) then
      cflux(1,2,1) =                (corners(0,np  ) - corners(1,np  ))     
      cflux(1,2,1) = cflux(1,2,1) + (corners(0,np+1) - corners(1,np  )) / 2.0d0
      cflux(1,2,1) = cflux(1,2,1) + (corners(0,np  ) - corners(1,np+1)) / 2.0d0
 
      cflux(1,2,2) =                (corners(1,np+1) - corners(1,np  ))     
      cflux(1,2,2) = cflux(1,2,2) + (corners(0,np+1) - corners(1,np  )) / 2.0d0
      cflux(1,2,2) = cflux(1,2,2) + (corners(1,np+1) - corners(0,np  )) / 2.0d0
    else
      cflux(1,2,1) =                (corners(0,np  ) - corners(1,np  ))     
      cflux(1,2,2) =                (corners(1,np+1) - corners(1,np  ))     
    endif
 
    if (getmapP(swest+3*max_corner_elem) /= -1) then
      cflux(2,2,1) =                (corners(np+1,np  ) - corners(np,np  ))     
      cflux(2,2,1) = cflux(2,2,1) + (corners(np+1,np+1) - corners(np,np  )) / 2.0d0
      cflux(2,2,1) = cflux(2,2,1) + (corners(np+1,np  ) - corners(np,np+1)) / 2.0d0
 
      cflux(2,2,2) =                (corners(np  ,np+1) - corners(np,np  ))     
      cflux(2,2,2) = cflux(2,2,2) + (corners(np+1,np+1) - corners(np,np  )) / 2.0d0
      cflux(2,2,2) = cflux(2,2,2) + (corners(np  ,np+1) - corners(np+1,np)) / 2.0d0
    else
      cflux(2,2,1) =                (corners(np+1,np  ) - corners(np,np  ))     
      cflux(2,2,2) =                (corners(np  ,np+1) - corners(np,np  ))     
    endif
  end subroutine



  subroutine smooth_phis(phis,elem,hybrid,deriv,nets,nete,minf,numcycle)
  use dimensions_mod, only : np, np, nlev
  use control_mod, only : smooth_phis_nudt, hypervis_scaling
  use hybrid_mod, only : hybrid_t
  use edge_mod, only : edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
  use edgetype_mod, only : EdgeBuffer_t
  use bndry_mod, only : bndry_exchangev
  use element_mod, only : element_t
  use derivative_mod, only : derivative_t 
  use derivative_mod, only : laplace_sphere_wk_routine
  use time_mod, only : TimeLevel_t
  implicit none

  integer :: nets,nete
  real (kind=real_kind), dimension(np,np,nets:nete), intent(inout)   :: phis
  type (hybrid_t)      , intent(in) :: hybrid
  type (element_t)     , intent(inout), target :: elem(:)
  type (derivative_t)  , intent(in) :: deriv
  real (kind=real_kind), intent(in)   :: minf
  integer,               intent(in) :: numcycle

  ! local
  real (kind=real_kind), dimension(np,np,nets:nete) :: pstens
  real (kind=real_kind), dimension(nets:nete) :: pmin,pmax
  real (kind=real_kind) :: mx,mn
  integer :: nt,ie,ic,i,j,order,order_max, iuse
  logical :: use_var_coef


  ! compute local element neighbor min/max
  do ie=nets,nete
     pstens(:,:,ie)=minval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,ie)
  enddo
  call t_startf('bndry_exchangeV.edge3p1', t_detail_medium)
  call bndry_exchangeV(hybrid,edge3p1,location='smooth_phis1')
  call t_stopf('bndry_exchangeV.edge3p1', t_detail_medium)

  do ie=nets,nete
     call edgeVunpackMin(edge3p1, pstens(:,:,ie), 1, 0, ie)
     pmin(ie)=minval(pstens(:,:,ie))
  enddo
  do ie=nets,nete
     pstens(:,:,ie)=maxval(phis(:,:,ie))
     call edgeVpack(edge3p1,pstens(:,:,ie),1,0,ie)
  enddo
  call t_startf('bndry_exchangeV.edge3p1', t_detail_medium)
  call bndry_exchangeV(hybrid,edge3p1,location='smooth_phis2')
  call t_stopf('bndry_exchangeV.edge3p1', t_detail_medium)
  do ie=nets,nete
     call edgeVunpackMax(edge3p1, pstens(:,:,ie), 1, 0, ie)
     pmax(ie)=maxval(pstens(:,:,ie))
  enddo

  ! order = 1   grad^2 laplacian
  ! order = 2   grad^4 (need to add a negative sign)
  ! order = 3   grad^6
  ! order = 4   grad^8 (need to add a negative sign)
  order_max = 1


  use_var_coef=.true.
  if (hypervis_scaling/=0) then
     ! for tensorHV option, we turn off the tensor except for *last* laplace operator
     use_var_coef=.false.
     if (hypervis_scaling>=3) then
        ! with a 3.2 or 4 scaling, assume hyperviscosity
        order_max = 2
     endif
  endif


  do ic=1,numcycle
     pstens=phis

     do order=1,order_max-1

        do ie=nets,nete
           call laplace_sphere_wk_routine(pstens(:,:,ie),deriv,elem(ie),var_coef=use_var_coef,laplace=pstens(:,:,ie))
           call edgeVpack(edge3p1,pstens(:,:,ie),1,0,ie)
        enddo
        call t_startf('bndry_exchangeV.edge3p1', t_detail_medium)
        call bndry_exchangeV(hybrid,edge3p1,location='smooth_phis3')
        call t_stopf('bndry_exchangeV.edge3p1', t_detail_medium)
        do ie=nets,nete
           call edgeVunpack(edge3p1, pstens(:,:,ie), 1, 0, ie)
           pstens(:,:,ie)=pstens(:,:,ie)*elem(ie)%rspheremp(:,:)
        enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
     enddo
     do ie=nets,nete
       call laplace_sphere_wk_routine(pstens(:,:,ie),deriv,elem(ie),var_coef=.true.,laplace=pstens(:,:,ie))
     enddo
     if (mod(order_max,2)==0) pstens=-pstens

     do ie=nets,nete
        !  ps(t+1) = ps(t) + Minv * DSS * M * RHS
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  RHS ]
        ! but output of biharminc_wk is of the form M*RHS.  rewrite as:
        !  ps(t+1) = Minv * DSS * M [ ps(t) +  M*RHS/M ]
        ! so we can apply limiter to ps(t) +  (M*RHS)/M
#if 1
        mn=pmin(ie)
        mx=pmax(ie)
        iuse = numcycle+1  ! always apply min/max limiter
#endif
        phis(:,:,ie)=phis(:,:,ie) + &
           smooth_phis_nudt*pstens(:,:,ie)/elem(ie)%spheremp(:,:)


        ! remove new extrema.  could use conservative reconstruction from advection
        ! but no reason to conserve mean PHI.
        if (ic < iuse) then
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)>mx) phis(i,j,ie)=mx
           if (phis(i,j,ie)<mn) phis(i,j,ie)=mn
        enddo
        enddo
        endif


        ! user specified minimum
        do i=1,np
        do j=1,np
           if (phis(i,j,ie)<minf) phis(i,j,ie)=minf
        enddo
        enddo

        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%spheremp(:,:)
        call edgeVpack(edge3p1,phis(:,:,ie),1,0,ie)

     enddo
     call t_startf('bndry_exchangeV.edge3p1', t_detail_medium)
     call bndry_exchangeV(hybrid,edge3p1,location='smooth_phis4')
     call t_stopf('bndry_exchangeV.edge3p1', t_detail_medium)
     do ie=nets,nete
        call edgeVunpack(edge3p1, phis(:,:,ie), 1, 0, ie)
        phis(:,:,ie)=phis(:,:,ie)*elem(ie)%rspheremp(:,:)
     enddo
#ifdef DEBUGOMP
!$OMP BARRIER
#endif
  enddo
  end subroutine smooth_phis

  subroutine overwrite_SEdensity(elem, fvm, dt_q, hybrid,nets,nete, np1)
    use fvm_reconstruction_mod, only: reconstruction_gradient, recons_val_cart
    use dimensions_mod, only : np, nlev, nc,nhe
    use hybrid_mod, only : hybrid_t
    use edge_mod, only : edgevpack, edgevunpack, edgevunpackmax, edgevunpackmin
    use bndry_mod, only : bndry_exchangev
    use element_mod, only : element_t
    use derivative_mod, only : derivative_t
    use time_mod, only : TimeLevel_t
    use fvm_control_volume_mod, only : fvm_struct
    use spelt_mod, only : spelt_struct


    type (element_t) , intent(inout)        :: elem(:)

#if defined(_SPELT)
      type(spelt_struct), intent(inout) :: fvm(:)
#else
      type(fvm_struct), intent(inout) :: fvm(:)
#endif
    type (hybrid_t), intent(in)           :: hybrid  ! distributed parallel structure (shared)

    integer, intent(in)                     :: nets  ! starting thread element number (private)
    integer, intent(in)                     :: nete  ! ending thread element number   (private)
    integer, intent(in)                     :: np1
    integer :: ie, k

    real (kind=real_kind)             :: xp,yp, tmpval, dt_q
    integer                           :: i, j,ix, jy, starti,endi,tmpi

    real (kind=real_kind), dimension(1-nhe:nc+nhe,1-nhe:nc+nhe,6)      :: recons

    if ((nc .ne. 4) .or. (np .ne. 4)) then
      if(hybrid%masterthread) then
        print *,"You are in OVERWRITE SE AIR DENSITY MODE"
        print *,"This only works for nc=4 and np=4"
        print *,"Write a new search algorithm or pay $10000!"
      endif
      stop
    endif
#if defined(_FVM)
    do ie=nets,nete
      call reconstruction_gradient(fvm(ie)%psc, fvm(ie),recons,.false.)
      do j=1,np
        do i=1,np
          xp=tan(elem(ie)%cartp(i,j)%x)
          yp=tan(elem(ie)%cartp(i,j)%y)
          ix=i
          jy=j
          ! Search index along "x"  (bisection method)
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (xp  >  fvm%acartx(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           ix = starti
!
!         ! Search index along "y"
!           starti = 1
!           endi = nc+1
!           do
!              if  ((endi-starti) <=  1)  exit
!              tmpi = (endi + starti)/2
!              if (yp  >  fvm%acarty(tmpi)) then
!                 starti = tmpi
!              else
!                 endi = tmpi
!              endif
!           enddo
!           jy = starti

          call recons_val_cart(fvm(ie)%psc, xp,yp,fvm(ie)%spherecentroid, fvm%recons_metrics, &
               recons,ix,jy,tmpval)
          elem(ie)%state%ps(i,j,np1)= elem(ie)%state%ps(i,j,np1) +&
               dt_q*(tmpval - elem(ie)%state%ps(i,j,np1) )/(7*24*60*60)
        end do
      end do
      elem(ie)%state%ps(:,:,np1)=elem(ie)%state%ps(:,:,np1)*elem(ie)%spheremp(:,:)
     call edgeVpack(edge3p1,elem(ie)%state%ps(:,:,np1),1,0,ie)
  enddo
  call t_startf('bndry_exchangeV.edge3p1', t_detail_medium)
  call bndry_exchangeV(hybrid,edge3p1,location='prim_advance_mod:4199')
  call t_stopf('bndry_exchangeV.edge3p1', t_detail_medium)
  do ie=nets,nete
     call edgeVunpack(edge3p1, elem(ie)%state%ps(:,:,np1), 1, 0, ie)
     elem(ie)%state%ps(:,:,np1)=elem(ie)%state%ps(:,:,np1)*elem(ie)%rspheremp(:,:)
  enddo
#endif
  end subroutine overwrite_SEdensity

#ifdef CAM
!-----------------------------------------------------------------------
!
! this is a dynamics decomposition version of calc_tot_energy in cam_diagnostics
!
  subroutine calc_tot_energy_dynamics(elem,nets,nete,tl,tl_qdp,outfld_name_suffix)
    use dimensions_mod, only: npsq,qsize,nlev,np,nelemd
    use dimensions_mod, only: ldry_mass_vertical_coordinates, qsize_condensate_loading
    use physconst     , only: gravit, cpair
    use control_mod   , only: rsplit
    use element_mod   , only: element_t
    use cam_history   , only: outfld, hist_fld_active
    use constituents  , only: cnst_get_ind
    use hycoef,         only: hyai, hybi, ps0
    use cam_abortutils, only: endrun

    !------------------------------Arguments--------------------------------
    
    type (element_t) , intent(in) :: elem(:)
    integer          , intent(in) :: tl, tl_qdp,nets,nete
    
    character*(*),intent(in) :: outfld_name_suffix ! suffix for "outfld" names
    
    !---------------------------Local storage-------------------------------
    
    real(kind=real_kind) :: se(npsq)                          ! Dry Static energy (J/m2)
    real(kind=real_kind) :: ke(npsq)                          ! kinetic energy    (J/m2)
    real(kind=real_kind) :: wv(npsq)                          ! column integrated vapor       (kg/m2)
    real(kind=real_kind) :: wl(npsq)                          ! column integrated liquid      (kg/m2)
    real(kind=real_kind) :: wi(npsq)                          ! column integrated ice         (kg/m2)
    real(kind=real_kind) :: tt(npsq)                          ! column integrated test tracer (kg/m2)
    real(kind=real_kind) :: se_tmp
    real(kind=real_kind) :: ke_tmp
    real(kind=real_kind) :: wv_tmp
    real(kind=real_kind) :: wl_tmp
    real(kind=real_kind) :: wi_tmp
    real(kind=real_kind) :: tt_tmp
    real(kind=real_kind) :: ps(np,np)
    real(kind=real_kind) :: pdel
    
    integer  ie,i,j,k,ic                                  ! column, level indices
    integer :: ixcldice, ixcldliq, ixtt              ! CLDICE, CLDLIQ and test tracer indices
    character(len=16) :: name_out                  ! output field name
    character(len=16) :: name_out1,name_out2,name_out3,name_out4,name_out5,name_out6
    
    !-----------------------------------------------------------------------

    name_out1 = 'SE_'   //trim(outfld_name_suffix)
    name_out2 = 'KE_'   //trim(outfld_name_suffix)
    name_out3 = 'WV_'   //trim(outfld_name_suffix)
    name_out4 = 'WL_'   //trim(outfld_name_suffix)
    name_out5 = 'WI_'   //trim(outfld_name_suffix)
    name_out6 = 'TT_'   //trim(outfld_name_suffix)

    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2).or.hist_fld_active(name_out3).or.&
         hist_fld_active(name_out4).or.hist_fld_active(name_out5).or.hist_fld_active(name_out6)) then
    
      call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
      call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
      call cnst_get_ind('TT_UN' , ixtt    , abort=.false.)
      
      ! Compute frozen static energy in 3 parts:  KE, SE, and energy associated with vapor and liquid
      
      do ie=nets,nete!1,nelemd
        se    = 0.0D0
        ke    = 0.0D0
        wv    = 0.0D0
        wl    = 0.0D0
        wi    = 0.0D0
        tt    = 0.0D0
        
!          ps(:,:)    = sum(elem(ie)%state%dp3d(:,:,:,tl),3) + hyai(1)*ps0
        ps(:,:)    = hyai(1)*ps0
        do k = 1, nlev
          do j=1,np
            do i = 1, np
              if (ldry_mass_vertical_coordinates) then
                pdel     = elem(ie)%state%dp3d(i,j,k,tl)+&
                     SUM(elem(ie)%state%qdp(i,j,k,1:qsize_condensate_loading,tl_qdp))
!*(1.0D0+SUM(elem(ie)%state%q(i,j,k,1:qsize_condensate_loading)))
              else
                pdel     = elem(ie)%state%dp3d(i,j,k,tl)
              end if
              ps(i,j)  = ps(i,j)+pdel
              
              ke_tmp   = 0.5D0*(elem(ie)%state%v(i,j,1,k,tl)**2+ elem(ie)%state%v(i,j,2,k,tl)**2)*pdel/gravit
              se_tmp   = cpair*elem(ie)%state%T(i,j,k,tl)*pdel/gravit
              wv_tmp   =  elem(ie)%state%qdp(i,j,k,1,tl_qdp)/gravit
              
              se   (i+(j-1)*np) = se   (i+(j-1)*np) + se_tmp
              ke   (i+(j-1)*np) = ke   (i+(j-1)*np) + ke_tmp
              wv   (i+(j-1)*np) = wv   (i+(j-1)*np) + wv_tmp
            end do
          end do
        end do
        do j=1,np
          do i = 1, np
            se(i+(j-1)*np) = se(i+(j-1)*np) + elem(ie)%state%phis(i,j)*ps(i,j)/gravit
          end do
        end do

        ! Don't require cloud liq/ice to be present.  Allows for adiabatic/ideal phys.
        
        if (ixcldliq > 1) then
          do k = 1, nlev
            do j = 1, np
              do i = 1, np
                wl_tmp   = elem(ie)%state%qdp(i,j,k,ixcldliq,tl_qdp)/gravit
                wl   (i+(j-1)*np) = wl(i+(j-1)*np) + wl_tmp
              end do
            end do
          end do
        end if
        
        if (ixcldice > 1) then
          do k = 1, nlev
            do j = 1, np
              do i = 1, np
                wi_tmp   = elem(ie)%state%qdp(i,j,k,ixcldice,tl_qdp)/gravit
                wi(i+(j-1)*np)    = wi(i+(j-1)*np) + wi_tmp
              end do
            end do
          end do
        end if
        
        if (ixtt > 1) then
          do k = 1, nlev
            do j = 1, np
              do i = 1, np
                tt_tmp   = elem(ie)%state%qdp(i,j,k,ixtt,tl_qdp)/gravit
                tt   (i+(j-1)*np) = tt(i+(j-1)*np) + tt_tmp
              end do
            end do
          end do
        end if
        
        ! Output energy diagnostics
        call outfld(name_out1  ,se       ,npsq,ie)
        call outfld(name_out2  ,ke       ,npsq,ie)
        call outfld(name_out3  ,wv       ,npsq,ie)
        call outfld(name_out4  ,wl       ,npsq,ie)
        call outfld(name_out5  ,wi       ,npsq,ie)
        call outfld(name_out6  ,tt       ,npsq,ie)
      end do
    end if
  end subroutine calc_tot_energy_dynamics

  subroutine output_qdp_var_dynamics(qdp,nets,nete,outfld_name)
    use dimensions_mod, only: npsq,qsize,nlev,np,nelemd      
    use physconst     , only: gravit
    use cam_history   , only: outfld, hist_fld_active
    use constituents  , only: cnst_get_ind
    use cam_abortutils, only: endrun
    !------------------------------Arguments--------------------------------
    
    real(kind=real_kind) :: qdp(np,np,nlev,qsize,nelemd)
    character*(*),intent(in) :: outfld_name
    integer      ,intent(in) :: nets,nete
    
    !---------------------------Local storage-------------------------------
    
    real(kind=real_kind) :: qdp1(npsq),qdp2(npsq),qdp3(npsq),qdp4(npsq)
    real(kind=real_kind) :: qdp_tmp
    
    integer :: i,j,k,ie
    integer :: ixcldice, ixcldliq, ixtt
    character(len=16) :: name_out1,name_out2,name_out3,name_out4
    
    !-----------------------------------------------------------------------

    name_out1 = 'WV_'       //trim(outfld_name)
    name_out2 = 'WI_'   //trim(outfld_name)
    name_out3 = 'WL_'   //trim(outfld_name)
    name_out4 = 'TT_'   //trim(outfld_name)

    if ( hist_fld_active(name_out1).or.hist_fld_active(name_out2).or.hist_fld_active(name_out3).or.&
         hist_fld_active(name_out4)) then
      
      call cnst_get_ind('CLDICE', ixcldice, abort=.false.)
      call cnst_get_ind('CLDLIQ', ixcldliq, abort=.false.)
      call cnst_get_ind('TT_UN' , ixtt    , abort=.false.)

      do ie=nets,nete
        qdp1 = 0.0D0
        qdp2 = 0.0D0
        qdp3 = 0.0D0
        qdp4 = 0.0D0
        
        do k = 1, nlev
          do j = 1, np
            do i = 1, np
              qdp_tmp   = qdp(i,j,k,1,ie)/gravit
              qdp1   (i+(j-1)*np) = qdp1(i+(j-1)*np) + qdp_tmp
            end do
          end do
        end do

        if (ixcldice > 0) then
          do k = 1, nlev
            do j = 1, np
              do i = 1, np
                qdp_tmp   = qdp(i,j,k,ixcldice,ie)/gravit
                qdp2   (i+(j-1)*np) = qdp2(i+(j-1)*np) + qdp_tmp
              end do
            end do
          end do
        end if

        if (ixcldliq > 0) then
          do k = 1, nlev
            do j = 1, np
              do i = 1, np
                qdp_tmp   = qdp(i,j,k,ixcldliq,ie)/gravit
                qdp3   (i+(j-1)*np) = qdp3(i+(j-1)*np) + qdp_tmp
              end do
            end do
          end do
        end if

        if (ixtt > 0) then
          do k = 1, nlev
            do j = 1, np
              do i = 1, np
                qdp_tmp   = qdp(i,j,k,ixtt,ie)/gravit
                qdp4   (i+(j-1)*np) = qdp4(i+(j-1)*np) + qdp_tmp
              end do
            end do
          end do
        end if
        
        call outfld(name_out1  ,qdp1       ,npsq,ie)
        call outfld(name_out2  ,qdp2       ,npsq,ie)
        call outfld(name_out3  ,qdp3       ,npsq,ie)
        call outfld(name_out4  ,qdp4       ,npsq,ie)
      end do
    end if
  end subroutine output_qdp_var_dynamics

#endif

end module prim_advance_mod

