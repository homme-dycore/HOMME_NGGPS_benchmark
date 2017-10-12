!
!  This module contains code for setting initial-conditions
!  and prescribed-winds for some of the DCMIP 2012 test cases.
!
!  https://www.earthsystemcog.org/projects/dcmip-2012/test_cases
!
!  For questions, please contact: david.hall@colorado.edu
!_______________________________________________________________________

#include "config.h"

module dcmip_wrapper_mod

use kinds,          only: real_kind
use dimensions_mod, only: np, nlev, qsize, nlevp, qsize_d, nelemd
use control_mod,    only: test_case
use hybrid_mod,     only: hybrid_t
use hybvcoord_mod,  only: hvcoord_t
use time_mod,       only: timelevel_t
use element_mod,    only: element_t, elem_state_t, nt=>timelevels, derived_state_t
use dcmip_123_mod,  only: test1_advection_deformation, test1_advection_hadley
use physical_constants, only: p0, g, Rd=>Rgas

implicit none

! parameters for constant-temperature atmosphere cases

real(real_kind), parameter :: T0 = 300.d0    ! temperature (K)
real(real_kind), parameter :: H  = Rd*T0/g   ! scale height (m)

! arrays for accumulation of point-wise values

real(real_kind) :: v_m(np,np,2,nlev)         ! horizontal velocity at layer midpoints
real(real_kind) :: w_i(np,np,nlevp)          ! vertical velocity at layer interfaces
real(real_kind) :: T_m(np,np,nlev)           ! temperature at layer midpoints
real(real_kind) :: q_m(np,np,nlev,qsize_d)   ! tracer mixing ratios at midpoints
real(real_kind) :: p_i(np,np,nlevp)          ! pressure at layer interfaces
real(real_kind) :: p_m(np,np,nlev)           ! pressure at layer midpoints
real(real_kind) :: phi_s(np,np)              ! geopotential at the surface
real(real_kind) :: z_m(np,np,nlev)           ! z at layer midpoints
real(real_kind) :: omega_m(np,np,nlev)       ! vertical pressure-velocity at midpoints
real(real_kind) :: eta_dot_dpdn(np,np,nlevp) ! vertical flux at layer interfaces

! make module-level variables private to each thread
!$OMP THREADPRIVATE(v_m,T_m,q_m,p_i,p_m,phi_s,z_m,omega_m,eta_dot_dpdn)
  
  public :: set_dcmip_1_1_fields
  public :: set_dcmip_1_2_fields

CONTAINS

!_______________________________________________________________________
subroutine init_dcmip_1_1(elem, hybrid, hvcoord, nets, nete, tl)
  type(element_t),    intent(inout), target :: elem(:)
  type(hybrid_t),     intent(in)            :: hybrid
  type(hvcoord_t),    intent(inout)         :: hvcoord
  integer,            intent(in)            :: nets,nete
  type(timelevel_t),  intent(in)            :: tl

  integer:: ie
  real(real_kind) :: vflux(np,np,nlevp)

  do ie = nets,nete
    call set_dcmip_1_1_fields(elem(ie),vflux, hybrid, hvcoord, tl%n0, 0.0d0)
    elem(ie)%derived%eta_dot_dpdn= vflux
  enddo

end subroutine

!_______________________________________________________________________
subroutine set_dcmip_1_1_fields(elem, vflux, hybrid, hvcoord, n0, time)

  ! Wrapper for DCMIP 1-1: 3D Deformational Flow

  type(element_t),    intent(inout), target :: elem
  real(real_kind),    intent(inout)         :: vflux(np,np,nlevp)
  type(hybrid_t),     intent(in)            :: hybrid
  type(hvcoord_t),    intent(inout)         :: hvcoord
  integer,            intent(in)            :: n0
  real(real_kind),    intent(in)            :: time                     ! simulation time in seconds

  real(real_kind) :: T,phis,ps,u,v,w,p,z,rho,q(qsize_d),lon,lat         ! point-wise values
  integer :: ie,i,j,k                                                   ! loop indices
  integer,  parameter :: zcoords  = 1                                   ! use z coordinate

  ! accumulate layer interfaces values
  do k=1,nlev
    z = H * log(1.0d0/hvcoord%etam(k))
    p = p0*hvcoord%etam(k)  ! for evenly spaced eta

    do j=1,np; do i=1,np
      lon=elem%spherep(i,j)%lon; lat=elem%spherep(i,j)%lat
      call test1_advection_deformation(time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2),q(3),q(4))
      call cache_midpoint_values(i,j,k,u,v,p,T,q,z)
    enddo; enddo!i,j
  enddo!k

  ! accumulate layer midpoint values
  do k=1,nlevp
    z = H * log(1.0d0/hvcoord%etai(k))
    p = p0*hvcoord%etam(k)  ! for evenly spaced eta

    do j=1,np; do i=1,np
      lon=elem%spherep(i,j)%lon; lat=elem%spherep(i,j)%lat
      call test1_advection_deformation(time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2),q(3),q(4))
      call cache_interface_values(i,j,k,p,phis,rho,w)
    enddo; enddo! i,j
  enddo!k

  ! fill element data-structure
  call set_element_state(elem,hvcoord,n0,time)

  vflux = eta_dot_dpdn

  ! fill unused tracers with a checkerboard pattern
  !if(time==0d0) then
  !   call set_extra_tracers(e,n0,5,qsize) ! fill tracers 5..qsize
  !endif

end subroutine

!_______________________________________________________________________
subroutine init_dcmip_1_2(elem, hybrid, hvcoord, nets, nete, tl)
  type(element_t),    intent(inout), target :: elem(:)
  type(hybrid_t),     intent(in)            :: hybrid
  type(hvcoord_t),    intent(inout)         :: hvcoord
  integer,            intent(in)            :: nets,nete
  type(timelevel_t),  intent(in)            :: tl

  integer:: ie
  real(real_kind) :: vflux(np,np,nlevp)

  do ie = nets,nete
    call set_dcmip_1_2_fields(elem(ie),vflux, hybrid, hvcoord, tl%n0, 0.0d0)
    elem(ie)%derived%eta_dot_dpdn= vflux
  enddo

end subroutine

!_______________________________________________________________________
subroutine set_dcmip_1_2_fields(elem, vflux, hybrid, hvcoord, n0, time)

  ! Wrapper for DCMIP 1-2: Hadley-like Flow

  type(element_t),    intent(inout), target :: elem
  real(real_kind),    intent(inout)         :: vflux(np,np,nlevp)
  type(hybrid_t),     intent(in)            :: hybrid
  type(hvcoord_t),    intent(inout)         :: hvcoord
  integer,            intent(in)            :: n0
  real(real_kind),    intent(in)            :: time                     ! simulation time in seconds

  real(real_kind) :: T,phis,ps,u,v,w,p,z,rho,q(qsize_d),lon,lat         ! point-wise values
  integer :: ie,i,j,k                                                   ! loop indices
  integer,  parameter :: zcoords  = 1                                   ! use z coordinate

  ! accumulate layer interfaces values
  do k=1,nlev
    z = H * log(1.0d0/hvcoord%etam(k))
    p = p0*hvcoord%etam(k)  ! for evenly spaced eta

    do j=1,np; do i=1,np
      lon=elem%spherep(i,j)%lon; lat=elem%spherep(i,j)%lat
      call test1_advection_hadley (time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2))
      call cache_midpoint_values(i,j,k,u,v,p,T,q,z)
    enddo; enddo!i,j
  enddo!k

  ! accumulate layer midpoint values
  do k=1,nlevp
    z = H * log(1.0d0/hvcoord%etai(k))
    p = p0*hvcoord%etam(k)  ! for evenly spaced eta

    do j=1,np; do i=1,np
      lon=elem%spherep(i,j)%lon; lat=elem%spherep(i,j)%lat
      call test1_advection_hadley (time,lon,lat,p,z,zcoords,u,v,w,T,phis,ps,rho,q(1),q(2))
      call cache_interface_values(i,j,k,p,phis,rho,w)
    enddo; enddo! i,j
  enddo!k

  ! fill element data-structure
  call set_element_state(elem,hvcoord,n0,time)

  vflux = eta_dot_dpdn

  ! fill unused tracers with a checkerboard pattern
  !if(time==0d0) then
  !   call set_extra_tracers(e,n0,1,1)      ! fill tracer  1
  !   call set_extra_tracers(e,n0,3,qsize)  ! fill tracers 3..qsize
  !endif

end subroutine

!_______________________________________________________________________
subroutine set_element_state(elem, hvcoord, tl, time)

  ! set element state from accumulated field values

  type (element_t), target, intent(inout) :: elem     ! current element
  type(hvcoord_t),          intent(in)    :: hvcoord  ! vertical coordinate
  integer,                  intent(in)    :: tl       ! time level
  real(real_kind),          intent(in)    :: time     ! simulation time in seconds

  integer :: qi,k

  type(elem_state_t),    pointer :: s
  type(derived_state_t), pointer :: d

  real(real_kind) :: dp(np,np,nlev), dn(nlev)
  real(real_kind) :: dpdn_m(np,np,nlev)

  ! set state variables

  dp = p_i(:,:,2:nlevp) - p_i(:,:,1:nlev)

  s => elem%state
  s%v(:,:,:,:,tl)   = v_m
  s%T(:,:,:,tl)     = T_m
  s%dp3d(:,:,:,tl)  = dp
  s%ps(:,:,tl)    = p_i(:,:,nlevp)
  s%lnps(:,:,tl)    = log(s%ps(:,:,tl))
  s%phis            = phi_s

  ! set derived variables

  d => elem%derived
  d%phi         = z_m * g             ! set geopotential height
  d%omega       = omega_m             ! set omega/p
  d%dp          = dp

  if(time==0.0d0) then

    ! Set additional initial conditions
    s%Q(:,:,:,:)      = q_m
    do qi=1,qsize_d
      s%Qdp(:,:,:,qi,1)  = q_m(:,:,:,qi)*s%dp3d(:,:,:,tl)
      s%Qdp(:,:,:,qi,2)  = q_m(:,:,:,qi)*s%dp3d(:,:,:,tl)
    enddo

  endif

end subroutine

!_______________________________________________________________________
subroutine set_extra_tracers(elem, tl, q1,q2)

  ! fill unusued tracers with a checkerboard pattern

  type (element_t), target, intent(inout):: elem  ! current element
  integer,                  intent(in)   :: tl    ! time level
  integer,                  intent(in)   :: q1    ! starting tracer index
  integer,                  intent(in)   :: q2    ! end tracer index

  integer :: i,j,qi,k
  real(real_kind) :: term

  do qi=q1,q2
    do j=1,np
      do i=1,np
         term = sin(9.*elem%spherep(i,j)%lon)*sin(9.*elem%spherep(i,j)%lat)
         if ( term < 0. ) then
            elem%state%Q(i,j,:,qi) = 0
         else
            elem%state%Q(i,j,:,qi) = 1
         endif
      enddo
    enddo

    elem%state%Qdp(:,:,:,qi,1)  = elem%state%Q(:,:,:,qi)*elem%state%dp3d(:,:,:,tl)
    elem%state%Qdp(:,:,:,qi,2)  = elem%state%Q(:,:,:,qi)*elem%state%dp3d(:,:,:,tl)
  enddo

end subroutine

!_______________________________________________________________________
subroutine cache_midpoint_values(i,j,k,u,v,p,T,q,z)
  integer   :: i,j,k
  real(real_kind)  :: T,u,v,w,p,z,rho,q(qsize_d)

  v_m(i,j,1,k)  = u;  v_m(i,j,2,k)  = v
  p_m(i,j,k)    = p;  T_m(i,j,k)    = T
  q_m(i,j,k,:)  = q;  z_m(i,j,k)    = z
  omega_m(i,j,k)= -g*rho*w

end subroutine

!_______________________________________________________________________
subroutine cache_interface_values(i,j,k,p,phis,rho,w)
  integer   :: i,j,k
  real(real_kind)  :: p,phis,rho,w

  p_i(i,j,k) = p
  phi_s(i,j) = phis
  w_i(i,j,k) = w
  eta_dot_dpdn(i,j,k) = -g*rho*w

end subroutine

!_______________________________________________________________________
subroutine write_level_files()

  ! write hyai and hybi levels to file for dcmip1-1 test

  real(real_kind), parameter :: z_top = 12000.d0 ! top of atmosphere (m)
  real(real_kind), parameter :: c     = 2.0d0

  real(real_kind) :: zi(nlevp), eta_i(nlevp)
  real(real_kind) :: Am(nlev), Bm(nlev), Ai(nlevp), Bi(nlevp)

  integer :: k

  ! get evenly spaced z
  forall(k=1:nlevp) zi(k) = z_top - z_top*(k-1)/nlev

  ! get eta coord from z for constant temperature atmosphere
  eta_i = exp(-zi/H)

  ! get hybrid coordinates at interfaces
  Bi = ((eta_i - eta_i(1))/(1.0-eta_i(1)))**c
  Ai = eta_i - Bi

  ! get hybrid coordinate at midpoints by averaging interfaces
  Bm = 0.5d0*(Bi(2:nlevp) + Bi(1:nlev))
  Am = 0.5d0*(Ai(2:nlevp) + Ai(1:nlev))

  ! write interface-level file

  open(UNIT=12, FILE="dcmip11i-64.ascii", ACTION="write", STATUS="replace")
  write(12,*) nlevp, " ! hyai";
  do k=1,nlevp; write(12,*) Ai(k); enddo
  write(12,*) nlevp, " ! hybi";
  do k=1,nlevp; write(12,*) Bi(k); enddo
  close(12)

  ! write midpoint-level file

  open(UNIT=12, FILE="dcmip11m-64.ascii", ACTION="write", STATUS="replace")
  write(12,*) nlev, " ! hyam";
  do k=1,nlev; write(12,*) Am(k); enddo
  write(12,*) nlev, " ! hybm";
  do k=1,nlev; write(12,*) Bm(k); enddo
  close(12)

end subroutine

end module
