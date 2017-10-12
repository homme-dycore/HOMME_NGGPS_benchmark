#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
! This module should 'use' only module 'kinds'
!
module physical_constants
  ! ------------------------------
  use kinds, only : real_kind, longdouble_kind
#ifdef CAM
  use physconst, only : pi, & ! _EXTERNAL
		        g => gravit, &
                        rearth, &
                        omega, &
                        Rgas => rair, &
                        rh2o, &
                        Cp => cpair, &
                        p0 => pstd, &
                        MWDAIR => mwdry, &
                        Rwater_vapor => rh2o, &
                        Cpwater_vapor => cpwv, &
                        Cp_liq => cpliq,&
                        Cp_ice => cpice,&
                        kappa => cappa, &
                        Rd_on_Rv => epsilo, &
                        Cpd_on_Cpv, &
                        rrearth => ra
#endif
  ! -----------------------------
  implicit none

  private
#ifdef CAM
  real (kind=real_kind), public, parameter :: DD_PI = pi
  real (kind=longdouble_kind), public, parameter :: QQ_PI = 3.141592653589793238462643383279_longdouble_kind
  public                                   :: g              ! m s^-2
  public                                   :: rearth         ! m
  public                                   :: omega          ! s^-1
  public                                   :: Rgas
  public                                   :: Rh2o
  public                                   :: Cp
  public                                   :: p0             ! Pa
  public                                   :: cp_ice
  public                                   :: cp_liq
  public                                   :: MWDAIR
  public                                   :: Rwater_vapor
  public                                   :: Cpwater_vapor
  public                                   :: kappa
  public                                   :: Rd_on_Rv
  public                                   :: Cpd_on_Cpv
  public                                   :: rrearth         ! m
#else
  real (kind=real_kind), public, parameter :: DD_PI = 3.141592653589793238462643383279_real_kind
  real (kind=longdouble_kind), public, parameter :: QQ_PI = 3.141592653589793238462643383279_longdouble_kind

  ! to run Held-Suarez in JUPITER mode:
  ! make rearth,g,omega namelist variables and set to:
  ! omega = 2e-4, rearth=7e7, g=23.0
  real (kind=real_kind), public, parameter :: rearth       = 6.376D6         ! m
  real (kind=real_kind), public, parameter :: g            = 9.80616D0         ! m s^-2
  real (kind=real_kind), public            :: omega        = 7.292D-5        ! radians/s  = 2*pi/86164.
  real (kind=real_kind), public, parameter :: Rgas         = 287.04D0
  real (kind=real_kind), public, parameter :: Rh2o           = 6.02214D26*1.38065D-23/18.016D0 !CAM values
!  real (kind=real_kind), public, parameter :: Cp           = 1.00464D3 !CAM value
  real (kind=real_kind), public, parameter :: Cp           = 1005.0D0   !old HOMME value - not correct
  real (kind=real_kind), public, parameter :: p0           = 100000.0D0        ! surface pressure (mbar)
  real (kind=real_kind), public, parameter :: MWDAIR       = 28.966D0
  real (kind=real_kind), public, parameter :: Rwater_vapor = 461.50D0
  real (kind=real_kind), public, parameter :: Cpwater_vapor= 1.810D3
  real (kind=real_kind), public, parameter :: Cp_liq       = 4.188D3    !cp liquid water
  real (kind=real_kind), public, parameter :: Cp_ice       = 2.11727D3  !cp for ice
  real (kind=real_kind), public, parameter :: kappa        = Rgas/Cp !2.0D0/7.0D0 !Rgas/Cp
  real (kind=real_kind), public, parameter :: Rd_on_Rv     = Rgas/Rwater_vapor
  real (kind=real_kind), public, parameter :: Cpd_on_Cpv     = Cp/Cpwater_vapor
  real (kind=real_kind), public, parameter :: rrearth      = 1.0_real_kind/rearth         ! m
  real (kind=real_kind), public, parameter :: Lc           = 2.5D+6 ! multicloud J/Kg
#endif

end module physical_constants
