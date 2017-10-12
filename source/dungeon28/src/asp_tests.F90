#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

!
!  This file contains the initial condititions for the ASP 2008 summer school test cases
!  and some tests from DCMIP 2012 (2nd ASP summer school)
MODULE jw 
  use kinds, only : real_kind
  IMPLICIT NONE
  !
  !  Functions for setting up initial conditions for the Jablonowski-Williamson test case.
  !
  !  Given longitude (radians), latitude (radians), eta (pressure) and rotation_angle (degrees)
  !  the functions will return temperature, surface geopotential, zonal and meridional wind
  !  components, respectively.
  !
  !  lperturb=.FALSE. result in initial conditions for the steady-state test case.
  !  lperturb=.TRUE.  result in initial conditions for the baroclinic wave test case.
  !
  !     T   : FUNCTION temperature         (lon,lat,eta,rotation_angle)
  !     PHIS: FUNCTION surface_geopotential(lon,lat,rotation_angle)
  !     U   : FUNCTION u_wind              (lon,lat,eta,lperturb,rotation_angle)
  !     V   : FUNCTION v_wind              (lon,lat,eta,lperturb,rotation_angle)
  !
  !  The non-rotated (rotation_angle=0) version of the test cases is described in: 
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A baroclinic instability 
  !                 test case for atmospheric model dynamical cores. 
  !                 Quart. J. Roy. Meteor. Soc., 132, 2943-2975.
  !
  !                 Jablonowski, C., and D. L. Williamson, 2006: A Baroclinic Wave Test Case 
  !                 for Dynamical Cores of General Circulation Models: Model Intercomparisons, 
  !                 NCAR Technical Note, NCAR/TN-469+STR, 89 pp. 
  !  
  !  The rotated version simply rotates the initial conditions so that the spherical coordinate
  !  poles do not conicide with the earth's rotation axis. Thereby the Coriolis parameter is
  !  a function of latitude and longitude:
  !
  !      f = 2*Omega*(-cos(lon)*cos(lat)*sin(rotation_angle)+sin(lat)*cos(rotation_angle))
  !
  !  where Omega = 7.292 x 10E-5/s and rotation_angle is the angle between the flow direction
  !  and equator.
  !
  !  Author: Peter Hjort Lauritzen (NCAR, pel@ucar.edu)
  !

  REAL(kind=real_kind), PARAMETER ::       &
       eta_strato = 0.2d0     ,&! tropopause level
       u0         = 35.d0     ,&! 35 m/s
       T0         = 288.d0    ,&! global mean T at surface
       p0         = 100000.d0 ,&! global mean surface pressure
       eta0       = 0.252d0   ,&! center of jets (hybrid)
       !
       radius                 = 10.d0,             & ! radius of the perturbation
       perturbation_amplitude =  1.d0,             & ! amplitude of u perturbation 1 m/s
       perturbation_longitude = 20.d0,             & ! longitudinal position, 20E
       perturbation_latitude  = 40.d0,             & ! latitudinal position, 40N
       !
       !
       !
       Rd         = 287.d0    ,&! gas constant J/(K kg)
       g          = 9.80616d0 ,&! gravitational acceleration (m/s^2)
       a          = 6371229.d0,&! Earth's radius in m
       omega      = 7.29212d-5,&! angular velocity 1/s
       gamma      = 0.005d0   ,&! lapse rate
       pi         = 3.14159265358979323846D0,&  ! pi
       deg2rad    = pi/180.d0, &
       perturbation_latitude_tracer = 55.d0

CONTAINS
  !
  !********************************************************************
  !
  ! Temperature (equation (6) in Jablonowski and Williamson, 2006)
  !
  !********************************************************************
  !
  REAL(kind=real_kind) FUNCTION temperature(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(kind=real_kind)             :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF

    temperature  = t_mean(eta) + t_deviation(rot_lon,rot_lat,eta)
  END FUNCTION temperature
  !
  ! Horizontally averaged temperature (equation (4) and (5) in Jablonowski and Williamson (2006))
  !
  REAL(kind=real_kind) FUNCTION t_mean(eta)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: eta
    REAL(kind=real_kind)             :: exponent, delta_T
    
    exponent = Rd*gamma/g
    delta_T  = 480000.d0  ! in K, for T mean calculation
    IF (eta.gt.(eta_strato)) THEN
       t_mean = T0*eta**exponent                                ! mean temperature at each level (troposphere)
    ELSE
       t_mean = T0*eta**exponent + delta_T*(eta_strato-eta)**5  ! mean temperature at each level (stratosphere)
    ENDIF
  END FUNCTION t_mean
  !
  ! Temperature deviation from the horizontal mean 
  ! (equation (6) minus horizontally averaged temperature)
  !
  REAL(kind=real_kind) FUNCTION t_deviation(lon,lat,eta)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: eta, lon, lat
    REAL(kind=real_kind)             :: factor, phi_vertical, rot_lon, rot_lat, a_omega

    factor       = eta*pi*u0/Rd             
    phi_vertical = (eta - eta0) * 0.5d0*pi
    a_omega      = a*omega

    rot_lon = lon
    rot_lat = lat

    t_deviation = factor * 1.5d0 * SIN(phi_vertical) * (cos(phi_vertical))**0.5d0 *                        &
                  ((-2.d0*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*              &
                  u0 * (COS(phi_vertical))**1.5d0  +                                                       &
                  (8.d0/5.d0*(COS(rot_lat))**3 * ((SIN(rot_lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega*0.5d0 )
  END FUNCTION t_deviation
  !
  !**************************************************************************
  !
  ! Surface geopotential (equaiton (7) in Jablonowski and Williamson, 2006)
  !
  !**************************************************************************
  !  
  REAL(kind=real_kind) FUNCTION surface_geopotential(lon,lat,rotation_angle)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: lon, lat, rotation_angle
    REAL(kind=real_kind)             :: cos_tmp, eta_sfc, rot_lon, rot_lat, a_omega

    IF (ABS(rotation_angle)<1.0D-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF

    eta_sfc    = 1.d0
    cos_tmp    = u0 * (cos((eta_sfc-eta0)*pi*0.5d0))**1.5d0
    a_omega    = a*omega

    surface_geopotential = ((-2.d0*(SIN(rot_lat))**6 * ((COS(rot_lat))**2 + 1.d0/3.d0) + 10.d0/63.d0)*COS_tmp   &
                 + (8.d0/5.d0*(COS(rot_lat))**3 * ((SIN(rot_lat))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega)*COS_tmp
  END FUNCTION surface_geopotential
  !
  !********************************************************************
  !
  ! wind components (equation 2 in Jablonowski and Williamson, 2006)
  !
  !********************************************************************
  !  
  REAL(kind=real_kind) FUNCTION u_wind(lon,lat,eta,lperturb,rotation_angle)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: lon,lat,eta,rotation_angle
    LOGICAL, INTENT(IN)  :: lperturb
    REAL(kind=real_kind) :: cos_lat, u_lat, phi_vertical, rot_lon, rot_lat, sin_tmp, cos_tmp, r, u_perturb, v_lat
    REAL(kind=real_kind) :: perturb_lon, perturb_lat, v_tmp

    perturb_lon = perturbation_longitude*deg2rad
    perturb_lat = perturbation_latitude*deg2rad

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF

    phi_vertical = (eta - eta0) *0.5d0*pi
    u_lat = (COS(phi_vertical))**1.5d0 * 4.d0 * u0 * (sin(rot_lat))**2 * (cos(rot_lat))**2
    u_wind = u_lat

    IF (lperturb) THEN

       sin_tmp = SIN(perturb_lat)*SIN(rot_lat)
       cos_tmp = COS(perturb_lat)*COS(rot_lat)
                  
       r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-perturb_lon) )    ! great circle distance
       u_perturb = perturbation_amplitude*EXP(- (r*radius)**2 )
       u_lat     = u_perturb + u_lat
    ENDIF
    IF (ABS(rotation_angle)<1.0E-8) THEN
       u_wind = u_lat
    ELSE
       v_lat = 0.0d0
       !
       ! rotate wind components
       !
       CALL turnwi(u_lat,v_lat, u_wind,v_tmp,lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,-1)
       IF (ABS(u_wind)<1.0E-10) u_wind=0.0d0
    ENDIF
  END FUNCTION u_wind

  REAL(kind=real_kind) FUNCTION v_wind(lon,lat,eta,lperturb,rotation_angle)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: lon,lat,eta,rotation_angle
    LOGICAL, INTENT(IN)  :: lperturb
    REAL(kind=real_kind) :: cos_lat, u_lat, phi_vertical, rot_lon, rot_lat, sin_tmp, cos_tmp, r, u_perturb, v_lat
    REAL(kind=real_kind) :: perturb_lon, perturb_lat, u_tmp

    perturb_lon = perturbation_longitude*deg2rad
    perturb_lat = perturbation_latitude*deg2rad

    IF (ABS(rotation_angle)<1.0E-8) THEN
       v_wind = 0.0d0
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)


       phi_vertical = (eta - eta0) *0.5d0*pi
       u_lat = (COS(phi_vertical))**1.5d0 * 4.d0 * u0 * (sin(rot_lat))**2 * (cos(rot_lat))**2
 
       IF (lperturb) THEN
          
          sin_tmp = SIN(perturb_lat)*SIN(rot_lat)
          cos_tmp = COS(perturb_lat)*COS(rot_lat)
          
          r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-perturb_lon) )    ! great circle distance
          u_perturb = perturbation_amplitude*EXP(- (r*radius)**2 )
          u_lat     = u_perturb + u_lat
       ENDIF

       v_lat = 0.0d0
       !
       ! pole point velocities are not well-defined
       !
       IF (ABS(pi*0.5d0-lat)<1.0E-8.OR.ABS(pi*0.5d0+lat)<1.0E-8) THEN
          v_wind = 0.0d0
       ELSE
          !
          ! rotate wind components
          !
          CALL turnwi(u_lat,v_lat, u_tmp,v_wind,lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,-1)
       ENDIF
    ENDIF
  END FUNCTION v_wind




!-----------------------------------------------------------------------
! Tracer q1 and q2
!-----------------------------------------------------------------------
  REAL(kind=real_kind) FUNCTION tracer_q1_q2(lon,lat,eta,rotation_angle, eta_c)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: eta, lon, lat, rotation_angle, eta_c
    REAL(kind=real_kind) :: rot_lon, rot_lat, sin_tmp, cos_tmp, r
    REAL(kind=real_kind) :: rot_perturb_lon, rot_perturb_lat, tmp
    
    rot_perturb_lon = perturbation_longitude*deg2rad
    rot_perturb_lat = perturbation_latitude_tracer *deg2rad

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    sin_tmp = SIN(rot_perturb_lat)*SIN(rot_lat)
    cos_tmp = COS(rot_perturb_lat)*COS(rot_lat)
    r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-rot_perturb_lon) )    ! great circle distance

    tmp = EXP(- ((r*radius)**2 + ((eta-eta_c)/0.1D0)**2))
    IF (ABS(tmp)<1.0E-8) tmp = 0.0
    tracer_q1_q2 = tmp !trunk version
!phl    tracer_q1_q2 = tmp+1.0D0
  END FUNCTION tracer_q1_q2
  
!-----------------------------------------------------------------------
! Tracer q3
!-----------------------------------------------------------------------
  REAL(kind=real_kind) FUNCTION tracer_q3(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(kind=real_kind) :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    tracer_q3 = 0.5D0 * ( tanh( 3.D0*abs(rot_lat)-pi ) + 1.D0)

  END FUNCTION tracer_q3

!-----------------------------------------------------------------------
! Tracer q, absolute value of the relative vorticity of the unperturbed initial state
!           multiplied by 10^5
!-----------------------------------------------------------------------
  REAL(kind=real_kind) FUNCTION tracer_q(lon,lat,eta,rotation_angle)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: eta, lon, lat, rotation_angle
    REAL(kind=real_kind) :: rot_lon, rot_lat

    IF (ABS(rotation_angle)<1.0E-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    tracer_q = abs(-4.D0 * u0/a * (cos((eta-eta0)*pi*0.5D0))**1.5D0 * sin(rot_lat) * &
               cos(rot_lat) * (2.D0-5.D0*(sin(rot_lat))**2)) * 1.D5
    if (tracer_q < 1.D-9) tracer_q = 0.D0  !  otherwise error in netcdf file
    tracer_q=tracer_q+1.0D0
  END FUNCTION tracer_q


!----------------------------------------------------------------------- 
! Tracer slotted cylinder
!-----------------------------------------------------------------------                                                   
  REAL(kind=real_kind) FUNCTION tracer_slotted_cylinder(lon,lat,eta,rotation_angle, eta_c)
    IMPLICIT NONE
    REAL(kind=real_kind), INTENT(IN) :: eta, lon, lat, rotation_angle, eta_c
    REAL(kind=real_kind) :: rot_lon, rot_lat, sin_tmp, cos_tmp, r
    REAL(kind=real_kind) :: rot_perturb_lon, rot_perturb_lat, tmp, R0

    rot_perturb_lon = perturbation_longitude*deg2rad
    rot_perturb_lat = perturbation_latitude_tracer *deg2rad

    IF (ABS(rotation_angle)<1.0D-8) THEN
       rot_lon = lon
       rot_lat = lat
    ELSE
       CALL regrot(lon,lat,rot_lon,rot_lat,0.0d0,-0.5d0*pi+rotation_angle*deg2rad,1)
    ENDIF
    sin_tmp = SIN(rot_perturb_lat)*SIN(rot_lat)
    cos_tmp = COS(rot_perturb_lat)*COS(rot_lat)
    r = ACOS( sin_tmp + cos_tmp*COS(rot_lon-rot_perturb_lon) )    ! great circle distance
    R0=0.25D0

    if ((r .le. R0) .AND. (abs(rot_lon-rot_perturb_lon).ge. R0/6.0D0)) then
       tracer_slotted_cylinder = 1.0D0
    elseif ((r .le. R0) .AND. (abs(rot_lon-rot_perturb_lon) < R0/6.0D0) &
         .AND. (rot_lat-rot_perturb_lat < -5.0D0*R0/12.0D0)) then
       tracer_slotted_cylinder = 1.0D0
    else
       tracer_slotted_cylinder = 0.0D0
    endif
    tracer_slotted_cylinder=tracer_slotted_cylinder+1.0D0
  END FUNCTION tracer_slotted_cylinder


  !
  !******************************************************************************
  !
  ! Subroutines for rotation
  !
  !******************************************************************************
  !

  SUBROUTINE regrot(pxreg,pyreg,pxrot,pyrot,pxcen,pycen,kcall)
    IMPLICIT NONE
!
!----------------------------------------------------------------------
!
!*    conversion between regular and rotated spherical coordinates.
!*
!*    pxreg     longitudes of the regular coordinates
!*    pyreg     latitudes of the regular coordinates
!*    pxrot     longitudes of the rotated coordinates
!*    pyrot     latitudes of the rotated coordinates
!*              all coordinates given in degrees n (negative for s)
!*              and degrees e (negative values for w)
!*    pxcen     regular longitude of the south pole of the rotated grid
!*    pycen     regular latitude of the south pole of the rotated grid
!*
!*    kcall=-1: find regular as functions of rotated coordinates.
!*    kcall= 1: find rotated as functions of regular coordinates.
!
!-----------------------------------------------------------------------
!
      integer kxdim,kydim,kx,ky,kcall
      real(kind=real_kind) :: pxreg,pyreg,&
                  pxrot,pyrot,&
                  pxcen,pycen
!
!-----------------------------------------------------------------------
!
      real(kind=real_kind) zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg, &
               zsyrot,zcyrot,zcxrot,zsxrot,zpi,zpih
      integer jy,jx

      zpih = pi*0.5d0
!
      !----------------------------------------------------------------------
!
      zsycen = SIN((pycen+zpih))
      zcycen = COS((pycen+zpih))
!
      IF (kcall.eq.1) then
!
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
         zsyrot = max(zsyrot,-1.0D0)
         zsyrot = min(zsyrot,+1.0D0)
         !
         pyrot = ASIN(zsyrot)
         !
         zcyrot = COS(pyrot)
         zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot
         zcxrot = max(zcxrot,-1.0D0)
         zcxrot = min(zcxrot,+1.0D0)
         zsxrot = zcyreg*zsxmxc/zcyrot
         !
         pxrot = ACOS(zcxrot)
         !
         IF (zsxrot<0.0) pxrot = -pxrot
               !
      ELSEIF (kcall.eq.-1) then
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
         zsyreg = max(zsyreg,-1.0D0)
         zsyreg = min(zsyreg,+1.0D0)
         !
         pyreg = ASIN(zsyreg)
         !
         zcyreg = COS(pyreg)
         zcxmxc = (zcycen*zcyrot*zcxrot -&
              zsycen*zsyrot)/zcyreg
         zcxmxc = max(zcxmxc,-1.0D0)
         zcxmxc = min(zcxmxc,+1.0D0)
         zsxmxc = zcyrot*zsxrot/zcyreg
         zxmxc  = ACOS(zcxmxc)
         IF (zsxmxc<0.0) zxmxc = -zxmxc
         !
         pxreg = zxmxc + pxcen
         !
      ELSE
         WRITE(6,'(1x,''invalid kcall in regrot'')')
         STOP
      ENDIF
    END SUBROUTINE regrot

    SUBROUTINE turnwi(puarg,pvarg,pures,pvres,         &
                      pxreg,pyreg,pxrot,pyrot,   &
                      pxcen,pycen,kcall)
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!*    turn horizontal velocity components between regular and
!*    rotated spherical coordinates.
!
!*    puarg : input u components
!*    pvarg : input v components
!*    pures : output u components
!*    pvres : output v components
!*    pa    : transformation coefficients
!*    pb    :    -"-
!*    pc    :    -"-
!*    pd    :    -"-
!*    pxreg : regular longitudes
!*    pyreg : regular latitudes
!*    pxrot : rotated longitudes
!*    pyrot : rotated latitudes
!*    kxdim              : dimension in the x (longitude) direction
!*    kydim              : dimension in the y (latitude) direction
!*    kx                 : number of gridpoints in the x direction
!*    ky                 : number of gridpoints in the y direction
!*    pxcen              : regular longitude of the south pole of the
!*                         transformed grid
!*    pycen              : regular latitude of the south pole of the
!*                         transformed grid
!*
!*    kcall < 0          : find wind components in regular coordinates
!*                         from wind components in rotated coordinates
!*    kcall > 0          : find wind components in rotated coordinates
!*                         from wind components in regular coordinates
!*    note that all coordinates are given in degrees n and degrees e.
!*       (negative values for s and w)
!
!-----------------------------------------------------------------------

      integer kxdim,kydim,kx,ky,kcall
      real(kind=real_kind) puarg,pvarg,    &
               pures,pvres,    &
               pa,   pb,       &
               pc,   pd,       &
               pxreg,pyreg,    &
               pxrot,pyrot
      real(kind=real_kind) pxcen,pycen
!
!-----------------------------------------------------------------------
!
      integer jy,jx
      real(kind=real_kind) zpih,zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,&
               zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot
!
!-----------------------------------------------------------------------
!
      IF (kcall.eq.1) then
         zpih = pi*0.5d0
         zsyc = SIN(pycen+zpih)
         zcyc = COS(pycen+zpih)
         !
         zsxreg = SIN(pxreg)
         zcxreg = COS(pxreg)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         !
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         !
         pa = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
         pb = zcyc*zcxmxc*zsyreg*zsxrot - zsyc*zcyreg*zsxrot - &
              zsxmxc*zsyreg*zcxrot
         pc = zsyc*zsxmxc/zcyrot
         pd = (zsyc*zcxmxc*zsyreg + zcyc*zcyreg)/zcyrot
         !
         pures = pa*puarg + pb*pvarg
         pvres = pc*puarg + pd*pvarg
      ELSEIF (kcall.eq.-1) then
         zpih = pi*0.5d0
         zsyc = SIN(pycen+zpih)
         zcyc = COS(pycen+zpih)
         !
         zsxreg = SIN(pxreg)
         zcxreg = COS(pxreg)
         zsyreg = SIN(pyreg)
         zcyreg = COS(pyreg)
         !
         zxmxc  = pxreg - pxcen
         zsxmxc = SIN(zxmxc)
         zcxmxc = COS(zxmxc)
         !
         zsxrot = SIN(pxrot)
         zcxrot = COS(pxrot)
         zsyrot = SIN(pyrot)
         zcyrot = COS(pyrot)
         !
         pa = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot
         pb = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot -&
              zcxmxc*zsxrot*zsyrot
         pc =-zsyc*zsxrot/zcyreg
         pd = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
         !
         pures = pa*puarg + pb*pvarg
         pvres = pc*puarg + pd*pvarg
      ELSE
         write(6,'(1x,''invalid kcall in turnwi'')')
         STOP
      ENDIF
    END SUBROUTINE turnwi  

END MODULE jw








MODULE testcases_3_4_5_6
  use kinds, only : real_kind
  !=======================================================================
  !
  !  Functions for setting up initial conditions for the dynamical core test cases 3-6
  !  The input parameters depend on the test case (see comments for each test case below).
  !
  !  Given longitude (radians), latitude (radians), eta (pressure) and rotation_angle (degrees)
  !  the functions will return temperature, surface geopotential, zonal and meridional wind
  !  components, respectively.
  !
  !  Author: Christiane Jablonowski (University of Michigan, cjablono@umich.edu)
  !
  !  May/5/2008
  !
  !=======================================================================
  
implicit none
!=======================================================================
!  physical constants
!=======================================================================



  real(kind=real_kind), parameter ::           &
       Rd         = 287.04D0,     &             ! gas constant J/(K kg)
       cp         = 1004.64D0,    &             ! specific heat at constant pressure J/(K kg)
       kappa      = Rd/cp,         &             ! kappa = 2/7
       g          = 9.80616D0,    &             ! gravitational acceleration (m/s^2)
       a          = 6371229.D0,   &             ! Earth's radius in m
       pi         = 3.14159265358979323846D0,&  ! pi
       omega      = 2.D0*pi/86164.D0, &        ! Earth's angular velocity 1/s
       pih        = pi*0.5D0,     &             ! pi/2
       deg2rad    = pi/180.D0, &
       tau     = 4*86400,       &                ! period: 4 days expressed in s (used to be 3 days)
       omega_0 = pi*40000.D0 / tau              ! 0.4848 Pa/s

      real(kind=real_kind),parameter :: p0      = 100000.D0  ! reference pressure 

CONTAINS
  SUBROUTINE advection_vertical(time,eta,eta_top,eta_dot)
  real(kind=real_kind), intent(in)  :: time,eta,eta_top
  real(kind=real_kind), intent(out)  :: eta_dot

  ! local
  real(kind=real_kind) :: s

  s = 2*sqrt(  sin( pi* (eta-eta_top)/(1-eta_top) )  )
  if (s>1) s=1
     
  eta_dot = omega_0*cos(2*pi*time/tau)*sin( s*pi/2)/p0


end subroutine



  SUBROUTINE advection (tracer_variant, lon, lat, height, rotation_angle,  &
                        u_wind, v_wind, temperature, surface_geopotential, &
                        surface_pressure, q5, q6)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      character*2, intent(in) :: tracer_variant                          ! identifies test variant 'yy', here tracers
                                                                         ! e.g. 0 : no tracer, set to zero
                                                                         !      5 : tracer q5 only
                                                                         !     56 : both tracers q5 and q6
      real(kind=real_kind), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               height,                                  & ! height of the level in m
                               rotation_angle                             ! alpha in degrees
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(kind=real_kind), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure,                        & ! surface pressure in Pa
                               q5,                                      & ! tracer q5
                               q6                                         ! tracer q6
!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
 real(kind=real_kind),parameter  ::      u0      = (2.D0*pi*a)/(12.D0*86400.D0), & ! circumference / 12 days
                             T0      = 300.D0,                         & ! constant temperature
                             H       = Rd * T0 / g,                     & ! scale height
                             RR      = 1/3.D0,                         & ! horizontal half width divided by 'a'
                             ZZ      = 1000.D0,                        & ! vertical half width
                             z0      = 4500.D0,                        & ! center point in z
                             lambda0 = 1.5D0*pi,                       & ! center point in longitudes
                             phi0    = 0.D0,                           & ! center point in latitudes
                             slot    = 1.D0/8.D0                        ! half width of the slot in radians
!----------------------------------------------------------------------- 
!     local variables
!-----------------------------------------------------------------------                             
      real(kind=real_kind) :: alpha
      real(kind=real_kind) :: sin_tmp, cos_tmp
      real(kind=real_kind) :: d1, d2, s, r
      
      alpha = rotation_angle*deg2rad
!-----------------------------------------------------------------------
!    initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * (cos(lat)*cos(alpha) + sin(lat)*cos(lon)*sin(alpha))
      v_wind = -u0 *  sin(lon) * sin(alpha)
!-----------------------------------------------------------------------
!     initialization of the vertical velocity: 
!     must be implemented in the dynamical core
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     initialize T (temperature)
!-----------------------------------------------------------------------
      temperature = T0
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0.D0
!-----------------------------------------------------------------------
!     initialize PS (surface pressure)
!-----------------------------------------------------------------------
      surface_pressure = p0
!-----------------------------------------------------------------------
!     Tracer variables
!-----------------------------------------------------------------------
      q5 = 0.D0   ! default
      q6 = 0.D0   ! default
!-----------------------------------------------------------------------
!     tracer q5
!-----------------------------------------------------------------------
      if (tracer_variant(1:1) == '5' .or. tracer_variant(2:2) == '5') then
        sin_tmp = sin(lat) * sin(phi0)
        cos_tmp = cos(lat) * cos(phi0)
        r  = ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))       ! great circle distance without 'a'
        d1 = min( 1.D0, (r/RR)**2 + ((height-z0)/ZZ)**2 )
        q5 = 0.5D0 * (1.D0 + cos(pi*d1))
      endif
!-----------------------------------------------------------------------
!     tracer q6
!-----------------------------------------------------------------------
      if (tracer_variant(1:1) == '6' .or. tracer_variant(2:2) == '6') then
        sin_tmp = sin(lat) * sin(phi0)
        cos_tmp = cos(lat) * cos(phi0)
        r  = ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))       ! great circle distance without 'a'
        d2 = (r/RR)**2 + ((height-z0)/ZZ)**2
        if (d2 <= 1.D0) then
          q6 = 1.D0
        else
          q6 = 0.D0
        endif
        if ((height > z0) .and. ((phi0-slot) < lat .and. lat < (phi0+slot)) ) q6 = 0.D0   ! slotted ellipse               
      endif
  end subroutine advection
 
!==========================================================================================
! Rossby_Haurwitz wave, wavenumber 4
!==========================================================================================
  SUBROUTINE Rossby_Haurwitz (lon, lat, pressure,                                &
                              u_wind, v_wind, temperature, surface_geopotential, &
                              surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(kind=real_kind), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               pressure                                   ! pressure at full model level in Pa
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(kind=real_kind), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure                           ! surface pressure in Pa

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
      real(kind=real_kind),parameter :: u0      = 50.D0,                          &   ! reference wind
                            T0      = 288.D0,                         &   ! reference temperature
                            n       = 4.D0,                           &   ! wavenumber
                            MM      = u0/(n*a),                        &   ! parameter M and p_ref=95500.
                            KK      = u0/(n*a),                        &   ! parameter K
                            gamma   = 0.0065D0,                       &   ! lapse rate in K/m
                            p_ref   = 95500.D0                            ! reference pressure
                            
!----------------------------------------------------------------------- 
!     local
!----------------------------------------------------------------------- 
      real(kind=real_kind) :: tmp1, tmp2, tmp3
      real(kind=real_kind) :: sin_lat, cos_lat, sin_slat, cos_slat
      real(kind=real_kind) :: exponent_1, exponent_2
      real(kind=real_kind) :: AA, BB, CC
      real(kind=real_kind) :: phis_perturb
      
!-----------------------------------------------------------------------
!     initialize the wind components
!-----------------------------------------------------------------------
      cos_lat = cos(lat)
      sin_lat = sin(lat)
      tmp1 = a * MM * cos_lat
      tmp2 = a * KK * cos_lat**(n-1.D0)*(n*sin_lat**2 - cos_lat**2)
      tmp3 = -a * KK * n * cos_lat**(n-1.D0) * sin_lat
      u_wind = tmp1 + tmp2 * cos(n*lon)
      v_wind = tmp3 * sin(n*lon)
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0.D0     
!-----------------------------------------------------------------------
!     initialize surface pressure and temperature
!-----------------------------------------------------------------------
      tmp1       = gamma/(g*T0)
      tmp2       = a*a
      exponent_1 = g/(gamma*Rd)
      exponent_2 = (gamma*Rd)/g
      
      cos_lat = cos(lat)
      AA = tmp2 * (0.5D0 * MM*(2.D0*omega+MM) * cos_lat**2 + 0.25D0 * KK**2 * cos_lat**(2.D0*n) * &
                  ( (n+1.D0)*cos_lat**2 + (2.D0*n*n - n - 2.D0)) - 0.5D0*n*n*KK**2 * cos_lat**(2.D0*(n-1)))
      BB = tmp2 * (2.D0*(omega+MM)*KK/((n+1.D0)*(n+2.D0)) * cos_lat**n * &
                   ( (n*n + 2.D0*n +2.D0) - (n+1.D0)**2 * cos_lat**2 ))
      CC = tmp2 * (0.25D0 * KK**2 * cos_lat**(2.D0*n) * ( (n+1.D0)*cos_lat**2 - (n+2.D0)))
      phis_perturb = AA + BB * cos(n*lon) + CC * cos(2.D0*n*lon)
      surface_pressure = p_ref * (1.D0 + tmp1*phis_perturb)**exponent_1   ! surface pressure
      temperature      = T0 * (pressure/p_ref)**exponent_2                 ! temperature
      
  end subroutine Rossby_Haurwitz

!==========================================================================================
! Mountain induced Rossby wave
!==========================================================================================
  SUBROUTINE mountain_Rossby (lon, lat, pressure,                                &
                              u_wind, v_wind, temperature, surface_geopotential, &
                              surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(kind=real_kind), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               pressure                                   ! pressure at full model level in Pa

!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(kind=real_kind), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure                           ! surface pressure in Pa
!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
      real(kind=real_kind),parameter :: u0      = 20.D0,                          &   ! 20 m/s
                            T0      = 288.D0,                         &   ! temperature
                            N2      = g*g/(cp*T0),                     &   ! squared Brunt Vaisala frequency N^2
                            h0      = 2000.D0,                        &   ! amplitude of the mountain, 2km
                            d       = 1500.D3,                      &   ! half width 1500 km
                            lambda0 = 0.5D0*pi,                       &   ! center point in longitudes
                            phi0    = pi/6.D0,                        &   ! center point in latitudes
                            p_sp    = 93000.D0                            ! pressure at the South Pole in Pa
!-----------------------------------------------------------------------
!   local variables
!-----------------------------------------------------------------------
      real(kind=real_kind) :: sin_tmp, cos_tmp
      real(kind=real_kind) :: tmp1, tmp2, tmp3
      real(kind=real_kind) :: r

!-----------------------------------------------------------------------
!    initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * cos(lat)
      v_wind = 0.D0
!-----------------------------------------------------------------------
!     initialize T (temperature)
!-----------------------------------------------------------------------
      temperature = T0
!-----------------------------------------------------------------------
!     initialize surface geopotential and surface pressure
!-----------------------------------------------------------------------
      tmp1 = (a * N2 * u0)/(2.D0 * g*g * kappa) * (u0/a + 2.D0 * omega)
      tmp2 = N2 / (g*g * kappa)
      sin_tmp = sin(lat) * sin(phi0)
      cos_tmp = cos(lat) * cos(phi0)
      tmp3 = tmp1*((sin(lat))**2 - 1.D0)
      r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))   ! great circle distance with 'a'
      surface_geopotential = g*h0 * exp(-(r/d)**2)        ! Gaussian profile of the mountain
      surface_pressure     = p_sp * exp( -tmp3 - tmp2*surface_geopotential)

  end subroutine mountain_Rossby
  
!==========================================================================================
! gravity waves
!==========================================================================================
  SUBROUTINE gravity_wave (choice, lon, lat, height,                          &
                           u_wind, v_wind, temperature, surface_geopotential, &
                           surface_pressure)
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      integer, intent(in)   :: choice                                     ! identifies test variant 'x'
                                                                          ! e.g. 0 : no Coriolis, N=0.01 1/s, u0=0  m/s
                                                                          !      1 : no Coriolis, N=0.01 1/s, u0=40 m/s
      real(kind=real_kind), intent(in)  :: lon,                                     & ! longitude in radians
                               lat,                                     & ! latitude in radians
                               height                                     ! height of the level in m
!-----------------------------------------------------------------------
!     input parameters
!-----------------------------------------------------------------------
      real(kind=real_kind), intent(out) :: u_wind,                                  & ! zonal wind in m/s
                               v_wind,                                  & ! meridional wind in m/s
                               temperature,                             & ! temperature in K
                               surface_geopotential,                    & ! surface geopotential in m^2/s^2
                               surface_pressure

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
    real(kind=real_kind),parameter ::      T0      = 300.D0,                         & ! reference temperature
                            RR      = a/3.D0,                         & ! half width   
                            Lz      = 20.D3,                        & ! vertical wave length, 20 km 
                            delta_theta = 10.D0                         ! potential temperature perturbation amplitude                   

!----------------------------------------------------------------------- 
!     local variables
!----------------------------------------------------------------------- 
      real(kind=real_kind) :: sin_tmp, cos_tmp
      real(kind=real_kind) :: tmp1, tmp2, tmp3
      real(kind=real_kind) :: theta                                                    ! potential temperature
      real(kind=real_kind) :: theta_mean
      real(kind=real_kind) :: pres                                                     ! pressure
      real(kind=real_kind) :: r                                                        ! great circle distance

!-----------------------------------------------------------------------
!     more test case parameters
!----------------------------------------------------------------------- 
      real(kind=real_kind) :: N2,                                                  &   ! squared Brunt-Vaisala frequency
                  S,                                                   &   ! parameter
                  u0,                                                  &   ! background wind speed
                  lambda0,                                             &   ! center point in longitudes 
                  phi0,                                                &   ! center point in latitudes
                  gw_omega,                                            &   ! rotation
                  p_eq
                  
!-----------------------------------------------------------------------
!    initialize parameters
!-----------------------------------------------------------------------
      lambda0 = pi                                                         ! center point in longitudes
      p_eq    = p0                                                         ! surface pressure at the equator

      select case (choice)
      case (0) 
        N2 = 1.D-4                                                      ! squared Brunt Vaisala frequency N^2
        u0 = 0.D0                                                         ! background wind speed
        phi0 = 0.D0                                                       ! center point in latitudes (0 deg)
        gw_omega = 0.D0                                                   ! no rotation
      case (1) 
        N2 = (g*g)/(cp*T0)                                                 ! squared Brunt Vaisala frequency N^2
        u0 = 0.D0                                                         ! background wind speed
        phi0 = 0.D0                                                       ! center point in latitudes (0 deg)
        gw_omega = 0.D0                                                   ! no rotation
      case (2) 
        N2 = (g*g)/(cp*T0)                                                 ! squared Brunt Vaisala frequency N^2
        u0 = 40.D0                                                        ! background wind speed
        phi0 = 0.D0                                                       ! center point in latitudes (0 deg)
        gw_omega = 0.D0                                                   ! no rotation
      case (3) 
        N2 = (g*g)/(cp*T0)                                                 ! squared  Brunt Vaisala frequency N^2
        u0 = 0.D0                                                         ! background wind speed
        phi0 = pi/4.D0                                                    ! center point in latitudes (45 deg)
        gw_omega = omega                                                   ! Earth's rotation
      end select
      S = g*g/(cp*N2)                                                      ! parameter

!-----------------------------------------------------------------------
!     initialize the wind components
!-----------------------------------------------------------------------
      u_wind = u0 * cos(lat)
      v_wind = 0.D0
!-----------------------------------------------------------------------
!     initialize surface geopotential
!-----------------------------------------------------------------------
      surface_geopotential = 0.D0
!-----------------------------------------------------------------------
!     initialize surface pressure
!-----------------------------------------------------------------------
      tmp1 = (a * N2 * u0)/(2.D0 * g*g * kappa) * (u0/a + 2.D0 * gw_omega)
      surface_pressure = p_eq * exp( - tmp1*((sin(lat))**2) )
!-----------------------------------------------------------------------
!     initialize temperature
!-----------------------------------------------------------------------                 
      pres    = p0 * ( (1.D0 - S/T0) + S/T0 * exp(- (N2*height)/g) )**(cp/Rd)
      sin_tmp = sin(lat) * sin(phi0)
      cos_tmp = cos(lat) * cos(phi0)
      theta_mean = T0 /( T0/S * ((pres/p0)**kappa - 1.D0) + 1.D0 )
      r = a * ACOS (sin_tmp + cos_tmp*cos(lon-lambda0))     ! great circle distance with radius
      if (r < RR) then
        tmp1 = 0.5D0 * (1.D0 + cos(pi*r/RR))
      else
        tmp1 = 0.D0
      endif
      theta = theta_mean + delta_theta * tmp1 * sin(2.D0*pi*height/Lz)
      temperature = theta * (pres/p0)**kappa

  end subroutine gravity_wave


END MODULE testcases_3_4_5_6












module asp_tests
!
!  This module contains the initial condititions for the baroclinic
!  instability probelms in Jablonowski and Williamson, QJR (2006) 132 
!
use element_mod, only : element_t, timelevels
use fvm_control_volume_mod, only : fvm_struct
use spelt_mod, only : spelt_struct
use hybrid_mod, only : hybrid_t
use hybvcoord_mod, only : hvcoord_t 
use kinds, only : real_kind

use physical_constants, only : p0, g
use dimensions_mod, only : nlev,np, qsize,nc,ntrac, nep
use control_mod, only : test_case, u_perturb
use cube_mod, only : rotate_grid
use jw, only : u_wind, v_wind, temperature, surface_geopotential, tracer_q1_q2,&
     tracer_q3, perturbation_longitude, perturbation_latitude, deg2rad ! _EXTERNAL
use parallel_mod, only : abortmp

implicit none



!-----------------------------------------------------------------------
!     Physical Parameters
!-----------------------------------------------------------------------

!!!! replace a with earth rad

	real(kind=real_kind), parameter ::	a	= 6371220.0d0,	&	! Earth's Radius (m)
				Rd 	= 287.0d0,	&	! Ideal gas const dry air (J kg^-1 K^1)
				cp	= 1004.5d0,	&	! Specific heat capacity (J kg^-1 K^1)
				pi	= 4.d0*atan(1.d0)       ! pi





private

public :: asp_baroclinic
public :: asp_tracer, asp_advection_vertical, asp_gravity_wave, asp_rossby, asp_mountain

public :: dcmip2_schar

contains

subroutine dcmip2_schar(elem,hybrid,hvcoord,nets,nete)
!=======================================================================================================!
!  
! new version of JW Baroclinic test case which also allows for rotation
! Code from Jablonowski and Lauritzen
! 
!=======================================================================================================!
    use prim_si_mod, only : preq_hydrostatic

    type(element_t), intent(inout) :: elem(:)

    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete

!   local
    real (kind=real_kind)  :: lat,lon,z,p, u, v, w, t, phis, ps, rho, q 
    integer :: i,j,k,ie


    do ie=nets,nete
    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
 


             !elem(ie)%state%T(i,j,k,:)  = temperature(lon,lat,hvcoord%etam(k),rotate_grid)

             call test2_steady_state_mountain (lon,lat,p,z,0,.true.,hvcoord%hyam(k),hvcoord%hybm(k),u,v,w,t,phis,ps,rho,q)

             elem(ie)%state%T(i,j,k,:)  = t
             elem(ie)%state%v(i,j,1,k,:)  = u
             elem(ie)%state%v(i,j,2,k,:)  = v
             !what to do about w?

             elem(ie)%state%phis(i,j) = phis
             elem(ie)%state%pswet(i,j) = ps ! is ps_v surf pressure in homme? yes

!p, z, w, rho --- not used

             if (qsize>=1) then
                elem(ie)%state%Q(i,j,k,1) = q
             endif
          enddo

       enddo
    enddo
    enddo

end subroutine dcmip2_schar


!=========================================================================
! Test 2-0:  Steady-State Atmosphere at Rest in the Presence of Orography
!=========================================================================
SUBROUTINE test2_steady_state_mountain (lon,lat,p,z,zcoords,hybrid_eta,hyam,hybm,u,v,w,t,phis,ps,rho,q)

IMPLICIT NONE
!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------

	real(8), intent(in)  :: lon, &		! Longitude (radians)
				lat, &		! Latitude (radians)
				z, &		! Height (m)
				hyam, &		! A coefficient for hybrid-eta coordinate, at model level midpoint
				hybm		! B coefficient for hybrid-eta coordinate, at model level midpoint

	logical, intent(in)  :: hybrid_eta      ! flag to indicate whether the hybrid sigma-p (eta) coordinate is used
                                                ! if set to .true., then the pressure will be computed via the 
                                                !    hybrid coefficients hyam and hybm, they need to be initialized
                                                ! if set to .false. (for pressure-based models): the pressure is already pre-computed
                                                !    and is an input value for this routine 
                                                ! for height-based models: pressure will always be computed based on the height and
                                                !    hybrid_eta is not used

	real(8), intent(inout) :: p		! Pressure  (Pa)
				
	integer,  intent(in) :: zcoords 	! 0 or 1 see below

	real(8), intent(out) :: u, & 		! Zonal wind (m s^-1)
				v, &		! Meridional wind (m s^-1)
				w, &		! Vertical Velocity (m s^-1)
				t, & 		! Temperature (K)
				phis, & 	! Surface Geopotential (m^2 s^-2)
				ps, & 		! Surface Pressure (Pa)
				rho, & 		! density (kg m^-3)
				q 		! Specific Humidity (kg/kg)

	! if zcoords = 1, then we use z and output p
	! if zcoords = 0, then we compute or use p
        !
	! In hybrid-eta coords: p = hyam p0 + hybm ps
        !
        ! The grid-point based initial data are computed in this routine. 

!-----------------------------------------------------------------------
!     test case parameters
!----------------------------------------------------------------------- 
	real(8), parameter :: 	T0      = 300.d0,		&	! temperature (K)
                                gamma   = 0.0065d0,             &       ! temperature lapse rate (K/m)      
                            	lambdam = 3.d0*pi/2.d0,		&	! mountain longitude center point (radians)   
                            	phim    = 0.d0,			&	! mountain latitude center point (radians)    
                            	h0      = 2000.d0,		&	! peak height of the mountain range (m)
                            	Rm      = 3.d0*pi/4.d0,		&	! mountain radius (radians)
                            	zetam   = pi/16.d0,		&	! mountain oscillation half-width (radians) 
                            	ztop    = 12000.d0			! model top (m)         
                            
      real(8) :: height							! Model level heights (m)
      real(8) :: r							! Great circle distance (radians)
      real(8) :: zs							! Surface elevation (m)
      real(8) :: exponent                                               ! exponent: g/(Rd * gamma)
      real(8) :: exponent_rev                                           ! reversed exponent


!-----------------------------------------------------------------------
!    compute exponents 
!-----------------------------------------------------------------------
     exponent     = g/(Rd*gamma)
     exponent_rev = 1.d0/exponent

!-----------------------------------------------------------------------
!    PHIS (surface geopotential)
!-----------------------------------------------------------------------
      
	r = acos( sin(phim)*sin(lat) + cos(phim)*cos(lat)*cos(lon - lambdam) )

	if (r .lt. Rm) then

		zs = (h0/2.d0)*(1.d0+cos(pi*r/Rm))*cos(pi*r/zetam)**2.d0   ! mountain height

	else

		zs = 0.d0

	endif

	phis = g*zs

!-----------------------------------------------------------------------
!    PS (surface pressure)
!-----------------------------------------------------------------------

	ps = p0 * (1.d0 - gamma/T0*zs)**exponent


!-----------------------------------------------------------------------
!    HEIGHT AND PRESSURE
!-----------------------------------------------------------------------

	! Height and pressure are aligned (p = p0 * (1.d0 - gamma/T0*z)**exponent)

	if (zcoords .eq. 1) then

		height = z
		p = p0 * (1.d0 - gamma/T0*z)**exponent

	else

                if (hybrid_eta) p = hyam*p0 + hybm*ps              ! compute the pressure based on the surface pressure and hybrid coefficients
		height = T0/gamma * (1.d0 - (p/p0)**exponent_rev)  ! compute the height at this pressure

	endif

!-----------------------------------------------------------------------
!    THE VELOCITIES ARE ZERO (STATE AT REST)
!-----------------------------------------------------------------------

	! Zonal Velocity

	u = 0.d0

	! Meridional Velocity

	v = 0.d0

        ! Vertical Velocity

        w = 0.d0

!-----------------------------------------------------------------------
!    TEMPERATURE WITH CONSTANT LAPSE RATE
!-----------------------------------------------------------------------

	t = T0 - gamma*height

!-----------------------------------------------------------------------
!    RHO (density)
!-----------------------------------------------------------------------

	rho = p/(Rd*t)

!-----------------------------------------------------------------------
!     initialize Q, set to zero 
!-----------------------------------------------------------------------

	q = 0.d0

END SUBROUTINE test2_steady_state_mountain









subroutine asp_baroclinic(elem,hybrid,hvcoord,nets,nete, fvm)
!=======================================================================================================!
!  
! new version of JW Baroclinic test case which also allows for rotation
! Code from Jablonowski and Lauritzen
! 
!=======================================================================================================!
    use prim_si_mod, only : preq_hydrostatic
    use jw, only : tracer_slotted_cylinder
    type(element_t), intent(inout) :: elem(:)
#if defined(_SPELT)
      type(spelt_struct), optional, intent(inout) :: fvm(:)
#else
      type(fvm_struct), optional, intent(inout) :: fvm(:)
#endif
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev)
    integer :: i,j,k,ie,idex,t
    logical :: perturbation 
    real(kind=real_kind) :: eta_c


!og, added slotted ellipse from asp advection test as tracer 5
    real (kind=real_kind)  :: height, q6, d2, cos_tmp, sin_tmp, perturb_lat, perturb_lon, r
    real (kind=real_kind)  :: p(np,np,nlev),dp(np,np,nlev)
    real (kind=real_kind)  :: ph(np,np,nlev+1)
    real (kind=real_kind)  :: phi(np,np,nlev)
    real (kind=real_kind),parameter  ::    RR      = 1.0d0/3.0d0,                         & ! horizontal half width divided by 'a'
                             ZZ      = 1000.0d0,                        & ! vertical half width
                             z0      = 4500.0d0,                        & ! center point in z
                             slot    = 1.0d0/8.0d0                        ! half width of the slot in radians




    perturbation = (u_perturb/=0)

    if (perturbation) then
       if (u_perturb /= 1) then
          call abortmp("Error: baroclinc test only supports u_perturb = 0 or 1")
       endif
    endif

    !print *,'JW initial condition  u_perturb, rotation = ',u_perturb,rotate_grid,perturbation
    do ie=nets,nete
    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             elem(ie)%state%v(i,j,1,k,:) = u_wind(lon,lat,hvcoord%etam(k),perturbation,rotate_grid)
             elem(ie)%state%v(i,j,2,k,:) = v_wind(lon,lat,hvcoord%etam(k),perturbation,rotate_grid)

             elem(ie)%state%T(i,j,k,:)  = temperature(lon,lat,hvcoord%etam(k),rotate_grid)
          enddo
          elem(ie)%state%phis(i,j) = surface_geopotential(lon,lat,rotate_grid)
          elem(ie)%state%pswet(i,j) = p0
       enddo
    enddo
    enddo


! initialize passive tracers.  just for testing.  
! first tracer:  temperature at T=0
    
    if (qsize>=1) then
       idex=1
       eta_c = 0.6
       do ie=nets,nete
          do j=1,np
             do i=1,np
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                do k=1,nlev
                   elem(ie)%state%Q(i,j,k,idex) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
!                   elem(ie)%state%Q(i,j,k,idex)=0.0
                enddo
             enddo
          enddo
       enddo
    endif
    
    if (qsize>=2) then
       idex=2
       eta_c = 1
       do ie=nets,nete
          do j=1,np
             do i=1,np
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                do k=1,nlev
                   elem(ie)%state%Q(i,j,k,idex) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
                enddo
             enddo
          enddo
       enddo
    endif
    
    if (qsize>=3) then
       idex=3
       do ie=nets,nete
          do j=1,np
             do i=1,np
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                do k=1,nlev
                   elem(ie)%state%Q(i,j,k,idex) = tracer_q3(lon,lat,hvcoord%etam(k),rotate_grid)
                enddo
             enddo
          enddo
       enddo
    endif
    
    
    if (qsize>=4) then
       idex=4
       do ie=nets,nete
          elem(ie)%state%Q(:,:,:,idex) = 1 !phl trunk code
          !add phl start
          !      do j=1,np
          !      do i=1,np
          !         lon = elem(ie)%spherep(i,j)%lon
          !         lat = elem(ie)%spherep(i,j)%lat
          !         do k=1,nlev
          !           elem(ie)%state%Q(i,j,k,idex) = tracer_slotted_cylinder(lon,lat,hvcoord%etam(k),rotate_grid,eta_c)
          !         enddo
          !      end do
          !      end do
          !add phl end
       enddo
    endif
    
    
    if (qsize>=5) then
      do idex=5,qsize
         do ie=nets,nete
            do j=1,np
               do i=1,np
                  lon = elem(ie)%spherep(i,j)%lon
                  lat = elem(ie)%spherep(i,j)%lat
                  do k=1,nlev
                     elem(ie)%state%Q(i,j,k,idex) = tracer_slotted_cylinder(lon,lat,hvcoord%etam(k),rotate_grid,eta_c)
                  enddo
               enddo
            enddo
         enddo
      end do
   endif

!if (qsize>=5) then
!   idex=5
!
!   perturb_lon = perturbation_longitude*deg2rad
!   perturb_lat = perturbation_latitude*deg2rad
!
!   ! now compute PHI, needed to init tracers:	
!   do ie=nets,nete
!      ! compute height
!      do k=1,nlev+1
!         do j=1,np
!            do i=1,np
!               ph(i,j,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%ps_v(i,j,1)
!            end do
!         end do
!      end do
!      do k=1,nlev
!         do j=1,np
!            do i=1,np
!               p(i,j,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps_v(i,j,1)
!               dp(i,j,k)  = ph(i,j,k+1) - ph(i,j,k)
!            end do
!         end do
!      end do
!      call preq_hydrostatic(phi,elem(ie)%state%phis,elem(ie)%state%T(:,:,:,1),p,dp)
!
!      ! init slotted ellipse 
!      do j=1,np
!      do i=1,np
!         lon = elem(ie)%spherep(i,j)%lon
!         lat = elem(ie)%spherep(i,j)%lat
!         do k=1,nlev
!
!	    height=phi(i,j,k)/g
!            sin_tmp = sin(lat) * sin(perturb_lat)
!            cos_tmp = cos(lat) * cos(perturb_lat)
!            r  = ACOS (sin_tmp + cos_tmp*cos(lon-perturb_lon))       ! great circle distance without 'a'
!            d2 = (r/RR)**2 + ((height-z0)/ZZ)**2
!            if (d2 <= 1.0d0) then
!               q6 = 1.0d0
!            else
!               q6 = 0.0d0
!            endif
!            if ((height > z0) .and. ((perturb_lat-slot) < lat .and. lat < (perturb_lat+slot)) ) q6 = 0.0d0   ! slotted ellipse     
!
!            elem(ie)%state%Q(i,j,k,idex) = q6
!
!         enddo
!      enddo
!      enddo
!   enddo
!
!endif

    
    if (qsize>=6) then
       ! set the rest of the tracers to the temperature
       do idex=6,qsize
          do ie=nets,nete
             do j=1,np
                do i=1,np
                   lon = elem(ie)%spherep(i,j)%lon
                   lat = elem(ie)%spherep(i,j)%lat
                   do k=1,nlev
                      elem(ie)%state%Q(:,:,:,idex) = temperature(lon,lat,hvcoord%etam(k),rotate_grid)
                   enddo
                enddo
             enddo
          enddo
       enddo
    endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! fvm tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (present(fvm)) then
#if defined(_SPELT)
  if (ntrac>=1) then
     idex=1
     ! First tracer will be airdensity
     do ie=nets,nete
        do j=1,nep
          do i=1,nep
            fvm(ie)%c(i,j,:,idex,:) = 1.0D0
          enddo
        enddo
     enddo
  endif

  if (ntrac>=2) then
     idex=2
     do ie=nets,nete
        do j=1,nep
        do i=1,nep
           lon = fvm(ie)%asphere(i,j)%lon
           lat = fvm(ie)%asphere(i,j)%lat
           do k=1,nlev
              fvm(ie)%c(i,j,k,idex,:) = tracer_q3(lon,lat,hvcoord%etam(k),rotate_grid)
           enddo
        enddo
        enddo
     enddo
  endif

  if (ntrac>=3) then
     idex=3
     do ie=nets,nete
        do j=1,nep
        do i=1,nep
           lon = fvm(ie)%asphere(i,j)%lon
           lat = fvm(ie)%asphere(i,j)%lat
           do k=1,nlev
              fvm(ie)%c(i,j,k,idex,:) = tracer_q3(lon,lat,hvcoord%etam(k),rotate_grid)
           enddo
        enddo
        enddo
     enddo
  endif
  if (ntrac>=4) then
     idex=4
     eta_c = 0.6
     do ie=nets,nete
        do j=1,nep
        do i=1,nep
           lon = fvm(ie)%asphere(i,j)%lon
           lat = fvm(ie)%asphere(i,j)%lat
           do k=1,nlev
              fvm(ie)%c(i,j,k,idex,:) = 1.0D0 !tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
           enddo
        enddo
        enddo
     enddo
  endif

  if (ntrac>=5) then
     do idex=5,ntrac
     do ie=nets,nete
        do k=1,nlev
        do t=1,timelevels
          do j=1,nep
            do i=1,nep
              fvm(ie)%c(i,j,k,idex,t) = 1.0D0
            enddo
          enddo
        enddo
        enddo
     enddo
     enddo
  endif
#else
  !
  ! CSLAM tracers
  !
  if (ntrac>=1) then
     idex=2
     eta_c = 0.6
     do ie=nets,nete
        do j=1,nc
           do i=1,nc
              lon = fvm(ie)%centersphere(i,j)%lon
              lat = fvm(ie)%centersphere(i,j)%lat
              do k=1,nlev
                 fvm(ie)%c(i,j,k,idex,:) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
              enddo
           enddo
        enddo
     enddo
  endif

  if (ntrac>=2) then
     idex=2
     eta_c = 1.0
     do ie=nets,nete
        do j=1,nc
           do i=1,nc
              lon = fvm(ie)%centersphere(i,j)%lon
              lat = fvm(ie)%centersphere(i,j)%lat
              do k=1,nlev
                 fvm(ie)%c(i,j,k,idex,:) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
              enddo
           enddo
        enddo
     enddo
  endif

  if (ntrac>=3) then
     idex=3
     do ie=nets,nete
        do j=1,nc
        do i=1,nc
           lon = fvm(ie)%centersphere(i,j)%lon
           lat = fvm(ie)%centersphere(i,j)%lat
           do k=1,nlev
              fvm(ie)%c(i,j,k,idex,:) = tracer_q3(lon,lat,hvcoord%etam(k),rotate_grid)
           enddo
        enddo
        enddo
     enddo
  endif
  if (ntrac>=4) then
     idex=4
     eta_c = 0.6
     do ie=nets,nete
        do j=1,nc
        do i=1,nc
           lon = fvm(ie)%centersphere(i,j)%lon
           lat = fvm(ie)%centersphere(i,j)%lat
           do k=1,nlev
              fvm(ie)%c(i,j,k,idex,:) = 1.0D0 !
           enddo
        enddo
        enddo
     enddo
  endif

  if (ntrac>=5) then
     idex=5
     do ie=nets,nete
        do k=1,nlev
           do j=1,nc
              do i=1,nc
                 lon = fvm(ie)%centersphere(i,j)%lon
                 lat = fvm(ie)%centersphere(i,j)%lat
                 fvm(ie)%c(i,j,k,idex,:) =  tracer_slotted_cylinder(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
              enddo
           enddo
        enddo
     enddo    
  endif

  if (ntrac>=6) then
     do idex=5,ntrac
     do ie=nets,nete
        do k=1,nlev
          do j=1,nc
            do i=1,nc
              lon = fvm(ie)%centersphere(i,j)%lon
              lat = fvm(ie)%centersphere(i,j)%lat
              fvm(ie)%c(i,j,k,idex,:) =  temperature(lon,lat,hvcoord%etam(k),rotate_grid)
            enddo
          enddo
        enddo
     enddo
     enddo
  endif

#endif

endif





end subroutine asp_baroclinic



!==========================================================================================
! pure 3D advection, time-dependent
!==========================================================================================

subroutine asp_advection_vertical(time,hvcoord,eta_dot_dpdn)
    use dimensions_mod, only : np,nlevp
    use testcases_3_4_5_6, only : advection_vertical ! _EXTERNAL
    real (kind=real_kind)  :: time
    real (kind=real_kind)  :: eta_dot_dpdn(np,np,nlevp)
    type (hvcoord_t)       :: hvcoord

    ! local	
    real (kind=real_kind)  :: dpdn(np,np,nlevp)
    real (kind=real_kind)  :: eta_dot,p1,p2
    integer i,j,k

    ! compute dp/dn at interfaces	
    ! note that for this problem: ps = ps0 = 1000 mb.  
    ! thus  p = A p0 + B ps = (A + B) p0 = eta p0
    ! dp/dn = p0
    eta_dot_dpdn(:,:,1)=0
    eta_dot_dpdn(:,:,nlevp)=0
    do k=2,nlevp-1
       call advection_vertical(time,hvcoord%etai(k),hvcoord%etai(1),eta_dot)
       ! we need eta_dot_dpdn
       eta_dot_dpdn(:,:,k) = eta_dot*p0   ! dpdn(i,j,k)
    enddo
     
end subroutine

subroutine asp_tracer(elem,hybrid,hvcoord,nets,nete)
    use dimensions_mod, only : np,nlev
    use cube_mod, only : rotate_grid
    use prim_si_mod, only : preq_hydrostatic
    use testcases_3_4_5_6, only : advection,g ! _EXTERNAL

    type(element_t), intent(inout) :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev),height,q5,q6
    real (kind=real_kind)  :: u_wind,v_wind,temperature
    real (kind=real_kind)  :: surface_geopotential, surface_pressure
    real (kind=real_kind)  :: p(np,np,nlev),dp(np,np,nlev)
    real (kind=real_kind)  :: ph(np,np,nlev+1)
    real (kind=real_kind)  :: phi(np,np,nlev)
    integer :: i,j,k,ie,idex


    ! initialize phis, ps
    height=0   ! not needed for ps, phis
    do ie=nets,nete
    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             call advection("  ", lon, lat, height, rotate_grid,  &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure, q5, q6)

             elem(ie)%state%phis(i,j) = surface_geopotential
             elem(ie)%state%pswet(i,j) = surface_pressure
             elem(ie)%state%T(i,j,k,:)  = temperature
             elem(ie)%state%v(i,j,1,k,:) = u_wind
             elem(ie)%state%v(i,j,2,k,:) = v_wind
          enddo
       enddo
    enddo
    enddo

    ! now compute PHI, needed to init tracers:	
    do ie=nets,nete
       ! compute height
       do k=1,nlev+1
          do j=1,np
             do i=1,np
                ph(i,j,k)   = hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*elem(ie)%state%pswet(i,j)
             end do
          end do
       end do
       do k=1,nlev
          do j=1,np
             do i=1,np
                p(i,j,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%pswet(i,j)
                dp(i,j,k)  = ph(i,j,k+1) - ph(i,j,k)
             end do
          end do
       end do
       call preq_hydrostatic(phi,elem(ie)%state%phis,elem(ie)%state%T(:,:,:,1),p,dp)

       ! init tracers
       do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             call advection("56", lon, lat, phi(i,j,k)/g, rotate_grid,  &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure, q5, q6)


             if (qsize>=2) then
                idex=1
                elem(ie)%state%Q(i,j,k,idex) = q5
                idex=2
                elem(ie)%state%Q(i,j,k,idex) = q6
             endif
          enddo
       enddo
       enddo
    enddo


end subroutine




subroutine asp_rossby(elem,hybrid,hvcoord,nets,nete)
    use dimensions_mod, only : np,nlev
    use prim_si_mod, only : preq_hydrostatic
    use testcases_3_4_5_6, only : rossby_haurwitz ! _EXTERNAL

    type(element_t), intent(inout) :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev),height,q5,q6
    real (kind=real_kind)  :: u_wind,v_wind,temperature
    real (kind=real_kind)  :: surface_geopotential, surface_pressure
    real (kind=real_kind)  :: p(np,np,nlev),dp(np,np,nlev)
    real (kind=real_kind)  :: ph(np,np,nlev+1)
    real (kind=real_kind)  :: phi(np,np,nlev)
    integer :: i,j,k,ie,idex


    ! initialize surface pressure 
    p=0	
    do ie=nets,nete
    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          call Rossby_Haurwitz(lon, lat, p(i,j,1), &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure)
          elem(ie)%state%phis(i,j) = surface_geopotential
          elem(ie)%state%pswet(i,j) = surface_pressure
       enddo
    enddo
    enddo

    ! now we can compute the pressure
    do ie=nets,nete
       do k=1,nlev
       do j=1,np
       do i=1,np
          p(i,j,k)   = hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%pswet(i,j)
       end do
       end do
       end do

       do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             call Rossby_Haurwitz(lon, lat, p(i,j,k), &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure)
             elem(ie)%state%v(i,j,1,k,:) = u_wind
             elem(ie)%state%v(i,j,2,k,:) = v_wind
             elem(ie)%state%T(i,j,k,:)  = temperature
          enddo
       enddo
    enddo
    enddo


end subroutine



subroutine asp_mountain(elem,hybrid,hvcoord,nets,nete)
    use dimensions_mod, only : np,nlev
    use prim_si_mod, only : preq_hydrostatic
    use testcases_3_4_5_6, only : mountain_rossby ! _EXTERNAL

    type(element_t), intent(inout) :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev),height,q5,q6
    real (kind=real_kind)  :: u_wind,v_wind,temperature
    real (kind=real_kind)  :: surface_geopotential, surface_pressure, p
    integer :: i,j,k,ie,idex


    ! initialize surface pressure 
    p=0	 ! not used by initial condition
    do ie=nets,nete
       do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat
          do k=1,nlev
             call mountain_Rossby(lon, lat, p, &
                 u_wind, v_wind, temperature, surface_geopotential, &
                  surface_pressure)
             if (k==1) then
                elem(ie)%state%phis(i,j) = surface_geopotential
                elem(ie)%state%pswet(i,j) = surface_pressure
             endif
             elem(ie)%state%v(i,j,1,k,:) = u_wind
             elem(ie)%state%v(i,j,2,k,:) = v_wind
             elem(ie)%state%T(i,j,k,:)  = temperature
          enddo
       enddo
       enddo

    enddo


end subroutine



subroutine asp_gravity_wave(elem,hybrid,hvcoord,nets,nete,choice)
    use dimensions_mod, only : np,nlev
    use cube_mod, only : rotate_grid
    use prim_si_mod, only : preq_hydrostatic
    use testcases_3_4_5_6, only : gravity_wave ! _EXTERNAL

    type(element_t), intent(inout) :: elem(:)
    type (hvcoord_t)                  :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer :: nets,nete, choice

!   local
    real (kind=real_kind)  :: lat,lon,eta(nlev),height,dz,q5,q6
    real (kind=real_kind)  :: u_wind,v_wind,temperature
    real (kind=real_kind)  :: surface_geopotential, surface_pressure
    real (kind=real_kind)  :: p(np,np,nlev),dp(np,np,nlev)
    real (kind=real_kind)  :: ph(np,np,nlev+1)
    real (kind=real_kind)  :: phi(np,np,nlev)
    integer :: i,j,k,ie,idex


    ! initialize surface pressure, geopotential
    do ie=nets,nete
    do j=1,np
       do i=1,np
          lon = elem(ie)%spherep(i,j)%lon
          lat = elem(ie)%spherep(i,j)%lat

          do k=1,nlev
             dz = 500
             height = (nlev-k)*dz + dz/2
             call gravity_wave (choice, lon, lat, height,                          &
                           u_wind, v_wind, temperature, surface_geopotential, &
                           surface_pressure)
             if (k.eq.1) then
                elem(ie)%state%phis(i,j) = surface_geopotential
                elem(ie)%state%pswet(i,j) = surface_pressure
             endif
             elem(ie)%state%v(i,j,1,k,:) = u_wind
             elem(ie)%state%v(i,j,2,k,:) = v_wind
             elem(ie)%state%T(i,j,k,:)  = temperature
          enddo

       enddo
    enddo
    enddo


end subroutine
end module asp_tests


!================================================================================================
! This is the "toy" chemistry module.
!================================================================================================

module chemistry  
  use kinds, only : real_kind
  implicit none
  private
  save

  public :: initial_value        ! initialize x and x2
  public :: compute_chemistry_FQ 


  real(kind=real_kind), parameter :: xt_constant = 4.D-6
  real(kind=real_kind), parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164D0
  real(kind=real_kind), parameter :: half_pi = pi*0.5D0
  real(kind=real_kind), parameter :: k1_lat_center =   20.d0*pi/180.0D0
  real(kind=real_kind), parameter :: k1_lon_center =  300.d0*pi/180.0D0

contains

  !===============================================================================
  !  Solar photolysis rate and recombination rate
  !===============================================================================

  subroutine k_vals( lat, lon, k1, k2 )

    !-----------------------------------------------------------------------
    ! Arguments:
    !-----------------------------------------------------------------------
    real(kind=real_kind), intent(in)    :: lat, lon  ! latitude and longitude, radians
    real(kind=real_kind), intent(out)   :: k1, k2    ! reaction rates

    k1 = 1.0D0*max(0.d0,sin(lat)*sin(k1_lat_center) + cos(lat)*cos(k1_lat_center)*cos(lon-k1_lon_center))
    k2 = 1.D0

    return

  end subroutine k_vals

  !===============================================================================
  !  Tendencies of x and x2
  !===============================================================================

  subroutine tendency( lat, lon, x, x2, dt, x_f, x2_f ,dp)

    !-----------------------------------------------------------------------
    ! Arguments:
    !-----------------------------------------------------------------------

    real(kind=real_kind), intent(in)    :: lat, lon  ! latitude and longitude, radians
    real(kind=real_kind), intent(in)    :: x, x2     ! molar mixing ratio of x and x2
    real(kind=real_kind), intent(in)    :: dt        ! size of physics time step
    real(kind=real_kind), intent(in)    :: dp        ! 

    real(kind=real_kind), intent(out)   :: x_f, x2_f  ! time rate of change of x and x2

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    real(kind=real_kind) :: r, det, expdt, el ! useful algebraic quantities used in the computation
    real(kind=real_kind) :: k1, k2            ! reaction rates
    real(kind=real_kind) :: xt                ! quantity that should be conserved

    call k_vals( lat, lon, k1, k2 )

    r = k1 / (4.D0*k2)
    xt = x + 2.D0* x2

    det = sqrt( r*r + 2.D0*r*xt )
    expdt = exp( -4.D0*k2*det*dt )

    if ( abs(det * k2 * dt) .gt. 1D-16 ) then
       el = (1.D0 - expdt) /det /dt
    else
       el = 4.D0*k2
    endif

    x_f  = -el * (x - det + r)*(x + det + r) / (1.D0 + expdt + dt*el*(x + r))
    x2_f = -x_f / 2.D0
!    x_f=x_f*dp
!    x2_f=x2_f*dp
    !
    ! make sure mixing ratio's are not going negative and make sure linear correlations are not
    ! violated in physics dynamics coupling
    !
!    if (x<1.0D-12.or.x2<1.0D-12) then
!       x_f = 0.0D0
!       x2_f= 0.0D0
!    else
!       if (x+dt*dp*x_f<0.0D0) then
!          x_f  = -x/dt
!          x2_f = -0.5D0*x_f
!       end if
!       if (x2+dt*dp*x2_f<0.0D0) then
!          x2_f = -x2/dt
!          x_f  = -2.D0*x2_f
!       end if
!    end if

    return

  end subroutine tendency

  !===============================================================================
  !  Compute initial values
  !===============================================================================

  subroutine initial_value( lat, lon, x, x2 )

    !-----------------------------------------------------------------------
    ! Arguments:
    !-----------------------------------------------------------------------

    real(kind=real_kind), intent(in)  :: lat, lon  ! latitude and longitude, radians
    real(kind=real_kind), intent(out) :: x, x2   ! molar mixing ratio of x and x2

    !-----------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------

    real(kind=real_kind) :: r, det  ! useful algebraic forms
    real(kind=real_kind) :: k1, k2  ! reaction rates

    call k_vals( lat, lon, k1, k2 )

    r = k1 / (4.D0*k2)
    det = sqrt(r*r + 2.D0*xt_constant*r)

    x  = (det-r)
    x2 = xt_constant/2.D0 - (det-r)/2.D0

    return

  end subroutine initial_value

  subroutine compute_chemistry_FQ(elem,hybrid,hvcoord,nets,nete, fvm,np0,n_fvm)
    use element_mod, only : element_t, timelevels
    use fvm_control_volume_mod, only : fvm_struct
    use hybrid_mod, only : hybrid_t
    use hybvcoord_mod, only : hvcoord_t 
    use kinds, only : real_kind  
    use dimensions_mod, only : nlev,np, qsize,nc,ntrac
    use parallel_mod, only : abortmp
    use time_mod, only : dt_phys
    implicit none
    
    type(element_t)           , intent(inout) :: elem(:)
    type(fvm_struct), optional, intent(inout) :: fvm (:)
    
    type (hvcoord_t)            :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in)         :: nets,nete,np0,n_fvm
    
    !   local
    real (kind=real_kind)  :: lat,lon,Cl,Cl2,Cl_f,Cl2_f,dp
    integer                :: i,j,k,ie,idex,t    
    if (qsize<=2.and.ntrac<=2) then
       call abortmp("qsize or ntrac must be >2 for terminator chemistry - ABORT")
    end if

    if (qsize>=3) then
       do ie=nets,nete
          do j=1,np
             do i=1,np
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                do k=1,nlev
                   Cl  = elem(ie)%state%Q(i,j,k,2)
                   Cl2 = elem(ie)%state%Q(i,j,k,3)

                   dp = ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                        ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(i,j,np0)
                  call tendency(lat, lon, Cl, Cl2, dt_phys, Cl_f, Cl2_f,dp)
                   elem(ie)%derived%FQ(i,j,k,2,1)= Cl_f*dp
                   elem(ie)%derived%FQ(i,j,k,3,1)= Cl2_f*dp
                end do
             end do
          end do
       end do
    end if
    if (qsize>=4) then
       do ie=nets,nete
          do j=1,np
             do i=1,np
                do k=1,nlev
                   elem(ie)%derived%FQ(i,j,k,4:qsize,1)= 0.0D0
                end do
             end do
          end do
       end do
    end if
    if (ntrac>=3) then
       do ie=nets,nete
          do j=1,nc
             do i=1,nc
                lon = fvm(ie)%centersphere(i,j)%lon
                lat = fvm(ie)%centersphere(i,j)%lat
                do k=1,nlev
                   Cl  = fvm(ie)%c(i,j,k,2,n_fvm)
                   Cl2 = fvm(ie)%c(i,j,k,3,n_fvm)
                   dp = fvm(ie)%dp_fvm(i,j,k,n_fvm)
                   call tendency(lat, lon, Cl, Cl2, dt_phys, Cl_f, Cl2_f,dp)
                   fvm(ie)%fc(i,j,k,2) = Cl_f
                   fvm(ie)%fc(i,j,k,3) = Cl2_f
                end do
             end do
          end do
       end do
    end if
    if (ntrac>=4) then
       do ie=nets,nete
          do j=1,nc
             do i=1,nc
                do k=1,nlev
                   fvm(ie)%fc(i,j,k,4:ntrac) = 0.0D0
                end do
             end do
          end do
       end do
    end if
  end subroutine compute_chemistry_FQ


end module chemistry





!===============================================================================
!  Test program
!===============================================================================

!program test!
!
!use chemistry, only: initial_value, tendency!
!
!    integer, parameter :: r8 = selected_real_kind (12)
!    real(kind=real_kind) :: x, x2, xf, x2f, lat, lon, dt
! 
!    lat = -45.00  ! latitude
!    lon =   8.65  ! longitude  (near terminator)
!    dt = 1800.00  ! seconds
!
!    call initial_value(lat, lon, x, x2)
!    print *, " "
!    print *, " Initial Values"
!    print *, " Latitude: ",lat," Longitude: ",lon
!    print *, " x: ",x," x2: ",x2
!
!    call tendency(lat, lon, x, x2, dt, xf, x2f)
!    print *, " "
!    print *, " Tendency should be zero for the initial values"
!!    print *, " x rate ",xf," x2 rate: ",x2f
 !   print *, " "
!
!end




MODULE umjs_baroclinic_wave_mod
  use kinds, only : real_kind
  use physical_constants, only: p0, Rgas, cp, kappa, g, rearth, DD_pi, omega, rh2o
  use parallel_mod, only : abortmp
  IMPLICIT NONE
  !=======================================================================
  !
  !  Date:  July 29, 2015
  !
  !  Functions for setting up idealized initial conditions for the
  !  Ullrich, Melvin, Jablonowski and Staniforth (QJ, 2014) baroclinic instability.
  !
  !  SUBROUTINE baroclinic_wave_sample(
  !    deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,w,t,phis,ps,rho,q)
  !
  !  Options:
  !     deep    deep atmosphere (1 = yes or 0 = no)
  !    moist    include moisture (1 = yes or 0 = no)
  !    pertt    type of perturbation (0 = exponential, 1 = stream function)
  !        X    Earth scaling factor
  !
  !  Given a point specified by: 
  !      lon    longitude (radians) 
  !      lat    latitude (radians) 
  !      p/z    pressure (Pa) / height (m)
  !  zcoords    1 if z is specified, 0 if p is specified
  !
  !  the functions will return:
  !        p    pressure if z is specified and zcoords = 1 (Pa)
  !        u    zonal wind (m s^-1)
  !        v    meridional wind (m s^-1)
  !        t    temperature (K)
  !   thetav    virtual potential temperature (K)
  !     phis    surface geopotential (m^2 s^-2)
  !       ps    surface pressure (Pa)
  !      rho    density (kj m^-3)
  !        q    water vapor mixing ratio (kg/kg)
  !
  !
  !  Author: Paul Ullrich
  !          University of California, Davis
  !          Email: paullrich@ucdavis.edu
  !
  !=======================================================================
  ! use phyical constants
  !=======================================================================
  real(kind=real_kind), parameter ::           &
       Rd         = Rgas,          &
       Rv         = Rh2o,          &
       a          = rearth,        &
       pi         = DD_pi,         &
       pih        = 0.5D0*pi,      &
       Mvap       = 0.608D0         ! Ratio of molar mass dry air/water vapor
  
  
  !=======================================================================
  !    Test case parameters
  !=======================================================================
  REAL(8), PARAMETER ::               &
       T0E        = 310.d0     ,      & ! temperature at equatorial surface (K)
       T0P        = 240.d0     ,      & ! temperature at polar surface (K)
       B          = 2.d0       ,      & ! jet half-width parameter
       KK         = 3.d0       ,      & ! jet width parameter
       lapse      = 0.005d0             ! lapse rate parameter
  
  REAL(8), PARAMETER ::               &
       pertu0     = 0.5d0      ,      & ! SF Perturbation wind velocity (m/s)
       pertr      = 1.d0/6.d0  ,      & ! SF Perturbation radius (Earth radii)
       pertup     = 1.0d0      ,      & ! Exp. perturbation wind velocity (m/s)
       pertexpr   = 0.1d0      ,      & ! Exp. perturbation radius (Earth radii)
       pertlon    = pi/9.d0    ,      & ! Perturbation longitude
       pertlat    = 2.d0*pi/9.d0,     & ! Perturbation latitude
       pertz      = 15000.d0   ,      & ! Perturbation height cap
       dxepsilon  = 1.d-5               ! Small value for numerical derivatives
  
  REAL(8), PARAMETER ::               &
       moistqlat  = 2.d0*pi/9.d0,     & ! Humidity latitudinal width
       moistqp    = 34000.d0,         & ! Humidity vertical pressure width
       moistq0    = 0.018d0             ! Maximum specific humidity


  INTEGER, PARAMETER  :: deep = 0! Deep (1) or Shallow (0) test case
  INTEGER, PARAMETER  :: pertt= 0!! 0: exponential, 1: streamfunction
  REAL(8), PARAMETER  :: bigx    = 1.0  ! factor for a reduced size earth

  REAL(8) :: lat,lon,ztop


  REAL(8), SAVE :: cly_mass_exact, cly_mass_exact_fvm
  
  private
  public :: umjs_baroclinic, prim_printstate_par_terminator
CONTAINS
  
  !=======================================================================
  !    Generate the baroclinic instability initial conditions
  !=======================================================================
  SUBROUTINE baroclinic_wave_test(moist,p,ptop,z,&
       u,v,temp,thetav,phis,ps,rho,q)
    use dimensions_mod, only : ldry_mass_vertical_coordinates    
    IMPLICIT NONE

    
    !-----------------------------------------------------------------------
    !     input/output params parameters at given location
    !-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: &
         moist        ! Moist (1) or Dry (0) test case
    
    
    REAL(8), INTENT(IN) :: &
         p            ,&! Pressure at the full model level (Pa)
         ptop           !

    
    REAL(8), INTENT(OUT) :: &
         u,          & ! Zonal wind (m s^-1)
         v,          & ! Meridional wind (m s^-1)
         temp,       & ! Temperature (K)
         thetav,     & ! Virtual potential temperature (K)
         phis,       & ! Surface Geopotential (m^2 s^-2)
         ps,         & ! Surface Pressure (Pa)
         rho,        & ! density (kg m^-3)
         q,          & ! water vapor mixing ratio (kg/kg)
         z               ! Altitude (m)
    z = iterate_z_given_pressure(p,ldry_mass_vertical_coordinates,ptop) !Find height of pressure surface
    ps = p0
    call uv_given_z(z,u,v)
    temp = Tv_given_z(z)
    phis = 0.d0
    if (moist .eq. 1) then
       q = qv_given_moist_pressure(moist_pressure_given_z(z))
    else
       q = 0.d0                  ! dry
    end if

    
    ! Convert virtual temperature to temperature
    temp = temp / (1.d0 + Mvap * q)
    rho = p / (Rd * temp * (1.d0 + 0.61d0 * q))
    thetav = temp * (1.d0 + 0.61d0 * q) * (p0 / p)**(Rd / cp)
    if (ldry_mass_vertical_coordinates) then
       q=q/(1-q)
    end if
  END SUBROUTINE baroclinic_wave_test
  
  REAL(kind=real_kind) FUNCTION qv_given_moist_pressure(pwet)
    implicit none
    REAL(8), INTENT(IN)  :: pwet

    REAL(8)  :: eta
    eta = pwet/p0 
    if (eta.gt.0.1D0) then  ! intialize q if p > 100 hPa
       qv_given_moist_pressure = moistq0 * exp(- (lat/moistqlat)**4)          & 
            * exp(- ((eta-1.d0)*p0/moistqp)**2)
    else
       qv_given_moist_pressure = 1.d-12              ! above 100 hPa set q to 1e-12 to avoid supersaturation
    endif
  END FUNCTION qv_given_moist_pressure

  REAL(kind=real_kind) FUNCTION weight_of_water_vapor_given_z(z,ptop)
    implicit none
    REAL(8), INTENT(IN)  :: z,ptop
    real (8)  :: dx,xm,xr,gaussw(10),gaussx(10),integral, tmp1, tmp2
    real(8)   :: temp, rho, qv, pressure, z1, z2, Tv,pwet, ztmp
    integer   :: jgw
    SAVE gaussw,gaussx    
    DATA gaussw/0.1527533871307258D0,0.1491729864726037D0,0.1420961093183820D0,0.1316886384491766D0,0.1181945319615184D0,&
         0.1019301198172404D0,0.0832767415767048D0,0.0626720483341091D0,0.0406014298003869D0,0.0176140071391521D0/
    DATA gaussx/0.0765265211334973D0,0.2277858511416451D0,0.3737060887154195D0,0.5108670019508271D0,0.6360536807265150D0,&
         0.7463319064601508D0,0.8391169718222188D0,0.9122344282513259D0,0.9639719272779138D0,0.9931285991850949D0/
    
    z1=z
    z2=ztop
    xm=0.5*(z1+z2)
    xr=0.5*(z2-z1)
    integral=0 
    do jgw=1,10 
       dx=xr*gaussx(jgw)
       ztmp=xm+dx
       pwet = moist_pressure_given_z(ztmp); qv= qv_given_moist_pressure(pwet);Tv= Tv_given_z(ztmp)
       tmp1=g*pwet*qv/(Rd*Tv)
       
       ztmp=xm-dx
       pwet = moist_pressure_given_z(ztmp); qv= qv_given_moist_pressure(pwet);Tv= Tv_given_z(ztmp)
       tmp2=g*pwet*qv/(Rd*Tv)
       integral=integral+gaussw(jgw)*(tmp1+tmp2)
    enddo
    integral=xr*integral    ! Scale the answer to the range of integration.    
    
    weight_of_water_vapor_given_z = integral
  end FUNCTION weight_of_water_vapor_given_z

  REAL(kind=real_kind) FUNCTION weight_of_dry_air_given_z(z,ptop)
    implicit none
    REAL(8), INTENT(IN)  :: z,ptop
    real (8)  :: dx,xm,xr,gaussw(10),gaussx(10),integral, tmp1, tmp2
    real(8)   :: temp, rho, qv, pressure, z1, z2, Tv,pwet, ztmp
    integer   :: jgw
    SAVE gaussw,gaussx    
    DATA gaussw/0.1527533871307258D0,0.1491729864726037D0,0.1420961093183820D0,0.1316886384491766D0,0.1181945319615184D0,&
         0.1019301198172404D0,0.0832767415767048D0,0.0626720483341091D0,0.0406014298003869D0,0.0176140071391521D0/
    DATA gaussx/0.0765265211334973D0,0.2277858511416451D0,0.3737060887154195D0,0.5108670019508271D0,0.6360536807265150D0,&
         0.7463319064601508D0,0.8391169718222188D0,0.9122344282513259D0,0.9639719272779138D0,0.9931285991850949D0/
    
    z1=z
    z2=ztop
    xm=0.5*(z1+z2)
    xr=0.5*(z2-z1)
    integral=0 
    do jgw=1,10 
       dx=xr*gaussx(jgw)
       ztmp=xm+dx
       pwet = moist_pressure_given_z(ztmp); qv= qv_given_moist_pressure(pwet);Tv= Tv_given_z(ztmp)
       tmp1=g*pwet*(1-qv)/(Rd*Tv)
       
       ztmp=xm-dx
       pwet = moist_pressure_given_z(ztmp); qv= qv_given_moist_pressure(pwet);Tv= Tv_given_z(ztmp)
       tmp2=g*pwet*(1-qv)/(Rd*Tv)
       integral=integral+gaussw(jgw)*(tmp1+tmp2)
    enddo
    integral=xr*integral    ! Scale the answer to the range of integration.    
    
    weight_of_dry_air_given_z = integral+ptop
  end FUNCTION weight_of_dry_air_given_z

  REAL(kind=real_kind) FUNCTION iterate_z_given_pressure(p,ldry,ptop)
    implicit none
    REAL(8), INTENT(IN)  :: &
         p,              &! Pressure (Pa)
         ptop             ! Pressure (Pa)
    LOGICAL, INTENT(IN)  :: ldry
    
    INTEGER :: ix
    
    REAL(8) :: z0, z1, z2
    REAL(8) :: p0, p1, p2
    z0 = 0.d0
    z1 = 10000.d0

    if (ldry) then
       p0 = weight_of_dry_air_given_z(z0,ptop)
       p1 = weight_of_dry_air_given_z(z1,ptop)
    else
       p0 =  moist_pressure_given_z(z0)
       p1 =  moist_pressure_given_z(z1)
    endif

    DO ix = 1, 100
       z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)
       if (ldry) then
          p2 = weight_of_dry_air_given_z(z2,ptop)
       else
          p2 = moist_pressure_given_z(z2)   
       end if
       
       IF (ABS((p2 - p)/p) .lt. 1.0d-13) THEN
          EXIT
       END IF
       
       z0 = z1
       p0 = p1
       
       z1 = z2
       p1 = p2
    END DO
    if (ix==101) then
       write(*,*) "iteration did not converge"
    end if    
    iterate_z_given_pressure = z2

  END FUNCTION iterate_z_given_pressure


  REAL(kind=real_kind) FUNCTION moist_pressure_given_z(z)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: z
    REAL(8) :: aref, omegaref
    REAL(8) :: T0, constA, constB, constC, constH, scaledZ
    REAL(8) :: tau1, tau2, inttau1, inttau2
    REAL(8) :: rratio, inttermT,pwet,wvp
    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = a / bigX
    omegaref = omega * bigX
    
    T0 = 0.5d0 * (T0E + T0P)
    constA = 1.d0 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5d0 * (KK + 2.d0) * (T0E - T0P) / (T0E * T0P)
    constH = Rd * T0 / g
    
    scaledZ = z / (B * constH)
    
    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    
    inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
         + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)
    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
       rratio = 1.d0
    else
       rratio = (z + aref) / aref;
    end if
    
    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**KK &
         - KK / (KK + 2.d0) * (rratio * cos(lat))**(KK + 2.d0)

    !--------------------------------------------
    !    hydrostatic pressure
    !--------------------------------------------
    moist_pressure_given_z = p0 * exp(- g / Rd * (inttau1 - inttau2 * inttermT))
  END FUNCTION moist_pressure_given_z

  REAL(kind=real_kind) FUNCTION Tv_given_z(z)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: z
    REAL(8) :: aref, omegaref
    REAL(8) :: T0, constA, constB, constC, constH, scaledZ
    REAL(8) :: tau1, tau2, inttau1, inttau2
    REAL(8) :: rratio, inttermT
    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = a / bigX
    omegaref = omega * bigX
    
    T0 = 0.5d0 * (T0E + T0P)
    constA = 1.d0 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5d0 * (KK + 2.d0) * (T0E - T0P) / (T0E * T0P)
    constH = Rd * T0 / g
    
    scaledZ = z / (B * constH)
    
    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    
    inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
         + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)
    
    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
       rratio = 1.d0
    else
       rratio = (z + aref) / aref;
    end if
    
    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**KK &
         - KK / (KK + 2.d0) * (rratio * cos(lat))**(KK + 2.d0)
    
    !--------------------------------------------
    !    temperature
    !--------------------------------------------
    Tv_given_z = 1.d0 / (rratio**2 * (tau1 - tau2 * inttermT))
  END FUNCTION Tv_given_z

  SUBROUTINE uv_given_z(z,u,v)
    IMPLICIT NONE
    REAL(8), INTENT(IN)  :: z
    REAL(8), INTENT(OUT) :: u,v
    REAL(8) :: aref, omegaref
    REAL(8) :: T0, constH, constC, scaledZ, inttau2, rratio
    REAL(8) :: inttermU, bigU, rcoslat, omegarcoslat
    !------------------------------------------------
    !   Compute test case constants
    !------------------------------------------------
    aref = a / bigx
    omegaref = omega * bigx
    
    T0 = 0.5d0 * (T0E + T0P)
    
    constH = Rd * T0 / g
    
    constC = 0.5d0 * (KK + 2.d0) * (T0E - T0P) / (T0E * T0P)
    
    scaledZ = z / (B * constH)
    
    inttau2 = constC * z * exp(- scaledZ**2)
    
    ! radius ratio
    if (deep .eq. 0) then
       rratio = 1.d0
    else
       rratio = (z + aref) / aref;
    end if
    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    inttermU = (rratio * cos(lat))**(KK - 1.d0) - (rratio * cos(lat))**(KK + 1.d0)
    bigU = g / aref * KK * inttau2 * inttermU * Tv_given_z(z)
    if (deep .eq. 0) then
       rcoslat = aref * cos(lat)
    else
       rcoslat = (z + aref) * cos(lat)
    end if
    
    omegarcoslat = omegaref * rcoslat
    
    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0.d0
    
    !-----------------------------------------------------
    !   Add perturbation to the velocity field
    !-----------------------------------------------------
!    if (.false.) then !xxxx
    ! Exponential type
    if (pertt .eq. 0) then
       u = u + evaluate_exponential(z)
       
       ! Stream function type
    elseif (pertt .eq. 1) then
       u = u - 1.d0 / (2.d0 * dxepsilon) *                       &
            ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
            - evaluate_streamfunction(lon, lat - dxepsilon, z))
       
       v = v + 1.d0 / (2.d0 * dxepsilon * cos(lat)) *            &
            ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
            - evaluate_streamfunction(lon - dxepsilon, lat, z))
    end if
!    endif!xxx
  END SUBROUTINE uv_given_z

  !-----------------------------------------------------------------------
  !    Exponential perturbation function
  !-----------------------------------------------------------------------
  REAL(8) FUNCTION evaluate_exponential(z)
    REAL(8), INTENT(IN)  :: &
         z             ! Altitude (meters)
    
    REAL(8) :: greatcircler, perttaper
    
    ! Great circle distance
    greatcircler = 1.d0 / pertexpr &
         * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))
    
    ! Vertical tapering of stream function
    if (z < pertz) then
       perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
       perttaper = 0.d0
    end if
    
    ! Zonal velocity perturbation
    if (greatcircler < 1.d0) then
       evaluate_exponential = pertup * perttaper * exp(- greatcircler**2)
    else
       evaluate_exponential = 0.d0
    end if
    
  END FUNCTION evaluate_exponential
  
  !-----------------------------------------------------------------------
  !    Stream function perturbation function
  !-----------------------------------------------------------------------
  REAL(8) FUNCTION evaluate_streamfunction(z,lon_local,lat_local)
    
    REAL(8), INTENT(IN)  :: &
         lon_local, lat_local,&
         z             ! Altitude (meters)
    
    REAL(8) :: greatcircler, perttaper, cospert
    
    ! Great circle distance
    greatcircler = 1.d0 / pertr &
         * acos(sin(pertlat) * sin(lat_local) + cos(pertlat) * cos(lat_local) * cos(lon_local - pertlon))
    
    ! Vertical tapering of stream function
    if (z < pertz) then
       perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
       perttaper = 0.d0
    end if
    
    ! Horizontal tapering of stream function
    if (greatcircler .lt. 1.d0) then
       cospert = cos(0.5d0 * pi * greatcircler)
    else
       cospert = 0.d0
    end if
    
    evaluate_streamfunction = &
         (- pertu0 * pertr * perttaper * cospert**4)
    
  END FUNCTION evaluate_streamfunction
  
  subroutine umjs_baroclinic(elem,hybrid,hvcoord,nets,nete, fvm, ldry)
    use element_mod, only : element_t, timelevels
    use fvm_control_volume_mod, only : fvm_struct
    use hybrid_mod, only : hybrid_t
    use hybvcoord_mod, only : hvcoord_t 
    use kinds, only : real_kind  
    use dimensions_mod, only : nlev,np, qsize,nc,ntrac,ldry_mass_vertical_coordinates
    use parallel_mod, only : abortmp
    use chemistry, only : initial_value
    !
    ! for tracers
    !
    use cube_mod, only : rotate_grid
    use jw, only : tracer_q1_q2,tracer_q3,tracer_slotted_cylinder
    implicit none
    
    type(element_t), intent(inout) :: elem(:)
    type(fvm_struct), optional, intent(inout) :: fvm(:)
    
    type (hvcoord_t)            :: hvcoord
    type (hybrid_t), intent(in) :: hybrid
    integer, intent(in)         :: nets,nete
    logical, intent(in)         :: ldry
    
    !   local
    integer                :: i,j,k,ie,idex,t,jgw
    
    real(kind=real_kind) :: u_out,v_out,T_out,p_out,phis_out,ps_out,q_out, p_in
    integer              :: dry                                           ! 0: dry, 1: moist
    real(kind=real_kind) :: zdummy, rho_dummy, thetav_dummy, Cl, Cl2 ! dummy vars
    real(kind=real_kind) :: eta_c, ptop
    real(kind=real_kind) :: dp_dry, dp_wet, psdry, pwet, pdry, z, pressure, temp, rho, qv, z1,z2

    real(kind=real_kind) :: zdry(nlev),wvp, pdry_half(nlev+1), pwet_half(nlev+1)

    ! !  print*,'Creating Initial Data for Test 70'
    if (ldry) then
       dry     = 0    ! 0: dry, 1:moist 
    else
       dry     = 1    ! 0: dry, 1:moist 
    end if


!    call check_balances(elem,ptop)
!    write(*,*) "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx"
    ptop = hvcoord%hyai(1)*hvcoord%ps0
    ztop = iterate_z_given_pressure(ptop,.false.,ptop)
!    write(*,*) "ztop",ztop
    if (.not.ldry_mass_vertical_coordinates) then
       do ie=nets,nete
          do j=1,np
             do i=1,np             
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                do k=1,nlev
                   p_in =  hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*p0
                   call baroclinic_wave_test(dry, p_in,ptop,zdummy,&
                        u_out,v_out,T_out,thetav_dummy,phis_out,ps_out,rho_dummy,q_out)
                   elem(ie)%state%v(i,j,1,k,:) = u_out
                   elem(ie)%state%v(i,j,2,k,:) = v_out
                   elem(ie)%state%T(i,j  ,k,:) = T_out
                   if (qsize>=1) elem(ie)%state%Q(i,j,k,1) = q_out
                enddo
                elem(ie)%state%phis(i,j ) = phis_out
                elem(ie)%state%pswet(i,j) = ps_out
                elem(ie)%state%lnps(i,j,:) = LOG(ps_out)
             end do
          enddo
       enddo
    else
!       write(*,*) "model top is",ztop
!       call check_balances_dry(elem,hvcoord,ztop,ptop)
       do ie=nets,nete
          do j=1,np
             do i=1,np             
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat

                !
                ! compute water vapor pressure at the surface
                !               
                wvp = weight_of_water_vapor_given_z(0.0D0,ptop)
                elem(ie)%state%ps(i,j,:) = p0-wvp

!                write(*,*) "water vapor pressure",wvp
                do k=1,nlev
                   p_in =  hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*elem(ie)%state%ps(i,j,1)
                   call baroclinic_wave_test(dry, p_in,ptop,zdummy,&
                        u_out,v_out,T_out,thetav_dummy,phis_out,ps_out,rho_dummy,q_out)
!                   zdry(k) = zdummy

!                   if (qsize>=1) elem(ie)%state%Q(i,j,k,1) = q_out
                   elem(ie)%state%v(i,j,1,k,:) = u_out
                   elem(ie)%state%v(i,j,2,k,:) = v_out
                   elem(ie)%state%T(i,j  ,k,:) = T_out
                end do
                do k=1,nlev+1
                   pdry_half(k) =  hvcoord%hyai(k)*hvcoord%ps0 + hvcoord%hybi(k)*p0
                   zdummy = iterate_z_given_pressure(pdry_half(k),ldry_mass_vertical_coordinates,ptop) !Find height of pressure surface
                   pwet_half(k) = moist_pressure_given_z(zdummy)
                end do
                do k=1,nlev
                   if (qsize>=1)&
                        elem(ie)%state%Q(i,j,k,1) = &
                        ((pwet_half(k+1)-pwet_half(k))/(pdry_half(k+1)-pdry_half(k)))-1.0D0
                end do


!                elem(ie)%state%phis(i,j ) = 0.0D0!phis_out
                elem(ie)%state%pswet(i,j) = p0
!                elem(ie)%state%ps(i,j,:) = psdry
!                elem(ie)%state%lnps(i,j,:) = LOG(ps_out)
             end do
          enddo
       enddo
    end if

    if (qsize>=3) then
       do ie=nets,nete
          do j=1,np
             do i=1,np
                !
                ! terminator chemistry initial conditions
                !
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                call initial_value( lat, lon, Cl, Cl2 )
                elem(ie)%state%Q(i,j,:,2) = Cl
                elem(ie)%state%Q(i,j,:,3:qsize) = Cl2
             end do
          end do
       end do
    end if
    !
    ! inert DCMIP tracers
    !   
    if (qsize>=4) then
       idex=4
       eta_c = 0.6
       do ie=nets,nete
          do j=1,np
             do i=1,np
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                do k=1,nlev
                   elem(ie)%state%Q(i,j,k,idex) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
                enddo
             enddo
          enddo
       enddo
    endif
    
    if (qsize>=5) then
       idex=5
       eta_c = 1
       do ie=nets,nete
          do j=1,np
             do i=1,np
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                do k=1,nlev
                   elem(ie)%state%Q(i,j,k,idex) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
                enddo
             enddo
          enddo
       enddo
    endif

    if (qsize>=6) then
       idex=6
       do ie=nets,nete
          do j=1,np
             do i=1,np
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                do k=1,nlev
                   elem(ie)%state%Q(i,j,k,idex) = tracer_q3(lon,lat,hvcoord%etam(k),rotate_grid)
                enddo
             enddo
          enddo
       enddo
    endif
    
    if (qsize>=7) then
       idex=7
       eta_c=1.0D0
       do ie=nets,nete
          do j=1,np
             do i=1,np
                lon = elem(ie)%spherep(i,j)%lon
                lat = elem(ie)%spherep(i,j)%lat
                do k=1,nlev
                   elem(ie)%state%Q(i,j,k,idex) = tracer_slotted_cylinder(lon,lat,hvcoord%etam(k),rotate_grid,eta_c)
                enddo
             end do
          end do
       enddo
    endif
    


    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! fvm tracers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (present(fvm)) then
       !
       ! CSLAM tracers
       !
       if (ntrac>=1) then
          do ie=nets,nete
             do j=1,nc
                do i=1,nc
                   lon = fvm(ie)%centersphere(i,j)%lon
                   lat = fvm(ie)%centersphere(i,j)%lat
                   do k=1,nlev
                      p_in =  hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*p0
                      call baroclinic_wave_test(dry, p_in,ptop,zdummy,&
                           u_out,v_out,T_out,thetav_dummy,phis_out,ps_out,rho_dummy,q_out)
                      fvm(ie)%c(i,j,k,1,:) = q_out
!                      if (ie==1.and.i==1.and.j==1.and.k==1) write(*,*) "not initialized correctly !!!!"
                   enddo
                enddo
             enddo
          enddo
       endif

       if (ntrac>=3) then
          do ie=nets,nete
             do j=1,nc
                do i=1,nc
                   lon = fvm(ie)%centersphere(i,j)%lon
                   lat = fvm(ie)%centersphere(i,j)%lat
                   call initial_value( lat, lon, Cl, Cl2 )
                   fvm(ie)%c(i,j,:,2      ,:) = Cl
                   fvm(ie)%c(i,j,:,3:ntrac,:) = Cl2
                enddo
             enddo
          enddo
       endif

    !
    ! inert DCMIP tracers
    !   
    if (ntrac>=4) then
       idex=4
       eta_c = 0.6
       do ie=nets,nete
          do j=1,nc
             do i=1,nc
                lon = fvm(ie)%centersphere(i,j)%lon
                lat = fvm(ie)%centersphere(i,j)%lat
                do k=1,nlev
                  fvm(ie)%c(i,j,k,idex,:) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
                enddo
             enddo
          enddo
       enddo
    endif
    
    if (ntrac>=5) then
       idex=5
       eta_c = 1
       do ie=nets,nete
          do j=1,nc
             do i=1,nc
                lon = fvm(ie)%centersphere(i,j)%lon
                lat = fvm(ie)%centersphere(i,j)%lat
                do k=1,nlev
                  fvm(ie)%c(i,j,k,idex,:) = tracer_q1_q2(lon,lat,hvcoord%etam(k),rotate_grid, eta_c)
                enddo
             enddo
          enddo
       enddo
    endif

    if (ntrac>=6) then
       idex=6
       do ie=nets,nete
          do j=1,nc
             do i=1,nc
                lon = fvm(ie)%centersphere(i,j)%lon
                lat = fvm(ie)%centersphere(i,j)%lat
                do k=1,nlev
                  fvm(ie)%c(i,j,k,idex,:) = tracer_q3(lon,lat,hvcoord%etam(k),rotate_grid)
                enddo
             enddo
          enddo
       enddo
    endif
    
    if (ntrac>=7) then
       idex=7
       eta_c=1.0D0
       do ie=nets,nete
          do j=1,nc
             do i=1,nc
                lon = fvm(ie)%centersphere(i,j)%lon
                lat = fvm(ie)%centersphere(i,j)%lat
                do k=1,nlev                   
                   fvm(ie)%c(i,j,k,idex,:) = tracer_slotted_cylinder(lon,lat,hvcoord%etam(k),rotate_grid,eta_c)
                enddo
             end do
          end do
       enddo
    endif
 endif
  end subroutine umjs_baroclinic

  subroutine prim_printstate_par_terminator(elem, tl,hybrid,hvcoord,nets,nete, fvm,cly_mass_init,lcompute_cly_mass_init)
    ! ------------------------------
    use kinds, only : real_kind, iulog
    ! ------------------------------
    use dimensions_mod, only : nlev, np, nc, qsize_d, qsize, nelemd, ntrac, ntrac_d
    ! ------------------------------
    use parallel_mod, only :  iam, ordered, parallel_t, syncmp
    use parallel_mod, only: global_shared_buf, global_shared_sum
    ! ------------------------------
    use global_norms_mod, only: wrap_repro_sum
    ! ------------------------------
    use hybrid_mod, only : hybrid_t
    ! ------------------------------
    use time_mod, only : tstep, secpday, timelevel_t, TimeLevel_Qdp, time_at
    use fvm_control_volume_mod, only : n0_fvm
    ! ------------------------------
!    use control_mod, only : integration, test_case, runtype, moisture, &
!         tstep_type,energy_fixer, qsplit, ftype, use_cpstar, rsplit
    ! ------------------------------
    use hybvcoord_mod, only : hvcoord_t 
    ! ------------------------------
    use global_norms_mod, only : global_integral, linf_snorm, l1_snorm, l2_snorm
    ! ------------------------------
    use element_mod, only : element_t
    ! ------------------------------
    use fvm_control_volume_mod, only : fvm_struct, n0_fvm
    ! ------------------------------
    use reduction_mod, only : parallelmax,parallelmin
    ! ------------------------------


    type (element_t), intent(in) :: elem(:)
    type (TimeLevel_t), target, intent(in) :: tl
    type (hybrid_t),intent(in)     :: hybrid
    type (hvcoord_t), intent(in)   :: hvcoord
    type(fvm_struct), optional, intent(in) :: fvm(:)
    integer,intent(in)             :: nets,nete
    character(len=*), parameter    :: fstub = "state_norms"
    integer	                   :: simday
!    type(parallel_t)               :: par

    real (kind=real_kind), intent(inout) :: cly_mass_init(2)
    logical              , intent(in) :: lcompute_cly_mass_init

    real (kind=real_kind) :: l1,l2,linf,day, cly_mass, cly_mass_fvm, linf_fvm, l2_cly_fvm

    real (kind=real_kind)  :: Cly(np,np,nets:nete), Cly_exact(np,np,nets:nete)
    real (kind=real_kind)  :: Cly_fvm(nc,nc,nets:nete), Cly_exact_fvm(nc,nc,nets:nete), Cly_constant
    real (kind=real_kind)  :: dp_fvm_init(nc,nc,nlev)
    real (kind=real_kind)  :: l2denom_Cly_fvm(nc,nc,nets:nete), l2nom_Cly_fvm(nc,nc,nets:nete)
    real (kind=real_kind)  :: tmp(nets:nete),dp3d_tmp(1:np,1:np,nlev),dp3d_init(1:np,1:np,nlev)

    integer               :: n0,i,j,k,ie,npts

    cly_constant = 4.0D-6
    n0=tl%n0
    day=Time_at(tl%nstep)/(24*3600)
    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       print *, "Computing terminator chemistry diagnostics at time step ",tl%nstep
       print *,"===================================================================="
    end if

    if (qsize>0) then
       if (qsize<3) call abortmp("Error: cannot do terminator chemistry diagnostics with less then 3 tracers")
       npts=np

       do ie=nets,nete
         do k=1,nlev
            dp3d_tmp(:,:,k)=&
                 ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                 ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(:,:,n0)
         end do
         
         do j=1,np
            do i=1,np
               Cly      (i,j,ie) = SUM(dp3d_tmp(i,j,:)*&
                    (elem(ie)%state%Q(i,j,:,2)+2.0D0*elem(ie)%state%Q(i,j,:,3)))!/&
            end do
         end do
      end do

      if (lcompute_cly_mass_init) then
         cly_mass_init(1)  = global_integral(elem, Cly(:,:,nets:nete),hybrid,npts,nets,nete)
      end if
      cly_mass_exact = cly_mass_init(1)
       !
       ! Cly holds column integrated Cly mass: compute global Cly mass integrals
       !
       cly_mass      =global_integral(elem, Cly      (:,:,nets:nete),hybrid,npts,nets,nete)
       !
       ! Convert Cly to mixing ratio
       !
       do ie=nets,nete
          do k=1,nlev
             dp3d_tmp(:,:,k)=&
                  ( hvcoord%hyai(k+1) - hvcoord%hyai(k) )*hvcoord%ps0 + &
                  ( hvcoord%hybi(k+1) - hvcoord%hybi(k) )*elem(ie)%state%ps(:,:,n0)
          end do
          do j=1,np
             do i=1,np
                Cly      (i,j,ie) = Cly      (i,j,ie)/SUM(dp3d_tmp(i,j,:))
                Cly_exact(i,j,ie) = Cly_constant
             end do
          end do
       end do
! parallel read file with state vector in unformatted blocks as written above                                                 
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER                                                                                                             
#endif
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER                                                                                                             
#endif
       npts=np          
       l2   = l2_snorm  (elem,Cly(:,:,nets:nete), Cly_exact(:,:,nets:nete),hybrid,npts,nets,nete)
       linf = linf_snorm(     Cly(:,:,nets:nete), Cly_exact(:,:,nets:nete),hybrid,npts,nets,nete)          
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
        print *, "Relative Cly mass change=",day,(cly_mass-cly_mass_exact)/cly_mass_exact
!          print *, "L1  =",day,l1
          print *, "L2  =",day,l2
          print *, "Linf=",day,linf
       end if
    end if
#if (defined HORIZ_OPENMP)
    !$OMP BARRIER                                                                                                             
#endif

    if (ntrac>0) then
       if (ntrac<3) call abortmp("Error: cannot do terminator chenistry diagnostics with less then 3 fvm tracers")
       do ie=nets,nete
          do j=1,nc
             do i=1,nc
                !
                ! cly_fvm holds Cly mass
                !
                cly_fvm(i,j,ie) = SUM(fvm(ie)%dp_fvm(i,j,:,n0_fvm)*(&
                     fvm(ie)%c(i,j,:,2,n0_fvm)+2.0D0*fvm(ie)%c(i,j,:,3,n0_fvm)))
                !
                ! norms are based on Cly mixing ratio
                !
                l2nom_cly_fvm(i,j,ie)   = (cly_fvm(i,j,ie)/SUM(fvm(ie)%dp_fvm(i,j,:,n0_fvm))-&
                     cly_constant)**2
                l2denom_cly_fvm(i,j,ie) = cly_constant**2 
             end do
          end do
          global_shared_buf(ie,1) = SUM(cly_fvm(:,:,ie)*fvm(ie)%area_sphere(1:nc,1:nc))
          global_shared_buf(ie,2) = SUM(l2nom_cly_fvm  (:,:,ie)*fvm(ie)%area_sphere(1:nc,1:nc))
          global_shared_buf(ie,3) = SUM(l2denom_cly_fvm(:,:,ie)*fvm(ie)%area_sphere(1:nc,1:nc))
          !
          ! convert to mixing ratio
          ! 
          do j=1,nc
             do i=1,nc
                cly_fvm(i,j,ie)=cly_fvm(i,j,ie)/SUM(fvm(ie)%dp_fvm(i,j,:,n0_fvm))
             end do
          end do
          tmp(ie) = MAXVAL(ABS(cly_fvm(:,:,ie)-cly_constant)) !for linf_fvm
       end do
       call wrap_repro_sum(nvars=3, comm=hybrid%par%comm)
       cly_mass_fvm = global_shared_sum(1)                
       if (lcompute_cly_mass_init) then
          cly_mass_exact_fvm = cly_mass_fvm
       else
          cly_mass_exact = cly_mass_init(2)
       end if

       l2_cly_fvm = SQRT(global_shared_sum(2))/SQRT(global_shared_sum(3))
       linf_fvm  = ParallelMax(tmp,hybrid)/cly_constant
      
       if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
          print *, "Relative Cly fvm mass change =",day,(cly_mass_fvm-cly_mass_exact_fvm)/cly_mass_exact_fvm
          print *,"l2   fvm =",day,l2_cly_fvm
          print *,"linf fvm =",day,linf_fvm
       end if
    endif
    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       print *, " "
       print *,"===================================================================="
    end if
  end subroutine prim_printstate_par_terminator

  subroutine check_balances(elem,ptop)
    use element_mod, only : element_t
    implicit none
    type(element_t), intent(in) :: elem(:)
    real (kind=real_kind)  :: lambda(-1:1), theta(-1:1)
    real (kind=real_kind)  :: dlambda, dtheta

    real (kind=real_kind), dimension(-1:1,-1:1)  :: u,v,temp,z,q, rho
    real (kind=real_kind), dimension(-1:1,-1:1)  :: thetav, phis, ps, p_out, p_in, utmp, vtmp

    real (kind=real_kind)  :: grad_p(2), grad_phi(2), vorticity, vort_term(2), grad_usq(2),ptop

    integer :: i,j,ie,ic,jc

    ie=572
    ic=2
    jc=2

    lambda(0)  = elem(ie)%spherep(ic,jc)%lon
    theta(0)   = elem(ie)%spherep(ic,jc)%lat
    dlambda = 0.5D0*pi/180D0
    dtheta  = 0.5D0*pi/180D0

    p_in = 85953.476652503072D0

    lambda(-1) = lambda(0)-dlambda;lambda(1) = lambda(0)+dlambda;
    theta (-1) = theta (0)-dtheta ;theta (1) = theta (0)+dtheta;
    
    do j=-1,1
       do i=-1,1
          lon = lambda(i)
          lat = theta (j)
          call baroclinic_wave_test(1, p_in(i,j),ptop,z(i,j),&
               u(i,j),v(i,j),temp(i,j),thetav(i,j),phis(i,j),ps(i,j),rho(i,j),q(i,j))
       end do
    end do

    grad_p(1) = (1.0D0/rho(0,0))*(1.0D0/cos(theta(0)))*(p_in(1,0)-p_in(-1,0))/(2.0D0*dlambda)
    grad_p(2) = (1.0D0/rho(0,0))*     (p_in(0,1)-p_in( 0,-1))/(2.0D0*dtheta)

    grad_p=grad_p/rearth

    grad_phi(1) = g*(1.0D0/cos(theta(0)))*(z(1,0)-z(-1,0))/(2.0D0*dlambda)
    grad_phi(2) =                       g*(z(0,1)-z( 0,-1))/(2.0D0*dtheta)


    grad_phi=grad_phi/rearth

    utmp = 0.5D0*(u*u+v*v)

    grad_usq(1) =   (1.0D0/cos(theta(0)))*(utmp(1,0)-utmp(-1,0))/(2.0D0*dlambda)
    grad_usq(2) =                         (utmp(0,1)-utmp( 0,-1))/(2.0D0*dtheta)

    grad_usq=grad_usq/rearth


    vorticity = (1.0D0/cos(theta(0)))*((v(1,0)-v(-1,0))/(2.0D0*dlambda)-&
                                       (cos(theta(1))*u(0,1)-cos(theta(-1))*u(0,-1))/(2.0D0*dtheta))
    vorticity = vorticity/rearth

    vort_term(1) = -v(0,0)*(vorticity+elem(ie)%fcor(ic,jc))
    vort_term(2) =  u(0,0)*(vorticity+elem(ie)%fcor(ic,jc))

    write(*,*) "grad_p(1)",grad_p(1)
    write(*,*) "grad_p(2)",grad_p(2)
    write(*,*) "grad_phi(1)",grad_phi(1)
    write(*,*) "grad_phi(2)",grad_phi(2)
    write(*,*) "grad_usq(1)",grad_usq(1)
    write(*,*) "grad_usq(2)",grad_usq(2)
    write(*,*) "vtemp1",grad_phi(1)+grad_usq(1)
    write(*,*) "vtemp2",grad_phi(2)+grad_usq(2)
    write(*,*) "vorticity  ",vorticity
    write(*,*) "vort_term1",vort_term(1)
    write(*,*) "vort_term2",vort_term(2)
    write(*,*) "=========="
    write(*,*) "Total1",vort_term(1)+grad_phi(1)+grad_usq(1)
    write(*,*) "Total2",vort_term(2)+grad_phi(2)+grad_usq(2)

  end subroutine check_balances

  subroutine check_balances_dry(elem,hvcoord,ztop,ptop)
    use element_mod, only : element_t
    use hybvcoord_mod, only : hvcoord_t 
    implicit none
    type(element_t), intent(in) :: elem(:)
    type (hvcoord_t)            :: hvcoord
    real (kind=real_kind)  :: lambda(-1:1), theta(-1:1)
    real (kind=real_kind)  :: dlambda, dtheta, ztop,ptop

    real (kind=real_kind), dimension(-1:1,-1:1)  :: u,v,temp,z,q, rho,pwet
    real (kind=real_kind), dimension(-1:1,-1:1)  :: thetav, phis, ps, p_out, p_in, utmp, vtmp

    real (kind=real_kind)  :: grad_p(2), grad_phi(2), vorticity, vort_term(2), grad_usq(2),wvp

    integer :: i,j,ie,ic,jc,k

    ie=572
    ic=2
    jc=2

    lambda(0)  = elem(ie)%spherep(ic,jc)%lon
    theta(0)   = elem(ie)%spherep(ic,jc)%lat
    dlambda = 0.5D0*pi/180D0
    dtheta  = 0.5D0*pi/180D0

    lambda(-1) = lambda(0)-dlambda;lambda(1) = lambda(0)+dlambda;
    theta (-1) = theta (0)-dtheta ;theta (1) = theta (0)+dtheta;

    k=24
    do j=-1,1
       do i=-1,1
          lon = lambda(i)
          lat = theta (j)          
          wvp = weight_of_water_vapor_given_z(0.0D0,ptop)
          p_in(i,j) =  hvcoord%hyam(k)*hvcoord%ps0 + hvcoord%hybm(k)*(p0-wvp)
          write(*,*) "p_in",p0,p_in(i,j)
       end do
    end do


    
    do j=-1,1
       do i=-1,1
          lon = lambda(i)
          lat = theta (j)
          call baroclinic_wave_test(1, p_in(i,j),ptop,z(i,j),&
               u(i,j),v(i,j),temp(i,j),thetav(i,j),phis(i,j),ps(i,j),rho(i,j),q(i,j))
          pwet(i,j)=p_in(i,j)!*(1+q(i,j)*1.608D0)
       end do
    end do



    grad_p(1) = (1.0D0/rho(0,0))*(1.0D0/cos(theta(0)))*(pwet(1,0)-pwet(-1,0))/(2.0D0*dlambda)
    grad_p(2) = (1.0D0/rho(0,0))*     (pwet(0,1)-pwet( 0,-1))/(2.0D0*dtheta)

    grad_p=grad_p/rearth

    grad_phi(1) = g*(1.0D0/cos(theta(0)))*(z(1,0)-z(-1,0))/(2.0D0*dlambda)
    grad_phi(2) =                       g*(z(0,1)-z( 0,-1))/(2.0D0*dtheta)


    grad_phi=grad_phi/rearth

    utmp = 0.5D0*(u*u+v*v)

    grad_usq(1) =   (1.0D0/cos(theta(0)))*(utmp(1,0)-utmp(-1,0))/(2.0D0*dlambda)
    grad_usq(2) =                         (utmp(0,1)-utmp( 0,-1))/(2.0D0*dtheta)

    grad_usq=grad_usq/rearth


    vorticity = (1.0D0/cos(theta(0)))*((v(1,0)-v(-1,0))/(2.0D0*dlambda)-&
                                       (cos(theta(1))*u(0,1)-cos(theta(-1))*u(0,-1))/(2.0D0*dtheta))
    vorticity = vorticity/rearth

    vort_term(1) = -v(0,0)*(vorticity+elem(ie)%fcor(ic,jc))
    vort_term(2) =  u(0,0)*(vorticity+elem(ie)%fcor(ic,jc))

    write(*,*) "grad_p(1)",grad_p(1)
    write(*,*) "grad_p(2)",grad_p(2)
    write(*,*) "grad_phi(1)",grad_phi(1)
    write(*,*) "grad_phi(2)",grad_phi(2)
    write(*,*) "grad_usq(1)",grad_usq(1)
    write(*,*) "grad_usq(2)",grad_usq(2)
    write(*,*) "vtemp1",grad_phi(1)+grad_usq(1)
    write(*,*) "vtemp2",grad_phi(2)+grad_usq(2)
    write(*,*) "vorticity  ",vorticity
    write(*,*) "vort_term1",vort_term(1)
    write(*,*) "vort_term2",vort_term(2)
    write(*,*) "=========="
    write(*,*) "Total1",vort_term(1)+grad_phi(1)+grad_usq(1)
    write(*,*) "Total2",vort_term(2)+grad_phi(2)+grad_usq(2)

  end subroutine check_balances_dry


END MODULE umjs_baroclinic_wave_mod
