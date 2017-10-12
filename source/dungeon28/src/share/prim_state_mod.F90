#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module prim_state_mod
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
  use hybrid_mod, only : hybrid_t, PrintHybrid
  ! ------------------------------
  use physical_constants, only : p0,Cp,g
  ! ------------------------------
  use time_mod, only : tstep, secpday, timelevel_t, TimeLevel_Qdp, time_at
  ! ------------------------------
  use control_mod, only : integration, test_case, runtype, moisture, &
       tstep_type,energy_fixer, qsplit, ftype, use_cpstar, rsplit
  ! ------------------------------
  use hybvcoord_mod, only : hvcoord_t 
  ! ------------------------------
  use global_norms_mod, only : global_integral, linf_snorm, l1_snorm, l2_snorm
  ! ------------------------------
  use element_mod, only : element_t
  ! ------------------------------
  use fvm_control_volume_mod, only : fvm_struct
  use spelt_mod, only : spelt_struct
  ! ------------------------------
  use viscosity_mod, only : compute_zeta_C0
  ! ------------------------------
  use reduction_mod, only : parallelmax,parallelmin
  ! ------------------------------
#if defined(_REFSOLNW) || defined(_REFSOLNR) || defined(_REFSOLN)
  use ref_state_mod, only : ref_state_read, ref_state_write
#endif

implicit none
private
  character(len=*), private, parameter :: massfname = "mass.out"

  integer, public :: naccum


  public :: prim_printstate
  public :: prim_printstate_par
!  public :: prim_printstate_par_terminator
  public :: prim_printstate_init
contains
!=======================================================================================================! 


  subroutine prim_printstate_init(par)
    type (parallel_t) :: par

    real (kind=real_kind) :: time
    integer               :: c0

    if (par%masterproc) then
       time=0.0D0
       c0  =0.0D0

#ifndef CAM
#ifndef _BGL
       open(unit=10,file=massfname,form="formatted",status="replace")
!       write(10,*)time,c0
       close(10)
#endif
#endif
    end if

  end subroutine prim_printstate_init
!=======================================================================================================! 

  subroutine prim_printstate(elem, tl,hybrid,hvcoord,nets,nete, fvm)
    use physical_constants    , only : dd_pi
    use control_mod           , only : tracer_transport_type
    use fvm_control_volume_mod, only : n0_fvm, np1_fvm
    use dimensions_mod        , only : qsize_condensate_loading, ldry_mass_vertical_coordinates

    type (element_t), intent(inout) :: elem(:)
    
#if defined(_SPELT)
      type(spelt_struct), optional, intent(in) :: fvm(:)
#else
      type(fvm_struct), optional, intent(inout) :: fvm(:)
#endif
    type (TimeLevel_t), target, intent(in) :: tl
    type (hybrid_t),intent(in)     :: hybrid
    type (hvcoord_t), intent(in)   :: hvcoord
    integer,intent(in)             :: nets,nete

    ! Local variables...
    integer :: i,j,k,ie
    integer,parameter  :: type=ORDERED

    real (kind=real_kind)  :: Mass2,Mass,Mass2_moist, Mass_moist
    real (kind=real_kind),save  :: time0

    real (kind=real_kind)  :: tmp(np,np,nets:nete)
    real (kind=real_kind),allocatable  :: tmp3d(:,:,:,:)
    real (kind=real_kind)  :: tmp1(nets:nete)
    real (kind=real_kind)  :: ps(np,np)
    real (kind=real_kind)  :: dp(np,np)
    !    real (kind=real_kind)  :: E(np,np)

    real (kind=real_kind) :: umin_local(nets:nete), umax_local(nets:nete), usum_local(nets:nete), &
         vmin_local(nets:nete), vmax_local(nets:nete), vsum_local(nets:nete), &
         tmin_local(nets:nete), tmax_local(nets:nete), tsum_local(nets:nete), &
         psmin_local(nets:nete),psmax_local(nets:nete),pssum_local(nets:nete), &
         fumin_local(nets:nete),fumax_local(nets:nete),fusum_local(nets:nete), &
         fvmin_local(nets:nete),fvmax_local(nets:nete),fvsum_local(nets:nete), &
         ftmin_local(nets:nete),ftmax_local(nets:nete),ftsum_local(nets:nete), &
         fqmin_local(nets:nete),fqmax_local(nets:nete),fqsum_local(nets:nete), &
         omegamin_local(nets:nete),omegamax_local(nets:nete),omegasum_local(nets:nete),&
         dpmin_local(nets:nete), dpmax_local(nets:nete), dpsum_local(nets:nete)


    real (kind=real_kind) :: umin_p, vmin_p, tmin_p, qvmin_p(qsize_d), cmin(ntrac_d),&
         psmin_p, dpmin_p


    real (kind=real_kind) :: umax_p, vmax_p, tmax_p, qvmax_p(qsize_d), cmax(ntrac_d),&
         psmax_p, dpmax_p

    real (kind=real_kind) :: usum_p, vsum_p, tsum_p, qvsum_p(qsize_d), csum(ntrac_d),&
         pssum_p, dpsum_p, relative_mass_change(ntrac+3),relative_mass_change_se(qsize)

    !
    ! for fvm diagnostics
    !
    real (kind=real_kind) :: psc_mass, psc_min, psc_max,dp_fvm_mass, dp_fvm_min, dp_fvm_max 


    real (kind=real_kind) :: fusum_p, fvsum_p, ftsum_p, fqsum_p
    real (kind=real_kind) :: fumin_p, fvmin_p, ftmin_p, fqmin_p
    real (kind=real_kind) :: fumax_p, fvmax_p, ftmax_p, fqmax_p
    real (kind=real_kind) :: omegamax, omegamin, omegasum


    real(kind=real_kind) :: vsum_t(1), relvort
    real(kind=real_kind) :: v1, v2, vco(np,np,2,nlev)

    real (kind=real_kind) :: time, time2,time1, scale, dt, dt_split
    real (kind=real_kind) :: KEvert,IEvert,T1,T2,T2_s,T2_m,S1,S2,S1_wet
    real (kind=real_kind) :: KEhorz,IEhorz,IEhorz_wet,IEvert_wet
    real (kind=real_kind) :: ddt_tot,ddt_diss
    integer               :: n0, nm1, pnm1, np1, n0_qdp, np1_qdp
    integer               :: npts,n,q
    
    if (hybrid%masterthread) then 
       write(iulog,*) "nstep=",tl%nstep," time=",Time_at(tl%nstep)/(24*3600)," [day]"
    end if
    if (.not. present(fvm) .and. ntrac>0) then
       print *,'ERROR: prim_state_mod.F90: optional fvm argument required if ntrac>0'
    endif

    if (tstep_type==0) return !BUG model crashes in 

    ! dynamics timelevels
    n0=tl%n0
    nm1=tl%nm1
    np1=tl%np1

    call TimeLevel_Qdp( tl, qsplit, n0_qdp, np1_qdp)

    ! forcing timelevel
#ifdef CAM
    pnm1 = 1
#else
    pnm1=tl%nm1  
#endif

    dt=tstep*qsplit
    if (rsplit>0) dt = tstep*qsplit*rsplit  ! vertical REMAP timestep 
    !
    !   dynamics variables in n0 are at time =  'time' 
    !
    !   Diagnostics computed a 4 different levels during one compute REMAP step
    ! in RK code:
    !   E(:,:,2)-E(:,:,1)   change due to dynamics step  from time-dt to time
    !   E(:,:,3)-E(:,:,2)   change due to energy fixer   
    !   E(:,:,1)-E(:,:,4)   impact of forcing
    !
    ! Dissipation rates were computed during the first dynamics timstep, and represent
    ! the change going from 'time-dt' to 'time-dt+tstep' (one dynamics step)
    !
    time=tl%nstep*tstep
    time2 = time 
    time1 = time - dt


    vsum_t(1) = 0.0D0

    ! npts = np
    npts=SIZE(elem(1)%state%lnps(:,:,n0),1)
    
    do q=1,qsize
      do ie=nets,nete
        tmp1(ie) = MINVAL(elem(ie)%state%Q(:,:,:,q))
      enddo
      qvmin_p(q) = ParallelMin(tmp1,hybrid)
      do ie=nets,nete
        tmp1(ie) = MAXVAL(elem(ie)%state%Q(:,:,:,q))
      enddo

      qvmax_p(q) = ParallelMax(tmp1,hybrid)
      
      do ie=nets,nete
        do j=1,np
          do i=1,np
            tmp(i,j,ie)=SUM(elem(ie)%state%Qdp(i,j,:,q,n0_qdp))
          end do
        enddo
      enddo
      !
      ! qvsum is total mass of tracer
      !
      qvsum_p(q) = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    end do
    if (tl%nstep==0) then
      relative_mass_change_se(:) = 0.0D0
      elem(nets)%derived%mass(1:qsize)   = qvsum_p(1:qsize)
    else
      do q=1,qsize
        if (ABS(elem(nets)%derived%mass(q))<1.0D-12) then
          relative_mass_change_se(q) = qvsum_p(q) - elem(nets)%derived%mass(q)
        else
          relative_mass_change_se(q) = (qvsum_p(q) - elem(nets)%derived%mass(q))/elem(nets)%derived%mass(q)
        end if
      end do
    end if

    !

    do ie=nets,nete
      tmp(:,:,ie)=elem(ie)%state%ps(:,:,n0)
      !       tmp(:,:,ie)=EXP(elem(ie)%state%lnps(:,:,n0))
      
      
      !======================================================  
      umax_local(ie)    = MAXVAL(elem(ie)%state%v(:,:,1,:,n0))
      vmax_local(ie)    = MAXVAL(elem(ie)%state%v(:,:,2,:,n0))
      
      fumax_local(ie)    = MAXVAL(elem(ie)%derived%FM(:,:,1,:,pnm1))
      fvmax_local(ie)    = MAXVAL(elem(ie)%derived%FM(:,:,2,:,pnm1))
      
      tmax_local(ie)    = MAXVAL(elem(ie)%state%T(:,:,:,n0))
      
      if (rsplit>0) &
           dpmax_local(ie)    = MAXVAL(elem(ie)%state%dp3d(:,:,:,n0))
      
      psmax_local(ie) = MAXVAL(tmp(:,:,ie))
      ftmax_local(ie)    = MAXVAL(elem(ie)%derived%FT(:,:,:,pnm1))
      fqmax_local(ie)    = MAXVAL(elem(ie)%derived%FQ(:,:,:,1,pnm1))
      omegamax_local(ie)    = MAXVAL(elem(ie)%derived%Omega(:,:,:))
      !======================================================
      
      umin_local(ie)    = MINVAL(elem(ie)%state%v(:,:,1,:,n0))
      vmin_local(ie)    = MINVAL(elem(ie)%state%v(:,:,2,:,n0))
      
      Fumin_local(ie)    = MINVAL(elem(ie)%derived%FM(:,:,1,:,pnm1))
      Fvmin_local(ie)    = MINVAL(elem(ie)%derived%FM(:,:,2,:,pnm1))
      
      tmin_local(ie)    = MINVAL(elem(ie)%state%T(:,:,:,n0))
      
      if (rsplit>0) &
           dpmin_local(ie)    = MINVAL(elem(ie)%state%dp3d(:,:,:,n0))
      
      Ftmin_local(ie)    = MINVAL(elem(ie)%derived%FT(:,:,:,pnm1))
      Fqmin_local(ie) = MINVAL(elem(ie)%derived%FQ(:,:,:,1,pnm1))
      Omegamin_local(ie) = MINVAL(elem(ie)%derived%Omega(:,:,:))
      
      
      psmin_local(ie) = MINVAL(tmp(:,:,ie))
      !======================================================
      
      usum_local(ie)    = SUM(elem(ie)%state%v(:,:,1,:,n0))
      vsum_local(ie)    = SUM(elem(ie)%state%v(:,:,2,:,n0))
      Fusum_local(ie)    = SUM(elem(ie)%derived%FM(:,:,1,:,pnm1))
      Fvsum_local(ie)    = SUM(elem(ie)%derived%FM(:,:,2,:,pnm1))
      
      tsum_local(ie)    = SUM(elem(ie)%state%t(:,:,:,n0))
      if (rsplit>0) then
        dpsum_local(ie)    = SUM(elem(ie)%state%dp3d(:,:,:,n0))
      else
        ! Make sure to initialize this to prevent possible
        ! floating point exceptions.
        dpsum_local(ie)    = 0.0D0
      end if
      
      Ftsum_local(ie)    = SUM(elem(ie)%derived%FT(:,:,:,pnm1))
      FQsum_local(ie) = SUM(elem(ie)%derived%FQ(:,:,:,1,pnm1))
      Omegasum_local(ie) = SUM(elem(ie)%derived%Omega(:,:,:))
      
      pssum_local(ie) = SUM(tmp(:,:,ie))
      !======================================================
      
      global_shared_buf(ie,1) = usum_local(ie)
      global_shared_buf(ie,2) = vsum_local(ie)
      global_shared_buf(ie,3) = tsum_local(ie)
      global_shared_buf(ie,4) = pssum_local(ie)
      global_shared_buf(ie,5) = FUsum_local(ie)
      global_shared_buf(ie,6) = FVsum_local(ie)
      global_shared_buf(ie,7) = FTsum_local(ie)
      global_shared_buf(ie,8) = FQsum_local(ie)
      global_shared_buf(ie,9) = Omegasum_local(ie)
      global_shared_buf(ie,10) = dpsum_local(ie)
    end do
    
    !JMD This is a Thread Safe Reduction 
    umin_p = ParallelMin(umin_local,hybrid)
    umax_p = ParallelMax(umax_local,hybrid)
    
    vmin_p = ParallelMin(vmin_local,hybrid)
    vmax_p = ParallelMax(vmax_local,hybrid)
    
    tmin_p = ParallelMin(tmin_local,hybrid)
    tmax_p = ParallelMax(tmax_local,hybrid)
    
    if (rsplit>0)then
      dpmin_p = ParallelMin(dpmin_local,hybrid)
      dpmax_p = ParallelMax(dpmax_local,hybrid)
    endif
    
    psmin_p = ParallelMin(psmin_local,hybrid)
    psmax_p = ParallelMax(psmax_local,hybrid)
    
    FUmin_p = ParallelMin(FUmin_local,hybrid)
    FUmax_p = ParallelMax(FUmax_local,hybrid)
    
    FVmin_p = ParallelMin(FVmin_local,hybrid)
    FVmax_p = ParallelMax(FVmax_local,hybrid)
    
    FTmin_p = ParallelMin(FTmin_local,hybrid)
    FTmax_p = ParallelMax(FTmax_local,hybrid)
    
    FQmin_p = ParallelMin(FQmin_local,hybrid)
    FQmax_p = ParallelMax(FQmax_local,hybrid)
    
    Omegamin = ParallelMin(Omegamin_local,hybrid)
    Omegamax = ParallelMax(Omegamax_local,hybrid)
    
    call wrap_repro_sum(nvars=10, comm=hybrid%par%comm)
    usum_p = global_shared_sum(1)
    vsum_p = global_shared_sum(2)
    tsum_p = global_shared_sum(3)
    pssum_p = global_shared_sum(4)
    FUsum_p = global_shared_sum(5)
    FVsum_p = global_shared_sum(6)
    FTsum_p = global_shared_sum(7)
    FQsum_p = global_shared_sum(8)
    Omegasum = global_shared_sum(9)
    dpsum_p = global_shared_sum(10)
    
    scale=1/g                                  ! assume code is using Pa
    if (hvcoord%ps0 <  2000 ) scale=100*scale  ! code is using mb
    ! after scaling, Energy is in J/m**2,  Mass kg/m**2
    ! so rate of changes are W/m**2
    
    
    !   mass = integral( ps-p(top) )
    do ie=nets,nete
      tmp(:,:,ie)=elem(ie)%state%ps(:,:,n0) 
    enddo
    Mass2 = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    do ie=nets,nete
      tmp(:,:,ie)=elem(ie)%state%pswet(:,:) 
    enddo
    Mass2_moist = global_integral(elem, tmp(:,:,nets:nete),hybrid,npts,nets,nete)
    
    
    
    !       tmp(:,:,ie)=elem(ie)%state%pswet(:,:) 
    
    !
    !   ptop =  hvcoord%hyai(1)*hvcoord%ps0)  + hvcoord%hybi(1)*ps(i,j)
    !   but we assume hybi(1) is zero at top of atmosphere (pure pressure coordinates)
    !    Mass = (Mass2-(hvcoord%hyai(1)*hvcoord%ps0) )*scale  ! this correction is a constant,
    !                                                         ! ~20 kg/m^2 (effects 4th digit of Mass)
    !   BUT: CAM EUL defines mass as integral( ps ), so to be consistent, ignore ptop contribution; 
    Mass       = Mass2*scale
    Mass_moist = Mass2_moist*scale
    
    !
    ! fvm diagnostics
    !
    if (ntrac>0) then
      do q=1,ntrac
        do ie=nets,nete
          tmp1(ie) = MINVAL(fvm(ie)%c(1:nc,1:nc,:,q,n0_fvm)) 
        enddo
        cmin(q) = ParallelMin(tmp1,hybrid)
        do ie=nets,nete
          tmp1(ie) = MAXVAL(fvm(ie)%c(1:nc,1:nc,:,q,n0_fvm))
        enddo
        cmax(q) = ParallelMax(tmp1,hybrid)
        !
        ! compute total tracer mass
        !
        global_shared_buf(:,1) = 0.0D0
        do k=1,nlev
          do ie=nets,nete
            global_shared_buf(ie,1) = global_shared_buf(ie,1)+&
                 SUM(fvm(ie)%c(1:nc,1:nc,k,q,n0_fvm)*&
                 fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)*&
                 fvm(ie)%area_sphere(1:nc,1:nc))
          end do
        enddo
        call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
        csum(q) = global_shared_sum(1)
      enddo
      if (tl%nstep==0) then
        relative_mass_change(:) = 0.0D0
        fvm(nets)%mass(1:ntrac)       = csum(1:ntrac)
        fvm(nets)%mass(ntrac+3)       = mass !se mass
      else
        !
        ! SE mass change
        !
        relative_mass_change(ntrac+3) = (mass-fvm(nets)%mass(ntrac+3))/fvm(nets)%mass(ntrac+3)
        !
        ! fvm tracer mass change
        !
        do q=1,ntrac
          if (ABS(fvm(nets)%mass(q))<1.0D-12) then
            relative_mass_change(q) = csum(q) - fvm(nets)%mass(q)
          else
            relative_mass_change(q) = (csum(q) - fvm(nets)%mass(q))/fvm(nets)%mass(q)
          end if
        end do
      end if
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
        global_shared_buf(ie,1) = 0.0D0
        do j=1,nc
          do i=1,nc
            global_shared_buf(ie,1) = global_shared_buf(ie,1)+fvm(ie)%psc(i,j)*fvm(ie)%area_sphere(i,j)
          end do
        end do
      enddo
      call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
      psc_mass = global_shared_sum(1)
      
      if (tl%nstep==0) then
        fvm(nets)%mass(ntrac+1) = psc_mass
      else
        relative_mass_change(ntrac+1) = (psc_mass - fvm(nets)%mass(ntrac+1))/fvm(nets)%mass(ntrac+1)
      end if
      
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
               SUM(fvm(ie)%dp_fvm(1:nc,1:nc,k,n0_fvm)*fvm(ie)%area_sphere(1:nc,1:nc))!+&
          !                  hvcoord%hyai(1)*hvcoord%ps0
        end do
      enddo
      call wrap_repro_sum(nvars=1, comm=hybrid%par%comm)
      dp_fvm_mass = global_shared_sum(1)
      if (tl%nstep==0) then
        fvm(nets)%mass(ntrac+2) = dp_fvm_mass
      else
        relative_mass_change(ntrac+2) = (dp_fvm_mass - fvm(nets)%mass(ntrac+2))/fvm(nets)%mass(ntrac+2)
      end if
    end if
    
    
    if(hybrid%masterthread) then
      write(iulog,100) "u     = ",umin_p,umax_p,usum_p
      write(iulog,100) "v     = ",vmin_p,vmax_p,vsum_p
      write(iulog,100) "omega = ",omegamin,omegamax,omegasum
      write(iulog,100) "t     = ",tmin_p,tmax_p,tsum_p
      if (rsplit>0) &
           write(iulog,100) "dp    = ",dpmin_p,dpmax_p,dpsum_p
      
      if (tstep_type>0) then  !no longer support tracer advection with tstep_type = 0
        do q=1,qsize
          write(iulog,'(A3,I4,4(E23.15))') "qv ",q,qvmin_p(q), qvmax_p(q), qvsum_p(q),relative_mass_change_se(q)
        enddo
      endif
      if (ldry_mass_vertical_coordinates) then
        write(iulog,100) "psdry= ",psmin_p,psmax_p,pssum_p
        write(iulog,'(a,E23.15,a,E23.15,a)') "  dry M = ",Mass,&
             ' kg/m^2',Mass2,      ' mb     '
      else
        write(iulog,100) "pswet= ",psmin_p,psmax_p,pssum_p
        write(iulog,'(a,E23.15,a,E23.15,a)') "  dry M = ",Mass-SUM(qvsum_p(1:qsize_condensate_loading))*scale,&
             ' kg/m^2',Mass2-SUM(qvsum_p(1:qsize_condensate_loading)),      ' mb     '
      end if
      write(iulog,'(a,E23.15,a,E23.15,a)') "  wet M = ",Mass_moist,' kg/m^2',Mass2_moist,' mb     '
      
      write(iulog,100) "fu = ",fumin_p,fumax_p,fusum_p
      write(iulog,100) "fv = ",fvmin_p,fvmax_p,fvsum_p
      write(iulog,100) "ft = ",ftmin_p,ftmax_p,ftsum_p
      write(iulog,100) "fq = ",fqmin_p, fqmax_p, fqsum_p
      
      !
      ! fvm diagnostics
      !
      if (ntrac>0) then
        write(iulog,'(A36)') "-----------------------------------"
        write(iulog,'(A36)') "fvm diagnostics                    "
        write(iulog,'(A36)') "-----------------------------------"
        do q=1,ntrac
          write(iulog,'(A29,I1,4(E23.15))')&
               "#c, min, max, ave, change = ",q,cmin(q), cmax(q), csum(q)/(4.0D0*DD_PI),relative_mass_change(q)
        enddo
        write(iulog,'(A37,3(E23.15))')&
             "   min(dp_), max(dp_), M loss(dp_) =  ",dp_fvm_min, dp_fvm_max, relative_mass_change(ntrac+2)
        write(iulog,'(A37,3(E23.15))')&
             "   min(psC), max(psC), M loss(psC) =  ",psc_min, psc_max,  relative_mass_change(ntrac+1)
        write(iulog,'(A37,E23.15)') "relative M change SE = ",relative_mass_change(ntrac+3)
        write(iulog,'(a,E23.15,a,E23.15,a)') "  M fvm = ",(psC_mass)*scale/(4.0D0*DD_PI),' kg/m^2',&
             (psC_mass)/(4.0D0*DD_PI),' mb     '
        write(iulog,'(A36)') "                                   "
        
      end if
    end if
100 format (A10,4(E23.15))
    
    if ( test_case(1:10) == "baroclinic" ) then
      ! zeta does not need to be made continious, but we  
      allocate(tmp3d(np,np,nlev,nets:nete))
      call compute_zeta_C0(tmp3d,elem,hybrid,nets,nete,n0)
      tmp=tmp3d(:,:,nlev,:)**2
      relvort = SQRT(global_integral(elem, tmp,hybrid,npts,nets,nete))
      deallocate(tmp3d)
    endif
    
    if(hybrid%par%masterproc .and. hybrid%ithr==0) then 
      if ( test_case(1:10) == "baroclinic" ) then
        write(iulog,101) "2-norm relative vorticity = ",relvort
101     format (A30,E24.15)
        
      else
#ifndef CAM
#ifndef _BGL
        open(unit=10,file=massfname,form="formatted",position="append")
        write(10,*)time/secpday,(Mass2-hvcoord%ps0)/hvcoord%ps0
        close(10)
#endif
#endif
      end if
    endif
  end subroutine prim_printstate
  
  
  subroutine prim_printstate_par(elem, tl,hybrid,hvcoord,nets,nete, par)
    type (element_t), intent(in) :: elem(:)
    type (TimeLevel_t), target, intent(in) :: tl
    type (hybrid_t),intent(in)     :: hybrid
    type (hvcoord_t), intent(in)   :: hvcoord
    integer,intent(in)             :: nets,nete
    character(len=*), parameter    :: fstub = "state_norms"
    integer	                   :: simday
    type(parallel_t)               :: par
    
    real (kind=real_kind)  :: v(np,np,2,nlev,nets:nete)
    real (kind=real_kind)  :: t(np,np,nlev,nets:nete)
    real (kind=real_kind)  :: ps(np,np,nets:nete)
    real (kind=real_kind)  :: vp(np,np,2,nlev,nets:nete)
    real (kind=real_kind)  :: tp(np,np,nlev,nets:nete)
    real (kind=real_kind)  :: psp(np,np,nets:nete)
    real (kind=real_kind) :: l1,l2,linf
    integer               :: n0,i,j,k,ie,npts

    npts=SIZE(elem(1)%state%lnps(:,:,n0),1)
    n0=tl%n0
    do ie=nets,nete
       v(:,:,:,:,ie)=elem(ie)%state%v(:,:,:,:,n0) 
       T(:,:,:,ie)=elem(ie)%state%T(:,:,:,n0) 
       ps(:,:,ie)=elem(ie)%state%ps(:,:,n0) 
    enddo
    simday = 0

#ifdef _REFSOLNW
! parallel write file with state vector in unformatted blocks for later calculation of norms
    call ref_state_write(v(:,:,:,:,nets:nete),T(:,:,:,nets:nete),ps(:,:,nets:nete), & 
	fstub,simday,nets,nete,par)
    do ie=nets,nete
       vp(:,:,:,:,ie)=v(:,:,:,:,ie)
       Tp(:,:,:,ie)=T(:,:,:,ie)
      psp(:,:,ie)=ps(:,:,ie)
    end do
#endif

#ifdef _REFSOLNR
! parallel read file with state vector in unformatted blocks as written above
    !$OMP BARRIER
!  Parallel version of ref_state, comment out if writing above
    call ref_state_read(vp(:,:,:,:,nets:nete),Tp(:,:,:,nets:nete),psp(:,:,nets:nete), & 
	fstub,simday,nets,nete,par)
    !$OMP BARRIER

    npts=np

    l1   = l1_snorm(elem,ps(:,:,nets:nete),  psp(:,:,nets:nete),hybrid,npts,nets,nete)
    l2   = l2_snorm(elem,ps(:,:,nets:nete),  psp(:,:,nets:nete),hybrid,npts,nets,nete)
    linf = linf_snorm(ps(:,:,nets:nete),psp(:,:,nets:nete),hybrid,npts,nets,nete)

    if (hybrid%par%masterproc .and. (hybrid%ithr==0)) then
       print *,simday, "L1=",l1
       print *,simday, "L2=",l2
       print *,simday, "Linf=",linf
    end if
    !$OMP BARRIER
#endif

  end subroutine prim_printstate_par
!=======================================================================================================! 
end module prim_state_mod
