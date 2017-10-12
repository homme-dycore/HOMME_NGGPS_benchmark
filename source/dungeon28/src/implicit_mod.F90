module implicit_mod
  
  use ,intrinsic :: iso_c_binding 
  use perf_mod, only : t_startf, t_stopf ! EXTERNAL
  use perf_utils, only : t_detail_low, t_detail_medium, t_detail_high, t_detail_max ! EXTERNAL
  use derivative_mod, only : derivative_t
  use derivative_mod, only : gradient_sphere_routine, divergence_sphere_routine
  use derivative_mod, only : vorticity_sphere_routine, vorticity_sphere_diag_function
  
#ifdef SWIM_USE_OPENACC
  use kinds          , only: real_kind
  
  implicit none
  private
  real(kind=real_kind), allocatable :: vco(:,:,:,:,:),vvtemp(:,:,:,:),vtemp(:,:,:,:),v1(:,:,:,:),v2(:,:,:,:)
  real(kind=real_kind), allocatable :: metdet(:,:,:),rmetdet(:,:,:),Dvv(:,:),Dinv(:,:,:,:,:),D(:,:,:,:,:)
  real(kind=real_kind), allocatable :: vn0(:,:,:,:,:),ps0(:,:,:),pn0(:,:,:,:),spheremp(:,:,:),fcor(:,:,:),gv(:,:,:,:,:)
  real(kind=real_kind), allocatable :: E(:,:,:,:),E_n(:,:,:,:),pv(:,:,:,:,:),pv_n(:,:,:,:,:),px(:,:,:,:),up(:,:,:,:,:)
  real(kind=real_kind), allocatable :: vtens(:,:,:,:,:),vtens_n(:,:,:,:,:),ptens(:,:,:,:),ptens_n(:,:,:,:)
  real(kind=real_kind), allocatable :: grade(:,:,:,:,:),grade_n(:,:,:,:,:),zeta(:,:,:,:),zeta_n(:,:,:,:),div(:,:,:,:),div_n(:,:,:,:)
  
  public :: acc_imp_int, acc_imp_fin, advance_imp_nonstag

#endif
contains
  
  subroutine advance_imp_nonstag(elem, edge1, edge2, edge3, red, deriv, cg,   &
       hybrid, blkjac, lambdasq, dt, pmean, tl, nets, nete, xstate)
    
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar
    use physical_constants, only : g
    use element_mod, only : element_t
    use parallel_mod, only : parallel_t, abortmp
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use time_mod, only : timelevel_t, smooth
    use control_mod, only :  topology, test_case
    use shallow_water_mod, only : tc1_velocity
    use bndry_mod, only : bndry_exchangev
    use cg_mod, only : cg_t
    use precon_mod, only : pcg_presolver_nonstag
    use solver_mod, only : blkjac_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use derived_type_mod ,only : derived_type, initialize
    use precon_type_mod ,only : precon_type, init_precon
    use ,intrinsic :: iso_c_binding 
    use global_norms_mod

    implicit none

    type (element_t) ,intent(inout)   :: elem(:)
    type (EdgeBuffer_t)  , intent(in) :: edge1
    type (EdgeBuffer_t)  , intent(in) :: edge2
    type (EdgeBuffer_t)  , intent(in) :: edge3
    type (reductionbuffer_ordered_1d_t)  ,intent(in) :: red
    type (derivative_t)  , intent(in) :: deriv
    type (cg_t)                       :: cg
    type (hybrid_t)      , intent(in) :: hybrid
    type (blkjac_t)      , intent(inout) :: blkjac(:)
    real (kind=real_kind), intent(inout) :: lambdasq(:)
    real (kind=real_kind), intent(in) :: dt
    real (kind=real_kind), intent(in) :: pmean
    type (timeLevel_t)   , intent(in) :: tl
    integer              , intent(in) :: nets
    integer              , intent(in) :: nete

    ! =================
    ! Local
    ! =================

    ! Thread private working set ...

    type (parallel_t)          :: par
    real (c_double) ,allocatable ,dimension(:) :: xstate

    type (element_t)  :: pc_elem(size(elem))
    type (element_t)  :: jac_elem(size(elem))

!    type (element_t)  :: pc_elem(nets-nete+1)
!    type (element_t)  :: jac_elem(nets-nete+1)

    integer    :: i,j,k,n,ie,shiftv,shiftp
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lenx, lx
    integer    :: comm
    integer    :: state_size

! state_object is a derived data type passed from doloca thru cpp to calc_f
! in the form of a pointer
    type(derived_type) ,target         :: state_object
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object

    !type(c_funptr)                     :: res_calc
    !type(precon_type) ,target          :: pre_object
    !type(precon_type) ,pointer         :: pptr=>NULL()

    type(derived_type) ,target          :: pre_object
    type(derived_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                         :: c_ptr_to_pre

    type(derived_type) ,target          :: ajac_object
    type(derived_type) ,pointer         :: jptr=>NULL()
    type(c_ptr)                         :: c_ptr_to_jac
    
    integer(c_int) :: ierr = 0
    
    !    type(c_funptr)                     :: precon
    
    interface 
       subroutine noxsolve(vectorSize, vector, v_container, p_container, j_container, ierr) &
            bind(C,name='noxsolve')
         use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
         integer(c_int)                :: vectorSize
         real(c_double)  ,dimension(*) :: vector
         type(c_ptr)                   :: v_container 
         type(c_ptr)                   :: p_container  ! precon ptr
         type(c_ptr)                   :: j_container  ! jac ptr
         integer(c_int)                :: ierr         ! error flag
       end subroutine noxsolve
    end interface
    
    nm1   = tl%nm1
    n0    = tl%n0
    np1   = tl%np1
    nstep = tl%nstep

    comm = par%comm
! n+1 guess is result at n - need this to match up with time_mod bootstrap style
    call t_startf('implicit header', t_detail_high)
    call t_startf('implicit_pre_noxsolve', t_detail_high)

    shiftv = np*np*nlev*(nete-nets+1)
    shiftp = 2*np*np*nlev*(nete-nets+1)
    lx = 0
    do ie=nets,nete
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                xstate(lx) = elem(ie)%state%v(i,j,1,k,n0)
                xstate(lx+shiftv) = elem(ie)%state%v(i,j,2,k,n0)
                xstate(lx+shiftp) = elem(ie)%state%p(i,j,k,n0)
#ifdef SWIM_USE_OPENACC
     vn0(i,j,1,k,ie) = elem(ie)%state%v(i,j,1,k,n0)
     vn0(i,j,2,k,ie) = elem(ie)%state%v(i,j,2,k,n0)
     pn0(i,j,k,ie) = elem(ie)%state%p(i,j,k,n0)
     ps0(i,j,ie) = elem(ie)%state%ps(i,j)
#endif
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie
 
#ifdef SWIM_USE_OPENACC
!$omp barrier
!$omp master
 
!$acc update device(vn0,pn0,ps0)
 
!$omp end master
!$omp barrier
#endif
    
    call t_stopf('implicit_pre_noxsolve', t_detail_high)
    
    call t_startf('implicit_init', t_detail_high)

    pc_elem=elem
    jac_elem=elem
     
    call initialize(state_object, lenx, elem, pmean, edge1,edge2,edge3, &
         hybrid, deriv, dt, tl, nets, nete)
    
    call initialize(pre_object, lenx, pc_elem, pmean, edge1,edge2,edge3, &
         hybrid, deriv, dt, tl, nets, nete)
    
    call initialize(ajac_object, lenx, jac_elem, pmean, edge1,edge2,edge3, &
         hybrid, deriv, dt, tl, nets, nete)

    
   ! call init_precon(pre_object, lenx, elem, blkjac, edge1, edge2, edge3, &
   !     red, deriv, cg, lambdasq, dt, pmean, tl, nets, nete)

    fptr => state_object
    c_ptr_to_object =  c_loc(fptr)

    pptr => pre_object
    c_ptr_to_pre =  c_loc(pptr)

    jptr => ajac_object
    c_ptr_to_jac =  c_loc(jptr)
    call t_stopf('implicit_init', t_detail_high)
 

! ForTrilinos interface to use nox and loca, and returns xstate(n+1)

    call t_startf('noxsolve', t_detail_high)
    call noxsolve(size(xstate), xstate, c_ptr_to_object, c_ptr_to_pre, c_ptr_to_jac, ierr)   
    if (ierr /= 0) call abortmp('Error in noxsolve: Newton failed to converge')
    call t_stopf('noxsolve', t_detail_high)
    call t_startf('implicit_post_noxsolve', t_detail_high)

!Moved conditionals out of do loops
      lx = 0
      do ie=nets,nete
        do k=1,nlev
          do j=1,np
            do i=1,np
              lx = lx+1
              elem(ie)%state%v(i,j,1,k,np1) = xstate(lx)
              elem(ie)%state%v(i,j,2,k,np1) = xstate(lx+shiftv)
              elem(ie)%state%p(i,j,k,np1)   = xstate(lx+shiftp)
            end do  !np
          end do  !np
        end do  !nlev
      end do !ie
      
     call t_stopf('implicit_post_noxsolve', t_detail_high)

       ! ==========================================
       ! If SW test case 1, velocity is prescribed.
       ! reset v back to initial values
       ! ==========================================

       ! TODO update with vortex and swirl possibly using set_prescribed_velocity
       do ie=nets,nete
       if (topology == "cube" .and. test_case=="swtc1") then
          do k=1,nlev
             elem(ie)%state%v(:,:,:,k,np1)=tc1_velocity(elem(ie)%spherep,elem(ie)%Dinv)
             elem(ie)%state%v(:,:,:,k,n0 )=elem(ie)%state%v(:,:,:,k,np1)
             elem(ie)%state%v(:,:,:,k,nm1)=elem(ie)%state%v(:,:,:,k,np1)
          end do
       end if
       end do !ie

    call t_stopf('implicit header', t_detail_high)

    !$OMP BARRIER

  end subroutine advance_imp_nonstag
  
#ifdef SWIM_USE_OPENACC
  subroutine acc_imp_int(elem, deriv, nets,nete)
    use element_mod, only : element_t
    use dimensions_mod, only : np,nlev,nelem
 
 type (element_t), intent(in)    :: elem(:)
 type (derivative_t), intent(in) :: deriv 
 integer, intent(in)             :: nets  
 integer, intent(in)             :: nete  
 integer  :: i,j,ie
 
 allocate(Dvv(np,np))
 allocate(vco(np,np,2,nlev,nelem))
 allocate(vtemp(np,np,nlev,nelem))
 allocate(vvtemp(np,np,nlev,nelem))
 allocate(v1(np,np,nlev,nelem))
 allocate(v2(np,np,nlev,nelem))
 allocate(rmetdet(np,np,nelem))
 allocate(metdet(np,np,nelem))
 allocate(Dinv(np,np,2,2,nelem))
 allocate(D(np,np,2,2,nelem))
 allocate(vn0(np,np,2,nlev,nelem))
 allocate(ps0(np,np,nelem))
 allocate(pn0(np,np,nlev,nelem))
 allocate(spheremp(np,np,nelem))
 allocate(fcor(np,np,nelem))
 allocate(gv(np,np,2,nlev,nelem))
 allocate(E(np,np,nlev,nelem))
 allocate(E_n(np,np,nlev,nelem))
 allocate(pv(np,np,2,nlev,nelem))
 allocate(pv_n(np,np,2,nlev,nelem))
 allocate(up(np,np,2,nlev,nelem))
 allocate(px(np,np,nlev,nelem))
 allocate(vtens(np,np,2,nlev,nelem))
 allocate(vtens_n(np,np,2,nlev,nelem))
 allocate(ptens(np,np,nlev,nelem))
 allocate(ptens_n(np,np,nlev,nelem))
 allocate(grade(np,np,2,nlev,nelem))
 allocate(grade_n(np,np,2,nlev,nelem))
 allocate(zeta(np,np,nlev,nelem))
 allocate(zeta_n(np,np,nlev,nelem))
 allocate(div(np,np,nlev,nelem))
 allocate(div_n(np,np,nlev,nelem))

 
 do j=1,np
  do i=1,np
   Dvv(i,j) = deriv%Dvv(i,j)
   do ie=nets,nete
    rmetdet(i,j,ie) = elem(ie)%rmetdet(i,j)
    metdet(i,j,ie) = elem(ie)%metdet(i,j)
    spheremp(i,j,ie) = elem(ie)%spheremp(i,j)
    fcor(i,j,ie) = elem(ie)%fcor(i,j)
      
    Dinv(i,j,1,1,ie) = elem(ie)%Dinv(i,j,1,1)
    Dinv(i,j,2,1,ie) = elem(ie)%Dinv(i,j,2,1)
    Dinv(i,j,1,2,ie) = elem(ie)%Dinv(i,j,1,2)
    Dinv(i,j,2,2,ie) = elem(ie)%Dinv(i,j,2,2)
   
    D(i,j,1,1,ie) = elem(ie)%D(i,j,1,1)
    D(i,j,2,1,ie) = elem(ie)%D(i,j,2,1)
    D(i,j,1,2,ie) = elem(ie)%D(i,j,1,2)
    D(i,j,2,2,ie) = elem(ie)%D(i,j,2,2)
      
   end do
  end do
 end do
 
!$omp barrier
!$omp master
 
!$acc enter data create(vvtemp,vtemp,vco,v1,v2,gv,E,E_n,pv,pv_n,vn0,pn0,ps0,up,px)
!$acc enter data create(vtens,vtens_n,ptens,ptens_n,grade,grade_n,zeta,zeta_n,div,div_n)
!$acc enter data copyin(Dvv,rmetdet,metdet,Dinv,D,spheremp,fcor)
 
!$omp end master
!$omp barrier 

  end subroutine acc_imp_int
  
  subroutine acc_imp_fin()

 deallocate(Dvv,vco,vtemp,vvtemp,v1,v2,rmetdet,metdet,Dinv,D,vn0,ps0,pn0,spheremp,fcor,gv,E,E_n,pv,pv_n,vn0,pn0,ps0,up,px)
 deallocate(vtens,vtens_n,ptens,ptens_n,grade,grade_n,zeta,zeta_n,div,div_n)
  end subroutine acc_imp_fin
   
  subroutine residual(xstate, fx, nelemd, c_ptr_to_object) bind(C,name='calc_f')

    use ,intrinsic :: iso_c_binding 
 use physical_constants, only : rrearth 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only: EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case, tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type

    implicit none

    real (c_double) ,intent(in)        :: xstate(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne
 
    real(kind=real_kind) ::  dsdx00, dsdy00, dudy00, dvdx00, dudx00,dvdy00
 
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam

    integer    :: i,j,k,l,n,ie,shiftv,shiftp
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lx

    call t_startf('FE implicit_acc', t_detail_high)
    call t_startf('FE_implicit_acc_init', t_detail_high)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

!Moved conditionals out of do loops
       shiftv = np*np*nlev*(ne-ns+1)
       shiftp = 2*np*np*nlev*(ne-ns+1)

       call t_stopf('FE_implicit_acc_init', t_detail_high)
       call t_startf('FE_implicit_acc_KE_resid_calc', t_detail_high)

  
        do ie=ns,ne
         do k=1,nlev
   do j=1,np
    do i=1,np
     ! set du
     lx = i+(j-1)*np+(k-1)*np*np+(ie-ns)*np*np*nlev
     up(i,j,1,k,ie)      = xstate(lx)   !
     up(i,j,2,k,ie)      = xstate(lx+shiftv)   ! 
     px(i,j,k,ie)        = xstate(lx+shiftp)
    end do
   end do
    end do
    end do
    
       !$omp barrier
       !$omp master

    !$acc update device(up,px)
    !$acc parallel loop gang vector collapse(4) present(E,E_n,pv,pv_n,up,px,ps0,vn0,pn0) 
       do ie=ns,ne
         do k=1,nlev 
   do j=1,np
    do i=1,np
     E(i,j,k,ie)   = 0.5D0*(up(i,j,1,k,ie)**2 + up(i,j,2,k,ie)**2)  + px(i,j,k,ie) + ps0(i,j,ie)
     !if (tstep_type==13) then  !compute extra terms needed for CN
      E_n(i,j,k,ie) = 0.5D0*(vn0(i,j,1,k,ie)**2 + vn0(i,j,2,k,ie)**2)  + pn0(i,j,k,ie) + ps0(i,j,ie)
      pv_n(i,j,1,k,ie) = vn0(i,j,1,k,ie)*(pmean + pn0(i,j,k,ie) )
      pv_n(i,j,2,k,ie) = vn0(i,j,2,k,ie)*(pmean + pn0(i,j,k,ie) )
     !endif

     pv(i,j,1,k,ie) = up(i,j,1,k,ie)*( pmean + px(i,j,k,ie) )
     pv(i,j,2,k,ie) = up(i,j,2,k,ie)*( pmean + px(i,j,k,ie) )
    end do
   end do
   end do
       end do
    
    !$acc parallel loop gang vector collapse(4) present(Dvv,E,v1,v2) private(dsdx00,dsdy00)
       do ie=ns,ne
        do k=1,nlev 
   do j=1,np
      do l=1,np
      dsdx00=0.0d0
      dsdy00=0.0d0
      do i=1,np
      dsdx00 = dsdx00 + Dvv(i,l)*E(i,j,k,ie)
      dsdy00 = dsdy00 + Dvv(i,l)*E(j,i,k,ie)
      end do
      v1(l,j,k,ie) = dsdx00*rrearth
      v2(j,l,k,ie) = dsdy00*rrearth
      end do
   end do
    end do
      end do
      
   !$acc parallel loop gang vector collapse(4) present(grade,Dinv,v1,v2) 
      do ie=ns,ne
       do k=1,nlev 
   do j=1,np
      do i=1,np
      grade(i,j,1,k,ie)=Dinv(i,j,1,1,ie)*v1(i,j,k,ie) + Dinv(i,j,2,1,ie)*v2(i,j,k,ie)
      grade(i,j,2,k,ie)=Dinv(i,j,1,2,ie)*v1(i,j,k,ie) + Dinv(i,j,2,2,ie)*v2(i,j,k,ie)
      enddo
   enddo
         end do
      end do
   
    !$acc parallel loop gang vector collapse(4) present(vco,D,up) 
      do ie=ns,ne
       do k=1,nlev 
   do j=1,np
      do i=1,np
      vco(i,j,1,k,ie)=(D(i,j,1,1,ie)*up(i,j,1,k,ie) + D(i,j,2,1,ie)*up(i,j,2,k,ie))
      vco(i,j,2,k,ie)=(D(i,j,1,2,ie)*up(i,j,1,k,ie) + D(i,j,2,2,ie)*up(i,j,2,k,ie))
      enddo
   enddo
    end do
  end do
 
 
  !$acc parallel loop gang vector collapse(4) present(Dvv,vco,zeta,vtemp) private(dudy00,dvdx00)
  do ie=ns,ne
   do k=1,nlev   
   do j=1,np
      do l=1,np

      dudy00=0.0d0
      dvdx00=0.0d0

      do i=1,np
      dvdx00 = dvdx00 + Dvv(i,l)*vco(i,j,2,k,ie)
      dudy00 = dudy00 + Dvv(i,l)*vco(j,i,1,k,ie)
      enddo
 
      zeta(l,j,k,ie) = dvdx00
      vtemp(j,l,k,ie) = dudy00
      enddo
   enddo
   end do
    end do


    !$acc parallel loop gang vector collapse(4) present(zeta,vtemp,rmetdet)
    do ie=ns,ne
     do k=1,nlev   
   do j=1,np
      do i=1,np
      zeta(i,j,k,ie)=(zeta(i,j,k,ie)-vtemp(i,j,k,ie))*(rmetdet(i,j,ie)*rrearth)
      end do
   end do
   
   end do
    end do


    !$acc parallel loop gang vector collapse(4) present(gv,pv,Dinv,metdet) 
    do ie=ns,ne
     do k=1,nlev
   
   do j=1,np
      do i=1,np
      gv(i,j,1,k,ie)=metdet(i,j,ie)*(Dinv(i,j,1,1,ie)*pv(i,j,1,k,ie) + Dinv(i,j,1,2,ie)*pv(i,j,2,k,ie))
      gv(i,j,2,k,ie)=metdet(i,j,ie)*(Dinv(i,j,2,1,ie)*pv(i,j,1,k,ie) + Dinv(i,j,2,2,ie)*pv(i,j,2,k,ie))
      enddo
   enddo

   end do
    end do


    !$acc parallel loop gang vector collapse(4) present(Dvv,gv,div,vvtemp) private(dudx00,dvdy00)
    do ie=ns,ne
     do k=1,nlev
   do j=1,np
      do l=1,np
      dudx00=0.0d0
      dvdy00=0.0d0
      do i=1,np
      dudx00 = dudx00 + Dvv(i,l)*gv(i,j,1,k,ie)
      dvdy00 = dvdy00 + Dvv(i,l)*gv(j,i,2,k,ie)
      end do
      div(l,j,k,ie) = dudx00
      vvtemp(j,l,k,ie) = dvdy00
      end do
   end do
   end do
    end do


    !$acc parallel loop gang vector collapse(4) present(div,vvtemp,rmetdet) 
    do ie=ns,ne
     do k=1,nlev   
   do j=1,np
      do i=1,np
      div(i,j,k,ie)=(div(i,j,k,ie)+vvtemp(i,j,k,ie))*(rmetdet(i,j,ie)*rrearth)
      end do
   end do
   end do
    end do
 
 !$acc parallel loop gang vector collapse(4) present(Dvv,E_n,v1,v2) private(dsdx00,dsdy00)
    do ie=ns,ne
     do k=1,nlev
   do j=1,np
      do l=1,np
      dsdx00=0.0d0
      dsdy00=0.0d0
      do i=1,np
      dsdx00 = dsdx00 + Dvv(i,l)*E_n(i,j,k,ie)
      dsdy00 = dsdy00 + Dvv(i,l)*E_n(j,i,k,ie)
      end do
      v1(l,j,k,ie) = dsdx00*rrearth
      v2(j,l,k,ie) = dsdy00*rrearth
      end do
   end do
   end do
    end do


    !$acc parallel loop gang vector collapse(4) present(grade_n,Dinv,v1,v2) 
    do ie=ns,ne
     do k=1,nlev
   do j=1,np
      do i=1,np
      grade_n(i,j,1,k,ie)=Dinv(i,j,1,1,ie)*v1(i,j,k,ie) + Dinv(i,j,2,1,ie)*v2(i,j,k,ie)
      grade_n(i,j,2,k,ie)=Dinv(i,j,1,2,ie)*v1(i,j,k,ie) + Dinv(i,j,2,2,ie)*v2(i,j,k,ie)
      enddo
   enddo
   end do
    end do


 !$acc parallel loop gang vector collapse(4) present(vco,D,vn0) 
    do ie=ns,ne
     do k=1,nlev
   do j=1,np
      do i=1,np
      vco(i,j,1,k,ie)=(D(i,j,1,1,ie)*vn0(i,j,1,k,ie) + D(i,j,2,1,ie)*vn0(i,j,2,k,ie))
      vco(i,j,2,k,ie)=(D(i,j,1,2,ie)*vn0(i,j,1,k,ie) + D(i,j,2,2,ie)*vn0(i,j,2,k,ie))
      enddo
   enddo
   end do
    end do


 !$acc parallel loop gang vector collapse(4) present(Dvv,vco,zeta_n,vtemp) private(dudy00,dvdx00)
    do ie=ns,ne
     do k=1,nlev   
   do j=1,np
      do l=1,np

      dudy00=0.0d0
      dvdx00=0.0d0

      do i=1,np
      dvdx00 = dvdx00 +Dvv(i,l)*vco(i,j,2,k,ie)
      dudy00 = dudy00 +Dvv(i,l)*vco(j,i,1,k,ie)
      enddo

      zeta_n(l,j,k,ie) = dvdx00
      vtemp(j,l,k,ie) = dudy00
      enddo
   enddo
   end do
    end do


    !$acc parallel loop gang vector collapse(4) present(zeta_n,vtemp,rmetdet)
    do ie=ns,ne
     do k=1,nlev   

   do j=1,np
      do i=1,np
      zeta_n(i,j,k,ie)=(zeta_n(i,j,k,ie)-vtemp(i,j,k,ie))*(rmetdet(i,j,ie)*rrearth)
      end do
   end do
   end do
    end do


    !$acc parallel loop gang vector collapse(4) present(gv,Dinv,pv_n,metdet)
    do ie=ns,ne
     do k=1,nlev   

   do j=1,np
      do i=1,np
      gv(i,j,1,k,ie)=metdet(i,j,ie)*(Dinv(i,j,1,1,ie)*pv_n(i,j,1,k,ie) + Dinv(i,j,1,2,ie)*pv_n(i,j,2,k,ie))
      gv(i,j,2,k,ie)=metdet(i,j,ie)*(Dinv(i,j,2,1,ie)*pv_n(i,j,1,k,ie) + Dinv(i,j,2,2,ie)*pv_n(i,j,2,k,ie))
      enddo
   enddo

   end do
    end do


 !$acc parallel loop gang vector collapse(4) present(Dvv,gv,div_n,vvtemp) private(dudx00,dvdy00)
    do ie=ns,ne
     do k=1,nlev     
   do j=1,np
      do l=1,np
      dudx00=0.0d0
      dvdy00=0.0d0
      do i=1,np
      dudx00 = dudx00 + Dvv(i,l)*gv(i,j,1,k,ie)
      dvdy00 = dvdy00 + Dvv(i,l)*gv(j,i,2,k,ie)
      end do
      div_n(l,j,k,ie) = dudx00
      vvtemp(j,l,k,ie) = dvdy00
      end do
   end do
   end do
    end do


 !$acc parallel loop gang vector collapse(4) present(div_n,vvtemp,rmetdet)
    do ie=ns,ne
     do k=1,nlev   
   do j=1,np
      do i=1,np
      div_n(i,j,k,ie)=(div_n(i,j,k,ie)+vvtemp(i,j,k,ie))*(rmetdet(i,j,ie)*rrearth)
      end do
   end do
   end do
    end do

    
 !$acc parallel loop gang vector collapse(4) present(vtens,vtens_n,ptens,ptens_n,grade,div,zeta,grade_n,div_n,zeta_n,spheremp,fcor,vn0,up,px,pn0)
    do ie=ns,ne
     do k=1,nlev   
   do j=1,np
    do i=1,np
       vtens_n(i,j,1,k,ie)=spheremp(i,j,ie)*(vn0(i,j,2,k,ie)*(fcor(i,j,ie) + zeta_n(i,j,k,ie)) - grade_n(i,j,1,k,ie))
       vtens_n(i,j,2,k,ie)=spheremp(i,j,ie)*(-vn0(i,j,1,k,ie)*(fcor(i,j,ie) + zeta_n(i,j,k,ie)) - grade_n(i,j,2,k,ie))
       ptens_n(i,j,k,ie) =  -spheremp(i,j,ie)*div_n(i,j,k,ie)
       vtens(i,j,1,k,ie)=spheremp(i,j,ie)*(up(i,j,2,k,ie)*(fcor(i,j,ie) + zeta(i,j,k,ie)) - grade(i,j,1,k,ie))
       vtens(i,j,2,k,ie)=spheremp(i,j,ie)*(-up(i,j,1,k,ie)*(fcor(i,j,ie) + zeta(i,j,k,ie)) - grade(i,j,2,k,ie))

       ptens(i,j,k,ie) =  -spheremp(i,j,ie)*div(i,j,k,ie)
       
       ptens(i,j,k,ie) = spheremp(i,j,ie)*(px(i,j,k,ie)-pn0(i,j,k,ie))*dti - 0.5*ptens(i,j,k,ie) - 0.5*ptens_n(i,j,k,ie)

       vtens(i,j,1,k,ie) = spheremp(i,j,ie)*(up(i,j,1,k,ie)-vn0(i,j,1,k,ie))*dti - 0.5*vtens(i,j,1,k,ie) - 0.5*vtens_n(i,j,1,k,ie)

       vtens(i,j,2,k,ie) = spheremp(i,j,ie)*(up(i,j,2,k,ie)-vn0(i,j,2,k,ie))*dti - 0.5*vtens(i,j,2,k,ie) - 0.5*vtens_n(i,j,2,k,ie)
        
    end do
   end do
  end do
    end do
 
    !$acc update host(ptens,vtens)
    !$omp end master
    !$omp barrier

    do ie=ns,ne
  ! ===================================================
  ! Pack cube edges of tendencies, rotate velocities
  ! ===================================================
  kptr=0
  call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,ie)
  kptr=nlev
  call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
 end do

   !pw++
   call t_stopf('FE_implicit_acc_KE_resid_calc', t_detail_high)
   !pw--

   !$OMP BARRIER
   !pw++
   call t_startf('FE_implicit_acc_bndry_ex', t_detail_high)
   !pw--
   !$OMP BARRIER
   call bndry_exchangeV(fptr%hybrid,fptr%edge3)
   !$OMP BARRIER
   !pw++
   call t_stopf('FE_implicit_acc_bndry_ex', t_detail_high)
   !pw--
   !$OMP BARRIER

   !pw++
   call t_startf('FE_implicit_acc_bndry_unpack', t_detail_high)
   !pw--
   

   do ie=ns,ne

    ! ===========================================================
    ! Unpack the edges for vgradp and vtens
    ! ===========================================================
    kptr=0
    call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, ie) 

    kptr=nlev
    call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

   end do !ie

   call t_stopf('FE_implicit_acc_bndry_unpack', t_detail_high)

   ! ===========================================================
   ! Compute velocity and pressure tendencies for all levels
   ! ===========================================================

   call t_startf('FE_implicit_acc_vel_pres', t_detail_high)
 

   !Moved conditionals out of do loops
   lx = 0
   if (topology == "cube" .and. test_case=="swtc1") then
    do ie=ns,ne
     do k=1,nlev
      do j=1,np
       do i=1,np
        lx = lx+1
        fx(lx) = 0.0
        fx(lx+shiftv) = 0.0
        fx(lx+shiftp) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
       end do
      end do
     end do
    end do !ie
   else
    do ie=ns,ne
     do k=1,nlev
      do j=1,np
       do i=1,np
        lx = lx+1
        fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
        fx(lx+shiftv)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
        fx(lx+shiftp) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
       end do
      end do
     end do
    end do !ie
   end if

     call t_stopf('FE_implicit_acc_vel_pres', t_detail_high)
     call t_stopf('FE implicit_acc', t_detail_high)

  end subroutine residual
#else
  subroutine residual(xstate, fx, nelemd, c_ptr_to_object) bind(C,name='calc_f')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only: EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case, tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type

    implicit none

    real (c_double) ,intent(in)        :: xstate(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,ps,px
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens_n
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens_n ! at t=n

    real (kind=real_kind), dimension(np,np,2)  :: grade,grade_n ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: pv,pv_n     ! p*v lat-lon at t=n+1,n
    real (kind=real_kind), dimension(np,np)    :: E,E_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta,zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: div,div_n   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up, um

    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam

    integer    :: i,j,k,n,ie,shiftv,shiftp
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: nstep
    integer    :: lx

    call t_startf('FE implicit', t_detail_high)
    call t_startf('FE_implicit_init', t_detail_high)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0         = fptr%tl%n0
    np1        = fptr%tl%np1
    nm1        = fptr%tl%nm1
    nstep      = fptr%tl%nstep
    dti        = 1.0d0/fptr%dt
    ns         = fptr%nets 
    ne         = fptr%nete 
    pmean      = fptr%pmean
    gam        = 0.5

!Moved conditionals out of do loops
       shiftv = np*np*nlev*(ne-ns+1)
       shiftp = 2*np*np*nlev*(ne-ns+1)

       call t_stopf('FE_implicit_init', t_detail_high)
       call t_startf('FE_implicit_KE_resid_calc', t_detail_high)
        lx = 0
        do ie=ns,ne
             do j=1,np
                do i=1,np
                   fcor(i,j)        = fptr%base(ie)%fcor(i,j)
                   spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
                   ps(i,j)          = fptr%base(ie)%state%ps(i,j)
                end do
             end do
        do k=1,nlev

! ==============================================
! Compute kinetic energy term at each time level
! ==============================================

   do j=1,np
    do i=1,np

     ! set u
     ulatlon(i,j,1)       = fptr%base(ie)%state%v(i,j,1,k,n0)   !
     ulatlon(i,j,2)       = fptr%base(ie)%state%v(i,j,2,k,n0)   !

     ! set du
     lx = lx+1
     up(i,j,1)      = xstate(lx)   !
     up(i,j,2)      = xstate(lx+shiftv)   ! 
     px(i,j)        = xstate(lx+shiftp)
     ! set un-1
     um(i,j,1)      = fptr%base(ie)%state%v(i,j,1,k,nm1)   ! 
     um(i,j,2)      = fptr%base(ie)%state%v(i,j,2,k,nm1)   !  


     E(i,j)   = 0.5D0*(up(i,j,1)**2 + up(i,j,2)**2)  + &
     px(i,j) + ps(i,j)
     if (tstep_type==13) then  !compute extra terms needed for CN
      E_n(i,j) = 0.5D0*(ulatlon(i,j,1)**2 + ulatlon(i,j,2)**2)  + &
      fptr%base(ie)%state%p(i,j,k,n0) + ps(i,j)
      pv_n(i,j,1) = ulatlon(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
      pv_n(i,j,2) = ulatlon(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
     endif

     pv(i,j,1) = up(i,j,1)*( fptr%pmean + px(i,j) )
     pv(i,j,2) = up(i,j,2)*( fptr%pmean + px(i,j) )

    end do
   end do

   !call flush(6)
   !stop
   call gradient_sphere_routine(E,fptr%deriv,fptr%base(ie)%Dinv,grade)  ! scalar -> latlon vector
   call vorticity_sphere_routine(up,fptr%deriv,fptr%base(ie),zeta) ! latlon vector -> scalar 
   call divergence_sphere_routine(pv,fptr%deriv,fptr%base(ie),div) ! latlon vector -> scalar 
   ! residual time level n
                        if (tstep_type==13) then  
    !compute extra terms needed for CN
    call gradient_sphere_routine(E_n,fptr%deriv,fptr%base(ie)%Dinv,grade_n)  ! scalar -> latlon vector
    call vorticity_sphere_routine(ulatlon,fptr%deriv,fptr%base(ie),zeta_n) ! latlon vector -> scalar 
    call divergence_sphere_routine(pv_n,fptr%deriv,fptr%base(ie),div_n) ! latlon vector -> scalar 
    do j=1,np
     do i=1,np
        vtens_n(i,j,1,k,ie)=spheremp(i,j)* &
        (ulatlon(i,j,2)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,1))
        vtens_n(i,j,2,k,ie)=spheremp(i,j)* &
        (-ulatlon(i,j,1)*(fcor(i,j) + zeta_n(i,j)) - grade_n(i,j,2))
        ptens_n(i,j,k,ie) =  -spheremp(i,j)*div_n(i,j)
        vtens(i,j,1,k,ie)=spheremp(i,j)* &
        (up(i,j,2)*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
        vtens(i,j,2,k,ie)=spheremp(i,j)* &
        (-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))

        ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)
        
        ptens(i,j,k,ie) = spheremp(i,j)*(px(i,j)- &
        fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
        - 0.5*ptens_n(i,j,k,ie)

        vtens(i,j,1,k,ie) = spheremp(i,j)* &
        (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
        - 0.5*vtens_n(i,j,1,k,ie)

        vtens(i,j,2,k,ie) = spheremp(i,j)* &
        (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
        - 0.5*vtens_n(i,j,2,k,ie)
         
     end do
    end do
                        else if (tstep_type==12) then !BDF2 2nd order
    do j=1,np
     do i=1,np
        vtens(i,j,1,k,ie)=spheremp(i,j)* &
        (up(i,j,2)*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
        vtens(i,j,2,k,ie)=spheremp(i,j)* &
        (-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))

        ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)
        if (nstep==0) then ! BE bootstrap
         ptens(i,j,k,ie)   = spheremp(i,j)*(px(i,j)- &
         fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)

         vtens(i,j,1,k,ie) = spheremp(i,j)* &
         (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)

         vtens(i,j,2,k,ie) = spheremp(i,j)* &
         (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)

         !      if (nstep==0) then ! CN bootstrap
         !       ptens(i,j,k,ie) = spheremp(i,j)*(px(i,j)- &
         !          fptr%base(ie)%state%p(i,j,k,n0))*dti - 0.5*ptens(i,j,k,ie) &
         !          - 0.5*ptens_n(i,j,k,ie)
         !
         !       vtens(i,j,1,k,ie) = spheremp(i,j)* &
         !         (up(i,j,1)-ulatlon(i,j,1))*dti - 0.5*vtens(i,j,1,k,ie) &
         !          - 0.5*vtens_n(i,j,1,k,ie)
         !
         !       vtens(i,j,2,k,ie) = spheremp(i,j)* &
         !         (up(i,j,2)-ulatlon(i,j,2))*dti - 0.5*vtens(i,j,2,k,ie) &
         !          - 0.5*vtens_n(i,j,2,k,ie)

        else 

         ptens(i,j,k,ie) = &
         (1+gam)*spheremp(i,j)*(px(i,j) - &
         fptr%base(ie)%state%p(i,j,k,n0))*dti - &
         gam*spheremp(i,j)*(fptr%base(ie)%state%p(i,j,k,n0) - &
         fptr%base(ie)%state%p(i,j,k,nm1))*dti - &
         ptens(i,j,k,ie) 

         vtens(i,j,1,k,ie) = &
         (1+gam)*spheremp(i,j)*(up(i,j,1)-ulatlon(i,j,1))*dti - &
         gam*spheremp(i,j)*(ulatlon(i,j,1)-um(i,j,1))*dti - &
         vtens(i,j,1,k,ie)

         vtens(i,j,2,k,ie) = &
         (1+gam)*spheremp(i,j)*(up(i,j,2)-ulatlon(i,j,2))*dti - &
         gam*spheremp(i,j)*(ulatlon(i,j,2)-um(i,j,2))*dti - &
         vtens(i,j,2,k,ie)
        end if
     end do
    end do
     
   else ! Backward Euler 1st order method   
    do j=1,np
     do i=1,np
        vtens(i,j,1,k,ie)=spheremp(i,j)* &
        (up(i,j,2)*(fcor(i,j) + zeta(i,j)) - grade(i,j,1))
        vtens(i,j,2,k,ie)=spheremp(i,j)* &
        (-up(i,j,1)*(fcor(i,j) + zeta(i,j)) - grade(i,j,2))

        ptens(i,j,k,ie) =  -spheremp(i,j)*div(i,j)
        
        ptens(i,j,k,ie)   = spheremp(i,j)*(px(i,j)- &
        fptr%base(ie)%state%p(i,j,k,n0))*dti - ptens(i,j,k,ie)

        vtens(i,j,1,k,ie) = spheremp(i,j)* &
        (up(i,j,1)-ulatlon(i,j,1))*dti - vtens(i,j,1,k,ie)

        vtens(i,j,2,k,ie) = spheremp(i,j)* &
        (up(i,j,2)-ulatlon(i,j,2))*dti - vtens(i,j,2,k,ie)
     end do
    end do
   endif 
  end do

  ! ===================================================
  ! Pack cube edges of tendencies, rotate velocities
  ! ===================================================
  kptr=0
  call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,ie)
  kptr=nlev
  call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
 end do

   !pw++
   call t_stopf('FE_implicit_KE_resid_calc', t_detail_high)
   !pw--

   !$OMP BARRIER
   !pw++
   call t_startf('FE_implicit_bndry_ex', t_detail_high)
   !pw--
   !$OMP BARRIER
   call bndry_exchangeV(fptr%hybrid,fptr%edge3)
   !$OMP BARRIER
   !pw++
   call t_stopf('FE_implicit_bndry_ex', t_detail_high)
   !pw--
   !$OMP BARRIER

   !pw++
   call t_startf('FE_implicit_bndry_unpack', t_detail_high)
   !pw--
   

   do ie=ns,ne

    ! ===========================================================
    ! Unpack the edges for vgradp and vtens
    ! ===========================================================
    kptr=0
    call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, ie) 

    kptr=nlev
    call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

   end do !ie

   call t_stopf('FE_implicit_bndry_unpack', t_detail_high)

   ! ===========================================================
   ! Compute velocity and pressure tendencies for all levels
   ! ===========================================================

   call t_startf('FE_implicit_vel_pres', t_detail_high)
 

   !Moved conditionals out of do loops
   lx = 0
   if (topology == "cube" .and. test_case=="swtc1") then
    do ie=ns,ne
     do k=1,nlev
      do j=1,np
       do i=1,np
        lx = lx+1
        fx(lx) = 0.0
        fx(lx+shiftv) = 0.0
        fx(lx+shiftp) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
       end do
      end do
     end do
    end do !ie
   else
    do ie=ns,ne
     do k=1,nlev
      do j=1,np
       do i=1,np
        lx = lx+1
        fx(lx)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
        fx(lx+shiftv)=fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
        fx(lx+shiftp) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
       end do
      end do
     end do
    end do !ie
   end if

     call t_stopf('FE_implicit_vel_pres', t_detail_high)
     call t_stopf('FE implicit', t_detail_high)

  end subroutine residual
  
#endif
  
  subroutine precon_gmres(vv, z, nelemd, xstate, c_ptr_to_object, c_ptr_to_pre) &
   bind(C,name='precon')

    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use physical_constants, only : g
    use element_mod, only : element_t
!    use parallel_mod, only : parallel_t
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
!    use hybrid_mod, only : hybrid_t
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use shallow_water_mod, only : tc1_velocity
    use bndry_mod, only : bndry_exchangev
    use cg_mod, only : cg_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use derived_type_mod ,only : derived_type
    use precon_type_mod ,only : precon_type
!    use global_norms_mod

    implicit none

    ! =================
    ! Local
    ! =================

    ! Thread private working set ...

    integer    :: i,j,k,n,ie
    integer    :: nm1,n0,np1
    integer     :: lx, nn
    integer     :: ns, ne

    real (c_double) ,intent(out)       :: vv(nelemd) ! combo of dv and dp 
    real (c_double) ,intent(in)        :: z(nelemd)  ! rhs from lin solver
    real (c_double) ,intent(in)        :: xstate(nelemd)
    integer(c_int)  ,intent(in) ,value :: nelemd
    type(derived_type) ,pointer         :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    type(precon_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre

    real (c_double)       :: zt(nelemd) ! tmp vv to prevent overwrite
    real (c_double)       :: xt(nelemd) ! current xstate to prevent overwrite

  interface 
  
!   subroutine precon_solve(vectorSize,vector,rhs,v_container,p_container) &
!     bind(C,name='precon_solve')
!    use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
!      integer(c_int)                :: vectorSize
!      real(c_double)  ,dimension(*) :: vector
!      real(c_double)  ,dimension(*) :: rhs
!      type(c_ptr)                   :: v_container 
!      type(c_ptr)                   :: p_container 
!    end subroutine precon_solve

!   subroutine precon_solve_wgmres(vectorSize,vector,rhs,v_container,p_container) &
!     bind(C,name='precon_wgmres')
!    use ,intrinsic :: iso_c_binding ,only : c_double ,c_int ,c_ptr
!      integer(c_int)                :: vectorSize
!      real(c_double)  ,dimension(*) :: vector
!      real(c_double)  ,dimension(*) :: rhs
!      type(c_ptr)                   :: v_container 
!      type(c_ptr)                   :: p_container 
!    end subroutine precon_solve

  end interface

    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
    call c_f_pointer(c_ptr_to_pre, pptr) ! convert C ptr to F ptr

    call t_startf('precon_gmres', t_detail_high)

    zt = z
    xt = xstate

! ForTrilinos interface to use applyJacobianInverse with F_l(x)=residual_lin(x)
!    call precon_solve(size(xstate), xt, zt, c_ptr_to_object, c_ptr_to_pre)
! This will point to a GMRES solve by trilinos
!    call precon_precon_wgmres(size(xstate), xt, zt, c_ptr_to_object, c_ptr_to_pre)

    vv = zt

    call t_stopf('precon_gmres', t_detail_high)

    !$OMP BARRIER

  end subroutine precon_gmres




! precon_si() is a preconditioner for the fully implicit solver based on 
! the semi-implicit solver advance_si_nonstag

  subroutine precon_si(vv, z, nelemd, xstate, c_ptr_to_object, c_ptr_to_pre) &
   bind(C,name='precon_si')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use physical_constants, only : rearth, g
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use reduction_mod, only : reductionbuffer_ordered_1d_t
    use time_mod, only : timelevel_t
    use control_mod, only : topology, test_case
    use shallow_water_mod, only : tc1_velocity
    use cg_mod, only : cg_t
    use solver_mod, only : blkjac_t
    use precon_mod, only : pcg_presolver_nonstag
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod ,only : derived_type
    use precon_type_mod ,only : precon_type

    implicit none

    real (c_double) ,intent(out)       :: vv(nelemd) ! combo of dv and dp 
    real (c_double) ,intent(in)        :: z(nelemd)  ! rhs from lin solver
    real (c_double) ,intent(in)        :: xstate(nelemd)
    integer(c_int)  ,intent(in) ,value :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    type(precon_type) ,pointer         :: pptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_pre

    ! =================
    ! Local Variables
    ! =================

    real (kind=real_kind), dimension(np,np)     :: mp
    real (kind=real_kind), dimension(np,np)     :: metdetp
    real (kind=real_kind), dimension(np,np)     :: rmv
    real (kind=real_kind), dimension(np,np)     :: mv
    real (kind=real_kind), dimension(np,np)     :: rspheremp
    real (kind=real_kind), dimension(np,np)     :: spheremp
    real (kind=real_kind), dimension(np,np)     :: metdet
    real (kind=real_kind), dimension(np,np,2,2) :: met
    real (kind=real_kind), dimension(np,np,2,2) :: metinv

    ! Thread private working set ...

    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: dv, zv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: Ru
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: grad_dp
    real (kind=real_kind), dimension(np,np,nlev,nelem)    :: vgradp 
    real (kind=real_kind), dimension(np,np,nlev,nelem)    :: Rs   
    real (kind=real_kind), dimension(np,np,nlev,nelem)    :: dp, zp 
    real (kind=real_kind), dimension(np,np,2) :: gradp,grappm1 ! weak p grad, time level n,n-1
    real (kind=real_kind), dimension(np,np,2) :: gv,gvm1 ! metdet*v(n,n-1), v contravar
    real (kind=real_kind), dimension(np,np,2) :: ulatlon     !(cov) latlon velocity
    real (kind=real_kind), dimension(np,np)   :: E       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)   :: div,divm1 ! tstep n,n-1 velocity div

    real (kind=real_kind) ::  gradp1,gradp2
    real (kind=real_kind) ::  grad_dp1,grad_dp2
    real (kind=real_kind) ::  Ru1,Ru2
 
    real (kind=real_kind) ::  lenscale
    real (kind=real_kind) ::  pmean
    real (kind=real_kind) ::  time_adv

    real*8  :: et,st
    integer i,j,k,ie,shiftv,shiftp
    integer kptr
    integer point
    integer iptr, ieptr
    integer nm1,n0,np1
    integer ns, ne
    integer nn, lx


    call c_f_pointer(c_ptr_to_object, fptr) ! convert C ptr to F ptr
    call c_f_pointer(c_ptr_to_pre,pptr) ! convert C ptr to F ptr

    n0         = pptr%tl%n0
    np1        = pptr%tl%np1
    nm1        = pptr%tl%nm1
    ns         = pptr%nets 
    ne         = pptr%nete 
    pmean      = pptr%pmean

    lenscale=rearth
    call t_startf('precon_si', t_detail_high)

 shiftv = np*np*nlev*(ne-ns+1)
 shiftp = 2*np*np*nlev*(ne-ns+1)
 lx = 0
 do ie=ns,ne
  do k=1,nlev
   do j=1,np
    do i=1,np
     lx = lx+1
     zv(i,j,1,k,ie) = z(lx)
     zv(i,j,2,k,ie) = z(lx+shiftv)
     zp(i,j,k,ie)   = z(lx+shiftp)
    end do  !np
   end do  !np
  end do  !nlev
 end do !ie

    !$OMP BARRIER

    do ie=ns,ne


        rspheremp(:,:)         = pptr%base(ie)%rspheremp(:,:)

       do k=1,nlev

          ! =========================================================
          ! Scale velocity tendency and v.grad(p) term by inverse mass
          ! matrix, scale velocity by metric g factor.
          !     11 np*np Flops
          ! =========================================================

          do j=1,np
             do i=1,np


             gvm1(i,j,1) = pptr%base(ie)%D(i,j,1,1)*zv(i,j,1,k,ie) + pptr%base(ie)%D(i,j,1,2)*zv(i,j,2,k,ie)
             gvm1(i,j,2) = pptr%base(ie)%D(i,j,2,1)*zv(i,j,1,k,ie) + pptr%base(ie)%D(i,j,2,2)*zv(i,j,2,k,ie)

             end do
          end do

          ! ==========================================================
          ! Compute divergence 
          ! ==========================================================

           call divergence_sphere_routine(gvm1,pptr%deriv,pptr%base(ie),divm1)


          do j=1,np
             do i=1,np
                Rs(i,j,k,ie) =        &
                 +  pptr%base(ie)%spheremp(i,j)*(zp(i,j,k,ie)  &
                -  pptr%dt*pmean*divm1(i,j) )
             end do
          end do

       end do
       !DBG print *,'advance_si: point #14'
       kptr=0
       call edgeVpack(pptr%edge1, Rs(1,1,1,ie),nlev,kptr,ie)

    end do

    call bndry_exchangeV(pptr%cg%hybrid,pptr%edge1)
    !$OMP BARRIER

    do ie=ns,ne
       kptr=0
       call edgeVunpack(pptr%edge1, Rs(1,1,1,ie), nlev, kptr, ie)
       do k=1,nlev
          do j=1,np
             do i=1,np
                Rs(i,j,k,ie) = pptr%base(ie)%rspheremp(i,j)*Rs(i,j,k,ie)
             enddo
          enddo
       enddo
    enddo

    ! ======================================================
    ! Invoke the solver!
    ! ======================================================

    !DBG print *,'advance_si: before call to pcg_presolver'
    !$OMP BARRIER

    dp(:,:,:,ns:ne) = pcg_presolver_nonstag(pptr, &
         Rs(:,:,:,ns:ne) )     ! rhs of Helmholtz problem

    do ie=ns,ne

       metinv(:,:,:,:)  = pptr%base(ie)%metinv(:,:,:,:)
       spheremp(:,:)  = pptr%base(ie)%spheremp(:,:)
       rspheremp(:,:)  = pptr%base(ie)%rspheremp(:,:)

       do k=1,nlev

          ! compute grad dp needed to back substitute for du

          call gradient_sphere_routine(dp(:,:,k,ie),pptr%deriv,pptr%base(ie)%Dinv,grad_dp(:,:,:,k,ie))

          ! ==================================================
          ! Rotate grad_dp to form contravariant object:
          !        6 np*np Flops
          ! ==================================================

          do j=1,np
             do i=1,np
                grad_dp(i,j,1,k,ie) = spheremp(i,j)*grad_dp(i,j,1,k,ie) 
                grad_dp(i,j,2,k,ie) = spheremp(i,j)*grad_dp(i,j,2,k,ie) 
             end do
          end do
       enddo

       kptr=0
       call edgeVpack(pptr%edge2, grad_dp(1,1,1,1,ie),2*nlev,kptr,ie)
    end do

    !$OMP BARRIER
    call bndry_exchangeV(pptr%cg%hybrid,pptr%edge2)
    !$OMP BARRIER
    do ie=ns,ne

       kptr=0      
       call edgeVunpack(pptr%edge2, grad_dp(1,1,1,1,ie), 2*nlev, kptr, ie)

       do k=1,nlev

          ! ==============================
          ! Update velocity
          !    16 np*np Flops
          ! ==============================

          do j=1,np
             do i=1,np
               grad_dp1 = pptr%base(ie)%rspheremp(i,j)*grad_dp(i,j,1,k,ie)
               grad_dp2 = pptr%base(ie)%rspheremp(i,j)*grad_dp(i,j,2,k,ie)
           grad_dp(i,j,1,k,ie)   = pptr%base(ie)%Dinv(i,j,1,1)*grad_dp1 + &
                      pptr%base(ie)%Dinv(i,j,1,2)*grad_dp2
           grad_dp(i,j,2,k,ie)   = pptr%base(ie)%Dinv(i,j,2,1)*grad_dp1 + &
                      pptr%base(ie)%Dinv(i,j,2,2)*grad_dp2

               dv(i,j,1,k,ie) = zv(i,j,1,k,ie) - pptr%dt*grad_dp(i,j,1,k,ie)
               dv(i,j,2,k,ie) = zv(i,j,2,k,ie) - pptr%dt*grad_dp(i,j,2,k,ie)  


             end do
          end do
       end do

    end do ! ie loop

    lx = 0
 if (topology == "cube" .and. test_case=="swtc1") then
  do ie=ns,ne
   do k=1,nlev
    do j=1,np
     do i=1,np
      lx = lx+1
      vv(lx) = zv(i,j,1,k,ie)
      vv(lx+shiftv) = zv(i,j,2,k,ie)
      vv(lx+shiftp) = dp(i,j,k,ie)
     end do  !np
    end do  !np
   end do  !nlev
  end do !ie  
 else
  do ie=ns,ne
   do k=1,nlev
    do j=1,np
     do i=1,np
      lx = lx+1
      vv(lx) = dv(i,j,1,k,ie)
      vv(lx+shiftv) = dv(i,j,2,k,ie)
      vv(lx+shiftp) = dp(i,j,k,ie)
     end do  !np
    end do  !np
   end do  !nlev
  end do !ie
   end if
    !$OMP BARRIER
    call t_stopf('precon_si', t_detail_high)

  end subroutine precon_si



  subroutine sw_picard_schur_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_schur')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case,tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use physical_constants, only : g

    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    real (c_double)       :: fxtemp(2*nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1

    real (kind=real_kind), dimension(np,np,2)  :: grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term

    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv

    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam

    integer    :: i,j,k,n,ie
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx

    call t_startf('Precon Schur', t_detail_high)


    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr


    n0     = fptr%tl%n0
    np1    = fptr%tl%np1
    nm1    = fptr%tl%nm1
    dti    = 1.0d0/fptr%dt
    dt     = fptr%dt
    ns     = fptr%nets 
    ne     = fptr%nete 
    pmean  = fptr%pmean
    gam    = 0.5

    if (tstep_type==13) then ! CN
       dti=2.0d0*dti ! CN has a factor of 3/2 in the time step coefficient
    endif
    
    if (tstep_type==12) then ! BDF2 
       dti=(3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 in the time step coefficient
    endif

    lx = 0
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%p(i,j,k,np1) = xs(lx)!p =g*h in notes
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie

    do ie=ns,ne

       do j=1,np
          do i=1,np
             mv(i,j)        = fptr%base(ie)%mp(i,j)
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
             ps(i,j)        = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
          end do !np
       end do !np

       do k=1,nlev

          ! scalar -> latlon vector        (grad)del phi
          call gradient_sphere_routine(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv,grade_n3)

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                vtens(i,j,1,k,ie)= spheremp(i,j)*( &
                     grade_n3(i,j,1)) ! + d(delp)/dx

                vtens(i,j,2,k,ie)= spheremp(i,j)*( &
                     grade_n3(i,j,2))!d(delp)/dy

             end do !np
          end do !np
       end do !nlev


       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

    end do !ie


    !!apply rspheremp
    !
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                vtens(i,j,1,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do !ie
    ! end of Bt cal


    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================
    !F_diag_inv(B'p) (B'p is in vtens and matches along element interfaces so F_diag_inv can be applied locally )

    !  ApplyFDiagInv 
    lx = 0

    if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fptr%base(ie)%state%v(i,j,1,k,np1) = 0.0
                   fptr%base(ie)%state%v(i,j,2,k,np1) = 0.0
                end do
             end do
          end do
       end do
    else
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fptr%base(ie)%state%v(i,j,1,k,np1) = dt*vtens(i,j,1,k,ie)
                   fptr%base(ie)%state%v(i,j,2,k,np1) = dt*vtens(i,j,2,k,ie)
                end do
             end do
          end do
       end do
    end if

    ! we now need to apply B to vtens, as well as G to dh (a.k.a. fptr%base(ie)%state%p(:,:,k,np1), and sum these two quantities
    ! now compute Gdp -BF_diag_inv(B'p) 
    !
    do ie=ns,ne
       do j=1,np
          do i=1,np
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
          end do !np
       end do !np

       do k=1,nlev
          !
          !         ! ==============================================
          !         ! Compute kinetic energy term at each time level
          !         ! ==============================================
          !
          do j=1,np
             do i=1,np


                ! set u
                ulatlon(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
                ulatlon(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,n0)   !  

                ! set du
                up(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
                up(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,np1)   !  

                dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
                dh(i,j)   = fptr%base(ie)%state%p(i,j,k,np1)

                pv1(i,j,1) = up(i,j,1)
                pv1(i,j,2) = 0.d00
                pv2(i,j,1) = 0.d00
                pv2(i,j,2) = up(i,j,2)

                dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )  
                dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )   

                dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
                dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )

             end do !np
          end do !np

          !! create operators for full residual calc
          call divergence_sphere_routine(pv1,fptr%deriv,fptr%base(ie),divpv1) ! latlon vector -> scalar 
          call divergence_sphere_routine(pv2,fptr%deriv,fptr%base(ie),divpv2) ! latlon vector -> scalar 
          call divergence_sphere_routine(dpv,fptr%deriv,fptr%base(ie),divdpv) ! latlon vector -> scalar 
          !
          !
          !          ! ==============================================
          !          ! Compute residual terms
          !          ! ==============================================
          !
          !
          !
          do j=1,np
             do i=1,np
                !
                ptens(i,j,k,ie) = spheremp(i,j)*(&
                                            !G
                     dh(i,j)*dti+&                         !delp/dt
                     divdpv(i,j) - &
                                !-BdiagFinvB'
                     (  ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
                     (divpv1(i,j)+divpv2(i,j))))

             end do !np
          end do !np
       end do !nlev

       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge1, ptens(1,1,1,ie),nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge1)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge1, ptens(1,1,1,ie), nlev, kptr, ie)
    end do !ie

    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0
    !    do n=1,nvar
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                !              if (n==1) fx(lx) = 0.0d0
                !              if (n==2) fx(lx) = 0.0d0
                fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
             end do
          end do
       end do
    end do !ie
    !     end do
    
    call t_stopf('Precon Schur', t_detail_high)

  end subroutine sw_picard_schur_op



  subroutine sw_picard_DFinvBt_op(xs, nelemd,fx,nret, c_ptr_to_object) bind(C,name='sw_picard_DFinvBt')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use physical_constants, only : g

    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    integer(c_int) ,intent(in) ,value  :: nret
    real (c_double) ,intent(out)       :: fx(nret)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: spheremp,rspheremp,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens

    real (kind=real_kind), dimension(np,np,2)  :: grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np)    :: zeta_n ! relative vorticity

    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv

    real (kind=real_kind) ::  dt,dti, pmean
    real (kind=real_kind) ::  gam, tot_nv

    integer    :: i,j,k,n,ie, shiftv
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx

    call t_startf('Precon DFinvBt', t_detail_high)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0    = fptr%tl%n0
    np1   = fptr%tl%np1
    nm1   = fptr%tl%nm1
    dti   = 1.0d0/fptr%dt
    dt    = fptr%dt
    ns    = fptr%nets 
    ne    = fptr%nete 
    pmean = fptr%pmean
    gam   = 0.5

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    !     do n=1,nvar
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%p(i,j,k,np1) = xs(lx)!p =g*h in notes
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie

    do ie=ns,ne

       do j=1,np
          do i=1,np
             mv(i,j)          = fptr%base(ie)%mp(i,j)
             spheremp(i,j)    = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j)   = fptr%base(ie)%rspheremp(i,j)
          end do !np
       end do !np

       do k=1,nlev

          call gradient_sphere_routine(fptr%base(ie)%state%p(:,:,k,np1),fptr%deriv,fptr%base(ie)%Dinv,grade_n3)  ! scalar -> latlon vector        (grad)del phi

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                vtens(i,j,1,k,ie)= spheremp(i,j)*( &
                     grade_n3(i,j,1)) ! + d(delp)/dx

                vtens(i,j,2,k,ie)= spheremp(i,j)*( &
                     grade_n3(i,j,2))!d(delp)/dy

             end do !np
          end do !np
       end do !nlev


       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge2, vtens(1,1,1,1,ie),2*nlev,kptr, ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge2)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge2, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

    end do !ie

    ! apply rspheremp
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                vtens(i,j,1,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                vtens(i,j,2,k,ie) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
             end do
          end do
       end do
    end do !ie
    ! end of Bt cal


    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0

    if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx) = 0.0
                   fx(lx+shiftv) = 0.0
                end do
             end do
          end do
       end do !ie
    else
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx)        = dt*vtens(i,j,1,k,ie)
                   fx(lx+shiftv) = dt*vtens(i,j,2,k,ie)
                end do
             end do
          end do
       end do !ie 
    end if

    call t_stopf('Precon DFinvBt', t_detail_high)

  end subroutine sw_picard_DFinvBt_op



  subroutine sw_picard_block_11(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_picard_block_11')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case, tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use physical_constants, only : g

    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1

    real (kind=real_kind), dimension(np,np)    :: zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity

    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2

    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, tot_nv

    integer    :: i,j,k,n,ie,shiftv
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx

    call t_startf('precon 11', t_detail_high)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0    = fptr%tl%n0
    np1   = fptr%tl%np1
    nm1   = fptr%tl%nm1
    dti   = 1.0d0/fptr%dt
    ns    = fptr%nets 
    ne    = fptr%nete 
    pmean = fptr%pmean
    gam   = 0.5

    if (tstep_type==13) then ! CN
       dti=2.0d0*dti ! CN has a factor of 3/2 in the time step coefficient
    endif
    
    if (tstep_type==12) then ! BDF2 
       dti=(3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 in the time step coefficient
    endif

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)

    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx+shiftv)
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie

    do ie=ns,ne

       do j=1,np
          do i=1,np
             fcor(i,j)      = fptr%base(ie)%fcor(i,j)
             mv(i,j)        = fptr%base(ie)%mp(i,j)
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
          end do !np
       end do !np

       do k=1,nlev

          ! ==============================================
          ! Compute kinetic energy term at each time level
          ! ==============================================

          do j=1,np
             do i=1,np

                ! set u
                ulatlon(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
                ulatlon(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,n0)   !  
          
                ! set du
                up(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
                up(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,np1)   !  

                delv1(i,j,1) = up(i,j,1)
                delv1(i,j,2) = 0.0d0
                delv2(i,j,1) = 0.0d0
                delv2(i,j,2) = up(i,j,2)

             end do !np
          end do !np

          ! create operators for full residual calc
          call vorticity_sphere_routine(ulatlon,fptr%deriv,fptr%base(ie),zeta_n) ! latlon vector -> scalar  w
          call vorticity_sphere_routine(delv2,fptr%deriv,fptr%base(ie),newton_vortdelv2) ! latlon vector -> scalar  w
          call vorticity_sphere_routine(delv1,fptr%deriv,fptr%base(ie),newton_vortdelv1) ! latlon vector -> scalar  w

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                !picard linearization of 11 block 
                !note this 11 block is representing a 2x2 system for up1 and up2 

                vtens(i,j,1,k,ie) = spheremp(i,j)*( &
                     up(i,j,1)*dti - &                      ! delv1/dt
                     ulatlon(i,j,2)*newton_vortdelv1(i,j)-& ! u2*d(delv1)/dy & dx
                     up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-&    ! delv2*(w+f)
                     ulatlon(i,j,2)*newton_vortdelv2(i,j))  ! v2*d(delv2)/dx &dy

                vtens(i,j,2,k,ie) = spheremp(i,j)*( &
                     up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+&    ! delu1*(w+f)
                     ulatlon(i,j,1)*newton_vortdelv1(i,j)+& ! u1*d(delu1)/dx & dy
                     up(i,j,2)*dti +&                       ! delu2/dt
                     ulatlon(i,j,1)*newton_vortdelv2(i,j))  ! u1*d(delu2)/dy & dx

             end do !np
          end do !np
       end do !nlev

       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,ie)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, ie)
       kptr=nlev
       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

    end do !ie
    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0

    if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx) = 0.0
                   fx(lx+shiftv) = 0.0
                end do
             end do
          end do
       end do
    else
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx)        = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                   fx(lx+shiftv) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
                end do
             end do
          end do
       end do
    end if

    call t_stopf('precon 11', t_detail_high)

  end subroutine sw_picard_block_11



  subroutine sw_picard_block_21(xs, nelemd,fx, nret, c_ptr_to_object) bind(C,name='sw_picard_block_21')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use physical_constants, only : g

    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    integer(c_int) ,intent(in) ,value  :: nret
    real (c_double) ,intent(out)       :: fx(nret)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne

    real (kind=real_kind), dimension(np,np) :: spheremp,rspheremp,mv
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1

    real (kind=real_kind), dimension(np,np)    :: zeta_n ! relative vorticity

    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv

    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, tot_nv

    integer    :: i,j,k,n,ie, shiftv
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx

    call t_startf('picard 21', t_detail_high)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0    = fptr%tl%n0
    np1   = fptr%tl%np1
    nm1   = fptr%tl%nm1
    dti   = 1.0d0/fptr%dt
    ns    = fptr%nets 
    ne    = fptr%nete 
    pmean = fptr%pmean
    gam   = 0.5

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx+shiftv)
             end do  !np
          end do  !np
       end do  !nlev
    end do !ie

    do ie=ns,ne

       do j=1,np
          do i=1,np
             mv(i,j)        = fptr%base(ie)%mp(i,j)
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
          end do !np
       end do !np

       do k=1,nlev

          ! ==============================================
          ! Compute kinetic energy term at each time level
          ! ==============================================

          do j=1,np
             do i=1,np

                ! set u
                ulatlon(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,n0)   ! 
                ulatlon(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,n0)   !  

                ! set du
                up(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,np1)   ! 
                up(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,np1)   !  

                ! set dh_n=p_notes=g(h+hmean)=phi
                !          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
                !          dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)

                pv1(i,j,1) = up(i,j,1)
                pv1(i,j,2) = 0.d00
                pv2(i,j,1) = 0.d00
                pv2(i,j,2) = up(i,j,2)

             end do !np
          end do !np

          ! create operators for full residual calc
          call divergence_sphere_routine(pv1,fptr%deriv,fptr%base(ie),divpv1) ! latlon vector -> scalar 
          call divergence_sphere_routine(pv2,fptr%deriv,fptr%base(ie),divpv2) ! latlon vector -> scalar 

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                ptens(i,j,k,ie) = spheremp(i,j)*(&
                     ( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )*&
                     (divpv1(i,j)+divpv2(i,j)))

             end do !np
          end do !np
       end do !nlev


       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, ie)

    end do !ie
    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fx(lx) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie) 
             end do
          end do
       end do
    end do !ie

    call t_stopf('picard 21', t_detail_high)

  end subroutine sw_picard_block_21
  
  
  
  subroutine sw_jacobian_op(xs, nelemd,fx, c_ptr_to_object) bind(C,name='sw_jacobian')
    
    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use edge_mod, only : edgevpack, edgerotate, edgevunpack
    use edgetype_mod, only : EdgeBuffer_t
    use hybrid_mod, only : hybrid_t
    use time_mod, only : timelevel_t
    use control_mod, only :  topology, test_case, tstep_type
    use bndry_mod, only : bndry_exchangev
    use derived_type_mod, only : derived_type
    use physical_constants, only : g
    
    implicit none
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    !    real (c_double)        :: base_vector(nelemd)
    integer              :: ns
    integer              :: ne
    
    real (kind=real_kind), dimension(np,np) :: fcor,spheremp,rspheremp,ps,mv
    real (kind=real_kind), dimension(np,np,2,nlev,nelem)  :: vtens
    real (kind=real_kind), dimension(np,np,nlev,nelem) :: ptens ! at t=n+1
    
    real (kind=real_kind), dimension(np,np,2)  :: grade_n1,grade_n2,grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np,2)  :: newton_grade,newton_grade_n1,newton_grade_n2,newton_grade_n3 ! strong KE gradient
    real (kind=real_kind), dimension(np,np)    :: dh,dh_n       ! kinetic energy term
    real (kind=real_kind), dimension(np,np)    :: zeta_n ! relative vorticity
    real (kind=real_kind), dimension(np,np)    :: newton_vortdelv1,newton_vortdelv2 ! relative vorticity
    
    real (kind=real_kind), dimension(np,np)    :: divpv1,divpv2,divdpv
    real (kind=real_kind), dimension(np,np)    :: div_h,div_hn   ! at t=n+1,n  
    real (kind=real_kind), dimension(np,np,2)  :: ulatlon, up
    real (kind=real_kind), dimension(np,np,2)  :: pv1,pv2,dpv,delv1,delv2
    !    real (kind=real_kind), dimension(np,np)     :: metdet
    
    real (kind=real_kind) ::  dti, pmean
    real (kind=real_kind) ::  gam, tot_nv
    
    integer    :: i,j,k,n,ie, shiftv, shiftp
    integer    :: kptr
    integer    :: nm1,n0,np1
    integer    :: lx
    
    !call flush(6)
    !write(6,*)'sw startf'
    !call flush(6)
    
    call t_startf('sw jacobian op', t_detail_high)
  
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    
    n0    = fptr%tl%n0
    np1   = fptr%tl%np1
    nm1   = fptr%tl%nm1
    dti   = 1.0d0/fptr%dt
    ns   = fptr%nets 
    ne   = fptr%nete 
    pmean = fptr%pmean
    gam   = 0.5
    
    if (tstep_type==13) then ! CN
       dti=2.0d0*dti ! CN has a factor of 3/2 in the time step coefficient
    endif
    
    if (tstep_type==12) then ! BDF2 
       dti=(3.0d0/2.0d0)*dti ! BDF2 has a factor of 3/2 in the time step coefficient
    endif

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    shiftp = 2*np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,np1) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,np1) = xs(lx+shiftv)
                fptr%base(ie)%state%p(i,j,k,np1)   = xs(lx+shiftp)
             end do !np
          end do !np
       end do !nlev
    end do !ie
    
    do ie=ns,ne
       
       do j=1,np
          do i=1,np
             fcor(i,j)      = fptr%base(ie)%fcor(i,j)
             mv(i,j)        = fptr%base(ie)%mp(i,j)
             spheremp(i,j)  = fptr%base(ie)%spheremp(i,j)
             rspheremp(i,j) = fptr%base(ie)%rspheremp(i,j)
             ps(i,j)        = fptr%base(ie)%state%ps(i,j)! surface height g*h_s in notes
          end do !np
       end do !np
       
       do k=1,nlev
          
          ! ==============================================
          ! Compute kinetic energy term at each time level
          ! ==============================================
          
          do j=1,np
             do i=1,np
                
                ! set u
                ulatlon(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,n0)
                ulatlon(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,n0)
                
                ! set du
                up(i,j,1) = fptr%base(ie)%state%v(i,j,1,k,np1)
                up(i,j,2) = fptr%base(ie)%state%v(i,j,2,k,np1)
                
                delv1(i,j,1) = up(i,j,1)
                delv1(i,j,2) = 0.0d0
                delv2(i,j,1) = 0.0d0
                delv2(i,j,2) = up(i,j,2)
                
                !test
                !         up(i,j,1)=v1p
                !          up(i,j,2)=v2p
                !end test
                
                ! where is ps???? shouldn't ps be included in dp? ps=g*h_s in comment above, 
                ! ANSWER: NO, also no pmean, otherwise we would be adding these contributions at each nonlinear update 
                !what is p 
                !what is pmean
                !ps= surface height g*h_s in notes
                
                ! set dh_n=p_notes=g(h+hmean)=phi
                dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)+pmean +ps(i,j)
                !test unit vector
                !          dh_n(i,j) = fptr%base(ie)%state%p(i,j,k,n0)
                
                ! set dp
                !dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)+ps(i,j)+pmean
                !dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)+pmean
                dh(i,j) = fptr%base(ie)%state%p(i,j,k,np1)
                
                pv1(i,j,1) = up(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
                pv1(i,j,2) = 0.d00
                
                pv2(i,j,1) = 0.d00
                pv2(i,j,2) = up(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,n0) )
                
                dpv(i,j,1) = ulatlon(i,j,1)*( fptr%base(ie)%state%p(i,j,k,np1) )
                dpv(i,j,2) = ulatlon(i,j,2)*( fptr%base(ie)%state%p(i,j,k,np1) )
                
                !dpv(i,j,1) = ulatlon(i,j,1)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1))
                !dpv(i,j,2) = ulatlon(i,j,2)*( fptr%pmean + fptr%base(ie)%state%p(i,j,k,np1))
                
             end do !np
          end do !np
          
          ! create operators for full residual calc
          call gradient_sphere_routine(up(:,:,1),fptr%deriv,fptr%base(ie)%Dinv,grade_n1)  ! scalar -> latlon vector (grad)delu_1
          call gradient_sphere_routine(up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv,grade_n2)  ! scalar -> latlon vector (grad)delu_2
          call gradient_sphere_routine(dh,fptr%deriv,fptr%base(ie)%Dinv,grade_n3)  ! scalar -> latlon vector        (grad)del phi
                    
          call gradient_sphere_routine(ulatlon(:,:,1)*up(:,:,1) ,fptr%deriv,fptr%base(ie)%Dinv,newton_grade_n1)  ! scalar -> latlon vector (grad)u_1
          call gradient_sphere_routine(ulatlon(:,:,2)*up(:,:,2),fptr%deriv,fptr%base(ie)%Dinv,newton_grade_n2)  ! scalar -> latlon vector (grad)u_2
          call gradient_sphere_routine(dh_n,fptr%deriv,fptr%base(ie)%Dinv,newton_grade_n3)  ! scalar -> latlon vector        (grad)phi
          
          call vorticity_sphere_routine(ulatlon,fptr%deriv,fptr%base(ie),zeta_n) ! latlon vector -> scalar  w
          call vorticity_sphere_routine(delv2,fptr%deriv,fptr%base(ie),newton_vortdelv2) ! latlon vector -> scalar  w
          call vorticity_sphere_routine(delv1,fptr%deriv,fptr%base(ie),newton_vortdelv1) ! latlon vector -> scalar  w
          
          call divergence_sphere_routine(pv1,fptr%deriv,fptr%base(ie),divpv1) ! latlon vector -> scalar 
          call divergence_sphere_routine(pv2,fptr%deriv,fptr%base(ie),divpv2) ! latlon vector -> scalar 
          call divergence_sphere_routine(dpv,fptr%deriv,fptr%base(ie),divdpv) ! latlon vector -> scalar 

          ! ==============================================
          ! Compute residual terms
          ! ==============================================

          do j=1,np
             do i=1,np

                ! Jacobian
                vtens(i,j,1,k,ie) = spheremp(i,j)*( &
                     up(i,j,1)*dti - &                      ! delv1/dt
                     ulatlon(i,j,2)*newton_vortdelv1(i,j)+& ! u2*d(delv1)/dy & dx
                     newton_grade_n1(i,j,1)- &              ! delv1*d(v1)/dx
                     up(i,j,2)*(zeta_n(i,j)+fcor(i,j))-&    ! delv2*(w+f)
                     ulatlon(i,j,2)*newton_vortdelv2(i,j)+& ! v2*d(delv2)/dx &dy
                     newton_grade_n2(i,j,1)+ &              ! -delv2*d(v2)/dx
                     grade_n3(i,j,1))                       ! + d(delp)/dx

                vtens(i,j,2,k,ie) = spheremp(i,j)*( &
                     up(i,j,1)*(zeta_n(i,j)+fcor(i,j))+&    ! delu1*(w+f)
                     ulatlon(i,j,1)*newton_vortdelv1(i,j)+& ! u1*d(delu1)/dx & dy
                     newton_grade_n2(i,j,2)+&               ! delu1*d(u2)/dy
                     up(i,j,2)*dti +&                       ! delu2/dt
                     ulatlon(i,j,1)*newton_vortdelv2(i,j)+& ! u1*d(delu2)/dy & dx
                     newton_grade_n1(i,j,2)+ &              ! delu1*d(u1)/dy
                     grade_n3(i,j,2))                       ! d(delp)/dy

                ptens(i,j,k,ie) = spheremp(i,j)*(&
                     dh(i,j)*dti+&                          ! delp/dt
                     divpv1(i,j)+divpv2(i,j)+divdpv(i,j))
                
             end do !np
          end do !np
       end do !nlev

       ! ===================================================
       ! Pack cube edges of tendencies, rotate velocities
       ! ===================================================
       kptr=0
       call edgeVpack(fptr%edge3, ptens(1,1,1,ie),nlev,kptr,ie)
       kptr=nlev
       call edgeVpack(fptr%edge3, vtens(1,1,1,1,ie),2*nlev,kptr,ie)
    end do !ie

    !$OMP BARRIER
    call bndry_exchangeV(fptr%hybrid,fptr%edge3)
    !$OMP BARRIER

    do ie=ns,ne

       ! ===========================================================
       ! Unpack the edges for vgradp and vtens
       ! ===========================================================
       kptr=0
       call edgeVunpack(fptr%edge3, ptens(1,1,1,ie), nlev, kptr, ie)
       kptr=nlev
       call edgeVunpack(fptr%edge3, vtens(1,1,1,1,ie), 2*nlev, kptr, ie)

    end do !ie
    ! ===========================================================
    ! Compute velocity and pressure tendencies for all levels
    ! ===========================================================

    lx = 0
    if (topology == "cube" .and. test_case=="swtc1") then
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx) = 0.0
                   fx(lx+shiftv) = 0.0
                end do
             end do
          end do
       end do
    else
       do ie=ns,ne
          do k=1,nlev
             do j=1,np
                do i=1,np
                   lx = lx+1
                   fx(lx)        = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,1,k,ie)
                   fx(lx+shiftv) = fptr%base(ie)%rspheremp(i,j)*vtens(i,j,2,k,ie)
                   fx(lx+shiftp) = fptr%base(ie)%rspheremp(i,j)*ptens(i,j,k,ie)
                end do
             end do
          end do
       end do
    end if

    call t_stopf('sw jacobian op', t_detail_high)

  end subroutine sw_jacobian_op



  subroutine update_jac_state(xs, nelemd, c_ptr_to_object) bind(C,name='update_jac_state')

    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use derived_type_mod, only : derived_type

    implicit none

    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer                            :: ns
    integer                            :: ne

    integer :: i, j, k, n, ie, shiftv, shiftp
    integer :: nm1, n0, np1
    integer :: lx

    ! write(6,*)'update jac state\n'
    ! call t_startf('update jac', t_detail_high)

    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr

    n0  = fptr%tl%n0
    np1 = fptr%tl%np1
    nm1 = fptr%tl%nm1
    ns = fptr%nets 
    ne = fptr%nete 

    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    shiftp = 2*np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,n0) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,n0) = xs(lx+shiftv)
                fptr%base(ie)%state%p(i,j,k,n0)   = xs(lx+shiftp)
             end do ! np
          end do ! np
       end do ! nlev
    end do !ie

    ! call t_stopf('update jac', t_detail_high)

  end subroutine update_jac_state



  subroutine update_prec_state(xs, nelemd, c_ptr_to_object) bind(C,name='update_prec_state')
    
    use ,intrinsic :: iso_c_binding 
    use kinds, only : real_kind
    use dimensions_mod, only : np, nlev, nvar, nelem
    use element_mod, only : element_t
    use derived_type_mod, only : derived_type
    
    implicit none
    
    real (c_double) ,intent(in)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne
    
    integer    :: i,j,k,n,ie, shiftv, shiftp
    integer    :: nm1,n0,np1
    integer    :: lx

    !write(6,*)'update prec state\n'
    !    call t_startf('residual lin', t_detail_high)
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    
    n0  = fptr%tl%n0
    np1 = fptr%tl%np1
    nm1 = fptr%tl%nm1
    ns  = fptr%nets 
    ne  = fptr%nete 
    
    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    shiftp = 2*np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                fptr%base(ie)%state%v(i,j,1,k,n0) = xs(lx)
                fptr%base(ie)%state%v(i,j,2,k,n0) = xs(lx+shiftv)
                fptr%base(ie)%state%p(i,j,k,n0)   = xs(lx+shiftp)
             end do !np
          end do !np
       end do !nlev
    end do !ie
    
  end subroutine update_prec_state
  
  
  
  subroutine get_jac_vector(xs, nelemd, c_ptr_to_object) bind(C,name='get_jac_vector')
    
    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use derived_type_mod, only : derived_type
    
    implicit none
    
    real (c_double)        :: xs(nelemd)
    integer(c_int) ,intent(in) ,value  :: nelemd
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    integer              :: ns
    integer              :: ne
    
    integer    :: i,j,k,n,ie,shiftv,shiftp
    integer    :: nm1,n0,np1
    integer    :: lx
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    
    n0 = fptr%tl%n0
    ns = fptr%nets
    ne = fptr%nete
    
    lx = 0
    shiftv = np*np*nlev*(ne-ns+1)
    shiftp = 2*np*np*nlev*(ne-ns+1)
    do ie=ns,ne
       do k=1,nlev
          do j=1,np
             do i=1,np
                lx = lx+1
                xs(lx)        = fptr%base(ie)%state%v(i,j,1,k,n0)
                xs(lx+shiftv) = fptr%base(ie)%state%v(i,j,2,k,n0)
                xs(lx+shiftp) = fptr%base(ie)%state%p(i,j,k,n0) 
             end do !np
          end do !np
       end do !nlev
    end do !ie
    
  end subroutine get_jac_vector
  

  
  subroutine test_id(xs, nelemd, fx, c_ptr_to_object) bind(C,name='test_id')
    use ,intrinsic :: iso_c_binding 
    
    implicit none 
    
    integer(c_int) ,intent(in) ,value  :: nelemd
    real (c_double) ,intent(in)        :: xs(nelemd)
    real (c_double) ,intent(out)       :: fx(nelemd)
    type(c_ptr)                        :: c_ptr_to_object
    
    fx=xs
    
  end subroutine test_id
  
  
  
  subroutine get_discrete_params(cnets,cnete,cnp,cnlev, c_ptr_to_object) bind(C,name='get_discrete_params')
    
    use ,intrinsic :: iso_c_binding
    use kinds, only : real_kind
    use dimensions_mod, only :  np, nlev, nvar, nelem
    use element_mod, only : element_t
    use derived_type_mod, only : derived_type
    
    implicit none
    
    integer (c_int)        :: cnets
    integer (c_int)        :: cnete
    integer (c_int)        :: cnp
    integer (c_int)        :: cnlev
    type(derived_type) ,pointer        :: fptr=>NULL()
    type(c_ptr)                        :: c_ptr_to_object
    
    call c_f_pointer(c_ptr_to_object,fptr) ! convert C ptr to F ptr
    
    cnets         = fptr%nets
    cnete         = fptr%nete
    
    cnp=np
    cnlev=nlev
    
  end subroutine get_discrete_params
  
end module implicit_mod
