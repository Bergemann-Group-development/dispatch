! ---------------------------------------------------------------
!  MAG_UNSPLIT Unsplit second order Godunov integrator for
!              polytropic magnetized gas dynamics using 
!              MUSCL-HANCOCK scheme
!              with various Riemann solvers and slope limiters.
!              The sheme follows closely the paper by
!              Londrillo & Del Zanna ApJ 2000, 530, 508, 
!
!  inputs/outputs
!  uin         => (const)  input state
!  gravin      => (const)  input gravitational acceleration
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,    
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
!
!  uin = (\rho, \rho u, \rho v, \rho w, Etot, A, B, C)
!  the hydro variable are cell-centered 
!  whereas the magnetic field B=(A,B,C) are face-centered.
!  Note that here we have 3 components for v and B whatever ndim.
!
!  This routine was written by Sebastien Fromang and Patrick Hennebelle
! ----------------------------------------------------------------
subroutine mag_unsplit(uin,gravin,flux,emfx,emfy,emfz,tmp,mark_low_order,mark_low_order_edge,dx,dy,dz,dt,ngrid,ind_grid)
  use amr_parameters
  use amr_commons, only: dbg_grid, print_id, xg, local_debug, myid
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  integer ::ind_grid(nvector)
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin 

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp 

  ! Output electromotive force
  REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3)::emfx
  REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3)::emfy
  REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3)::emfz

  ! Markers that tell if we should default to llf
  logical ,dimension(1:nvector,twotondim,twondim):: mark_low_order      ! This is for marking cell interfaces
  logical ,dimension(1:nvector,1:3,1:3,1:3,3)    :: mark_low_order_edge ! This is for marking edges

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+1)::qin 
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  

  ! Cell-centered slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::dS

  ! Corner-centered slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq2

  ! Face-centered slopes
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf

  ! Face-averaged left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp
  
  ! Edge-averaged left-right and top-bottom state arrays
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLB

  ! Intermediate fluxes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::fx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2   )::tx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)       ::emf

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi
  real(dp)::dtdx,dtdy,dtdz,a2,ctot,idim
  logical in_debug_region
  character(len=80):: id='umuscl.f90 $Id$'
  logical, save:: first_time=.true.

  call print_id(id)
!  call trace_mpi('mag_unsplit','enter',2)

!  if (local_debug) then
!     call dump_debug_cube('umuscl','uin',uin(1:ngrid,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar), &
!         iu2-iu1+1,ju2-ju1+1,ku2-ku1+1,ngrid,nvar,ndim,ind_grid)
!  end if

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  dtdx=dt/dx
  dtdy=dt/dy
  dtdz=dt/dz

  ! Translate to primitive variables, compute sound speeds  
  call ctoprim(uin,qin,bf,gravin,dt,ngrid)

  ! Compute TVD slopes
  call uslope(bf,qin,dq,dbf,dx,dt,nvector,slope_type)
  if (slope_type==3.0) then
     dq2=dq
  else
     call uslope(bf,qin,dq2,dbf,dx,dt,nvector,3.0_dp) ! Use slope_type=3 for corner states
  endif

  ! Compute peak conserving slopes of entropy
  call aslope(qin(:,:,:,:,nvar+1),dS,nvector,slope_type)

  ! Compute 3D traced-states in all three directions
#if NDIM==1
     call trace1d(qin   ,dq    ,qm,qp                ,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call trace2d(qin,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call trace3d(qin,bf,dq,dq2,dbf,dS,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dz,dt,nvector,ind_grid)
#endif

  ! Pick up the advected entropy value, for testing in god_fine1
  if (entropy_fix) then
     uin(:,:,:,:,5) = qin(:,:,:,:,5)
  end if

  ! Solve for 1D flux in X direction
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,6,7,8, &
       &       mark_low_order,fx,tx,ngrid,ind_grid)
!  if (local_debug) then
!     call dump_debug_cube('umuscl','x-flux',fx(1:ngrid,if1:if2,jlo:jhi,klo:khi,1:nvar), &
!         if2-if1+1,jhi-jlo+1,khi-klo+1,ngrid,nvar,ndim,ind_grid)
!  end if
  ! Save flux in output array
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,1)=fx(l,i,j,k,ivar)*dtdx
           !if (in_debug_region(ind_grid(l)).and.ivar==5) then
           !  print '(a,i8,5i3,1p,g12.4)','flx1',ind_grid(l),l,i,j,k,1,flux(l,i,j,k,ivar,1)
           !end if
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,1)=tx(l,i,j,k,ivar)*dtdx
        end do
     end do
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 3,2,4,7,6,8, &
       &       mark_low_order,fx,tx,ngrid,ind_grid)
  if (local_debug) then
     call dump_debug_cube('umuscl','y-flux',fx(1:ngrid,ilo:ihi,jf1:jf2,klo:khi,1:nvar), &
         ihi-ilo+1,jf2-jf1+1,khi-klo+1,ngrid,nvar,ndim,ind_grid)
  end if
  ! Save flux in output array
  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,2)=fx(l,i,j,k,ivar)*dtdy
           !if (in_debug_region(ind_grid(l)).and.ivar==5) then
           !  print '(a,i8,5i3,1p,g12.4)','flx2',ind_grid(l),l,i,j,k,2,flux(l,i,j,k,ivar,2)
           !end if
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,2)=tx(l,i,j,k,ivar)*dtdy
        end do
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
#if NDIM==3
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,8,6,7, &
       &       mark_low_order,fx,tx,ngrid,ind_grid)
  if (local_debug) then
     call dump_debug_cube('umuscl','z-flux',fx(1:ngrid,ilo:ihi,jlo:jhi,kf1:kf2,1:nvar), &
         ihi-ilo+1,jhi-jlo+1,kf2-kf1+1,ngrid,nvar,ndim,ind_grid)
  end if
  ! Save flux in output array
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,3)=fx(l,i,j,k,ivar)*dtdz
           !if (in_debug_region(ind_grid(l)).and.ivar==5) then
           !  print '(a,i8,5i3,1p,g12.4)','flx3',ind_grid(l),l,i,j,k,3,flux(l,i,j,k,ivar,3)
           !end if
        end do
     end do
     do ivar=1,2
        do l=1,ngrid
           tmp (l,i,j,k,ivar,3)=tx(l,i,j,k,ivar)*dtdz
        end do
     end do
  end do
  end do
  end do
#endif

#if NDIM>1
  ! Solve for EMF in Z direction
  CALL cmp_mag_flx(qRT,iu1+1,iu2+1,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qRB,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &           qLT,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qLB,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &               if1  ,if2  ,jf1  ,jf2  ,klo  ,khi  , 2,3,4,6,7,8, &
       &           mark_low_order_edge(:,:,:,:,3),emf,ngrid)
 ! Save vector in output array
  do k=klo,khi
  do j=jf1,jf2
  do i=if1,if2
     do l=1,ngrid
        emfz(l,i,j,k)=emf(l,i,j,k)*dtdx
     end do
  end do
  end do
  end do
#if NDIM==2
  do k=2,2
  do j=jf1,jf2
  do i=if1,if2
     do l=1,ngrid
        emfz(l,i,j,k)=emf(l,i,j,k-1)*dtdx
     end do
  end do
  end do
  end do
#endif
#endif

#if NDIM>2
  ! Solve for EMF in Y direction
  CALL cmp_mag_flx(qRT,iu1+1,iu2+1,ju1  ,ju2  ,ku1+1,ku2+1, &
       &           qLT,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &           qRB,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &           qLB,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &               if1  ,if2  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,8,6,7, &
       &           mark_low_order_edge(:,:,:,:,2),emf,ngrid)
  ! Save vector in output array
  do k=kf1,kf2
  do j=jlo,jhi
  do i=if1,if2
     do l=1,ngrid
        emfy(l,i,j,k)=emf(l,i,j,k)*dtdx
     end do
  end do
  end do
  end do
  ! Solve for EMF in X direction
  CALL cmp_mag_flx(qRT,iu1  ,iu2  ,ju1+1,ju2+1,ku1+1,ku2+1, &
       &           qRB,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &           qLT,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &           qLB,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &               ilo,ihi,jf1  ,jf2  ,kf1  ,kf2  , 3,4,2,7,8,6, &
       &           mark_low_order_edge(:,:,:,:,1),emf,ngrid)
  ! Save vector in output array
  do k=kf1,kf2
  do j=jf1,jf2
  do i=ilo,ihi
     do l=1,ngrid
        emfx(l,i,j,k)=emf(l,i,j,k)*dtdx
     end do
  end do
  end do
  end do
#else
  ! Make sure they are zeroed.
  do k=1,3

  do j=jf1,jf2
  do i=if1,if2
     do l=1,ngrid
        emfy(l,i,j,k)=0.0_dp
     end do
  end do
  end do

  do j=jf1,jf2
  do i=if1,if2
     do l=1,ngrid
        emfx(l,i,j,k)=0.0_dp
     end do
  end do
  end do

  end do
#endif

!  call trace_mpi('mag_unsplit','exit',2)
end subroutine mag_unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==1
SUBROUTINE  trace1d(q,dq,qm,qp,dx,dt,ngrid)
  USE amr_parameters
  USE hydro_parameters
  USE const
  IMPLICIT NONE

  INTEGER ::ngrid
  REAL(dp)::dx,dt

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  ! declare local variables
  INTEGER ::i, j, k, l, n
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  INTEGER ::ir, iu, iv, iw, ip, iA, iB, iC
  REAL(dp)::dtdx
  REAL(dp)::r, u, v, w, p, A, B, C
  REAL(dp)::drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  REAL(dp)::sr0, su0, sv0, sw0, sp0, sA0, sB0, sC0
#ifdef KROME
  REAL(dp)::eff_gamma
#endif

!  call trace_mpi('trace1d','enter',2)

  su0=0.0_dp; sv0=0.0_dp; sw0=0.0_dp
  
  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir = 1; iu = 2; iv = 3 ; iw = 4 ; ip = 5; iA = 6; iB = 7; iC = 8

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi
           DO l = 1, ngrid

              ! Cell centered values
              r = q(l,i,j,k,ir)
              u = q(l,i,j,k,iu)
              v = q(l,i,j,k,iv)
              w = q(l,i,j,k,iw)
              p = q(l,i,j,k,ip)
              A = q(l,i,j,k,iA)
              B = q(l,i,j,k,iB)
              C = q(l,i,j,k,iC)
#ifdef KROME
              eff_gamma = merge(q(l,i,j,k,ichem),gamma,chemistry)
#endif

              ! TVD slopes in X direction
              drx = half*dq(l,i,j,k,ir,1)
              dux = half*dq(l,i,j,k,iu,1)
              dvx = half*dq(l,i,j,k,iv,1)
              dwx = half*dq(l,i,j,k,iw,1)
              dpx = half*dq(l,i,j,k,ip,1)
              dBx = half*dq(l,i,j,k,iB,1)
              dCx = half*dq(l,i,j,k,iC,1)

              ! Source terms (including transverse derivatives)
              sr0 = -u*drx-r*dux
              if(scheme.ne.'induction')then
              su0 = -u*dux-(dpx+B*dBx+C*dCx)/r  
              sv0 = -u*dvx+(A*dBx)/r
              sw0 = -u*dwx+(A*dCx)/r
              endif
#ifdef KROME
              sp0 = -u*dpx-eff_gamma*p*dux
#else
              sp0 = -u*dpx-gamma*p*dux
#endif
              sB0 = -u*dBx+A*dvx-B*dux
              sC0 = -u*dCx+A*dwx-C*dux

              ! Cell-centered predicted states
              r = r + sr0*dtdx
              u = u + su0*dtdx
              v = v + sv0*dtdx
              w = w + sw0*dtdx
              p = p + sp0*dtdx
              B = B + sB0*dtdx
              C = C + sC0*dtdx

              ! Right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              qp(l,i,j,k,iA,1) = A
              qp(l,i,j,k,iB,1) = B - dBx
              qp(l,i,j,k,iC,1) = C - dCx
              qp(l,i,j,k,ir,1) = MAX(smallr, qp(l,i,j,k,ir,1))

              ! Left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              qm(l,i,j,k,iA,1) = A
              qm(l,i,j,k,iB,1) = B + dBx
              qm(l,i,j,k,iC,1) = C + dCx
              qm(l,i,j,k,ir,1) = MAX(smallr, qm(l,i,j,k,ir,1))
           END DO
        END DO
     END DO
  END DO
  
  ! passive scalars
  DO n = 9, nvar
     DO k = klo, khi
        DO j = jlo, jhi
           DO i = ilo, ihi
              DO l = 1, ngrid
                 a   = q(l,i,j,k,n )           ! Cell centered values
                 u   = q(l,i,j,k,iu)  
                 dax = half * dq(l,i,j,k,n,1)  ! TVD slope
                 sa0 = -u*dax                  ! Source terms
                 a   = a + sa0*dtdx            ! Predicted state
                 qp(l,i,j,k,n,1) = a - dax     ! Right state
                 qm(l,i,j,k,n,1) = a + dax     ! Left state
              END DO
           END DO
        END DO
     END DO
  END DO


!  call trace_mpi('trace1d','exit',2)
END SUBROUTINE trace1d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==2
SUBROUTINE trace2d(q,bf,dq,dbf,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dt,ngrid)
  USE amr_parameters
  USE hydro_parameters
  USE const
  IMPLICIT NONE

  INTEGER ::ngrid
  REAL(dp)::dx, dy, dt

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  REAL(dp),DIMENSION(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLB

  ! Declare local variables
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Ez
  INTEGER ::i, j, k, l, n
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  INTEGER ::ir, iu, iv, iw, ip, iA, iB, iC 
  REAL(dp)::dtdx, dtdy, smallp
  REAL(dp)::r, u, v, w, p, A, B, C
  REAL(dp)::ELL, ELR, ERL, ERR
  REAL(dp)::drx, dux, dvx, dwx, dpx, dAx, dBx, dCx
  REAL(dp)::dry, duy, dvy, dwy, dpy, dAy, dBy, dCy
  REAL(dp)::sr0, su0, sv0, sw0, sp0, sA0, sB0, sC0
  REAL(dp)::AL, AR, BL, BR
  REAL(dp)::dALy, dARy, dBLx, dBRx
  REAL(DP)::sAL0, sAR0, sBL0, sBR0
#ifdef KROME
  REAL(DP)::eff_gamma
#endif

!  call trace_mpi('trace2d','enter',2)

  su0=0.0_dp; sv0=0.0_dp; sw0=0.0_dp
  
  dtdx = dt/dx
  dtdy = dt/dy
  smallp = smallr*smallc**2/gamma

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; iw=4; ip=5; ia=6; ib=7; ic=8 

  DO k = klo, ku2
     DO j = jlo, ju2
        DO i = ilo, iu2
           DO l = 1, ngrid
              u = 0.25*(q(l,i-1,j-1,k,iu)+q(l,i-1,j,k,iu)+q(l,i,j-1,k,iu)+q(l,i,j,k,iu))
              v = 0.25*(q(l,i-1,j-1,k,iv)+q(l,i-1,j,k,iv)+q(l,i,j-1,k,iv)+q(l,i,j,k,iv))
              A = 0.5*(bf(l,i,j-1,k,1)+bf(l,i,j,k,1))
              B = 0.5*(bf(l,i-1,j,k,2)+bf(l,i,j,k,2))
              Ez(l,i,j,k)=u*B-v*A
           END DO
        END DO
     END DO
  END DO

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi
           DO l = 1, ngrid

              ! Cell centered values
              r =    q(l,i,j,k,ir)
              u =    q(l,i,j,k,iu)
              v =    q(l,i,j,k,iv)
              w =    q(l,i,j,k,iw)
              p =    q(l,i,j,k,ip)
              A =    q(l,i,j,k,iA)
              B =    q(l,i,j,k,iB)
              C =    q(l,i,j,k,iC)
#ifdef KROME
              eff_gamma = merge(q(l,i,j,k,ichem),gamma,chemistry)
#endif

              ! Face centered variables
              AL =  bf(l,i  ,j  ,k,1)
              AR =  bf(l,i+1,j  ,k,1)
              BL =  bf(l,i  ,j  ,k,2)
              BR =  bf(l,i  ,j+1,k,2)

              ! Cell centered TVD slopes in X direction
              drx = half * dq(l,i,j,k,ir,1)
              dux = half * dq(l,i,j,k,iu,1)
              dvx = half * dq(l,i,j,k,iv,1)
              dwx = half * dq(l,i,j,k,iw,1)
              dpx = half * dq(l,i,j,k,ip,1)
              dBx = half * dq(l,i,j,k,iB,1)
              dCx = half * dq(l,i,j,k,iC,1)

              ! Cell centered TVD slopes in Y direction
              dry = half * dq(l,i,j,k,ir,2)
              duy = half * dq(l,i,j,k,iu,2)
              dvy = half * dq(l,i,j,k,iv,2)
              dwy = half * dq(l,i,j,k,iw,2)
              dpy = half * dq(l,i,j,k,ip,2)
              dAy = half * dq(l,i,j,k,iA,2)
              dCy = half * dq(l,i,j,k,iC,2)

              ! Face centered TVD slopes in transverse direction
              dALy = half * dbf(l,i  ,j  ,k,1,1)
              dARy = half * dbf(l,i+1,j  ,k,1,1)
              dBLx = half * dbf(l,i  ,j  ,k,2,1)
              dBRx = half * dbf(l,i  ,j+1,k,2,1)

              ! Edge centered electric field Ez = uB-vA
              ELL = Ez(l,i  ,j  ,k)
              ELR = Ez(l,i  ,j+1,k)
              ERL = Ez(l,i+1,j  ,k)
              ERR = Ez(l,i+1,j+1,k)

              ! Face-centered predicted states
              sAL0 = +(ELR-ELL)*dtdy*half
              sAR0 = +(ERR-ERL)*dtdy*half
              sBL0 = -(ERL-ELL)*dtdx*half
              sBR0 = -(ERR-ELR)*dtdx*half

              AL = AL + sAL0
              AR = AR + sAR0
              BL = BL + sBL0
              BR = BR + sBR0
              
              ! Source terms (including transverse derivatives)
              sr0 = (-u*drx-dux*r)*dtdx + (-v*dry-dvy*r)*dtdy
              if(scheme.ne.'induction')then
              su0 = (-u*dux-(dpx+B*dBx+C*dCx)/r)*dtdx + (-v*duy+B*dAy/r)*dtdy 
              sv0 = (-u*dvx+A*dBx/r)*dtdx + (-v*dvy-(dpy+A*dAy+C*dCy)/r)*dtdy
              sw0 = (-u*dwx+A*dCx/r)*dtdx + (-v*dwy+B*dCy/r)*dtdy
              endif
#ifdef KROME
              sp0 = (-u*dpx-dux*eff_gamma*p)*dtdx + (-v*dpy-dvy*eff_gamma*p)*dtdy
#else
              sp0 = (-u*dpx-dux*gamma*p)*dtdx + (-v*dpy-dvy*gamma*p)*dtdy
#endif
              sC0 = (-u*dCx-C*dux+A*dwx)*dtdx + (-v*dCy-C*dvy+B*dwy)*dtdy

              ! Cell-centered predicted states
              r = r + sr0
              r = max(r,smallr)
              u = u + su0
              v = v + sv0
              w = w + sw0
              p = p + sp0
              p = max(p,smallp)
              C = C + sC0
              A = 0.5*(AL+AR)
              B = 0.5*(BL+BR)

              if (debug) then
                 if (count(r==smallr)>0) print*,'MK: ',myid,count(r==smallr)
                 if (count(p==smallp)>0) print*,'MK: ',myid,count(p==smallp)
              end if

              ! Face averaged right state at left interface
              qp(l,i,j,k,ir,1) = r - drx
              qp(l,i,j,k,iu,1) = u - dux
              qp(l,i,j,k,iv,1) = v - dvx
              qp(l,i,j,k,iw,1) = w - dwx
              qp(l,i,j,k,ip,1) = p - dpx
              qp(l,i,j,k,iA,1) = AL     
              qp(l,i,j,k,iB,1) = B - dBx
              qp(l,i,j,k,iC,1) = C - dCx
              qp(l,i,j,k,ir,1) = MAX(smallr, qp(l,i,j,k,ir,1))
              qp(l,i,j,k,ip,1) = MAX(smallp, qp(l,i,j,k,ip,1))

              ! Face averaged left state at right interface
              qm(l,i,j,k,ir,1) = r + drx
              qm(l,i,j,k,iu,1) = u + dux
              qm(l,i,j,k,iv,1) = v + dvx
              qm(l,i,j,k,iw,1) = w + dwx
              qm(l,i,j,k,ip,1) = p + dpx
              qm(l,i,j,k,iA,1) = AR     
              qm(l,i,j,k,iB,1) = B + dBx
              qm(l,i,j,k,iC,1) = C + dCx
              qm(l,i,j,k,ir,1) = MAX(smallr, qm(l,i,j,k,ir,1))
              qm(l,i,j,k,ip,1) = MAX(smallp, qm(l,i,j,k,ip,1))

              ! Face averaged top state at bottom interface
              qp(l,i,j,k,ir,2) = r - dry
              qp(l,i,j,k,iu,2) = u - duy
              qp(l,i,j,k,iv,2) = v - dvy
              qp(l,i,j,k,iw,2) = w - dwy
              qp(l,i,j,k,ip,2) = p - dpy
              qp(l,i,j,k,iA,2) = A - dAy
              qp(l,i,j,k,iB,2) = BL     
              qp(l,i,j,k,iC,2) = C - dCy
              qp(l,i,j,k,ir,2) = MAX(smallr, qp(l,i,j,k,ir,2))
              qp(l,i,j,k,ip,2) = MAX(smallp, qp(l,i,j,k,ip,2))

              ! Face averaged bottom state at top interface
              qm(l,i,j,k,ir,2) = r + dry
              qm(l,i,j,k,iu,2) = u + duy
              qm(l,i,j,k,iv,2) = v + dvy
              qm(l,i,j,k,iw,2) = w + dwy
              qm(l,i,j,k,ip,2) = p + dpy
              qm(l,i,j,k,iA,2) = A + dAy
              qm(l,i,j,k,iB,2) = BR     
              qm(l,i,j,k,iC,2) = C + dCy
              qm(l,i,j,k,ir,2) = MAX(smallr, qm(l,i,j,k,ir,2))
              qm(l,i,j,k,ip,2) = MAX(smallp, qm(l,i,j,k,ip,2))

              ! Edge averaged right-top corner state (RT->LL)
              qRT(l,i,j,k,ir,3) = r + (+drx+dry)
              qRT(l,i,j,k,iu,3) = u + (+dux+duy)
              qRT(l,i,j,k,iv,3) = v + (+dvx+dvy)
              qRT(l,i,j,k,iw,3) = w + (+dwx+dwy)
              qRT(l,i,j,k,ip,3) = p + (+dpx+dpy)
              qRT(l,i,j,k,iC,3) = C + (+dCx+dCy)
              qRT(l,i,j,k,iA,3) = AR+ (   +dARy)
              qRT(l,i,j,k,iB,3) = BR+ (+dBRx   )
              qRT(l,i,j,k,ir,3) = MAX(smallr, qRT(l,i,j,k,ir,3))
              qRT(l,i,j,k,ip,3) = MAX(smallp, qRT(l,i,j,k,ip,3))

              ! Edge averaged right-bottom corner state (RB->LR)
              qRB(l,i,j,k,ir,3) = r + (+drx-dry)
              qRB(l,i,j,k,iu,3) = u + (+dux-duy)
              qRB(l,i,j,k,iv,3) = v + (+dvx-dvy)
              qRB(l,i,j,k,iw,3) = w + (+dwx-dwy)
              qRB(l,i,j,k,ip,3) = p + (+dpx-dpy)
              qRB(l,i,j,k,iC,3) = C + (+dCx-dCy)
              qRB(l,i,j,k,iA,3) = AR+ (   -dARy)
              qRB(l,i,j,k,iB,3) = BL+ (+dBLx   )
              qRB(l,i,j,k,ir,3) = MAX(smallr, qRB(l,i,j,k,ir,3))
              qRB(l,i,j,k,ip,3) = MAX(smallp, qRB(l,i,j,k,ip,3))

              ! Edge averaged left-top corner state (LT->RL)
              qLT(l,i,j,k,ir,3) = r + (-drx+dry)
              qLT(l,i,j,k,iu,3) = u + (-dux+duy)
              qLT(l,i,j,k,iv,3) = v + (-dvx+dvy)
              qLT(l,i,j,k,iw,3) = w + (-dwx+dwy)
              qLT(l,i,j,k,ip,3) = p + (-dpx+dpy)
              qLT(l,i,j,k,iC,3) = C + (-dCx+dCy)
              qLT(l,i,j,k,iA,3) = AL+ (   +dALy)
              qLT(l,i,j,k,iB,3) = BR+ (-dBRx   )
              qLT(l,i,j,k,ir,3) = MAX(smallr, qLT(l,i,j,k,ir,3))
              qLT(l,i,j,k,ip,3) = MAX(smallp, qLT(l,i,j,k,ip,3))

              ! Edge averaged left-bottom corner state (LB->RR)
              qLB(l,i,j,k,ir,3) = r + (-drx-dry)
              qLB(l,i,j,k,iu,3) = u + (-dux-duy)
              qLB(l,i,j,k,iv,3) = v + (-dvx-dvy)
              qLB(l,i,j,k,iw,3) = w + (-dwx-dwy)
              qLB(l,i,j,k,ip,3) = p + (-dpx-dpy)
              qLB(l,i,j,k,iC,3) = C + (-dCx-dCy)
              qLB(l,i,j,k,iA,3) = AL+ (   -dALy)
              qLB(l,i,j,k,iB,3) = BL+ (-dBLx   )
              qLB(l,i,j,k,ir,3) = MAX(smallr, qLB(l,i,j,k,ir,3))
              qLB(l,i,j,k,ip,3) = MAX(smallp, qLB(l,i,j,k,ip,3))

           END DO
        END DO
     END DO
  END DO

  ! Passive scalars
  DO n = 9, nvar
     DO k = klo, khi
        DO j = jlo, jhi
           DO i = ilo, ihi
              DO l = 1, ngrid
                 r   = q(l,i,j,k,n )              ! Cell centered values
                 u   = q(l,i,j,k,iu)
                 v   = q(l,i,j,k,iv)
                 drx = half * dq(l,i,j,k,n,1)     ! TVD slopes
                 dry = half * dq(l,i,j,k,n,2)
                 sr0 = -u*drx*dtdx -v*dry*dtdy    ! Source terms
                 r   = r + sr0                    ! Predicted state
                 qp(l,i,j,k,n,1) = r - drx        ! Right state
                 qm(l,i,j,k,n,1) = r + drx        ! Left state
                 qp(l,i,j,k,n,2) = r - dry        ! Top state
                 qm(l,i,j,k,n,2) = r + dry        ! Bottom state
              END DO
           END DO
        END DO
     END DO
  END DO

!  call trace_mpi('trace2d','exit',2)
END SUBROUTINE trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM==3
SUBROUTINE trace3d(q,bf,dq,dq2,dbf,dS,qm,qp,qRT,qRB,qLT,qLB,dx,dy,dz,dt,ngrid,ind_grid)
  USE amr_parameters
  USE amr_commons, only: dbg_grid, dbg_ind, myid, xg, local_debug
  USE hydro_parameters
  USE const
  IMPLICIT NONE

  INTEGER ::ngrid
  INTEGER ::ind_grid(ngrid), ig
  REAL(dp)::dx, dy, dz, dt

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+1)::q  
  REAL(dp),DIMENSION(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq,dq2
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::dS
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm 
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp 

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:3)::qLB

  ! Declare local variables
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Ex
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Ey
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::Ez

  INTEGER ::i, j, k, l, n
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  INTEGER, parameter ::ir=1, iu=2, iv=3, iw=4, ip=5, ia=6, ib=7, ic=8 
  REAL(dp)::dtdx, dtdy, dtdz, smallp
  REAL(dp), dimension(nvector) :: invr, r, u, v, w, p, A, B, C
  REAL(dp), dimension(nvector) :: S, dSx, dSy, dSz, sS0
  REAL(dp), dimension(nvector) :: ELL, ELR, ERL, ERR
  REAL(dp), dimension(nvector) :: FLL, FLR, FRL, FRR
  REAL(dp), dimension(nvector) :: GLL, GLR, GRL, GRR
  REAL(dp), dimension(nvector) :: drx, drx2, dux, dvx, dwx, dpx, dpx2, dAx, dBx, dCx
  REAL(dp), dimension(nvector) :: dry, dry2, duy, dvy, dwy, dpy, dpy2, dAy, dbY, dCy
  REAL(dp), dimension(nvector) :: drz, drz2, duz, dvz, dwz, dpz, dpz2, dAz, dBz, dCz
  REAL(dp), dimension(nvector) :: t1, t2, t3, t4, t5, t6
  REAL(dp), dimension(nvector) :: sr0, su0, sv0, sw0, sp0, sA0, sB0, sC0
  REAL(dp), dimension(nvector) :: AL, AR, BL, BR, CL, CR
  REAL(dp), dimension(nvector) :: dALy, dARy, dALz, dARz
  REAL(dp), dimension(nvector) :: dBLx, dBRx, dBLz, dBRz
  REAL(dp), dimension(nvector) :: dCLx, dCRx, dCLy, dCRy
  REAL(dp), dimension(nvector) :: sAL0, sAR0, sBL0, sBR0, sCL0, sCR0
  REAL(dp), DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2):: rho_t, p_t
  REAL(dp), dimension(nvector) :: eff_gamma
#ifdef KROME
  INTEGER, save ::nprintX=100
  INTEGER ::idir,
  REAL(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,6) :: totmass ! Total mass in passive scalars
  REAL(dp), dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,ichem+1:nvar):: chem_t
#endif

  LOGICAL ::noinduc 
  LOGICAL ::in_debug_region
  INTEGER ::iS
  
!  call trace_mpi('trace3d','enter',2)

  iS = nvar+1
  su0=0.0_dp; sv0=0.0_dp; sw0=0.0_dp

  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  smallp = smallr*smallc**2/gamma
  noinduc = scheme.ne.'induction'

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  DO k = klo, ku2
     DO j = jlo, ju2
        DO i = ilo, iu2
              v = 0.25*(q(:,i,j-1,k-1,iv)+q(:,i,j-1,k,iv)+q(:,i,j,k-1,iv)+q(:,i,j,k,iv))
              w = 0.25*(q(:,i,j-1,k-1,iw)+q(:,i,j-1,k,iw)+q(:,i,j,k-1,iw)+q(:,i,j,k,iw))
              B = 0.5*(bf(:,i,j,k-1,2)+bf(:,i,j,k,2))
              C = 0.5*(bf(:,i,j-1,k,3)+bf(:,i,j,k,3))
              Ex(:,i,j,k)=v*C-w*B

              u = 0.25*(q(:,i-1,j,k-1,iu)+q(:,i-1,j,k,iu)+q(:,i,j,k-1,iu)+q(:,i,j,k,iu))
              w = 0.25*(q(:,i-1,j,k-1,iw)+q(:,i-1,j,k,iw)+q(:,i,j,k-1,iw)+q(:,i,j,k,iw))
              A = 0.5*(bf(:,i,j,k-1,1)+bf(:,i,j,k,1))
              C = 0.5*(bf(:,i-1,j,k,3)+bf(:,i,j,k,3))
              Ey(:,i,j,k)=w*A-u*C

              u = 0.25*(q(:,i-1,j-1,k,iu)+q(:,i-1,j,k,iu)+q(:,i,j-1,k,iu)+q(:,i,j,k,iu))
              v = 0.25*(q(:,i-1,j-1,k,iv)+q(:,i-1,j,k,iv)+q(:,i,j-1,k,iv)+q(:,i,j,k,iv))
              A = 0.5*(bf(:,i,j-1,k,1)+bf(:,i,j,k,1))
              B = 0.5*(bf(:,i-1,j,k,2)+bf(:,i,j,k,2))
              Ez(:,i,j,k)=u*B-v*A
        END DO
     END DO
  END DO

  ! fill positive definite variables with reasonable values in "ghost-zones"
  rho_t(:,:,:,:) = q(:,:,:,:,ir)
  p_t(:,:,:,:)   = q(:,:,:,:,ip)
#ifdef KROME
  chem_t(:,:,:,:,:) = q(:,:,:,:,ichem+1:nvar)
#endif

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi

              ! Cell centered values
              r =    q(:,i,j,k,ir)
              u =    q(:,i,j,k,iu)
              v =    q(:,i,j,k,iv)
              w =    q(:,i,j,k,iw)            
              p =    q(:,i,j,k,ip)
              A =    q(:,i,j,k,iA)
              B =    q(:,i,j,k,iB)
              C =    q(:,i,j,k,iC)            
              S =    q(:,i,j,k,iS)

#ifdef KROME
              eff_gamma = merge(q(:,i,j,k,ichem),gamma,chemistry)
#else
              eff_gamma = gamma
#endif

              ! Face centered variables
              AL =  bf(:,i  ,j  ,k  ,1)
              AR =  bf(:,i+1,j  ,k  ,1)
              BL =  bf(:,i  ,j  ,k  ,2)
              BR =  bf(:,i  ,j+1,k  ,2)
              CL =  bf(:,i  ,j  ,k  ,3)
              CR =  bf(:,i  ,j  ,k+1,3)

              ! Cell centered TVD slopes in X, Y and Z directions
              drx = half * dq(:,i,j,k,ir,1)
              drx2= half *dq2(:,i,j,k,ir,1)
              dux = half * dq(:,i,j,k,iu,1)
              dvx = half * dq(:,i,j,k,iv,1)
              dwx = half * dq(:,i,j,k,iw,1)
              dpx = half * dq(:,i,j,k,ip,1)
              dpx2= half *dq2(:,i,j,k,ip,1)
              dBx = half * dq(:,i,j,k,iB,1)
              dCx = half * dq(:,i,j,k,iC,1)
              dSx = half * dS(:,i,j,k,1)

              dry = half * dq(:,i,j,k,ir,2)
              dry2= half *dq2(:,i,j,k,ir,2)
              duy = half * dq(:,i,j,k,iu,2)
              dvy = half * dq(:,i,j,k,iv,2)
              dwy = half * dq(:,i,j,k,iw,2)
              dpy = half * dq(:,i,j,k,ip,2)
              dpy2= half *dq2(:,i,j,k,ip,2)
              duy = half * dq(:,i,j,k,iu,2)
              dAy = half * dq(:,i,j,k,iA,2)
              dCy = half * dq(:,i,j,k,iC,2)
              dSy = half * dS(:,i,j,k,2)

              drz = half * dq(:,i,j,k,ir,3)
              drz2= half *dq2(:,i,j,k,ir,3)
              duz = half * dq(:,i,j,k,iu,3)
              dvz = half * dq(:,i,j,k,iv,3)
              dwz = half * dq(:,i,j,k,iw,3)
              dpz = half * dq(:,i,j,k,ip,3)
              dpz2= half *dq2(:,i,j,k,ip,3)
              dAz = half * dq(:,i,j,k,iA,3)
              dBz = half * dq(:,i,j,k,iB,3)
              dSz = half * dS(:,i,j,k,3)

              ! Face centered TVD slopes in transverse directions
              dALy = half * dbf(:,i  ,j  ,k  ,1,1)
              dARy = half * dbf(:,i+1,j  ,k  ,1,1)
              dALz = half * dbf(:,i  ,j  ,k  ,1,2)
              dARz = half * dbf(:,i+1,j  ,k  ,1,2)

              dBLx = half * dbf(:,i  ,j  ,k  ,2,1)
              dBRx = half * dbf(:,i  ,j+1,k  ,2,1)
              dBLz = half * dbf(:,i  ,j  ,k  ,2,2)
              dBRz = half * dbf(:,i  ,j+1,k  ,2,2)

              dCLx = half * dbf(:,i  ,j  ,k  ,3,1)
              dCRx = half * dbf(:,i  ,j  ,k+1,3,1)
              dCLy = half * dbf(:,i  ,j  ,k  ,3,2)
              dCRy = half * dbf(:,i  ,j  ,k+1,3,2)

              ! Edge centered electric field in X, Y and Z directions
              ELL = Ex(:,i,j  ,k  )
              ELR = Ex(:,i,j  ,k+1)
              ERL = Ex(:,i,j+1,k  )
              ERR = Ex(:,i,j+1,k+1)

              FLL = Ey(:,i  ,j,k  )
              FLR = Ey(:,i  ,j,k+1)
              FRL = Ey(:,i+1,j,k  )
              FRR = Ey(:,i+1,j,k+1)

              GLL = Ez(:,i  ,j  ,k)
              GLR = Ez(:,i  ,j+1,k)
              GRL = Ez(:,i+1,j  ,k)
              GRR = Ez(:,i+1,j+1,k)

              ! Face-centered predicted states
              sAL0 = +(GLR-GLL)*dtdy*half -(FLR-FLL)*dtdz*half
              sAR0 = +(GRR-GRL)*dtdy*half -(FRR-FRL)*dtdz*half
              sBL0 = -(GRL-GLL)*dtdx*half +(ELR-ELL)*dtdz*half
              sBR0 = -(GRR-GLR)*dtdx*half +(ERR-ERL)*dtdz*half
              sCL0 = +(FRL-FLL)*dtdx*half -(ERL-ELL)*dtdy*half
              sCR0 = +(FRR-FLR)*dtdx*half -(ERR-ELR)*dtdy*half
              
              AL = AL + sAL0
              AR = AR + sAR0
              BL = BL + sBL0
              BR = BR + sBR0
              CL = CL + sCL0
              CR = CR + sCR0

              ! Source terms (including transverse derivatives)
              invr = one / r
              sr0 = (-u*drx-dux*r)*dtdx + (-v*dry-dvy*r)*dtdy + (-w*drz-dwz*r)*dtdz 
              su0 = merge((-u*dux-(dpx+B*dBx+C*dCx)*invr)*dtdx + (-v*duy+B*dAy*invr)*dtdy + (-w*duz+C*dAz*invr)*dtdz,zero,noinduc)
              sv0 = merge((-u*dvx+A*dBx*invr)*dtdx + (-v*dvy-(dpy+A*dAy+C*dCy)*invr)*dtdy + (-w*dvz+C*dBz*invr)*dtdz,zero,noinduc)
              sw0 = merge((-u*dwx+A*dCx*invr)*dtdx + (-v*dwy+B*dCy*invr)*dtdy + (-w*dwz-(dpz+A*dAz+B*dBz)*invr)*dtdz,zero,noinduc)

              ! Evolve entropy instead of pressure
              sp0 = (-u*dpx-dux*eff_gamma*p)*dtdx + (-v*dpy-dvy*eff_gamma*p)*dtdy + (-w*dpz-dwz*eff_gamma*p)*dtdz
              sS0 = (-u*dSx)*dtdx + (-v*dSy)*dtdy + (-w*dSz)*dtdz

              ! Save predicted cell-centered pseudo-entropy; we are borrowing space in q(:,:,:,:,ip),
              ! which is not used below this point.  In the calling routine this is again copied over
              ! into uin, and is thus made available to godfine1, through that array.
              if (entropy_fix) then
                 if ((i==1.or.i==2).and.(j==1.or.j==2).and.(k==1.or.k==2)) then
                    q(:,i,j,k,ip) = S + 2.0*sS0
                    if (local_debug) then
                       do ig=1,ngrid
                          if (in_debug_region(ind_grid(ig))) then
                             print'(3hMKE,i8,i3,3f10.6,1p,12e11.3)', &
                               ind_grid(ig),(i-1)+(j-1)*2+(k-1)*4, &
                               xg(ind_grid(ig),1)+(i-1.5)/32.0, &
                               xg(ind_grid(ig),2)+(j-1.5)/32.0, &
                               xg(ind_grid(ig),3)+(k-1.5)/32.0, &
                               log(p(ig)-dpx(ig)) - gamma*log(r(ig)-drx(ig)), &
                               log(p(ig))-gamma*log(r(ig)), &
                               log(p(ig)+dpx(ig)) - gamma*log(r(ig)+drx(ig)), &
                               q(ig,i,j,k,ip), &
                               dSx(ig),dSy(ig),dSz(ig), &
                               u(ig),v(ig),w(ig),dtdx
                          end if
                       end do
                    end if
                 end if
              end if

              ! Cell-centered predicted states
              r = r + sr0
              r = max(r,smallr)
              rho_t(:,i,j,k) = r
              u = u + su0
              v = v + sv0
              w = w + sw0
              p = p + sp0
              p = max(p,smallp)
              p_t(:,i,j,k) = p
              A = half*(AL+AR)
              B = half*(BL+BR)
              C = half*(CL+CR)


              ! Face averaged right state at left interface
              qp(:,i,j,k,ir,1) = r - drx
              qp(:,i,j,k,iu,1) = u - dux
              qp(:,i,j,k,iv,1) = v - dvx
              qp(:,i,j,k,iw,1) = w - dwx
              qp(:,i,j,k,ip,1) = p - dpx
              qp(:,i,j,k,iA,1) = AL
              qp(:,i,j,k,iB,1) = B - dBx
              qp(:,i,j,k,iC,1) = C - dCx

              ! Face averaged left state at right interface
              qm(:,i,j,k,ir,1) = r + drx
              qm(:,i,j,k,iu,1) = u + dux
              qm(:,i,j,k,iv,1) = v + dvx
              qm(:,i,j,k,iw,1) = w + dwx
              qm(:,i,j,k,ip,1) = p + dpx
              qm(:,i,j,k,iA,1) = AR
              qm(:,i,j,k,iB,1) = B + dBx
              qm(:,i,j,k,iC,1) = C + dCx

              ! Face averaged top state at bottom interface
              qp(:,i,j,k,ir,2) = r - dry
              qp(:,i,j,k,iu,2) = u - duy
              qp(:,i,j,k,iv,2) = v - dvy
              qp(:,i,j,k,iw,2) = w - dwy
              qp(:,i,j,k,ip,2) = p - dpy
              qp(:,i,j,k,iA,2) = A - dAy
              qp(:,i,j,k,iB,2) = BL
              qp(:,i,j,k,iC,2) = C - dCy

              ! Face averaged bottom state at top interface
              qm(:,i,j,k,ir,2) = r + dry
              qm(:,i,j,k,iu,2) = u + duy
              qm(:,i,j,k,iv,2) = v + dvy
              qm(:,i,j,k,iw,2) = w + dwy
              qm(:,i,j,k,ip,2) = p + dpy
              qm(:,i,j,k,iA,2) = A + dAy
              qm(:,i,j,k,iB,2) = BR
              qm(:,i,j,k,iC,2) = C + dCy

              ! Face averaged front state at back interface
              qp(:,i,j,k,ir,3) = r - drz
              qp(:,i,j,k,iu,3) = u - duz
              qp(:,i,j,k,iv,3) = v - dvz
              qp(:,i,j,k,iw,3) = w - dwz
              qp(:,i,j,k,ip,3) = p - dpz
              qp(:,i,j,k,iA,3) = A - dAz
              qp(:,i,j,k,iB,3) = B - dBz
              qp(:,i,j,k,iC,3) = CL

              ! Face averaged back state at front interface
              qm(:,i,j,k,ir,3) = r + drz
              qm(:,i,j,k,iu,3) = u + duz
              qm(:,i,j,k,iv,3) = v + dvz
              qm(:,i,j,k,iw,3) = w + dwz
              qm(:,i,j,k,ip,3) = p + dpz
              qm(:,i,j,k,iA,3) = A + dAz
              qm(:,i,j,k,iB,3) = B + dBz
              qm(:,i,j,k,iC,3) = CR

              ! X-edge averaged right-top corner state (RT->LL)
              qRT(:,i,j,k,ir,1) = r + (+dry2+drz2)
              qRT(:,i,j,k,iu,1) = u + (+duy+duz)
              qRT(:,i,j,k,iv,1) = v + (+dvy+dvz)
              qRT(:,i,j,k,iw,1) = w + (+dwy+dwz)
              qRT(:,i,j,k,ip,1) = p + (+dpy2+dpz2)
              qRT(:,i,j,k,iA,1) = A + (+dAy+dAz)
              qRT(:,i,j,k,iB,1) = BR+ (   +dBRz)
              qRT(:,i,j,k,iC,1) = CR+ (+dCRy   )

              ! X-edge averaged right-bottom corner state (RB->LR)
              qRB(:,i,j,k,ir,1) = r + (+dry2-drz2)
              qRB(:,i,j,k,iu,1) = u + (+duy-duz)
              qRB(:,i,j,k,iv,1) = v + (+dvy-dvz)
              qRB(:,i,j,k,iw,1) = w + (+dwy-dwz)
              qRB(:,i,j,k,ip,1) = p + (+dpy2-dpz2)
              qRB(:,i,j,k,iA,1) = A + (+dAy-dAz)
              qRB(:,i,j,k,iB,1) = BR+ (   -dBRz)
              qRB(:,i,j,k,iC,1) = CL+ (+dCLy   )

              ! X-edge averaged left-top corner state (LT->RL)
              qLT(:,i,j,k,ir,1) = r + (-dry2+drz2)
              qLT(:,i,j,k,iu,1) = u + (-duy+duz)
              qLT(:,i,j,k,iv,1) = v + (-dvy+dvz)
              qLT(:,i,j,k,iw,1) = w + (-dwy+dwz)
              qLT(:,i,j,k,ip,1) = p + (-dpy2+dpz2)
              qLT(:,i,j,k,iA,1) = A + (-dAy+dAz)
              qLT(:,i,j,k,iB,1) = BL+ (   +dBLz)
              qLT(:,i,j,k,iC,1) = CR+ (-dCRy   )

              ! X-edge averaged left-bottom corner state (LB->RR)
              qLB(:,i,j,k,ir,1) = r + (-dry2-drz2)
              qLB(:,i,j,k,iu,1) = u + (-duy-duz)
              qLB(:,i,j,k,iv,1) = v + (-dvy-dvz)
              qLB(:,i,j,k,iw,1) = w + (-dwy-dwz)
              qLB(:,i,j,k,ip,1) = p + (-dpy2-dpz2)
              qLB(:,i,j,k,iA,1) = A + (-dAy-dAz)
              qLB(:,i,j,k,iB,1) = BL+ (   -dBLz)
              qLB(:,i,j,k,iC,1) = CL+ (-dCLy   )

              ! Y-edge averaged right-top corner state (RT->LL)
              qRT(:,i,j,k,ir,2) = r + (+drx2+drz2)
              qRT(:,i,j,k,iu,2) = u + (+dux+duz)
              qRT(:,i,j,k,iv,2) = v + (+dvx+dvz)
              qRT(:,i,j,k,iw,2) = w + (+dwx+dwz)
              qRT(:,i,j,k,ip,2) = p + (+dpx2+dpz2)
              qRT(:,i,j,k,iA,2) = AR+ (   +dARz)
              qRT(:,i,j,k,iB,2) = B + (+dBx+dBz)
              qRT(:,i,j,k,iC,2) = CR+ (+dCRx   )

              ! Y-edge averaged right-bottom corner state (RB->LR)
              qRB(:,i,j,k,ir,2) = r + (+drx2-drz2)
              qRB(:,i,j,k,iu,2) = u + (+dux-duz)
              qRB(:,i,j,k,iv,2) = v + (+dvx-dvz)
              qRB(:,i,j,k,iw,2) = w + (+dwx-dwz)
              qRB(:,i,j,k,ip,2) = p + (+dpx2-dpz2)
              qRB(:,i,j,k,iA,2) = AR+ (   -dARz)
              qRB(:,i,j,k,iB,2) = B + (+dBx-dBz)
              qRB(:,i,j,k,iC,2) = CL+ (+dCLx   )

              ! Y-edge averaged left-top corner state (LT->RL)
              qLT(:,i,j,k,ir,2) = r + (-drx2+drz2)
              qLT(:,i,j,k,iu,2) = u + (-dux+duz)
              qLT(:,i,j,k,iv,2) = v + (-dvx+dvz)
              qLT(:,i,j,k,iw,2) = w + (-dwx+dwz)
              qLT(:,i,j,k,ip,2) = p + (-dpx2+dpz2)
              qLT(:,i,j,k,iA,2) = AL+ (   +dALz)
              qLT(:,i,j,k,iB,2) = B + (-dBx+dBz)
              qLT(:,i,j,k,iC,2) = CR+ (-dCRx   )

              ! Y-edge averaged left-bottom corner state (LB->RR)
              qLB(:,i,j,k,ir,2) = r + (-drx2-drz2)
              qLB(:,i,j,k,iu,2) = u + (-dux-duz)
              qLB(:,i,j,k,iv,2) = v + (-dvx-dvz)
              qLB(:,i,j,k,iw,2) = w + (-dwx-dwz)
              qLB(:,i,j,k,ip,2) = p + (-dpx2-dpz2)
              qLB(:,i,j,k,iA,2) = AL+ (   -dALz)
              qLB(:,i,j,k,iB,2) = B + (-dBx-dBz)
              qLB(:,i,j,k,iC,2) = CL+ (-dCLx   )

              ! Z-edge averaged right-top corner state (RT->LL)
              qRT(:,i,j,k,ir,3) = r + (+drx2+dry2)
              qRT(:,i,j,k,iu,3) = u + (+dux+duy)
              qRT(:,i,j,k,iv,3) = v + (+dvx+dvy)
              qRT(:,i,j,k,iw,3) = w + (+dwx+dwy)
              qRT(:,i,j,k,ip,3) = p + (+dpx2+dpy2)
              qRT(:,i,j,k,iA,3) = AR+ (   +dARy)
              qRT(:,i,j,k,iB,3) = BR+ (+dBRx   )
              qRT(:,i,j,k,iC,3) = C + (+dCx+dCy)

              ! Z-edge averaged right-bottom corner state (RB->LR)
              qRB(:,i,j,k,ir,3) = r + (+drx2-dry2)
              qRB(:,i,j,k,iu,3) = u + (+dux-duy)
              qRB(:,i,j,k,iv,3) = v + (+dvx-dvy)
              qRB(:,i,j,k,iw,3) = w + (+dwx-dwy)
              qRB(:,i,j,k,ip,3) = p + (+dpx2-dpy2)
              qRB(:,i,j,k,iA,3) = AR+ (   -dARy)
              qRB(:,i,j,k,iB,3) = BL+ (+dBLx   )
              qRB(:,i,j,k,iC,3) = C + (+dCx-dCy)

              ! Z-edge averaged left-top corner state (LT->RL)
              qLT(:,i,j,k,ir,3) = r + (-drx2+dry2)
              qLT(:,i,j,k,iu,3) = u + (-dux+duy)
              qLT(:,i,j,k,iv,3) = v + (-dvx+dvy)
              qLT(:,i,j,k,iw,3) = w + (-dwx+dwy)
              qLT(:,i,j,k,ip,3) = p + (-dpx2+dpy2)
              qLT(:,i,j,k,iA,3) = AL+ (   +dALy)
              qLT(:,i,j,k,iB,3) = BR+ (-dBRx   )
              qLT(:,i,j,k,iC,3) = C + (-dCx+dCy)

              ! Z-edge averaged left-bottom corner state (LB->RR)
              qLB(:,i,j,k,ir,3) = r + (-drx2-dry2)
              qLB(:,i,j,k,iu,3) = u + (-dux-duy)
              qLB(:,i,j,k,iv,3) = v + (-dvx-dvy)
              qLB(:,i,j,k,iw,3) = w + (-dwx-dwy)
              qLB(:,i,j,k,ip,3) = p + (-dpx2-dpy2)
              qLB(:,i,j,k,iA,3) = AL+ (   -dALy)
              qLB(:,i,j,k,iB,3) = BL+ (-dBLx   )
              qLB(:,i,j,k,iC,3) = C + (-dCx-dCy)

              ! Check that density and pressure are positive definite at cell interface
              call check_negative(qm,i,j,k,drx,dpx,dry,dpy,drz,dpz,'qm',ngrid)
              t1 = -drx; t2 = -dpx; t3 = -dry; t4 = -dpy; t5 = -drz; t6 = -dpz
              call check_negative(qp,i,j,k,t1,t2,t3,t4,t5,t6,'qp',ngrid)
        END DO
     END DO
  END DO

  t2 = 0.; t4 = 0.; t6 = 0.
  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi
              r = rho_t(:,i,j,k)
              p = p_t(:,i,j,k)
              if (debug) then
                 do l=1,ngrid
                    if (rho_t(l,i,j,k) <= 0.) write(*,'(a,1p,3e11.4)') &
                     'density <= 0.:', rho_t(l,i,j,k), q(l,i,j,k,ir), rho_t(l,i,j,k)-q(l,i,j,k,ir)
                    if (p_t(l,i,j,k) <= 0.) write(*,'(a,1p,3e11.4)') &
                     'pressure <= 0.:', p_t(l,i,j,k), q(l,i,j,k,ip), p_t(l,i,j,k)-q(l,i,j,k,ip)
                 enddo
              endif

              ! X-edge averaged right-top corner state (RT->LL)
              !qRT(:,i,j,k,ir,1) = r + (+dry+drz)
              drx = sqrt(sqrt(r*rho_t(:,i,j+1,k)*rho_t(:,i,j,k+1)*rho_t(:,i,j+1,k+1))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i,j+1,k)*  p_t(:,i,j,k+1)*  p_t(:,i,j+1,k+1))) - p
              dry = qRT(:,i,j,k,ir,1) - r
              dpy = qRT(:,i,j,k,ip,1) - p
              qRT(:,i,j,k,ir,1) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qRT(:,i,j,k,ip,1) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! X-edge averaged right-bottom corner state (RB->LR)
              !qRB(:,i,j,k,ir,1) = r + (+dry-drz)
              drx = sqrt(sqrt(r*rho_t(:,i,j+1,k)*rho_t(:,i,j,k-1)*rho_t(:,i,j+1,k-1))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i,j+1,k)*  p_t(:,i,j,k-1)*  p_t(:,i,j+1,k-1))) - p
              dry = qRB(:,i,j,k,ir,1) - r
              dpy = qRB(:,i,j,k,ip,1) - p
              qRB(:,i,j,k,ir,1) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qRB(:,i,j,k,ip,1) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! X-edge averaged left-top corner state (LT->RL)
              !qLT(:,i,j,k,ir,1) = r + (-dry+drz)
              drx = sqrt(sqrt(r*rho_t(:,i,j-1,k)*rho_t(:,i,j,k+1)*rho_t(:,i,j-1,k+1))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i,j-1,k)*  p_t(:,i,j,k+1)*  p_t(:,i,j-1,k+1))) - p
              dry = qLT(:,i,j,k,ir,1) - r
              dpy = qLT(:,i,j,k,ip,1) - p
              qLT(:,i,j,k,ir,1) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qLT(:,i,j,k,ip,1) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! X-edge averaged left-bottom corner state (LB->RR)
              !qLB(:,i,j,k,ir,1) = r + (-dry-drz)
              drx = sqrt(sqrt(r*rho_t(:,i,j-1,k)*rho_t(:,i,j,k-1)*rho_t(:,i,j-1,k-1))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i,j-1,k)*  p_t(:,i,j,k-1)*  p_t(:,i,j-1,k-1))) - p
              dry = qLB(:,i,j,k,ir,1) - r
              dpy = qLB(:,i,j,k,ip,1) - p
              qLB(:,i,j,k,ir,1) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qLB(:,i,j,k,ip,1) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! Y-edge averaged right-top corner state (RT->LL)
              !qRT(:,i,j,k,ir,2) = r + (+drx+drz)
              drx = sqrt(sqrt(r*rho_t(:,i+1,j,k)*rho_t(:,i,j,k+1)*rho_t(:,i+1,j,k+1))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i+1,j,k)*  p_t(:,i,j,k+1)*  p_t(:,i+1,j,k+1))) - p
              dry = qRT(:,i,j,k,ir,2) - r
              dpy = qRT(:,i,j,k,ip,2) - p
              qRT(:,i,j,k,ir,2) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qRT(:,i,j,k,ip,2) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! Y-edge averaged right-bottom corner state (RB->LR)
              !qRB(:,i,j,k,ir,2) = r + (+drx-drz)
              drx = sqrt(sqrt(r*rho_t(:,i+1,j,k)*rho_t(:,i,j,k-1)*rho_t(:,i+1,j,k-1))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i+1,j,k)*  p_t(:,i,j,k-1)*  p_t(:,i+1,j,k-1))) - p
              dry = qRB(:,i,j,k,ir,2) - r
              dpy = qRB(:,i,j,k,ip,2) - p
              qRB(:,i,j,k,ir,2) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qRB(:,i,j,k,ip,2) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! Y-edge averaged left-top corner state (LT->RL)
              !qLT(:,i,j,k,ir,2) = r + (-drx+drz)
              drx = sqrt(sqrt(r*rho_t(:,i-1,j,k)*rho_t(:,i,j,k+1)*rho_t(:,i-1,j,k+1))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i-1,j,k)*  p_t(:,i,j,k+1)*  p_t(:,i-1,j,k+1))) - p
              dry = qLT(:,i,j,k,ir,2) - r
              dpy = qLT(:,i,j,k,ip,2) - p
              qLT(:,i,j,k,ir,2) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qLT(:,i,j,k,ip,2) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! Y-edge averaged left-bottom corner state (LB->RR)
              !qLB(:,i,j,k,ir,2) = r + (-drx-drz)
              drx = sqrt(sqrt(r*rho_t(:,i-1,j,k)*rho_t(:,i,j,k-1)*rho_t(:,i-1,j,k-1))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i-1,j,k)*  p_t(:,i,j,k-1)*  p_t(:,i-1,j,k-1))) - p
              dry = qLB(:,i,j,k,ir,2) - r
              dpy = qLB(:,i,j,k,ip,2) - p
              qLB(:,i,j,k,ir,2) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qLB(:,i,j,k,ip,2) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! Z-edge averaged right-top corner state (RT->LL)
              !qRT(:,i,j,k,ir,3) = r + (+drx+dry)
              drx = sqrt(sqrt(r*rho_t(:,i+1,j,k)*rho_t(:,i,j+1,k)*rho_t(:,i+1,j+1,k))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i+1,j,k)*  p_t(:,i,j+1,k)*  p_t(:,i+1,j+1,k))) - p
              dry = qRT(:,i,j,k,ir,3) - r
              dpy = qRT(:,i,j,k,ip,3) - p
              qRT(:,i,j,k,ir,3) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qRT(:,i,j,k,ip,3) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! Z-edge averaged right-bottom corner state (RB->LR)
              !qRB(:,i,j,k,ir,3) = r + (+drx-dry)
              drx = sqrt(sqrt(r*rho_t(:,i+1,j,k)*rho_t(:,i,j-1,k)*rho_t(:,i+1,j-1,k))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i+1,j,k)*  p_t(:,i,j-1,k)*  p_t(:,i+1,j-1,k))) - p
              dry = qRB(:,i,j,k,ir,3) - r
              dpy = qRB(:,i,j,k,ip,3) - p
              qRB(:,i,j,k,ir,3) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qRB(:,i,j,k,ip,3) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! Z-edge averaged left-top corner state (LT->RL)
              !qLT(:,i,j,k,ir,3) = r + (-drx+dry)
              drx = sqrt(sqrt(r*rho_t(:,i-1,j,k)*rho_t(:,i,j+1,k)*rho_t(:,i-1,j+1,k))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i-1,j,k)*  p_t(:,i,j+1,k)*  p_t(:,i-1,j+1,k))) - p
              dry = qLT(:,i,j,k,ir,3) - r
              dpy = qLT(:,i,j,k,ip,3) - p
              qLT(:,i,j,k,ir,3) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qLT(:,i,j,k,ip,3) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              ! Z-edge averaged left-bottom corner state (LB->RR)
              !qLB(:,i,j,k,ir,3) = r + (-drx-dry)
              drx = sqrt(sqrt(r*rho_t(:,i-1,j,k)*rho_t(:,i,j-1,k)*rho_t(:,i-1,j-1,k))) - r
              dpx = sqrt(sqrt(p*  p_t(:,i-1,j,k)*  p_t(:,i,j-1,k)*  p_t(:,i-1,j-1,k))) - p
              dry = qLB(:,i,j,k,ir,3) - r
              dpy = qLB(:,i,j,k,ip,3) - p
              qLB(:,i,j,k,ir,3) = merge(drx, dry, abs(drx) <= abs(dry) .or. dry+r<=0.) + r
              qLB(:,i,j,k,ip,3) = merge(dpx, dpy, abs(dpx) <= abs(dpy) .or. dpy+p<=0.) + p

              t1 = qRT(:,i,j,k,ir,1)-r; t3 = qRT(:,i,j,k,ir,2)-r; t5 = qRT(:,i,j,k,ir,3)-r
              t2 = qRT(:,i,j,k,ip,1)-p; t4 = qRT(:,i,j,k,ip,2)-p; t6 = qRT(:,i,j,k,ip,3)-p
              call check_negative(qRT,i,j,k,t1,t2,t3,t4,t5,t6,'qRT',ngrid)
              t1 = qRB(:,i,j,k,ir,1)-r; t3 = qRB(:,i,j,k,ir,2)-r; t5 = qRB(:,i,j,k,ir,3)-r
              t2 = qRB(:,i,j,k,ip,1)-p; t4 = qRB(:,i,j,k,ip,2)-p; t6 = qRB(:,i,j,k,ip,3)-p
              call check_negative(qRB,i,j,k,t1,t2,t3,t4,t5,t6,'qRB',ngrid)
              t1 = qLT(:,i,j,k,ir,1)-r; t3 = qLT(:,i,j,k,ir,2)-r; t5 = qLT(:,i,j,k,ir,3)-r
              t2 = qLT(:,i,j,k,ip,1)-p; t4 = qLT(:,i,j,k,ip,2)-p; t6 = qLT(:,i,j,k,ip,3)-p
              call check_negative(qLT,i,j,k,t1,t2,t3,t4,t5,t6,'qLT',ngrid)
              t1 = qLB(:,i,j,k,ir,1)-r; t3 = qLB(:,i,j,k,ir,2)-r; t5 = qLB(:,i,j,k,ir,3)-r
              t2 = qLB(:,i,j,k,ip,1)-p; t4 = qLB(:,i,j,k,ip,2)-p; t6 = qLB(:,i,j,k,ip,3)-p
              call check_negative(qLB,i,j,k,t1,t2,t3,t4,t5,t6,'qLB',ngrid)
        END DO
     END DO
  END DO

  ! Passive scalars
  DO n = 9, nvar
     DO k = klo, khi
        DO j = jlo, jhi
           DO i = ilo, ihi
                 r   = q(:,i,j,k,n )            ! Cell centered values
                 u   = q(:,i,j,k,iu)
                 v   = q(:,i,j,k,iv)
                 w   = q(:,i,j,k,iw)
                 drx = half * dq(:,i,j,k,n,1)   ! TVD slopes
                 dry = half * dq(:,i,j,k,n,2)
                 drz = half * dq(:,i,j,k,n,3)

                 sr0 = -u*drx*dtdx -v*dry*dtdy -w*drz*dtdz   ! Source terms
#ifdef KROME
                 if (n .ne. ichem) &            ! No time evolution for adiabatic index
#endif
                    r   = r + sr0               ! Predicted state

                 qp(:,i,j,k,n,1) = r - drx      ! Right state
                 qm(:,i,j,k,n,1) = r + drx      ! Left state
                 qp(:,i,j,k,n,2) = r - dry      ! Top state
                 qm(:,i,j,k,n,2) = r + dry      ! Bottom state
                 qp(:,i,j,k,n,3) = r - drz      ! Front state
                 qm(:,i,j,k,n,3) = r + drz      ! Back state

#ifdef KROME
                 ! Make sure interpolated abundances are positive definite, and skip time evolution if not the case
                 if (n .ne. ichem) then
                    r = chem_t(:,i,j,k,n)
                    if (debug) then
                       do l=1,ngrid
                          if (chem_t(l,i,j,k,n) < 0.) write(*,'(a,i3,1p,3e11.4)') &
                           'Chemical density < 0.:', n, chem_t(l,i,j,k,n), q(l,i,j,k,n), chem_t(l,i,j,k,n)-q(l,i,j,k,n)
                       enddo
                    endif

                    !qp(:,i,j,k,n,1) = r - drx      ! Right state
                    qp(:,i,j,k,n,1) = merge(sqrt(r*chem_t(:,i-1,j,k,n)), qp(:,i,j,k,n,1), qp(:,i,j,k,n,1)<0.)

                    !qm(:,i,j,k,n,1) = r + drx      ! Left state
                    qm(:,i,j,k,n,1) = merge(sqrt(r*chem_t(:,i+1,j,k,n)), qm(:,i,j,k,n,1), qm(:,i,j,k,n,1)<0.)

                    !qp(:,i,j,k,n,2) = r - dry      ! Top state
                    qp(:,i,j,k,n,2) = merge(sqrt(r*chem_t(:,i,j-1,k,n)), qp(:,i,j,k,n,2), qp(:,i,j,k,n,2)<0.)

                    !qm(:,i,j,k,n,2) = r + dry      ! Bottom state
                    qm(:,i,j,k,n,2) = merge(sqrt(r*chem_t(:,i,j+1,k,n)), qm(:,i,j,k,n,2), qm(:,i,j,k,n,2)<0.)

                    !qp(:,i,j,k,n,3) = r - drz      ! Front state
                    qp(:,i,j,k,n,3) = merge(sqrt(r*chem_t(:,i,j,k-1,n)), qp(:,i,j,k,n,3), qp(:,i,j,k,n,3)<0.)

                    !qm(:,i,j,k,n,3) = r + drz      ! Back state
                    qm(:,i,j,k,n,3) = merge(sqrt(r*chem_t(:,i,j,k+1,n)), qm(:,i,j,k,n,3), qm(:,i,j,k,n,3)<0.)
                 end if

                 ! Find adiabatic index at corner states
                 if (chemistry .and. n==ichem) then ! Adiabatic index
                    qRT(:,i,j,k,ichem,1) = r + (+dry+drz) ! X-edge averaged right-top corner state (RT->LL)
                    qRB(:,i,j,k,ichem,1) = r + (+dry-drz) ! X-edge averaged right-bottom corner state (RB->LR)
                    qLT(:,i,j,k,ichem,1) = r + (-dry+drz) ! X-edge averaged left-top corner state (LT->RL)
                    qLB(:,i,j,k,ichem,1) = r + (-dry-drz) ! X-edge averaged left-bottom corner state (LB->RR)
                    qRT(:,i,j,k,ichem,2) = r + (+drx+drz) ! Y-edge averaged right-top corner state (RT->LL)
                    qRB(:,i,j,k,ichem,2) = r + (+drx-drz) ! Y-edge averaged right-bottom corner state (RB->LR)
                    qLT(:,i,j,k,ichem,2) = r + (-drx+drz) ! Y-edge averaged left-top corner state (LT->RL)
                    qLB(:,i,j,k,ichem,2) = r + (-drx-drz) ! Y-edge averaged left-bottom corner state (LB->RR)
                    qRT(:,i,j,k,ichem,3) = r + (+drx+dry) ! Z-edge averaged right-top corner state (RT->LL)
                    qRB(:,i,j,k,ichem,3) = r + (+drx-dry) ! Z-edge averaged right-bottom corner state (RB->LR)
                    qLT(:,i,j,k,ichem,3) = r + (-drx+dry) ! Z-edge averaged left-top corner state (LT->RL)
                    qLB(:,i,j,k,ichem,3) = r + (-drx-dry) ! Z-edge averaged left-bottom corner state (LB->RR)
                 endif
#endif
           END DO
        END DO
     END DO
  END DO

#ifdef KROME
  ! Check that abundances have been consistently slope interpolated,
  ! in the sense that they sum to 1 at the cell interface -- q is a primitive variable.
  ! We use the rescaling formalism of Plewa and Muller A&A 342, 179 (1999)

  ! First find the total abundances, according to the chemical species
  totmass = 0.0_dp
  DO idir = 1, 3
  DO n = ichem+1, nvar
  DO k = klo, khi
  DO j = jlo, jhi
  DO i = ilo, ihi
     totmass(:,i,j,k,idir  ) = totmass(:,i,j,k,idir  ) + qm(:,i,j,k,n,idir)
     totmass(:,i,j,k,idir+3) = totmass(:,i,j,k,idir+3) + qp(:,i,j,k,n,idir)
  END DO
  END DO
  END DO
  END DO
  END DO

  ! Compute rescaling factor to make abundances sum to 1
  DO idir = 1, 3
  DO k = klo, khi
  DO j = jlo, jhi
  DO i = ilo, ihi
    totmass(:,i,j,k,idir)   = 1.0_dp / totmass(:,i,j,k,idir  )
    totmass(:,i,j,k,idir+3) = 1.0_dp / totmass(:,i,j,k,idir+3)
  END DO
  END DO
  END DO
  END DO

  if (debug .and. nprintX>0) then
     do l=1,ngrid
     do idir = 1, 3
     do k = klo, khi
     do j = jlo, jhi
     do i = ilo, ihi
        if (nprintX<=0) cycle
        if (totmass(l,i,j,k,idir) < 0.9 .or. totmass(l,i,j,k,idir) > 1.1) then
           !$omp critical
           if (nprintX>0) then
              nprintX = nprintX - 1
              print '(a,5i4,4e12.2)', 'Large rescaling of abundances QM l,i,j,k,dir,totX,min(X),max(X):', l, i, j, k, idir,  &
                 totmass(l,i,j,k,idir), minval(qm(l,i,j,k,ichem+1:nvar,idir)), maxval(qm(l,i,j,k,ichem+1:nvar,idir))
           end if
           !$omp end critical
        end if
        if (totmass(l,i,j,k,idir+3) < 0.9 .or. totmass(l,i,j,k,idir+3) > 1.1) then
           !$omp critical
           if (nprintX>0) then
              nprintX = nprintX - 1
              print '(a,5i4,4e12.2)', 'Large rescaling of abundances QP l,i,j,k,dir,totX,rho,min(X),max(X):', l, i, j, k, idir,  &
                 totmass(l,i,j,k,idir+3), minval(qp(l,i,j,k,ichem+1:nvar,idir)), maxval(qp(l,i,j,k,ichem+1:nvar,idir))
           end if
           !$omp end critical
        end if
     end do
     end do
     end do
     end do
     end do
  end if

  ! Apply rescaling factor
  DO idir = 1, 3
  DO n = ichem+1, nvar
  DO k = klo, khi
  DO j = jlo, jhi
  DO i = ilo, ihi
     qm(1:ngrid,i,j,k,n,idir) = totmass(1:ngrid,i,j,k,idir  ) * qm(1:ngrid,i,j,k,n,idir)
     qp(1:ngrid,i,j,k,n,idir) = totmass(1:ngrid,i,j,k,idir+3) * qp(1:ngrid,i,j,k,n,idir)
  END DO
  END DO
  END DO
  END DO
  END DO
#endif

!  call trace_mpi('trace3d','exit',2)

CONTAINS
  subroutine check_negative(q,i,j,k,rx,px,ry,py,rz,pz,label,ngrid)
    real(dp), dimension(nvector,iu1:iu2,ju1:ju2,ku1:ku2,nvar,3), intent(inout) :: q
    real(dp), dimension(nvector),        intent(in)    :: rx,ry,rz,px,py,pz
    real(dp), dimension(nvector)                       :: r,p
    integer                                            :: l,idim,ngrid
    character(len=*)                                   :: label
    integer, intent(in) :: i, j, k
    integer, save:: nprint=20
!    call trace_mpi('check_negative','enter',3)
    if (debug .and. nprint>0) then
       !$omp atomic
       nprint=nprint-1
       do idim=1,3
          if (idim==1) then; r=rx; p=px; endif
          if (idim==2) then; r=ry; p=py; endif
          if (idim==3) then; r=rz; p=pz; endif
          do l=1,ngrid
             if (q(l,i,j,k,1,idim) <= smallr) then
                write(*,2) label//' Small density found!', q(l,i,j,k,1,idim)-r(l),r(l),smallr,q(l,i,j,k,1,idim),idim,i,j,k
              2 format(1x,a,1p,4e14.6,4i3) 
                q(l,i,j,k,1,idim) = smallr 
                !stop
             end if
             if (q(l,i,j,k,5,idim) <= smallp) then
                write(*,2) label//' Small pressure found!',q(l,i,j,k,5,idim)-p(l),p(l),smallp,q(l,i,j,k,5,idim),idim,i,j,k
                q(l,i,j,k,5,idim) = smallp
                !stop
             end if
          end do
       end do
    else
       do idim=1,3
          q(:,i,j,k,1,idim) = max(smallr,q(:,i,j,k,1,idim))
          q(:,i,j,k,5,idim) = max(smallp,q(:,i,j,k,5,idim))
       end do
    end if
!    call trace_mpi('check_negative','exit',3)
  end subroutine check_negative
END SUBROUTINE trace3d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, &
     &                ln ,lt1,lt2,bn ,bt1,bt2, &
     &                mark_low_order,flx,tmp,ngrid,ind_grid)
  use amr_parameters
  use amr_commons, only: local_debug, myid, xg
  use hydro_parameters
  use hydro_commons, only: courant_number
  use const
  !use mpi_base, only: clean_stop
  implicit none

  integer ::ngrid
  integer ::ind_grid(nvector)
  integer ::ln,lt1,lt2,bn,bt1,bt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp 
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar)::flx
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:2)::tmp
  logical ,dimension(1:nvector,twotondim,twondim)::mark_low_order
  
  ! local variables
  integer ::i, j, k, n, l, idim, xdim
  real(dp),dimension(1:nvector,1:nvar),target ::qleft,qright,fgdnv
  real(dp),dimension(1:nvar)::qlefte,qrighte,fgdnve
  real(dp),dimension(1:nvector)::bn_mean
  logical ,dimension(1:nvector)::mark_low_order_vec
  real(dp)::zero_flux, entho, max_speed

  real(dp)::xx=0.0
  logical:: dbg=.false.
  logical in_debug_region


!  call trace_mpi('cmpflxm','enter',2)

  xdim=ln-1
  entho=one/(gamma-one)
  max_speed=0.0_dp

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi

            ! Enforce continuity for normal magnetic field
            bn_mean(1:ngrid) = half*(qm(1:ngrid,i,j,k,bn,xdim)+qp(1:ngrid,i,j,k,bn,xdim))

            ! Left state
            qleft (1:ngrid,1) = qm(1:ngrid,i,j,k,1  ,xdim) ! Mass density
            qleft (1:ngrid,2) = qm(1:ngrid,i,j,k,5  ,xdim) ! Pressure
            qleft (1:ngrid,3) = qm(1:ngrid,i,j,k,ln ,xdim) ! Normal velocity
            qleft (1:ngrid,4) = bn_mean(1:ngrid)              ! Normal magnetic field
            qleft (1:ngrid,5) = qm(1:ngrid,i,j,k,lt1,xdim) ! Tangential velocity 1
            qleft (1:ngrid,6) = qm(1:ngrid,i,j,k,bt1,xdim) ! Tangential magnetic field 1
            qleft (1:ngrid,7) = qm(1:ngrid,i,j,k,lt2,xdim) ! Tangential velocity 2
            qleft (1:ngrid,8) = qm(1:ngrid,i,j,k,bt2,xdim) ! Tangential magnetic field 2

            ! Right state
            qright(1:ngrid,1) = qp(1:ngrid,i,j,k,1  ,xdim) ! Mass density
            qright(1:ngrid,2) = qp(1:ngrid,i,j,k,5  ,xdim) ! Pressure
            qright(1:ngrid,3) = qp(1:ngrid,i,j,k,ln ,xdim) ! Normal velocity
            qright(1:ngrid,4) = bn_mean(1:ngrid)           ! Normal magnetic field
            qright(1:ngrid,5) = qp(1:ngrid,i,j,k,lt1,xdim) ! Tangential velocity 1
            qright(1:ngrid,6) = qp(1:ngrid,i,j,k,bt1,xdim) ! Tangential magnetic field 1
            qright(1:ngrid,7) = qp(1:ngrid,i,j,k,lt2,xdim) ! Tangential velocity 2
            qright(1:ngrid,8) = qp(1:ngrid,i,j,k,bt2,xdim) ! Tangential magnetic field 2

            ! Vector marking for low order
            if (cbmax > 0) then
               select case (ln)
               case(2) ! X-direction
                  if (i == 1) mark_low_order_vec(1:ngrid) = mark_low_order(1:ngrid,1+2*(j-1+2*(k-1)),1)
                  if (i == 2) mark_low_order_vec(1:ngrid) = mark_low_order(1:ngrid,1+2*(j-1+2*(k-1)),2)
                  if (i == 3) mark_low_order_vec(1:ngrid) = mark_low_order(1:ngrid,2+2*(j-1+2*(k-1)),2)
#if NDIM>1
               case(3) ! Y-direction                                        
                  if (j == 1) mark_low_order_vec(1:ngrid) = mark_low_order(1:ngrid,i+2*(0+2*(k-1)),3)
                  if (j == 2) mark_low_order_vec(1:ngrid) = mark_low_order(1:ngrid,i+2*(0+2*(k-1)),4)
                  if (j == 3) mark_low_order_vec(1:ngrid) = mark_low_order(1:ngrid,i+2*(1+2*(k-1)),4)
#endif
#if NDIM>2
               case(4) ! Z-direction                                        
                  if (k == 1) mark_low_order_vec(1:ngrid) = mark_low_order(1:ngrid,i+2*(j-1+2*0),5)
                  if (k == 2) mark_low_order_vec(1:ngrid) = mark_low_order(1:ngrid,i+2*(j-1+2*0),6)
                  if (k == 3) mark_low_order_vec(1:ngrid) = mark_low_order(1:ngrid,i+2*(j-1+2*1),6)
#endif
               end select
            endif

            ! Other advected quantities
            do n = 9, nvar
               qleft (1:ngrid,n) = qm(1:ngrid,i,j,k,n,xdim)
               qright(1:ngrid,n) = qp(1:ngrid,i,j,k,n,xdim)    
            end do
          
            ! Solve 1D Riemann problem
            zero_flux = one
            IF(scheme.NE.'induction')THEN
            SELECT CASE (riemann)
            CASE ('roe')
               do l = 1, ngrid
                  qlefte = qleft(l,:); qrighte=qright(l,:)
                  CALL athena_roe     (qlefte,qrighte,fgdnve,zero_flux)
                  fgdnv(l,:) = fgdnve
               end do
            CASE ('llf')
               CALL lax_friedrich_vec  (qleft,qright,fgdnv,zero_flux)
            CASE ('hll')
               do l = 1, ngrid
                  qlefte = qleft(l,:); qrighte=qright(l,:)
                  CALL hll            (qlefte,qrighte,fgdnve)
                  fgdnv(l,:) = fgdnve
               end do
            CASE ('hlld')
               do l = 1, ngrid
                  qlefte = qleft(l,:); qrighte=qright(l,:)
                  !xx = xg(ind_grid(l),1) + (i-1.5)/32.
                  !dbg = local_debug.and.xdim==1.and.in_debug_region(ind_grid(l)).and.j==1.and.k==1
                  CALL hlld           (qlefte,qrighte,fgdnve,mark_low_order_vec(l),dbg,xx)
                  fgdnv(l,:) = fgdnve
                  if (dbg) then
                     print'(a,1p,12(e13.6,1h,))',' qlefte=',qlefte
                     print'(a,1p,12(e13.6,1h,))','qrighte=',qrighte
                     print'(a,1p,12(e13.6,1h,))','  fgdnv=',fgdnve
                  end if
               end do
            CASE ('hll3r')
               CALL hll3r(qleft,qright,fgdnv,ngrid,max_speed)
            CASE ('hll5r')
               CALL hll5r(qleft,qright,fgdnv,ngrid,max_speed)
            CASE ('upwind')
               CALL lax_friedrich_vec (qleft,qright,fgdnv,zero_flux)
            CASE ('hydro')
               do l = 1, ngrid
                  qlefte = qleft(l,:); qrighte=qright(l,:)
                  CALL hydro_acoustic (qlefte,qrighte,fgdnve)
                  fgdnv(l,:) = fgdnve
               end do
            CASE DEFAULT
               write(*,*)'unknown riemann solver'
               stop
            END SELECT
            ELSE
               do l = 1, ngrid
                  qlefte = qleft(l,:); qrighte=qright(l,:)
                  CALL upwind(qlefte,qrighte,fgdnve,zero_flux)
                  fgdnv(l,:) = fgdnve
               end do
            ENDIF
           
            ! Output fluxes
            flx(1:ngrid,i,j,k,1  ) = fgdnv(1:ngrid,1)  ! Mass density
            flx(1:ngrid,i,j,k,5  ) = fgdnv(1:ngrid,2)  ! Total energy
            flx(1:ngrid,i,j,k,ln ) = fgdnv(1:ngrid,3)  ! Normal momentum
            flx(1:ngrid,i,j,k,bn ) = fgdnv(1:ngrid,4)  ! Normal magnetic field
            flx(1:ngrid,i,j,k,lt1) = fgdnv(1:ngrid,5)  ! Transverse momentum 1
            flx(1:ngrid,i,j,k,bt1) = fgdnv(1:ngrid,6)  ! Transverse magnetic field 1
            flx(1:ngrid,i,j,k,lt2) = fgdnv(1:ngrid,7)  ! Transverse momentum 2
            flx(1:ngrid,i,j,k,bt2) = fgdnv(1:ngrid,8)  ! Transverse magnetic field 2

            ! Other advected quantities
            do n = 9, nvar
               flx(1:ngrid,i,j,k,n) = fgdnv(1:ngrid,n)
            end do

#ifdef KROME
            if (debug) then
               do l=1,ngrid
                  if (fgdnv(l,1) .ne. 0.) then
                     if (abs(sum(fgdnv(l,ichem+1:nvar)) / fgdnv(l,1) -1.) > 1e-6_dp) print '(a,i5,3e12.3)', &
                        'cmpflxm Abundance problem :', myid, sum(fgdnv(l,ichem+1:nvar)), fgdnv(l,1), &
                        sum(fgdnv(l,ichem+1:nvar)) / fgdnv(l,1) -1.
                  end if
               end do
            end if
#endif
            
            ! Normal velocity estimate
            tmp(1:ngrid,i,j,k,1) = half*(qleft(1:ngrid,3)+qright(1:ngrid,3))
            ! Internal energy flux
#ifdef KROME
            if (chemistry) then
            tmp(1:ngrid,i,j,k,2) = merge(qleft (1:ngrid,2)/max(qleft (1:ngrid,1),smallr)/(qleft (1:ngrid,ichem)-one)*fgdnv(1:ngrid,1), &
                                         qright(1:ngrid,2)/max(qright(1:ngrid,1),smallr)/(qright(1:ngrid,ichem)-one)*fgdnv(1:ngrid,1), &
                                         fgdnv(1:ngrid,1)>zero)
            else
#endif
            tmp(1:ngrid,i,j,k,2) = merge(qleft (1:ngrid,2)/max(qleft (1:ngrid,1),smallr)*entho*fgdnv(1:ngrid,1), &
                                         qright(1:ngrid,2)/max(qright(1:ngrid,1),smallr)*entho*fgdnv(1:ngrid,1), &
                                         fgdnv(1:ngrid,1)>zero)
#ifdef KROME
            endif
#endif

        end do
     end do
  end do

  if (trim(riemann)=='hll3r' .or. trim(riemann)=='hll5r') then
    !$omp atomic
    courant_number%speed = max(courant_number%speed,max_speed)
  endif
!  call trace_mpi('cmpflxm','exit',2)
end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cmp_mag_flx(qRT,irt1,irt2,jrt1,jrt2,krt1,krt2, &
       &               qRB,irb1,irb2,jrb1,jrb2,krb1,krb2, &
       &               qLT,ilt1,ilt2,jlt1,jlt2,klt1,klt2, &
       &               qLB,ilb1,ilb2,jlb1,jlb2,klb1,klb2, &
       &                   ilo ,ihi ,jlo ,jhi ,klo ,khi , &
       &                   lp1 ,lp2 ,lor ,bp1 ,bp2 ,bor , &
       &               mark_edge,emf,ngrid)
  ! 2D Riemann solver to compute EMF at cell edges
  USE amr_parameters
  USE amr_commons, only : myid
  USE hydro_parameters
  USE const
  !USE mpi_base, only: clean_stop
  IMPLICIT NONE

  INTEGER ::ngrid
  ! indices of the 2 planar velocity lp1 and lp2 and the orthogonal one,
  ! lor and idem for the magnetic field
  INTEGER ::lp1,lp2,lor,bp1,bp2,bor
  INTEGER ::irt1,irt2,jrt1,jrt2,krt1,krt2
  INTEGER ::irb1,irb2,jrb1,jrb2,krb1,krb2
  INTEGER ::ilt1,ilt2,jlt1,jlt2,klt1,klt2
  INTEGER ::ilb1,ilb2,jlb1,jlb2,klb1,klb2
  INTEGER ::ilo,ihi,jlo,jhi,klo,khi
  REAL(dp),DIMENSION(1:nvector,irt1:irt2,jrt1:jrt2,krt1:krt2,1:nvar,1:3)::qRT
  REAL(dp),DIMENSION(1:nvector,irb1:irb2,jrb1:jrb2,krb1:krb2,1:nvar,1:3)::qRB
  REAL(dp),DIMENSION(1:nvector,ilt1:ilt2,jlt1:jlt2,klt1:klt2,1:nvar,1:3)::qLT
  REAL(dp),DIMENSION(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2,1:nvar,1:3)::qLB

  LOGICAL, DIMENSION(1:nvector,3,3,3)              :: mark_edge
  REAL(dp),DIMENSION(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2):: emf

  ! local variables
  INTEGER ::i, j, k, n, l, idim, xdim, m
  REAL(dp),DIMENSION(1:nvector,1:nvar)::qLL,qRL,qLR,qRR
  REAL(dp),DIMENSION(1:nvector,1:nvar)::qleftv,qrightv,fmean_xv,fmean_yv
  REAL(dp),DIMENSION(1:nvector)       ::Ev
  REAL(dp),DIMENSION(1:nvar),target::qleft,qright,fmean_x,fmean_y,qtmp
  REAL(dp) :: ELL,ERL,ELR,ERR,S,SL,SR,SB,ST,SAL,SAR,SAT,SAB
  REAL(dp) :: zero_flux,E,alpha
  REAL(dp) :: cLLx,cRLx,cLRx,cRRx,cLLy,cRLy,cLRy,cRRy
  REAL(dp) :: cfastLLx,cfastRLx,cfastLRx,cfastRRx,cfastLLy,cfastRLy,cfastLRy,cfastRRy
  REAL(dp) :: calfvenR,calfvenL,calfvenT,calfvenB
  REAL(dp) :: vLLx,vRLx,vLRx,vRRx,vLLy,vRLy,vLRy,vRRy
  REAL(dp) :: rLL,rLR,rRL,rRR,pLL,pLR,pRL,pRR,uLL,uLR,uRL,uRR,vLL,vLR,vRL,vRR
  REAL(dp) :: ALL,ALR,ARL,ARR,BLL,BLR,BRL,BRR,CLL,CLR,CRL,CRR
  REAL(dp) :: PtotLL,PtotLR,PtotRL,PtotRR,rcLLx,rcLRx,rcRLx,rcRRx,rcLLy,rcLRy,rcRLy,rcRRy
  REAL(dp) :: ustar,vstar,rstarLLx,rstarLRx,rstarRLx,rstarRRx,rstarLLy,rstarLRy,rstarRLy,rstarRRy
  REAL(dp) :: rstarLL,rstarLR,rstarRL,rstarRR,AstarLL,AstarLR,AstarRL,AstarRR,BstarLL,BstarLR,BstarRL,BstarRR
  REAL(dp) :: EstarLLx,EstarLRx,EstarRLx,EstarRRx,EstarLLy,EstarLRy,EstarRLy,EstarRRy,EstarLL,EstarLR,EstarRL,EstarRR
  REAL(dp) :: AstarT,AstarB,BstarR,BstarL,lmLx,lmRx,lmLy,lmRy,lpLx,lpRx,lpLy,lpRy
  CHARACTER(LEN=1), DIMENSION(3) :: dir=(/'x','y','z'/)
  xdim = lor - 1

!  call trace_mpi('cmp_mag_flx','enter',2)

  DO k = klo, khi
     DO j = jlo, jhi
        DO i = ilo, ihi

           ! Density
!           DO l = 1, ngrid
              qLL (1:ngrid,1) = qRT(1:ngrid,i,j,k,1,xdim)
              qRL (1:ngrid,1) = qLT(1:ngrid,i,j,k,1,xdim)
              qLR (1:ngrid,1) = qRB(1:ngrid,i,j,k,1,xdim)
              qRR (1:ngrid,1) = qLB(1:ngrid,i,j,k,1,xdim)
!           END DO           

           ! Pressure 
!           DO l = 1, ngrid
              qLL (1:ngrid,2) = qRT(1:ngrid,i,j,k,5,xdim)
              qRL (1:ngrid,2) = qLT(1:ngrid,i,j,k,5,xdim)
              qLR (1:ngrid,2) = qRB(1:ngrid,i,j,k,5,xdim)
              qRR (1:ngrid,2) = qLB(1:ngrid,i,j,k,5,xdim)
!           END DO

           ! First parallel velocity 
!           DO l = 1, ngrid
              qLL (1:ngrid,3) = qRT(1:ngrid,i,j,k,lp1,xdim)
              qRL (1:ngrid,3) = qLT(1:ngrid,i,j,k,lp1,xdim)
              qLR (1:ngrid,3) = qRB(1:ngrid,i,j,k,lp1,xdim)
              qRR (1:ngrid,3) = qLB(1:ngrid,i,j,k,lp1,xdim)
!           END DO

           ! Second parallel velocity 
!           DO l = 1, ngrid
              qLL (1:ngrid,4) = qRT(1:ngrid,i,j,k,lp2,xdim)
              qRL (1:ngrid,4) = qLT(1:ngrid,i,j,k,lp2,xdim)
              qLR (1:ngrid,4) = qRB(1:ngrid,i,j,k,lp2,xdim)
              qRR (1:ngrid,4) = qLB(1:ngrid,i,j,k,lp2,xdim)
!           END DO

           ! First parallel magnetic field (enforce continuity)
!           DO l = 1, ngrid
              qLL (1:ngrid,6) = half*(qRT(1:ngrid,i,j,k,bp1,xdim)+qLT(1:ngrid,i,j,k,bp1,xdim))
              qRL (1:ngrid,6) = half*(qRT(1:ngrid,i,j,k,bp1,xdim)+qLT(1:ngrid,i,j,k,bp1,xdim))
              qLR (1:ngrid,6) = half*(qRB(1:ngrid,i,j,k,bp1,xdim)+qLB(1:ngrid,i,j,k,bp1,xdim))
              qRR (1:ngrid,6) = half*(qRB(1:ngrid,i,j,k,bp1,xdim)+qLB(1:ngrid,i,j,k,bp1,xdim))
!           END DO

           ! Second parallel magnetic field (enforce continuity)
!           DO l = 1, ngrid
              qLL (1:ngrid,7) = half*(qRT(1:ngrid,i,j,k,bp2,xdim)+qRB(1:ngrid,i,j,k,bp2,xdim))
              qRL (1:ngrid,7) = half*(qLT(1:ngrid,i,j,k,bp2,xdim)+qLB(1:ngrid,i,j,k,bp2,xdim))
              qLR (1:ngrid,7) = half*(qRT(1:ngrid,i,j,k,bp2,xdim)+qRB(1:ngrid,i,j,k,bp2,xdim))
              qRR (1:ngrid,7) = half*(qLT(1:ngrid,i,j,k,bp2,xdim)+qLB(1:ngrid,i,j,k,bp2,xdim))
!           END DO

           ! Orthogonal velocity 
!           DO l = 1, ngrid
              qLL (1:ngrid,5) = qRT(1:ngrid,i,j,k,lor,xdim)
              qRL (1:ngrid,5) = qLT(1:ngrid,i,j,k,lor,xdim)
              qLR (1:ngrid,5) = qRB(1:ngrid,i,j,k,lor,xdim)
              qRR (1:ngrid,5) = qLB(1:ngrid,i,j,k,lor,xdim)
!           END DO

           ! Orthogonal magnetic Field
!           DO l = 1, ngrid
              qLL (1:ngrid,8) = qRT(1:ngrid,i,j,k,bor,xdim)
              qRL (1:ngrid,8) = qLT(1:ngrid,i,j,k,bor,xdim)
              qLR (1:ngrid,8) = qRB(1:ngrid,i,j,k,bor,xdim)
              qRR (1:ngrid,8) = qLB(1:ngrid,i,j,k,bor,xdim)
!           END DO

#ifdef KROME
           ! Effective adiabatic index
           if (chemistry) then
              qLL (1:ngrid,ichem) = qRT(1:ngrid,i,j,k,ichem,xdim)
              qRL (1:ngrid,ichem) = qLT(1:ngrid,i,j,k,ichem,xdim)
              qLR (1:ngrid,ichem) = qRB(1:ngrid,i,j,k,ichem,xdim)
              qRR (1:ngrid,ichem) = qLB(1:ngrid,i,j,k,ichem,xdim)
           end if
#endif

           ! Compute final fluxes
           if (riemann2d=='llf') then
           ! Compute final fluxes
              DO l = 1, ngrid 

                 ! vx*by - vy*bx at the four edge centers
                 ELL = qLL(l,3)*qLL(l,7) - qLL(l,4)*qLL(l,6)
                 ERL = qRL(l,3)*qRL(l,7) - qRL(l,4)*qRL(l,6)
                 ELR = qLR(l,3)*qLR(l,7) - qLR(l,4)*qLR(l,6)  
                 ERR = qRR(l,3)*qRR(l,7) - qRR(l,4)*qRR(l,6) 

                 ! find the average value of E
                 Ev(l) = forth*(ELL+ERL+ELR+ERR)
              ENDDO
                
              ! call the first solver in the x direction
              ! density
              qleftv (:,1) = half*(qLL(:,1)+qLR(:,1))
              qrightv(:,1) = half*(qRR(:,1)+qRL(:,1))
              
              ! pressure
              qleftv (:,2) = half*(qLL(:,2)+qLR(:,2))
              qrightv(:,2) = half*(qRR(:,2)+qRL(:,2))
              
              ! vt1 becomes normal velocity
              qleftv (:,3) = half*(qLL(:,3)+qLR(:,3))
              qrightv(:,3) = half*(qRR(:,3)+qRL(:,3))
              
              ! bt1 becomes normal magnetic field
              qleftv (:,4) = half*(qLL(:,6)+qLR(:,6))
              qrightv(:,4) = half*(qRR(:,6)+qRL(:,6))
              
              ! vt2 becomes transverse velocity field
              qleftv (:,5) = half*(qLL(:,4)+qLR(:,4))
              qrightv(:,5) = half*(qRR(:,4)+qRL(:,4))
              
              ! bt2 becomes transverse magnetic field 
              qleftv (:,6) = half*(qLL(:,7)+qLR(:,7))
              qrightv(:,6) = half*(qRR(:,7)+qRL(:,7))
              
              ! velocity component perp. to the plane is now transverse
              qleftv (:,7) = half*(qLL(:,5)+qLR(:,5))
              qrightv(:,7) = half*(qRR(:,5)+qRL(:,5))
              
              ! magnetic field component perp. to the plane is now transverse
              qleftv (:,8) = half*(qLL(:,8)+qLR(:,8))
              qrightv(:,8) = half*(qRR(:,8)+qRL(:,8))
              zero_flux = 0.0

#ifdef KROME
              ! Effective adiabatic index
              if (chemistry) then
                 qleftv (:,ichem) = half*(qLL(:,ichem)+qLR(:,ichem))
                 qrightv(:,ichem) = half*(qRR(:,ichem)+qRL(:,ichem))
              end if
#endif
              CALL lax_friedrich_vec(qleftv,qrightv,fmean_xv,zero_flux)

              ! call the second solver in the y direction
              ! density
              qleftv (:,1) = half*(qLL(:,1)+qRL(:,1))
              qrightv(:,1) = half*(qRR(:,1)+qLR(:,1))
            
              ! pressure
              qleftv (:,2) = half*(qLL(:,2)+qRL(:,2))
              qrightv(:,2) = half*(qRR(:,2)+qLR(:,2))
            
              ! vt2 becomes normal velocity
              qleftv (:,3) = half*(qLL(:,4)+qRL(:,4))
              qrightv(:,3) = half*(qRR(:,4)+qLR(:,4))
              
              ! bt2 becomes normal magnetic field
              qleftv (:,4) = half*(qLL(:,7)+qRL(:,7))
              qrightv(:,4) = half*(qRR(:,7)+qLR(:,7))
              
              ! vt1 becomes transverse velocity field 
              qleftv (:,5) = half*(qLL(:,3)+qRL(:,3))
              qrightv(:,5) = half*(qRR(:,3)+qLR(:,3))
              
              ! bt1 becomes transverse magnetic field 
              qleftv (:,6) = half*(qLL(:,6)+qRL(:,6))
              qrightv(:,6) = half*(qRR(:,6)+qLR(:,6))
              
              ! velocity component perp. to the plane is now transverse
              qleftv (:,7) = half*(qLL(:,5)+qRL(:,5))
              qrightv(:,7) = half*(qRR(:,5)+qLR(:,5))
              
              ! magnetic field component perp. to the plane is now transverse
              qleftv (:,8) = half*(qLL(:,8)+qRL(:,8))
              qrightv(:,8) = half*(qRR(:,8)+qLR(:,8))
              
              zero_flux = 0.
#ifdef KROME
              ! Effective adiabatic index
              if (chemistry) then
                 qleftv (:,ichem) = half*(qLL(:,ichem)+qRL(:,ichem))
                 qrightv(:,ichem) = half*(qRR(:,ichem)+qLR(:,ichem))
              end if
#endif
              CALL lax_friedrich_vec(qleftv,qrightv,fmean_yv,zero_flux)
              
              ! compute the final value of E including the 2D diffusive
              ! terms that ensure stability 
              emf(:,i,j,k) = Ev + (fmean_xv(:,6) - fmean_yv(:,6))
           else if (riemann2d=='hlld') then

            DO l = 1, ngrid 
               if (cbmax > 0 .and. mark_edge(l,i,j,k)) then
                  if (h_verbose >  3) write(*,'(a,3i3,a,i5)') &
                     ' HLLD: Switched to low order llf solver for edge along the '//dir(lor-1)//'-direction'// &
                     ' in edge ', i, j, k, ' on thread ', myid
                  call cmp_mag_flx_llf(qLL,qRL,qLR,qRR,emf,l,i,j,k,ilb1,ilb2,jlb1,jlb2,klb1,klb2)
               else
                  call cmp_mag_flx_hlld(qLL,qRL,qLR,qRR,emf,l,i,j,k,ilb1,ilb2,jlb1,jlb2,klb1,klb2)
               endif
            END DO
           
           else

            DO l = 1, ngrid 

               ! vx*by - vy*bx at the four edge centers
               ELL = qLL(l,3)*qLL(l,7) - qLL(l,4)*qLL(l,6)
               ERL = qRL(l,3)*qRL(l,7) - qRL(l,4)*qRL(l,6)
               ELR = qLR(l,3)*qLR(l,7) - qLR(l,4)*qLR(l,6)  
               ERR = qRR(l,3)*qRR(l,7) - qRR(l,4)*qRR(l,6) 

               if(riemann2d=='hll')then

                  ! Compute 4 fast magnetosonic velocity relative to x direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
#ifdef KROME
                  if (chemistry) qtmp(ichem) = qLL(l,ichem)  ! Effective adiabatic index
#endif
                  vLLx=qtmp(3); call find_speed_fast(qtmp,cLLx)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
#ifdef KROME
                  if (chemistry) qtmp(ichem) = qLR(l,ichem)  ! Effective adiabatic index
#endif
                  vLRx=qtmp(3); call find_speed_fast(qtmp,cLRx)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
#ifdef KROME
                  if (chemistry) qtmp(ichem) = qRL(l,ichem)  ! Effective adiabatic index
#endif
                  vRLx=qtmp(3); call find_speed_fast(qtmp,cRLx)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
#ifdef KROME
                  if (chemistry) qtmp(ichem) = qRR(l,ichem)  ! Effective adiabatic index
#endif
                  vRRx=qtmp(3); call find_speed_fast(qtmp,cRRx)

                  ! Compute 4 fast magnetosonic velocity relative to y direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
#ifdef KROME
                  if (chemistry) qtmp(ichem) = qLL(l,ichem)  ! Effective adiabatic index
#endif
                  vLLy=qtmp(3); call find_speed_fast(qtmp,cLLy)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
#ifdef KROME
                  if (chemistry) qtmp(ichem) = qLR(l,ichem)  ! Effective adiabatic index
#endif
                  vLRy=qtmp(3); call find_speed_fast(qtmp,cLRy)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
#ifdef KROME
                  if (chemistry) qtmp(ichem) = qRL(l,ichem)  ! Effective adiabatic index
#endif
                  vRLy=qtmp(3); call find_speed_fast(qtmp,cRLy)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
#ifdef KROME
                  if (chemistry) qtmp(ichem) = qRR(l,ichem)  ! Effective adiabatic index
#endif
                  vRRy=qtmp(3); call find_speed_fast(qtmp,cRRy)

                  SL=min(min(vLLx,vLRx,VRLx,vRRx)-max(cLLx,cLRx,cRLx,cRRx),zero)
                  SR=max(max(vLLx,vLRx,VRLx,vRRx)+max(cLLx,cLRx,cRLx,cRRx),zero)
                  SB=min(min(vLLy,vLRy,VRLy,vRRy)-max(cLLy,cLRy,cRLy,cRRy),zero)
                  ST=max(max(vLLy,vLRy,VRLy,vRRy)+max(cLLy,cLRy,cRLy,cRRy),zero)

                  emf(l,i,j,k) = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
                       -ST*SB/(ST-SB)*(qRR(l,6)-qLL(l,6)) &
                       +SR*SL/(SR-SL)*(qRR(l,7)-qLL(l,7))

               else if(riemann2d=='hll3r')then
                  ! This uses the 2D riemann flux described in Balsara 0911.1613
                  ! Corresponding to considering HLL3R a HLLE solver.

                  ! Compute 4 fast magnetosonic velocity relative to x direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
                  vLLx=qtmp(3); call find_speed_fast(qtmp,cLLx)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
                  vLRx=qtmp(3); call find_speed_fast(qtmp,cLRx)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
                  vRLx=qtmp(3); call find_speed_fast(qtmp,cRLx)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
                  vRRx=qtmp(3); call find_speed_fast(qtmp,cRRx)

                  ! Compute 4 fast magnetosonic velocity relative to y direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
                  vLLy=qtmp(3); call find_speed_fast(qtmp,cLLy)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
                  vLRy=qtmp(3); call find_speed_fast(qtmp,cLRy)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
                  vRLy=qtmp(3); call find_speed_fast(qtmp,cRLy)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
                  vRRy=qtmp(3); call find_speed_fast(qtmp,cRRy)

                  ! Compute max wave velocites in between the LX and RX in the x-direction 
                  qleft (1) = qLL(l,1); qright(1) = qRL(l,1) ! density
                  qleft (2) = qLL(l,2); qright(2) = qRL(l,2) ! pressure
                  qleft (3) = qLL(l,3); qright(3) = qRL(l,3) ! vt1 becomes normal velocity
                  qleft (4) = qLL(l,6); qright(4) = qRL(l,6) ! bt1 becomes normal magnetic field
                  qleft (5) = qLL(l,4); qright(5) = qRL(l,4) ! vt2 becomes transverse velocity field
                  qleft (6) = qLL(l,7); qright(6) = qRL(l,7) ! bt2 becomes transverse magnetic field 
                  qleft (7) = qLL(l,5); qright(7) = qRL(l,5) ! velocity component perp. to the plane is now transverse
                  qleft (8) = qLL(l,8); qright(8) = qRL(l,8) ! magnetic field component perp. to the plane is now transverse
                  call hll3r_speeds(qleft,qright,lmLx,lpLx)
                  qleft (1) = qLR(l,1); qright(1) = qRR(l,1) ! density
                  qleft (2) = qLR(l,2); qright(2) = qRR(l,2) ! pressure
                  qleft (3) = qLR(l,3); qright(3) = qRR(l,3) ! vt1 becomes normal velocity
                  qleft (4) = qLR(l,6); qright(4) = qRR(l,6) ! bt1 becomes normal magnetic field
                  qleft (5) = qLR(l,4); qright(5) = qRR(l,4) ! vt2 becomes transverse velocity field
                  qleft (6) = qLR(l,7); qright(6) = qRR(l,7) ! bt2 becomes transverse magnetic field 
                  qleft (7) = qLR(l,5); qright(7) = qRR(l,5) ! velocity component perp. to the plane is now transverse
                  qleft (8) = qLR(l,8); qright(8) = qRR(l,8) ! magnetic field component perp. to the plane is now transverse
                  call hll3r_speeds(qleft,qright,lmRx,lpRx)
                  ! Compute max wave velocites in between the XL and XR in the y-direction 
                  qleft (1) = qLL(l,1); qright(1) = qLR(l,1) ! density
                  qleft (2) = qLL(l,2); qright(2) = qLR(l,2) ! pressure
                  qleft (3) = qLL(l,4); qright(3) = qLR(l,4) ! vt2 becomes transverse velocity field
                  qleft (4) = qLL(l,7); qright(4) = qLR(l,7) ! bt2 becomes transverse magnetic field 
                  qleft (5) = qLL(l,3); qright(5) = qLR(l,3) ! vt1 becomes normal velocity
                  qleft (6) = qLL(l,6); qright(6) = qLR(l,6) ! bt1 becomes normal magnetic field
                  qleft (7) = qLL(l,5); qright(7) = qLR(l,5) ! velocity component perp. to the plane is now transverse
                  qleft (8) = qLL(l,8); qright(8) = qLR(l,8) ! magnetic field component perp. to the plane is now transverse
                  call hll3r_speeds(qleft,qright,lmLy,lpLy)
                  qleft (1) = qRL(l,1); qright(1) = qRR(l,1) ! density
                  qleft (2) = qRL(l,2); qright(2) = qRR(l,2) ! pressure
                  qleft (3) = qRL(l,4); qright(3) = qRR(l,4) ! vt2 becomes transverse velocity field
                  qleft (4) = qRL(l,7); qright(4) = qRR(l,7) ! bt2 becomes transverse magnetic field 
                  qleft (5) = qRL(l,3); qright(5) = qRR(l,3) ! vt1 becomes normal velocity
                  qleft (6) = qRL(l,6); qright(6) = qRR(l,6) ! bt1 becomes normal magnetic field
                  qleft (7) = qRL(l,5); qright(7) = qRR(l,5) ! velocity component perp. to the plane is now transverse
                  qleft (8) = qRL(l,8); qright(8) = qRR(l,8) ! magnetic field component perp. to the plane is now transverse
                  call hll3r_speeds(qleft,qright,lmRy,lpRy)

                  SL=min(vLLx-cLLx,vLRx-cLRx,lmLx,lmRx,zero)
                  SR=max(vRLx+cRLx,vRRx+cRRx,lpLx,lpRx,zero)
                  SB=min(vLLy-cLLy,vRLy-cRLy,lmLy,lmRy,zero)
                  ST=max(vLRy+cLRy,vRRy+cRRy,lpLy,lpRy,zero)

                  S = max(abs(SL), abs(SR), abs(SB), abs(ST))
                  alpha = 0.5_dp

                  emf(l,i,j,k) = - (( half*((ELR + ERL)*(SR*ST + SL*SB) - (ERR + ELL)*(SR*SB + SL*ST)) &
                                   - (one - alpha)*(SR*SL*(ST*(qRR(l,7)-qLR(l,7)) - SB*(qRL(l,7)*qLL(l,7))) &
                                                  - ST*SB*(SR*(qRR(l,6)-qRL(l,6)) - SL*(qLR(l,6)*qLL(l,6))) &
                                                   ) &
                                 ) / ((SR-SL)*(ST-SB)) &
                                 + alpha*S*forth*( qRR(l,7) - qLR(l,7) + qRL(l,7) - qLL(l,7) &
                                                - (qRR(l,6) - qRL(l,6) + qLR(l,6) - qLL(l,6)) ))

               else if (riemann2d=='hlla')then

                  ! Compute 4 Alfven velocity relative to x direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
                  vLLx=qtmp(3); call find_speed_alfven(qtmp,cLLx)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
                  vLRx=qtmp(3); call find_speed_alfven(qtmp,cLRx)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
                  vRLx=qtmp(3); call find_speed_alfven(qtmp,cRLx)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
                  vRRx=qtmp(3); call find_speed_alfven(qtmp,cRRx)

                  ! Compute 4 Alfven relative to y direction
                  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
                  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
                  vLLy=qtmp(3); call find_speed_alfven(qtmp,cLLy)
                  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
                  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
                  vLRy=qtmp(3); call find_speed_alfven(qtmp,cLRy)
                  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
                  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
                  vRLy=qtmp(3); call find_speed_alfven(qtmp,cRLy)
                  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
                  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
                  vRRy=qtmp(3); call find_speed_alfven(qtmp,cRRy)

                  SL=min(min(vLLx,vLRx,VRLx,vRRx)-max(cLLx,cLRx,cRLx,cRRx),zero)
                  SR=max(max(vLLx,vLRx,VRLx,vRRx)+max(cLLx,cLRx,cRLx,cRRx),zero)
                  SB=min(min(vLLy,vLRy,VRLy,vRRy)-max(cLLy,cLRy,cRLy,cRRy),zero)
                  ST=max(max(vLLy,vLRy,VRLy,vRRy)+max(cLLy,cLRy,cRLy,cRRy),zero)

                  emf(l,i,j,k) = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
                       -ST*SB/(ST-SB)*(qRR(l,6)-qLL(l,6)) &
                       +SR*SL/(SR-SL)*(qRR(l,7)-qLL(l,7))

               else

                  ! find the average value of E
                  E = forth*(ELL+ERL+ELR+ERR)
                  
                  ! call the first solver in the x direction
                  ! density
                  qleft (1) = half*(qLL(l,1)+qLR(l,1))
                  qright(1) = half*(qRR(l,1)+qRL(l,1))

                  ! pressure
                  qleft (2) = half*(qLL(l,2)+qLR(l,2))
                  qright(2) = half*(qRR(l,2)+qRL(l,2))
                  
                  ! vt1 becomes normal velocity
                  qleft (3) = half*(qLL(l,3)+qLR(l,3))
                  qright(3) = half*(qRR(l,3)+qRL(l,3))
                  
                  ! bt1 becomes normal magnetic field
                  qleft (4) = half*(qLL(l,6)+qLR(l,6))
                  qright(4) = half*(qRR(l,6)+qRL(l,6))
                  
                  ! vt2 becomes transverse velocity field
                  qleft (5) = half*(qLL(l,4)+qLR(l,4))
                  qright(5) = half*(qRR(l,4)+qRL(l,4))
                  
                  ! bt2 becomes transverse magnetic field 
                  qleft (6) = half*(qLL(l,7)+qLR(l,7))
                  qright(6) = half*(qRR(l,7)+qRL(l,7))
                  
                  ! velocity component perp. to the plane is now transverse
                  qleft (7) = half*(qLL(l,5)+qLR(l,5))
                  qright(7) = half*(qRR(l,5)+qRL(l,5))
                  
                  ! magnetic field component perp. to the plane is now transverse
                  qleft (8) = half*(qLL(l,8)+qLR(l,8))
                  qright(8) = half*(qRR(l,8)+qRL(l,8))
                  
                  zero_flux = 0.0
                  SELECT CASE (riemann2d)
                  CASE ('roe')
                     CALL athena_roe   (qleft,qright,fmean_x,zero_flux)
                  CASE ('llf')
                     CALL lax_friedrich(qleft,qright,fmean_x,zero_flux)
                  CASE ('upwind')
                     CALL upwind       (qleft,qright,fmean_x,zero_flux)
                  CASE DEFAULT
                     write(*,*)'unknown 2D riemann solver'
                     stop
                  END SELECT

                  ! call the second solver in the y direction
                  ! density
                  qleft (1) = half*(qLL(l,1)+qRL(l,1))
                  qright(1) = half*(qRR(l,1)+qLR(l,1))
                  
                  ! pressure
                  qleft (2) = half*(qLL(l,2)+qRL(l,2))
                  qright(2) = half*(qRR(l,2)+qLR(l,2))
                  
                  ! vt2 becomes normal velocity
                  qleft (3) = half*(qLL(l,4)+qRL(l,4))
                  qright(3) = half*(qRR(l,4)+qLR(l,4))
                  
                  ! bt2 becomes normal magnetic field
                  qleft (4) = half*(qLL(l,7)+qRL(l,7))
                  qright(4) = half*(qRR(l,7)+qLR(l,7))
                  
                  ! vt1 becomes transverse velocity field 
                  qleft (5) = half*(qLL(l,3)+qRL(l,3))
                  qright(5) = half*(qRR(l,3)+qLR(l,3))
                  
                  ! bt1 becomes transverse magnetic field 
                  qleft (6) = half*(qLL(l,6)+qRL(l,6))
                  qright(6) = half*(qRR(l,6)+qLR(l,6))
                  
                  ! velocity component perp. to the plane is now transverse
                  qleft (7) = half*(qLL(l,5)+qRL(l,5))
                  qright(7) = half*(qRR(l,5)+qLR(l,5))
                  
                  ! magnetic field component perp. to the plane is now transverse
                  qleft (8) = half*(qLL(l,8)+qRL(l,8))
                  qright(8) = half*(qRR(l,8)+qLR(l,8))
                  
                  zero_flux = 0.
                  SELECT CASE (riemann2d)
                  CASE ('roe')
                     CALL athena_roe   (qleft,qright,fmean_y,zero_flux)
                  CASE ('llf')
                     CALL lax_friedrich(qleft,qright,fmean_y,zero_flux)
                  CASE ('upwind')
                     CALL upwind       (qleft,qright,fmean_y,zero_flux)
                  CASE DEFAULT
                     write(*,*)'unknown 2D riemann solver'
                     stop
                  END SELECT
                  
                  ! compute the final value of E including the 2D diffusive
                  ! terms that ensure stability 
                  emf(l,i,j,k) = E + (fmean_x(6) - fmean_y(6))
                  
               endif

            END DO
           endif
        END DO
     END DO
  END DO

!  call trace_mpi('cmp_mag_flx','exit',2)
END SUBROUTINE cmp_mag_flx
!
! Compute hlld 2d fluxes
SUBROUTINE cmp_mag_flx_hlld(qLL,qRL,qLR,qRR,emf,l,i,j,k, &
                           ilb1,ilb2,jlb1,jlb2,klb1,klb2)
  USE amr_parameters,   only : dp, nvector, chemistry, debug
  USE hydro_parameters, only : nvar, ichem, smallc
  USE const,            only : half, zero
  implicit none
  INTEGER, INTENT(IN) :: l,i,j,k,ilb1,ilb2,jlb1,jlb2,klb1,klb2
  REAL(dp),DIMENSION(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2):: emf
  REAL(dp),DIMENSION(1:nvector,1:nvar),INTENT(IN)::qLL,qRL,qLR,qRR
  ! local variables
  REAL(dp) :: E,ELL,ERL,ELR,ERR,S,SL,SR,SB,ST,SAL,SAR,SAT,SAB
  REAL(dp) :: cfastLLx,cfastRLx,cfastLRx,cfastRRx,cfastLLy,cfastRLy,cfastLRy,cfastRRy
  REAL(dp) :: calfvenR,calfvenL,calfvenT,calfvenB
  REAL(dp) :: vLLx,vRLx,vLRx,vRRx,vLLy,vRLy,vLRy,vRRy
  REAL(dp) :: PtotLL,PtotLR,PtotRL,PtotRR,rcLLx,rcLRx,rcRLx,rcRRx,rcLLy,rcLRy,rcRLy,rcRRy
  REAL(dp) :: ustar,vstar,rstarLLx,rstarLRx,rstarRLx,rstarRRx,rstarLLy,rstarLRy,rstarRLy,rstarRRy
  REAL(dp) :: rstarLL,rstarLR,rstarRL,rstarRR,AstarLL,AstarLR,AstarRL,AstarRR,BstarLL,BstarLR,BstarRL,BstarRR
  REAL(dp) :: EstarLLx,EstarLRx,EstarRLx,EstarRRx,EstarLLy,EstarLRy,EstarRLy,EstarRRy,EstarLL,EstarLR,EstarRL,EstarRR
  REAL(dp) :: AstarT,AstarB,BstarR,BstarL

  REAL(dp) :: rLL,rLR,rRL,rRR,pLL,pLR,pRL,pRR,uLL,uLR,uRL,uRR,vLL,vLR,vRL,vRR, &
                       ALL,ALR,ARL,ARR,BLL,BLR,BRL,BRR,CLL,CLR,CRL,CRR
  REAL(dp),DIMENSION(1:nvar) :: qtmp

!  call trace_mpi('cmp_mag_flx_hlld','enter',2)

  rLL=qLL(l,1); pLL=qLL(l,2); uLL=qLL(l,3); vLL=qLL(l,4); ALL=qLL(l,6); BLL=qLL(l,7) ; CLL=qLL(l,8) 
  rLR=qLR(l,1); pLR=qLR(l,2); uLR=qLR(l,3); vLR=qLR(l,4); ALR=qLR(l,6); BLR=qLR(l,7) ; CLR=qLR(l,8) 
  rRL=qRL(l,1); pRL=qRL(l,2); uRL=qRL(l,3); vRL=qRL(l,4); ARL=qRL(l,6); BRL=qRL(l,7) ; CRL=qRL(l,8) 
  rRR=qRR(l,1); pRR=qRR(l,2); uRR=qRR(l,3); vRR=qRR(l,4); ARR=qRR(l,6); BRR=qRR(l,7) ; CRR=qRR(l,8) 

  ! Compute 4 fast magnetosonic velocity relative to x direction
  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
  qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
#ifdef KROME
  if (chemistry) qtmp(ichem) = qLL(l,ichem)  ! Effective adiabatic index
#endif
  call find_speed_fast(qtmp,cfastLLx)
  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
  qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
#ifdef KROME
  if (chemistry) qtmp(ichem) = qLR(l,ichem)  ! Effective adiabatic index
#endif
  call find_speed_fast(qtmp,cfastLRx)
  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
  qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
#ifdef KROME
  if (chemistry) qtmp(ichem) = qRL(l,ichem)  ! Effective adiabatic index
#endif
  call find_speed_fast(qtmp,cfastRLx)
  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
  qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
#ifdef KROME
  if (chemistry) qtmp(ichem) = qRR(l,ichem)  ! Effective adiabatic index
#endif
  call find_speed_fast(qtmp,cfastRRx)

  ! Compute 4 fast magnetosonic velocity relative to y direction
  qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
  qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
#ifdef KROME
  if (chemistry) qtmp(ichem) = qLL(l,ichem)  ! Effective adiabatic index
#endif
  call find_speed_fast(qtmp,cfastLLy)
  qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
  qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
#ifdef KROME
  if (chemistry) qtmp(ichem) = qLR(l,ichem)  ! Effective adiabatic index
#endif
  call find_speed_fast(qtmp,cfastLRy)
  qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
  qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
#ifdef KROME
  if (chemistry) qtmp(ichem) = qRL(l,ichem)  ! Effective adiabatic index
#endif
  call find_speed_fast(qtmp,cfastRLy)
  qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
  qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
#ifdef KROME
  if (chemistry) qtmp(ichem) = qRR(l,ichem)  ! Effective adiabatic index
#endif
  call find_speed_fast(qtmp,cfastRRy)

  SL=min(uLL,uLR,uRL,uRR)-max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
  SR=max(uLL,uLR,uRL,uRR)+max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
  SB=min(vLL,vLR,vRL,vRR)-max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
  ST=max(vLL,vLR,vRL,vRR)+max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)

  if(SB>0d0)then
     if(SL>0d0)then
        ELL=uLL*BLL-vLL*ALL

        E=ELL
     else if(SR<0d0)then
        ERL=uRL*BRL-vRL*ARL

        E=ERL
     else
        PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
        PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
        PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
        PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
        
        rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL) 
        rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
        rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR) 
        rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
        
        ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
        vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
        
        rstarLLx=rLL*(SL-uLL)/(SL-ustar)
        rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
        AstarLL=ALL*(SB-vLL)/(SB-vstar)
      
        rstarLRx=rLR*(SL-uLR)/(SL-ustar)
        rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
        AstarLR=ALR*(ST-vLR)/(ST-vstar)
      
        rstarRLx=rRL*(SR-uRL)/(SR-ustar)
        rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
        AstarRL=ARL*(SB-vRL)/(SB-vstar)        
      
        rstarRRx=rRR*(SR-uRR)/(SR-ustar)
        rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
        AstarRR=ARR*(ST-vRR)/(ST-vstar)

        calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                     abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
        calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                     abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)

        SAL=min(ustar-calfvenL,zero)
        SAR=max(ustar+calfvenR,zero)

        BstarLL=BLL*(SL-uLL)/(SL-ustar)
        BstarRL=BRL*(SR-uRL)/(SR-ustar)

        EstarLLx=ustar*BstarLL-vLL*ALL
        EstarRLx=ustar*BstarRL-vRL*ARL
        
        E=(SAR*EstarLLx-SAL*EstarRLx+SAR*SAL*(BRL-BLL))/(SAR-SAL)
     endif
  else if (ST<0d0)then
     if(SL>0d0)then
        ELR=uLR*BLR-vLR*ALR

        E=ELR
     else if(SR<0d0)then
        ERR=uRR*BRR-vRR*ARR

        E=ERR
     else
        PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
        PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
        PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
        PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
        
        rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL) 
        rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
        rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR) 
        rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
        
        ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
        vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
        
        rstarLLx=rLL*(SL-uLL)/(SL-ustar)
        rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
        AstarLL=ALL*(SB-vLL)/(SB-vstar)
      
        rstarLRx=rLR*(SL-uLR)/(SL-ustar)
        rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
        AstarLR=ALR*(ST-vLR)/(ST-vstar)
      
        rstarRLx=rRL*(SR-uRL)/(SR-ustar)
        rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
        AstarRL=ARL*(SB-vRL)/(SB-vstar)        
      
        rstarRRx=rRR*(SR-uRR)/(SR-ustar)
        rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
        AstarRR=ARR*(ST-vRR)/(ST-vstar)

        calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                     abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
        calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                     abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)

        SAL=min(ustar-calfvenL,zero)
        SAR=max(ustar+calfvenR,zero)

        BstarLR=BLR*(SL-uLR)/(SL-ustar)
        BstarRR=BRR*(SR-uRR)/(SR-ustar)

        EstarLRx=ustar*BstarLR-vLR*ALR
        EstarRRx=ustar*BstarRR-vRR*ARR

        E=(SAR*EstarLRx-SAL*EstarRRx+SAR*SAL*(BRR-BLR))/(SAR-SAL)
     endif
  else if(SL>0d0)then
     PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
     PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
     PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
     PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
      
     rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL) 
     rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
     rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR) 
     rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
      
     ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
     vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
      
     rstarLLy=rLL*(SB-vLL)/(SB-vstar)
     rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
     BstarLL=BLL*(SL-uLL)/(SL-ustar)
      
     rstarLRy=rLR*(ST-vLR)/(ST-vstar)
     rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
     BstarLR=BLR*(SL-uLR)/(SL-ustar)
      
     rstarRLy=rRL*(SB-vRL)/(SB-vstar)
     rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
     BstarRL=BRL*(SR-uRL)/(SR-ustar)      
      
     rstarRRy=rRR*(ST-vRR)/(ST-vstar)
     rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
     BstarRR=BRR*(SR-uRR)/(SR-ustar)
    
     calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                  abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
     calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                  abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)

     SAB=min(vstar-calfvenB,zero)
     SAT=max(vstar+calfvenT,zero)

     AstarLL=ALL*(SB-vLL)/(SB-vstar)
     AstarLR=ALR*(ST-vLR)/(ST-vstar)
  
     EstarLLy=uLL*BLL-vstar*AstarLL
     EstarLRy=uLR*BLR-vstar*AstarLR

     E=(SAT*EstarLLy-SAB*EstarLRy-SAT*SAB*(ALR-ALL))/(SAT-SAB)
  else if(SR<0d0)then
     PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
     PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
     PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
     PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
      
     rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL) 
     rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
     rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR) 
     rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
      
     ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
     vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
      
     rstarLLy=rLL*(SB-vLL)/(SB-vstar)
     rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
     BstarLL=BLL*(SL-uLL)/(SL-ustar)
      
     rstarLRy=rLR*(ST-vLR)/(ST-vstar)
     rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
     BstarLR=BLR*(SL-uLR)/(SL-ustar)
      
     rstarRLy=rRL*(SB-vRL)/(SB-vstar)
     rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
     BstarRL=BRL*(SR-uRL)/(SR-ustar)      
      
     rstarRRy=rRR*(ST-vRR)/(ST-vstar)
     rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
     BstarRR=BRR*(SR-uRR)/(SR-ustar)
 
     calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                  abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
     calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                  abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)

     SAB=min(vstar-calfvenB,zero)
     SAT=max(vstar+calfvenT,zero)

     AstarRL=ARL*(SB-vRL)/(SB-vstar)
     AstarRR=ARR*(ST-vRR)/(ST-vstar)

     EstarRLy=uRL*BRL-vstar*AstarRL
     EstarRRy=uRR*BRR-vstar*AstarRR

     E=(SAT*EstarRLy-SAB*EstarRRy-SAT*SAB*(ARR-ARL))/(SAT-SAB)
  else
     PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
     PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
     PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
     PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
      
     rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL) 
     rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
     rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR) 
     rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
      
     ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
     vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
      
     rstarLLx=rLL*(SL-uLL)/(SL-ustar); BstarLL=BLL*(SL-uLL)/(SL-ustar)
     rstarLLy=rLL*(SB-vLL)/(SB-vstar); AstarLL=ALL*(SB-vLL)/(SB-vstar)
     rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
      
     rstarLRx=rLR*(SL-uLR)/(SL-ustar); BstarLR=BLR*(SL-uLR)/(SL-ustar)
     rstarLRy=rLR*(ST-vLR)/(ST-vstar); AstarLR=ALR*(ST-vLR)/(ST-vstar)
     rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
      
     rstarRLx=rRL*(SR-uRL)/(SR-ustar); BstarRL=BRL*(SR-uRL)/(SR-ustar)
     rstarRLy=rRL*(SB-vRL)/(SB-vstar); AstarRL=ARL*(SB-vRL)/(SB-vstar)
     rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
      
     rstarRRx=rRR*(SR-uRR)/(SR-ustar); BstarRR=BRR*(SR-uRR)/(SR-ustar)
     rstarRRy=rRR*(ST-vRR)/(ST-vstar); AstarRR=ARR*(ST-vRR)/(ST-vstar)
     rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)

     calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                  abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
     calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                  abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
     calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                  abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
     calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                  abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
     SAL=min(ustar-calfvenL,zero)
     SAR=max(ustar+calfvenR,zero)
     SAB=min(vstar-calfvenB,zero)
     SAT=max(vstar+calfvenT,zero)
     EstarLL =ustar*BstarLL-vstar*AstarLL
     EstarRL =ustar*BstarRL-vstar*AstarRL
     EstarLR =ustar*BstarLR-vstar*AstarLR
     EstarRR =ustar*BstarRR-vstar*AstarRR

     AstarT=(SAR*AstarRR-SAL*AstarLR)/(SAR-SAL); AstarB=(SAR*AstarRL-SAL*AstarLL)/(SAR-SAL)
     BstarR=(SAT*BstarRR-SAB*BstarRL)/(SAT-SAB); BstarL=(SAT*BstarLR-SAB*BstarLL)/(SAT-SAB)

     E=(SAL*SAB*EstarRR-SAL*SAT*EstarRL-SAR*SAB*EstarLR+SAR*SAT*EstarLL)/(SAR-SAL)/(SAT-SAB) &
          & -SAT*SAB/(SAT-SAB)*(AstarT-AstarB)+SAR*SAL/(SAR-SAL)*(BstarR-BstarL)
  endif

  if (debug) then
  if (abs(E) > 1e12) then
    print 1, 'S  B,  T,  L,  R      :', SB, ST, SL, SR
    print 1, 'S AB, AT, AL, AR      :', SAB, SAT, SAL, SAR
    print 1, 'E LL, LR, RL, RR      :', ELL, ELR, ERL, ERR
    print 1, 'A LL, LR, RL, RR      :', ALL, ALR, ARL, ARR
    print 1, 'B LL, LR, RL, RR      :', BLL, BLR, BRL, BRR
    print 1, 'Estar  LL, LR, RL, RR :', EstarLL, EstarLR, EstarRL, EstarRR 
    print 1, 'Estarx LL, LR, RL, RR :', EstarLLx, EstarLRx, EstarRLx, EstarRRx 
    print 1, 'Estary LL, LR, RL, RR :', EstarLLy, EstarLRy, EstarRLy, EstarRRy 
    print 1, 'Astar  T, B           :', AstarT, AstarB
    print 1, 'Bstar  L, R           :', BstarL, BstarR
    print 1, 'ustar, vstar          :', ustar, vstar  
    print 1, 'calfven, B, T, L, R   :', calfvenB, calfvenT, calfvenL, calfvenR
    print 1, 'rstarx LL, LR, RL, RR :', rstarLLx, rstarLRx, rstarRLx, rstarRRx 
    print 1, 'rstary LL, LR, RL, RR :', rstarLLy, rstarLRy, rstarRLy, rstarRRy 
    print 1, 'qLL                   :', qLL(l,:)
    print 1, 'qLR                   :', qLR(l,:)
    print 1, 'qRL                   :', qRL(l,:)
    print 1, 'qRR                   :', qRR(l,:)
    print 1, 'E ', E
  1 format(1x,a,1p,13e14.5)
  endif
  endif
  
  emf(l,i,j,k) = E

!  call trace_mpi('cmp_mag_flx_hlld','exit',2)
END SUBROUTINE cmp_mag_flx_hlld
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! Compute llf 2d fluxes
SUBROUTINE cmp_mag_flx_llf(qLL,qRL,qLR,qRR,emf,l,i,j,k, &
                           ilb1,ilb2,jlb1,jlb2,klb1,klb2)
  USE amr_parameters,   only : dp, nvector, chemistry
  USE hydro_parameters, only : nvar, ichem
  USE const,            only : half,forth
  implicit none
  INTEGER, INTENT(IN) :: l,i,j,k,ilb1,ilb2,jlb1,jlb2,klb1,klb2
  REAL(dp),DIMENSION(1:nvector,ilb1:ilb2,jlb1:jlb2,klb1:klb2):: emf
  REAL(dp),DIMENSION(1:nvector,1:nvar),INTENT(IN)::qLL,qRL,qLR,qRR
  ! local variables
  REAL(dp)                            ::ELL,ERL,ELR,ERR,Ev,zero_flux
  REAL(dp),DIMENSION(1:nvar), target  ::qleft,qright,fmean_x,fmean_y
  
!  call trace_mpi('cmp_mag_flx_llf','enter',2)

  ! vx*by - vy*bx at the four edge centers
  ELL = qLL(l,3)*qLL(l,7) - qLL(l,4)*qLL(l,6)
  ERL = qRL(l,3)*qRL(l,7) - qRL(l,4)*qRL(l,6)
  ELR = qLR(l,3)*qLR(l,7) - qLR(l,4)*qLR(l,6)  
  ERR = qRR(l,3)*qRR(l,7) - qRR(l,4)*qRR(l,6) 

  ! find the average value of E
  Ev = forth*(ELL+ERL+ELR+ERR)
    
  ! call the first solver in the x direction
  ! density
  qleft (1) = half*(qLL(l,1)+qLR(l,1))
  qright(1) = half*(qRR(l,1)+qRL(l,1))
  
  ! pressure
  qleft (2) = half*(qLL(l,2)+qLR(l,2))
  qright(2) = half*(qRR(l,2)+qRL(l,2))
  
  ! vt1 becomes normal velocity
  qleft (3) = half*(qLL(l,3)+qLR(l,3))
  qright(3) = half*(qRR(l,3)+qRL(l,3))
  
  ! bt1 becomes normal magnetic field
  qleft (4) = half*(qLL(l,6)+qLR(l,6))
  qright(4) = half*(qRR(l,6)+qRL(l,6))
  
  ! vt2 becomes transverse velocity field
  qleft (5) = half*(qLL(l,4)+qLR(l,4))
  qright(5) = half*(qRR(l,4)+qRL(l,4))
  
  ! bt2 becomes transverse magnetic field 
  qleft (6) = half*(qLL(l,7)+qLR(l,7))
  qright(6) = half*(qRR(l,7)+qRL(l,7))
  
  ! velocity component perp. to the plane is now transverse
  qleft (7) = half*(qLL(l,5)+qLR(l,5))
  qright(7) = half*(qRR(l,5)+qRL(l,5))
  
  ! magnetic field component perp. to the plane is now transverse
  qleft (8) = half*(qLL(l,8)+qLR(l,8))
  qright(8) = half*(qRR(l,8)+qRL(l,8))
  zero_flux = 0.0
#ifdef KROME
  ! Effective adiabatic index
  if (chemistry) then
     qleft (ichem) = half*(qLL(l,ichem)+qLR(l,ichem))
     qright(ichem) = half*(qRR(l,ichem)+qRL(l,ichem))
  end if
#endif
  CALL lax_friedrich(qleft,qright,fmean_x,zero_flux)

  ! call the second solver in the y direction
  ! density
  qleft (1) = half*(qLL(l,1)+qRL(l,1))
  qright(1) = half*(qRR(l,1)+qLR(l,1))

  ! pressure
  qleft (2) = half*(qLL(l,2)+qRL(l,2))
  qright(2) = half*(qRR(l,2)+qLR(l,2))

  ! vt2 becomes normal velocity
  qleft (3) = half*(qLL(l,4)+qRL(l,4))
  qright(3) = half*(qRR(l,4)+qLR(l,4))
  
  ! bt2 becomes normal magnetic field
  qleft (4) = half*(qLL(l,7)+qRL(l,7))
  qright(4) = half*(qRR(l,7)+qLR(l,7))
  
  ! vt1 becomes transverse velocity field 
  qleft (5) = half*(qLL(l,3)+qRL(l,3))
  qright(5) = half*(qRR(l,3)+qLR(l,3))
  
  ! bt1 becomes transverse magnetic field 
  qleft (6) = half*(qLL(l,6)+qRL(l,6))
  qright(6) = half*(qRR(l,6)+qLR(l,6))
  
  ! velocity component perp. to the plane is now transverse
  qleft (7) = half*(qLL(l,5)+qRL(l,5))
  qright(7) = half*(qRR(l,5)+qLR(l,5))
  
  ! magnetic field component perp. to the plane is now transverse
  qleft (8) = half*(qLL(l,8)+qRL(l,8))
  qright(8) = half*(qRR(l,8)+qLR(l,8))
  
  zero_flux = 0.0_dp
#ifdef KROME
  ! Effective adiabatic index
  if (chemistry) then
     qleft (ichem) = half*(qLL(l,ichem)+qRL(l,ichem))
     qright(ichem) = half*(qRR(l,ichem)+qLR(l,ichem))
  end if
#endif
  CALL lax_friedrich(qleft,qright,fmean_y,zero_flux)
  
  ! compute the final value of E including the 2D diffusive
  ! terms that ensure stability 
  emf(l,i,j,k) = Ev + (fmean_x(6) - fmean_y(6))
!  call trace_mpi('cmp_mag_flx_llf','exit',2)
END SUBROUTINE cmp_mag_flx_llf
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,bf,gravin,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  use hydro_commons, only : Grho, ichem
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+1)::q  
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  

  integer ::i, j, k, l, n, idim
  real(dp)::eint, smalle, smallp, etot
  real(dp),dimension(1:nvector)::eken,emag

!  call trace_mpi('ctoprim','enter',2)

  smalle = smallc**2/gamma/(gamma-one)
  smallp = smallr*smallc**2/gamma

  ! Store face centered magnetic field
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2+1
           if(i<=iu2)then
              DO l = 1, ngrid
                 bf(l,i,j,k,1) = uin(l,i,j,k,6)
              END DO
           else
              DO l = 1, ngrid
                 bf(l,i,j,k,1) = uin(l,i-1,j,k,nvar+1)
              END DO
           endif
        end do
     end do
  end do
  do k = ku1, ku2
     do j = ju1, ju2+1
        do i = iu1, iu2
           if(j<=ju2)then
              DO l = 1, ngrid
                 bf(l,i,j,k,2) = uin(l,i,j,k,7)
              END DO
           else
              DO l = 1, ngrid
                 bf(l,i,j,k,2) = uin(l,i,j-1,k,nvar+2)
              END DO
           endif
        end do
     end do
  end do
  do k = ku1, ku2+1
     do j = ju1, ju2
        do i = iu1, iu2
           if(k<=ku2)then
              DO l = 1, ngrid
                 bf(l,i,j,k,3) = uin(l,i,j,k,8)
              END DO
           else
              DO l = 1, ngrid
                 bf(l,i,j,k,3) = uin(l,i,j,k-1,nvar+3)
              END DO
           endif
        end do
     end do
  end do
  if (ngrid<nvector) bf(ngrid+1:nvector,:,:,:,:) = 0.0

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2

           ! Compute density
           do l = 1, ngrid
              q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)
           end do
           do l = ngrid+1, nvector
              q(l,i,j,k,1) = smallr
           end do

           ! Compute velocities
           do l = 1, ngrid
              q(l,i,j,k,2) = uin(l,i,j,k,2)/uin(l,i,j,k,1)
              q(l,i,j,k,3) = uin(l,i,j,k,3)/uin(l,i,j,k,1)
              q(l,i,j,k,4) = uin(l,i,j,k,4)/uin(l,i,j,k,1)
           end do
           do l = ngrid+1, nvector
              q(l,i,j,k,2) = 0.0
              q(l,i,j,k,3) = 0.0
              q(l,i,j,k,4) = 0.0
           end do

           ! Compute cell centered magnetic field
           DO l = 1, ngrid
              q(l,i,j,k,6) = (uin(l,i,j,k,6)+uin(l,i,j,k,nvar+1))*half
              q(l,i,j,k,7) = (uin(l,i,j,k,7)+uin(l,i,j,k,nvar+2))*half
              q(l,i,j,k,8) = (uin(l,i,j,k,8)+uin(l,i,j,k,nvar+3))*half
           END DO
           do l = ngrid+1, nvector
              q(l,i,j,k,6) = 0.0
              q(l,i,j,k,7) = 0.0
              q(l,i,j,k,8) = 0.0
           end do

           ! Compute specific kinetic energy and magnetic energy
           do l = 1, ngrid
              eken(l) = half*(q(l,i,j,k,2)**2+q(l,i,j,k,3)**2+q(l,i,j,k,4)**2)
              emag(l) = half*(q(l,i,j,k,6)**2+q(l,i,j,k,7)**2+q(l,i,j,k,8)**2)
           end do

           ! Compute thermal pressure through EOS
#ifdef KROME
           do l = 1, ngrid
              etot = uin(l,i,j,k,5) - emag(l)
              eint = etot/uin(l,i,j,k,1)-eken(l)
              if (chemistry) then
                 q(l,i,j,k,5)=MAX((uin(l,i,j,k,ichem)-one)*q(l,i,j,k,1)*eint,smallp)
              else
                 q(l,i,j,k,5)=MAX((gamma-one)*q(l,i,j,k,1)*eint,smallp)
              endif
              q(l,i,j,k,nvar+1) = log(q(l,i,j,k,5)) - uin(l,i,j,k,ichem)*log(q(l,i,j,k,1))
           end do
#else
           do l = 1, ngrid
              etot = uin(l,i,j,k,5) - emag(l)
              eint = etot/uin(l,i,j,k,1)-eken(l)
              q(l,i,j,k,5)=MAX((gamma-one)*q(l,i,j,k,1)*eint,smallp)
              q(l,i,j,k,nvar+1) = log(q(l,i,j,k,5)) - gamma*log(q(l,i,j,k,1))
           end do
#endif
           do l = ngrid+1, nvector
              q(l,i,j,k,5) = smallp
              q(l,i,j,k,nvar+1) = 0.0
           end do

           ! Gravity predictor step
           do idim = 1, ndim
              do l = 1, ngrid
                 q(l,i,j,k,idim+1) = q(l,i,j,k,idim+1) + Grho*gravin(l,i,j,k,idim)*dt*half
              end do
              do l = ngrid+1, nvector
                 q(l,i,j,k,idim+1) = 0.0
              end do
           end do

        end do
     end do
  end do

  ! Passive scalar
  do n = 9, nvar
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 q(l,i,j,k,n) = uin(l,i,j,k,n)/uin(l,i,j,k,1)
#ifdef KROME
                 if (chemistry .and. n==ichem) q(l,i,j,k,n) = uin(l,i,j,k,n) ! adiabatic index should just be copied
#endif
              end do
              do l = ngrid+1, nvector
                 q(l,i,j,k,n) = smallr
              end do
           end do
        end do
     end do
  end do
 
!  call trace_mpi('ctoprim','exit',2)
end subroutine ctoprim
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope(bf,q,dq,dbf,dx,dt,ngrid,rslope_type)
  use amr_parameters
  use hydro_parameters, only : iu1,iu2,ju1,ju2,ku1,ku2,nvar,ichem, &
                               limit_slope,gamma,slope_type_global=>slope_type
  use const
  implicit none

  integer::ngrid, slope_type
  real(dp)::dx,dt,rslope_type
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+1)::q 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:ndim)::dbf

  ! local arrays
  integer::i, j, k, l, n, idir
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop, alphav, malphav, palphav
  real(dp), dimension(1:nvector)::vdsgn, vdlim, vdcen, vdlft, vdrgt, vslop, vlimd, vlimp
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff,idff,xslope_type
  integer::ilo,ihi,jlo,jhi,klo,khi
  
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  slope_type=rslope_type

  if(slope_type==0)then
     dq=zero
     dbf=zero
     return
  end if
!  call trace_mpi('uslope','enter',2)

#if NDIM==1
  do n = 1, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              if(slope_type==1.or.slope_type==2)then  ! minmod or average==MC-limiter
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else
                 write(*,*)'Unknown slope type'
                 stop
              end if
           end do
        end do
     end do     
  end do
#endif

#if NDIM==2              
  if(slope_type==1.or.slope_type==2)then  ! minmod or average==MC-limiter
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 2d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dfll = q(l,i-1,j-1,k,n)-q(l,i,j,k,n)
                    dflm = q(l,i-1,j  ,k,n)-q(l,i,j,k,n)
                    dflr = q(l,i-1,j+1,k,n)-q(l,i,j,k,n)
                    dfml = q(l,i  ,j-1,k,n)-q(l,i,j,k,n)
                    dfmm = q(l,i  ,j  ,k,n)-q(l,i,j,k,n)
                    dfmr = q(l,i  ,j+1,k,n)-q(l,i,j,k,n)
                    dfrl = q(l,i+1,j-1,k,n)-q(l,i,j,k,n)
                    dfrm = q(l,i+1,j  ,k,n)-q(l,i,j,k,n)
                    dfrr = q(l,i+1,j+1,k,n)-q(l,i,j,k,n)
                    
                    vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    
                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dff  = half*(abs(dfx)+abs(dfy))
                    
                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif
                    
                    dlim = slop
                    
                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy

                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif
  ! 1D transverse TVD slopes for face-centered magnetic fields
  ! Bx along direction Y 
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi+1 ! WARNING HERE
           do l = 1, ngrid
              dlft = slope_type*(bf(l,i,j  ,k,1) - bf(l,i,j-1,k,1))
              drgt = slope_type*(bf(l,i,j+1,k,1) - bf(l,i,j  ,k,1))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(l,i,j,k,1,1) = dsgn*min(dlim,abs(dcen))
           end do
        enddo
     end do
  end do
  ! By along direction X
  do k = klo, khi
     do j = jlo, jhi+1 ! WARNING HERE
        do i = ilo, ihi
           do l = 1, ngrid
              dlft = slope_type*(bf(l,i  ,j,k,2) - bf(l,i-1,j,k,2))
              drgt = slope_type*(bf(l,i+1,j,k,2) - bf(l,i  ,j,k,2))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dbf(l,i,j,k,2,1) = dsgn*min(dlim,abs(dcen))
           end do
        enddo
     end do
  end do
#endif

#if NDIM==3
  if(slope_type==1.or.slope_type==2)then  ! minmod or average==MC-limiter
     xslope_type = slope_type
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 vdlft = (q(:,i  ,j,k,n) - q(:,i-1,j,k,n))        ! = dq
                 vdrgt = (q(:,i+1,j,k,n) - q(:,i  ,j,k,n))        ! = dq
                 vdcen = half*(vdlft+vdrgt)                       ! = dq
                 vdsgn = sign(one, vdcen)                         ! sign(slope)
                 vslop = xslope_type*min(abs(vdlft),abs(vdrgt))   ! slope_type=1: vslope=dq, slope_type=2: vslope=2*dq 
                 vslop = merge(vslop,zero,(vdlft*vdrgt)>zero)     ! if both are non-zero, use the smaller one, or times 2
                 dq(:,i,j,k,n,1) = vdsgn*min(vslop,abs(vdcen))    ! the smaller of the twice the smaller slope, and the central one
                 ! slopes in second coordinate direction
                 vdlft = (q(:,i,j  ,k,n) - q(:,i,j-1,k,n))
                 vdrgt = (q(:,i,j+1,k,n) - q(:,i,j  ,k,n))
                 vdcen = half*(vdlft+vdrgt)
                 vdsgn = sign(one, vdcen)
                 vslop = xslope_type*min(abs(vdlft),abs(vdrgt))
                 vslop = merge(vslop,zero,(vdlft*vdrgt)>zero)
                 dq(:,i,j,k,n,2) = vdsgn*min(vslop,abs(vdcen))
                 ! slopes in third coordinate direction
                 vdlft = (q(:,i,j,k  ,n) - q(:,i,j,k-1,n))
                 vdrgt = (q(:,i,j,k+1,n) - q(:,i,j,k  ,n))
                 vdcen = half*(vdlft+vdrgt)
                 vdsgn = sign(one, vdcen)
                 vslop = xslope_type*min(abs(vdlft),abs(vdrgt))
                 vslop = merge(vslop,zero,(vdlft*vdrgt)>zero)
                 dq(:,i,j,k,n,3) = vdsgn*min(vslop,abs(vdcen))
              end do
           end do
        end do
     end do
  else if(slope_type==5.or.slope_type==6)then ! van Albada 1
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 vdlft = (q(:,i  ,j,k,n) - q(:,i-1,j,k,n))
                 vdrgt = (q(:,i+1,j,k,n) - q(:,i  ,j,k,n))
                 vslop = (vdlft**2+vdlft*vdrgt)/(vdlft**2+vdrgt**2+1e-30)
                 vslop = vdlft*(1.-vslop)+vdrgt*vslop
                 dq(:,i,j,k,n,1) = merge(vslop,zero,(vdlft*vdrgt)>zero)
                 ! slopes in second coordinate direction
                 vdlft = (q(:,i,j  ,k,n) - q(:,i,j-1,k,n))
                 vdrgt = (q(:,i,j+1,k,n) - q(:,i,j  ,k,n))
                 vslop = (vdlft**2+vdlft*vdrgt)/(vdlft**2+vdrgt**2+1e-30)
                 vslop = vdlft*(1.-vslop)+vdrgt*vslop
                 dq(:,i,j,k,n,2) = merge(vslop,zero,(vdlft*vdrgt)>zero)
                 ! slopes in third coordinate direction
                 vdlft = (q(:,i,j,k  ,n) - q(:,i,j,k-1,n))
                 vdrgt = (q(:,i,j,k+1,n) - q(:,i,j,k  ,n))
                 vslop = (vdlft**2+vdlft*vdrgt)/(vdlft**2+vdrgt**2+1e-30)
                 vslop = vdlft*(1.-vslop)+vdrgt*vslop
                 dq(:,i,j,k,n,3) = merge(vslop,zero,(vdlft*vdrgt)>zero)
              end do
           end do
        end do
     end do
  else if(slope_type==3.or.slope_type==4)then ! positivity preserving 3d unsplit slope
     xslope_type = real(rslope_type - 2_dp,kind=dp)
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dflll = q(l,i-1,j-1,k-1,n)-q(l,i,j,k,n)
                    dflml = q(l,i-1,j  ,k-1,n)-q(l,i,j,k,n)
                    dflrl = q(l,i-1,j+1,k-1,n)-q(l,i,j,k,n)
                    dfmll = q(l,i  ,j-1,k-1,n)-q(l,i,j,k,n)
                    dfmml = q(l,i  ,j  ,k-1,n)-q(l,i,j,k,n)
                    dfmrl = q(l,i  ,j+1,k-1,n)-q(l,i,j,k,n)
                    dfrll = q(l,i+1,j-1,k-1,n)-q(l,i,j,k,n)
                    dfrml = q(l,i+1,j  ,k-1,n)-q(l,i,j,k,n)
                    dfrrl = q(l,i+1,j+1,k-1,n)-q(l,i,j,k,n)

                    dfllm = q(l,i-1,j-1,k  ,n)-q(l,i,j,k,n)
                    dflmm = q(l,i-1,j  ,k  ,n)-q(l,i,j,k,n)
                    dflrm = q(l,i-1,j+1,k  ,n)-q(l,i,j,k,n)
                    dfmlm = q(l,i  ,j-1,k  ,n)-q(l,i,j,k,n)
                    dfmmm = q(l,i  ,j  ,k  ,n)-q(l,i,j,k,n)
                    dfmrm = q(l,i  ,j+1,k  ,n)-q(l,i,j,k,n)
                    dfrlm = q(l,i+1,j-1,k  ,n)-q(l,i,j,k,n)
                    dfrmm = q(l,i+1,j  ,k  ,n)-q(l,i,j,k,n)
                    dfrrm = q(l,i+1,j+1,k  ,n)-q(l,i,j,k,n)

                    dfllr = q(l,i-1,j-1,k+1,n)-q(l,i,j,k,n)
                    dflmr = q(l,i-1,j  ,k+1,n)-q(l,i,j,k,n)
                    dflrr = q(l,i-1,j+1,k+1,n)-q(l,i,j,k,n)
                    dfmlr = q(l,i  ,j-1,k+1,n)-q(l,i,j,k,n)
                    dfmmr = q(l,i  ,j  ,k+1,n)-q(l,i,j,k,n)
                    dfmrr = q(l,i  ,j+1,k+1,n)-q(l,i,j,k,n)
                    dfrlr = q(l,i+1,j-1,k+1,n)-q(l,i,j,k,n)
                    dfrmr = q(l,i+1,j  ,k+1,n)-q(l,i,j,k,n)
                    dfrrr = q(l,i+1,j+1,k+1,n)-q(l,i,j,k,n)

                    vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                    vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)

                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dfz  = half*(q(l,i,j,k+1,n)-q(l,i,j,k-1,n))
                    idff  = xslope_type / sqrt(dfx**2+dfy**2+dfz**2+tiny(1.0_dp)) ! xslope_type = 1 or 2 dep. on slope_type

                    ! vmin and vmax are the lower and upper limit of one-cell differences
                    ! so a constant slope with step dq will give vmin=-dq and vmax=+dq,
                    ! while dfx (dff for any direction) will also be = dq, unless we keep
                    ! the 'half' in front of the sqrt.  With that kept, the dlim factor
                    ! will remain unity until the slopes differ with a factor of 2. Then,
                    ! if any of the 26 slopes is steeper than twice the gradient slope, 
                    ! the gradient will be diminished.
                    
                    ! idff > 0, so 1/dff = idff < 0.001 * huge is equivalent to dff > 0
                    slop = merge(min(one,min(abs(vmin),abs(vmax))*idff), one, idff<0.001_dp * huge(1.0_dp))

                    dlim = slop

                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy
                    dq(l,i,j,k,n,3) = dlim*dfz

                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type'
     stop
  endif     
  ! 2D transverse TVD slopes for face-centered magnetic fields
  ! Bx along direction Y and Z
  xslope_type=MIN(slope_type,2)
  if(slope_type>=1.and.slope_type<=5)then  ! minmod or average==MC-limiter
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi+1 ! WARNING HERE
              ! slopes in first coordinate direction
              vdlft = (bf(:,i,j  ,k,1) - bf(:,i,j-1,k,1))
              vdrgt = (bf(:,i,j+1,k,1) - bf(:,i,j  ,k,1))
              vdcen = half*(vdlft+vdrgt)
              vdsgn = sign(one, vdcen)
              vslop = xslope_type*min(abs(vdlft),abs(vdrgt))
              vslop = merge(vslop,zero,(vdlft*vdrgt)>zero)
              dbf(:,i,j,k,1,1) = vdsgn*min(vslop,abs(vdcen))
              ! slopes in second coordinate direction
              vdlft = (bf(:,i,j,k  ,1) - bf(:,i,j,k-1,1))
              vdrgt = (bf(:,i,j,k+1,1) - bf(:,i,j,k  ,1))
              vdcen = half*(vdlft+vdrgt)
              vdsgn = sign(one, vdcen)
              vslop = xslope_type*min(abs(vdlft),abs(vdrgt))
              vslop = merge(vslop,zero,(vdlft*vdrgt)>zero)
              dbf(:,i,j,k,1,2) = vdsgn*min(vslop,abs(vdcen))
           end do
        end do
     end do
  elseif (slope_type==6) then
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi+1 ! WARNING HERE
              ! slopes in first coordinate direction
              vdlft = (bf(:,i,j  ,k,1) - bf(:,i,j-1,k,1))
              vdrgt = (bf(:,i,j+1,k,1) - bf(:,i,j  ,k,1))
              vslop = (vdlft**2+vdlft*vdrgt)/(vdlft**2+vdrgt**2+1e-30)
              vslop = vdlft*(1.-vslop)+vdrgt*vslop
              dbf(:,i,j,k,1,1) = merge(vslop,zero,(vdlft*vdrgt)>zero)
              ! slopes in second coordinate direction
              vdlft = (bf(:,i,j,k  ,1) - bf(:,i,j,k-1,1))
              vdrgt = (bf(:,i,j,k+1,1) - bf(:,i,j,k  ,1))
              vslop = (vdlft**2+vdlft*vdrgt)/(vdlft**2+vdrgt**2+1e-30)
              vslop = vdlft*(1.-vslop)+vdrgt*vslop
              dbf(:,i,j,k,1,2) = merge(vslop,zero,(vdlft*vdrgt)>zero)
           end do
        end do
     end do
  endif
  ! By along direction X and Z
  if(slope_type>=1.and.slope_type<=5)then  ! minmod or average==MC-limiter
     do k = klo, khi
        do j = jlo, jhi+1 ! WARNING HERE
           do i = ilo, ihi
              ! slopes in first coordinate direction
              vdlft = (bf(:,i  ,j,k,2) - bf(:,i-1,j,k,2))
              vdrgt = (bf(:,i+1,j,k,2) - bf(:,i  ,j,k,2))
              vdcen = half*(vdlft+vdrgt)
              vdsgn = sign(one, vdcen)
              vslop = xslope_type*min(abs(vdlft),abs(vdrgt))
              vslop = merge(vslop,zero,(vdlft*vdrgt)>zero)
              dbf(:,i,j,k,2,1) = vdsgn*min(vslop,abs(vdcen))
              ! slopes in second coordinate direction
              vdlft = (bf(:,i,j,k  ,2) - bf(:,i,j,k-1,2))
              vdrgt = (bf(:,i,j,k+1,2) - bf(:,i,j,k  ,2))
              vdcen = half*(vdlft+vdrgt)
              vdsgn = sign(one, vdcen)
              vslop = xslope_type*min(abs(vdlft),abs(vdrgt))
              vslop = merge(vslop,zero,(vdlft*vdrgt)>zero)
              dbf(:,i,j,k,2,2) = vdsgn*min(vslop,abs(vdcen))
           end do
        end do
     end do
  elseif (slope_type==6) then
     do k = klo, khi
        do j = jlo, jhi+1 ! WARNING HERE
           do i = ilo, ihi
              ! slopes in first coordinate direction
              vdlft = (bf(:,i  ,j,k,2) - bf(:,i-1,j,k,2))
              vdrgt = (bf(:,i+1,j,k,2) - bf(:,i  ,j,k,2))
              vslop = (vdlft**2+vdlft*vdrgt)/(vdlft**2+vdrgt**2+1e-30)
              vslop = vdlft*(1.-vslop)+vdrgt*vslop
              dbf(:,i,j,k,2,1) = merge(vslop,zero,(vdlft*vdrgt)>zero)
              ! slopes in second coordinate direction
              vdlft = (bf(:,i,j,k  ,2) - bf(:,i,j,k-1,2))
              vdrgt = (bf(:,i,j,k+1,2) - bf(:,i,j,k  ,2))
              vslop = (vdlft**2+vdlft*vdrgt)/(vdlft**2+vdrgt**2+1e-30)
              vslop = vdlft*(1.-vslop)+vdrgt*vslop
              dbf(:,i,j,k,2,2) = merge(vslop,zero,(vdlft*vdrgt)>zero)
           end do
        end do
     end do
  endif
  ! Bz along direction X and Y
  if(slope_type>=1.and.slope_type<=5)then  ! minmod or average==MC-limiter
     do k = klo, khi+1 ! WARNING HERE
        do j = jlo, jhi
           do i = ilo, ihi
              ! slopes in first coordinate direction
              vdlft = (bf(:,i  ,j,k,3) - bf(:,i-1,j,k,3))
              vdrgt = (bf(:,i+1,j,k,3) - bf(:,i  ,j,k,3))
              vdcen = half*(vdlft+vdrgt)
              vdsgn = sign(one, vdcen)
              vslop = xslope_type*min(abs(vdlft),abs(vdrgt))
              vslop = merge(vslop,zero,(vdlft*vdrgt)>zero)
              dbf(:,i,j,k,3,1) = vdsgn*min(vslop,abs(vdcen))
              ! slopes in second coordinate direction
              vdlft = (bf(:,i,j  ,k,3) - bf(:,i,j-1,k,3))
              vdrgt = (bf(:,i,j+1,k,3) - bf(:,i,j  ,k,3))
              vdcen = half*(vdlft+vdrgt)
              vdsgn = sign(one, vdcen)
              vslop = xslope_type*min(abs(vdlft),abs(vdrgt))
              vslop = merge(vslop,zero,(vdlft*vdrgt)>zero)
              dbf(:,i,j,k,3,2) = vdsgn*min(vslop,abs(vdcen))
           end do
        end do
     end do
  elseif (slope_type==6) then
     do k = klo, khi+1 ! WARNING HERE
        do j = jlo, jhi
           do i = ilo, ihi
              ! slopes in first coordinate direction
              vdlft = (bf(:,i  ,j,k,3) - bf(:,i-1,j,k,3))
              vdrgt = (bf(:,i+1,j,k,3) - bf(:,i  ,j,k,3))
              vslop = (vdlft**2+vdlft*vdrgt)/(vdlft**2+vdrgt**2+1e-30)
              vslop = vdlft*(1.-vslop)+vdrgt*vslop
              dbf(:,i,j,k,3,1) = merge(vslop,zero,(vdlft*vdrgt)>zero)
              ! slopes in second coordinate direction
              vdlft = (bf(:,i,j  ,k,3) - bf(:,i,j-1,k,3))
              vdrgt = (bf(:,i,j+1,k,3) - bf(:,i,j  ,k,3))
              vslop = (vdlft**2+vdlft*vdrgt)/(vdlft**2+vdrgt**2+1e-30)
              vslop = vdlft*(1.-vslop)+vdrgt*vslop
              dbf(:,i,j,k,3,2) = merge(vslop,zero,(vdlft*vdrgt)>zero)
           end do
        end do
     end do
  endif
  if(limit_slope>=1)then                      ! use additional limiter by Knut Waagan J. Comp Phys 228 (2009) 8609
      alphav = 0.8_dp / (1.0_dp + gamma)/limit_slope      ! to ensure positive density and pressure. Inspired by his implementation in Flash v2.5
     palphav = alphav * dx/dt                                                   ! Limit dvpar to 80/(1+gamma)% compared to dx/dt
     malphav = - palphav
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              vlimd = 0.40_dp*q(:,i,j,k,1)/limit_slope                          ! Limit drho to 40%
              vlimp = alphav*q(:,i,j,k,5)                                       ! Limit dpressure to 80/(1+gamma)%
              do idir=1,3
                 dq(:,i,j,k,1,     idir) = max( -vlimd,min(dq(:,i,j,k,1,     idir),vlimd))
                 dq(:,i,j,k,1+idir,idir) = max(malphav,min(dq(:,i,j,k,1+idir,idir),palphav)) 
                 dq(:,i,j,k,5,     idir) = max( -vlimp,min(dq(:,i,j,k,5,     idir),vlimp))
              end do
           end do
        end do
     end do
     ! Check that corner states are ok too :
      alphav = 0.9_dp / (1.0_dp + gamma)/limit_slope                            ! Limit dvpar to 90/(1+gamma)% compared to dx/dt
     palphav = alphav * dx/dt           
     malphav = - palphav
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              vlimd = 0.45_dp*q(:,i,j,k,1)/limit_slope                          ! Limit drho to 45%
              vlimp = alphav*q(:,i,j,k,5)                                       ! Limit dpressure to 90/(1+gamma)%
              do l=1,ngrid
                 if ( abs(dq(l,i,j,k,1,2))+abs(dq(l,i,j,k,1,3)) > vlimd(l) .or. &
                      abs(dq(l,i,j,k,1,3))+abs(dq(l,i,j,k,1,1)) > vlimd(l) .or. &
                      abs(dq(l,i,j,k,1,1))+abs(dq(l,i,j,k,1,2)) > vlimd(l) .or. &
                      abs(dq(l,i,j,k,5,2))+abs(dq(l,i,j,k,5,3)) > vlimp(l) .or. &
                      abs(dq(l,i,j,k,5,3))+abs(dq(l,i,j,k,5,1)) > vlimp(l) .or. &
                      abs(dq(l,i,j,k,5,1))+abs(dq(l,i,j,k,5,2)) > vlimp(l) .or. &
                      abs(dq(l,i,j,k,3,2))+abs(dq(l,i,j,k,4,3)) > palphav .or. &
                      abs(dq(l,i,j,k,4,3))+abs(dq(l,i,j,k,2,1)) > palphav .or. &
                      abs(dq(l,i,j,k,2,1))+abs(dq(l,i,j,k,3,2)) > palphav ) then
                    do idir=1,3
                       dq(l,i,j,k,1,     idir) = max( -half*vlimd(l),min(dq(l,i,j,k,1,     idir),half*vlimd(l)))
                       dq(l,i,j,k,1+idir,idir) = max(  half*malphav, min(dq(l,i,j,k,1+idir,idir),half*palphav )) 
                       dq(l,i,j,k,5,     idir) = max( -half*vlimp(l),min(dq(l,i,j,k,5,     idir),half*vlimp(l)))
                    end do
                 endif
              end do
           end do
        end do
     end do
  endif
#endif
  
!  call trace_mpi('uslope','exit',2)
end subroutine uslope
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine aslope(q,dq,ngrid,rslope_type)
  use amr_parameters
  use hydro_parameters, only : iu1,iu2,ju1,ju2,ku1,ku2
  implicit none

  integer::ngrid
  real(dp)::rslope_type
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::q 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::dq

  ! local variables
  integer::i,j,k,l
  integer::ilo,ihi,jlo,jhi,klo,khi
  real(dp)::dfx,dfy,dfz,fmin,fmax,qmin,qmax,dlim
  real(dp),dimension(-1:1,-1:1,-1:1)::f
  
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  do k = klo, khi
    do j = jlo, jhi
      do i = ilo, ihi
        do l = 1, ngrid

          ! Gradient components
          dfx  = 0.5_dp*(q(l,i+1,j,k)-q(l,i-1,j,k))
          dfy  = 0.5_dp*(q(l,i,j+1,k)-q(l,i,j-1,k))
          dfz  = 0.5_dp*(q(l,i,j,k+1)-q(l,i,j,k-1))

          ! Difference values based on the gradient
          f(-1, 0,-1) = -dfx    -dfz
          f( 0, 0,-1) =         -dfz
          f(+1, 0,-1) = +dfx    -dfz
          f(-1,-1,-1) = -dfx-dfy-dfz
          f( 0,-1,-1) =     -dfy-dfz
          f(+1,-1,-1) = +dfx-dfy-dfz
          f(-1,+1,-1) = -dfx+dfy-dfz
          f( 0,+1,-1) =     +dfy-dfz
          f(+1,+1,-1) = +dfx+dfy-dfz
          f(-1, 0, 0) = -dfx
          f( 0, 0, 0) =  0.0
          f(+1, 0, 0) = +dfx
          f(-1,-1, 0) = -dfx-dfy
          f( 0,-1, 0) =     -dfy
          f(+1,-1, 0) = +dfx-dfy
          f(-1,+1, 0) = -dfx+dfy
          f( 0,+1, 0) =     +dfy
          f(+1,+1, 0) = +dfx+dfy
          f(-1, 0,+1) = -dfx    +dfz
          f( 0, 0,+1) =         +dfz
          f(+1, 0,+1) = +dfx    +dfz
          f(-1,-1,+1) = -dfx-dfy+dfz
          f( 0,-1,+1) =     -dfy+dfz
          f(+1,-1,+1) = +dfx-dfy+dfz
          f(-1,+1,+1) = -dfx+dfy+dfz
          f( 0,+1,+1) =     +dfy+dfz
          f(+1,+1,+1) = +dfx+dfy+dfz

          ! Max positive and negative changes, separately
          fmax = abs(maxval(f))
          fmin = abs(minval(f))
          qmin = abs(minval(min(q(l,i-1:i+1,j-1:j+1,k-1:k+1)-q(l,i,j,k),0.0)))
          qmax = abs(maxval(max(q(l,i-1:i+1,j-1:j+1,k-1:k+1)-q(l,i,j,k),0.0)))

          ! If any of the function values exceed the predictions, reduce
          ! the gradient correspondingly.  This is quite conservative, since
          ! the actual change is only dt/ds*u*gradient, which is significantly
          ! smaller than the full change. Before possibly adopting this for
          ! use in the general case one should consider how to parametrize 
          ! being more or less aggressive.
          dlim = min(qmax/(fmax+tiny(1d0)),qmin/(fmin+tiny(1d0)),1.0)

          dq(l,i,j,k,1) = dlim*dfx
          dq(l,i,j,k,2) = dlim*dfy
          dq(l,i,j,k,3) = dlim*dfz
        end do
      end do
    end do
  end do
end subroutine aslope
