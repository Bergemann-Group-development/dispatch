!define DBG
!===============================================================================
!> Riemann solvers for HD
!===============================================================================
MODULE riemann_mod
  implicit none
  ! private
  ! type, public:: riemann_t
  !   real(8)          :: g1, g2
  !   character(len=4) :: method='hll'
  ! contains
  !   procedure        :: solve
  ! end type
  ! type(riemann_t), public:: riemann
CONTAINS

!===============================================================================
!> Interface to choose method
!===============================================================================
SUBROUTINE solve (self_g1, self_g2, dir, l, u, left, rght, flux)
  !class(riemann_t):: self
  real(8) :: self_g1, self_g2
  integer:: dir, l(3), u(3)
  real, dimension(:,:,:,:), contiguous, intent(in) :: left
  real, dimension(:,:,:,:),             intent(in) :: rght
  real, dimension(:,:,:,:), contiguous, intent(out):: flux
  !$omp declare target
  !.............................................................................
  ! select case (trim(self%method))
  ! case ('hll')
    call hll  (self_g1, self_g2, dir, l, u, left, rght, flux)
  ! case ('hllc')
  !   call hllc (self, dir, l, u, left, rght, flux)
  ! end select
END SUBROUTINE solve

!===============================================================================
!> Solve for the HLL flux give two HD state vectors (left and right of the
!> interface) stored in s
!===============================================================================
SUBROUTINE hll (self_g1, self_g2, dir, l, u, left, rght, flux)
  !class(riemann_t):: self
  real(8) :: self_g1, self_g2
  integer:: dir, l(3), u(3)
  real, dimension(:,:,:,:), contiguous, intent(in) :: left
  real, dimension(:,:,:,:),             intent(in) :: rght
  real, dimension(:,:,:,:), contiguous, intent(out):: flux
  !$omp declare target
  !.............................................................................
  real:: &
    Dl, Ul, Pl, El, Vl, Wl, Tl, cl, fDl, fTl, fUl, fVl, fWl, Sl, &
    Dr, Ur, Pr, Er, Vr, Wr, Tr, cr, fDr, fTr, fUr, fVr, fWr, Sr, &
    cmax
  integer:: i1, i2, i3, j1, j2, j3, v1, v2, v3
  integer, save:: itimer=0

  !-----------------------------------------------------------------------------
  !  Parallel and perpendicular velocity components
  !-----------------------------------------------------------------------------
  v1 = 2 + mod(dir-1,3)
  v2 = 2 + mod(dir  ,3)
  v3 = 2 + mod(dir+1,3)
  !-----------------------------------------------------------------------------
  ! For the left states, we need to refer to index-1 in the flux direction
  !-----------------------------------------------------------------------------
  if (dir==1) then
    j1=1; j2=0; j3=0
  else if (dir==2) then
    j1=0; j2=1; j3=0
  else
    j1=0; j2=0; j3=1
  end if
  !-----------------------------------------------------------------------------

  !$omp teams distribute parallel do schedule(static,1) shared(self_g1, self_g2, dir, l, u, left, rght, flux, v1, v2, v3, j1, j2, j3) private(Sr, Sl, fWr, fWl, fVr, fVl, fTr, fTl, fUr, fUl, fDr, fDl, Tr, Tl, cmax, cr, cl, Pr, Pl, Wr, Vr, Ur, Er, Dr, Wl, Dl, El, Ul, Vl) default(none) collapse(3)
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    Dl = exp(&
         left(i1-j1,i2-j2,i3-j3,1))
    El = exp(&
         left(i1-j1,i2-j2,i3-j3,5))
    Ul = left(i1-j1,i2-j2,i3-j3,v1)
    Vl = left(i1-j1,i2-j2,i3-j3,v2)
    Wl = left(i1-j1,i2-j2,i3-j3,v3)
    Dr = exp(&
         rght(i1,i2,i3,1))
    Er = exp(&
         rght(i1,i2,i3,5))
    Ur = rght(i1,i2,i3,v1)
    Vr = rght(i1,i2,i3,v2)
    Wr = rght(i1,i2,i3,v3)
    !---------------------------------------------------------------------------
    ! Pressure and sound speed for each side of interface (l==left, r==right)
    !---------------------------------------------------------------------------
    Pl = Dl*El*self_g1
    Pr = Dr*Er*self_g1
    cl = sqrt(self_g2*El)
    cr = sqrt(self_g2*Er)
    cmax = max(cl,cr)
    !---------------------------------------------------------------------------
    ! Total energy per unit volume
    !---------------------------------------------------------------------------
    Tl = Dl*(El + 0.5*(Ul**2 + Vl**2 + Wl**2))
    Tr = Dr*(Er + 0.5*(Ur**2 + Vr**2 + Wr**2))
    !---------------------------------------------------------------------------
    ! Hydro fluxes
    !---------------------------------------------------------------------------
    fDl = Dl*Ul
    fDr = Dr*Ur
    fUl = Dl*Ul*Ul + Pl
    fUr = Dr*Ur*Ur + Pr
    fTl = Ul*(Tl + Pl)
    fTr = Ur*(Tr + Pr)
    fVl = Dl*Ul*Vl
    fVr = Dr*Ur*Vr
    fWl = Dl*Ul*Wl
    fWr = Dr*Ur*Wr
    !---------------------------------------------------------------------------
    ! maximum wave speeds to the l and r (guaranteed to have correct sign)
    !---------------------------------------------------------------------------
    SL = min(min(Ul,Ur)-cmax,0.0)
    SR = max(max(Ul,Ur)+cmax,0.0)
    !---------------------------------------------------------------------------
    ! HLL flux based on wavespeeds.  The general form is 
    !    (SR*F_l - SL*F_r + SL*SR *(u_r - u_l))/(SR-SL)
    ! where u = (rho, rho*v, E_tot) are the conserved variables
    !---------------------------------------------------------------------------
    flux(i1,i2,i3, 1) = (SR * fDl - SL*fDr + SL*SR*(Dr    - Dl   )) / (SR - SL)
    flux(i1,i2,i3,v1) = (SR * fUl - SL*fUr + SL*SR*(Dr*Ur - Dl*Ul)) / (SR - SL)
    flux(i1,i2,i3,v2) = (SR * fVl - SL*fVr + SL*SR*(Dr*Vr - Dl*Vl)) / (SR - SL)
    flux(i1,i2,i3,v3) = (SR * fWl - SL*fWr + SL*SR*(Dr*Wr - Dl*Wl)) / (SR - SL)
    flux(i1,i2,i3, 5) = (SR * fTl - SL*fTr + SL*SR*(Tr    - Tl   )) / (SR - SL)
#ifdef DBG
    if (dir==1 .and. i2==l(2) .and. i3==l(3)) then
      print '(i3,3(2x,5f10.4))', &
        i1, Dl,ul,el,Pl,SL, Dr,ur,er,Pr,SR, flux(i1,i2,i3,1:5)
      end if
#endif
  end do
  end do
  end do
 
  !$omp end teams distribute parallel do
  
  
 END SUBROUTINE hll

! !===============================================================================
! !> Solve for the HLLC flux give two HD state vectors (left and right of the
! !> interface) stored in fgdnv
! !===============================================================================
! SUBROUTINE hllc (self, dir, l, u, qleft, qrght, fgdnv)
!   implicit none
!   class(riemann_t):: self
!   !-----------------------------------------------------------------------------
!   ! Parameters
!   !-----------------------------------------------------------------------------
!   integer, parameter:: sp=4, dp=8
!   real, parameter:: half=0.5, zero=0.0, one=1.0
!   !--- single precision
!   REAL(sp),dimension(:,:,:,:):: qleft,qrght,fgdnv
!   REAL(sp):: rl,pl,ul,ecinl,etotl,ptotl,el,vl,wl
!   REAL(sp):: rr,pr,ur,ecinr,etotr,ptotr,er,vr,wr
!   REAL(sp):: rstarl,etotstarl
!   REAL(sp):: rstarr,etotstarr
!   REAL(sp):: ro,uo,ptoto,etoto
!   REAL(sp):: ptotstar
!   !--- double precision
!   REAL(dp):: cfastl,SL,rcl
!   REAL(dp):: cfastr,SR,rcr
!   REAL(dp):: ustar
!   !---
!   INTEGER:: l(3),u(3),i2,i3,dir,n_cheap,n_expensive
!   INTEGER:: i,j1,j2,j3,v1,v2,v3
!   INTEGER, SAVE:: itimer=0
!   !-----------------------------------------------------------------------------
!   ! Choose variable indices; v1=parallel, v2-v3: perpendicular
!   !-----------------------------------------------------------------------------
!   v1 = 2 + mod(dir-1,3)
!   v2 = 2 + mod(dir  ,3)
!   v3 = 2 + mod(dir+1,3)
!   n_cheap = 0
!   n_expensive = 0
!   !-----------------------------------------------------------------------------
!   ! For the left states, we need to refer to index-1 in the flux direction
!   !-----------------------------------------------------------------------------
!   if (dir==1) then
!     j1=1; j2=0; j3=0
!   else if (dir==2) then
!     j1=0; j2=1; j3=0
!   else
!     j1=0; j2=0; j3=1
!   end if
!   !---------------------------------------------------------------------------
!   ! Loops over 2nd and 3rd index
!   !---------------------------------------------------------------------------
!   do i3=l(3),u(3)
!   do i2=l(2),u(2)
!     !---------------------------------------------------------------------------
!     ! Loop over 1st index
!     !---------------------------------------------------------------------------
!     do i=l(1),u(1)
!       !
!       ! Left variables
!       rl = exp(qleft(i-j1,i2-j2,i3-j3,1))
!       el = exp(qleft(i-j1,i2-j2,i3-j3,5))
!       ul=      qleft(i-j1,i2-j2,i3-j3,v1)
!       vl=      qleft(i-j1,i2-j2,i3-j3,v2)
!       wl=      qleft(i-j1,i2-j2,i3-j3,v3)
!       !
!       ! Right variables
!       rr = exp(qrght(i,i2,i3,1))
!       er = exp(qrght(i,i2,i3,5))
!       ur=      qrght(i,i2,i3,v1)
!       vr=      qrght(i,i2,i3,v2)
!       wr=      qrght(i,i2,i3,v3)
!       !
!       ! Find the largest eigenvalues in the normal direction to the interface
!       cfastl = sqrt(self%g2*el)
!       cfastr = sqrt(self%g2*er)
!       !
!       ! Compute HLL wave speed; SL positive when positive supersonic, SR negative when negative ss
!       SL = min(ul,ur) - max(cfastl,cfastr)
!       SR = max(ul,ur) + max(cfastl,cfastr)
!       if (SL > 0.0) then
!         !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'SL case',i
!         ro = rl
!         uo = ul
!         Ptoto = rl*el*self%g1
!         etoto = rl*el + half*rl*(ul**2+vl**2+wl**2)
!         n_cheap = n_cheap+1
!       else if (SR < 0.0) then
!         !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'SR case',i
!         ro = rr
!         uo = ur
!         Ptoto = rr*er*self%g1
!         etoto = rr*er + half*rr*(ur**2+vr**2+wr**2)
!         n_cheap = n_cheap+1
!       else
!         n_expensive = n_expensive+1
!         !
!         ! Compute acoustic star state and lagrangian sound speed
!         rcl = rl*(ul-SL)
!         rcr = rr*(SR-ur)
!         Ptotl = rl*el*self%g1
!         Ptotr = rr*er*self%g1
!         ustar =    (rcr*ur    + rcl*ul    +   (Ptotl-Ptotr))/(rcr+rcl)
!         Ptotstar = (rcr*Ptotl + rcl*Ptotr + rcl*rcr*(ul-ur))/(rcr+rcl)
!         !
!         ! Sample the solution at x/t=0
!         if (ustar > 0.0)then
!           !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'rstarl case',i
!           ! Left star region variables
!           ro = rl*(SL-ul)/(SL-ustar)  ! = rstarl
!           uo = ustar
!           Ptoto = Ptotstar
!           etotl = rl*el + half*rl*(ul**2+vl**2+wl**2)
!           etoto = ((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar)/(SL-ustar) ! = etotstarl
!         else
!           ! Right star region variables
!           !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'rstarr case',i
!           ro = rr*(SR-ur)/(SR-ustar)  ! = rstarr
!           uo = ustar
!           Ptoto = Ptotstar
!           etotr = rr*er + half*rr*(ur**2+vr**2+wr**2)
!           etoto = ((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar) ! = etotstarr
!         end if
!       end if
!       !
!       ! Compute the Godunov flux
!       fgdnv(i,i2,i3, 1) = ro*uo
!       fgdnv(i,i2,i3, 5) = (etoto+Ptoto)*uo
!       fgdnv(i,i2,i3,v1) = ro*uo*uo+Ptoto
!       fgdnv(i,i2,i3,v2) = fgdnv(i,i2,i3,1)*merge(vl,vr,uo>0.0)
!       fgdnv(i,i2,i3,v3) = fgdnv(i,i2,i3,1)*merge(wl,wr,uo>0.0)
! #ifdef DBG
!       if (dir==1 .and. i2==l(2) .and. i3==l(3)) then
!         Pl = rl*el*self%g1
!         Pr = rr*er*self%g1
!         print '(i3,3(2x,5f10.4))', &
!           i, rl,ul,el,Pl,SL, rr,ur,er,Pr,SR, fgdnv(i,i2,i3,1:5)
!       end if
! #endif
!     end do
!   end do
!   end do
!   !-----------------------------------------------------------------------------
!   ! Update solver statistics
!   !-----------------------------------------------------------------------------
! END SUBROUTINE hllc

END MODULE riemann_mod
