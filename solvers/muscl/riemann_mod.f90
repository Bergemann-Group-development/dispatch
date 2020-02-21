!===============================================================================
!> Riemann solvers for HD
!===============================================================================
MODULE riemann_mod
  implicit none
  private
  type, public:: riemann_t
    real(8)     :: gamma
  contains
    procedure   :: hll
  end type
CONTAINS

!===============================================================================
!> Solve for the HLL flux give two HD state vectors (left and right of the
!> interface) stored in s
!===============================================================================
SUBROUTINE hll (self, dir, l, u, left, rght, flux)
  class(riemann_t):: self
  integer:: dir, l(3), u(3)
  !real, dimension(:,:,:,:), allocatable, target:: left, rght
  !real, dimension(:,:,:,:), contiguous         :: flux
  real, dimension(:,:,:,:), contiguous, intent(in) :: left
  real, dimension(:,:,:,:),             intent(in) :: rght
  real, dimension(:,:,:,:), contiguous, intent(out):: flux
  !.............................................................................
  real:: &
    Dl, Ul, Pl, El, Vl, Wl, Tl, cl, fDl, fTl, fUl, fVl, fWl, Sl, &
    Dr, Ur, Pr, Er, Vr, Wr, Tr, cr, fDr, fTr, fUr, fVr, fWr, Sr, &
    cmax
  integer:: i1, i2, i3, j1, j2, j3, v1, v2, v3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  v1 = 2 + mod(dir-1,3)
  v2 = 2 + mod(dir  ,3)
  v3 = 2 + mod(dir+1,3)
  !-----------------------------------------------------------------------------
  ! For the right states, we need to refer to ix-1 in the flux direction
  !-----------------------------------------------------------------------------
  if (dir==1) then
    j1=1; j2=0; j3=0
  else if (dir==2) then
    j1=0; j2=1; j3=0
  else
    j1=0; j2=0; j3=1
  end if
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    Dl = exp(&
         left(i1-j1,i2-j2,i3-j3,1))
    El = left(i1-j1,i2-j2,i3-j3,5)
    Ul = left(i1-j1,i2-j2,i3-j3,v1)
    Vl = left(i1-j1,i2-j2,i3-j3,v2)
    Wl = left(i1-j1,i2-j2,i3-j3,v3)
    Dr = exp(&
         rght(i1,i2,i3,1))
    Er = rght(i1,i2,i3,5)
    Ur = rght(i1,i2,i3,v1)
    Vr = rght(i1,i2,i3,v2)
    Wr = rght(i1,i2,i3,v3)
    !---------------------------------------------------------------------------
    ! Pressure and sound speed for each side of interface (l==left, r==right)
    !---------------------------------------------------------------------------
    Pl = Dl*El*(self%gamma-1d0)
    Pr = Dr*Er*(self%gamma-1d0)
    cl = sqrt(self%gamma*Pl/Dl)
    cr = sqrt(self%gamma*Pr/Dr)
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
  end do
  end do
  end do
END SUBROUTINE hll

END MODULE riemann_mod
