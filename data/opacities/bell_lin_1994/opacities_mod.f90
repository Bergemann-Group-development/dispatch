!===============================================================================
!> Opacity module for the piece-wise analytic opacities of Bell & Lin
!> (1994, ApJ, 427, 987).
!===============================================================================
MODULE opacities_mod
  USE trace_mod
  USE scaling_mod
  implicit none

  type, public ::  opacity_t
    real(8) :: onethird = 1.0/3.0
    real(8) :: twothird = 2.0/3.0
    real(8) :: beta = 2.0
    real(8) :: sc_temp, sc_rho, sc_kap
    real(8), dimension(8) :: ki, aexp, bexp
  contains
    procedure :: init   => init_bell_lin
    procedure :: solve  => BellLinOpacity
  end type
CONTAINS

SUBROUTINE init_bell_lin(self)
  class(opacity_t) :: self
  ! --------------------------------------------------------------------
  call trace_begin('opacities_mod::init_bell_lin',3)
  self%ki   = (/2.0d-4, 2.0d16, 0.1_8, 2.0d81, 1.0d-8, 1.0d-36, 1.5d20, 0.348_8 /)
  self%aexp = (/ 0d0, 0d0, 0d0, 1.0_8, self%twothird, self%onethird, 1.0_8, 0.0_8 /)
  self%bexp = (/ 2.0_8, -7.0_8, 0.5_8, -24.0_8, 3.0_8, 10.0_8, -2.5_8, 0.0_8 /)
  self%sc_rho  = 1_8/scaling%d
  self%sc_temp = 1_8/scaling%temp
  self%sc_kap  = scaling%m/scaling%l**2
  call trace_end
END SUBROUTINE init_bell_lin

!===============================================================================
!> Calculate a temperature and density dependent opacity following the piece-wise
!> analytic functions of Bell & Lin (1994, ApJ, 427, 987).
!> Near the location of transitions between pieces, the opacity is smoothed with
!> what is roughly a generalised power mean.
!> The inputs are assumed to have units of g cm**-3 and K, respectively.
!>
!> NOTES:
!> - The opacity is generally a much stronger function of temperature than density.
!> - Thanks to Kees Dullemond for explaining the smoothing procedure (which wasn't
!>   document anywhere).
!> - The smoothing procedure is specific to these opacities.  If a different
!>   opacity formula is used, the smoothing will need to be re-derived.  See JPR's
!>   notes titled "Lin & Papaloizou analytic opacity function" from April 2013.
!===============================================================================
SUBROUTINE BellLinOpacity (self,rho0, temp0, nwb, kappa)
  class(opacity_t) :: self
  real(4), dimension(:,:,:), pointer :: rho0, temp0
  real(4), dimension(:,:,:,:), pointer :: kappa
  real(8), dimension(size(rho0,1),size(rho0,2),size(rho0,3)) :: k123, kappad, rho, temp
  real(8), dimension(size(rho0,1),size(rho0,2),size(rho0,3), 8) :: k
  integer :: i, nwb
  !real (8) :: sc_kap , sc_rho, wc_t
  !.............................................................................
  call trace_begin('opacity_mod::BellLinOpacity',3)

  ! each piece of the opacity function is, in general, a power-law in density
  ! and temperature. It just might be that, for a particular piece, one of the
  ! exponents is zero.
  associate (beta => self%beta)
  
  rho = rho0 * self%sc_rho
  temp = temp0 * self%sc_temp
  do i=1,8
    k(:,:,:,i) = self%ki(i) * rho**self%aexp(i) * temp**self%bexp(i)
  end do
  
  ! piece them all together and smooth out the transition points.
  k123 = ( ( k(:,:,:,1)**(-2) + k(:,:,:,2)**(-2) )**(-2) &
         + ( k(:,:,:,3) / ( 1.0 + 1.0d22 * temp**(-10) ) )**4 )**0.25
  kappad = ( ( ( k123**(-beta) +     k(:,:,:,4)**(-beta) )**(-1) &
              + k(:,:,:,5)**  beta + k(:,:,:,6)**beta )   **(-1) &
              + k(:,:,:,7)**(-beta) )**(-1.0/beta)
  kappad = ( kappad**beta &
           + ( k(:,:,:,8) / ( 1.0 + 1.0d41 * temp**(-10) ) )**beta )**(1.0/beta)
  kappad = kappad * self%sc_kap
  kappa(:,:,:,1) = real(kappad, kind=4)
  
  if (nwb>1) then
    do i= 2, nwb
      kappa(:,:,:,i) = kappa(:,:,:,1)
    end do
  end if
  end associate
  call trace_end
END SUBROUTINE BellLinOpacity

END MODULE opacities_mod
