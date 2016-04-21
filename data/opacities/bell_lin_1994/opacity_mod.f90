!> =============================================================================
!> Opacity module for the piece-wise analytic opacities of Bell & Lin
!> (1994, ApJ, 427, 987).
!> =============================================================================
MODULE opacity_mod
  USE trace_mod
  implicit none
  private
    
  interface mean_opacity
    module procedure :: BellLinOpacity
  end interface
  public mean_opacity

  real(8), parameter :: onethird = 1.0/3.0, twothird = 2.0/3.0, beta = 2.0
  ! data for the indivdual pieces of the opacity
  real(8), dimension(8) :: ki, aexp, bexp
  data  ki  &
        / 2.0d-4, 2.0d16, 0.1, 2.0d81, 1.0d-8, 1.0d-36, 1.5d20, 0.348 /
  data  aexp  &
        / 0.0, 0.0, 0.0, 1.0, twothird, onethird, 1.0, 0.0 /
  data  bexp  &
        / 2.0, -7.0, 0.5, -24.0, 3.0, 10.0, -2.5, 0.0 /

CONTAINS

!> =============================================================================
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
!> =============================================================================
FUNCTION BellLinOpacity (rho, temp) RESULT (kappa)
  real(4), intent(in) :: rho, temp
  real(4) :: kappa
  integer :: i
  real(8) :: k123, k(8), kappad
  !.............................................................................
  call trace_begin('opacity_mod::BellLinOpacity')

  ! each piece of the opacity function is, in general, a power-law in density
  ! and temperature. It just might be that, for a particular piece, one of the
  ! exponents is zero.
  do i=1,8
    k(i) = ki(i) * rho**aexp(i) * temp**bexp(i)
  end do

  ! piece them all together and smooth out the transition points.
  k123 = ( ( k(1)**(-2) + k(2)**(-2) )**(-2) &
         + ( k(3) / ( 1.0 + 1.0d22 * temp**(-10) ) )**4 )**0.25
  kappad = ( ( ( k123**(-beta) + k(4)**(-beta) )**(-1) &
              + k(5)**beta + k(6)**beta )**(-1) + k(7)**(-beta) )**(-1.0/beta)
  kappad = ( kappad**beta &
           + ( k(8) / ( 1.0 + 1.0d41 * temp**(-10) ) )**beta )**(1.0/beta)
  kappa = real(kappad, kind=4)

  call trace_end
END FUNCTION BellLinOpacity

END MODULE opacity_mod
