MODULE boundary_wavekill_mod
  USE kinds_mod
  implicit none
  private
  
PUBLIC wave_killing
CONTAINS
!===============================================================================
!> Following Section 3.2 of de Val-Borro et al. (2006, MNRAS, 370, 529).
!> These "wave killing" conditions are designed to reduce wave reflection at
!> solid boundaries.
!>
!> In this approach, the following ODE is applied at each application of the
!> boundary conditions:
!>    dq = - q - q0
!>    --     ------ * R(x)
!>    dt      tau
!> where `q` is the quantity that you wish to damp towards value `q0`, `tau` is
!> the timescale over which `q` is damped, `xbndy` is the location of the actual
!> boundary, `xkill` is the location at which to begin the damping procedure and
!> `x` is the location of the current grid point.
!>
!> `R(x)` is a simple quadratic weighting function that is 1 at the actual boundary
!> (`xbndy`) and 0 at `xkill`.
!>
!> The function returns the *rate* at which quantity `q` is damped; the quantity
!> should then be updated along the lines of:
!>    qold = qnew + dt * wave_killing(qold, q0, tau, xkill, xbndy, x)
!>
!===============================================================================
FUNCTION wave_killing (q, q0, tau, xkill, xbndy, x) RESULT (rate)
  real(kind=KindScalarVar), intent(in) :: q, q0
  real(kind=KindScalarVar), intent(in) :: tau, x, xkill, xbndy
  real(kind=KindScalarVar) :: rate, quadfunc
  !-----------------------------------------------------------------------------

  quadfunc = ((x - xkill) / (xbndy - xkill))**2
  rate = -(q - q0) * quadfunc / tau

END FUNCTION wave_killing

END MODULE boundary_wavekill_mod
