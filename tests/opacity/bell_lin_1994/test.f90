!===============================================================================
!> A program to test the functionality of the Bell & Lin opacity function.
!===============================================================================
PROGRAM test
  USE opacity_mod
  implicit none

  real(4) :: rho, T, kappa
  !.............................................................................
  rho = 1.0e-5  ! g/cm**3
  T = 1.0e2  ! K
  kappa = mean_opacity(rho, T) ! cm**2 /g
  print*,kappa

  rho = 1.0e-9 ! g/cm**3
  T = 1.0e4    ! K
  kappa = mean_opacity(rho, T) ! cm**2 /g
  print*,kappa

  rho = 1.0e-7 ! g/cm**3
  T = 4.0e3    ! K
  kappa = mean_opacity(rho, T) ! cm**2 /g
  print*,kappa

END PROGRAM
