PROGRAM test
  USE opacity_mod
  USE scaling_mod
  implicit none
  integer, parameter:: n=10
  real:: rho(n), T(n), kap(n)
  !-----------------------------------------------------------------------------
  rho = 1.
  T = 100.
  kap = opacity%mean(rho,T)
  print *,kap
  T = 200.
  kap = opacity%mean(rho,T)
  print *,kap
END PROGRAM test
