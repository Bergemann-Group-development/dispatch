!> =============================================================================
!> Program to test functionality and speed of particles, when using an array
!> of structures
!> =============================================================================

!> =============================================================================
PROGRAM test
  USE p_arrays_mod
  implicit none
  integer:: n, np, niter, i1, i2, i3
  type(p_arrays_t):: pa
  !-----------------------------------------------------------------------------
  n = 16
  np = 10
  niter = 4
  call pa%test (np, n, niter)
END PROGRAM
