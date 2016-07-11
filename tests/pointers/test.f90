!> =============================================================================
!> Test showing that Fortran pointers are just pointers, not pointers to pointers
!> =============================================================================
PROGRAM test_pointers
  real, pointer:: a,b,c,d
  allocate (a,b)
  a=1.
  b=2.
  c=>a
  d=>c
  print*,d
  c=>b
  print*,d
END PROGRAM
