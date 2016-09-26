!> =============================================================================
!> Test showing that Fortran pointers are just pointers, not pointers to pointers
!> =============================================================================
PROGRAM test_pointers
  type:: a_t
    real:: a
  end type
  type, extends(a_t):: b_t
    real:: b
  end type
  class(a_t), dimension(:), pointer:: aa
  real, pointer:: a,b,c,d
  allocate (a,b)
  allocate (a_t::aa(3))
  a=1.
  b=2.
  c=>a
  d=>c
  print*,d
  c=>b
  print*,d
END PROGRAM
