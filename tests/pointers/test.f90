!> =============================================================================
!> Test showing that Fortran pointers are just pointers, not pointers to pointers,
!> and demonstrating that one can allocate on dimensioned class pointers.
!> =============================================================================
PROGRAM test_pointers
  type:: a_t
    real:: a
  end type
  type, extends(a_t):: b_t
    real:: b
  end type
  class(a_t), dimension(:), pointer:: aa, bb
  type(b_t), dimension(:), pointer:: bp
  real, pointer:: a,b,c,d
  !.............................................................................
  allocate (a,b)
  a=1.
  b=2.
  c=>a
  d=>c
  print*,d
  c=>b
  print*,d
  !-----------------------------------------------------------------------------
  ! Test allocating on arrays of class pointers
  !-----------------------------------------------------------------------------
  allocate (a_t::aa(3))
  aa(:)%a = 1.
  print *, aa%a
  !-----------------------------------------------------------------------------
  ! Need to select type to access the %b part of pointers to class(a_t)
  !-----------------------------------------------------------------------------
  allocate (b_t::bb(5))
  select type (bb)
    type is (b_t)
    bb(:)%a = 2.
    bb(:)%b = 3.
    bp => bb
  end select
  print *, bb%a
  print *, bp%b
END PROGRAM
