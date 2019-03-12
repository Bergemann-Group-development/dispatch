!===============================================================================
!> Test showing that Fortran pointers are just pointers, not pointers to pointers,
!> demonstrating that one can allocate on dimensioned class pointers, and showing
!> that "target, allocatable" is essentially equivalent to "pointer".
!===============================================================================
PROGRAM test_pointers
  type:: a_t
    real:: a
  end type
  type, extends(a_t):: b_t
    real:: b
  end type
  class(a_t), dimension(:), pointer:: aa, bb
  class(a_t), dimension(:), target, allocatable:: cc
  type(b_t), dimension(:), pointer:: bp, cp
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
  allocate (b_t::bb(5),cc(4))
  select type (bb)
  type is (b_t)
    bb(:)%a = 2.
    bb(:)%b = 3.
    bp => bb
  end select
  !-----------------------------------------------------------------------------
  ! A class array of allocatable with target attribute works the same as an array 
  ! pointers to class.
  !-----------------------------------------------------------------------------
  select type (cc)
  type is (b_t)
    cc(:)%a = 4.
    cc(:)%b = 5.
    cp => cc
  end select
  print *, cc%a
  print *, cp%b
END PROGRAM
