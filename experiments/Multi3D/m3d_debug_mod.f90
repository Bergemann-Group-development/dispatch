module m3d_debug_mod
  implicit none

  real(4), public, dimension(6) :: memcheck
  integer, public :: pcheck(3)=[1,1,1], fcheck=1, linecheck=3, lvlcheck=1
  logical, public :: pdebug=.false.
  !$omp threadprivate(memcheck, pdebug)
  
contains
  
end module m3d_debug_mod