!===============================================================================
!> Module to handle running means
!===============================================================================
MODULE running_mean_mod
  USE trace_mod
  USE io_mod
  implicit none
  private
  type, public:: running_mean_t
    integer:: verbose=0                 ! verbosity
    integer:: current                   ! current index
    integer:: length = 10               ! running length
    real, pointer:: ring(:) => null()   ! values
    real:: average                      ! value value
  contains
    procedure:: init
    procedure:: value
    procedure:: update
    procedure:: test
    procedure:: dealloc
  end type
CONTAINS

!===============================================================================
!> Initialize a running mean
!===============================================================================
SUBROUTINE init (self, length, initial_value, verbose)
  class(running_mean_t):: self
  integer:: length
  real:: initial_value
  integer, optional:: verbose
  !.............................................................................
  call trace%begin ('running_mean_t%init')
  if (present(verbose)) &
    self%verbose = verbose
  self%length = length
  allocate (self%ring(length))
  if (self%verbose > 0) &
    write (io_unit%log,'(a,i8,1p,g14.5,l3)') &
     'running_mean_t%init: length, initial_value, associated', &
      length, initial_value, associated(self%ring)
  self%ring = initial_value
  self%average = initial_value
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Return the running mean
!===============================================================================
FUNCTION value (self, increment)
  class(running_mean_t):: self
  real:: increment, value
  !.............................................................................
  call trace%begin ('running_mean_t%value')
  call self%update (increment)
  value = self%average
  call trace%end ()
END FUNCTION value

!===============================================================================
!> Update the running mean
!===============================================================================
SUBROUTINE update (self, increment, all)
  class(running_mean_t):: self
  real:: increment
  logical, optional:: all
  !.............................................................................
  call trace%begin ('running_mean_t%update')
  if (.not. associated(self%ring)) then
     call self%init (self%length, increment)
  else if (present(all)) then
    self%ring = increment
    if (self%verbose > 1) &
      write (io_unit%log, *) 'running_mean: update(all) value =', increment
  else
    self%current = 1 + mod(self%current,self%length)
    self%ring (self%current) = increment
  end if
  self%average = sum(self%ring)/self%length
  if (self%verbose > 2) &
   write (io_unit%log, *) 'running_mean: increment, value =', &
     increment, self%average
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Create a running mean, starting at zero, and add the value 20.0, until the
!> average is 15.0.  Then add the value 0.0 until the value is 5.0, and return
!===============================================================================
SUBROUTINE test (self)
  class(running_mean_t):: self
  !-----------------------------------------------------------------------------
  print *, 'init'
  self%verbose = 1
  call self%init (10,1.0)
  print *, self%value(0.0)
  do while (self%value (20.0) < 15.0)
  end do
  do while (self%value (0.0) > 5.0)
  end do
  call self%dealloc
  self%verbose = 0
END SUBROUTINE test

!===============================================================================
!> Deallocate the ring
!===============================================================================
SUBROUTINE dealloc (self)
  class(running_mean_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('running_mean_t%update')
  print *, associated(self%ring)
  deallocate (self%ring)
  call trace%end ()
END SUBROUTINE dealloc
END MODULE running_mean_mod
