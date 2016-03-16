! $Id$
!> =============================================================================
MODULE omp_lock_mod
  USE omp_mod
  USE omp_lib
  implicit none
  private
  type, public:: lock_t
    integer:: id=0
    integer(kind=omp_lock_kind):: lock
    character(len=4):: kind
  contains
    procedure:: init
    procedure:: set
    procedure:: test
    procedure:: unset
    procedure:: destroy
  end type
  integer:: id=0
  logical:: debug=.false.
CONTAINS

!> =============================================================================
SUBROUTINE init (self, kind)
  class(lock_t):: self
  character(len=4):: kind
  !-----------------------------------------------------------------------------
  if (self%id > 0) return
  call omp_init_lock (self%lock)
  !$omp atomic
  id = id+1
  self%id = id
  self%kind = kind
  if (debug) &
    print*, omp_mythread(),' initialized lock', self%id, ' '//kind
END SUBROUTINE init

!> =============================================================================
SUBROUTINE set (self, label)
  class(lock_t):: self
  character(len=*):: label
  !-----------------------------------------------------------------------------
  if (omp_test_lock (self%lock)) then
    if (debug) &
      print '(i12,a,i14,1x,a)', omp_mythread(),'  set lock', self%id, self%kind//': '//trim(label)
  else
    if (debug) &
      print '(i12,a,i14,1x,a)', omp_mythread(),' wait lock', self%id, self%kind//': '//trim(label)
    call omp_set_lock (self%lock)
    if (debug) &
      print '(i12,a,i14,1x,a)', omp_mythread(),'  set lock', self%id, self%kind//': '//trim(label)
  end if
END SUBROUTINE set

!> =============================================================================
FUNCTION test (self) RESULT (out)
  class(lock_t):: self
  logical:: out
  !-----------------------------------------------------------------------------
  if (debug) &
    print*, omp_mythread(),'  test lock', self%id, self%kind
  out =  omp_test_lock (self%lock)
END FUNCTION test

!> =============================================================================
SUBROUTINE unset (self, label)
  class(lock_t):: self
  character(len=*):: label
  !-----------------------------------------------------------------------------
  if (debug) &
    print*, omp_mythread(),' unset lock', self%id, self%kind//': '//trim(label)
  call omp_unset_lock (self%lock)
END SUBROUTINE unset

!> =============================================================================
SUBROUTINE destroy (self)
  class(lock_t):: self
  !-----------------------------------------------------------------------------
  if (debug) &
    print*, omp_mythread(),' destr lock', self%id, self%kind
  call omp_destroy_lock (self%lock)
END SUBROUTINE destroy

!> =============================================================================
END MODULE omp_lock_mod
