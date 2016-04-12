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
  integer,save:: id=0
  logical,save:: debug=.false.
CONTAINS

!> =============================================================================
!> Initialize a lock, taking care not to do it more than once
!> =============================================================================
SUBROUTINE init (self, kind)
  class(lock_t):: self
  character(len=4):: kind
  logical omp_in_parallel
  !-----------------------------------------------------------------------------
  !$omp critical
  if (self%id == 0) then
    call omp_init_lock (self%lock)
    id = id+1
    self%id = id
    self%kind = kind
    if (debug) &
      print*, omp_mythread(),' initialized lock:', self%id, ' '//kind
  end if
  !$omp end critical
END SUBROUTINE init

!> =============================================================================
SUBROUTINE set (self, label)
  class(lock_t):: self
  character(len=*):: label
  !-----------------------------------------------------------------------------
  if (omp_test_lock (self%lock)) then
    if (debug) then
      print '(i12,a,i14,1x,a)', omp_mythread(),'  set lock', self%id, self%kind//': '//trim(label)
      flush(6)
    end if
  else
    if (debug) then
      print '(i12,a,i14,1x,a)', omp_mythread(),' wait lock', self%id, self%kind//': '//trim(label)
      flush(6)
    end if
    call omp_set_lock (self%lock)
    if (debug) then
      print '(i12,a,i14,1x,a)', omp_mythread(),'  set lock', self%id, self%kind//': '//trim(label)
      flush(6)
    end if
  end if
END SUBROUTINE set

!> =============================================================================
FUNCTION test (self) RESULT (out)
  class(lock_t):: self
  logical:: out
  !-----------------------------------------------------------------------------
  if (debug) then
    print*, omp_mythread(),'  test lock', self%id, self%kind
    flush(6)
  end if
  out =  omp_test_lock (self%lock)
END FUNCTION test

!> =============================================================================
SUBROUTINE unset (self, label)
  class(lock_t):: self
  character(len=*):: label
  !-----------------------------------------------------------------------------
  if (debug) then
    print*, omp_mythread(),' unset lock', self%id, self%kind//': '//trim(label)
    flush(6)
  end if
  call omp_unset_lock (self%lock)
END SUBROUTINE unset

!> =============================================================================
SUBROUTINE destroy (self)
  class(lock_t):: self
  !-----------------------------------------------------------------------------
  if (debug) then
    print*, omp_mythread(),' destr lock', self%id, self%kind
    flush(6)
  end if
  call omp_destroy_lock (self%lock)
END SUBROUTINE destroy

!> =============================================================================
END MODULE omp_lock_mod
