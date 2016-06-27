! $Id$
!> =============================================================================
MODULE omp_lock_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_lib
  USE omp_timer_mod
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
  logical:: debug=.true.
  real(8):: wc0=0d0
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
  end if
  !$omp end critical
END SUBROUTINE init

!> =============================================================================
SUBROUTINE set (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  character(len=32):: filename
  !-----------------------------------------------------------------------------
  if (omp_test_lock (self%lock)) then
    if (debug) then
      if (present(label)) then
        write (log_unit,'(f15.7,2x,a,i14,1x,a)') wallclock()-wc0, '  set lock', self%id, self%kind//': '//trim(label)
        flush (log_unit)
      else
        write (log_unit,'(f15.7,2x,a,i14,1x,a)') wallclock()-wc0, '  set lock', self%id, self%kind
        flush (log_unit)
      end if
      flush(6)
    end if
  else
    if (debug) then
      if (present(label)) then
        write (log_unit,'(f15.7,2x,a,i14,1x,a)') wallclock()-wc0, ' wait lock', self%id, self%kind//': '//trim(label)
        flush (log_unit)
      else
        write (log_unit,'(f15.7,2x,a,i14,1x,a)') wallclock()-wc0, ' wait lock', self%id, self%kind
        flush (log_unit)
      end if
      flush(6)
    end if
    call omp_set_lock (self%lock)
    if (debug) then
      if (present(label)) then
        write (log_unit,'(f15.7,2x,a,i14,1x,a)') wallclock()-wc0, '  set lock', self%id, self%kind//': '//trim(label)
        flush (log_unit)
      else
        write (log_unit,'(f15.7,2x,a,i14,1x,a)') wallclock()-wc0, '  set lock', self%id, self%kind
        flush (log_unit)
      end if
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
    write (log_unit,'(f15.7,2x,a,i14,1x,a)') wallclock()-wc0, ' test lock', self%id, self%kind
    flush (log_unit)
  end if
  out =  omp_test_lock (self%lock)
END FUNCTION test

!> =============================================================================
SUBROUTINE unset (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (debug) then
    if (present(label)) then
      write (log_unit,'(f15.7,2x,a,i14,1x,a)') wallclock()-wc0, 'unset lock', self%id, self%kind//': '//trim(label)
      flush (log_unit)
    else
      write (log_unit,'(f15.7,2x,a,i14,1x,a)') wallclock()-wc0, 'unset lock', self%id, self%kind
      flush (log_unit)
    end if
  end if
  call omp_unset_lock (self%lock)
END SUBROUTINE unset

!> =============================================================================
SUBROUTINE destroy (self)
  class(lock_t):: self
  !-----------------------------------------------------------------------------
  if (debug) then
    print*, omp_mythread,' destr lock', self%id, self%kind
    flush(6)
  end if
  call omp_destroy_lock (self%lock)
END SUBROUTINE destroy

!> =============================================================================
END MODULE omp_lock_mod
