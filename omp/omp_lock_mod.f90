! $Id$
!> =============================================================================
MODULE omp_lock_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lib
  !USE trace_mod
  implicit none
  private
  type, public:: lock_t
    integer:: id=-1
    integer(kind=omp_lock_kind):: lock
    integer:: thread=-1
    character(len=4):: kind
    integer:: level=0
  contains
    procedure:: init
    procedure:: set
    procedure:: test
    procedure:: check
    procedure:: unset
    procedure:: destroy
  end type
  type(lock_t), public:: lock
  real(8):: lock_verbose = -1d0
  integer, save:: id_save=-1
CONTAINS

!> =============================================================================
!> Initialize a lock, taking care not to do it more than once
!> =============================================================================
SUBROUTINE init (self, kind, id)
  class(lock_t):: self
  integer, optional:: id
  character(len=4), optional:: kind
  logical omp_in_parallel
  !-----------------------------------------------------------------------------
  if (self%id == -1) then
    !$omp critical (lock_cr)
    if (self%id == -1) then
      call omp_init_lock (self%lock)
      if (present(id)) then
        self%id = id
        id_save = max(id_save,id)
      else
        id_save = id_save+1
        self%id = id_save
      end if
      if (present(kind)) self%kind = kind
      !print *, omp_get_thread_num(), 'initialized lock id', self%id
    end if
    !$omp end critical (lock_cr)
  end if
END SUBROUTINE init

!> =============================================================================
!> Return with the lock, waiting if necessary, and creating the lock if needed
!> =============================================================================
SUBROUTINE set (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  character(len=32):: filename
  integer:: thread
  real(8):: wc
  !-----------------------------------------------------------------------------
  call self%init
  !call trace%begin('lock_t%set')
  thread = omp_get_thread_num()
  !print *, thread, ' set begin id, thread', self%id, self%thread
  if (self%thread == thread .and. self%level > 0) then
    self%level = self%level+1
  else
    if (.not. self%test()) then
      call omp_set_lock (self%lock)
    end if
    !$omp flush
    if (self%level > 0 .or. self%thread > -1) then
      !write (io_unit%output,*) 'ERROR: lock_t%set', self%level, self%thread
      !self%level = 0
      !call trace%back
    end if
    self%thread = thread
    self%level = self%level+1
    wc = wallclock()
    if (wc < lock_verbose .and. present(label)) &
      write (io_unit%log,'(f10.6,i4,a,i6,i3,a)') wc, thread, ' lock id', self%id, self%level, '   set at '//label
    !print *, thread, ' set end id, thread', self%id, self%thread
  end if
  !call trace%end()
END SUBROUTINE set

!> =============================================================================
!> Return true if we already have the lock, or if we get it now
!> =============================================================================
FUNCTION test (self) RESULT (out)
  class(lock_t):: self
  logical:: out
  integer:: thread
  real(8):: wc
  !-----------------------------------------------------------------------------
  !call trace%begin('lock_t%test')
  thread = omp_get_thread_num()
  !print *, thread, ' set lock id, thread', self%id, self%thread
  !$omp flush
  !print *, thread, ' testing lock id', self%id, self%thread
  if (thread == self%thread) then
    out = .true.
    !print *, thread, ' already has lock id', self%id
  else
    !print *, thread, ' omp_test_lock for id', self%id
    out =  omp_test_lock (self%lock)
    if (out) then
      self%thread = thread
      self%level = self%level+1
      wc = wallclock()
      if (wc < lock_verbose) &
        write (io_unit%log,'(f10.6,i4,a,i6,i3,a)') wc, thread, ' lock id', self%id, self%level, '  test'
      !print *, thread, ' test has lock id', self%id, self%thread
    else
      !print *, thread, ' failed to acquire lock id', self%id, self%thread
    end if
  end if
  !$omp flush
  !call trace%end()
END FUNCTION test

!> =============================================================================
!> Return true only if we have rightful access to the locked item
!> =============================================================================
FUNCTION check (self) RESULT (out)
  class(lock_t):: self
  logical:: out
  integer:: thread
  !-----------------------------------------------------------------------------
  thread = omp_get_thread_num()
  out = (thread == self%thread)
END FUNCTION check

!> =============================================================================
!> Unset the associated thread number and relase the lock (in that order!)
!> =============================================================================
SUBROUTINE unset (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  integer:: thread
  real(8):: wc
  !-----------------------------------------------------------------------------
  !call trace%begin('lock_t%unset')
  thread = omp_get_thread_num()
  !$omp flush
  if (self%thread == thread) then
    !print *, self%thread, ' releases lock id', self%id
    self%level = self%level-1
    if (self%level==0) then
      call omp_unset_lock (self%lock)
      self%thread = -1
      wc = wallclock()
      if (wc < lock_verbose .and. present(label)) &
        write (io_unit%log,'(f10.6,i4,a,i6,i3,a)') wc, thread, ' lock id', self%id, self%level, ' unset at '//label
    end if
  else if (self%thread == -1) then
    print *, thread, 'ERROR: trying to release unlocked lock id', self%id
    !call trace%back
  else
    print *, thread, 'ERROR: trying to release lock of thread', self%thread, self%id
    !call trace%back
  end if
  !$omp flush
  !call trace%end()
END SUBROUTINE unset

!> =============================================================================
SUBROUTINE destroy (self)
  class(lock_t):: self
  !-----------------------------------------------------------------------------
  call omp_destroy_lock (self%lock)
  self%thread = -1
END SUBROUTINE destroy

!> =============================================================================
END MODULE omp_lock_mod
