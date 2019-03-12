! $Id$
!===============================================================================
!> Using trace_mod causes a dependency loop, and can only be compiled if/when
!> the code is already successfully compiled.  It should only be used to trace
!> problems with locks.
!===============================================================================
MODULE omp_lock_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lib
!define TRACE
#ifdef TRACE 
  USE trace_mod
#endif
  implicit none
  private
  type, public:: lock_t
    integer:: id=-1
    integer(kind=omp_nest_lock_kind):: lock
    integer:: thread=-1
    character(len=4):: kind = 'void'
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
  real:: verbose = 0.0
  integer, save:: id_save=-1
CONTAINS

!===============================================================================
!> Initialize a lock, taking care not to do it more than once
!===============================================================================
SUBROUTINE init (self, kind, id)
  class(lock_t):: self
  integer, optional:: id
  character(len=4), optional:: kind
  logical omp_in_parallel
  integer:: iostat
  logical:: set_id
  logical:: first_time=.true.
  namelist /lock_params/ verbose
  !-----------------------------------------------------------------------------
  if (self%id == -1) then
#ifdef TRACE 
    call trace%begin('lock_t%init')
#endif
    !$omp critical (lock_cr)
    if (self%id == -1) then
      if (first_time) then
        first_time = .false.
        rewind (io_unit%input)
        read (io_unit%input, lock_params, iostat=iostat)
        write (io_unit%output, lock_params)
      end if
      call omp_init_nest_lock (self%lock)
      if (present(id)) then
        self%id = id
        id_save = max(id_save,id)
        set_id = .true.
      else
        id_save = id_save+1
        self%id = id_save
        set_id = .false.
      end if
      if (present(kind)) self%kind = kind
      if (verbose > 0.0) &
        print *, omp_get_thread_num(), 'initialized lock kind, id ', &
          self%kind, self%id, set_id
    end if
    !$omp end critical (lock_cr)
#ifdef TRACE 
    call trace%end()
#endif
  end if
END SUBROUTINE init

!===============================================================================
!> Return with the lock, waiting if necessary, and creating the lock if needed
!===============================================================================
SUBROUTINE set (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  character(len=32):: filename
  integer:: thread, othread, level
  real(8):: wc
  !-----------------------------------------------------------------------------
  call self%init
#ifdef TRACE 
  call trace%begin('lock_t%set')
#endif
  thread = omp_get_thread_num()
  call omp_set_nest_lock (self%lock)
  self%thread = thread
  self%level = self%level+1
  !write(io_unit%output,*) &
  !  'lock%  set: kind, id, thread, level = ', self%kind, self%id, self%thread, self%level
  !$omp flush
  wc = wallclock()
  if (wc < verbose) then
    if (present(label)) then
      write (io_unit%output,'(f10.6,i4,a,a4,i6,i3,a)') wc, thread, &
        ' lock kind, id, level = ', self%kind, self%id, self%level, '   set at '//label
    else
      write (io_unit%output,'(f10.6,i4,a,a4,i6,i3,a)') wc, thread, &
        ' lock kind, id, level = ', self%kind, self%id, self%level, '   set'
    end if
  end if
#ifdef TRACE 
  call trace%end()
#endif
END SUBROUTINE set

!===============================================================================
!> Return true if we already have the lock, or if we get it now
!===============================================================================
FUNCTION test (self) RESULT (out)
  class(lock_t):: self
  logical:: out
  integer:: thread
  real(8):: wc
  !-----------------------------------------------------------------------------
#ifdef TRACE 
  call trace%begin('lock_t%test')
#endif
  thread = omp_get_thread_num()
  !print *, thread, ' set lock id, thread', self%id, self%thread
  !$omp flush
  !print *, thread, ' testing lock id', self%id, self%thread
  if (thread == self%thread) then
    out = .true.
    !print *, thread, ' already has lock id', self%id
  else
    !print *, thread, ' omp_test_lock for id', self%id
    out =  omp_test_nest_lock (self%lock)
    if (out) then
      self%thread = thread
      self%level = self%level+1
      wc = wallclock()
      if (wc < verbose) &
        write (io_unit%log,'(f10.6,i4,a,i6,i3,a)') wc, thread, ' lock id', self%id, self%level, '  test'
      !print *, thread, ' test has lock id', self%id, self%thread
    else
      !print *, thread, ' failed to acquire lock id', self%id, self%thread
    end if
  end if
  !$omp flush
#ifdef TRACE 
  call trace%end()
#endif
END FUNCTION test

!===============================================================================
!> Return true only if we have rightful access to the locked item
!===============================================================================
FUNCTION check (self) RESULT (out)
  class(lock_t):: self
  logical:: out
  integer:: thread
  !-----------------------------------------------------------------------------
  thread = omp_get_thread_num()
  out = (thread == self%thread)
END FUNCTION check

!===============================================================================
!> Unset the associated thread number and relase the lock (in that order!)
!===============================================================================
SUBROUTINE unset (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  integer:: thread
  real(8):: wc
  !-----------------------------------------------------------------------------
#ifdef TRACE 
  call trace%begin('lock_t%unset')
#endif
  thread = omp_get_thread_num()
  self%level = self%level-1
  if (self%level==0) then
    self%thread = -1
  end if
  !$omp flush
  call omp_unset_nest_lock (self%lock)
  !write(io_unit%output,*) &
  !  'lock%unset: kind, id, thread, level = ', self%kind, self%id, self%thread, self%level
  wc = wallclock()
  if (wc < verbose) then
    if (present(label)) then
      write (io_unit%output,'(f10.6,i4,a,a4,i6,i3,a)') wc, thread, &
        ' lock kind, id, level = ', self%kind, self%id, self%level, ' unset at '//label
    else
      write (io_unit%output,'(f10.6,i4,a,a4,i6,i3,a)') wc, thread, &
        ' lock kind, id, level = ', self%kind, self%id, self%level, ' unset'
    end if
  end if
#ifdef TRACE 
  call trace%end()
#endif
END SUBROUTINE unset

!===============================================================================
SUBROUTINE destroy (self)
  class(lock_t):: self
  !-----------------------------------------------------------------------------
  call omp_destroy_nest_lock (self%lock)
  self%thread = -1
END SUBROUTINE destroy

!===============================================================================
END MODULE omp_lock_mod
