! $Id$
!===============================================================================
!> This module implements a multi-thread-read / singl-thread-write lock.
!> It uses a single nested lock and an atomically updated count of the number
!> of readers.  A write request starts with acquiring the lock, which blocks
!> new read requests, and then spin-waits for the count of readers to reach
!> zero.  A read request engages the lock briefly, allowing the request to be
!> blocked during writes, and then increments the reader count.
!===============================================================================
MODULE omp_lock_mod
  USE omp_lib
  USE omp_mod
  USE omp_timer_mod
  USE io_unit_mod
  USE org_lock_mod
  implicit none
  private
  type, public, extends(org_lock_t):: lock_t
  contains
    procedure:: init
    procedure:: set
    procedure:: set_read
    procedure:: unset
    procedure:: unset_read
    procedure:: destroy
  end type
  integer:: verbose=0
  type(lock_t), save, public:: omp_lock
CONTAINS

!===============================================================================
!> Initialize read and write locks
!===============================================================================
SUBROUTINE init (self, kind, id)
  class(lock_t), target:: self
  character(len=*), optional:: kind
  integer, optional:: id
  !.............................................................................
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Make sure the first lock initialization reads parameters
  !-----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      call self%org_lock_t%init (kind, id)
      call self%org_lock_t%init_unit ()
      first_time = .false.
    end if
    !$omp end critical (input_cr)
  else
    call self%org_lock_t%init (kind, id)
  end if
END SUBROUTINE init

!===============================================================================
!> To get read access, acquire the write lock briefly and increase the reader
!> count.
!===============================================================================
SUBROUTINE set_read (self, label, itimer)
  class(lock_t):: self
  character(len=*), optional:: label
  integer, optional:: itimer
  integer:: writer
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  if (self%fall_through) then
    call self%org_lock_t%set (label, itimer)
    return
  end if
  !-----------------------------------------------------------------------------
  ! Use the access to increment the number of readers; if a writer is active
  ! the access lock is being held until the write is done
  !-----------------------------------------------------------------------------
  call self%org_lock_t%set ('R'//label)
  !$omp atomic update
  self%readers = self%readers + 1
  call self%org_lock_t%unset ('R'//label)
  !-----------------------------------------------------------------------------
END SUBROUTINE set_read

!===============================================================================
!> Unset the associated thread number and release the lock (in that order!)
!===============================================================================
SUBROUTINE unset_read (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  if (self%fall_through) then
    call self%org_lock_t%unset (label)
    return
  end if
  !-----------------------------------------------------------------------------
  !$omp atomic update
  self%readers = self%readers - 1
  call self%unset_log (label)
  !-----------------------------------------------------------------------------
END SUBROUTINE unset_read

!===============================================================================
!> Hold on to a write_lock while a thread has write access, with the self%write
!> setting signalling to readers that they must wait for the lock before reads
!===============================================================================
SUBROUTINE set (self, label, itimer)
  class(lock_t):: self
  character(len=*), optional:: label
  integer, optional:: itimer
  integer:: readers, oreaders
  real(8):: last, now
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  if (self%fall_through) then
    call self%org_lock_t%set (label, itimer)
    return
  end if
  !-----------------------------------------------------------------------------
  ! Acquire the write lock, which will hold back new read attempts
  !-----------------------------------------------------------------------------
  call self%org_lock_t%set (label)
  !-----------------------------------------------------------------------------
  ! Loop until self%readers is zero, engaging the access lock while testing
  !-----------------------------------------------------------------------------
  readers = 1
  oreaders = 0
  last = omp_get_wtime()
  do while (readers > 0)
    !$omp atomic read
    readers = self%readers
    if (omp%nthreads==1 .and. readers>0) then
      write (*,*) 'lock_t%set_write: inside own read lock'
      stop
    end if
    if (readers /= oreaders) then
      write (self%unit,'(a,i6,2i4)') &
        'set: id,thread,readers=', self%id, omp%thread, readers
      oreaders = readers
    end if
    now = omp_get_wtime()
    if (now > last + 1d0) then
      write (self%unit,'(a,i6,2i4)') &
        'set: time-out, id,thread,readers=', self%id, omp%thread, readers
      last = now
    end if
    if (readers == 0) &
      exit
  end do
  !-----------------------------------------------------------------------------
END SUBROUTINE set

!===============================================================================
!> Unset the write lock, which has held back read requests in the meantime
!===============================================================================
SUBROUTINE unset (self, label)
  class(lock_t):: self
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  if (self%fall_through) then
    call self%org_lock_t%unset (label)
    return
  end if
  !-----------------------------------------------------------------------------
  call self%org_lock_t%unset (label)
END SUBROUTINE unset

!===============================================================================
SUBROUTINE destroy (self)
  class(lock_t):: self
  !-----------------------------------------------------------------------------
  if (self%fall_through) then
    call self%org_lock_t%destroy
    return
  end if
  if (self%initialized) then
    call self%org_lock_t%destroy ()
    self%initialized = .false.
  end if
END SUBROUTINE destroy

!===============================================================================
END MODULE omp_lock_mod
