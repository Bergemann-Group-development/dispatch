! $Id$
!===============================================================================
!===============================================================================
MODULE read_lock_mod
  USE omp_lib
  USE omp_mod
  USE omp_timer_mod
  USE io_unit_mod
  USE omp_lock_mod
  implicit none
  private
  type, public:: read_lock_t
    type(lock_t), pointer:: access => null()
    type(lock_t), pointer:: write  => null()
    logical:: on=.false.
    integer:: writer=0
    integer:: readers=0
  contains
    procedure:: init
    procedure:: set
    procedure:: unset
    procedure:: dealloc
  end type
CONTAINS

!===============================================================================
!> Initialize read and write locks
!===============================================================================
SUBROUTINE init (self, kind)
  class(read_lock_t), target:: self
  character(len=*), optional:: kind
  !.............................................................................
  if (.not.associated(self%access)) then
    allocate (self%access, self%write)
  end if
  call self%access%init (kind)
  call self%write%init (kind)
  !-----------------------------------------------------------------------------
  call self%access%set (kind)
  self%on = .true.
  call self%access%unset (kind)
END SUBROUTINE init

!===============================================================================
!> Return with the lock, waiting for writing to finish if necessary, and then
!> increment the readers counter
!===============================================================================
SUBROUTINE set (self, label)
  class(read_lock_t):: self
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  !-----------------------------------------------------------------------------
  call self%access%set (label)
  if (self%writer > 0) then
    call self%write%set (label)
    call self%write%unset (label)
  end if
  self%readers = self%readers + 1
  call self%access%unset (label)
  !-----------------------------------------------------------------------------
END SUBROUTINE set

!===============================================================================
!> Unset the associated thread number and release the lock (in that order!)
!===============================================================================
SUBROUTINE unset (self, label)
  class(read_lock_t):: self
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  !-----------------------------------------------------------------------------
  call self%access%set (label)
  self%readers = max (self%readers-1, 0)
  call self%access%unset (label)
  !-----------------------------------------------------------------------------
END SUBROUTINE unset

!===============================================================================
!> Hold on to a write_lock while a thread has write access, with the self%write
!> setting signalling to readers that they must wait for the lock before reads
!===============================================================================
SUBROUTINE set_write (self, label)
  class(read_lock_t):: self
  character(len=*), optional:: label
  integer:: readers
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  !-----------------------------------------------------------------------------
  ! Loop until the self%readers is zero, engaging the read_lock while testing
  !-----------------------------------------------------------------------------
  readers = 1
  do while (readers > 0)
    call self%access%set (label)
    readers = self%readers
    if (readers == 0) then
      call self%write%set (label)
      self%writer = 1
    end if
    call self%access%unset (label)
  end do
  !-----------------------------------------------------------------------------
END SUBROUTINE set_write

!===============================================================================
!> Engage the access lock and unset the write state
!===============================================================================
SUBROUTINE unset_write (self, label)
  class(read_lock_t):: self
  character(len=*), optional:: label
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  !-----------------------------------------------------------------------------
  call self%access%set (label)
  self%writer = 0
  call self%write%unset (label)
  call self%access%unset (label)
  !-----------------------------------------------------------------------------
END SUBROUTINE unset_write

!===============================================================================
SUBROUTINE dealloc (self)
  class(read_lock_t):: self
  !-----------------------------------------------------------------------------
  call self%access%dealloc (self%access)
  call self%write%dealloc (self%write)
END SUBROUTINE dealloc

!===============================================================================
END MODULE read_lock_mod
