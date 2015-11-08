!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  implicit none
  integer, parameter:: ntime=7
  !-----------------------------------------------------------------------------
  ! Basic task info
  !-----------------------------------------------------------------------------
  type, abstract:: task_t
    integer:: id = 0
    integer:: rank = 0
    real(8):: time = 0d0
    real(8):: dtime = 0d0
    real(8):: position(3) = 0d0
    real(8):: velocity(3) = 0d0
    real(8):: size(3) = 0d0
    integer:: status = 0
  end type
  !-----------------------------------------------------------------------------
  ! Timeslot info added
  !-----------------------------------------------------------------------------
  type, abstract, extends(task_t):: tasks_t
    real(8):: t(ntime)
    real(8):: dt(ntime)
    integer:: it = 0
    contains
    procedure(init)      , deferred:: init
    procedure(allocate)  , deferred:: allocate
    procedure(deallocate), deferred:: deallocate
    procedure(update)    , deferred:: update
    procedure(print)     , deferred:: print
  end type
  abstract interface 
    SUBROUTINE init (self)
      import tasks_t
      class(tasks_t):: self
    END SUBROUTINE
    SUBROUTINE allocate (self)
      import tasks_t
      class(tasks_t):: self
    END SUBROUTINE
    SUBROUTINE deallocate (self)
      import tasks_t
      class(tasks_t):: self
    END SUBROUTINE
    SUBROUTINE update (self)
      import tasks_t
      class(tasks_t):: self
    END SUBROUTINE
    SUBROUTINE print (self)
      import tasks_t
      class(tasks_t):: self
    END SUBROUTINE
  end interface
PRIVATE
PUBLIC ntime, task_t, tasks_t
END MODULE
