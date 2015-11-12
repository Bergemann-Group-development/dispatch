!> =============================================================================
!> $Id$
!> This object extends a task with procedures to pack the task into a bufer and
!> send it to another rank.  At that rank, other procedures handle the unpacking
!> and further processing.
!> =============================================================================
MODULE message_mod
  USE task_mod
  USE list_mod
  type, extends(tasks_t):: message_t
    integer(1), dimension(:), pointer:: buffer
    integer:: req
    integer:: msg
    contains
    procedure:: pack
    procedure:: send_to
    procedure:: probe
    procedure:: receive
    procedure:: move_to
    procedure:: unpack
    procedure:: allocate   => void
    procedure:: deallocate => void
    procedure:: print      => void
    procedure:: update     => void
    procedure:: init       => void
    procedure:: check      => void
  end type
CONTAINS

SUBROUTINE pack (self)
  class (message_t):: self
  class (task_t), pointer:: task
END SUBROUTINE

SUBROUTINE send_to (self, rank)
  class (task_t), pointer:: task
  class (message_t):: self
  class (list_t), pointer:: list
  integer:: rank
END SUBROUTINE

SUBROUTINE probe (self, list)
  class (message_t):: self
  class (list_t), pointer:: list
END SUBROUTINE

SUBROUTINE receive (self, list)
  class (message_t):: self
  class (list_t), pointer:: list
END SUBROUTINE

SUBROUTINE move_to (self, list1, list2)
  class (message_t):: self
  class (list_t), pointer:: list1, list2
END SUBROUTINE

SUBROUTINE unpack (self, list)
  class (message_t):: self
  class (list_t), pointer:: list
END SUBROUTINE

SUBROUTINE allocate (self)
  class (message_t):: self
END SUBROUTINE

SUBROUTINE print (self)
  class (message_t):: self
END SUBROUTINE

SUBROUTINE void (self)
  class (message_t):: self
END SUBROUTINE

END MODULE
