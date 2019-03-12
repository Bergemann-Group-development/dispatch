!===============================================================================
!> $Id$
!> A message list type has procedures that only need to handle the generic
!> part of task structures; the parameters, memory, and mesh structures.
!> However, it needs to be able add task to task lists, which means it must
!> be able to deal with the same tasks as the task_list_t points to.
!===============================================================================
MODULE mesg_list_mod
  USE list_mod
  USE task_list_mod
  implicit none
  type, extends(list_t):: mesg_list_t
  contains
    procedure:: probe
    procedure:: receive
    procedure:: unpack
  end type
CONTAINS

!===============================================================================
!> Probe for new MPI message and collect them on the recv_list
!===============================================================================
SUBROUTINE probe (recv_list)
  class(mesg_list_t):: recv_list
END SUBROUTINE probe

!===============================================================================
!> Check for complete messages in recv_list and move them to mesg_list
!===============================================================================
SUBROUTINE receive (recv_list, unpack_list)
  class(mesg_list_t)::  recv_list, unpack_list
END SUBROUTINE receive

!===============================================================================
!> Unpack messages in the unpack_list to the appropriate virtual patch
!===============================================================================
SUBROUTINE unpack (unpack_list, task_list)
  class(mesg_list_t):: unpack_list
  type(task_list_t):: task_list
END SUBROUTINE unpack

END MODULE mesg_list_mod
