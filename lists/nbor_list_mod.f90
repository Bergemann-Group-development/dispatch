!> =============================================================================
!> Handle neighbor lists for tasks
!> =============================================================================
MODULE nbor_list_mod
  USE bits_mod
  USE mpi_mod
  USE item_mod
  USE task_mod
  USE list_mod
  USE simple_list_mod
  implicit none
  type, extends(simple_list_t):: nbor_list_t
  contains
    procedure:: init
  end type
PRIVATE
PUBLIC  nbor_list_t
CONTAINS

!> =============================================================================
!> Construct the nbor_list for a task, by running through a task_list.  Note that
!> making this a task procedure rather than an nbor list procedure does not work,
!> since an task definition cannot use a task_list
!> =============================================================================
SUBROUTINE init (self, name, list, item)
  class(nbor_list_t):: self
  character(len=*):: name
  class(list_t), pointer, optional:: list
  class(item_t), pointer, optional:: item
  type(link_t), pointer:: link
  !.............................................................................
  call self%simple_list_t%init (name)             ! this purges the list
  link => list%head                               ! start on the list
  do while (associated(link))                     ! until end
    item => link%item                             ! get the task point
    select type (item)                            ! determine the type
    class is (task_t)                             ! select only items
      if (item%overlaps (item)) then
        if (item%rank == mpi%rank) then
          call item%clear (bits%virtual)          ! item is not virtual
          call item%set (bits%internal)           ! task is internal
        else
          call item%set (bits%virtual)            ! item is virtual
          call item%set (bits%boundary)           ! task is boundary
        end if
!! FIXME        call self%append (item)                   ! append to self%nbor_list
      end if
    end select
    link => link%next                             ! next list link
  end do
END SUBROUTINE init

END MODULE nbor_list_mod
