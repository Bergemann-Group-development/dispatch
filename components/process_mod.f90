!===============================================================================
!> Data type to keep and maintain information on MPI processes.  The process
!> data type maintains a list of nbor MPI ranks, which is used to send and recv.
!>
!> To send a buffer to all nbors, repeatedly do:
!>
!> call process%send (buffer)
!>
!> To process incoming buffers (with "call unpack(buffer)"), repeatedly do
!>
!> call process%recv (buffer, unpack)
!>
!===============================================================================
MODULE nbor_list_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE mpi_comm_mod
  USE link_mod
  USE task_list_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! process nbor
  !-----------------------------------------------------------------------------
  type, public:: nbor_t
    type(nbor_t), pointer:: next => null()
    type(mpi_comm_t):: out, in
    integer:: rank
    logical:: present
  end type
  !-----------------------------------------------------------------------------
  ! process nbor list
  !-----------------------------------------------------------------------------
  type, public:: nbor_list_t
    type(nbor_t), pointer:: head => null()
    type(nbor_t), pointer:: tail => null()
  contains
    procedure:: update
    procedure:: add_by_rank
  end type
CONTAINS

!===============================================================================
!> Initialize a neighbor list
!===============================================================================
SUBROUTINE update (self, task_list)
  class(nbor_list_t):: self
  type(task_list_t):: task_list
  !.............................................................................
  type(link_t), pointer:: link
  type(nbor_t), pointer:: nbor, next, prev
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Clear the present flag
  !-----------------------------------------------------------------------------
  call trace%begin ('nbor_list_t%update')
  nbor => self%head
  do while (associated(nbor))
    nbor%present = .false.
    nbor => nbor%next
  end do
  !-----------------------------------------------------------------------------
  ! Add to the nbor_list and set present flags
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    if (link%task%rank /= mpi%rank) &
      call self%add_by_rank (link%task%rank)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Prune the nbor_list = remove ranks that are no longer present
  !-----------------------------------------------------------------------------
  nullify (prev)
  nbor => self%head
  if (first_time.and.io%master) &
    print *, 'nbors to rank', mpi%rank, ':'
  do while (associated(nbor))
    next => nbor%next
    if (.not.nbor%present) then
      if (associated(prev)) then
        prev%next => next
      else
        self%head => next
      end if
      deallocate (nbor)
    end if
    if (first_time.and.io%master) &
      print *, 'nbor rank:', nbor%rank
    nbor => next
  end do
  first_time = .false.
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Add an nbor to the list, keeping it sorted by rank
!===============================================================================
SUBROUTINE add_by_rank (self, rank)
  class(nbor_list_t):: self
  integer:: rank
  !.............................................................................
  type(nbor_t), pointer:: this, nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('nbor_list_t%add_by_rank')
  this => self%head
  if (associated(this)) then
    !---------------------------------------------------------------------------
    ! If nbor%rank is less than self%head%rank, place it at the head
    !---------------------------------------------------------------------------
    if (rank < this%rank) then
      nbor => new_nbor (rank)
      nbor%next => this
      self%head => nbor
      call trace%end()
      return
    else
      do while (associated(this))
        !-----------------------------------------------------------------------
        ! If not at the tail, and nbor%rankd is in between this and the next
        !-----------------------------------------------------------------------
        if (rank == this%rank) then
          this%present = .true.
          return
        else if (associated(this%next)) then
          if (rank > this%rank .and. rank < this%next%rank) then
            nbor => new_nbor (rank)
            nbor%next => this%next
            this%next => nbor
            call trace%end()
            return
          end if
        !-----------------------------------------------------------------------
        ! At the tail
        !-----------------------------------------------------------------------
        else if (rank > this%rank) then
          nbor => new_nbor (rank)
          this%next => nbor
          self%tail => nbor
          call trace%end()
          return
        end if
        this => this%next
        call trace%end()
        return
      end do
    end if
  !-----------------------------------------------------------------------------
  ! If self%head is not associated, nbor is the first list elemend to be added
  !-----------------------------------------------------------------------------
  else
    nbor => new_nbor (rank)
    self%head => nbor
    self%tail => nbor
  end if
  call trace%end()
contains
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  function new_nbor (rank)
    integer:: rank
    type(nbor_t), pointer:: new_nbor
    !---------------------------------------------------------------------------
    allocate (new_nbor)
    new_nbor%rank = rank
    new_nbor%present = .true.
  end function
END SUBROUTINE add_by_rank

END MODULE nbor_list_mod

!===============================================================================
!> Data type to keep and maintain information on MPI processes
!===============================================================================
MODULE process_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE link_mod
  USE task_list_mod
  USE nbor_list_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! process
  !-----------------------------------------------------------------------------
  type, public:: process_t
    integer:: rank
    type(nbor_list_t):: nbor_list
  contains
    procedure:: update
    procedure:: send
    procedure:: recv
  end type
CONTAINS

!===============================================================================
!> Initialize a process_t instance
!===============================================================================
SUBROUTINE update (self, task_list)
  class(process_t):: self
  type(task_list_t):: task_list
  !.............................................................................
  call self%nbor_list%update (task_list)
  self%rank = mpi%rank
END SUBROUTINE update

!===============================================================================
!> Send a buffer to all nbors.  When a send is complete, renew the buffer, and
!> resend.
!===============================================================================
SUBROUTINE send (self, buffer)
  class(process_t):: self
  integer, pointer:: buffer(:)
  type(nbor_t), pointer:: nbor
  !.............................................................................
  nbor => self%nbor_list%head
  do while (associated(nbor))
    if (.not.allocated(nbor%out%buffer)) then
      allocate (nbor%out%buffer(size(buffer)))
      nbor%out%buffer = buffer
      call nbor%out%send (nbor%rank, 111)
    end if
    if (nbor%out%sent ()) then
      nbor%out%buffer = buffer
      call nbor%out%send (nbor%rank, 111)
    end if
    nbor => nbor%next
  end do
END SUBROUTINE send

!===============================================================================
!> Receive and process a buffer from all nbors
!===============================================================================
SUBROUTINE recv (self, buffer, process)
  class(process_t):: self
  integer, pointer:: buffer(:)
  procedure(template):: process
  !.............................................................................
  type(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  ! Loop over process nbors, asking for input packages, and process when recv'd
  !-----------------------------------------------------------------------------
  nbor => self%nbor_list%head
  do while (associated(nbor))
    if (.not.allocated(nbor%in%buffer)) then
      allocate (nbor%in%buffer(size(buffer)))
    end if
    if (nbor%in%recv (nbor%rank, 111)) then
      call process (nbor%in%buffer)
    end if
    nbor => nbor%next
  end do
END SUBROUTINE recv

!===============================================================================
!> Template processing procedure
!===============================================================================
SUBROUTINE template (buffer)
  integer:: buffer(:)
END SUBROUTINE template

END MODULE process_mod
