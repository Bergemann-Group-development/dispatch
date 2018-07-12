!> =============================================================================
!> Help keep track of when all patches for a snapshot have been written out,
!> by decrementing a counter, from io%ntask to 0.  Then print a message, and
!> reset the counter for that snapshot.  Use patch%iout to determine which
!> snapshot is being written out, adding an instance of a snapshot data type
!> to the list if it does not exist, and removing it from the list when done.
!> =============================================================================
MODULE snapshot_list_mod
  USE io_mod
  USE omp_mod
  USE mpi_mod
  USE trace_mod
  USE patch_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Snapshot data type, to keep track of output
  !-----------------------------------------------------------------------------
  type, public:: snapshot_t
    type(snapshot_t), pointer:: prev=>null(), next=>null()
    integer:: iout=-1, ibuf=0, nbuf=0
    logical:: input=.false.
  end type
  type, public:: snapshot_list_t
    type(snapshot_t), pointer:: head=>null(), tail=>null()
    integer:: n=0
  contains
    procedure:: count
    procedure:: append
    procedure:: remove
  end type
  integer, save:: verbose=1
  type(snapshot_list_t), public:: snapshot_list
CONTAINS

!> =============================================================================
!> Find the relevant snapshot instance and decrement counter
!> =============================================================================
SUBROUTINE count (self, patch, snapshot, input)
  class(snapshot_list_t):: self
  class(patch_t):: patch
  type(snapshot_t), pointer:: snapshot
  logical:: input
  !-----------------------------------------------------------------------------
  ! Search for the snapshot instance
  !-----------------------------------------------------------------------------
  !$omp critical (output_cr)
  snapshot => self%head
  do while (associated(snapshot))
    if (snapshot%iout==patch%iout .and. snapshot%input.eqv.input) then
      exit
    end if
    snapshot => snapshot%next
  end do
  !-----------------------------------------------------------------------------
  ! If a snapshot instance does not exist, append one to the list
  !-----------------------------------------------------------------------------
  if (.not.associated(snapshot)) then
    call self%append (snapshot, patch, input)
  end if
  !-----------------------------------------------------------------------------
  ! If output is finished, remove from list
  !-----------------------------------------------------------------------------
  snapshot%ibuf = snapshot%ibuf-1
  if (verbose>1) print *, 'snapshot_list_t%count:', patch%iout, patch%id, snapshot%ibuf, io%ntask
  if (snapshot%ibuf==0) then
    if (verbose>1) print *, 'count:', omp%thread, mpi%rank, mpi%master, io%master
    if (verbose>0 .and. mpi%master) print '(1x,a,2x,a,i5.5,5x,"time:",1p,2g14.5,2x,9g14.5)', &
      trim(io%method)//'_io:'//merge('input ','output',input), &
      'snapshot '//trim(io%outputname), patch%iout, patch%out_next, patch%time
  end if
  !$omp end critical (output_cr)
END SUBROUTINE count

!> =============================================================================
!> Append a new snapshot record counter to the list, and initialize it 
!> =============================================================================
SUBROUTINE append (self, snapshot, patch, input)
  class(snapshot_list_t):: self
  type(snapshot_t), pointer:: snapshot
  class(patch_t):: patch
  logical:: input
  !-----------------------------------------------------------------------------
  call trace%begin ('snapshot_list_t%append')
  allocate (snapshot)
  !-----------------------------------------------------------------------------
  ! initialize the new snapshot instance counters
  !-----------------------------------------------------------------------------
  snapshot%iout = patch%iout
  snapshot%nbuf = io%ntask
  snapshot%ibuf = io%ntask
  snapshot%input = input
  !-----------------------------------------------------------------------------
  ! append to the list
  !-----------------------------------------------------------------------------
  snapshot%prev => self%tail
  if (associated(self%tail)) then
    self%tail%next => snapshot
  end if
  self%tail => snapshot
  self%n = self%n+1
  if (.not.associated(self%head)) self%head => snapshot
  call trace%end
END SUBROUTINE append

!> =============================================================================
!> Remove the snapshot instance from the list
!> =============================================================================
SUBROUTINE remove (self, snapshot)
  class(snapshot_list_t):: self
  type(snapshot_t), pointer:: snapshot
  !-----------------------------------------------------------------------------
  if (associated(snapshot%prev)) then
    snapshot%prev%next => snapshot%next
  else
    snapshot_list%head => snapshot%next
  end if
  if (associated(snapshot%next)) then
    snapshot%next%prev => snapshot%prev
  else
    snapshot_list%tail => snapshot%prev
  end if
  deallocate(snapshot)
  nullify (snapshot)
  snapshot_list%n = snapshot_list%n-1
END SUBROUTINE remove

END MODULE snapshot_list_mod
