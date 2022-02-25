!===============================================================================
!> Module storing the filenames and unit numbers of binary and metadata files.
!>
!> Each snapshot instance keeps track of the snapshot number, the record number,
!> and the unit numbers of data and metadata files.  At most %max snapshots are
!> kept (set to 5 mainly to test the mechanism; both ifort and gfortran allow
!> at least 1000 newunit requests).
!===============================================================================
MODULE snapshots_mod
  USE io_mod
  USE trace_mod
  USE dll_mod
  implicit none
  private
  type, public:: snapshot_t
    logical:: initialized = .false.
    integer:: iout = -1
    integer:: record = 0
    integer:: data_unit = 0
    integer:: meta_unit = 0
  contains
  end type
  type, public, extends(dll_t):: snapshots_t
    integer:: max = 10
  contains
    procedure:: get
  end type
  integer, save:: verbose = 0
  type(snapshots_t), public:: snapshots
CONTAINS

!===============================================================================
!> Search for a snapshot, returning a non-null pointer if successful
!===============================================================================
FUNCTION get (self, iout) RESULT (found)
  class(snapshots_t)         :: self
  integer                    :: iout
  class(snapshot_t), pointer :: found
  !.............................................................................
  class(dll_node_t), pointer :: item
  class(*), pointer          :: ptr
  !-----------------------------------------------------------------------------! Initialize?
  ! Initialize?
  !-----------------------------------------------------------------------------
  call trace%begin ('snapshots_t%get')
  if (.not. self%initialized) then
    call io%lock%set
    if (.not. self%initialized) then
      call self%init
      self%initialized = .true.
    end if
    call io%lock%unset
  end if
  !-----------------------------------------------------------------------------
  ! Search for matching snapshot
  !-----------------------------------------------------------------------------
  nullify (found)
  item => self%head
  do while (associated(item))
    ptr => item%car
    select type (ptr)
    class is (snapshot_t)
      found => ptr
      if (found%iout == iout) then
        call trace%end ()
        return
      end if
    end select
    item => item%next
  end do
  !-----------------------------------------------------------------------------
  ! Not found, check first if it's time to remove an old one
  !-----------------------------------------------------------------------------
  if (self%n > self%max) then
    item => self%head
    ptr => item%car
    select type (ptr)
    class is (snapshot_t)
      found => ptr
      if (verbose > 0) &
        write (stdout,*) 'snapshots_t%get: deallocating snapshot', found%iout
      if (found%data_unit /= 0) close (found%data_unit)
      if (found%meta_unit /= 0) close (found%meta_unit)
      deallocate (found)
    end select
    call self%remove (item)
  end if
  !-----------------------------------------------------------------------------
  ! Allocate a new one
  !-----------------------------------------------------------------------------
  allocate (item, found)
  call self%append (item)
  item%car => found
  found%iout = iout
  if (verbose > 1) &
    write (stdout,*) 'snapshots_t%get: allocating snapshot', found%iout
  call trace%end ()
END FUNCTION get

END MODULE snapshots_mod
