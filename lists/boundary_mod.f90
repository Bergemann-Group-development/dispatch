!===============================================================================
!> Extend the boundary_t data type for guard zones with the AMR level
!===============================================================================
MODULE boundary_mod
  USE trace_mod
  USE io_mod
  USE task_mod
  USE link_mod
  implicit none
  private
  type, public, extends(nbor_t):: boundary_t
  contains
    procedure:: init
    procedure:: update
  end type
  integer:: verbose=0
CONTAINS

!===============================================================================
!> Register the source and target, and initialize the range
!===============================================================================
SUBROUTINE init (self)
  class(boundary_t):: self
  !.............................................................................
  integer:: iostat
  namelist /boundary_params/ verbose
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('boundary_t%init')
  call io%lock%set ('boundary_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, boundary_params, iostat=iostat)
    call io%namelist_warning ('boundary_params', iostat)
    write (io%namelist, boundary_params)
    first_time = .false.
  end if
  call io%lock%unset ('boundary_t%init')
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update the guard zones covered by this data type 
!===============================================================================
SUBROUTINE update (self, all, only)
  class(boundary_t):: self
  logical, optional:: all
  integer, optional:: only
  !-----------------------------------------------------------------------------
END SUBROUTINE update

END MODULE boundary_mod
