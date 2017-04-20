!> =============================================================================
!> The gpatch_t data type extends the patch_t with a generic guard zone dnload
!> procedure, which covers the need for all simple, Cartesian mesh based 
!> solvers, and eliminates the need for a dnload procedure in each one of them.
!> It also adds restart functionality when initializing patches.
!> =============================================================================
MODULE gpatch_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE download_mod
  USE initial_mod
  USE kinds_mod
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
    type(initial_t):: initial
  contains
    procedure:: init
    procedure:: dnload
  end type
  real(8), save:: restart_time = 0d0
CONTAINS

!> =============================================================================
!> Add restart functionality to patch initialization
!> =============================================================================
SUBROUTINE init (self)
  class(gpatch_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: ff
  logical:: ok
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%init')
  call self%patch_t%init
  call self%initial%init
  !-----------------------------------------------------------------------------
  ! Try reading a restart snapshot.
  !-----------------------------------------------------------------------------
  call self%input (ok)
  if (ok) then
    restart_time = self%time
  !-----------------------------------------------------------------------------
  ! If the restart parameter is not set, or the file does not exist, initialize
  ! the patch from scratch, resetting other necessary counters only if positive
  ! restart_time indicates that restart files were read for previous patches.
  !-----------------------------------------------------------------------------
  else
    ff => self%mem(:,:,:,:,1,1)
    call self%initial%condition (self%mesh, ff)
  end if
  !-----------------------------------------------------------------------------
  ! If restart succeeded for some patch, use it to set time and output cadence.
  ! FIXME: This may fail in some corner cases when running on multiple MPI ranks,
  ! if one rank happens to handle only patches that do not exist yet.
  !-----------------------------------------------------------------------------
  if (restart_time > 0d0) then
    self%time = restart_time
    self%t(self%it) = self%time
    self%out_next = (int(self%time/self%out_time)+1)*self%out_time
    self%iout = self%restart+1
  end if
  call trace%end
END SUBROUTINE

!> =============================================================================
!> Add default guard zone download functionality.  This may be overloaded in
!> experiment_mod, if more specialized guard zone loading is required.
!> =============================================================================
SUBROUTINE dnload (self)
  class(gpatch_t):: self
  call download%download_link (self%link, all_cells=.false.)
END SUBROUTINE

END MODULE gpatch_mod
