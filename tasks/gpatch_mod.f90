!> =============================================================================
!> The gpath_t layer now essentially only handles restarts
!> =============================================================================
MODULE gpatch_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE patch_mod
  USE initial_mod
  USE kinds_mod
  USE mesh_mod
  USE units_mod
  USE timer_mod
  USE omp_timer_mod
  USE bits_mod
  USE data_io_mod
  USE download_mod
  USE data_hub_mod
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
    type(initial_t)  :: initial
    type(data_hub_t) :: data_hub
    real(8):: d0
  contains
    procedure:: input
    procedure:: init
    procedure:: output
    procedure:: counter_update
    procedure:: dnload
  end type
  real(8), save :: restart_time = 0d0
  logical, save :: detailed_timer=.false.
  integer, save :: verbose=0
  integer, save :: order=1
CONTAINS

!> =============================================================================
!> Interface routine to the data_io input procedure, signalling if the patch
!> could be read, or not
!> =============================================================================
SUBROUTINE input (self, ok)
  class(gpatch_t):: self
  logical:: ok
  !-----------------------------------------------------------------------------
  call data_io%input (self, ok)
END SUBROUTINE input

!> =============================================================================
!> Interface routine to the data_io output procedure
!> =============================================================================
SUBROUTINE update (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call data_io%output (self)
END SUBROUTINE update

!> =============================================================================
!> Interface routine to the data_io output procedure
!> =============================================================================
SUBROUTINE output (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call data_io%output (self)
END SUBROUTINE output

!> =============================================================================
!> Add restart functionality to patch initialization.  This %init routine is
!> called for both internal and virtual patches, which both can be initialized
!> fully by the self%ininitial%condition() procedure. When instead calling the
!> self%input() procedure, we need to make a decision:  Either we must make sure
!> that the data_io%input() routine can handle off-rank reading, or we need to
!> nevertheless return ok=.true., and instead rely on a negative time to force
!> any task that relies on the patch for guard zone values to wait until data
!> has been recieved with MPI from the owner rank.
!> =============================================================================
SUBROUTINE init (self)
  class(gpatch_t):: self
  real(kind=KindScalarVar), pointer :: ff(:,:,:,:)
  character(len=64)                 :: filename
  logical                           :: ok
  real(kind=KindScalarVar), pointer :: d(:,:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%init')
  call self%initial%init (self%kind, real(self%gamma))
  call self%patch_t%init
  !-----------------------------------------------------------------------------
  ! Try reading a restart snapshot.
  !-----------------------------------------------------------------------------
  call self%input (ok)
  if (ok) then
    restart_time = self%time
  else 
    !---------------------------------------------------------------------------
    ! Make sure the default time for patches with IC values is zero, to ensure 
    ! boundary tasks are not considered ready for updating until their virtual
    ! task nbors have been updated.
    !---------------------------------------------------------------------------
    restart_time = 0.0
    self%time = 0.0
    self%mem = 0.0
    ff => self%mem(:,:,:,:,1,1)
    call self%initial%condition (self%mesh, ff, self%idx)
  end if
  !-----------------------------------------------------------------------------
  ! If restart succeeded for some patch, use it to set time and output cadence.
  ! FIXME: This may fail in some corner cases when running on multiple MPI ranks,
  ! if one rank happens to handle only patches that do not exist yet.
  !-----------------------------------------------------------------------------
  if (restart_time > 0d0) then
    self%time = restart_time
    self%t(self%it) = self%time
    self%out_next = (int(self%time/io%out_time+1.0e-4)+1)*io%out_time ! add an `eps` to avoid round-off errors
    self%iout = self%restart+1
    !$omp master
    write (filename,'(a,i5.5,"/")') trim(io%outputname), self%iout
    if (io%iodir/=self%iout) call os%mkdir (trim(filename))
    io%iodir = self%iout
    !$omp end master
  end if
  call trace%end()
END SUBROUTINE init

!> =============================================================================
!> Increment the cell udpate counter.  This may be overloaded in more complex
!> solvers, which prefer a different way of counting.
!> =============================================================================
SUBROUTINE counter_update (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%counter_update')
  !$omp atomic
  timer%n_update = timer%n_update + product(self%n)
  call trace%end
END SUBROUTINE counter_update

!> =============================================================================
!> Add default guard zone download functionality.  This may be overloaded in
!> experiment_mod, if more specialized guard zone loading is required.
!> =============================================================================
SUBROUTINE dnload (self, only)
  class(gpatch_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('extras_t%dnload')
  if (self%use_data_hub) then
    call self%data_hub%update (self%link, only=only)
  else
    call download%download_link (self%link, only=only)
  end if
  call trace%end()
END SUBROUTINE dnload

END MODULE gpatch_mod
