!===============================================================================
!> The gpath_t layer now essentially only handles restarts
!===============================================================================
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
  USE guard_zones_mod
  USE download_mod
  USE data_hub_mod
  USE validate_mod
  USE list_mod
  USE task_mod
  implicit none
  private
  type, public, extends(patch_t):: gpatch_t
    type(initial_t):: initial
    type(data_hub_t):: data_hub
    class(list_t), pointer:: task_list => null()
    real(8):: d0
  contains
    procedure, nopass:: cast2gpatch
    procedure:: init_task_list
    procedure:: update
    procedure:: input
    procedure:: output
    procedure:: init
    procedure:: counter_update
    procedure:: dnload
    procedure:: courant_condition
  end type
  real(8), save :: restart_time = 0d0
  logical, save :: detailed_timer=.false.
  integer, save :: verbose=0
  integer, save :: order=1
  type(gpatch_t), public:: gpatch
CONTAINS

!===============================================================================
!> Cast a generic task_t to patch_t
!===============================================================================
FUNCTION cast2gpatch (task) RESULT(gpatch)
  class(task_t), pointer:: task
  class(gpatch_t), pointer:: gpatch
  !.............................................................................
  select type (task)
  class is (gpatch_t)
  gpatch => task
  class default
  nullify(gpatch)
  call io%abort ('gpatch_t%cast: failed to cast a task to patch_t')
  end select
END FUNCTION cast2gpatch

!===============================================================================
!> Make a copy of the task list pointer
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(gpatch_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  call trace%begin ('gpatch_t%init_task_list')
  self%task_list => task_list
  call trace%end()
END SUBROUTINE init_task_list

!===============================================================================
!> Interface allowing intercept
!===============================================================================
SUBROUTINE update (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call self%patch_t%update
END SUBROUTINE update

!===============================================================================
!> Interface routine to the data_io input procedure, signalling if the patch
!> could be read, or not
!===============================================================================
SUBROUTINE input (self, ok)
  class(gpatch_t):: self
  logical:: ok
  !-----------------------------------------------------------------------------
  call data_io%input (self, ok)
END SUBROUTINE input

!===============================================================================
!> Interface routine to the data_io output procedure
!===============================================================================
SUBROUTINE output (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call data_io%output (self)
END SUBROUTINE output

!===============================================================================
!> Add restart functionality to patch initialization.  This %init routine is
!> called for both internal and virtual patches, which both can be initialized
!> fully by the self%ininitial%condition() procedure. When instead calling the
!> self%input() procedure, we need to make a decision:  Either we must make sure
!> that the data_io%input() routine can handle off-rank reading, or we need to
!> nevertheless return ok=.true., and instead rely on a negative time to force
!> any task that relies on the patch for guard zone values to wait until data
!> has been recieved with MPI from the owner rank.
!===============================================================================
SUBROUTINE init (self)
  class(gpatch_t):: self
  real(kind=KindScalarVar), pointer :: ff(:,:,:,:)
  character(len=64)                 :: filename
  logical                           :: ok
  real(kind=KindScalarVar), pointer :: d(:,:,:)
  !-----------------------------------------------------------------------------
  if (self%is_set(bits%frozen)) &
    return
  print *, 'gpath_t%init: time =', self%time
  call trace%begin ('gpatch_t%init')
  call self%initial%init (self%kind, real(self%gamma))
  call data_io%init (self)
  call guard_zones%init
  !-----------------------------------------------------------------------------
  ! Try reading a restart snapshot.
  !-----------------------------------------------------------------------------
  call self%input (ok)
  if (ok) then
    restart_time = self%time
    self%istep = 0
  else 
    !---------------------------------------------------------------------------
    ! Make sure the default time for patches with IC values is zero, to ensure 
    ! boundary tasks are not considered ready for updating until their virtual
    ! task nbors have been updated, in evolved cases.
    !---------------------------------------------------------------------------
    restart_time = 0.0
    self%mem = 0.0
    ff => self%mem(:,:,:,:,1,1)
    call self%initial%condition (self%mesh, ff, self%idx)
    !---------------------------------------------------------------------------
    ! patch_t%setup sets time = -1, to prevent downloading until nbors have been
    ! initialized, and their time has been set to 0.0 here
    !---------------------------------------------------------------------------
    self%time = 0.0
    self%t = 0.0
  end if
  !-----------------------------------------------------------------------------
  ! If restart succeeded for some patch, use it to set time and output cadence.
  ! FIXME: This may fail in some corner cases when running on multiple MPI ranks,
  ! if one rank happens to handle only patches that do not exist yet.
  !-----------------------------------------------------------------------------
  if (restart_time > 0d0) then
    call self%lock%set ('gpatch_t')
    self%time = restart_time
    self%t(self%it) = self%time
    self%out_next = (int(self%time/io%out_time+1.0e-4)+1)*io%out_time ! add an `eps` to avoid round-off errors
    self%iout = self%restart+1
    !$omp master
    write (filename,'(a,i5.5,"/")') trim(io%outputname), self%iout
    if (io%iodir/=self%iout) call os%mkdir (trim(filename))
    io%iodir = self%iout
    !$omp end master
    call self%lock%unset ('gpatch_t')
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Increment the cell udpate counter.  This may be overloaded in more complex
!> solvers, which prefer a different way of counting.
!===============================================================================
SUBROUTINE counter_update (self)
  class(gpatch_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%counter_update')
  !$omp atomic
  timer%n_update = timer%n_update + product(self%n)
  call trace%end()
END SUBROUTINE counter_update

!===============================================================================
!> Add default guard zone download functionality.  This may be overloaded in
!> experiment_mod, if more specialized guard zone loading is required.
!===============================================================================
SUBROUTINE dnload (self, only)
  class(gpatch_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%dnload')
!if (self%id==1) print'("dnload",1p,22e9.1)',self%mem(10,1:15,15,5,1,1)-self%mem(10,10,10,5,1,1)
  if (self%is_clear (bits%frozen)) then
    if (self%use_data_hub) then
      call self%data_hub%update (self%link, only=only)
    else
      call download%download_link (self%link, only=only)
    end if
  end if
!if (self%id==1) print'("dnload",1p,22e9.1)',self%mem(10,1:15,15,5,1,1)-self%mem(10,10,10,5,1,1)
  call trace%end()
END SUBROUTINE dnload

!===============================================================================
!> Intercept to make validate%check call
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(gpatch_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('gpatch_t%courant_condition')
  call validate%check (self%link, self%u_max, 'courant')
  call self%patch_t%courant_condition (detailed_timer)
  call trace%end()
END SUBROUTINE courant_condition

END MODULE gpatch_mod
