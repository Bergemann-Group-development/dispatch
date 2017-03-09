!> =============================================================================
!> This module contains all experiment specific information.  The experiment_t
!> object extends the generic MHD object mhd_t, which is again an extension of
!> the basic patch_t object for tasks (object task_t) with meshes.
!>
!> The MHD object is chosen with the SOLVER macro in the Makefile in the curren
!> directory, while the generic task_t and patch_t objects are defined in the
!> $(TOP)/tasks/ directory, in task_mod.f90 and patch_mod.f90.
!>
!> The init procedure here (experiment_t::init) is called from the cartesian_mod
!> module, which is specified as the patch distributor in the main program
!> dispatch.f90.   The update procedure here is called by the task_list_mod
!> update procedure in $(TOP)/lists/task_list_mod.f90.
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE rmhd_mod
  USE initial_mod
  USE boundary_mod
  implicit none
  private
  type, public, extends(rmhd_t):: experiment_t
    type(initial_t):: initial
    type(boundary_t) :: boundary
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%nv = 8
  self%mhd = .false.
  call self%rmhd_t%init()
  call self%set(bits%static)
  !self%BC_direction = 1
  if (self%position(1) - self%size(1) < 0.0) then
    call self%set(bits%has_BC)
  else if (self%position(1) + self%size(1) > self%box(1)) then
    call self%set(bits%has_BC)
  end if
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  call self%initial%init()
  call self%boundary%init()
  call self%initial%initial_conditions(self)
  call self%stats ('initial')
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).  The next step
!> is to call a timestep%udpate, whose actual implementation is in a file
!> $(TOP)/timestep/$(UPD)/timestep_mod.f90 -- currently an Adams-Bashforth
!> method.
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%rmhd_t%update
   if (io%verbose>2) &
    call self%check_density('before timestep')
  call self%check_density(' after timestep')

  call trace_end
END SUBROUTINE update

END MODULE
