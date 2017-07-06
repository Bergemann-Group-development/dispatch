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
  USE pic_mod
  USE mhd_mod
  USE initial_mod
  USE initial_pic_mod
  USE timestep_mod
  USE timer_mod
  USE link_mod
  USE task_mod
  implicit none
  private
  type, public, extends(pic_t):: experiment_t
      type(initial_pic_t):: initial
  contains
    procedure :: init
    procedure :: update
  end type
  logical, save:: first_run=.true.

CONTAINS

!> =============================================================================
!> Setup initial state
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  real, dimension(:,:,:,:), pointer:: ff
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  !----------------------------------------------------------------------------
  ! Initialize the patch
  !----------------------------------------------------------------------------
  self%nv = 8
  call self%print_patch_type(self%initial%pic_needed)
  call self%initial%init                              ! initialize the experiment
  call self%pic_t%init
  call self%set(bits%static)
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  ff => self%mem(:,:,:,:,1,1)
  call self%initial%init_mhd_vars(self%mesh, ff)
  self%mhd = self%initial%mhd
  self%time = self%initial%time
  self%t(self%it) = self%initial%time
  call trace_end (itimer=itimer)
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
  class(link_t), pointer:: link
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  !print*, "patch id =", self%id, "pic_needed =", self%initial%pic_needed
  if (first_run) then
     first_run = .false.
    if (.not.self%initial%pic_needed) call self%initial%cleanup
  end if
  call self%dnload
  call self%output
  if (self%initial%pic_needed) then
    call self%pic_t%update
  else
    call self%mhd_t%update
    call self%courant_condition
    call self%check_density('before timestep')
    call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
    call self%check_density(' after timestep')
  end if
  call trace_end
END SUBROUTINE update
END MODULE
