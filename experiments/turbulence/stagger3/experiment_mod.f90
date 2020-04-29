!===============================================================================
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
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE mhd_mod
  USE initial_mod
  USE timer_mod
  USE link_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  real, dimension(:,:,:,:), pointer:: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  call self%initial%init (self%kind)
  if (self%initial%mhd) then
    self%nv = 8
  else
    self%nv = 5
  end if
  self%mhd = self%initial%mhd
  call self%mhd_t%init
  call self%set (bits%static)
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  mem => self%mem(:,:,:,:,1,1)
  call self%initial%conditions (self%mesh, mem, self%idx)
  self%time       = self%initial%time
  self%t(self%it) = self%initial%time
  call self%periodic_grid
  call trace_end (itimer=itimer)
END SUBROUTINE init

!===============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    call self%initial%cleanup
  end if
  end if
  call self%mhd_t%update
  call trace_end
END SUBROUTINE update

END MODULE
