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
  USE rt_nbors_mod
  USE solver_mod
  USE list_mod
  USE download_mod
  USE patch_mod
  USE link_mod
  USE rt_mod
  USE mesh_mod
  USE task_mod
  USE eos_mod
  USE data_io_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE scaling_mod
  USE bits_mod
  USE omp_timer_mod
  implicit none
  private
  !.............................................................................
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: init_task_list
    procedure:: update
  end type
  !.............................................................................
  integer, save:: verbose=1
  type(experiment_t), public:: experiment
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%mhd = .true.
  call self%solver_t%init
  self%unsigned = .false.
  !----------------------------------------------------------------------------
  ! Boundary setup
  !----------------------------------------------------------------------------
  call self%init_bdries
  call self%boundaries%init (self)
  self%boundaries%id = self%id
END SUBROUTINE init

!===============================================================================
!> Intercept the init_task_list call to finalize RT setup
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(experiment_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  call trace%begin ('experiment_t%init_task_list')
  ! -- punt to the next level; picked up by extras_t ---------------------------
  call self%solver_t%init_task_list (task_list)
  ! -- for RT tasks, modify the nbor lists and add tasks to the task list-------
  if (associated (self%rt)) then
    if (self%rt%on) then
      call rt_nbors%init (self)
      call self%rt%add_to_task_list
    end if
  end if
  call trace%end()
END SUBROUTINE init_task_list

!===============================================================================
!> Experiment update
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  real, dimension(:,:,:), pointer :: mem
  integer:: iv
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  do iv=1, self%nv
    mem => self%mem(:,:,:,iv,self%it,1)
    call self%boundaries%extrapolate(mem, iv)
  end do
  call self%solver_t%update
  call trace%end()
END SUBROUTINE update

END MODULE experiment_mod
