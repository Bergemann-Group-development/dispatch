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
  USE kinds_mod
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
  USe vector_ops
  implicit none
  private
  !.............................................................................
  type, public, extends(solver_t):: experiment_t
    real:: d_bot  = 0.0
    real:: e_bot  = 0.0
    real:: ee_top = 0.0
  contains
    procedure:: init
    procedure:: update
    procedure:: bndries
    procedure:: clamp_ee
    procedure, nopass:: cast2experiment
  end type
  !.............................................................................
  integer, save:: verbose=1
  logical, save:: clamp_e = .false.
  real, save:: d_bot=0.0, e_bot=0.0
  real, save:: bottom_bc_rate=0.01
  integer, save:: divb_iter = 4
  type(experiment_t), public:: experiment
CONTAINS

!===============================================================================
!> Experiment setup
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: iostat
  real, pointer:: ff(:,:,:,:)
  logical, save:: first_time=.true.
  namelist /experiment_params/ clamp_e
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')
  !-----------------------------------------------------------------------------
  call io%lock%set ('experiment_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    write (io%namelist, experiment_params)
    !d_bot = d_bot/scaling%d
    !e_bot = e_bot/scaling%e
    first_time = .false.
  end if
  call io%lock%unset ('experiment_t%init')
  if (.not. io%master) verbose=0
  !-----------------------------------------------------------------------------
  ! Initialize the MHD task and the RT sub-tasks, and setup their nbor lists
  !-----------------------------------------------------------------------------
  self%mhd = .true.
  call self%solver_t%init
  ! ----------------------------------------------------------------------------
  ! Send additional things to boundaries
  ! ----------------------------------------------------------------------------
  self%boundaries%gravity = self%gravity%const_grav
  self%boundaries%mhd = self%mhd
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Experiment update
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(mesh_t), pointer:: m(:)
  real, pointer:: e(:,:,:), d(:,:,:)
  real, pointer:: e1(:,:,:), d1(:,:,:)
  real:: ee1, ee
  integer:: mlee(3)
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  m => self%mesh
  call self%solver_t%update

  if (clamp_e) call self%clamp_ee

  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:), pointer, contiguous:: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%bndries')
  mem => self%mem(:,:,:,:,self%it,1)
  call self%boundaries%primary (mem)
  call trace%end()
END SUBROUTINE bndries

!===============================================================================
!> Clamp e/d to a minimum table value
!===============================================================================
SUBROUTINE clamp_ee(self)
  class(experiment_t):: self
  real, pointer:: e(:,:,:), d(:,:,:)
  real::ee
  integer:: i,j,k, l(3), u(3)
  ! ----------------------------------------------------------------------------
  l = self%mesh%lb
  u = self%mesh%ub
  e => self%mem(:,:,:,self%idx%e,self%it,1)
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  !
  do k=l(3),u(3)
    do j=l(2),u(2)
      do i=l(1),u(1)
        ee = e(i,j,k)/d(i,j,k)
        e(i,j,k)= max(ee,eos%ee_min)*d(i,j,k)
      end do
    end do
  end do
END subroutine clamp_ee

!===============================================================================
!> Cast a generic task_t to patch_t
!===============================================================================
FUNCTION cast2experiment (task) RESULT(experiment)
  class(task_t), pointer:: task
  class(experiment_t), pointer:: experiment
  !.............................................................................
  select type (task)
  class is (experiment_t)
  experiment => task
  class default
  nullify(experiment)
  call io%abort ('experiment_t%cast2experiment: failed to cast a task to experiment_t')
  end select
END FUNCTION cast2experiment

END MODULE experiment_mod
