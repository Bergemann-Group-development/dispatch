!===============================================================================
!> Template version of a module intended for adding extra features, at a level
!> between the basic patch and the layer of solvers.
!>
!> To add extra features to the code, copy this file to experiments/your_dir/,
!> and add (only) the extra dependencies to the Makefile there.
!>
!>                    experiment_t
!>                     solver_t
!>                       mhd_t
!>                      extras_t
!>                      /    \
!>                     /      \
!>                extras_1_t extras_2_t
!>                       \    /
!>                      gpatch_t
!>                       patch_t
!>                       /
!>                    task_t
!>
!> The extra modules needed in an experiment should be chosen and coordinated
!> in this module, and the ones that are "below" (USEd by) this module should
!> generally refer to array memory in patch_t, and/or to local memory allocated
!> inside these extra modules.  Only key physics effects, such as forcing,
!> cooling or heating, ionization, etc., should be added to the patch_t module,
!>
!> Modules that prepare input to the solvers should generally be called from the
!> pre_upate() procedure, while modules that rely on results from the solver
!> should be called from the post_update() procedure.
!>
!> This module should generally NOT use the h5_mod (HDF5) module itself. It
!> should instead be USEd from the actual extra modules, to assist in saving
!> ad hoc results in a standardized format.
!===============================================================================

MODULE extras_mod
  !-----------------------------------------------------------------------------
  ! Generic modules, which are always included in the code in any case
  !-----------------------------------------------------------------------------
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE patch_mod
  USE task_mod
  USE list_mod
  USE link_mod
  USE bits_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.
  !-----------------------------------------------------------------------------
!  USE boundaries_mod
  USE selfgravity_mod                                           ! selfg
  USE particle_solver_mod                                       ! dust
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    class(gpatch_t), pointer:: gpatch => null()                 ! keep
    class(selfgravity_t), pointer:: selfgravity                 ! selfg
    class(particle_solver_t), pointer:: particle_solver         ! dust
    real(8), dimension(:,:,:,:),   contiguous, pointer:: r       => null()
    real(8), dimension(:,:,:,:,:), contiguous, pointer:: r_hat   => null()
    real(8), dimension(:,:,:,:),   contiguous, pointer:: phi_hat => null()
    real(8), dimension(:,:,:,:),   contiguous, pointer:: v_patch => null()
    real(8), dimension(:,:,:,:,:), contiguous, pointer:: f_grav  => null()
    real(8), dimension(:,:,:,:,:), contiguous, pointer:: r_inv   => null()
    real(8), dimension(:,:,:),     contiguous, pointer:: v_head  => null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_task_list
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
    procedure:: check_refine
    procedure:: courant_condition
    procedure:: statistics
    procedure:: is_ahead_of
  end type
  type(extras_t), public:: extras
CONTAINS

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE pre_init (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_init')
  self%gpatch => self                                             ! keep
  call self%gpatch_t%pre_init ()                                  ! keep
  allocate (self%selfgravity)                                     ! selfg
  call self%selfgravity%pre_init (self%gpatch)                    ! selfg
  self%particle_solver%patch => self                              ! dust
  call self%particle_solver%pre_init ()                           ! dust
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  call self%gpatch_t%init                                         ! keep
  call self%selfgravity%init                                      ! selfg
  call self%particle_solver%init                                  ! dust
  call trace%end()
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize pointers to the task list
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(extras_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  call trace%begin ('extras_t%init_task_list')
  call self%gpatch_t%init_task_list (task_list)
  call self%particle_solver%init_task_list (task_list)
  call trace%end()
END SUBROUTINE init_task_list

!===============================================================================
!> The pre_update hierarchy is called before the solver updates
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t), target   :: self
  class(gpatch_t), pointer  :: patch
  real, allocatable         :: dd(:,:,:,:)  ! dust density
  real, pointer, contiguous :: dg(:,:,:)    ! gas density
  real, pointer             :: pg(:,:,:,:)  ! gas momentum
  integer                   :: j
  integer, save             :: itimer=0
  !.............................................................................
  call trace%begin ('extras_t%pre_update', itimer=itimer)
  self%gpatch => self                                              ! keep
  call self%gpatch_t%pre_update ()                                 ! keep
  if (self%particle_solver%on) then
    if (self%selfgravity%on) then
      !---------------------------------------------------------------------------
      ! Add dust to gas density before computation of self-gravity
      !---------------------------------------------------------------------------
      patch => self%particle_solver%patch
      dg => patch%mem(:,:,:,patch%idx%d,patch%it,1)
      allocate(dd(size(dg,1),size(dg,2),size(dg,3),self%particle_solver%n_bin))
      call self%particle_solver%density (dd)
      dg(:,:,:) = dg(:,:,:) + sum (dd(:,:,:,:),dim=4)
      call self%selfgravity%pre_update                             ! selfg
      !---------------------------------------------------------------------------
      ! Subtract dust from gas density after computation of self-gravity
      !---------------------------------------------------------------------------
      dg(:,:,:) = dg(:,:,:) - sum (dd(:,:,:,:),dim=4)
      deallocate(dd)
      !FIXME: self%particle_solver%grad_phi => self%selfgravity%grad_phi
    end if
    if (self%time + self%dtime >= self%out_next) &
      call self%particle_solver%exports_info
    self%particle_solver%v_patch => self%v_patch
    call self%particle_solver%update                               ! dust
  end if
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  class(gpatch_t), pointer  :: patch
  real, allocatable         :: dd(:,:,:,:)  ! dust density
  real, pointer, contiguous :: dg(:,:,:)    ! gas density
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('extras_t%post_update', itimer=itimer)
  call self%gpatch_t%post_update ()                                ! keep
  !-----------------------------------------------------------------------------
  ! Add dust to gas density before computation of self-gravity
  !-----------------------------------------------------------------------------
  if (self%selfgravity%on) then
    patch => self%particle_solver%patch
    dg => patch%mem(:,:,:,patch%idx%d,patch%new,1)
    allocate(dd(size(dg,1),size(dg,2),size(dg,3),self%particle_solver%n_bin))
    call self%particle_solver%density (dd)
    dg(:,:,:) = dg(:,:,:) + sum (dd(:,:,:,:),dim=4)
    call self%selfgravity%post_update                                ! selfg
    !---------------------------------------------------------------------------
    ! Subtract dust from gas density after computation of self-gravity
    !---------------------------------------------------------------------------
    dg(:,:,:) = dg(:,:,:) - sum (dd(:,:,:,:),dim=4)
    deallocate(dd)
  end if
  !-----------------------------------------------------------------------------
  ! Generate new nbor lists when patches have moved enough (FIXME)
  !-----------------------------------------------------------------------------
  call self%set (bits%init_nbors)
  call trace%end (itimer)
END SUBROUTINE post_update

!===============================================================================
!> Optionally call specific output modules for the extra features
!===============================================================================
SUBROUTINE output (self)
  class(extras_t):: self
  !.............................................................................
  call self%gpatch_t%output                                        ! keep
END SUBROUTINE output

!===============================================================================
!> Interface to refinement.
!===============================================================================
INTEGER function check_refine (self)
  class(extras_t), target:: self
  !.............................................................................
  check_refine = -1
  !FIXME: check_refine = self%particle_solver%refine_Flash ()
END FUNCTION check_refine

!===============================================================================
!> Intercept courant_condition() calls; make sure the gpatch_t call is the last
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(extras_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('extras_t%courant_condition')
  call self%gpatch_t%courant_condition ()                          ! keep
  call self%particle_solver%courant_condition                      ! dust
  call trace%end()
END SUBROUTINE courant_condition

!===============================================================================
!> Called to allow statistics collection
!===============================================================================
SUBROUTINE statistics (self)
  class(extras_t):: self
  call self%particle_solver%statistics
END SUBROUTINE statistics

!===============================================================================
!> Modify grace, to force nbor tasks to be very close in time, before output.
!> This is called from task_list%check_ready(), with self as the nbor and
!> target as the task to check.
!===============================================================================
FUNCTION is_ahead_of (self, target)
  class(extras_t):: self
  class(task_t):: target
  logical:: is_ahead_of
  real:: grace
  !-----------------------------------------------------------------------------
  grace=self%grace
  if (target%out_next - target%time < 1.1*target%dtime) then
    grace = self%grace
    self%grace=0.0001
    is_ahead_of = self%task_t%is_ahead_of (target)
    self%grace = grace
  else
    is_ahead_of = self%task_t%is_ahead_of (target)
  end if
END FUNCTION is_ahead_of

END MODULE extras_mod
