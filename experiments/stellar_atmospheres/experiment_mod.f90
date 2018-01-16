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
  USE timestep_mod
  USE timer_mod
  USE link_mod
  implicit none
  private
  type, public, extends(rmhd_t):: experiment_t
  contains
    procedure init
    procedure update
    procedure boundary_conditions
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  call self%initial%init (self%kind)
  self%nv  = self%initial%nv
  self%mhd = self%initial%mhd
  self%periodic(2) = .false.                     ! non-periodic vertically
  self%initial%id = self%id
  call self%rmhd_t%init
  if (self%mhd) then
    self%idx%bx=6
    self%idx%by=7
    self%idx%bz=8
    self%idx%tt=9
  else
    self%idx%tt=6
  end if
  call self%set(bits%static)
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  self%time = self%initial%time
  self%t(:) = self%initial%time
  self%dtime = 1.0
  self%dt(:) = self%dtime
  call self%boundary_conditions
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
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  if (first_time) then
    first_time = .false.
    call self%initial%cleanup
  end if
  call self%rmhd_t%update
  call self%boundary_conditions
  call self%courant_condition
  if (io%verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
  call trace_end
END SUBROUTINE update

SUBROUTINE boundary_conditions (self)
  class(experiment_t):: self
  integer:: iy, li(3), ui(3)
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%boundary_conditions')
  li = self%mesh%li
  ui = self%mesh%ui
  if (self%mesh(2)%lower_boundary) then
    do iy=self%mesh(2)%lb,self%mesh(2)%lo
      self%mem(:,iy,:,:,self%it,1) = self%mem(:,self%mesh(2)%li,:,:,self%it,1)
    end do
  end if
  if (self%mesh(2)%upper_boundary) then
    do iy=self%mesh(2)%uo,self%mesh(2)%ub
      self%mem(:,iy,:,:,self%it,1) = self%mem(:,self%mesh(2)%ui,:,:,self%it,1)
    end do
  end if
  call trace_end
END SUBROUTINE boundary_conditions

END MODULE
