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
  USE pic_mod
  USE solver_mod
  USE initial_mod
  USE timestep_mod
  USE timer_mod
  USE link_mod
  USE task_mod
  USE boundary_mod
  USE driver_mod
  USE mesh_mod
  implicit none
  private
  type, public, extends(pic_t):: experiment_t
    type(initial_t):: initial
    logical:: boundary=.false.
    logical:: lower_boundary=.false.
    logical:: upper_boundary=.false.
  contains
    procedure :: init
    procedure :: update
    procedure :: boundary_conditions
  end type
  type, public, extends(task_t):: field
    type(field), pointer:: fields => null()
  end type
  logical, save:: first_run=.true.
  logical, save:: closed_upper=.true.
CONTAINS

!===============================================================================
!> Setup initial state
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  real, dimension(:,:,:,:), pointer:: ff
  logical, save:: first_time=.true.
  namelist /experiment_params/ closed_upper
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
  end if
  call io%lock%unset
  !----------------------------------------------------------------------------
  ! Initialize the patch
  !----------------------------------------------------------------------------  
  self%lower_boundary = (self%position(1) < self%size(1))
  self%upper_boundary = (self%position(1) + self%size(1) > self%box(1))
  self%boundary = self%lower_boundary .or. self%lower_boundary
  call self%initial%init                              ! initialize the experiment 
  if (self%lower_boundary) then
    self%nv = 9
  else
    self%nv = 8
  end if
  !call io%print_patch_type(self%initial%pic_needed)
  call self%solver_t%init
  if (self%lower_boundary) then
    print *,'extra slot for patch id', self%id, self%nv, size(self%mem,4)
  end if
  call self%set(bits%static)
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  ff => self%mem(:,:,:,:,1,1) 
  call self%initial%init_mhd_vars (self%mesh, ff)
  self%mhd = self%initial%mhd
  self%time = self%initial%time
  self%t(self%it) = self%initial%time
  call trace_end (itimer=itimer)
END SUBROUTINE init


!===============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).  The next step
!> is to call a timestep%udpate, whose actual implementation is in a file
!> $(TOP)/timestep/$(UPD)/timestep_mod.f90 -- currently an Adams-Bashforth
!> method.
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  
  !print*, "patch id =", self%id, "pic_needed =", self%initial%pic_needed
  if (first_run) then
    !$omp critical (first_cr)
    if (first_run) then
      first_run = .false.
      if (.not.self%initial%pic_needed) call self%initial%cleanup
    end if
    !$omp end critical (first_cr)
  end if
  if (self%initial%pic_needed) then
    call self%pic_t%update
  else
    !io%do_trace = .true.
    call self%boundary_conditions
    call self%solver_t%update
  end if
  call trace_end
END SUBROUTINE update

!===============================================================================
!===============================================================================
SUBROUTINE boundary_conditions (self)
  class(experiment_t):: self
  real, dimension(:,:,:), pointer:: d, e, px, ee, bx, by, bz, ux, uy, uz
  class(mesh_t), pointer, dimension(:):: m
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%boundary_conditions')
  associate (lb=>self%mesh(1)%lb, lo=>self%mesh(1)%lo, &
             uo=>self%mesh(1)%uo, ub=>self%mesh(1)%ub)
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  bx => self%mem(:,:,:,self%idx%bx,self%it,1)
  by => self%mem(:,:,:,self%idx%by,self%it,1)
  bz => self%mem(:,:,:,self%idx%bz,self%it,1)
  m  => self%mesh
  if (self%lower_boundary) then
    ! -- impose lower BCs here ---
    call log_extrapolate_x_lb_center (d, m)
    call extrapolate_x_lb_face (by, m)
    call extrapolate_x_lb_face (bz, m)
    call self%p2U
    call cylinder_driver (self,m)
    call self%U2p
    call antisymmetric_x_lb_face (px, m)
    allocate (ee(self%gn(1),self%gn(2),self%gn(3)))
    ee = e/d
    call symmetric_x_lb_center (ee, m)
    e(1:lo-1,:,:) = ee(1:lo-1,:,:)*d(1:lo-1,:,:)
    call divB_correction_x_lb (bx,by,bz,m)
    deallocate (ee)
  end if
  if (self%upper_boundary) then
    ! -- impose upper BCs here ---
    call log_extrapolate_x_ub_center (d,m)
    call extrapolate_x_ub_face (by, m)
    call extrapolate_x_ub_face (bz, m)
    if (closed_upper) then
      !-------------------------------------------------------------------------
      ! A closed upper boundary has Ux=0.0 => antisymmetric px
      !-------------------------------------------------------------------------
      call antisymmetric_x_ub_face (px, m)
    else
      !-------------------------------------------------------------------------
      ! An open upper boundary might be possible with vanishing U-derivatives
      !-------------------------------------------------------------------------
      call self%p2U
      ux => self%U(:,:,:,1)
      uy => self%U(:,:,:,2)
      uz => self%U(:,:,:,3)
      call symmetric_x_ub_face   (ux, m)
      call symmetric_x_ub_center (uy, m)
      call symmetric_x_ub_center (uz, m)
      call self%U2p
    end if
    allocate (ee(self%gn(1),self%gn(2),self%gn(3)))
    ee = e/d
    call symmetric_x_ub_center (ee, m)
    e(uo+1:ub,:,:) = ee(uo+1:ub,:,:)*d(uo+1:ub,:,:)
    call divB_correction_x_ub (bx,by,bz,m)
    deallocate (ee)
  end if
  end associate
  call trace_end
END SUBROUTINE boundary_conditions

END MODULE
