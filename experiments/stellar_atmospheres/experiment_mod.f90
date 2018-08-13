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
  USE timer_mod
  USE lagrange_mod
  USE mesh_mod
  USE gpatch_mod
  USE link_mod
  USE rt_solver_mod
  USE solver_mod
  implicit none
  private
  type, public, extends(rt_solver_t):: experiment_t
    real, dimension(:,:,:,:), pointer:: lower => null()
    real, dimension(:,:,:,:), pointer:: upper => null()
  contains
    procedure init
    procedure update
    procedure boundary_conditions
  end type
  real, save, allocatable, dimension(:,:,:,:):: buffer
  integer, save:: verbose=1
  integer, save:: dims(3)=63
  integer, save:: order=2
  real, save:: grav=2.75
CONTAINS

!> =============================================================================
!> Experiment setup
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  namelist /experiment_params/ grav, dims, order, verbose
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  !----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (init_cr)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, experiment_params, iostat=iostat)
      if (io%master) write (*, experiment_params)
    end if
    !$omp end critical (init_cr)
  end if
  if (.not. io%master) verbose=0
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  call self%rt_solver_t%init
  call self%set(bits%static)
  self%force_per_unit_mass => force
  !----------------------------------------------------------------------------
  ! Boundary setup
  !----------------------------------------------------------------------------
  self%periodic(3) = .false.                     ! non-periodic vertically
  call self%init_bdries
  if (verbose > 0) then
    if (self%mesh(3)%lower_boundary) print *, self%id, 'at lower B'
    if (self%mesh(3)%upper_boundary) print *, self%id, 'at upper B'
  end if
  call trace_end (itimer=itimer)
END SUBROUTINE init

!> =============================================================================
!> Experiment update
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  !call self%link%print_nbors
  if (first_time) then
    first_time = .false.
    if (allocated(buffer)) &
      deallocate (buffer)
      call io%bits_mem (-storage_size (buffer), product(shape(buffer)),'buffer')
  end if
  call self%test_d (50)
  call self%rt_solver_t%update
  call self%test_d (51)
  call self%boundary_conditions
  call self%test_d (52)
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> For now, uses initially frozen BCs, which will work for at least a while --
!> enough to develop and validate the RT
!> =============================================================================
SUBROUTINE boundary_conditions (self)
  class(experiment_t):: self
  integer, dimension(3):: lb, lo, li, ui, uo, ub
  integer:: iz, m(6)
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%boundary_conditions')
  lb = self%mesh%lb; lo = self%mesh%lo; li = self%mesh%li
  ub = self%mesh%ub; uo = self%mesh%uo; ui = self%mesh%ui
  m = shape(self%mem)
  m(3) = self%ng(3)
  if (self%mesh(3)%lower_boundary) then
    if (.not.associated(self%lower)) then
      allocate (self%lower(m(1),m(2),m(3),m(4)))
      call io%bits_mem (storage_size(self%lower), product(shape(self%lower)),'lower')
      self%lower = self%mem(:,:,lb(3):lo(3),:,self%it,1)
      if (verbose > 0) &
        print *,'LBC',self%id,self%position, &
        minval(self%lower(:,:,:,1)),maxval(self%lower(:,:,:,1))
    else
      self%mem(:,:,lb(3):lo(3),:,self%it,1) = self%lower
    end if
  end if
  if (self%mesh(3)%upper_boundary) then
    if (.not.associated(self%upper)) then
      allocate (self%upper(m(1),m(2),m(3),m(4)))
      call io%bits_mem (storage_size(self%upper), product(shape(self%upper)),'upper')
      self%upper = self%mem(:,:,uo(3):ub(3),:,self%it,1)
      if (verbose > 0) &
        print *,'UBC',self%id,self%position, &
        minval(self%upper(:,:,:,1)),maxval(self%upper(:,:,:,1))
    else
      self%mem(:,:,uo(3):ub(3),:,self%it,1) = self%upper
    end if
  end if
  call trace_end
END SUBROUTINE boundary_conditions

!> =============================================================================
!> Stellar atmoshpere with gravity pointing in the positive y-direction
!> =============================================================================
FUNCTION force (self) RESULT (ff)
  class(gpatch_t)   :: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3),3) :: ff
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%force')
  ff(:,:,:,1) = 0.0
  ff(:,:,:,2) = 0.0
  ff(:,:,:,3) = grav
  call trace_end
END FUNCTION force

END MODULE
