!> =============================================================================
!> $Id$
!> Advection experiment, runs runs the advection problem in separate patches.
!> =============================================================================
MODULE experiment_mod
  USE pde_mod                   ! patches/centered/advection_mod.f90
  USE patch_mod                 ! tasks/patch_mod.f90
  USE task_list_mod             ! lists/task_list_mod.f90
  USE io_mod                    ! io/io_mod.f90
  USE mpi_basic, only: mpi_t
  implicit none
  !-----------------------------------------------------------------------------
  ! Specialized PDE object for advection with constant velocity
  !-----------------------------------------------------------------------------
  type, extends(pde_t):: experiment_t
    real:: u(3), w, s
    integer:: n
  contains
    procedure:: init
  end type
  !-----------------------------------------------------------------------------
  ! Parameters for namelist
  !-----------------------------------------------------------------------------
  real:: velocity(3)   = [1.0, 0.5, 0.0]        ! velocity components
  real:: width         = 1.0                    ! Gaussian width
  real:: size          = 16.                    ! experiment size
  integer:: resolution = 64                     ! total number of points
  integer:: patch_resolution = 16               ! number of points per patch
  type(io_t):: io
  type(mpi_t):: mpi
CONTAINS

!> =============================================================================
!> Populate a task_list (complete list of all tasks relevant to this rank --
!> local tasks and virtual task = copies of tasks from neighbor ranks).
!> =============================================================================
SUBROUTINE init_experiment (task_list)
  type(task_list_t):: task_list
  type(template_t):: template
  type(experiment_t), pointer:: task
  class(item_t), pointer:: item
  integer:: ntask, itask, resolution
  namelist /parameters/ size, width, velocity, resolution, patch_resolution
  !.............................................................................
  call io%init
  call mpi%init
  rewind (io%input)
  read (io%input, parameters); if (io%master) write (io%output, parameters)
  template%box%size = size
  template%box%dim  = resolution
  template%task%dim = patch_resolution
  do itask=1,ntask
    template%itask = itask
    call template%init                          ! computes task position etc
    task%position = template%position           ! copy position
    ! now we can set task parameters from the template
    if (template%rank == mpi%rank) &            ! if on this rank
      call task%allocate                        ! allocate memory and mesh
      item => task                              ! generic task pointer
      call task_list%append (item)              ! append to task_list
  end do
END SUBROUTINE

!> =============================================================================
!> Set the initial conditions of a task.
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real, dimension(:,:,:), pointer:: rho
  integer:: i1, i2, i3
  real:: x1, x2, x3
  !.............................................................................
  rho => self%mem(:,:,:,1,self%it,1)
  do i3=1,self%m(3)%gn
  do i2=1,self%m(2)%gn
  do i1=1,self%m(1)%gn
    x1 = self%m(1)%p + self%m(1)%r(i1)
    x2 = self%m(2)%p + self%m(2)%r(i2)
    x3 = self%m(3)%p + self%m(3)%r(i3)
    rho(i1,i2,i3) = exp(-(x1**2+x2**2+x3**2)/self%w**2)
  end do
  end do
  end do
END SUBROUTINE

END MODULE experiment_mod
