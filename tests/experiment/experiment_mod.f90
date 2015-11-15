!> =============================================================================
!> $Id$
!> Skeleton experiment, which does nothing, but has the necessary procedures
!> =============================================================================
MODULE experiment_mod
  USE list_mod
  USE patch_mod
  USE pde_mod
  USE mhd_mod
  USE mpi_coord_mod
  USE io_mod
  USE task_list_mod
  type, extends(pde_t):: experiment_t
    real:: u(3), w, s
    integer:: n
  contains
    procedure:: init
  end type
  real:: velocity(3)   = [1.0, 0.5, 0.0]        ! velocity components
  real:: width         = 1.0                    ! Gaussian width
  real:: size          = 16.                    ! experiment size
  integer:: resolution = 64                     ! total number of points
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
  integer:: npatch, dims(3)
  namelist /parameters/ size, width, velocity, resolution
  !.............................................................................
  read (input, parameters); if (io%master) write (io%output, parameters)
  template%box%size = size
  template%box%dim  = resolution
  dims = resolution
  call mpi_coord%init (dims=dims)
  do i=1,ntask
    template%itask = i
    call template%init                          ! allocates and initializes
    ! now we can set task parameters from the template
    if (template%rank == mpi_rank) &            ! if on this rank
      allocate (task)
      item => task
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
  !.............................................................................
  rho => self%mem(:,:,:,1,self%it,1)
  !rho = gaussian(self)
  do i3=1,self%m(3)%gn
  do i2=1,self%m(2)%gn
  do i1=1,self%m(1)%gn
    x1 = self%m(1)%p + self%m(1)%r(i1)
    x2 = self%m(2)%p + self%m(2)%r(i2)
    x3 = self%m(3)%p + self%m(3)%r(i3)
    rho(i1,i2,i3) = exp(-(x**2+y**2+z**2)/self%w**2)
  end do
  end do
  end do
END SUBROUTINE

!> =============================================================================
!> The Gaussian gets the local task position via the task argument, and thus
!> returns the proper piece of the IC.
!> =============================================================================
FUNCTION gaussion (task) result (f)
  class(experiment_t), pointer:: task
  real, dimension(task%m(1)%gn, task%m(2)%gn, task%m(3)%gn):: f
  !.............................................................................
  do i3=1,task%m(3)%gn
  do i2=1,task%m(2)%gn
  do i1=1,task%m(1)%gn
    x1 = task%m(1)%p + task%m(1)%r(i1)
    x2 = task%m(2)%p + task%m(2)%r(i2)
    x3 = task%m(3)%p + task%m(3)%r(i3)
    f(i1,i2,i3) = exp(-(x**2+y**2+z**2)/task%w**2)
  end do
  end do
  end do
END FUNCTION

END MODULE experiment_mod
