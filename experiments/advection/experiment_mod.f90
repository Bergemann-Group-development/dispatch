!> =============================================================================
!> $Id$
!> Advection experiment, runs runs the advection problem in separate patches.
!> =============================================================================
MODULE experiment_mod
  USE pde_mod                   ! patches/centered/advection_mod.f90
  USE patch_mod                 ! tasks/patch_mod.f90
  USE task_list_mod             ! lists/task_list_mod.f90
  USE io_mod                    ! io/io_mod.f90
  USE mpi_coord_mod             ! use cartesian coordinates
  USE nbor_list_mod             !
  USE bits_mod
  USE mpi_mod
  USE item_mod
  USE list_mod
  USE unigrid_mod
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
  real:: end_time      = 1.0
  real:: velocity(3)   = [1.0, 0.5, 0.0]        ! velocity components
  real:: width         = 1.0                    ! Gaussian width
  real:: size          = 16.                    ! experiment size
  integer:: resolution = 64                     ! total number of points
  integer:: patch_resolution = 16               ! number of points per patch
CONTAINS

!> =============================================================================
!> Populate a task_list (complete list of all tasks relevant to this rank --
!> local tasks and virtual task = copies of tasks from neighbor ranks).
!> =============================================================================
SUBROUTINE init_experiment (task_list)
  implicit none
  type(task_list_t):: task_list
  type(box_t):: box, task_box, rank_box
  class(item_t), pointer:: item
  type(link_t), pointer:: link
  integer:: resolution, dims(3), patch_dims(3), rank_c(3), task_c(3)
  real:: size, width, velocity(3)
  type(experiment_t), pointer:: task
  integer:: mpi_x, mpi_y, mpi_z, patch_x, patch_y, patch_z, rank
  namelist /parameters/ end_time, size, width, velocity, resolution
  !.............................................................................
  rewind (io%input); read (io%input, parameters)  ! experiment parameters
  if (io%master) write (*, parameters)
  dims = resolution                               ! unigrid 3-D resolution
  patch_dims = patch_resolution                   ! unigrid patch resolution
  call unigrid%init                               ! get MPI geometry
  unigrid%size = size                             ! set size
  unigrid%resolution = resolution                 ! set resolution

  ! FIXME call unigrid%generate (task_list)               ! make a task list
  !-----------------------------------------------------------------------------
  ! Generate initial conditions, using the task init routine, which relies on
  ! using the task meshes, which were initialized when the task was allocated
  !-----------------------------------------------------------------------------
  link => task_list%head
  do rank=0,mpi%size-1
    rank_c = mpi_coord%rank_to_coords(rank)
    do patch_z=0,unigrid%task_box%num(3)-1
    do patch_y=0,unigrid%task_box%num(2)-1
    do patch_x=0,unigrid%task_box%num(1)-1
      box      = unigrid%box
      task_box = unigrid%task_box
      rank_box = unigrid%rank_box
      allocate (task)                             ! allocate new task
      task%rank = rank                            ! set its rank
      task%size = task_box%size                   ! set task size
      task_c = [patch_x,patch_y,patch_z]          ! integer patch position
      task%position = &                           ! set task position
        (rank_c*task_box%num + task_c + 0.5d0)*task_box%size - 0.5d0*box%size
      item => task                                ! generic task pointer
      call task_list%append (item)                ! append to task_list
      if (task%rank==mpi%rank) then
        call task%set (bits%internal)              ! owned by this rank
      else
        call task%set (bits%external)              ! owned by another rank
      end if
    end do
    end do
    end do
  end do
!  call task_list%init_nbors
  link => task_list%head
  do while (associated(link))
    item => link%item
    select type (item)
    class is (patch_t)
      if (item%is_set(bits%internal+bits%virtual)) then
        call item%init
      end if
    end select
    link => link%next
  end do
END SUBROUTINE

!> =============================================================================
!> Set the initial conditions for a single task.  This must reside in the same
!> module that reads in experiment parameters, so the place is here.  It could
!> in principle also contain a loop over tasks, but that is a more generic thing.
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real, dimension(:,:,:), pointer:: rho
  integer:: i1, i2, i3
  real:: x1, x2, x3
  !.............................................................................
  self%end_time = end_time
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
