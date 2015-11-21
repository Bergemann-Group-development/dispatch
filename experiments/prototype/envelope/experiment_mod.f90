!> =============================================================================
!> The duty of the experiment_mod module is to define a generic experiment_t
!> object, which has an init procedure that populates a task list.
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE task_list_mod
  implicit none
  private
  type, public:: experiment_t
    real:: size=1.0, density=1.0
    integer:: n_patch(3)=4
    type(task_list_t), pointer:: task_list
  contains
    procedure:: init
    procedure:: end
  end type
contains

!> =============================================================================
!> The init procedure reads parameters, allocates and adds patches to the task
!> list
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  class(task_t), pointer:: patch
  integer:: ix, iy, iz
  real:: size=1.0, density=1.0, dtime=0.1, end_time=1
  integer:: n_patch(3)=4
  namelist /envelope_params/ n_patch, size, density, end_time
  !-----------------------------------------------------------------------------
  ! Envelop section, with a number of MHD patches
  !-----------------------------------------------------------------------------
  call trace_begin ('envelop experiment')
  rewind (io%input); read (io%input, envelope_params)
  self%n_patch = n_patch
  self%size = size
  self%density = density
  allocate (self%task_list)
  self%task_list%name = 'envelope task list'
  !-----------------------------------------------------------------------------
  print *,'envelope/experiment_t%init: envelope experiment patches'
  do iz=1,self%n_patch(3)
  do iy=1,self%n_patch(2)
  do ix=1,self%n_patch(1)
    allocate (patch)
    patch%size = self%size/self%n_patch
    patch%dtime = dtime
    patch%end_time = end_time
    patch%position = ([ix,iy,iz]-0.5)*patch%size
    print '(3i4,4x,3f8.3)',ix,iy,iz,patch%position
    call patch%init
    call self%task_list%append (patch)
  end do
  end do
  end do
  call trace_end
END SUBROUTINE

!> =============================================================================
SUBROUTINE end (self)
  class(experiment_t):: self
END SUBROUTINE

END MODULE experiment_mod
