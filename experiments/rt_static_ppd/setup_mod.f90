!> =============================================================================
!> The duty of the `experiment_mod` module is to define a generic (or specific)
!> `experiment_t` object, which has an initialisation procedure and a task list.
!> The tasks can be of any type belonging to the task class.
!> =============================================================================
MODULE setup_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE units_mod
  USE task_mod
  USE task_list_mod
  USE mesh_mod
  USE experiment_mod
  USE scaling_mod
  implicit none
  private

  ! components of the simulation
  ! an `setup_t` object is required to tie everything together
  type, public:: setup_t
    type(task_list_t), pointer :: task_list
    character(len=64) :: experiment_name = 'static RT PPD benchmark'
  contains
    procedure:: init => init_setup
  end type
  type(setup_t), public :: setup

CONTAINS

!> =============================================================================
!> The init procedure just splices together the relevant task lists.  The
!> specific behavior of each component is described in its respective module.
!> =============================================================================
SUBROUTINE init_setup (self)
  class(setup_t):: self
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  logical, save:: first_time = .true.
  real, save:: r_min=1.0, r_max=1.0e3, z_min=-5.0e2, z_max=5.0e2
  integer, save:: nr=128, nz=180, npatches(3)=1
  integer:: ip1, ip2, ip3, n_per_patch(3)
  real:: zextent, rextent, phiextent
  namelist / domain_params / r_min, r_max, z_min, z_max, nr, nz, npatches
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')

  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, domain_params)
    if (io%master) write (*, domain_params)
    call scaling%init()
  end if

  ! Create an empty task list.
  allocate(self%task_list)

  zextent = (z_max - z_min) * cgs%AU / npatches(1)
  rextent = (r_max - r_min) * cgs%AU / npatches(2)
  phiextent = scaling%pi * 2.0  ! npatches(3) is ignored!
  n_per_patch = [nz / npatches(1), nr / npatches(2), 1]

  do ip3=1,npatches(3)
    do ip2=1,npatches(2)
      do ip1=1,npatches(1)
        allocate (experiment_t:: patch)
        patch%size = [zextent, rextent, phiextent]
        patch%position = patch%size * 0.5 ! FIXME: This is a hack because the STAGGER solver does not currently support cylindrical coords.
        patch%n = n_per_patch
        !patch%origin = 0.0  ! FIXME: introduced in curvilinear_coords branch
        !patch%llc = 0.0  ! FIXME: introduced in curvilinear_coords branch
        !patch%staggered = .true.  ! FIXME: introduced in curvilinear_coords branch
        !patch%mesh_type = mesh_types%Cartesian  ! FIXME: introduced in curvilinear_coords branch
        patch%box = [(z_max - z_min), (r_max - r_min), phiextent]

        call patch%init()
        task => patch
        call self%task_list%append(task)
        if (io%verbose > 0) &
          write (*, '(a,i6,2(1x,3e12.4))') 'init_setup: id, size, pos:', &
            patch%id, patch%size/cgs%AU, patch%position/cgs%AU
      end do
    end do
  end do

  ! Done. Perform final list initialisations.
  call self%task_list%init(trim(self%experiment_name)//' task list')

  call trace_end
END SUBROUTINE init_setup

END MODULE setup_mod
