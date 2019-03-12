!===============================================================================
!> Setup a simple envelope component, consisting of a Cartesian arrangement of
!> MHD patches.  The module is an extension of the generic component_t, since it
!> contains the extra parameters density and n_patch
!===============================================================================
MODULE envelope_mod
  USE task_list_mod
  USE list_mod
  USE task_mod
  USE mhd_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: envelope_t
    real:: density=1.0
    real:: size=1.0
    integer:: n_patch(3)=2
    type(task_list_t), pointer:: task_list
  contains
    procedure:: init
  end type
  type(envelope_t), pointer:: envelope
contains
!===============================================================================
!> The init procedure reads parameters, allocates and adds patches to the task
!> list
!===============================================================================
SUBROUTINE init (self, name)
  class(envelope_t):: self
  character(len=*):: name
  class(task_t), pointer:: patch
  integer:: ix, iy, iz
  real:: size=1.0, density=1.0
  integer:: n_patch(3)=4
  namelist /envelope_params/ n_patch, size, density
  !-----------------------------------------------------------------------------
  ! Envelop section, with a large number of MHD patches
  !-----------------------------------------------------------------------------
  call trace_begin ('envelope_t%init')
  rewind (io%input); read (io%input, envelope_params)
  self%n_patch = n_patch
  self%size = size
  self%density = density
  !-----------------------------------------------------------------------------
  allocate (self%task_list)
  call self%task_list%list_t%init ('envelope task list')
  do iz=1,self%n_patch(3)
  do iy=1,self%n_patch(2)
  do ix=1,self%n_patch(1)
    allocate (patch)
    call patch%init
    patch%size = self%size/self%n_patch
    patch%position = ([ix,iy,iz]-0.5)*patch%size
    patch%name => self%task_list%name
    call trace('append')
    call self%task_list%append (patch)
  end do
  end do
  end do
  call trace_end
END SUBROUTINE
END MODULE envelope_mod
