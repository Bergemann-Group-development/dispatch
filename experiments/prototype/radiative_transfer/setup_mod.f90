MODULE setup_mod
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a planet at a given position, and a set of Rubik's Cube patches
!> surrounding it, representing its extended atmosphere.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE bits_mod
  USE omp_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE planet_mod
  USE atmosphere_mod
  USE patch_mod
  USE pebbles_mod
  USE rubiks_mod
  USE radiative_transfer_mod
  USE list_mod
  implicit none
  private
  type, public, extends(rubiks_t):: setup_t
  contains
    procedure:: init
    procedure, nopass:: close
  end type
  type(setup_t), public:: setup
CONTAINS

!> =============================================================================
!> Setup a general Rubik hierarchy: NxNzN patches per level, and zoom factor M
!>
!> Rubik's cube has N=3, M=3
!> Rubik's revenge has N=4, M=2
!>
!> We set the position of the planet as the origin of the coordinate system,
!> which is natural when treating the gravitational influence of the Sun by going
!> into a comoving and rotating coordinate system.  The basic size unit is the
!> size of the central patch, which is taken to be 10% larger than the planet.
!> =============================================================================
SUBROUTINE init (self, label)
  class(setup_t):: self
  character(len=*), optional:: label
  !.............................................................................
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch
  integer:: iray
  !-----------------------------------------------------------------------------
  call trace_begin ('setup_t%init')
  call io%init (label)
  call scaling%init
  call planet%init
  call atmosphere%init
  self%size = planet%radius*2.2
  planet%position = 0.0
  self%position = planet%position
  call self%rubiks_t%init_hierarchy (label)
  link => self%task_list%head
  do while (associated(link))
    patch => task2patch (link%task)
    if (all(abs(patch%position) < patch%ds)) call patch%set (bits%has_BC)
    link => link%next
  end do
  !----------------------------------------------------------------------------
  ! Loop over the existing MHD patches and construct RT tasks for each
  !----------------------------------------------------------------------------
  link => self%task_list%head
  do while (associated(link))
    do iray=1,self%nray
      call self%rt_setup (link, iray)
    end do
    link => link%next
  end do
  call trace_end
END SUBROUTINE init

SUBROUTINE rt_setup (self, iray)
  class(setup_t):: self
  integer:: iray
  type(rt_t), pointer:: task_internal, task_left, task_right
  !-----------------------------------------------------------------------------
  ! Construct the RT task for the internal solution and add it to the task list
  ! giving it a nbor list consisting only of the MHD task
  !-----------------------------------------------------------------------------
  allocate (task_internal)
  call task_internal%init_internal (iray)
  allocate (link, nbor)
  link%task => task_internal
  nbor%task => self%task
  link%nbor => nbor
  call self%task_list%insert (link)
  !-----------------------------------------------------------------------------
  ! Construct the RT task for the left correction and add it to the task list,
  ! starting the task list with the MHD task, and adding selected MHD nbors
  !-----------------------------------------------------------------------------
  allocate (link)
  allocate (task_left)
  link%task => task_left
  call task_left%init_left
  allocate (nbor)
  link%nbor => self
  link%nbor%task => self%task
  link%nbor%next => self%nbor
  nbor => self%nbor
  prev => nbor
  do while associated(nbor)
    if (nbor%is_to_the_left_of(self)) then
      allocate (new)
      prev%next => new
      new%task => nbor%task
    end if
    nbor => nbor%next
  end do
  call self%task_list%insert (link)
  !> ---------------------------------------------------------------------------
  !> Construct the RT task for the right correction and add it to the task list
  !> ---------------------------------------------------------------------------
  allocate (link)
  allocate (task_right)
  link%task => task_right
  call task_left%init_right
  allocate (nbor)
  link%nbor => self
  link%nbor%task => self%task
  link%nbor%next => self%nbor
  nbor => self%nbor
  prev => nbor
  do while associated(nbor)
    if (nbor%is_to_the_right_of(self)) then
      allocate (new)
      prev%next => new
      new%task => nbor%task
    end if
    nbor => nbor%next
  end do
  call self%task_list%insert (link)
END SUBROUTINE rt_setup

SUBROUTINE close
END SUBROUTINE close

END MODULE setup_mod
