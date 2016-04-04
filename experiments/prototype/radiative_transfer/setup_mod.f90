MODULE setup_mod
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a planet at a given position, and a set of Rubik's Cube patches
!> surrounding it, representing its extended atmosphere.  Then add the 3*nray
!> radiative transfer tasks
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
  call io%init (label)                          ! make sure I/O is initialized
  call scaling%init                             ! initialize scaling to code units
  call planet%init                              ! initialize planet
  call atmosphere%init                          ! initialize its atmosphere
  self%size = planet%radius*2.2                 ! size of central patch
  planet%position = 0.0                         ! choose as the coordinate center
  self%position = planet%position               ! also the center of the setup
  call self%rubiks_t%init_hierarchy (label)     ! create a Rubiks cube hierarchy
  link => self%task_list%head                   ! loop variable
  do while (associated(link))                   ! find tasks containing planet
    patch => task2patch (link%task)             ! patch form of task
    if (all(abs(patch%position) < patch%ds)) &
      call patch%set (bits%has_BC)              ! mark it has having BCs
    link => link%next                           ! next link
  end do
  !----------------------------------------------------------------------------
  ! Loop over the existing MHD patches and construct RT tasks for each
  !----------------------------------------------------------------------------
  link => self%task_list%head                   ! loop variable
  do while (associated(link))                   ! task list loop
    do iray=1,self%nray                         ! loop over ray directions
      call self%rt_setup (link, iray)           ! setup for that direction
    end do
    link => link%next                           ! next task list link
  end do
  call trace_end
END SUBROUTINE init

SUBROUTINE rt_setup (self, iray)
  class(setup_t):: self
  integer:: iray
  class(link_t), pointer:: link, nbor
  type(rt_t), pointer:: task
  !-----------------------------------------------------------------------------
  ! Construct the RT task for the internal solution and add it to the task list,
  ! giving it a nbor list consisting of only the MHD task
  !-----------------------------------------------------------------------------
  allocate (task)                               ! the internal RT task
  call task%init_internal (iray)                ! initialize in ray direction
  allocate (link, nbor)                         ! the task link and a neighbor
  link%task => task                             ! connect to RT task
  nbor%task => self%task                        ! connect to MHD task
  link%nbor => nbor                             ! make MHD the (only) nbor
  call self%task_list%insert (link)             ! insert into the task list
  !> ---------------------------------------------------------------------------
  !> Construct the RT task for the upstream correction and add to the task list
  !> ---------------------------------------------------------------------------
  allocate (link)                               ! a new link
  allocate (task)                               ! its new task
  link%task => task                             ! connect
  call task%init_upstream (iray)                ! initialize as upstream
  allocate (nbor)                               ! new nbor link
  link%nbor => nbor                             ! is the 1st neighbor
  nbor%task => self%task                        ! hosting the MHD task
  nbor%next => self%nbor                        ! continuing with normal nbor
  nbor => self%nbor                             ! loop pointer
  prev => nbor                                  ! remember previous
  do while associated(nbor)                     ! loop over nbors
    if (nbor%is_upstream (self)) then           ! forward ray
      allocate (new)                            ! new nbor link
      new%task => nbor%task                     ! task is from the nbor list
      prev%next => new                          ! linked into the list
    end if
    nbor => nbor%next                           ! try the next nbor
  end do
  call self%task_list%insert (link)             ! insert into the task list
  !> ---------------------------------------------------------------------------
  !> Construct the RT task for the downstream correction and add to the task list
  !> ---------------------------------------------------------------------------
  allocate (link)                               ! a new link
  allocate (task)                               ! its new task
  link%task => task                             ! connect
  call task%init_downstream (iray)              ! initialize as downstream
  allocate (nbor)                               ! new nbor link
  link%nbor => nbor                             ! is the 1st neighbor
  nbor%task => self%task                        ! hosting the MHD task
  nbor%next => self%nbor                        ! continuing with normal nbor
  nbor => self%nbor                             ! loop pointer
  prev => nbor                                  ! remember previous
  do while associated(nbor)                     ! loop over nbors
    if (nbor%is_downstream (self)) then         ! backward ray
      allocate (new)                            ! new nbor link
      new%task => nbor%task                     ! task is from the nbor list
      prev%next => new                          ! linked into the list
    end if
    nbor => nbor%next                           ! try the next nbor
  end do
  call self%task_list%insert (link)             ! insert into the task list
END SUBROUTINE rt_setup

SUBROUTINE close
END SUBROUTINE close

END MODULE setup_mod
