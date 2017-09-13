!> =============================================================================
!> Particle patches are collections (lists) of particles histories, which may
!> be used for example to compute the forces at other particles at a particular
!> time.  They may also be evolved, by performing kick-drift-kick on each of its
!> members.
!>
!> Each particle patch has a pointer to list of nbor, whose particles are also
!> treated explicitly, and also a pointer to a list with all the other particle
!> patches, which are added as a collective force
!> =============================================================================
MODULE particle_patch_mod
  USE particle_mod
  USE particle_history_mod
  implicit none
  private
  type, public:: particle_patch_t
    integer:: n
    real:: moves
    real:: mass
    real(8):: position(3)
    type (particle_history_t), pointer :: head => null()
    type (particle_patch_t), pointer:: nbor  => null()
    type (particle_patch_t), pointer:: other => null()
  contains
    procedure :: init
    procedure :: force_on
    procedure :: add_patch_force_to
    procedure :: move_past
    procedure :: randomize
    procedure :: print
  end type
CONTAINS

!> =============================================================================
!> Initialize a list of particle histories
!> =============================================================================
SUBROUTINE init (self, n)
  class (particle_patch_t)                :: self
  integer                            :: n
  !.............................................................................
  type (particle_history_t), pointer :: history, prev
  integer                            :: i
  !-----------------------------------------------------------------------------
  allocate (history)
  call history%init
  self%head => history
  prev => self%head
  self%n = n
  do i=2,n
    allocate (history)
    call history%init
    prev%next => history
    prev => history
  end do
END SUBROUTINE init

!> =============================================================================
!> Add the force on a particle from particles in the local particle patch, and
!> those in the list of nbor particle patches
!> =============================================================================
SUBROUTINE force_on (self, particle)
  class (particle_patch_t)                :: self
  class (particle_t), pointer        :: particle
  type (particle_history_t), pointer :: history
  class (particle_patch_t), pointer       :: nbor, other
  real(8)                            :: mass
  !-----------------------------------------------------------------------------
  ! Add forces from particles in the current patch
  !-----------------------------------------------------------------------------
  history => self%head
  particle%a = 0.0
  mass = 0.0
  do while (associated(history))
    if (history%id /= particle%id) then
      call history%add_particle_force_to (particle)
    end if
    mass = mass + history%mass
    history => history%next
  end do
  !$omp atomic write
  self%mass = mass
  !-----------------------------------------------------------------------------
  ! Add forces from particles in nbor patches
  !-----------------------------------------------------------------------------
  nbor => self%nbor
  do while (associated(nbor))
    history => self%head
    do while (associated(history))
      call history%add_particle_force_to (particle)
      history => history%next
    end do
    nbor => nbor%nbor
  end do
  !-----------------------------------------------------------------------------
  ! Finally, add forces from other patches, further away
  !-----------------------------------------------------------------------------
  other => self%other
  do while (associated(other))
    call self%add_patch_force_to (other, particle)
    other => other%other
  end do
END SUBROUTINE force_on

!> =============================================================================
!> Update the positions of all particles, maintaining their histories
!> =============================================================================
SUBROUTINE add_patch_force_to (self, other, particle)
  class (particle_patch_t)          :: self
  class (particle_patch_t), pointer :: other
  class (particle_t), pointer  :: particle
  real                         :: r(3), ri3
  !-----------------------------------------------------------------------------
  r = particle%p + (self%position - other%position)
  ri3 = (1./sqrt(sum(r**2)))**3
  particle%a = particle%a - (self%mass*ri3)*r
END SUBROUTINE add_patch_force_to

!> =============================================================================
!> Update the positions of all particles, maintaining their histories
!> =============================================================================
SUBROUTINE move_past (self, time)
  class (particle_patch_t)       :: self
  real                               :: time
  class (particle_t), pointer        :: particle
  type (particle_history_t), pointer :: history
  integer                            :: n, nh
  !-----------------------------------------------------------------------------
  n = self%n
  history => self%head
  nh = history%n
  do while (associated(history))
    particle => history%particle(nh)
    !print *, 'move_past', particle%id, particle%t
    do while (particle%time < time)
      !print*,'force_on'
      call self%force_on (particle)
      !print*,'kick'
      call particle%kick
      !print*,'drift'
      call particle%drift
      call self%force_on (particle)
      call particle%kick
      call history%rotate
      !history%particle(n) => particle
      history%moves = history%moves+1
    end do
    history => history%next
  end do
END SUBROUTINE move_past

!> =============================================================================
!> Randomize the particle velocities
!> =============================================================================
SUBROUTINE randomize (self)
  class (particle_patch_t)                :: self
  type (particle_history_t), pointer :: history
  !-----------------------------------------------------------------------------
  history => self%head
  do while (associated(history))
    call history%randomize
    history => history%next
  end do
END SUBROUTINE randomize

!> =============================================================================
!> Print information about the particles in the particle patch
!> =============================================================================
SUBROUTINE print (self)
  class (particle_patch_t)                :: self
  type (particle_history_t), pointer :: history
  integer(8):: moves
  !-----------------------------------------------------------------------------
  moves = 0
  history => self%head
  do while (associated(history))
    !call history%print
    moves = moves+history%moves
    history => history%next
  end do
  print *,'total moves:', real(moves)
  self%moves = moves
END SUBROUTINE print

END MODULE particle_patch_mod
