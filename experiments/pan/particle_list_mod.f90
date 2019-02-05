!> =============================================================================
!> Particle list, extending a doubly-linked list.  Each particle maintains 
!> arrays with previous positions, velocities, and times, which are brought
!> along if/when it changes owner patch or rank.  The memory footprint of a
!> particle is 12 words for position, 12 words for velocity, and 8 words for
!> time, all together 32 words = 128 bytes (plus a few words for id and weight).
!> This could be reduced to half, by keeping only two previous positions and
!> velocities in the data type.
!> =============================================================================
MODULE particle_list_mod
  USE dll_mod
  USE omp_timer_mod
  implicit none
  private
  integer, parameter:: nt=4
  integer, save:: id=0
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: particle_t
    integer(8):: id=0
    real, allocatable:: r(:,:), dr(:,:)
    real(8), allocatable:: t(:)
    real:: w
    integer:: it, nt
  contains
    procedure:: init
  end type
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: particle_list_t
  contains
    procedure:: test
    procedure:: print
  end type
  !-----------------------------------------------------------------------------
CONTAINS

!> =============================================================================
!> Initialize a particle with a unique ID
!> =============================================================================
SUBROUTINE init (self)
  class (particle_t):: self
  !-----------------------------------------------------------------------------
  !$omp atomic
  id = id+1
  self%id = id
  self%it = 1
  self%nt = nt
  allocate (self%r(3,nt), self%dr(3,nt), self%t(nt))
END SUBROUTINE init

!> =============================================================================
!> Initialize a particle list
!> =============================================================================
SUBROUTINE test (self)
  class(particle_list_t):: self
  integer:: id, n
  class(particle_t), pointer:: p
  class(dll_node_t), pointer :: node, next, new
  real(8)::used
  logical:: flip
  !-----------------------------------------------------------------------------
  call self%init ('test')
  !-----------------------------------------------------------------------------
  ! Make a list with 100000 particles
  !-----------------------------------------------------------------------------
  n = 100000
  used = wallclock()
  do id=1,n
    allocate (p)
    call p%init
    p%r = 2.0
    p%dr = 1.0
    node => p
    call self%append (node)
  end do
  used = wallclock()-used
  print *, 'particle list n =', self%n, used/self%n
  !-----------------------------------------------------------------------------
  ! Remove every 2nd particle
  !-----------------------------------------------------------------------------
  flip = .false.
  used = wallclock()
  node => self%head
  do while (associated(node))
    next => node%next
    select type (node)
    type is (particle_t)
    p => node
    end select
    if (flip) then
      call self%remove (node)
      deallocate (node)
    end if
    flip = .not.flip
    node => next
  end do
  used = wallclock()-used
  print *, 'particle list n =', self%n, used/self%n
  !-----------------------------------------------------------------------------
  ! Add one new particle for every 2nd existing one
  !-----------------------------------------------------------------------------
  flip = .false.
  used = wallclock()
  node => self%head
  do while (associated(node))
    next => node%next
    allocate (p)
    new => p
    if (flip) then
      call self%insert_before (node, new)
    end if
    flip = .not.flip
    node => next
  end do
  used = wallclock()-used
  print *, 'particle list n =', self%n, used/self%n
END SUBROUTINE test

!> =============================================================================
!> Print a particle list
!> =============================================================================
SUBROUTINE print (self)
  class(particle_list_t):: self
  class(dll_node_t), pointer:: p
  p => self%head
  do while (associated(p))
    select type (p)
    class is (particle_t)
    print *, p%id, p%it, p%r(:,p%it)
    end select
    p => p%next
  end do
END SUBROUTINE print

END MODULE particle_list_mod
