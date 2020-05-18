!===============================================================================
!> Particle data types for dust-gas interaction:
!>
!>   particle_t         : data type for particles
!>   particles_t        : data type for particle lists
!===============================================================================
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE omp_lock_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Particle data type
  !-----------------------------------------------------------------------------
  type, public:: particle_t
    type(particle_t), pointer:: next => null()  ! next particle
    type(particle_t), pointer:: prev => null()  ! previous particle
    integer(8)      :: id = 0                   ! particle id
    real(8)         :: t = 0.0_8                ! time
    real(8)         :: dt                       ! time step
    real            :: w                        ! weight factor
    real            :: c                        ! drag coefficient
    real            :: v(3)                     ! velocity (cell units / time)
    real            :: dv(3)                    ! vgas - velocity
    real            :: a(3) = 0.0               ! acceleration
    real            :: d(3)                     ! drag
    real            :: r(3)                     ! real positions (cell units)
    integer(kind=2) :: q(3)                     ! integer position (cell units)
    integer(kind=2) :: bin=1                    ! size bin
    integer(kind=2) :: flag                     ! flag / bits
  contains
    procedure:: init
  end type
  !-----------------------------------------------------------------------------
  ! List of particles
  !-----------------------------------------------------------------------------
  type, public:: particles_t
    class(particle_t), pointer:: head => null() ! pointer to start of list
    class(particle_t), pointer:: tail => null() ! pointer to end of list
    integer(8):: n=0, max=0                     ! number of particles in list
    real(8):: dtime                             ! collective dtime
    type(lock_t), pointer:: lock => null()      ! OMP lock
  contains
    procedure:: init => init_particles
    procedure:: prepend
    procedure:: append
    procedure:: remove
    procedure:: delete
  end type
  integer:: verbose=0
  integer(8):: id=0
CONTAINS

!===============================================================================
!> Initialize a particle
!===============================================================================
SUBROUTINE init (self)
  class(particle_t):: self
  !.............................................................................
  self%id = 1 + mpi%rank + mpi%size*id
  !$omp atomic update
  id = id+1
END SUBROUTINE init

!===============================================================================
!> Initialize a particle list
!===============================================================================
SUBROUTINE init_particles (self)
  class(particles_t):: self
  !.............................................................................
  call trace%begin ('particles_t%init')
  nullify (self%head, self%tail)
  associate (n=>self%n)
  !$omp atomic write
  n = 0
  end associate
  allocate (self%lock)
  call self%lock%init ('particles')
  call trace%end ()
END SUBROUTINE init_particles

!===============================================================================
!> Add an particle at the tail
!===============================================================================
SUBROUTINE prepend (self, new)
  class(particles_t)              :: self
  class(particle_t), pointer :: new
  !.............................................................................
  nullify (new%prev)
  new%next => self%head
  if (.not.associated(self%tail)) self%tail => new
  if (associated(self%head)) self%head%next => new
  self%head => new
  associate (n=>self%n, max=>self%max)
  !$omp atomic
  n = n+1
  !$omp atomic write
  max = id
  end associate
END SUBROUTINE prepend

!===============================================================================
!> Add an particle at the tail
!===============================================================================
SUBROUTINE append (self, new, atomic)
  class(particles_t):: self
  class(particle_t), pointer :: new
  logical, optional:: atomic
  !.............................................................................
  nullify (new%next)
  new%prev => self%tail
  if (self%n > 0) then
    self%tail%next => new
  else
    self%head => new
  end if
  self%tail => new
  if (present(atomic)) then
    associate (n=>self%n, max=>self%max)
    !$omp atomic
    n = n+1
    !$omp atomic write
    max = id
    end associate
  else
    self%n = self%n+1
    self%max = id
  end if
END SUBROUTINE append

!===============================================================================
!> Remove a particle from the list (the particle is not deallocated)
!===============================================================================
SUBROUTINE remove (self, old)
  class(particles_t):: self
  class(particle_t), pointer:: old
  !.............................................................................
  call self%lock%set ('particles_t%remove')
  if (associated(old%prev)) old%prev%next => old%next
  if (associated(old%next)) old%next%prev => old%prev
  if (associated(old, self%head)) then
    self%head => old%next
  end if
  if (associated(old, self%tail)) then
    self%tail => old%prev
  end if
  nullify (old%prev, old%next)
  associate (n=>self%n)
  !$omp atomic
  n = n-1
  end associate
  call self%lock%unset ('particles_t%remove')
END SUBROUTINE remove

!===============================================================================
!> Delete a particle list, including particles
!===============================================================================
SUBROUTINE delete (self)
  class(particles_t):: self
  !.............................................................................
  class (particle_t), pointer:: item, next
  !-----------------------------------------------------------------------------
  item => self%head
  do while (associated(item))
    next => item%next
    deallocate (item)
    item => next
  end do
  nullify (self%head)
  nullify (self%tail)
  self%n = 0
  call self%init
END SUBROUTINE delete

END MODULE particles_mod
