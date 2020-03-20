!===============================================================================
!> Particle data types for dust-gas interaction:
!>
!>   particle_t         : attributes for particles
!>   particles_t        : attributes for particle lists
!===============================================================================
MODULE particles_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  integer, parameter:: nt=5
  !-----------------------------------------------------------------------------
  ! Particle data type
  !-----------------------------------------------------------------------------
  type, public:: particle_t
    type(particle_t), pointer:: next => null()  ! next particle
    type(particle_t), pointer:: prev => null()  ! previous particle
    real(8)         :: time = 0d0               ! current time
    real(8)         :: dt                       ! time step
    real(8)         :: t(nt)                    ! time slices
    real            :: w                        ! weight factor
    real            :: v(3)                     ! velocity (cell units / time)
    real            :: r(3,nt)                  ! real positions (cell units)
    integer(kind=2) :: q(3,nt)                  ! integer position (cell units)
    integer(kind=2) :: it=1                     ! time slice index
    integer(kind=2) :: bin=1                    ! size bin
  end type
  !-----------------------------------------------------------------------------
  ! List of particles
  !-----------------------------------------------------------------------------
  type, public:: particles_t
    class(particle_t), pointer:: head => null() ! pointer to start of list
    class(particle_t), pointer:: tail => null() ! pointer to end of list
    integer:: n=0                               ! number of particles in list
  contains
    procedure:: init
    procedure:: append
    procedure:: remove
    procedure:: delete
  end type
  integer:: verbose=0
CONTAINS

!===============================================================================
!> Initialize a list
!===============================================================================
SUBROUTINE init (self)
  class(particles_t):: self
  !.............................................................................
  nullify (self%head, self%tail)
  !$omp atomic write
  self%n = 0
END SUBROUTINE init

!===============================================================================
!> Add an item at the tail
!===============================================================================
SUBROUTINE append (self, new)
  class(particles_t)              :: self
  class(particle_t), pointer :: new
  !.............................................................................
  new%prev => self%tail
  nullify(new%next)
  if (.not.associated(self%head)) self%head => new
  if (associated(self%tail)) self%tail%next => new
  self%tail => new
  !$omp atomic
  self%n = self%n+1
END SUBROUTINE append

!===============================================================================
!> Remove an item from the list; note that the item is not deallocated!
!===============================================================================
SUBROUTINE remove (self, old)
  class(particles_t):: self
  class(particle_t), pointer:: old
  !.............................................................................
  if (associated(old%prev)) old%prev%next => old%next
  if (associated(old%next)) old%next%prev => old%prev
  if (associated(old, self%head)) then
    self%head => old%next
  end if
  if (associated(old, self%tail)) then
    self%tail => old%prev
  end if
  !$omp atomic
  self%n = self%n-1
END SUBROUTINE remove

!===============================================================================
!> Delete a double-linked list, including carried content
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
