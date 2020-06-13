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
  USE omp_mod
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
    integer(kind=2) :: flag=0                   ! flag / bits
  contains
    procedure:: init
  end type
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  type, public:: sorted_t
    class(particle_t), pointer:: head => null()
    class(particle_t), pointer:: tail => null()
    integer:: n=0
  end type
  !-----------------------------------------------------------------------------
  ! List of particles
  !-----------------------------------------------------------------------------
  type, public:: particles_t
    class(particle_t), pointer:: head => null() ! pointer to start of list
    class(particle_t), pointer:: tail => null() ! pointer to end of list
    class(sorted_t), pointer:: sorted(:,:,:,:) => null()
    integer(8):: n=0, max=0                     ! number of particles in list
    real(8):: dtime                             ! collective dtime
    type(lock_t), pointer:: lock => null()      ! OMP lock
  contains
    procedure:: init => init_particles
    procedure:: prepend
    procedure:: append
    procedure:: remove
    procedure:: dealloc
    procedure:: add_sorted
    procedure:: remove_sorted
  end type
  integer:: verbose=0
  integer(8):: id=0
  !$omp threadprivate (id)
CONTAINS

!===============================================================================
!> Initialize a particle
!===============================================================================
SUBROUTINE init (self)
  class(particle_t):: self
  !.............................................................................
  self%id = 1 + omp%thread + (mpi%rank + mpi%size*id)*omp%nthreads
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
SUBROUTINE dealloc (self)
  class(particles_t):: self
  !.............................................................................
  class (particle_t), pointer:: item, next
  integer:: n
  !-----------------------------------------------------------------------------
  item => self%head
  n = 0
  do while (associated(item))
    n = n+1
    next => item%next
    deallocate (item)
    item => next
  end do
  if (n > 0) then
    call io%warning ('particles_t%deallocate: particles remained')
  end if
  call omp_lock%destroy (self%lock)
call self%init
END SUBROUTINE dealloc

!===============================================================================
SUBROUTINE add_sorted (self, p, i1, i2, i3, i4)
  class(particles_t):: self
  class(particle_t), pointer:: p
  class(sorted_t), pointer:: sorted
  integer(kind=2):: i1, i2, i3, i4
  !-----------------------------------------------------------------------------
  if (.not. associated(self%sorted)) then
    write(stdout,*) 'add_sorted: not associated', i1, i2, i3, i4
  end if
  sorted => self%sorted(i1,i2,i3,i4)
  !-----------------------------------------------------------------------------
  ! If the sorted box is empty, the particle becomes the only one in the box,
  ! but since it was removed from the particle list it needs to also be linked
  ! into the list again
  !-----------------------------------------------------------------------------
  if (sorted%n == 0) then
    call self%append (p)
    sorted%tail => p
    sorted%head => p
  else
    if (associated(sorted%tail,self%tail)) then     ! tail of ALL particles
      nullify(p%next)                               ! maintain null at end
      self%tail => p                                ! new tail
    else
      p%next => sorted%tail%next                    ! connect to next cell
      p%next%prev => p                              ! point back past it
    end if
    p%prev => sorted%tail                           ! attach to the tail
    p%prev%next => p                                ! point fwd past it
    sorted%tail => p
    self%n = self%n + 1                             ! increment total count
  end if
  sorted%n = sorted%n + 1                           ! increment box count
END SUBROUTINE add_sorted

!===============================================================================
SUBROUTINE remove_sorted (self, p, i1, i2, i3, i4)
  class(particles_t):: self
  class(particle_t), pointer:: p
  class(sorted_t), pointer:: sorted
  integer(kind=2):: i1, i2, i3, i4
  !-----------------------------------------------------------------------------
  if (.not. associated(self%sorted)) then
    write(stdout,*) 'remove_sorted: not associated', i1, i2, i3, i4
  end if
  sorted => self%sorted(i1,i2,i3,i4)
  !--- if there is a particle after this one, change its %prev to p%prev
  if (associated(p%next)) then
    p%next%prev => p%prev
  end if
  !--- if there is a particle before this one, change its %next to p%next
  if (associated(p%prev)) then
    p%prev%next => p%next
  end if
  !--- if the remove particle is the tail in a cell, repoint the tail
  if (associated(p,sorted%tail)) then
    sorted%tail => p%prev
    !--- if the removed particle is alse the very last tail, reassign tail
    if (associated(p,self%tail)) then
      self%tail => p%prev
      !--- for good measure, make sure the tail has no next
      nullify (self%tail%next)
    end if
  end if
  !--- if the remove particle is the head in a cell, repoint the head
  if (associated(p,sorted%head)) then
    sorted%head => p%next
    !--- if the removed particle is alse the very first head, reassign head
    if (associated(p,self%head)) then
      self%head => p%next
      !--- for good measure, make sure the head has no prev
      nullify (self%head%prev)
    end if
  end if
  sorted%n = sorted%n - 1                                ! decrement box count
  if (sorted%n == 0) then                                ! if box empty
    nullify (sorted%head, sorted%tail)                   ! reset head & tail
  end if
  self%n =  self%n - 1                                   ! decrement total
END SUBROUTINE remove_sorted

END MODULE particles_mod
