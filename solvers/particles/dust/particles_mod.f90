!===============================================================================
!> Particle data types for dust-gas interaction:
!>
!>   particle_t         : data type for particles
!>   particles_t        : data type for particle lists
!>        bins_t        : Cartesian collection of particle lists
!>
!> Optimized for speed of list operations by adding permanent head and tail
!===============================================================================
MODULE particles_mod
  USE io_mod
  implicit none
  private
  !--- link structure for particles ---
  type, public:: particle_t
    type(particle_t), pointer:: next, prev  ! 2*8 = 16
    integer(8):: id = -1                    !        8
    real:: r(3), v(3), w, rpad              ! 7*4 = 28
    integer(2):: q(3), bin                  ! 4*2 =  8
    integer(2):: flag = 0, ipad             !        2
  end type                                  ! sum:  62  
  !--- collection of particles ---
  type, public:: particles_t
    type(particle_t), pointer:: head => null(), tail => null()
    integer:: n = 0
  contains
    procedure:: init
    procedure:: append
    !procedure:: prepend
    procedure:: remove
    !procedure:: transfer
    procedure:: pop
  end type
  !--- cartesian arrangement of particle bins ---
  type, public:: bins_t
    class(particles_t), pointer:: particles(:,:,:,:)
    integer:: n
  end type
CONTAINS

!===============================================================================
!> Initialize particle list
!===============================================================================
SUBROUTINE init (self)
  class(particles_t):: self
  allocate (self%head, self%tail)
  self%head%next => self%tail
  self%tail%prev => self%head
END SUBROUTINE init

!===============================================================================
!> Append to the list
!===============================================================================
SUBROUTINE append (self, p)
  class(particles_t):: self
  type(particle_t), pointer:: p
  !.............................................................................
  p%next => self%tail
  p%prev => self%tail%prev
  p%next%prev => p
  p%prev%next => p
  self%n = self%n + 1
END SUBROUTINE append

!===============================================================================
!> Prepend to the list
!===============================================================================
SUBROUTINE prepend (self, p)
  class(particles_t):: self
  type(particle_t), pointer:: p
  !.............................................................................
  p%prev => self%head
  p%next => self%head%next
  p%next%prev => p
  p%prev%next => p
  self%n = self%n + 1
END SUBROUTINE prepend

!===============================================================================
!> Remove given particle
!===============================================================================
SUBROUTINE remove (self, p)
  class(particles_t):: self
  type(particle_t), pointer:: p
  p%next%prev => p%prev
  p%prev%next => p%next
  self%n = self%n - 1
END SUBROUTINE remove

!===============================================================================
!> Remove the last particle appended
!===============================================================================
FUNCTION pop(self) RESULT (p)
  class(particles_t):: self
  type(particle_t), pointer:: p
  if (self%n > 0) then
    p => self%tail%prev
    self%tail%prev => p%prev
    self%n = self%n - 1
  else
    allocate (p)
    call io%bits_mem (storage_size(p), 1, 'p')
  end if
END FUNCTION pop

!===============================================================================
!> Trasnfer all particles to a pool
!===============================================================================
SUBROUTINE transfer (self, pool)
  class(particles_t):: self
  class(particles_t), pointer:: pool
  type(particle_t), pointer:: p, next
  integer:: i
  p => self%head%next
  do i=1,self%n
    next => p%next
    call self%remove (p)
    call pool%append (p)
    p => next
  end do
END SUBROUTINE transfer
END MODULE particles_mod
