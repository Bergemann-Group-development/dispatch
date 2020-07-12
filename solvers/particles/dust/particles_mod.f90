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
  USE omp_mod
  implicit none
  private
  !--- link structure for particles ---
  type, public:: particle_t
    type(particle_t), pointer:: next, prev  ! 2*8 = 16
    real(8):: t                             !        8
    integer(8):: id = -1                    !        8
    real:: r(3), v(3), w, rpad              ! 7*4 = 28
    integer(2):: q(3), bin                  ! 4*2 =  8
    integer(2):: flag = 0, ipad             !        2
  contains
    procedure:: init => init_particle
  end type                                  ! sum:  70  
  !--- collection of particles ---
  type, public:: particles_t
    type(particle_t), pointer:: head => null(), tail => null()
    integer:: n = 0
    integer:: id(5)
  contains
    procedure:: init
    procedure:: append
    !procedure:: prepend
    procedure:: remove
    procedure:: transfer
    procedure:: pop
    procedure:: validate
    procedure:: dump
  end type
  !--- cartesian arrangement of particle bins ---
  type, public:: bins_t
    class(particles_t), pointer:: particles(:,:,:,:)
    integer:: n
  end type
  integer(8), save:: id=0
  !$omp threadprivate(id)
CONTAINS

!===============================================================================
!> Initialize particle, with unique ID, based on the threadprivate global
!> thread number omp%thread_global
!===============================================================================
SUBROUTINE init_particle (self)
  class(particle_t):: self
  self%id = omp%thread_global + id*omp%nthreads_global
  id = id+1
END SUBROUTINE init_particle

!===============================================================================
!> Initialize particle list
!===============================================================================
SUBROUTINE init (self, id)
  class(particles_t):: self
  integer:: id(5)
  allocate (self%head, self%tail)
  self%id = id
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
!define DBG
#ifdef DBG
  if (p%id == -1 .or. .not.associated(p%next) &
                 .or. .not.associated(p%next%prev)) then
    write (stdlog,*)'remove error:', p%id
    call self%dump
    call io%warning ('not associated')
    write (stdlog,*) associated(p%next)
    write (stdlog,*) associated(p%next%prev)
  end if
  if (self%n == 0) then
    print *, self%id
    call io%warning ('self%n == 0')
  end if
#endif
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
    p%next%prev => p%prev
    p%prev%next => p%next
    self%n = self%n - 1
  else
    allocate (p)
    call io%bits_mem (storage_size(p), 1, 'p')
    call p%init
  end if
END FUNCTION pop

!===============================================================================
!> Trasnfer all particles to a pool
!===============================================================================
SUBROUTINE transfer (self, pool)
  class(particles_t):: self
  type(particles_t), pointer:: pool
  type(particle_t), pointer:: p, next
  integer:: j
  !.............................................................................
  p => self%head%next
  do j=1,self%n
    next => p%next
    call self%remove (p)
    call pool%append (p)
    p => next
  end do
END SUBROUTINE transfer

!===============================================================================
!> Validate particles_t
!===============================================================================
SUBROUTINE validate (self, label)
  class(particles_t):: self
  character(len=*):: label
  type(particle_t), pointer:: p
  integer:: j, k
  !.............................................................................
  p => self%head%next
  do j=1,self%n
    if (p%id < 0) then
      call self%dump
      call io%abort ('particles_t%validate')
    end if
    p => p%next
  end do
END SUBROUTINE validate

!===============================================================================
!> Dump particles_t
!===============================================================================
SUBROUTINE dump (self)
  class(particles_t):: self
  type(particle_t), pointer:: p
  integer:: j
  !.............................................................................
  write (stdlog,*) 'dump:', self%id
  p => self%head%next
  do j=1,self%n
    write (stdlog,'(i4,3i10)') j, p%id, p%prev%id, p%next%id
    p => p%next
  end do
  flush (stdlog)
END SUBROUTINE dump

END MODULE particles_mod
