!===============================================================================
!> Program to test functionality and speed of particle lists
!>
!> Data types:
!>  particle_t: add particle data
!> tparticle_t: add an 8-byte ID
!> particles_t: collection of class(plink_t)
!>      bins_t: Cartesian collection of class(particles_t)
!===============================================================================
MODULE bins_mod
  USE io_mod
  implicit none
  private
  !--- link structure for particles ---
  type, public:: particle_t
    type(particle_t), pointer:: next, prev
    real:: r(3), v(3), w
    integer:: q(3)
    integer(8):: id = 0
  end type
  !--- collection of particles ---
  type, public:: particles_t
    type(particle_t), pointer:: head => null(), tail => null()
    integer:: n = 0
  contains
    procedure:: init => init_particles
    procedure:: append
    procedure:: prepend
    procedure:: remove
  end type
  !--- cartesian arrangement of particle bins ---
  type, public:: bins_t
    class(particles_t), pointer:: particles(:,:,:,:)
  end type
CONTAINS

SUBROUTINE init_particles (self)
  class(particles_t):: self
  allocate (self%head, self%tail)
  self%head%next => self%tail
  self%tail%prev => self%head
END SUBROUTINE init_particles

SUBROUTINE append (self, p)
  class(particles_t):: self
  class(particle_t), pointer:: p
  !.............................................................................
  p%next => self%tail
  p%prev => self%tail%prev
  p%next%prev => p
  p%prev%next => p
  self%n = self%n + 1
END SUBROUTINE append

SUBROUTINE prepend (self, p)
  class(particles_t):: self
  class(particle_t), pointer:: p
  !.............................................................................
  p%prev => self%head
  p%next => self%head%next
  p%next%prev => p
  p%prev%next => p
  self%n = self%n + 1
END SUBROUTINE prepend

SUBROUTINE remove (self, p)
  class(particles_t):: self
  class(particle_t), pointer:: p
  p%next%prev => p%prev
  p%prev%next => p%next
  self%n = self%n - 1
END SUBROUTINE remove
END MODULE bins_mod

PROGRAM test
  USE io_mod
  USE bins_mod
  USE omp_lib
  implicit none
  integer:: i, n, bytes, s, iter    
  real(8):: wc
  real:: gb
  class(particle_t), pointer:: p, next
  type(particles_t):: particles
  type(bins_t):: bins
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call particles%init
  !-----------------------------------------------------------------------------
  ! Add & remove one GB of plinks
  !-----------------------------------------------------------------------------
  s = 10
  gb = 5.
  print *, 'check top! -- basic size'
  call sleep (s)
  allocate (p)
  bytes = storage_size(p)/8
  deallocate (p)
  n = gb*1024_8**3/bytes
  print 1, bytes, ' bytes', n*1e-6, 'million plinks/GB -- check top!'
  1 format(/,i6,a,f6.1,2x,a)
  do iter=1,2
    wc = omp_get_wtime()
    do i=1,n
      allocate(p)
      call particles%append (p)
      call io%bits_mem (storage_size(p), +1, 'p')
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' allocate + append'
    2 format(i10,9p,f8.1," ns/p -- ",a)
    call sleep (s)
    s = 0
    !---------------------------------------------------------------------------
    wc = omp_get_wtime()
    p => particles%head%next
    do i=1,n
      next => p%next
      p%id = i
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' empty loop'
    !---------------------------------------------------------------------------
    wc = omp_get_wtime()
    p => particles%head%next
    do i=1,n
      next => p%next
      select type (p)
      class is (particle_t)
        p%id = i
      end select
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' select'
    !---------------------------------------------------------------------------
    wc = omp_get_wtime()
    p => particles%head%next
    do i=1,n
      next => p%next
      if (associated(p)) p%id = i
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' associated test'
    !---------------------------------------------------------------------------
    wc = omp_get_wtime()
    p => particles%head%next
    do i=1,n
      next => p%next
      if (associated(p,particles%head)) p%id = i
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' 2-associated test'
    !---------------------------------------------------------------------------
    wc = omp_get_wtime()
    p => particles%head%next
    do i=1,n
      next => p%next
      call particles%remove (p)
      call particles%prepend (p)
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' remove + prepend'
    !---------------------------------------------------------------------------
    wc = omp_get_wtime()
    p => particles%head%next
    do i=1,n
      next => p%next
      call particles%remove (p)
      call io%bits_mem (storage_size(p), -1, 'p')
      deallocate(p)
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' remove + dealloc'
  end do
  !-----------------------------------------------------------------------------
  ! Add & remove one GB of plinks
  !-----------------------------------------------------------------------------
  allocate(p)
  bytes = storage_size(p)/8
  deallocate (p)
  n = gb*1024_8**3/bytes
  print 1, bytes, ' bytes', n*1e-6, 'million particles/GB'
  do iter=1,2
    wc = omp_get_wtime()
    do i=1,n
      allocate(p)
      call particles%append (p)
      call io%bits_mem (storage_size(p), +1, 'p')
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' allocate + append'
    call sleep (s)
    !---------------------------------------------------------------------------
    wc = omp_get_wtime()
    p => particles%head%next
    do i=1,n
      next => p%next
      call particles%remove (p)
      call io%bits_mem (storage_size(p), -1, 'p')
      deallocate(p)
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' remove + dealloc'
  end do
  !-----------------------------------------------------------------------------
  ! Add & remove one GB of plinks
  !-----------------------------------------------------------------------------
  allocate (p)
  bytes = storage_size(p)/8
  deallocate (p)
  n = gb*1024_8**3/bytes
  print 1, bytes, ' bytes', n*1e-6, 'million tparticles/GB -- with bits_mem'
  do iter=1,2
    wc = omp_get_wtime()
    do i=1,n
      allocate (p)
      call particles%append (p)
      call io%bits_mem (storage_size(p), +1, 'p')
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' allocate + append'
    call sleep (s)
    !---------------------------------------------------------------------------
    p => particles%head%next
    wc = omp_get_wtime()
    do i=1,n
      next => p%next
      call particles%remove (p)
      call io%bits_mem (storage_size(p), -1, 'p')
      deallocate(p)
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' remove + dealloc'
  end do
  !-----------------------------------------------------------------------------
  ! Add & remove one GB of plinks
  !-----------------------------------------------------------------------------
  allocate (p)
  bytes = storage_size(p)/8
  deallocate (p)
  n = gb*1024_8**3/bytes
  print 1, bytes, ' bytes', n*1e-6, 'million tparticles/GB -- no bits_mem'
  do iter=1,2
    wc = omp_get_wtime()
    do i=1,n
      allocate (p)
      call particles%append (p)
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' allocate + append'
    call sleep (s)
    !---------------------------------------------------------------------------
    p => particles%head%next
    wc = omp_get_wtime()
    do i=1,n
      next => p%next
      call particles%remove (p)
      deallocate(p)
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' remove + dealloc'
  end do
  !-----------------------------------------------------------------------------
  ! Add & remove one GB of plinks
  !-----------------------------------------------------------------------------
  allocate (p)
  bytes = storage_size(p)/8
  deallocate (p)
  n = gb*1024_8**3/bytes
  print 1, bytes, ' bytes', n*1e-6, 'million tparticles/GB'
  do iter=1,2
      wc = omp_get_wtime()
    do i=1,n
      allocate (p)
      call particles%append (p)
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' allocate + append'
    call sleep (s)
    !---------------------------------------------------------------------------
    p => particles%head%next
    wc = omp_get_wtime()
    do i=1,n
      next => p%next
      call particles%remove (p)
      p => next
    end do
    wc = omp_get_wtime() - wc
    print 2, particles%n, real(wc/n), ' remove'
  end do
END PROGRAM
