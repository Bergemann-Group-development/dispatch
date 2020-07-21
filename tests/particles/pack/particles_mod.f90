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
  USE mpi_buffer_mod
  implicit none
  private
  !--- link structure for particles ---
  type, public:: particle_t
    type(particle_t), pointer:: next, prev  ! 2*8 = 16
    real(8):: t                             !        8
    integer(2):: q(4)                       ! 4*2 =  8
    real:: r(3), v(3), w                    ! 7*4 = 28
  contains
    procedure:: init => init_particle
    procedure:: pack => pack_particle
  end type                                  ! sum:  60
  type, extends(particle_t):: dressed_particle_t
    integer(8):: id
    real:: vg(3), d
  end type
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
    procedure:: pack => pack_particles
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
END SUBROUTINE init_particle

!===============================================================================
SUBROUTINE pack_particle (self, mpi_buffer)
  class(particle_t):: self
  type(mpi_buffer_t):: mpi_buffer
  call mpi_buffer%append (self%r)
  call mpi_buffer%append (self%v)
  call mpi_buffer%append (self%q)
  call mpi_buffer%append (self%w)
  call mpi_buffer%append (self%t)
END SUBROUTINE pack_particle

!===============================================================================
!> Initialize particle list
!===============================================================================
SUBROUTINE init (self, id)
  class(particles_t):: self
  integer:: id(5)
  allocate (self%head, self%tail)
  self%id = id
  self%n = 0
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
      call self%dump
      call io%abort ('particles_t%validate')
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
    write (stdlog,'(i4,3i10)') j
    p => p%next
  end do
  flush (stdlog)
END SUBROUTINE dump

!===============================================================================
SUBROUTINE pack_particles (self, mpi_buffer)
  class(particles_t):: self
  class(particle_t), pointer:: p
  type(mpi_buffer_t):: mpi_buffer
  integer(2), allocatable:: qp(:,:), qs(:,:)
  integer(8), allocatable:: id(:)
  real, allocatable:: rp(:,:), vp(:,:), rs(:,:), vs(:,:), vg(:,:), ds(:)
  real, allocatable:: wp(:), ws(:)
  real(8), allocatable:: tp(:), ts(:)
  integer:: i, n, np, ns
  !-----------------------------------------------------------------------------
  n = self%n
  np = 0
  ns = 0
  p => self%head%next
  !-----------------------------------------------------------------------------
  ! Count categories
  !-----------------------------------------------------------------------------
  do i=1,n
    select type (p)
    type is (dressed_particle_t)
      ns = ns+1
    type is (particle_t)
      np = np+1
    end select
    p => p%next
  end do
  !-----------------------------------------------------------------------------
  ! Allocate and fill buffers
  !-----------------------------------------------------------------------------
  if (np > 0) &
    allocate (rp(3,np), vp(3,np), qp(4,np), tp(np), wp(np))
  if (ns > 0) &
    allocate (rs(3,ns), vs(3,ns), qs(4,ns), ts(ns), ws(ns), &
    vg(3,ns), ds(ns), id(ns))
  np = 0
  ns = 0
  p => self%head%next
  do i=1,n
    select type (p)
    type is (dressed_particle_t)
      ns = ns+1
      rs(:,ns) = p%r        ! float position
      qs(:,ns) = p%q        ! integer position
      vs(:,ns) = p%v        ! particle velocity
      vg(:,ns) = p%vg       ! gas velocity
      ds(ns)   = p%d        ! gas density
      ts(ns)   = p%t        ! particle time
      ws(ns)   = p%w        ! particle weight
      id(ns)   = p%id       ! particle id
    type is (particle_t)
      np = np+1
      rp(:,np) = p%r
      vp(:,np) = p%v
      qp(:,np) = p%q
      tp(np)   = p%t
      wp(np)   = p%w
    end select
    p => p%next
  end do
  !-----------------------------------------------------------------------------
  ! Pack into mpi_buffer
  !-----------------------------------------------------------------------------
  call mpi_buffer%append (np)         ! npp
  if (np > 0) then
    call mpi_buffer%append (rp)         ! 3*n
    call mpi_buffer%append (vp)         ! 3*n
    call mpi_buffer%append (qp)         ! 2*n
    call mpi_buffer%append (wp)         ! 1*n
    call mpi_buffer%append (tp)         ! 2*n = 11*n + npp
    deallocate (rp, vp, qp, tp, wp)
  end if
  call mpi_buffer%append (ns)
  if (ns > 0) then
    call mpi_buffer%append (rs)
    call mpi_buffer%append (vs)
    call mpi_buffer%append (qs)
    call mpi_buffer%append (ws)
    call mpi_buffer%append (ts)
    call mpi_buffer%append (ds)         ! 1*n
    call mpi_buffer%append (vg)         ! 3*n
    call mpi_buffer%append (id)         ! 2*n = 6*n
    deallocate (id,rs,vs,qs,ts,ws,vg,ds)
  end if
END SUBROUTINE pack_particles

END MODULE particles_mod
