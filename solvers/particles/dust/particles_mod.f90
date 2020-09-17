!===============================================================================
!> Particle data types for dust-gas interaction:
!>
!>   particle_t         : data type for particles
!>   vparticle_t        : data type for visualization particles
!>   particles_t        : data type for particle lists
!>        bins_t        : Cartesian collection of particle lists
!>
!> Optimized for speed of list operations by adding permanent head and tail
!>
!> To include particles in the MPI communication, one could send the exports
!> particle bins along with the patch data, but this would be problematic,
!> since the exports particle bins also contain particles that don't belong
!> to the MPI nbor, but to local (boundary or internal) nbors.  The most
!> common MPI nbor has only one face in common with the sender, and hence it
!> would be a 6 times overkill to send all particle, and delete or never
!> append the irrelevant ones.
!>
!> But the alternative means to send different buffers to different nbors,
!> which is also not a good alternative. So, at least initially we should
!> accept the simpler alternative
!===============================================================================
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE mpi_mod
  USE mpi_buffer_mod
  implicit none
  private
  !--- link structure for particles ---
  type, public:: particle_t
    class(particle_t), pointer:: next, prev ! 2*8 = 16
    real(8):: t                             !        8
    integer(2):: q(4)                       ! 4*2 =  8
    real:: r(3), v(3), w                    ! 7*4 = 28
  contains
    procedure:: init => init_particle
    procedure:: pack => pack_particle
  end type                                  ! sum:  60
  integer, parameter:: vparticle_key=1
  type, public, extends(particle_t):: vparticle_t
    integer(8):: id
    real:: vg(3), dg, dd
  contains
    procedure:: init => init_vparticle
    procedure:: pack => pack_vparticle
    procedure:: output => output_vparticle
  end type
  !--- collection of particles ---
  integer, parameter:: particles_key=2
  type, public:: particles_t
    class(particle_t), pointer:: head => null(), tail => null()
    integer:: n = 0
    integer:: type = 0
    integer:: id(5)
    real(8):: momentum(3)=0d0
  contains
    procedure:: init => init_particles
    procedure:: dealloc
    procedure:: append
    !procedure:: prepend
    procedure:: remove
    procedure:: transfer
    procedure:: pop
    procedure:: validate
    procedure:: dump
    procedure:: pack => pack_particles
    procedure:: statistics
    procedure:: total_momentum
    procedure:: density
    procedure:: dispersion
  end type
  !--- cartesian arrangement of particle bins ---
  type, public:: bins_t
    class(particles_t), pointer:: particles(:,:,:,:)
    integer:: n
  end type
  integer(8), save:: id=0
  !$omp threadprivate(id)
  type(particles_t), public :: pool, vpool    ! particle pool
  !$omp threadprivate (pool, vpool)           ! for each thread
CONTAINS

!===============================================================================
!> Initialize particle, with unique ID, based on the threadprivate global
!> thread number omp%thread_global
!===============================================================================
SUBROUTINE init_particle (self)
  class(particle_t):: self
END SUBROUTINE init_particle
SUBROUTINE init_vparticle (self)
  class(vparticle_t):: self
  self%id = 1 + omp%thread_global + id*omp%nthreads_global
  id = id+1
END SUBROUTINE init_vparticle

!===============================================================================
SUBROUTINE pack_particle (self, buffer)
  class(particle_t):: self
  type(mpi_buffer_t):: buffer
  call buffer%append (self%r)
  call buffer%append (self%v)
  call buffer%append (self%q)
  call buffer%append (self%w)
  call buffer%append (self%t)
END SUBROUTINE pack_particle

!===============================================================================
SUBROUTINE pack_vparticle (self, buffer)
  class(vparticle_t):: self
  type(mpi_buffer_t):: buffer
  call buffer%append (self%r)
  call buffer%append (self%v)
  call buffer%append (self%q)
  call buffer%append (self%w)
  call buffer%append (self%t)
  call buffer%append (self%dg)
  call buffer%append (self%vg)
  call buffer%append (self%id)
END SUBROUTINE pack_vparticle

!===============================================================================
SUBROUTINE output_vparticle (self)
  class(vparticle_t) :: self
  integer            :: n, j
  integer, save      :: itimer=0
  !-----------------------------------------------------------------------------
  write (io_unit%tmp) vparticle_key
  write (io_unit%tmp) self%q, self%r, self%v, self%vg, self%w, self%dg, &
    self%t, self%id
END SUBROUTINE output_vparticle

!===============================================================================
!> Initialize particle list
!===============================================================================
SUBROUTINE init_particles (self, id, type)
  class(particles_t):: self
  integer:: id(5), type
  character(len=120), save:: ids= &
  '$Id$
  !-----------------------------------------------------------------------------
  call trace%print_id (ids)
  allocate (self%head, self%tail)
  self%id = id
  self%type = type
  self%n = 0
  self%head%next => self%tail
  self%tail%prev => self%head
END SUBROUTINE init_particles

!===============================================================================
!> Deallocate particle list
!===============================================================================
SUBROUTINE dealloc (self)
  class(particles_t):: self
  !-----------------------------------------------------------------------------
  deallocate (self%head, self%tail)
END SUBROUTINE dealloc

!===============================================================================
!> Append to the list
!===============================================================================
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
  class(particle_t), pointer:: p
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
  class(particle_t), pointer:: p
  type(particle_t), pointer:: p0
  type(vparticle_t), pointer:: p1
  if (self%n > 0) then
    p => self%tail%prev
    p%next%prev => p%prev
    p%prev%next => p%next
    self%n = self%n - 1
  else
    if (self%type == 0) then
      allocate (p0)
      call io%bits_mem (storage_size(p0), 1, 'p0')
      call p0%init
      p => p0
    else
      if (mpi%size == 1) then
        write (stderr,*) 'ERROR:', self%id
        allocate (p0)
        call io%bits_mem (storage_size(p0), 1, 'p0')
        call p0%init
        p => p0
      else
        allocate (p1)
        call io%bits_mem (storage_size(p1), 1, 'p1')
        call p1%init
        p => p1
      end if
    end if
  end if
END FUNCTION pop

!===============================================================================
!> Trasnfer all particles to a pool
!===============================================================================
SUBROUTINE transfer (self, pool)
  class(particles_t):: self
  type(particles_t), pointer:: pool
  class(particle_t), pointer:: p, next
  integer:: j
  !.............................................................................
  p => self%head%next
  do j=1,self%n
    next => p%next
    call self%remove (p)
    call pool%append (p)
    if (abs(maxval(p%q)) > 100) &
      write (stdlog,*) 'transfer error', j, self%n, p%q
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
    p => p%next
  end do
  flush (stdlog)
END SUBROUTINE dump

!===============================================================================
SUBROUTINE pack_particles (self, mpi_buffer)
  class(particles_t):: self
  type(mpi_buffer_t):: mpi_buffer
  class(particle_t), pointer:: p
  integer(2), pointer, contiguous:: qp(:,:), qs(:,:)
  integer(8), pointer, contiguous:: id(:)
  real(8),    pointer, contiguous:: tp(:), ts(:)
  real,       pointer, contiguous:: rp(:,:), vp(:,:), wp(:)
  real,       pointer, contiguous:: rs(:,:), vs(:,:), ws(:)
  integer:: i, n, np, ns
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pack_particles')
  n = self%n
  np = 0
  ns = 0
  p => self%head%next
  !-----------------------------------------------------------------------------
  ! Count categories
  !-----------------------------------------------------------------------------
  do i=1,n
    select type (p)
    type is (vparticle_t)
      ns = ns+1
    type is (particle_t)
      np = np+1
    end select
    p => p%next
  end do
  !print *,'counts(1):',np,ns
  !-----------------------------------------------------------------------------
  ! Allocate and fill buffers
  !-----------------------------------------------------------------------------
  if (np > 0) &
    allocate (rp(3,np), vp(3,np), qp(4,np), tp(np), wp(np))
  if (ns > 0) &
    allocate (rs(3,ns), vs(3,ns), qs(4,ns), ts(ns), ws(ns), id(ns))
  np = 0
  ns = 0
  p => self%head%next
  do i=1,n
    select type (p)
    type is (vparticle_t)
      ns = ns+1
      rs(:,ns) = p%r        ! float position
      qs(:,ns) = p%q        ! integer position
      vs(:,ns) = p%v        ! particle velocity
      ts(ns)   = p%t        ! particle time
      ws(ns)   = p%w        ! particle weight
      id(ns)   = p%id       ! particle id
    type is (particle_t)
      np = np+1
      rp(:,np) = p%r
      qp(:,np) = p%q
      vp(:,np) = p%v
      tp(np)   = p%t
      wp(np)   = p%w
    end select
    p => p%next
  end do
  !print *,'counts(2):',np,ns
  !-----------------------------------------------------------------------------
  ! Pack into mpi_buffer
  !-----------------------------------------------------------------------------
  call mpi_buffer%append (1)            ! normal particles
  call mpi_buffer%append (np)           ! np
  if (np > 0) then
    call mpi_buffer%append (rp)         ! 3*n
    call mpi_buffer%append (vp)         ! 3*n
    call mpi_buffer%append (qp)         ! 2*n
    call mpi_buffer%append (wp)         ! 1*n
    call mpi_buffer%append (tp)         ! 2*n = 11*n + npp
    deallocate (rp, vp, qp, tp, wp)
  end if
  call mpi_buffer%append (2)            ! dressed particles
  call mpi_buffer%append (ns)
  if (ns > 0) then
    call mpi_buffer%append (rs)
    call mpi_buffer%append (vs)
    call mpi_buffer%append (qs)
    call mpi_buffer%append (ws)
    call mpi_buffer%append (ts)
    call mpi_buffer%append (id)         ! 2*n = 6*n
    deallocate (id,rs,vs,qs,ts,ws)
  end if
  call trace%end ()
END SUBROUTINE pack_particles

!===============================================================================
!> Unpack particle buffer, from MPI or I/O
!===============================================================================
SUBROUTINE unpack_particles (self, mpi_buffer)
  class(particles_t)           :: self
  type(mpi_buffer_t)           :: mpi_buffer
  class(particle_t),      pointer:: p
  integer(2), pointer, contiguous:: qp(:,:)
  integer(8), pointer, contiguous:: id(:)
  real(8)   , pointer, contiguous:: tp(:)
  real      , pointer, contiguous:: rp(:,:), vp(:,:), wp(:)
  integer:: j, n, np, ns, type
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%unpack_particles')
  !-----------------------------------------------------------------------------
  ! Read mpi_buffer into particle arrays
  !-----------------------------------------------------------------------------
  call mpi_buffer%read (type)        ! particle type?
  call mpi_buffer%read (np)
  allocate (rp(3,np), vp(3,np), qp(4,np), tp(np), wp(np))
  call mpi_buffer%read (rp)
  call mpi_buffer%read (vp)
  call mpi_buffer%read (qp)
  call mpi_buffer%read (wp)
  call mpi_buffer%read (tp)
  if (type==0) then
    do j=1,np
      p => pool%pop()
      p%r = rp(:,j)
      p%q = qp(:,j)
      p%v = vp(:,j)
      p%t = tp(j)
      p%w = wp(j)
      call self%append (p)
    end do
  else
    allocate (id(ns))
    call mpi_buffer%read (id)
    do j=1,np
      p => vpool%pop()
      select type (p)
      class is (vparticle_t)
      p%r = rp(:,j)
      p%q = qp(:,j)
      p%v = vp(:,j)
      p%t = tp(j)
      p%w = wp(j)
      p%id = id(j)
      end select
      call self%append (p)
    end do
    deallocate (id)
  end if
  deallocate (rp, vp, qp, tp, wp)
  !.............................................................................
  call trace%end ()
END SUBROUTINE unpack_particles

!===============================================================================
SUBROUTINE statistics (self)
  class(particles_t):: self
  write (stdlog,'(a,i12,2x,a,i5)') &
    'particles_t%statistics:',id,'particles on thread',omp%thread_global
END SUBROUTINE statistics

!===============================================================================
!> Compute total momentum of particles
!===============================================================================
FUNCTION total_momentum (self) RESULT (sum)
  class(particles_t):: self
  real:: sum(3)
  type(particle_t), pointer:: p
  integer:: j
  !-----------------------------------------------------------------------------
  sum = 0.0
  p => self%head%next
  do j=1,self%n
    sum = sum + p%v*p%w
    p => p%next
  end do
END FUNCTION total_momentum

!===============================================================================
!> Compute total momentum of particles
!===============================================================================
FUNCTION density (self) RESULT (sum)
  class(particles_t):: self
  real:: sum
  type(particle_t), pointer:: p
  integer:: j
  !-----------------------------------------------------------------------------
  sum = 0.0
  p => self%head%next
  do j=1,self%n
    sum = sum + p%w
    p => p%next
  end do
END FUNCTION density

!===============================================================================
!> Compute total momentum of particles
!===============================================================================
FUNCTION dispersion (self) RESULT (rms)
  class(particles_t):: self
  real:: sum, vav(3), rms
  type(particle_t), pointer:: p
  integer:: j
  !-----------------------------------------------------------------------------
  if (self%n > 0) then
    !---------------------------------------------------------------------------
    ! Compute mean velocity
    !---------------------------------------------------------------------------
    sum = 0.0
    vav = 0.0
    p => self%head%next
    do j=1,self%n
      sum = sum + p%w
      vav = vav + p%w*p%v
      p => p%next
    end do
    vav = vav/sum
    !---------------------------------------------------------------------------
    ! Compute RMS velocity dispersion
    !---------------------------------------------------------------------------
    rms = 0.0
    p => self%head%next
    do j=1,self%n
      rms = rms + p%w*((p%v(1)-vav(1))**2 &
                      +(p%v(2)-vav(2))**2 &
                      +(p%v(3)-vav(3))**2)
      p => p%next
    end do
    rms = sqrt(rms/sum)
  else
    rms = 0.0
  end if
END FUNCTION dispersion

END MODULE particles_mod
