!===============================================================================
!> Program to test functionality and speed of particles, when using an array
!> of structures
!===============================================================================
MODULE dust_mod
  type part_t
    real:: p(3)=0.0, v(3)=1.0
    real:: m=1.0, s=0.0, t=0.0, dt=1.0
    type(part_t), pointer:: next => null()
  end type
CONTAINS

!===============================================================================
SUBROUTINE copy_to_buffer (dust, p, v, m, s, t, dt, nbin, ndim)
  implicit none
  class(part_t), pointer:: dust
  integer:: nbin, ndim
  real, dimension(nbin,ndim) :: p, v, r, f, a
  real, dimension(nbin)      :: m, s, t, dt, ra
  real(8), dimension(ndim)   :: r0
  integer:: id, ib
  !.............................................................................
  do ib=1,nbin
    p(ib,:) = dust%p
    v(ib,:) = dust%v
    m(ib)   = dust%m
    t(ib)   = dust%t
    dt(ib)  = dust%dt
    dust => dust%next
  end do
END

!===============================================================================
SUBROUTINE copy_fr_buffer (dust, p, v, m, s, t, dt, nbin, ndim)
  implicit none
  class(part_t), pointer:: dust
  integer:: nbin, ndim
  real, dimension(nbin,ndim) :: p, v, r, f, a
  real, dimension(nbin)      :: m, s, t, dt, ra
  real(8), dimension(ndim)   :: r0
  integer:: id, ib
  !.............................................................................
  do ib=1,nbin
    dust%p  = p(ib,:)
    dust%v  = v(ib,:)
    dust%m  = m(ib)
    dust%t  = t(ib)
    dust%dt = dt(ib)
    dust => dust%next
  end do
END
END MODULE dust_mod

!===============================================================================
!===============================================================================
PROGRAM test
  USE list_mod
  USE task_mod
  USE mpi_timer
  USE dust_mod
  USE planet_mod
  USE star_mod
  implicit none
  class(star_t), pointer:: star
  class(planet_t), pointer:: planet
  class(part_t), pointer:: head, dust, prev => null(), next
  integer:: ndim=3, nbin=64, npart=100000, ip, ib, id
  real, dimension(:,:), pointer:: p, v
  real, dimension(:,:), pointer:: a, f, r
  real, dimension(:), pointer:: m, s, t, dt, ra
  !-----------------------------------------------------------------------------
  ! Allocate buffer variables
  !-----------------------------------------------------------------------------
  allocate (star, planet)
  allocate (p(nbin,ndim))
  allocate (v(nbin,ndim))
  allocate (a(nbin,ndim))
  allocate (f(nbin,ndim))
  allocate (r(nbin,ndim))
  allocate (m(nbin))
  allocate (s(nbin))
  allocate (t(nbin))
  allocate (dt(nbin))
  allocate (ra(nbin))
  !-----------------------------------------------------------------------------
  ! Allocate dust
  !-----------------------------------------------------------------------------
  allocate (head)
  prev => head
  do ip=2,npart*nbin
    allocate (dust)
    prev%next => dust
    prev => dust
  end do
  !-----------------------------------------------------------------------------
  ! Initialize positions
  !-----------------------------------------------------------------------------
  star%position = 1d0
  !-----------------------------------------------------------------------------
  ! Evolve with buffers
  !-----------------------------------------------------------------------------
  call tic
  dust => head
  do ip=1,npart
    prev => dust
    call copy_to_buffer (dust, p, v, m, s, t, dt, nbin, ndim)
    call move_particles (p, v, m, s, t, dt, ra, r, f, a, star%position, nbin, ndim)
    dust => prev
    call copy_fr_buffer (dust, p, v, m, s, t, dt, nbin, ndim)
  end do
  call timer%toc ("particle motion", npart*nbin)
END PROGRAM

!===============================================================================
SUBROUTINE move_particles (p, v, m, s, t, dt, ra, r, f, a, r0, nbin, ndim)
  implicit none
  integer:: nbin, ndim
  real, dimension(nbin,ndim) :: p, v, r, f, a
  real, dimension(nbin)      :: m, s, t, dt, ra
  real(8), dimension(ndim)   :: r0
  integer:: id, ib
  !.............................................................................
  ra = 0.0
  do id=1,ndim
  do ib=1,nbin
    r(ib,id) = p(ib,id) - r0(id)
    ra(ib) = ra(ib) + r(ib,id)**2
  end do
  end do
  ra = sqrt(ra)**3
  do id=1,ndim
  do ib=1,nbin
    f(ib,id) = -r(ib,id)/ra(ib)
    a(ib,id) =  f(ib,id)/m(ib)
  end do
  end do
  do id=1,ndim
  do ib=1,nbin
    v(ib,id) = v(ib,id) + a(ib,id)*dt(ib)
    p(ib,id) = p(ib,id) + v(ib,id)*dt(ib)
  end do
  end do
END SUBROUTINE
