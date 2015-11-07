!> =============================================================================
!> Program to test functionality and speed of particles
!> =============================================================================

PROGRAM test
  USE list_mod
  USE task_mod
  USE timer_mod
  implicit none
  class(star_t), pointer:: star
  class(planet_t), pointer:: planet
  class(particles_t), pointer:: dust
  integer:: ndim=3, nbin=64, npart=1000000, ip, ib, id
  real, dimension(:,:), pointer:: p, v, a, f, r
  real, dimension(:),   pointer:: m, s, t, dt, ra

  !-----------------------------------------------------------------------------
  ! Allocate variables
  !-----------------------------------------------------------------------------
  allocate (dust, star, planet)
  allocate (dust%particles(npart))
  allocate (f(nbin,3))
  allocate (a(nbin,3))
  allocate (r(nbin,3))
  allocate (m(nbin))
  allocate (s(nbin))
  allocate (ra(nbin))
  do ip=1,npart
    allocate (dust%particles(ip)%p(nbin,3))
    allocate (dust%particles(ip)%v(nbin,3))
    allocate (dust%particles(ip)%m(nbin))
    allocate (dust%particles(ip)%s(nbin))
    allocate (dust%particles(ip)%t(nbin))
    allocate (dust%particles(ip)%dt(nbin))
    dust%particles(ip)%p(:,:)  = 0.0
    dust%particles(ip)%v  = 0.0
    dust%particles(ip)%m  = 1.0
    dust%particles(ip)%s  = 1.0
    dust%particles(ip)%t  = 0.0
    dust%particles(ip)%dt = 1.0
  end do
  !-----------------------------------------------------------------------------
  ! Initialize positions
  !-----------------------------------------------------------------------------
  star%position = 1d0

  !-----------------------------------------------------------------------------
  ! Evolve with explicit variables
  !-----------------------------------------------------------------------------
  call tic
  do ip=1,npart
    p  => dust%particles(ip)%p
    v  => dust%particles(ip)%v
    m  => dust%particles(ip)%m
    s  => dust%particles(ip)%s
    t  => dust%particles(ip)%t
    dt => dust%particles(ip)%dt
    ra = 0.0
    do id=1,ndim
    do ib=1,nbin
      r(ib,id) = p(ib,id) - star%position(id)
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
  end do
  call toc ("particle motion", npart*nbin)

  !-----------------------------------------------------------------------------
  ! Evolve with subroutine call exposing arrays
  !-----------------------------------------------------------------------------
  call tic
  do ip=1,npart
    p  => dust%particles(ip)%p
    v  => dust%particles(ip)%v
    m  => dust%particles(ip)%m
    s  => dust%particles(ip)%s
    t  => dust%particles(ip)%t
    dt => dust%particles(ip)%dt
    call move_particles (p, v, m, s, t, dt, ra, r, f, a, star%position, nbin, ndim)
  end do
  call toc ("particle motion", npart*nbin)

  deallocate (dust)
END PROGRAM

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
