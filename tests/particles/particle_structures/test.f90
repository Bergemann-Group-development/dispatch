!===============================================================================
!> Program to test functionality and speed of particles, when using an array
!> of structures
!===============================================================================

PROGRAM test
  USE list_mod
  USE task_mod
  USE mpi_timer
  USE particle_mod
  USE planet_mod
  USE star_mod
  implicit none
  type part_t
    real:: p(3), v(3)
    real:: m, s, t, dt
  end type
  class(star_t), pointer:: star
  class(planet_t), pointer:: planet
  class(part_t), dimension(:), pointer:: dust
  integer:: ndim=3, nbin=64, npart=10000000, ip, ib, id
  real, dimension(:), pointer:: p, v
  real, dimension(3):: a, f, r
  real, pointer:: m, s, t, dt
  real:: ra

  !-----------------------------------------------------------------------------
  ! Allocate variables
  !-----------------------------------------------------------------------------
  allocate (star, planet)
  allocate (dust(npart))
  do ip=1,npart
    dust(ip)%p  = 0.0
    dust(ip)%v  = 0.0
    dust(ip)%m  = 1.0
    dust(ip)%s  = 1.0
    dust(ip)%t  = 0.0
    dust(ip)%dt = 1.0
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
    p  => dust(ip)%p
    v  => dust(ip)%v
    m  => dust(ip)%m
    s  => dust(ip)%s
    t  => dust(ip)%t
    dt => dust(ip)%dt
    ra = 0.0
    do id=1,3
      r(id) = p(id) - star%position(id)
      ra = ra + r(id)**2
    end do
    ra = sqrt(ra)**3
    do id=1,ndim
      f(id) = -r(id)/ra
      a(id) =  f(id)/m
    end do
    do id=1,ndim
      v(id) = v(id) + a(id)*dt
      p(id) = p(id) + v(id)*dt
    end do
  end do
  call timer%toc ("particle motion", npart)

  !-----------------------------------------------------------------------------
  ! Evolve with explicit variables
  !-----------------------------------------------------------------------------
  call tic
  do ip=1,npart
    ra = 0.0
    do id=1,3
      r(id) = dust(ip)%p(id) - star%position(id)
      ra = ra + r(id)**2
    end do
    ra = sqrt(ra)**3
    do id=1,3
      f(id) = -r(id)/ra
      a(id) =  f(id)/dust(ip)%m
    end do
    do id=1,3
      dust(ip)%v(id) = dust(ip)%v(id) + a(id)*dust(ip)%dt
      dust(ip)%p(id) = dust(ip)%p(id) + dust(ip)%v(id)*dust(ip)%dt
    end do
  end do
  call timer%toc ("particle motion", npart)

  deallocate (dust)
END PROGRAM
