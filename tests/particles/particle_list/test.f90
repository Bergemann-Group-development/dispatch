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
END MODULE dust_mod

!===============================================================================
PROGRAM test
  USE list_mod
  USE task_mod
  USE mpi_timer
  USE dust_mod
  USE star_mod
  USE planet_mod
  implicit none
  class(star_t), pointer:: star
  class(planet_t), pointer:: planet
  class(part_t), pointer:: head, dust, prev => null(), next
  integer:: ndim=3, nbin=64, npart=100000, ip, ib, id
  real, dimension(:), pointer:: p, v
  real, dimension(3):: a, f, r
  real, pointer:: m, s, t, dt
  real:: ra

  !-----------------------------------------------------------------------------
  ! Allocate variables
  !-----------------------------------------------------------------------------
  allocate (star, planet)
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
  ! Evolve with explicit variables
  !-----------------------------------------------------------------------------
  call tic
  dust => head
  do while(associated(dust))
    p  => dust%p
    v  => dust%v
    m  => dust%m
    s  => dust%s
    t  => dust%t
    dt => dust%dt
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
    dust => dust%next
  end do
  call timer%toc ("particle motion", npart*nbin)

  !-----------------------------------------------------------------------------
  ! Evolve with explicit variables
  !-----------------------------------------------------------------------------
  call tic
  dust => head
  do while(associated(dust))
    ra = 0.0
    do id=1,3
      r(id) = dust%p(id) - star%position(id)
      ra = ra + r(id)**2
    end do
    ra = sqrt(ra)**3
    do id=1,3
      f(id) = -r(id)/ra
      a(id) =  f(id)/dust%m
    end do
    do id=1,3
      dust%v(id) = dust%v(id) + a(id)*dust%dt
      dust%p(id) = dust%p(id) + dust%v(id)*dust%dt
    end do
    dust => dust%next
  end do
  call timer%toc ("particle motion", npart*nbin)

  !-----------------------------------------------------------------------------
  ! Evolve with explicit variables
  !-----------------------------------------------------------------------------
  call tic
  dust => head
  do while(associated(dust))
    ra = 0.0
    r = dust%p - star%position
    ra = sum(r**2)
    ra = sqrt(ra)**3
    f = -r/ra
    a =  f/dust%m
    dust%v = dust%v + a*dust%dt
    dust%p = dust%p + dust%v*dust%dt
    dust => dust%next
  end do
  call timer%toc ("particle motion", npart*nbin)

  !-----------------------------------------------------------------------------
  ! Evolve with explicit variables
  !-----------------------------------------------------------------------------
  call tic
  dust => head
  do while(associated(dust))
    ra = 0.0
    r = dust%p - star%position
    ra = sum(r**2)
    ra = sqrt(ra)**3
    f = -r/ra
    a =  f/dust%m
    dust%v = dust%v + a*dust%dt
    dust%p = dust%p + dust%v*dust%dt
    dust => dust%next
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
