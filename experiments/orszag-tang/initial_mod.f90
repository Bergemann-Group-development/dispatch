!> =============================================================================
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  implicit none
  private
  type, public :: initial_t
    integer :: idirect
    real(8) :: p0, v0, b0
    procedure(initial_conditions), pointer :: init_mhd_vars => null()
  contains
    procedure :: init
  end type

  real(8), parameter:: pi = asin(1.0d0) * 2.0d0
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(initial_t) :: self
  integer, save :: idirect
  real(8), save :: p0, v0, b0
  logical, save :: first_time=.true.
  namelist /initial_params/ idirect, b0, p0, v0
  !.............................................................................
  call trace_begin ('initial_t%init')

  ! default values for namelist parameters
  idirect = 3
  b0 = 0.5 / sqrt(pi)
  p0 = 5.0 / (12.0 * pi)
  v0 = 1.0

  !$omp critical (initial_t)
  if (first_time) then
    rewind (io%input); read (io%input,initial_params)
    if (io%master) write(*,initial_params)
    first_time = .false.
  end if
  !$omp end critical (initial_t)

  self%idirect = idirect
  self%p0 = p0
  self%v0 = v0
  self%b0 = b0

  self%init_mhd_vars => initial_conditions

  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Initial conditions, callable also from refine_mod.f90
!> =============================================================================
SUBROUTINE initial_conditions (self, m, gamma, ff)
  class(initial_t) :: self
  class(mesh_t), allocatable, dimension(:) :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,8) :: ff
  real, intent(in) :: gamma
  !.............................................................................
  integer:: ix, iy, iz                                       ! loop indices
  real:: entropy, twopi, xa, ya, za, xb, yb, zb, d0          ! local variables
  real,dimension(:,:,:),allocatable :: vecpot
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  !.............................................................................
  call trace_begin ('initial_t%init_mhd_vars')

  allocate(vecpot(m(1)%gn,m(2)%gn,m(3)%gn))

  twopi = 2.0 * pi
  d0 = gamma * self%p0

  ! "pre"-initialisation of all variables, independent of `idirect`.
  if (gamma==1.0) then
    entropy = 1.0
  else
    entropy = (log(self%p0)-log(d0)*gamma)/(gamma-1.0)
  end if
  do iz=m(3)%lb,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ff(ix,iy,iz,1) = d0
        ff(ix,iy,iz,2) = entropy*d0
        ff(ix,iy,iz,3) = 0.0
        ff(ix,iy,iz,4) = 0.0
        ff(ix,iy,iz,5) = 0.0
        ff(ix,iy,iz,6) = 0.0
        ff(ix,iy,iz,7) = 0.0
        ff(ix,iy,iz,8) = 0.0
      end do
    end do
  end do

  ! initialisation based on value of `idirect`
  ! ---- x1-x2 plane ----
  if (self%idirect == 3) then

    ! initialise vector potential
    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        ya = m(2)%p + r2(iy) - 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          xa = m(1)%p + r1(ix) - 0.5*m(1)%d
          vecpot(ix,iy,iz) = (self%b0 / twopi) * ( 0.5 * cos(2.0 * twopi * xa) + cos(twopi  * ya) )
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub-1
        yb = m(2)%p + r2(iy)
        print*,iy,yb
        do ix=m(1)%lb,m(1)%ub-1
          xb = m(1)%p + r1(ix)
          ff(ix,iy,iz,3) = -d0 * self%v0 * sin(twopi * yb)
          ff(ix,iy,iz,4) =  d0 * self%v0 * sin(twopi * xb)
          ff(ix,iy,iz,6) =  (vecpot(ix,iy+1,iz) - vecpot(ix,iy  ,iz) ) / m(2)%d
          ff(ix,iy,iz,7) = -(vecpot(ix+1,iy,iz) - vecpot(ix  ,iy,iz) ) / m(1)%d
        end do
      end do
    end do
    print*,'mk1', m(1:2)%p, minval(ff(:,:,:,3)), maxval(ff(:,:,:,3))
  end if

  ! ---- x2-x3 plane ----
  if (self%idirect == 1) then

    ! initialise vector potential
    do iz=m(3)%lb,m(3)%ub
      za = m(3)%p + r3(iz) - 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
        ya = m(2)%p + r2(iy) - 0.5*m(2)%d
        do ix=m(1)%lb,m(1)%ub
          vecpot(ix,iy,iz) = self%b0 / twopi * ( 0.5d0 * cos (2.0 * twopi * ya) + cos (twopi  * za) )
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub-1
      zb = m(3)%p + r3(iz)
      do iy=m(2)%lb,m(2)%ub-1
        yb = m(2)%p + r2(iy)
        do ix=m(1)%lb,m(1)%ub
          ff(ix,iy,iz,4) = -d0 * self%v0 * sin(twopi * zb)
          ff(ix,iy,iz,5) =  d0 * self%v0 * sin(twopi * yb)
          ff(ix,iy,iz,7) =  (vecpot(ix,iy,iz+1) - vecpot(ix,iy,iz  ) ) / m(3)%d
          ff(ix,iy,iz,8) = -(vecpot(ix,iy+1,iz) - vecpot(ix,iy  ,iz) ) / m(2)%d
        end do
      end do
    end do
  end if

  ! ---- x3-x1 plane ----
  if (self%idirect == 2) then

    ! initialise vector potential
    do iz=m(3)%lb,m(3)%ub
      za = m(3)%p + r3(iz) - 0.5*m(3)%d
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          xa = m(1)%p + r1(iy) - 0.5*m(1)%d
          vecpot(ix,iy,iz) = self%b0 / twopi * ( 0.5d0 * cos (2.0 * twopi * za) + cos (twopi  * xa) )
        end do
      end do
    end do

    do iz=m(3)%lb,m(3)%ub-1
      zb = m(3)%p + r3(iz)
      do iy=m(2)%lb,m(2)%ub
        do ix=m(1)%lb,m(1)%ub-1
          xb = m(1)%p + r1(ix)
          ff(ix,iy,iz,5) = -d0 * self%v0 * sin(twopi * xb)
          ff(ix,iy,iz,3) =  d0 * self%v0 * sin(twopi * zb)
          ff(ix,iy,iz,8) =  (vecpot(ix+1,iy,iz) - vecpot(ix+1,iy,iz) ) / m(1)%d
          ff(ix,iy,iz,6) = -(vecpot(ix,iy,iz+1) - vecpot(ix,iy,iz+1) ) / m(3)%d
        end do
      end do
    end do
  end if

  deallocate(vecpot)

  call trace_end
  end associate
END SUBROUTINE initial_conditions

END MODULE initial_mod
