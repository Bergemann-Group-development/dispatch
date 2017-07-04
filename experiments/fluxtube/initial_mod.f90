!> =============================================================================
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  implicit none
  private
  type, public :: initial_t
    real :: x0, y0, z0, r0, b0=0.1, temp=1.0, gamma=5./3., beta
    integer :: expe
    logical :: mhd
    real(8) :: time
    logical :: pic_needed
    logical :: first_pic
    procedure(initial_conditions), pointer :: init_mhd_vars => null()
  contains
    procedure :: init
    procedure, nopass :: cleanup
  end type
  real, save, pointer, dimension(:,:,:,:):: buffer => null()
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(initial_t) :: self
  real, save :: x0, y0, z0, r0, B0, temp, gamma, beta
  integer, save :: expe
  logical, save :: first_time=.true.
  logical, save :: mhd=.true.
  logical, save :: pic_needed = .false.
  logical, save :: first_pic = .true.
  namelist /experiment_params/ x0, y0, z0, r0, B0, temp, gamma, beta, expe, first_pic
  !.............................................................................
  call trace_begin ('initial_t%init')

  !$omp critical (initial_t)
  if (first_time) then
    rewind (io%input); read (io%input,experiment_params)
    if (io%master) write(*,experiment_params)
    first_time = .false.
  end if
  !$omp end critical (initial_t)

  self%x0 = x0
  self%y0 = y0
  self%z0 = z0
  self%r0 = r0
  self%B0 = B0
  self%temp = temp
  self%gamma = gamma
  self%mhd = mhd
  self%time = 0d0
  self%beta = beta
  self%expe = expe
  self%pic_needed = pic_needed
  self%first_pic = first_pic

  self%init_mhd_vars => initial_conditions

  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE cleanup
  if (associated(buffer)) deallocate(buffer)
END SUBROUTINE cleanup

!> =============================================================================
!> Initial conditions, callable also from refine_mod.f90
!> =============================================================================
SUBROUTINE initial_conditions (self, m, ff)
  class(initial_t) :: self
  class(mesh_t), pointer :: m(:)
  real, dimension(:,:,:,:), pointer :: ff
  !class(_t), pointer:: patch
  !.............................................................................
  integer:: ix, iy, iz                                       ! loop indices
  real:: x, y, z, dx, dy, dz
  real:: mu0 = 1.0, pg !, r0 = 10., B0 = 5.0
  real:: jx, jy, Flx, Fly, profile, g1                       ! local variables
  real, dimension(:,:,:), pointer:: d, e, px, py, pz, Bx, By, Bz, pp
  !.............................................................................
  call trace_begin ('initial_t%init_mhd_vars')

  dx = 0.5*m(1)%d
  dy = 0.5*m(2)%d
  dz = 0.5*m(3)%d
  d  => ff(:,:,:,1)
  e  => ff(:,:,:,2)
  px => ff(:,:,:,3)
  py => ff(:,:,:,4)
  pz => ff(:,:,:,5)
  bx => ff(:,:,:,6)
  by => ff(:,:,:,7)
  bz => ff(:,:,:,8)

  g1 = self%gamma-1.0

  select case(self%expe)

  ! The static equilibrium fluxtube experiment
  case(1)

    do iz=m(3)%lb,m(3)%ub
      do iy=m(2)%lb,m(2)%ub
        y = m(2)%p + m(2)%r(iy)
        do ix=m(1)%lb,m(1)%ub
          x = m(1)%p + m(1)%r(ix)
          profile = exp(-2*((x-self%x0)**2+(y-self%y0)**2)/self%r0**2)     ! Gaussian
          e(ix,iy,iz)  = self%B0**2/(2*mu0)*(1.+self%beta-profile)/g1      ! internal energy
          d(ix,iy,iz)  = e(ix,iy,iz)*g1/self%temp                          ! density
          px(ix,iy,iz) = 0.0                                               ! momentum
          py(ix,iy,iz) = 0.0
          pz(ix,iy,iz) = 0.0
          Bx(ix,iy,iz) = 0.0                                               ! Magnetic field
          By(ix,iy,iz) = 0.0
          Bz(ix,iy,iz) = self%B0 * exp(-((x+dx-self%x0)**2+(y+dx-self%y0)**2)/self%r0**2)
        end do
      end do
    end do

  case(2)
  ! Fluxtube experiment with a pertubation in the y-direction

    do iz=m(3)%lb,m(3)%ub
      z = m(3)%p + m(3)%r(iz)
      do iy=m(2)%lb,m(2)%ub
        y = m(2)%p + m(2)%r(iy)
        do ix=m(1)%lb,m(1)%ub
          x = m(1)%p + m(1)%r(ix)
          profile = exp(-2*((x-self%x0)**2+(y-self%y0)**2)/self%r0**2)     ! Gaussian
          e(ix,iy,iz)  = self%B0**2/(2*mu0)*(1.+self%beta-profile)/g1      ! Gas pressure
          d(ix,iy,iz)  = e(ix,iy,iz)*g1/self%temp                          ! density
          px(ix,iy,iz) = 0.0                                               ! momentum
          py(ix,iy,iz) = 0.0
          pz(ix,iy,iz) = 0.0
          Bx(ix,iy,iz) = 0.0                                               ! Magnetic field
          By(ix,iy,iz) = exp(-((x+dx-self%x0)**2+(z+dz-self%z0)**2)/self%r0**2)
          Bz(ix,iy,iz) = self%B0 * exp(-((x+dx-self%x0)**2+(y+dy-self%y0)**2)/self%r0**2)
        end do
      end do
    end do

  end select
  print *,'Bz', sum(bz**2), -3

  if (io%debug(2)) then
    iz = (m(3)%lb+m(3)%ub)/2
    do ix=m(1)%li,m(1)%ui
      print *,'initial B(x):',ix,iz,bz(ix,iz,iz)
    end do
    do iy=m(1)%li,m(1)%ui
      print *,'initial B(y):',iz,iy,bz(iz,iy,iz)
    end do
    print *,'minval,maxval(d)',minval(d),maxval(d)
  end if

  call trace_end
END SUBROUTINE initial_conditions

END MODULE initial_mod
