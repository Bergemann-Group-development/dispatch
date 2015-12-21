!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Radiative Transfer problem in 3-D
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pde_mod
  USE io_mod
  USE omp_timer_mod
  USE trace_mod
  USE patches_mod
  USE scalar_mod
  USE scaling_mod
  USE microphysics_mod
  USE rt_mod
  implicit none
  private
  integer, parameter:: i_d=1, i_s=2, i_tt=3
  type, public, extends(patch_t):: pde_t
    type(scalar_t):: d, s, dsdt, tt
    real:: w_gauss
    real:: p_gauss
    real:: a_gauss
    real:: amplitude
    real:: gamma
    real(kind=8):: tolerance
    type(rt_t)  :: rt
    type(microphysics_t) :: microphysics
  contains
    procedure:: init
    procedure:: solve
  end type
  character(len=12):: solver='integral'
CONTAINS

!> =============================================================================
!> Setup a patch with density and source function for the radiative transfer solver
!> =============================================================================
SUBROUTINE init (self)
  class(pde_t):: self
  integer:: i, j, k
  real:: x, y, z, r, density
  logical:: first_time = .true.
  real:: w_gauss=0.15, p_gauss=4.0, a_gauss=0.10
  real:: amplitude=25., temperature=100., gamma=1.4
  namelist /rt_params/ w_gauss, p_gauss, a_gauss, amplitude, temperature, gamma, solver
associate(gn=>self%gn, n=>self%n, l=>self%li, u=>self%ui)
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  experimentname = 'rt'
  call io%init (experimentname)
  call trace_begin ('hd_mod:pde_t%init', 2)
  if (first_time) then
    call scaling%init
    first_time = .false.
    rewind (io%input)
    read (io%input, rt_params)
    if (io%master) write (*,rt_params)
  end if
  self%size      = 1.
  self%position  = 0.
  self%w_gauss   = w_gauss
  self%p_gauss   = p_gauss
  self%a_gauss   = a_gauss
  self%gamma     = gamma
  self%amplitude = amplitude
  self%nv        = 3
  call self%patch_t%init
  self%nv        = 3
  !-----------------------------------------------------------------------------
  ! Setup microphysics for RT. One microphysics object per wavelength is needed
  !-----------------------------------------------------------------------------
  call self%microphysics%init
  !-----------------------------------------------------------------------------
  ! Memory mappinf og density, temperature, entropy, and cooling
  !-----------------------------------------------------------------------------
  call   self%d%init (gn, self%size/n, 'd')
  call   self%s%init (gn, self%size/n, 's')
  call  self%tt%init (gn, self%size/n, 'tt')
  self%d%f    => self%mem(:,:,:,i_d ,self%it,1)
  self%tt%f   => self%mem(:,:,:,i_tt,self%it,1)
  self%s%f    => self%mem(:,:,:,i_s ,self%it,1)
  self%dsdt%f => self%mem(:,:,:,i_s ,self%it,2)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Gaussian top hat with 1/r**2 wings in density,
  ! with constant temperature defining the entropy
  !-----------------------------------------------------------------------------
  if (self%time==0.0) then
    do k=1,gn(3)
    do j=1,gn(2)
    do i=1,gn(1)
      x = (i-(gn(1)+1)*0.5)*self%ds(1) + self%position(1)
      y = (j-(gn(2)+1)*0.5)*self%ds(2) + self%position(2)
      z = (k-(gn(3)+1)*0.5)*self%ds(3) + self%position(3)
      r = sqrt(x**2+y**2+z**2)/self%w_gauss
      r = min(r,0.5/w_gauss)
     density = &
        self%amplitude*exp(-r**self%p_gauss/(1.+r**self%p_gauss/log(r**2/self%a_gauss)))
      self%d%f(i,j,k)  = density
      self%tt%f(i,j,k) = temperature/scaling%temp
      self%s%f(i,j,k)  = log(self%tt%f(i,j,k))/(self%gamma-1.0) - log(density)
    end do
    end do
    end do
  end if
  print*,'mass:',sum(self%d%f(l(1):u(1),l(2):u(2),l(3):u(3)))*product(self%ds)
  !-----------------------------------------------------------------------------
  ! Set up the right dimensions from current patch, etc for radiative transfer
  !-----------------------------------------------------------------------------
  call self%rt%init (solver)
  !
  call trace_end
end associate
END SUBROUTINE init

!> =============================================================================
!> Compute the solution to the radiative transfer equation
!> =============================================================================
SUBROUTINE solve (self)
  class(pde_t):: self
  real, pointer, dimension(:,:) :: dtau, src, source, qcool, opacity, &
                                   density, temperature
  real, pointer, dimension(:)   :: ifw0, irv0, ifw1, irv1
  integer                       :: ix, iy, iz
  real(8)                       :: wc(6), cost(5)
  !.............................................................................
  call trace_begin ('radiation_mod:pde_t%update', 2)
associate(gn=>self%gn, n=>self%n, ntau=>self%n(3), l=>self%li, u=>self%ui)
  allocate (ifw0(n(1)), irv0(n(1)), ifw1(n(1)), irv1(n(1)))
  allocate (density    (l(1):u(1), l(3):u(3)))
  allocate (temperature(l(1):u(1), l(3):u(3)))
  allocate (opacity    (l(1):u(1), l(3):u(3)))
  allocate (dtau       (l(1):u(1), l(3):u(3)))
  allocate (source     (l(1):u(1), l(3):u(3)))
  allocate (qcool      (l(1):u(1), l(3):u(3)))
  !.............................................................................
  ifw0 = 0.0; irv1 = 0.0
  cost = 0d0
  !-----------------------------------------------------------------------------
  ! Loop over xz-slices, while accumulating the costs in wallclock seconds, and
  ! computing the time derivative dSdt = Q/T, where the heating per unit volume
  ! is Q = integral d\Omega opacity*(intensiy-source), and where we approximate
  ! the integral over space angle \Omega as just a product with 4\pi.
  !-----------------------------------------------------------------------------
  do iy=l(2),u(2)
      wc(1) = wallclock()
    density     = self%d%f (l(1):u(1), iy, l(3):u(3))
    temperature = self%tt%f(l(1):u(1), iy, l(3):u(3))
      wc(2) = wallclock()
    opacity = self%microphysics%opacity (density, temperature)
    dtau    = self%microphysics%dtau (self%ds(3), opacity)
      wc(3) = wallclock()
    source = self%microphysics%source (density, temperature)
      wc(4) = wallclock()
    call self%rt%solver (n(1), n(3), dtau, source, ifw0, ifw1, qcool, irv0, irv1)
      wc(5) = wallclock()
    self%dsdt%f(l(1):u(1),iy,l(3):u(3)) = 4.*scaling%pi*opacity*qcool/temperature
      wc(6) = wallclock()
      cost = cost + wc(2:6) - wc(1:5)
  end do
  !-----------------------------------------------------------------------------
  ! Convert the cost to nanoseconds per point
  !-----------------------------------------------------------------------------
  cost = cost*1e9/product(n)
  print'(a,6f6.1))', ' mus/pt: copy, dtau, source, solve, copy, total =', cost, sum(cost)
  print*, 'minmax(dSdt):', minval(self%dsdt%f(l(1):u(1),l(2):u(2),l(3):u(3))), &
                           maxval(self%dsdt%f(l(1):u(1),l(2):u(2),l(3):u(3)))
  deallocate (dtau, source, qcool, ifw0, irv0, ifw1, irv1)
  !-----------------------------------------------------------------------------
  ! Just to get dSdt included in the output file, use the temperature slot
  !-----------------------------------------------------------------------------
  self%tt%f = self%dsdt%f
end associate
  call trace_end
END SUBROUTINE solve
!> =============================================================================
END MODULE pde_mod
