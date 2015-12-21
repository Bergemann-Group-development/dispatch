!> =============================================================================
!> $Id$
!> Implements alternative ways to compute the Planck source funtion.
!> =============================================================================
MODULE planck_mod
  USE io_mod
  USE units_mod
  USE scaling_mod
  implicit none
  private
  type, public:: planck_t
    integer   :: n
    logical   :: debug
    real      :: cx
    real      :: xmin
    real      :: xmax
    real, dimension(:), pointer:: table
    procedure(function), pointer :: source => null()
  contains
    procedure :: init
    procedure :: lookup
    procedure :: function
    procedure :: stefan
    procedure :: test
  end type
  type(planck_t), public:: planck
CONTAINS

!> =============================================================================
!> Initialize table lookup of the Planck function.  Express everything in the
!> non-dimensional variable x = h*c/(lambda*k*T), which makes the Planck function
!> = c*T**4/((exp(x)-1.0)*x**4)
!> =============================================================================
SUBROUTINE init (self)
  class (planck_t):: self
  integer:: n=10000
  real(8):: xmin=1d-2, xmax=1d1
  character(len=12) :: source='function'        ! choice of source function
  integer:: i
  real:: x, lnx
  logical:: debug=.false.
  namelist /planck_params/ n, xmin, xmax, debug, source
  !.............................................................................
  rewind (io%input)
  read (io%input, planck_params)
  if (io%master) write (*, planck_params)
  select case(source)
  case('function')
    self%source => function
  case('lookup')
    self%source => lookup
  case('stefan')
    self%source => stefan
  end select
  self%n = n
  self%xmin = xmin
  self%xmax = xmax
  self%debug = debug
  self%cx = real(self%n,kind=8)/(xmax-xmin)
  allocate (self%table(0:self%n))
  if (self%debug) print *, 'test table values'
  do i=0,self%n
    x = self%xmin + i/self%cx
    self%table(i) = x**4/(exp(x)-1.0)
    if (self%debug) print *, i, x, self%table(i)
  end do
  call self%test
END SUBROUTINE init

!> =============================================================================
!> The planck function, normalized to give radiation per unit \ln\lambda,  for
!> a given temperature and wavelength. Input and output are in code units, but
!> we use lambda in micrometers, to avoid underflow and divide fault in the
!> expression when computing lambda**4.
!> =============================================================================
SUBROUTINE function (self, n, lambda, temperature, source)
  class(planck_t):: self
  integer:: n(2)
  real:: lambda, lambda_mum, lambda_4
  real, dimension(n(1),n(2)):: temperature, source
  integer :: i
  real(8) :: c1, c2
  real, dimension(n(1)):: e
  !.............................................................................
  c1 = 2d0*cgs%h_p*cgs%c**2/(cgs%mum**4)                ! 2*h*c**2/lambda**4
  c2 = cgs%h_p*cgs%c/(cgs%k_b*cgs%mum*scaling%temp)     ! h*c/(lambda*T)
  lambda_mum = lambda*scaling%l/cgs%mum                 ! lambda in micrometers
  lambda_4 = lambda_mum**4
  do i=1,n(2)
    e(:) = exp(-c2/(temperature(:,i)*lambda_mum))
    !e(:) = c2/(temperature(:,i)*lambda_mum)
    source(:,i) = c1*e(:)/((1.0-e(:))*lambda_4)
  end do
END SUBROUTINE function

!> =============================================================================
!> Linear interpolation in linear wavelength, for maximum speed.  To handle 
!> cases where some x are below xmin one could implement a test on the smallest
!> x, leading to a slower version, with an assumed x**3 behavior below xmin.
!> =============================================================================
SUBROUTINE lookup (self, n, wavelength, temperature, source)
  class (planck_t):: self
  integer                       :: n(2)
  real                          :: wavelength
  real,    dimension(n(1),n(2)) :: temperature
  real,    dimension(n(1),n(2)) :: source
  !.............................................................................
  real                          :: c1, c2, src
  integer                       :: i, j
  real,    dimension(n(1))      :: x
  integer, dimension(n(1))      :: k
  !-----------------------------------------------------------------------------
  c1 = cgs%h_p*cgs%c/(wavelength*scaling%l*cgs%k_b*scaling%temp)
  c2 = 2.*(cgs%k_b*scaling%temp)**4/(cgs%h_p**3*cgs%c**2)
  !-----------------------------------------------------------------------------
  ! Case with printout
  !-----------------------------------------------------------------------------
  if (self%debug) then
    print*,'table lookup', c1, c2, wavelength, temperature(1,1)
    call self%function (n, wavelength, temperature, source)
    do i=1,n(1)
    do j=1,n(2)
      x(i) = (c1/temperature(i,j) - self%xmin)*self%cx
      k(i)  = max (0, min (self%n-1, int(x(i))))
      x(i)  = x(i)-k(i)
      src   = source(i,j)
      source(i,j) = c2*temperature(i,j)**4*( &
              x(i) *self%table(k(i)+1) + &
         (1.0-x(i))*self%table(k(i)  ) )
      print *,i, j, x(i), k(i), src, source(i,j)
    end do
    end do
  !-----------------------------------------------------------------------------
  ! Optimized case
  !-----------------------------------------------------------------------------
  else
    do j=1,n(2)
      x = (c1/temperature(:,j) - self%xmin)*self%cx
      k  = max (0, min (self%n-1, int(x)))
      x  = x-k
      source(:,j) = c2*temperature(:,j)**4*(x*self%table(k+1) + (1.0-x)*self%table(k))
    end do
  end if
END SUBROUTINE lookup

!> =============================================================================
!> The integrated Planck function is proportional to Stefan's constant times the
!> fourth power of the temperature.  This is the relevant source function when
!> the opacity is constant as a function of wavelength, so one can integrate 
!> over wavelengthe before solving the radiative transfer equation
!> =============================================================================
SUBROUTINE stefan (self, n, lambda, temperature, source)
  class(planck_t):: self
  integer:: n(2)
  real, dimension(n(1),n(2)):: temperature, source
  real:: lambda, c
!...............................................................................
  c = cgs%stefan/cgs%pi*scaling%temp**4/(scaling%e/(scaling%l**2*scaling%t))
  source = c*temperature**4
END SUBROUTINE stefan

!> =============================================================================
!> Test the table lookup.  Wavelength is in micrometers, and the temperatures
!> are the ones the span over the x-range in the table.
!> =============================================================================
SUBROUTINE test (self)
  class(planck_t):: self
  integer:: m, n(2), i, j
  real:: wavelength, x
  real, dimension(:,:), pointer:: temperature, source1, source2
  !.............................................................................
  m = self%n*3
  allocate (temperature(m,1), source1(m,1), source2(m,1))
  wavelength = 100.*cgs%mum/scaling%l
  if (self%debug) print *, 'test temperatures'
  do i=1,m
    x = self%xmin + (i-1)*((self%xmax-self%xmin)/real(m))
    temperature(i,1) = cgs%h_p*cgs%c/(cgs%k_b*x*wavelength*scaling%l)
    if (self%debug) print *, i, x, temperature(i,1)
    temperature(i,1) = temperature(i,1)/scaling%temp
  end do
  n = [m,1]
  call self%lookup (n, wavelength, temperature, source1)
  call self%function (n, wavelength, temperature, source2)
  source1 = source1/source2-1.0
  print*,'planck_mod: min and max deviations =', minval(source1), maxval(source1)
  deallocate (temperature, source1, source2)
END SUBROUTINE test

!> =============================================================================
END MODULE planck_mod
