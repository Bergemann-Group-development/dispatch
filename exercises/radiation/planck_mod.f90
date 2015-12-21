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
!> = c*T**4*x**4/(exp(x)-1).  Interpolate in x/(exp(x)-1), which is smooth. The
!> limits xmin=0.001 and xmax=20 are perfect, since the error below 0.01 is less
!> 0.05%, and the value above x=20 is negligible.
!> =============================================================================
SUBROUTINE init (self)
  class (planck_t):: self
  integer:: n=1000
  real(8):: xmin=1d-3, xmax=2d1
  character(len=12) :: source='lookup'        ! choice of source function
  integer:: i
  real:: x, lnx, dx
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
  dx = (xmax-xmin)/real(self%n,kind=8)
  allocate (self%table(0:self%n))
  if (self%debug) print *, 'test table values'
  do i=0,self%n
    x = self%xmin + i*dx
    self%table(i) = x/(exp(x)-1.0)
    if (self%debug) print *, i, x, self%table(i)
  end do
  call self%test
  self%debug = .false.
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
    source(:,i) = c1*e(:)/((1.0-e(:))*lambda_4)
  end do
END SUBROUTINE function

!> =============================================================================
!> Linear interpolation in linear wavelength, for maximum speed.  To handle 
!> uncommon cases where some x are below xmin a test is made on x, before 
!> choosing how to interpolate.  This is only slightly slower than if one knows
!> everything is in range.
!> =============================================================================
SUBROUTINE lookup (self, n, wavelength, temperature, source)
  class (planck_t):: self
  integer                       :: n(2)
  real                          :: wavelength
  real,    dimension(n(1),n(2)) :: temperature
  real,    dimension(n(1),n(2)) :: source
  !.............................................................................
  real                          :: c1, c2, src, x1, f
  integer                       :: i, j
  real,    dimension(n(1))      :: p
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
    print*,self%table(0:1)
    do i=1,n(1)
    do j=1,n(2)
      x1 = c1/temperature(i,j)
      src   = source(i,j)
      if (x1 < self%xmin) then
        k(i) = 0
        source(i,j) = c2*temperature(i,j)**4*x1**3*self%table(0)
        f = source(i,j)/(c2*temperature(i,j)**4)
      else
        p(i) = (x1 - self%xmin)*self%cx
        k(i)  = min (self%n-1, int(p(i)))
        p(i)  = p(i)-k(i)
        f =   p(i) *self%table(k(i)+1) + &
         (1.0-p(i))*self%table(k(i)  )
        source(i,j) = c2*temperature(i,j)**4*f*x1**3
      end if
      print *,i, x1, p(i), k(i), src/(c2*temperature(i,j)**4), f
    end do
    end do
  !-----------------------------------------------------------------------------
  ! Optimized case
  !-----------------------------------------------------------------------------
  else
    do j=1,n(2)
      p = c1/temperature(:,j)
      if (minval(p) >= self%xmin) then
        source(:,j) = c2*temperature(:,j)**4*p**3
        p = (p-self%xmin)*self%cx
        k = min (self%n-1, int(p))
        p = p-k
        source(:,j) = source(:,j)* &
                      (p*self%table(k+1) + (1.0-p)*self%table(k))
      else
        where (p < self%xmin)
          source(:,j) = c2*temperature(:,j)**4*p**3*self%table(0)
        else where
          source(:,j) = c2*temperature(:,j)**4*p**3
          p = (p-self%xmin)*self%cx
          k  = min (self%n-1, int(p))
          p  = p-k
          source(:,j) = source(:,j)* &
                        (p*self%table(k+1) + (1.0-p)*self%table(k))
        end where
      end if
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
!> Test the table lookup.  Wavelength (in micrometers) is unimportant. The
!> temperatures spans over the x-range in the table + an order of magnitude below
!> =============================================================================
SUBROUTINE test (self)
  class(planck_t):: self
  integer:: m, n(2), i, j
  real:: wavelength, x, xmin
  real, dimension(:,:), pointer:: temperature, source1, source2
  !.............................................................................
  m = self%n*4
  allocate (temperature(m,1), source1(m,1), source2(m,1))
  wavelength = 100.*cgs%mum/scaling%l
  if (self%debug) print *, 'test temperatures'
  do i=1,m
    xmin = self%xmin*0.1
    x = xmin + (i-1)*((self%xmax-xmin)/real(m-1))
    temperature(i,1) = cgs%h_p*cgs%c/(cgs%k_b*x*wavelength*scaling%l)
    if (self%debug) print *, i, x, temperature(i,1)
    temperature(i,1) = temperature(i,1)/scaling%temp
  end do
  n = [m,1]
  call self%source (n, wavelength, temperature, source1)
  call self%function (n, wavelength, temperature, source2)
  source1 = source1/source2-1.0
  print*,'planck_mod: min and max deviations =', minval(source1), maxval(source1)
  deallocate (temperature, source1, source2)
END SUBROUTINE test

!> =============================================================================
END MODULE planck_mod
