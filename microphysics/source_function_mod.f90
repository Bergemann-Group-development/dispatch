MODULE source_function_mod
  USE scaling_mod
  USE units_mod
  USE io_mod
  USE trace_mod
  USE eos_mod
  implicit none
  private
  !  -------------------------------------------------------------------
  type :: planck_lookup_t
    real    :: xmin, xmax, cx
    integer :: n
    real, dimension(:), pointer:: table
  end type
  !  -------------------------------------------------------------------
  type, public :: source_function_t
    type(planck_lookup_t), pointer  :: lookup => null() 
    character(len=16)        :: source ='Planck' ! choice of source funcion
    procedure (source_Stefan), pointer :: solve => null()
  contains
    procedure :: init
    procedure :: flook
  end type
  !  -------------------------------------------------------------------
CONTAINS

!===============================================================================
!> Initialize a source function object, which includes a choice of source function
!===============================================================================
SUBROUTINE init (self, source)
  class(source_function_t) :: self
  character(len=16) :: source
  integer           :: i
  real              :: x, dx
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  call trace_begin('source_function_mod::Init',3)
  self%source = source  
  select case(source)
  case('Planck')
    allocate(self%lookup)
    self%lookup%xmin = 1d-3
    self%lookup%xmax = 2d1
    self%lookup%n    = 1000
    self%lookup%cx = real(self%lookup%n,kind=8)/(self%lookup%xmax - self%lookup%xmin)
    dx = (self%lookup%xmax - self%lookup%xmin)/real(self%lookup%n,kind=8)
    allocate (self%lookup%table(0:self%lookup%n))
    do i = 0, self%lookup%n
      x = self%lookup%xmin + i*dx
      self%lookup%table(i) = x/(exp(x)-1.0)
    end do
    self%solve => source_Planck
  case('Stefan')
    self%solve => source_Stefan
  case('eos_lookup')
    self%solve => source_eos_lookup
  case default
    print *, "source_function_mod::init: source function",source," not recognised"
    error stop "stopping"
  end SELECT
  call trace_end
END SUBROUTINE init

!=======================================================================
!> Source function is looked up from the EOS table. Suitable for Stagger
!> solver.
!=======================================================================
SUBROUTINE source_eos_lookup (self, rho, entr, wl, src)
  class(source_function_t) :: self
  real, dimension(:,:,:) :: rho
  real, dimension(:,:,:), optional :: entr
  real, dimension(:), optional :: wl
  real, dimension(:,:,:,:) :: src
!-----------------------------------------------------------------------
  call trace_begin('source_function_mod::source_eos_lookup',2)
  call eos%lookup(rho=rho, ee=entr,src=src)
  call trace_end
END SUBROUTINE source_eos_lookup

!=======================================================================
!> Simplest source function. The integrated Planck function is proportional
!> to Stefan's constant times the fourth power of the temperature, divided by pi.
!> This is the relevant source function when the opacity is constant as a
!> function of wavelength, so one can integrate over wavelength before solving
!> the radiative transfer equation.
!=======================================================================
SUBROUTINE source_Stefan (self, temp, dummy, wl, src)
  class(source_function_t) :: self
  real, dimension(:,:,:) :: temp
  real, dimension(:), optional     :: wl
  real, dimension(:,:,:), optional :: dummy
  real, dimension(:,:,:,:) :: src
  real (kind=8):: lambda, sc_stf
  integer:: iwl
  !.............................................................................
  call trace_begin('source_function_mod::source_Stefan',2)
  sc_stf = cgs%stefan/cgs%pi*scaling%temp**4/scaling%f
  src(:,:,:,1) = sc_stf*temp**4
  if (size(src,4)>1) then
    error stop "source_function_mod::source_Stefan: n_wl_bins >1, but Stefan source function is selected."
  end if
  call trace_end
END SUBROUTINE source_Stefan


!===============================================================================
!> Set the source function equal to a Planck function lookup table. Linear
!> interpolation in linear wavelength, for maximum speed.  To handle
!> uncommon cases where some x are below xmin a test is made on x, before
!> choosing how to interpolate.  This is only slightly slower than if one knows
!> everything is in range.
!===============================================================================
SUBROUTINE source_Planck (self, temp, dummy, wl, src)
  class(source_function_t) :: self
  real, dimension(:,:,:)   :: temp
  real, dimension(:), optional     :: wl
  real, dimension(:,:,:), optional :: dummy
  real, dimension(:,:,:,:) :: src
  !.....................................................................
  real,    dimension(size(temp,1),size(temp,2),size(temp,3)) :: p
  integer, dimension(size(temp,1),size(temp,2),size(temp,3)) :: k
  real (kind=8) :: c1, c2
  integer :: i
  !---------------------------------------------------------------------
  call trace_begin('source_function_mod::source_Planck',2)
  c2 = 2.*(cgs%k_b*scaling%temp)**4/(cgs%h_p**3*cgs%c**2)/scaling%f
  
  do i = 1, size(wl)
    c1 = cgs%h_p*cgs%c/(wl(i)*scaling%l*cgs%k_b*scaling%temp)
    p = c1/temp
    if (all(p.ge.self%lookup%xmin)) then
      src(:,:,:,i) = c2 * temp**4*p**3
      p = (p-self%lookup%xmin)*self%lookup%cx
      k = min (self%lookup%n-1, int(p))
      p = p-k
      src(:,:,:,i) = src(:,:,:,i)* &
                      (p*self%flook(k+1) + (1.0-p)*self%flook(k))
    else
      where (p < self%lookup%xmin)
        src(:,:,:,i) = c2*temp**4*p**3*self%lookup%table(0)
      else where
            src(:,:,:,i) = c2*temp**4*p**3
            p = (p-self%lookup%xmin)*self%lookup%cx
            k  = min (self%lookup%n-1, int(p))
            p  = p-k
            src(:,:,:,i) = src(:,:,:,i)* &
                        (p*self%flook(k+1) + (1.0-p)*self%flook(k))
      end where
    end if
  end do
  
  call trace_end  
END SUBROUTINE source_Planck

FUNCTION flook(self,k) result (pp)
  class(source_function_t)  :: self
  integer, dimension(:,:,:) :: k
  real,    dimension(size(k,1),size(k,2),size(k,3)) :: pp
  integer :: y, z
  ! --------------------------------------------------------------------
  do z = 1, size(k,3)
    do y  = 1, size(k,2)
      pp(:,y,z) = self%lookup%table(k(:,y,z))
    end do
  end do
END FUNCTION flook

!> OUTDATED:
!===============================================================================
!> The planck function, normalized to give radiation per unit \ln\lambda,  for
!> a given temperature and wavelength. Input and output are in code units, but
!> we use lambda in micrometers, to avoid underflow and divide fault in the
!> expression when computing lambda**4.
!===============================================================================
SUBROUTINE source_Pl_fn (self, n, lambda, temperature, source)
  class(source_function_t):: self
  integer:: n(2)
  real:: lambda, lambda_mum, lambda_4
  real, dimension(n(1),n(2)):: temperature, source
  integer :: i
  real(8) :: c1, c2
  real, dimension(n(1)):: e
  !.............................................................................
  c1 = 2d0*cgs%h_p*cgs%c**2/(cgs%mum**4)             ! 2*h*c**2/lambda**4
  c2 = cgs%h_p*cgs%c/(cgs%k_b*cgs%mum)              ! h*c/(lambda*T)
  lambda_mum = lambda*scaling%l/cgs%mum                 ! lambda in micrometers
  lambda_4 = lambda_mum**4
  do i=1,n(2)
    e(:) = exp(-c2/(temperature(:,i)*lambda_mum))
    source(:,i) = c1*e(:)/((1.0-e(:))*lambda_4)
  end do
END SUBROUTINE source_Pl_fn

END MODULE source_function_mod
