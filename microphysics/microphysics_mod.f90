MODULE microphysics_mod
  USE scaling_mod
  USE units_mod
  USE io_mod
  USE opacities_mod
  USE source_function_mod
  USE eos_mod
  USE trace_mod
  implicit none
  private
  type, public :: microphysics_t
    type (opacity_t) :: opacity 
    type (source_function_t) :: source_fn
    integer :: n_wl_bins
    real, dimension(:), pointer :: wl_bins                              ! (begin|end,2)
    real(kind=8)                :: cr_heating = 0                       ! cosmic ray heating
    ! add other types, e.g. chemistry when needed
  contains
    procedure :: init
    procedure :: set_wl_bins
    procedure :: get_data
  end type
  
  type(microphysics_t), public:: microphysics
  
CONTAINS

!> =============================================================================
!> Initialize a microphysics object, which includes a choice of source function
!> and opacity, and one or more methods to compute optical depth increments
!> =============================================================================
SUBROUTINE init (self)
  class(microphysics_t) :: self
  logical :: first_time=.true.                  ! namelist on first invocation
  character(len=16) :: source ='Planck'         ! choice of source funcion
  integer           :: i, n_wl_bins = 1, cr_heating = 0
  real              :: x, dx
  namelist /microphysics_params/ source, n_wl_bins, cr_heating
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  call trace_begin('microphysics_mod::Init',2)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, microphysics_params)
    if (io%master) write (*,microphysics_params)
    self%n_wl_bins = n_wl_bins
    allocate(self%wl_bins(n_wl_bins))
    self%cr_heating = cr_heating/cgs%m_u*scaling%d/(scaling%e/scaling%l**3/scaling%t)
    call self%opacity%init()
    call self%source_fn%init(source)
  end if
  call trace_end
END SUBROUTINE init

SUBROUTINE get_data(self, rho, entropy, source_fn, opacity, gam, cs)
  class(microphysics_t) :: self
  real, dimension(:,:,:), pointer :: rho, entropy, temp
  real, dimension(:,:,:,:), pointer, optional :: opacity, source_fn
  real, optional :: gam, cs
  integer :: it
  ! --------------------------------------------------------------------
  select case(eos%eos_type)
  case(1)  ! ideal gas eos
    allocate(temp(size(rho,1),size(rho,2),size(rho,3)))
    call eos%lookup(rho = rho, s=entropy, tt=temp, gam = gam, csound = cs)
    call self%opacity%solve(rho, temp, self%n_wl_bins, opacity)
    call self%source_fn%solve(temp=temp,wl=self%wl_bins,src=source_fn)
    
  case(2)  ! stagger eos
    call eos%lookup(rho=rho,s=entropy,src=source_fn,rk=opacity(:,:,:,1))
    if (size(opacity,4)>1) then
      do it =2, size(opacity,4)
        opacity(:,:,:,it) = opacity(:,:,:,it-1)*10.
      end do
    end if
    
  case default
    error stop "microphysics_mod::get_data: eos_type not set"
  end select

END SUBROUTINE get_data

SUBROUTINE set_wl_bins(self)
  class(microphysics_t) :: self
  integer :: bin
  !---------------------------------------------------------------------
  ! FIXME add appropriate subroutine for n-bins. Doesn't matter for 
  ! Marconi.
  !---------------------------------------------------------------------
  do bin =1, self%n_wl_bins
    self%wl_bins(bin) = 500.
  end do
END SUBROUTINE set_wl_bins

!> =============================================================================
END MODULE microphysics_mod
