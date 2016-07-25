MODULE microphysics_mod
  USE scaling_mod
  USE units_mod
  USE io_mod
  USE opacities_mod
  USE source_function_mod
  USE trace_mod
  implicit none
  private
  type, public :: microphysics_t
    type (opacity_t) :: opacity 
    type (source_function_t) :: source_fn
    integer :: n_wl_bins
    real, dimension(:), pointer :: wl_bins ! (begin|end,2)
    real(kind=8)                :: cr_heating = 0     ! cosmic ray heating
    ! add other types, e.g. chemistry when needed
  contains
    procedure :: init
    procedure :: dtau
    procedure :: set_wl_bins
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


!> =============================================================================
!> compute dtau given opacity and step ds along the ray
!> =============================================================================
FUNCTION dtau (self, ds, opacity, wlens)
  class(microphysics_t) :: self
  real, dimension(:,:)    :: ds
  real, dimension(:,:,:)  :: opacity
  real, dimension(:)   :: wlens
  real, dimension(size(opacity,1), size(opacity,2), size(wlens)):: dtau
  integer :: i, j
  !-----------------------------------------------------------------------------
  ! in DISPATCH ray geometry points are on the cell borders, so dtau doesn't
  ! have to be downstaggered and is siply opacity*optical path length
  !-----------------------------------------------------------------------------
  call trace_begin('microphysics_mod::dtau')
  do i = 1, size(wlens)
    dtau(:,:,i) = ds(:,:)*opacity(:,:,i)
  end do
  call trace_end
END FUNCTION dtau

SUBROUTINE set_wl_bins(self)
  class(microphysics_t) :: self
  !---------------------------------------------------------------------
  ! FIXME add appropriate subroutine for n-bins
  !---------------------------------------------------------------------
  if (self%n_wl_bins ==1) then
    self%wl_bins(1) = 500.
  end if
END SUBROUTINE set_wl_bins

!> =============================================================================
END MODULE microphysics_mod
