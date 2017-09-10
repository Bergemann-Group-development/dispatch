!> =============================================================================
!> Equation of state module for any sort of tables, provided by a reader
!> =====================================================================
!> Simple opacity, expressed as
!> kappa = (T/T0)^-3.5*(rho/rho0)^2
!> T0 = 290K, rho0 = 1 g cm-3
!> =====================================================================
!> =============================================================================
MODULE eos_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE scaling_mod
  USE units_mod
  USE kinds_mod
  implicit none
  private
  type, public:: eos_t
    real:: gamma = 1.6666667                  ! in case eos%init is never called
    integer :: eos_type = 0             ! 0 - ideal, x - any other
    real :: t0    = 290.
    real :: rho0 = 1.
    real :: alfa0 = 1.
    real(8) :: sc_stf
    integer :: n_lambda = 1
    real :: w_lambda = 1.
    logical :: do_rt = .false.
  contains
    procedure:: init
    procedure:: pressure
    procedure:: temperature
    procedure:: lookup
  end type
  type (eos_t), public:: eos
CONTAINS

!> =============================================================================
!> The eos reader is supposed to deliver table date in the eos_reader data
!> type, with pressure and internal energy in an (nT,nd,2) table, with NT
!> temperature value T(nT), and nd density values d(nd).
!> =============================================================================
SUBROUTINE init (self)
  class (eos_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  real, save:: gamma=1.66666667
  logical, save:: do_rt = .false.
  namelist /eos_params/ gamma, do_rt
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  call scaling%init
  !$omp critical (init_cr)
  if (first_time) then
  if (io%master) print *, &
     '------------------------ ideal equation-of-state -----------------------------'
    gamma = self%gamma
    rewind (io%input); read (io%input, eos_params, iostat=iostat)
    if (io%master) write (io%output, eos_params)
    self%gamma = gamma
    self%t0   = 1./(290./scaling%temp)
    self%rho0 = 1./(1.  /scaling%d)
    self%alfa0 = 1. * scaling%m/scaling%l**2
    first_time = .false.
    self%do_rt = do_rt
  end if
  !$omp end critical (init_cr)
  call trace%end
END SUBROUTINE init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION pressure (self, d, e) RESULT (pg)
  class (eos_t):: self
  real, dimension(:,:,:), pointer:: d, e
  real, dimension(size(d,1),size(d,2),size(d,3)):: pg
  !-----------------------------------------------------------------------------
  pg = (self%gamma-1.0)*d*e
END FUNCTION

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION temperature (self, d, e) RESULT (tt)
  class (eos_t):: self
  real, dimension(:,:,:), pointer:: d, e
  real, dimension(size(d,1),size(d,2),size(d,3)):: tt
  !-----------------------------------------------------------------------------
  tt = (self%gamma-1.0)*e
END FUNCTION

!> =============================================================================
!> General EOS lookup routine, returning values in optional arguments
!> =============================================================================
SUBROUTINE lookup (self, dim, lnx, x, lny, y, pg, tt, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, optional:: gamma
  real(kind=KindScalarVar), dimension(:,:,:), intent(in), pointer, optional :: lnx, x, lny, y
  real(kind=KindScalarVar), dimension(:,:,:),                      optional :: pg, tt
  real(kind=KindRTVar), dimension(:,:,:,:),                        optional :: src, rk
  real, dimension(:,:,:), pointer :: x_loc, y_loc
  real(kind=8), dimension(:,:,:), pointer :: tt_loc
  integer:: i
  integer, save:: itimer = 0
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t: lookup', itimer=itimer)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  !-----------------------------------------------------------------------------
  if (present(y)) then
    y_loc => y
  else if (present(lny)) then
    allocate (y_loc(dim(1),dim(2),dim(3)))
    y_loc = exp(lny)
  else
    call mpi%abort ('eos_t: missing density')
  end if
  !-----------------------------------------------------------------------------
  if (present(x)) then
    x_loc => x
  else if (present(lnx)) then
    allocate (x_loc(dim(1),dim(2),dim(3)))
    x_loc = exp(lnx)
  else
    call mpi%abort ('eos_t: missing internal energy')
  end if
  !-----------------------------------------------------------------------------
  allocate (tt_loc(dim(1),dim(2),dim(3)))
  tt_loc = y_loc**self%gamma*exp(x_loc/y_loc*(self%gamma-1))/y_loc
  if (present(tt)) then
    tt = tt_loc
  end if
  if (present(src)) then
    do i=1,size(src,4)
      src(:,:,:,i) = scaling%stefan*tt_loc**4
    end do
  end if
  if (present(pg)) then
    pg = y_loc**self%gamma*exp(x_loc/y_loc*(self%gamma-1))
  end if
  if (present(rk)) then
    do i=1,size(rk,4)
      rk(:,:,:,i) = self%alfa0*(tt_loc*self%t0)**(-3.5)*(y_loc*self%rho0)**2
    end do
  end if
  if (present(lnx)) deallocate(x_loc)
  if (present(lny)) deallocate( y_loc)
  deallocate(tt_loc)
    call trace%end(itimer)
END SUBROUTINE lookup

END MODULE eos_mod

