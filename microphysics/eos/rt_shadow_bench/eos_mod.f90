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
  implicit none
  private
  type, public:: eos_t
    real:: gamma = 1.4                  ! in case eos%init is never called
    integer :: eos_type = 0             ! 0 - ideal, x - any other
    real :: t0    = 290.
    real :: rho0 = 1.
    real :: sc_stf
    real :: n_lambda = 1.
    real :: w_lambda = 1.
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
  real, save:: gamma=1.4
  namelist /eos_params/ gamma
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  if (io%master) print *, &
 '------------------------ ideal equation-of-state -----------------------------'
  call scaling%init
  !$omp critical (init_cr)
  if (first_time) then
    gamma = self%gamma
    rewind (io%input); read (io%input, eos_params, iostat=iostat)
    if (io%master) write (io%output, eos_params)
    self%gamma = gamma
    self%t0   = 1./(290./scaling%temp)
    self%rho0 = 1./(1.  /scaling%d)
    print *, cgs%stefan, cgs%pi, scaling%temp, scaling%f
    self%sc_stf = cgs%stefan/cgs%pi*scaling%temp**4/scaling%f
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
SUBROUTINE lookup (self, dim, x, lny, lne, d, pg, tt, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, optional:: gamma
  real, dimension(:,:,:),   optional, pointer :: x, lny, d, lne, pg, tt
  real, dimension(:,:,:,:), optional, pointer :: src, rk
  real, dimension(:,:,:),             pointer :: ee_loc, d_loc
  integer:: i
  !-----------------------------------------------------------------------------
  if (present(gamma)) then
    self%gamma = gamma
  end if
  !-----------------------------------------------------------------------------
  if (present(d)) then
    d_loc => d
  else if (present(lny)) then
    allocate (d_loc(dim(1),dim(2),dim(3)))
    d_loc = exp(lny)
  else
    call mpi%abort ('eos_t: missing density')
  end if
  !-----------------------------------------------------------------------------
  if (present(x)) then
    ee_loc => x
  else if (present(lne)) then
    allocate (ee_loc(dim(1),dim(2),dim(3)))
    ee_loc = exp(lne)
  else
    call mpi%abort ('eos_t: missing internal energy')
  end if
  !-----------------------------------------------------------------------------
  if (present(tt)) then
    tt = (self%gamma-1.0)*ee_loc
  end if
  if (present(src)) then
    do i=1,size(src,4)
      src(:,:,:,i) = self%sc_stf*((self%gamma-1.0)*ee_loc)**4
    end do
  end if
  if (present(pg)) then
    pg = (self%gamma-1.0)*d_loc*ee_loc
  end if
  if (present(rk)) then
    do i=1,size(rk,4)
      rk(:,:,:,i) = (((self%gamma-1.0)*ee_loc*self%t0*scaling%temp)**(-3.5) * &
       (d_loc*self%rho0*scaling%d)**2)* scaling%m/scaling%l**2  
    end do
  end if
  if (present(lne)) deallocate(ee_loc)
  if (present(lny)) deallocate( d_loc)
END SUBROUTINE lookup

END MODULE eos_mod

