!> =============================================================================
!> Equation of state module for any sort of tables, provided by a reader
!> =============================================================================
MODULE eos_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  implicit none
  private
  type, public:: eos_t
    real:: gamma
    integer :: eos_type = 0                                             ! 0 - ideal, x - any other
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
  real, save:: gamma=5./3.
  namelist /eos_params/ gamma
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  if (io%master) print *, &
 '------------------------ ideal equation-of-state -----------------------------'
  rewind (io%input); read (io%input, eos_params, iostat=iostat)
  if (iostat == 0) then
    if (io%master) write (io%output, eos_params)
  end if
  self%gamma = gamma
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
SUBROUTINE lookup (self, ee, lnd, d, pg, tt, rk, src, gamma)
  class (eos_t):: self
  real, optional:: gamma
  real, dimension(:,:,:)            , pointer :: ee
  real, dimension(:,:,:),   optional, pointer :: lnd, d, pg, tt, rk
  real, dimension(:,:,:,:), optional, pointer :: src
  integer:: dim(3), i
  !-----------------------------------------------------------------------------
  dim = shape(ee)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  if (present(tt)) then
    tt = (self%gamma-1.0)*ee
  else if (present(src)) then
    do i=1,size(src,4)
      src(:,:,:,i) = self%temperature(lnd,ee)**4
      src(:,:,:,i) = ((self%gamma-1.0)*ee)**4
    end do
  else if (present(d)) then
    if (present(pg)) then
      pg = (self%gamma-1.0)*d*ee
    end if
    if (present(rk)) then
      rk = d
    end if
  else if (present (lnd)) then
    if (present(pg)) then
      pg = (self%gamma-1.0)*exp(lnd)*ee
    end if
    if (present(rk)) then
      rk = exp(lnd)
    end if
  else
    call mpi%abort ('missing d or lnd in eos_mod::lookup')
  end if
END SUBROUTINE lookup

END MODULE eos_mod

