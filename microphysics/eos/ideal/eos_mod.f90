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
    integer:: nd, ne, nv
    real:: gamma
    real(8):: lnd_min, lnd_max, lnd_stp
    real(8):: lne_min, lne_max, lne_stp
    real(8), dimension(:), pointer:: lnd, lne
    real(8), dimension(:,:,:), pointer:: table
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
  real(8), dimension(:), pointer:: x, y
  real(8):: Emin, Emax, E, T, P
  integer:: i, j, nd, ne, iostat
  real, save:: gamma=5./3.
  namelist /eos_params/ gamma
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
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
SUBROUTINE lookup (self, dim, d, e, pg, tt, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3), i
  real, optional:: gamma
  real, dimension(:,:,:),   optional, pointer :: d, e, pg, tt, rk
  real, dimension(:,:,:,:), optional, pointer :: src
  !-----------------------------------------------------------------------------
  if (present(gamma)) then
    self%gamma = gamma
  end if
  pg = (self%gamma-1.0)*d*e
  if (present(rk)) then
    rk = d
  end if
  if (present(src)) then
    do i=1,size(src,4)
      src(:,:,:,i) = self%temperature(d,e)**4
    end do
  end if
END SUBROUTINE lookup

END MODULE eos_mod

