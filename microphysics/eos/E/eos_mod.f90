!> =============================================================================
!> Equation of state module for any sort of tables, provided by a reader
!> =============================================================================
MODULE eos_mod
  USE io_mod
  USE trace_mod
  USE eos_reader_mod
  USE lagrange_mod
  implicit none
  private
  type, public:: eos_t
    integer:: nd, ne
    real(8):: lnd_min, lnd_max, lnd_stp
    real(8):: lne_min, lne_max, lne_stp
    real(8), dimension(:), pointer:: lnd, lne
    real(8), dimension(:,:,:), pointer:: table
  contains
    procedure:: init
    procedure:: pressure
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
  integer:: i, j, nd, ne
  associate (reader=>eos_reader)
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  call reader%init
  nd = reader%nd
  ne = reader%nT
  self%nd = nd
  self%ne = ne
  allocate (self%lnd(nd), self%lne(ne), self%table(ne,nd,2))
  !-----------------------------------------------------------------------------
  ! Fint extremum values of internal energy
  !-----------------------------------------------------------------------------
  Emin = reader%table(1,1,2)
  Emax = reader%table(1,1,2)
  do j=1,nd
  do i=1,ne
    Emin = min(Emin,reader%table(i,j,2))
    Emax = max(Emax,reader%table(i,j,2))
  end do
  end do
  print *, 'Emin, Emax =', Emin, Emax
  !-----------------------------------------------------------------------------
  ! Construct table x- and y-values
  !-----------------------------------------------------------------------------
  self%lne_min = log(Emin)
  self%lne_max = log(Emax)
  self%lne_stp = (self%lne_max-self%lne_min)/(nd-1.0)
  do i=1,ne
    self%lne(i) = self%lne_min + (i-1)*self%lne_stp
  end do
  self%lnd = log(reader%d)
  self%lnd_min = self%lnd(1)
  self%lnd_max = self%lnd(nd)
  self%lnd_stp = (self%lnd_max-self%lnd_min)/(nd-1.0)
  do j=1,nd
    self%lnd(j) = self%lnd_min + (j-1)*self%lnd_stp
  end do
  !-----------------------------------------------------------------------------
  ! For each density, interpolate in the table to find the T corrsponding to 
  ! the new table E, and then again, to find the pressure at that T
  !-----------------------------------------------------------------------------
  do j=1,nd
    do i=1,ne
      E = exp(self%lne(i))              ! new table E-value
      x => reader%table(:,j,2)          ! x-values = known E-values
      y => reader%T                     ! y-values = known T-values
      T = interpolate (x, y, E)         ! find the corresponding T
      x => reader%T                     ! x-values = table T-values
      y => reader%table(:,j,1)          ! y-values = table P-values
      P = interpolate (x, y, T)         ! find the pressure at T
      self%table(i,j,1) = P             ! store the pressure
      self%table(i,j,2) = T             ! store the temperature
      print '(2i5,1p,3e15.6)', i, j, T, P, E
    end do
  end do
  call trace%end
  end associate
END SUBROUTINE init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION interpolate (x, y, x0) RESULT (out)
  !USE lagrange_mod
  real(8), dimension(:):: x, y
  real(8):: x0, out, p
  integer:: i
  !-----------------------------------------------------------------------------
  !do i=2,size(x)-1
  !  if (x(i) > x0) exit
  !end do
  !p = (x0-x(i-1))/(x(i)-x(i-1))
  !out = (1.0-p)*y(i-1) + p*y(i)
  out = lagrange%sequence (x0, x, y)
END FUNCTION interpolate

FUNCTION pressure (self, lnd, lne) RESULT (out)
  class (eos_t):: self
  real:: lnd, lne, out
END FUNCTION pressure

END MODULE eos_mod

