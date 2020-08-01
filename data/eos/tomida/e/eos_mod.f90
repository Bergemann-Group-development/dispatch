!===============================================================================
!> Module using a table prepared by a Python procedure.  This module is
!> initialized while still in single-thread mode, so no OMP locks are not
!> needed. The input data table is in a custom format, which may be written by
!> a custom reader, written in Fortran or Python.
!===============================================================================
MODULE eos_mod
  USE scaling_mod
  USE trace_mod
  USE table_mod
  private
  type, public:: eos_t
    character(len=16) :: name='Tomida+Hori(e)'
    type(table_t)     :: table
    integer :: ne, nd, nv
    integer :: w_lambda = 1
    integer :: n_lambda = 1
    real, pointer:: lnd(:) => null()
    real :: gamma = 1.4
  contains
    procedure :: init
    procedure :: lookup
    procedure :: lookup_table
  end type
  type(eos_t), public:: eos
CONTAINS

!===============================================================================
! Initialize the eos table data, table axes and dimensions
!===============================================================================
SUBROUTINE init (self)
  class(eos_t):: self
  !-----------------------------------------------------------------------------
  call self%table%init ('/data/eos/Tomida+Hori_2016/e/eos.dat')
  self%ne = self%table%n1
  self%nd = self%table%n2
  self%nv = self%table%n3
  self%table%x_min = self%table%x_min - log(scaling%d)
  self%table%y_min = self%table%y_min - log(scaling%ee)
  allocate (self%lnd(self%nd))
  do i=1,self%nd
    self%lnd(i) = self%table%x_min + (i-1)*self%table%x_step
  end do
END SUBROUTINE init

!===============================================================================
!> Generic lookup interface, using x,y-names
!===============================================================================
SUBROUTINE lookup (self, dim, lnx, x, lny, y, pg, tt, ss, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:), intent(in), pointer, optional :: lnx, x, lny, y
  real, dimension(:,:,:),                      optional :: pg, tt, ss
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, optional:: gamma
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup_table')
  call self%table%lookup (dim, lnx, x, lny, y, pg, tt, ss, rk, src, gamma)
  call trace%end ()
END SUBROUTINE lookup

!===============================================================================
!> Lookup interface compatible with the stagger solver calls
!===============================================================================
SUBROUTINE lookup_table (self, dim, e, ee, d, lnd, lne, pg, tt, ne, rk, src, &
  ss, gamma, e_below)
class (eos_t):: self
integer:: dim(3)
real, dimension(:,:,:)  , pointer, optional :: lnd,d,lne,e,ee,tt,ne,ss
real, dimension(:,:,:,:)         , optional :: src,rk
real, dimension(:,:,:)           , optional :: pg
real, optional:: gamma
logical, optional:: e_below
!-----------------------------------------------------------------------------
call trace%begin ('eos_t%lookup_table')
call self%table%lookup (dim, lnd, d, lne, ee, pg, tt, ss, rk, src, gamma)
call trace%end ()
END SUBROUTINE lookup_table 

END MODULE eos_mod
