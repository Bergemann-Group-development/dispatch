!===============================================================================
!> Module using a table prepared by a Python procedure
!===============================================================================
MODULE eos_mod
  USE scaling_mod
  USE table_mod
  private
  type, public:: eos_t
    character(len=16) :: name='Tomida+Hori(s)'
    type(table_t)     :: table
  contains
    procedure         :: init
    procedure         :: lookup
  end type
  type(eos_t), public:: eos, eose
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(eos_t):: self
  !-----------------------------------------------------------------------------
  call self%table%init ('/data/eos/Tomida+Hori_2016/s/eos.dat')
  call eose%table%init ('/data/eos/Tomida+Hori_2016/e/eos.dat')
  self%table%x_min = self%table%x_min - log(scaling%d)
  eose%table%x_min = eose%table%x_min - log(scaling%d)
  eose%table%y_min = eose%table%y_min - log(scaling%ee)
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE lookup (self, dim, lnx, x, lny, y, pg, tt, ss, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:), intent(in), pointer, optional :: lnx, x, lny, y
  real, dimension(:,:,:),                      optional :: pg, tt, ss
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, optional:: gamma
  !-----------------------------------------------------------------------------
  call self%table%lookup (dim, lnx, x, lny, y, pg, tt, ss, rk, src, gamma)
END SUBROUTINE lookup

END MODULE eos_mod
