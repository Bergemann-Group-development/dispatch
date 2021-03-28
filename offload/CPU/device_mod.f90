!===============================================================================¨
!===============================================================================¨
MODULE device_mod
  implicit none
  private
  type, public:: device_t
    integer:: n
    character(16):: type='CPU'
  contains
    procedure:: init
    procedure:: update
  end type
  type(device_t), public:: device
CONTAINS

!===============================================================================¨
!===============================================================================¨
SUBROUTINE init(self)
  class(device_t):: self
END SUBROUTINE init

!===============================================================================¨
!===============================================================================¨
SUBROUTINE update (self, bunch)
  class(device_t):: self
  real, pointer:: bunch(:,:,:,:,:)
END SUBROUTINE update

END MODULE device_mod
