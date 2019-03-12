!===============================================================================
!===============================================================================
MODULE opacity_mod
  USE io_mod
  implicit none
  private
  type, public:: opacity_t
  contains
    procedure:: init
  end type
  type (opacity_t), public:: opacity
CONTAINS

!===============================================================================
!> Read in the data from a specific table file
!===============================================================================
SUBROUTINE init (self)
  class (opacity_t):: self
END SUBROUTINE init

END MODULE opacity_mod
