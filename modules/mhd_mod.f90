!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for MHD task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mhd_mod
  USE hd_mod
  implicit none
  type, extends(hd_t):: mhd_t
    real:: eta(4)
    contains
    procedure:: init       => init_mhd
    procedure:: allocate   => allocate_mhd
    procedure:: deallocate => deallocate_mhd
    procedure:: update     => update_mhd
    procedure:: print      => print_mhd
  end type
CONTAINS

SUBROUTINE init_mhd (self, id)
  class(mhd_t):: self
  integer:: id
END SUBROUTINE

SUBROUTINE update_mhd (self)
  class(mhd_t):: self
END SUBROUTINE

SUBROUTINE allocate_mhd (self)
  class(mhd_t):: self
END SUBROUTINE

SUBROUTINE deallocate_mhd (self)
  class(mhd_t):: self
END SUBROUTINE

SUBROUTINE print_mhd (self)
  class(mhd_t):: self
  print*,'print_mhd:', self%id, self%time
END SUBROUTINE

END MODULE
