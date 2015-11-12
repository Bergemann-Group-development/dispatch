!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for HD task.  Copy this to your experiment directory, or
!> to a specialized solver directory, and modify with relevant parameters and 
!> procedures..
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE hd_mod
  USE patch_mod
  implicit none
  type, extends(patch_t):: hd_t
    real:: nu(4)
    contains
    procedure:: init       => init_hd
    procedure:: allocate   => allocate_hd
    procedure:: deallocate => deallocate_hd
    procedure:: update     => update_hd
    procedure:: print      => print_hd
  end type
CONTAINS

SUBROUTINE init_hd (self)
  class(hd_t):: self
END SUBROUTINE

SUBROUTINE update_hd (self)
  class(hd_t):: self
END SUBROUTINE

SUBROUTINE allocate_hd (self)
  class(hd_t):: self
END SUBROUTINE

SUBROUTINE deallocate_hd (self)
  class(hd_t):: self
END SUBROUTINE

SUBROUTINE print_hd (self)
  class(hd_t):: self
  print*,'print_hd:', self%id, self%time
END SUBROUTINE

END MODULE
