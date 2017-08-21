!> =============================================================================
!> Demonstate that calling a procedure that only knows about the basic type in
!> a class, one can still reach the procedure of an extended type
!> =============================================================================
MODULE a
  type:: a_t
    real:: a=1.
  contains
    procedure:: info => info_a
  end type
CONTAINS
SUBROUTINE info_a (self)
  class(a_t):: self
  print *,'info a=', self%a
END SUBROUTINE info_a
END MODULE

!> =============================================================================
!> Extended object, in a separate module
!> =============================================================================
MODULE b
  USE a
  type, extends(a_t):: b_t
    real:: b=2
  contains
    procedure:: info => info_b
  end type
CONTAINS
SUBROUTINE info_b (self)
  class(b_t):: self
  print *,'info b=', self%b
END SUBROUTINE info_b
END MODULE

!> =============================================================================
!> Module that only knows about the basic class, while the caller may know more
!> =============================================================================
MODULE c
CONTAINS
SUBROUTINE info (f)
  USE a
  class(a_t):: f
  call f%info
END SUBROUTINE info
END MODULE c

!> =============================================================================
!> Watch what happens in the calls to info:
!> =============================================================================
PROGRAM test
  call test1
  call test2
END PROGRAM test

!> =============================================================================
!> Watch what happens in the calls to info:
!> =============================================================================
SUBROUTINE test1
USE a
USE b
USE c
  type(a_t):: aa
  type(b_t), target:: bb
  class(a_t), pointer:: p
  print *,'call from test1, which knows b_t'
  call aa%info
  call bb%info
  print *,'call from test1 via c, which does not know b_t'
  call info (aa)
  call info (bb)
  p => bb
  call info (p)
END SUBROUTINE

!> =============================================================================
!> Watch what happens in the calls to info:
!> =============================================================================
SUBROUTINE test2
USE a
USE c
  type(a_t), target:: aa
  class(a_t), pointer:: p
  print *,'call from test2, which does not know b_t'
  p => aa
  call info (p)
END SUBROUTINE
