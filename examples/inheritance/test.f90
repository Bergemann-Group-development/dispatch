!> =============================================================================
!> Example illustrating how procedures are shared between extended object. The
!> basic obect 'a' has an 'assert' procedure, which is inhereted by the extended
!> object 'a_x'.  The assert procedure may be called from any of the two objects,
!> as well as from other object procedures.
!> =============================================================================
MODULE mod
  type:: a_t
  contains
    procedure:: assert
  end type
  type, extends(a_t):: a_x_t
  contains
    procedure sub
  end type
  type(a_t):: a
  type(a_x_t):: a_x
CONTAINS
!> =============================================================================
SUBROUTINE assert (self, text, code)
  class(a_t):: self
  character(len=*):: text
  integer:: code
  print'(a10,i4)', text, code
END SUBROUTINE
!> =============================================================================
SUBROUTINE sub (self)
  class(a_x_t):: self
  call self%assert ('sub', 20)
END SUBROUTINE
END MODULE mod
!> =============================================================================
PROGRAM test
  use mod
  call a%assert ('a', 10)
  call a_x%sub
  call a_x%assert ('a_x', 30)
END PROGRAM test
