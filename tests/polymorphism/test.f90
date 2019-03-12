!===============================================================================
!> Demonstration of polymorphism, and the corresonding overloading of procedures.
!===============================================================================
MODULE mod
 private
 !------------------------------------------------------------------------------
 ! Object defines a procedure "a"
 !------------------------------------------------------------------------------
 type, public:: a_t
   real:: ar=0.0
 contains
   procedure:: a
 end type
 !------------------------------------------------------------------------------
 ! Object overrides procedure "a"
 !------------------------------------------------------------------------------
 type, public, extends(a_t):: b_t
   real:: br=1.0
 contains
   procedure:: a => b
 end type
 public c
CONTAINS
!> ------------------------------------------------------------------------------
SUBROUTINE a(self)
  class(a_t):: self
  print*,'a',self%ar
END SUBROUTINE
!> ------------------------------------------------------------------------------
SUBROUTINE b(self)
  class(b_t):: self
  print*,'b',self%br
END SUBROUTINE
!> ------------------------------------------------------------------------------
SUBROUTINE c(x)
  class(a_t):: x
  call x%a
END SUBROUTINE
END MODULE

!> ------------------------------------------------------------------------------
!> The explicit calls to a%a and b%a show the expected output, and indirect 
!> call via the c subroutine shows the same choice of procedures.
!> ------------------------------------------------------------------------------
PROGRAM test
  USE mod
  type(a_t):: a
  type(b_t):: b
  call a%a
  call b%a
  call c(a)
  call c(b)
END
