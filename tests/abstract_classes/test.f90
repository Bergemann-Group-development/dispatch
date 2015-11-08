!> =============================================================================
!> Template hieararchy, with the first abstract level minimal, the 2nd level
!> has an added procedure, and the 3rd and 4th non-abstract levels add one 
!> more procedure, which is overridden in the 4th level.  This is simlar to
!> what we want to do with the task structure levels
!> =============================================================================
MODULE mod
type, abstract:: a_t
  real:: a
end type
type, abstract, extends(a_t):: aa_t
  contains
  procedure:: sub_a
end type
type, extends(aa_t):: b_t
  real:: b
  contains
  procedure:: sub_b
end type
type, extends(b_t):: c_t
  real:: c
  contains
  procedure:: sub_b => sub_c
end type
!> =============================================================================
!> Note that the arguments below are declared with "class", and that the b_t
!> type argument being passed in the 2nd argument can be treated as a
!> class(a_t), as long as we only refer to the %a element, which is present in
!> all types in the class.
!> =============================================================================
CONTAINS
!> =============================================================================
SUBROUTINE sub(aa)
  class(a_t):: aa
  print*,'sub  ',aa%a
END SUBROUTINE
!> =============================================================================
SUBROUTINE sub_a (self, aa)
  class(aa_t):: self
  class(a_t):: aa
  print*,'sub_a',aa%a
END SUBROUTINE
!> =============================================================================
SUBROUTINE sub_b (self, aa)
  class(b_t):: self
  class(a_t):: aa
  print*,'sub_b',aa%a
END SUBROUTINE
!> =============================================================================
SUBROUTINE sub_c (self, aa)
  class(c_t):: self
  class(a_t):: aa
  print*,'sub_c', aa%a
END SUBROUTINE
END MODULE

!> =============================================================================
!> Note that the declarations in the main program should be with "type", not
!> with "class"!
!> =============================================================================
PROGRAM prog
  USE mod
  type(b_t):: bb
  type(c_t):: cc
  bb%a = 2.0
  bb%b = 1.0
  cc%c = 3.9
  call sub(bb)
  call bb%sub_a(bb)
  call bb%sub_b(bb)
  call cc%sub_a(bb)
  call cc%sub_b(bb)
END
