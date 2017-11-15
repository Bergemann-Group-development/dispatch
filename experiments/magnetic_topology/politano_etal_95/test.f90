MODULE mod
CONTAINS
SUBROUTINE sub (a)
  real, dimension(:):: a
  print *,a
END
END MODULE mod

PROGRAM test
  USE mod
  real, dimension(10), target:: a, b
  real, dimension(:), pointer:: pa, pb
  real:: c, d
  c = 0.0
  d = -c
  print *,c, d
  a = 0.0
  b = -a
  pa => a
  pb => b
  call sub (a)
  call sub (b)
  call sub (pa)
  pa => pb
  call sub (pa)
END PROGRAM test
