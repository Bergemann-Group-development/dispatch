MODULE mod
CONTAINS
subroutine sub(a)
  real, optional:: a
  if (present(a)) then
    print*,a
  else
    print*,'not present'
  end if
end subroutine
END MODULE mod

USE mod
  call sub
  call sub (a=2.)
END
