SUBROUTINE planck (n, lambda, temperature, source)
  integer:: n
  real, dimension(n,n,n):: temprature, source
  real, dimension(n):: e
  real:: lambda, cp, c2
!...............................................................................
  cp=1.191042870e7                                                              ! first radiation constant: 1.191042870e-16 W m^2 sr^-1
  c2=1.438777e4                                                                 ! second radiation constant [micro-meter Kelvin]
  do iz=1,n
  do iy=1,n
    e(:) = exp(-c2/(temprature(:,iy,iz)*lambda))
    source(:,iy,iz) = cp*e(:)/((1.0-e(:))*lambda**5)
   end do
   end do
END SUBROUTINE
