
real(8):: r(3), v(3), a(3), dt, ra, c1, c2, third, d1, d3
integer(8):: n
real(8), dimension(:), allocatable:: c, d
integer:: m

v0=0.4
orbits=100
print *,'orbits:', orbits

third=1d0/3d0
do m=1,4
print *, 'case:', m
do i=4,0,-1
  m1 = m
  if (m==1) m1=2
  allocate (c(m1),d(m1))
  if (m==1) then
    c = [0.5d0,0.5d0]
    d = [1.0d0,0.0d0]
  else if (m==2) then
    c = [0.0d0,1.0d0]
    d = [0.5d0,0.5d0]
  else if (m==3) then
    c = [1d0,-2d0/3d0,2d0/3d0]
    d = [-1d0/24d0,3d0/4d0,7d0/24d0]
  else
    c1 = 1d0/(2d0*(2d0-2d0**third))
    c2 = (1d0-2d0**third)*c1
    c = [c1,c2,c2,c1]
    d1 = 1d0/(2d0-2d0**third)
    d3 = -2d0**third*d1
    d = [d1,d3,d1,0d0]
  end if
  !print *,sum(c),sum(d)

  r=[0,1,0]
  v=[1,0,0]*v0

  dt = 1e-3*10.**(i*0.5)
  n = orbits/dt*2.*3.1415965

  e0 =  0.5*sum(v**2) - 1d0/sqrt(sum(r**2))
  emin = 0.0
  emax = -1e9

  do j=1,n
    do k=1,m1
      r = r + c(k)*v*dt
      ra = sqrt(sum(r**2))
      a = -r/ra**3
      v = v + d(k)*a*dt
    end do
    e1 =  0.5*sum(v**2) - 1d0/ra
    emin = min(emin,e1)
    emax = max(emax,e1)
    !print *,e1
  end do
  print 1,' order:', m,' dt:', dt,' e0:', e0,' emin:', emin-e0,'  emax:', emax-e0
  1 format(a,i5,3x,8(a,f12.4,4x))
  deallocate (c,d)
end do
end do

END
