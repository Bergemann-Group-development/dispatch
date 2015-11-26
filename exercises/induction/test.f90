PROGRAM test
  USE tensor_mod
  implicit none
  type(scalar_t):: a, c, e, f
  type(vector_t):: v, u, w
  type(tensor_t):: t
  real:: vel(3)
  real, dimension(:,:,:), pointer:: b
  integer:: n(3)
  !.............................................................................
  n = 5
  call set_dims (n)
  allocate (b(n(1),n(2),n(3)))
  b(:,:,:) = 1.1
  a = 1.2
  a = b
  c = a + 2.0
  e = a + 3
  f = a + e
  v = grad(1.)
  print*,2
  print*,a%f(:,1,1)
  v = grad(a)
  u = grad(a)
  t = grad(v)
  w = cross(u,v)
  vel = [1.,2.,3.]
  w = cross(vel,u)
  w = cross([1d0,2d0,3d0],u)
  w = curl(v,u)
  print*,b(:,1,1)
  print*,a%f(:,1,1)
  print*,c%f(:,1,1)
  print*,e%f(:,1,1)
  print*,f%f(:,1,1)
  print*,a%f(:,1,1)
END PROGRAM test
