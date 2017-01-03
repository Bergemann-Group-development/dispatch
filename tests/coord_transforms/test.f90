!> =============================================================================
!> Program to test coordinate transformation functions.
!> =============================================================================
PROGRAM test
  USE io_mod
  USE kinds_mod
  USE mesh_mod
  USE mhd_mod
  USE patch_mod
  implicit none

  integer:: size_seed
  integer, dimension(8):: values
  integer, dimension(:), allocatable:: seed
  real(kind=KindScalarVar):: pi = asin(1.0) * 2.0
  real(kind=KindScalarVar), dimension(3):: xp, R, x, x0, R2
  real(kind=KindScalarVar), dimension(:), allocatable:: vp, vCyl, vSph
  class(mhd_t), allocatable:: source
  class(patch_t), pointer:: target
  integer, dimension(3):: ii, jj
  !-----------------------------------------------------------------------------
  call io%init()

  ! initialise random number generator
  call random_seed(size=size_seed)
  allocate (seed(size_seed))
  call random_seed(get=seed)
  call date_and_time(values=values)
  seed = values(8)
  call random_seed(put=seed)
  deallocate(seed)

  !-- test conversion of positions ---------------------------------------------
  print*,"---- cyl => Cart => cyl ----"
  xp = random_coord()
  x0 = random_coord()
  ! rescale so the possible values are 0 -> 2pi
  xp(3) = xp(3) * 2.0 * pi
  print*, "Z,R,Phi :", xp
  print*, "x0,y0,z0:", x0
  x = CylindricalToCartesian(xp(1),xp(2),xp(3))
  x = x - x0
  print*, "x,y,z   :",x
  R = CartesianToCylindrical(x(1),x(2),x(3),OCart=-x0)
  print*, "Z,R,Phi :", R

  print*,"---- sph => Cart => sph ----"
  xp = random_coord()
  x0 = random_coord()
  ! rescale so the possible values are 0 -> pi
  xp(2) = xp(2) * pi
  ! rescale so the possible values are 0 -> 2pi
  xp(3) = xp(3) * 2.0 * pi
  print*, "r,theta,phi:", xp
  print*, "x0,y0,z0:", x0
  x = SphericalToCartesian(xp(1),xp(2),xp(3))
  x = x - x0
  print*, "x,y,z   :",x
  r = CartesianToSpherical(x(1),x(2),x(3),OCart=-x0)
  print*, "r,theta,phi:", r

  print*,"---- Cart => cyl => Cart ----"
  xp = random_coord()
  x0 = random_coord()
  print*, "x,y,z   :",xp
  print*, "x0,y0,z0:", x0
  xp = xp - x0
  R = CartesianToCylindrical(xp(1),xp(2),xp(3))
  print*, "Z,R,Phi :", R
  x = CylindricalToCartesian(R(1),R(2),R(3),OCart=-x0)
  print*, "x,y,z   :",x

  print*,"---- Cart => sph => Cart ----"
  xp = random_coord()
  x0 = random_coord()
  print*, "x,y,z   :",xp
  print*, "x0,y0,z0:", x0
  xp = xp - x0
  r = CartesianToSpherical(xp(1),xp(2),xp(3))
  print*, "r,theta,phi:", r
  x = SphericalToCartesian(r(1),r(2),r(3),OCart=-x0)
  print*, "x,y,z   :",x

  print*,"---- cyl => sph => cyl ----"
  xp = random_coord()
  x0 = random_coord()
  ! rescale so the possible values are 0 -> 2pi
  xp(3) = xp(3) * 2.0 * pi
  print*, "Z,R,Phi :", xp
  print*, "x0,y0,z0:", x0
  r = CylindricalToSpherical(xp(1),xp(2),xp(3))
  print*, "r,theta,phi:", r
  R2 = SphericalToCylindrical(r(1),r(2),r(3))
  print*, "Z,R,Phi :", R2

  !-- test conversion of velocities --------------------------------------------
  allocate (source)
  source%mesh_type = mesh_types%cylindrical
  source%size = [1.0d0, 2.0d0, 2.0d0 * pi]
  source%position = 0.0
  source%n = [10, 64, 180]
  source%llc_nat = [0.0d0, 0.1d0, 0.0d0]
  source%llc_cart = source%position + CurrentToCartesian(source%mesh_type, 0.0d0, 0.1d0, 0.0d0)
  source%nv = Bidxs(3)
  call source%patch_t%init()
  allocate (target, source=source)

  print*,"---- Vcyl => VCart => Vcyl ----"
  vp = random_vector(source%nv)
  allocate (vCyl(source%nv))
  ii = [5,9,60]
  print*, "v_Z,v_R,v_Phi:", vp(sidxs(1):sidxs(3))
  print*, "B_Z,B_R,B_Phi:", vp(Bidxs(1):Bidxs(3))
  vCyl = source%MapVectors(vp, target, ii, ii)
  print*, "v_Z,v_R,v_Phi:", vCyl(sidxs(1):sidxs(3))
  print*, "B_Z,B_R,B_Phi:", vCyl(Bidxs(1):Bidxs(3))

  deallocate (target, source, vCyl, vp)

  print*,"---- Vsph => VCart => Vsph ----"
  allocate (source)
  source%mesh_type = mesh_types%spherical
  source%size = [1.0d0, 0.5d0 * pi, 2.0d0 * pi]
  source%position = 0.0
  source%n = [10, 90, 180]
  source%llc_nat = [0.1d0, 0.0d0, 0.0d0]
  source%llc_cart = source%position + CurrentToCartesian(source%mesh_type, 0.1d0, 0.0d0, 0.0d0)
  source%nv = Bidxs(3)
  call source%patch_t%init()
  allocate (target, source=source)

  vp = random_vector(source%nv)
  allocate (vSph(source%nv))
  ii = [10,45,20]
  print*, "v_r,v_theta,v_phi:", vp(sidxs(1):sidxs(3))
  print*, "B_r,B_theta,B_phi:", vp(Bidxs(1):Bidxs(3))
  vSph = source%MapVectors(vp, target, ii, ii)
  print*, "v_r,v_theta,v_phi:", vSph(sidxs(1):sidxs(3))
  print*, "B_r,B_theta,B_phi:", vSph(Bidxs(1):Bidxs(3))

  deallocate (target, source, vSph, vp)

CONTAINS

!> =============================================================================
!> Generate a triplet of random numbers
!> =============================================================================
FUNCTION random_coord() result (vec)
  real(kind=KindScalarVar), dimension(3):: vec
  
  call random_number(vec(1))
  call random_number(vec(2))
  call random_number(vec(3))

END FUNCTION random_coord

!> =============================================================================
!> Generate `nvar` random numbers
!> =============================================================================
FUNCTION random_vector(nvar) result (vec)
  integer, intent(in):: nvar
  integer:: i, iend
  real(kind=KindScalarVar), dimension(:), allocatable:: vec

  allocate (vec(nvar))
  iend = ubound(vec,dim=1)
  do i=1,iend
    call random_number(vec(i))
  end do
END FUNCTION random_vector

END PROGRAM
