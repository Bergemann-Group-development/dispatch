!> =============================================================================
!> Program to test functionality of mesh generation and destruction.
!> =============================================================================
PROGRAM test
  USE io_mod
  USE mesh_mod
  USE patch_mod
  implicit none

  class(mesh_t), allocatable:: meshSingle
  class(mesh_t), allocatable, dimension(:):: meshCart, meshSph, meshCyl
  class(patch_t), pointer :: patch
  class(mesh_t), pointer :: m3
  class(mesh_t), pointer, dimension(:) :: m

  integer :: nTheta, i, l(3), u(3)
  real :: ds(3)
  real(8) :: z
  real(8), parameter :: pi = asin(1.0d0) * 2.0d0
  !-----------------------------------------------------------------------------
  io%verbose = 1

  !-----------------------------------------------------------------------------
  ! Test an individual meshes
  !-----------------------------------------------------------------------------
  print*,repeat("=",70)
  print*,"**Individual meshes**"
  print*,repeat("-",25)

  allocate (cylindrical_mesh:: meshSingle)
  meshSingle%type = 'cylindrical_mesh'
  meshSingle%idir = 3
  call meshSingle%init(10,2,1.0d0,0.5d0,0.0d0)
  call meshSingle%init_geometry(3,0.0d0)
  call meshSingle%print()
  print*,meshSingle%dvol3f(1)
  deallocate (meshSingle)

  allocate (spherical_mesh:: meshSingle)
  meshSingle%type = 'spherical_mesh'
  meshSingle%idir = 2
  call meshSingle%init(10,2,1.0d0,0.5d0,0.0d0)
  call meshSingle%init_geometry(meshSingle%idir,0.0d0)
  call meshSingle%print()
  print*,associated(meshSingle%dvol3f)," (should be F)"
  if (associated(meshSingle%dvol3f)) print*,meshSingle%dvol3f(1)
  print*,meshSingle%h32f(3)
  deallocate (meshSingle)

  !-----------------------------------------------------------------------------
  ! Test collections of meshes
  !-----------------------------------------------------------------------------
  print*,repeat("=",70)
  print*,"**Sets of meshes**"

  ! ---- Cartesian; all metric factors should be 1.0 ----
  call MeshFactory(mesh_types%Cartesian, meshCart, s=[1.0d0,1.0d0,1.0d0], &
                      p=[0.0d0,0.0d0,0.0d0])
  call meshCart(1)%print()

  print*,"h31c; result should be: T,F,F"
  do i=1,size(meshCart)
    select type (meshCart)
    class is (mesh_t)
      print*,associated(meshCart(i)%h31c)
    end select
    ! allocate and set staggering indices.
    allocate(meshCart(i)%h(8))
    meshCart(i)%h = 0.0
    meshCart(i)%h(i) = -0.5
  end do
  print*,meshCart%d
  print*,meshCart(1)%h

  print*,"Columns 2 and 3 should be identical:"
  associate (m1r => meshCart(1)%r)
  do i=meshCart(1)%lb,meshCart(1)%ub
    print*,i,meshCart(1)%r(i),m1r(i)
  end do
  end associate
  print*,"======================="

  ! ---- cylindrical ----
  call MeshFactory(mesh_types%cylindrical, meshCyl, n=[1,10,45], &
                      s=[1.0d0,0.5d0,2.0d0*pi], p=[0.0d0,0.0d0,0.0d0], &
                      origin=[0.0d0,0.1d0,0.5d0*pi],ng=[0,3,3])
  call meshCyl(1)%print()
  call meshCyl(2)%print()
  call meshCyl(3)%print()

  print*,"The third column should be the absolute value of the first:"
  select type (meshCyl)
  type is (cylindrical_mesh)
    associate(cmesh => meshCyl(2))
    if (associated(cmesh%h32c)) then
      do i=cmesh%li,cmesh%ui
        print*,i,cmesh%r(i),cmesh%h32c(i),cmesh%dar32c(i)
      end do
    endif
    end associate
  end select

  select type (meshCyl)
  type is (cylindrical_mesh)
    if (associated(meshCyl(1)%dvol1c)) print*,meshCyl(1)%d,meshCyl(1)%dvol1c(1),"(both should be 1.0)"
  end select

  select type (meshCyl)
  type is (cylindrical_mesh)
    associate (cmesh => meshCyl(3))
    print*,associated(cmesh%h32f)," (should be F)"
    print*,associated(cmesh%dvol3c)," (should be T)"
    print*,cmesh%zeta," (should be 0.5 * pi = 1.570796...)"
    end associate
  end select
  print*,"======================="

  ! ---- spherical ----
  nTheta = 10
  call MeshFactory(mesh_types%spherical, meshSph, n=[15,nTheta,90], &
                      s=[1.0d0,pi,2.0d0*pi], p=[0.0d0,0.0d0,0.0d0], origin=[1.0d0,0.0d0,0.0d0])
  call meshSph(1)%print()
  call meshSph(2)%print()
  select type (meshSph)
  type is (spherical_mesh)
    associate (cmesh2 => meshSph(2))
    print*,cmesh2%xi," (should be 0)"
    print*,allocated(cmesh2%zeta)," (should be F)"
    print*,cmesh2%xi," (should be 0)"
    print*,associated(cmesh2%h2c)," (should be F)"
    do i=1,nTheta
      print*,i,cmesh2%h32f(i),sin(cmesh2%p + cmesh2%r(i) - 0.5*cmesh2%d)
    end do
    end associate
  end select

  select type (meshSph)
  type is (spherical_mesh)
    associate(cmesh3 => meshSph(3))
    print*,cmesh3%zeta," (should be 0)"
    print*,allocated(cmesh3%xi)," (should be F)"
    end associate
  end select

  print*,"======================="

  call MeshRecycler(meshCart)
  call MeshRecycler(meshCyl)
  call MeshRecycler(meshSph)

  !-----------------------------------------------------------------------------
  ! Test integration with patches
  !-----------------------------------------------------------------------------
  print*,"**Test integration with patches**"
  print*,"======================="

  allocate(patch_t :: patch)
  call MeshFactory(mesh_types%cartesian, patch%mesh)
  call patch%mesh(1)%print()
  call patch%mesh(2)%print()

! The following works in ifort but not in gfortran (v5.3); as far as I can tell,
! it is valid Fortran 2003.
!  associate (m => patch%mesh))
!    print*,shape(m),size(m),'(should be 3 3)'
!  end associate

  call MeshRecycler(patch%mesh)
  deallocate(patch)

  allocate(patch_t :: patch)
  patch%n = [2,3,4]
  patch%size = [2.0,3.0,4.0]
  patch%position = [0.25,0.30,0.35]
  call MeshFactory(mesh_types%cartesian, patch%mesh, n=patch%n, s=patch%size, &
                   p=patch%position)
  call patch%mesh(1)%print()
  call patch%mesh(2)%print()

  print*,lbound(patch%mesh),ubound(patch%mesh)
  print*,"The following three lines should be identical:"
  l = patch%mesh%lb
  u = patch%mesh%ub
  print*,patch%mesh(3)%lb,patch%mesh(3)%ub
  print*,l(3),u(3)

  do i=1,3
    l(i) = patch%mesh(i)%lb
    u(i) = patch%mesh(i)%ub
  end do
  print*,l(3),u(3)

  print*,"Columns 2, 3, 4 and 5 should be identical:"
  m => patch%mesh
  m3 => patch%mesh(3)
  associate (m3r => patch%mesh(3)%r)
  do i=m3%lb,m3%ub
    print*,i,m(3)%r(i),m3%r(i),patch%mesh(3)%r(i),m3r(i)
  end do
  end associate

  select type(m3)
  type is (Cartesian_mesh)
    print*,'m3: Cartesian mesh recognised'
  type is (cylindrical_mesh)
    print*,'m3: Cylindrical mesh recognised'
  type is (spherical_mesh)
    print*,'m3: Spherical mesh recognised'
  class default
    print*,'m3: Unknown mesh type!'
  end select

  ! NOTE!: the following code causes an ICE in gfortran v5.4.1 and earlier.
!  associate (m2 => patch%mesh(2))
!  select type(m2)
!  type is (Cartesian_mesh)
!    print*,'m2: Cartesian mesh recognised'
!  type is (cylindrical_mesh)
!    print*,'m2: Cylindrical mesh recognised'
!  type is (spherical_mesh)
!    print*,'m2: Spherical mesh recognised'
!  class default
!    print*,'m3: Unknown mesh type!'
!  end select
!  end associate
  
  deallocate(patch)

  print*,"======================="

END PROGRAM
