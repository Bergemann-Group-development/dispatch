!> =============================================================================
!> Program to test functionality of mesh generation and destruction.
!> =============================================================================
PROGRAM test
  USE mesh_mod
  USE io_mod
  USE patch_mod
  implicit none

  class(mesh_t), allocatable, dimension(:) :: gridCart, gridSph, gridCyl
  class(patch_t), pointer :: patch
  class(mesh_t), pointer :: m3
  class(mesh_t), pointer, dimension(:) :: m

  integer :: nTheta, i, l(3), u(3)
  real :: ds(3)
  real(8) :: z
  real(8), parameter :: pi = 3.14159265359
  !.............................................................................
  io%verbose = 1

  !-----------------------------------------------------------------------------
  ! Test collections of meshes
  !-----------------------------------------------------------------------------
  print*,"**Sets of meshes**"
  print*,"======================="

  ! ---- Cartesian; all metric factors should be 1.0 ----
  call MeshFactory(mesh_types%Cartesian, gridCart, s=[1.0d0,1.0d0,1.0d0], &
                      p=[0.0d0,0.0d0,0.0d0])
  call gridCart(1)%print()

  print*,"h31c; result should be: T,F,F"
  do i=1,size(gridCart)
    select type (gridCart)
    class is (mesh_t)
      print*,associated(gridCart(i)%h31c)
    end select
    ! allocate and set staggering indices.
    allocate(gridCart(i)%h(8))
    gridCart(i)%h = 0.0
    gridCart(i)%h(i) = -0.5
  end do
  print*,gridCart%d
  print*,gridCart(1)%h

  print*,"Columns 2 and 3 should be identical:"
  associate (m1r => gridCart(1)%r)
  do i=gridCart(1)%lb,gridCart(1)%ub
    print*,i,gridCart(1)%r(i),m1r(i)
  end do
  end associate
  print*,"======================="

  ! ---- cylindrical ----
  call MeshFactory(mesh_types%cylindrical, gridCyl, n=[1,10,45], &
                      s=[1.0d0,0.5d0,2.0d0*pi], p=[0.0d0,0.0d0,0.0d0], &
                      origin=[0.0d0,0.1d0,0.5d0*pi],ng=[0,3,3])
  call gridCyl(1)%print()
  call gridCyl(2)%print()
  call gridCyl(3)%print()

  print*,"The third column should be the absolute value of the first:"
  select type (gridCyl)
  type is (cylindrical_mesh)
    associate(cmesh => gridCyl(2))
    if (associated(cmesh%h32c)) then
      do i=cmesh%li,cmesh%ui
        print*,i,cmesh%r(i),cmesh%h32c(i),cmesh%dar32c(i)
      end do
    endif
    end associate
  end select

  select type (gridCyl)
  type is (cylindrical_mesh)
    if (associated(gridCyl(1)%dvol1c)) print*,gridCyl(1)%d,gridCyl(1)%dvol1c(1),"(both should be 1.0)"
  end select

  select type (gridCyl)
  type is (cylindrical_mesh)
    associate (cmesh => gridCyl(3))
    print*,associated(cmesh%h32f)," (should be F)"
    print*,associated(cmesh%dvol3c)," (should be T)"
    print*,cmesh%zeta," (should be 0.5 * pi = 1.570796...)"
    end associate
  end select
  print*,"======================="

  ! ---- spherical ----
  nTheta = 10
  call MeshFactory(mesh_types%spherical, gridSph, n=[15,nTheta,90], &
                      s=[1.0d0,pi,2.0d0*pi], p=[0.0d0,0.0d0,0.0d0], origin=[1.0d0,0.0d0,0.0d0])
  call gridSph(1)%print()
  call gridSph(2)%print()
  select type (gridSph)
  type is (spherical_mesh)
    associate (cmesh2 => gridSph(2))
    print*,cmesh2%xi," (should be 0)"
    print*,allocated(cmesh2%zeta)," (should be F)"
    print*,cmesh2%xi," (should be 0)"
    print*,associated(cmesh2%h2c)," (should be F)"
    do i=1,nTheta
      print*,i,cmesh2%h32f(i),sin(cmesh2%p + cmesh2%r(i) - 0.5*cmesh2%d)
    end do
    end associate
  end select

  select type (gridSph)
  type is (spherical_mesh)
    associate(cmesh3 => gridSph(3))
    print*,cmesh3%zeta," (should be 0)"
    print*,allocated(cmesh3%xi)," (should be F)"
    end associate
  end select

  print*,"======================="

  call MeshRecycler(gridCart)
  call MeshRecycler(gridCyl)
  call MeshRecycler(gridSph)

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
  associate (m3r => patch%mesh(3)%r, mm => patch%mesh)
  do i=m3%lb,m3%ub
    print*,i,m(3)%r(i),m3%r(i),patch%mesh(3)%r(i),m3r(i)
  end do

  select type(m3)
  type is (Cartesian_mesh)
    print*,'Cartesian mesh recognised'
  type is (cylindrical_mesh)
    print*,'Cylindrical mesh recognised'
  type is (spherical_mesh)
    print*,'Spherical mesh recognised'
  class default
    print*,'Unknown mesh type!'
  end select
  end associate
  
  deallocate(patch)

  print*,"======================="

END PROGRAM
