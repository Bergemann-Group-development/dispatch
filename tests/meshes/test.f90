!===============================================================================
!> Program to test functionality of mesh generation and destruction.
!===============================================================================
PROGRAM test
  USE io_mod
  USE mesh_mod
  USE patch_mod
  implicit none

  class(mesh_t), pointer:: meshSingle => null()
  class(mesh_t), pointer, dimension(:):: meshCart => null(), meshSph => null(), meshCyl => null()
  type(Cartesian_mesh), pointer, dimension(:):: meshCart2 => null()
  type(Cartesian_mesh), dimension(3), target:: meshCart3
  class(patch_t), pointer :: patch => null()
  class(mesh_t), pointer, dimension(:) :: m => null()

  integer :: nTheta, i, l(3), u(3)
  real :: ds(3)
  real(8) :: z
  real(8), parameter :: pi = asin(1.0d0) * 2.0d0
  !-----------------------------------------------------------------------------
  io%verbose = 1

  !-----------------------------------------------------------------------------
  ! Tests of individual meshes
  !-----------------------------------------------------------------------------
  print*,repeat("=",70)
  print*,"**Individual meshes**"
  print*,repeat("-",20)

  allocate (cylindrical_mesh:: meshSingle)
  meshSingle%type = 'cylindrical_mesh'
  meshSingle%idir = 3
  call meshSingle%init(11,2,1.0d0,0.5d0,0.0d0)
  call meshSingle%init_geometry(3,0.0d0)
  call meshSingle%print()
  print'(a,1pg13.6E2,a)',"dvol3f(1) = ",meshSingle%dvol3f(1), " (should be 0.1)"

  print*,"Columns 2 and 3 should be identical"
  associate (msr => meshSingle%r)
  do i=meshSingle%lb,meshSingle%ub
    print*,i,meshSingle%r(i),msr(i)
  end do
  end associate
  deallocate (meshSingle)

  allocate (spherical_mesh:: meshSingle)
  meshSingle%type = 'spherical_mesh'
  meshSingle%idir = 2
  call meshSingle%init(11,s=0.5d0*pi,p=0.0d0)
  call meshSingle%init_geometry(meshSingle%idir,0.0d0)
  call meshSingle%print()
  print*,associated(meshSingle%dvol3f)," (should be F)"
  if (associated(meshSingle%dvol3f)) print*,meshSingle%dvol3f(1)
  print'(a,1pg13.6E2,a)',"h32f(3) = ",meshSingle%h32f(3)," (should be ?)"
  print'(a,1pg13.6E2,a)',"h32c(3) = ",meshSingle%h32c(3)," (should be ?)"
  deallocate (meshSingle)
  print*,''

  !-----------------------------------------------------------------------------
  ! Test collections of meshes
  !-----------------------------------------------------------------------------
  print*,repeat("=",70)
  print*,"**Sets of meshes**"
  print*,repeat("-",20)

  ! -------- Cartesian; all metric factors should be 1.0 --------
  call MeshFactory(mesh_types%Cartesian, meshCart, s=[1.0d0,1.0d0,1.0d0], &
                      p=[0.0d0,0.0d0,0.0d0])
  call meshCart(1)%print()
  call meshCart(2)%print()
  call meshCart(3)%print()

  print*,"is h31c associated in each mesh direction? the result should be: T,F,F"
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
  print'(3(1pg12.4E2),a)',meshCart%d," (should be 6.6667e-2 x3)"
  print'(8(1pg11.3))',meshCart(1)%h
  print*,"(the first element should be -0.5; the rest should be 0)"
  print*,''

  block
  class(mesh_t), pointer:: m1
  m1 => meshCart(1)
  ! here's another place where gfortran (v5.4.1 and earlier) chokes.
#ifdef __GFORTRAN__
  print*,"Columns 2, 3 and 4 should be identical:"
  associate (m1r => meshCart(1)%r)
  do i=meshCart(1)%lb,meshCart(1)%ub
    print*,i,meshCart(1)%r(i),m1r(i),m1%r(i)
  end do
  end associate
#else
  print*,"Columns 2 - 5 should be identical:"
  associate (m1r => meshCart(1)%r, m11=> meshCart(1))
  do i=meshCart(1)%lb,meshCart(1)%ub
    print*,i,meshCart(1)%r(i),m1r(i),m1%r(i),m11%r(i)
  end do
  end associate
#endif
  nullify (m1)
  end block
  print*,repeat("-",20)

  ! -------- cylindrical; the 1-direction is collapsed in this example --------
  call MeshFactory(mesh_types%cylindrical, meshCyl, n=[1,11,9], &
                      s=[1.0d0,0.5d0,2.0d0*pi], p=[0.0d0,0.0d0,0.0d0], &
                      llc_native=[0.0d0,0.1d0,0.5d0*pi],ng=[0,3,3])
  call meshCyl(1)%print()
  call meshCyl(2)%print()
  call meshCyl(3)%print()

  print*,"The 3rd column should be the absolute value of the 2nd;"
  print*,"the 4th column should be delta(r^2)*0.5, where r is the 2nd column."
  select type (meshCyl)
  type is (cylindrical_mesh)
    associate(cmesh => meshCyl(2))
    if (associated(cmesh%h32c)) then
      do i=cmesh%lo,cmesh%uo
        print'(i3.1,1x,3(1pg14.6E2))',i,cmesh%r(i)-cmesh%d*0.5+cmesh%centre_nat,cmesh%h32f(i),cmesh%dvol2c(i)
      end do
    endif
    end associate
  end select

  select type (meshCyl)
  type is (cylindrical_mesh)
    if (associated(meshCyl(1)%dvol1c)) &
      print'(2(1pg12.4E2),a)',meshCyl(1)%d,meshCyl(1)%dvol1c(1),"(both should be 1.0)"
  end select

  select type (meshCyl)
  type is (cylindrical_mesh)
    associate (cmesh => meshCyl(3))
    print*,associated(cmesh%h32f)," (should be F)"
    print*,associated(cmesh%dvol3c)," (should be T)"
    print'(1pg14.8,a)',cmesh%zeta," (should be 0.5 * pi = 1.570796...)"
    end associate
  end select
  print*,repeat("-",20)

  ! ---- spherical ----
  nTheta = 10
  call MeshFactory(mesh_types%spherical, meshSph, n=[5,nTheta,8], &
                      s=[1.0d0,pi,2.0d0*pi], p=[0.0d0,0.0d0,0.0d0], llc_native=[1.0d0,0.0d0,0.0d0])
  call meshSph(1)%print()
  call meshSph(2)%print()
  call meshSph(3)%print()

  select type (meshSph)
  type is (spherical_mesh)
    associate (cmesh2 => meshSph(2))
    print'(a,1pg12.6,a)','xi: ',cmesh2%xi," (should be 0 and equal to the value immediately above)"
    print*,'zeta: ',allocated(cmesh2%zeta)," (should be F)"
    print*,'h2c:  ',associated(cmesh2%h2c)," (should be F)"
    print*,"Column 2 should be the absolute value of column 3:"
    do i=1,nTheta
      print'(i4.1,1x,2(1pe14.6E2))',i,cmesh2%h32f(i),sin(cmesh2%centre_nat + cmesh2%r(i) - 0.5*cmesh2%d)
    end do
    end associate
  end select

  select type (meshSph)
  type is (spherical_mesh)
    associate(cmesh3 => meshSph(3))
    print'(a,1pg12.6,a)','zeta: ',cmesh3%zeta," (should be 0)"
    print*,'xi: ',allocated(cmesh3%xi)," (should be F)"
    end associate
  end select

  block
  class(mesh_t), pointer:: m1
  m1 => meshSph(1)
  ! here's another place where gfortran (v5.4.1 and earlier) chokes.
#ifdef __GFORTRAN__
  print*,"Columns 2, 3 and 4 should be identical:"
  associate (m1r => meshSph(1)%r)
  do i=meshSph(1)%lb,meshSph(1)%ub
    print*,i,meshSph(1)%r(i),m1r(i),m1%r(i)
  end do
  end associate
#else
  print*,"Columns 2 - 5 should be identical:"
  associate (m1r => meshSph(1)%r, m11=> meshSph(1))
  do i=meshSph(1)%lb,meshSph(1)%ub
    print*,i,meshSph(1)%r(i),m1r(i),m1%r(i),m11%r(i)
  end do
  end associate
#endif

  print*,"Columns 2, 3 and 4 should be identical:"
  select type (meshSph)
  type is (spherical_mesh)
    associate (m1r => meshSph(1)%r)
    do i=meshSph(1)%lb,meshSph(1)%ub
      print*,i,meshSph(1)%r(i),m1r(i),m1%r(i)
    end do
    end associate
  end select

  nullify (m1)
  end block

  call MeshRecycler(meshCart)
  call MeshRecycler(meshCyl)
  call MeshRecycler(meshSph)
  print*,repeat("-",20)
  print*,repeat("-",20)

  ! -------- Cartesian 2; all metric factors should be 1.0 --------
  allocate(meshCart2(3))
  do i=1,3
    meshCart2(i)%idir = i
    call meshCart2(i)%init(5,2,1.0d0,0.0d0,0.0d0)
  end do
  call meshCart2(1)%print()
  call meshCart2(2)%print()
  call meshCart2(3)%print()

  block
  class(mesh_t), pointer:: m1
  m1 => meshCart2(1)
  ! here's another place where gfortran (v5.4.1 and earlier) chokes.
#ifdef __GFORTRAN__
  print*,"Columns 2, 3 and 4 should be identical:"
  associate (m1r => meshCart2(1)%r)
  do i=meshCart2(1)%lb,meshCart2(1)%ub
    print*,i,meshCart2(1)%r(i),m1r(i),m1%r(i)
  end do
  end associate
#else
  print*,"Columns 2 - 5 should be identical:"
  associate (m1r => meshCart2(1)%r, m11=> meshCart2(1))
  do i=meshCart2(1)%lb,meshCart2(1)%ub
    print*,i,meshCart2(1)%r(i),m1r(i),m1%r(i),m11%r(i)
  end do
  end associate
#endif
  nullify (m1)
  end block
  deallocate(meshCart2)
  print*,repeat("-",20)

  ! -------- Cartesian 3; all metric factors should be 1.0 --------
  do i=1,3
    meshCart3(i)%idir = i
    call meshCart3(i)%init(6,2,1.0d0,0.5d0,0.0d0)
  end do
  call meshCart3(1)%print()
  call meshCart3(2)%print()
  call meshCart3(3)%print()

  block
  class(mesh_t), pointer:: m1
  m1 => meshCart3(1)
  ! here's another place where gfortran (v5.4.1 and earlier) chokes.
#ifdef __GFORTRAN__
  print*,"Columns 2, 3 and 4 should be identical:"
  associate (m1r => meshCart3(1)%r)
  do i=meshCart3(1)%lb,meshCart3(1)%ub
    print*,i,meshCart3(1)%r(i),m1r(i),m1%r(i)
  end do
  end associate
#else
  print*,"Columns 2 - 5 should be identical:"
  associate (m1r => meshCart3(1)%r, m11=> meshCart3(1))
  do i=meshCart3(1)%lb,meshCart3(1)%ub
    print*,i,meshCart3(1)%r(i),m1r(i),m1%r(i),m11%r(i)
  end do
  end associate
#endif
  nullify (m1)
  end block

  !-----------------------------------------------------------------------------
  ! Test integration with patches
  !-----------------------------------------------------------------------------
  print*,repeat("=",70)
  print*,"**Test integration with patches**"
  print*,repeat("-",20)

  ! ---- Take 1: Cartesian ----
  allocate(patch_t :: patch)
  call patch%idx%init (ienergy=2)
  call MeshFactory(mesh_types%cartesian, patch%mesh, n=[11,11,1], ng=[2,2,0])
  call patch%mesh(1)%print()
  call patch%mesh(2)%print()
  call patch%mesh(3)%print()

  ! The following works in ifort but not in gfortran (v5.4.1 and earlier).
#ifndef __GFORTRAN__
  associate (m => patch%mesh)
  print*,shape(m),'(should be 3)'
  end associate
#endif

  print'(3(1pg10.4),a)',patch%mesh(1)%h2c(3),patch%mesh(2)%h32c(3),patch%mesh(3)%dvol3c(1),&
    " (all should be 1.0)"

  call MeshRecycler(patch%mesh)
  deallocate(patch)
  print*,repeat("-",20)

  ! ---- Take 2: cylindrical ----
  allocate(patch_t :: patch)
  patch%n = [3,4,5]
  patch%size = [2.0,3.0,2.0*real(pi,kind=4)]
  patch%position = [0.25,0.30,0.35] ! in Cartesian coords.
  call MeshFactory(mesh_types%cylindrical, patch%mesh, n=patch%n, s=patch%size, &
                   p=patch%position)
  call patch%mesh(1)%print()
  call patch%mesh(2)%print()
  call patch%mesh(3)%print()

  print*,lbound(patch%mesh),ubound(patch%mesh)," (should be 1 and 3)"
  print*,"The following two lines should be identical (and equal to 1 and 9):"
  l = patch%mesh%lb
  u = patch%mesh%ub
  print*,patch%mesh(3)%lb,patch%mesh(3)%ub
  print*,l(3),u(3)

  block
  class(mesh_t), pointer :: m3
  m => patch%mesh
  m3 => patch%mesh(3)
  print*,"Columns 2, 3, 4 and 5 should be identical:"
  associate (m3r => patch%mesh(3)%r)
  do i=m3%lb,m3%ub
    print'(i4.1,4(1x,1pe14.6))',i,m(3)%r(i),m3%r(i),patch%mesh(3)%r(i),m3r(i)
  end do
  end associate

  print*,''
  print*,"Tests using `select type`:"
  ! NOTE!: the following does not work (it causes an ICE) in gfortran v5.4.1 and earlier.
#ifndef __GFORTRAN__
  associate (m2 => patch%mesh(2))
  select type(m2)
  type is (Cartesian_mesh)
    print*,'m2: Cartesian mesh recognised'
  type is (cylindrical_mesh)
    print*,'m2: Cylindrical mesh recognised'
  type is (spherical_mesh)
    print*,'m2: Spherical mesh recognised'
  class default
    print*,'m2: Unknown mesh type!'
  end select
  end associate
#endif

  ! NOTE!: The following does not seem to work with IFORT 17.0.2!!!
  print*,"The following is known to not work in IFORT v17.0.2, but it does work in v17.0.1!!"
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

  ! But this does work in ifort (but causes an ICE in gfortran v5.4.1)!
#ifndef __GFORTRAN__
  print*,"Workaround:"
  associate(m33 => m(3))
  select type(m33)
  type is (Cartesian_mesh)
    print*,'m33: Cartesian mesh recognised'
  type is (cylindrical_mesh)
    print*,'m33: Cylindrical mesh recognised'
  type is (spherical_mesh)
    print*,'m33: Spherical mesh recognised'
  class default
    print*,'m33: Unknown mesh type!'
  end select
  end associate
#endif

  print*,"Columns 2, 3, 4 and 5 should be identical:"
  select type(m)
  type is (cylindrical_mesh)
    associate (m3r => patch%mesh(3)%r)
    do i=m3%lb,m3%ub
      print'(i4.1,4(1x,1pe14.6))',i,m(3)%r(i),m3%r(i),patch%mesh(3)%r(i),m3r(i)
    end do
    end associate
  end select

  nullify (m, m3)
  end block
  deallocate(patch)
  print*,repeat("-",20)
  print*,"Finished."

END PROGRAM
