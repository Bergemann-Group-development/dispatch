!===============================================================================
MODULE operators_mod
  USE scalar_mod
  implicit none
CONTAINS
!===============================================================================
FUNCTION crossf (f, g) result (h)
  real, dimension(:,:,:,:), intent(in):: f, g
  real, dimension(size(f,1),size(f,2),size(f,3),size(f,4)):: h
  h(:,:,:,1) = f(:,:,:,2)*g(:,:,:,3)-f(:,:,:,3)*g(:,:,:,2)
  h(:,:,:,2) = f(:,:,:,3)*g(:,:,:,1)-f(:,:,:,1)*g(:,:,:,3)
  h(:,:,:,3) = f(:,:,:,1)*g(:,:,:,2)-f(:,:,:,2)*g(:,:,:,1)
END FUNCTION
!===============================================================================
SUBROUTINE cross1 (f, g, h)
  real, dimension(:,:,:,:), pointer:: f, g, h
  integer:: i1, i2, i3
  do i1=1,3
    i2=mod(i1,3)+1
    i3=mod(i2,3)+1
    h(:,:,:,i1) = f(:,:,:,i2)*g(:,:,:,i3)-f(:,:,:,i3)*g(:,:,:,i2)
  end do
END SUBROUTINE cross1
!===============================================================================
SUBROUTINE cross2 (f, g, h)
  real, dimension(:,:,:,:), pointer:: f, g, h
  integer:: i1,i2,i3
  do i1=1,3
    i2=mod(i1,3)+1
    i3=mod(i2,3)+1
    h(:,:,:,i1) = f(:,:,:,i2)*g(:,:,:,i3)-f(:,:,:,i3)*g(:,:,:,i2)
  end do
END SUBROUTINE cross2
!===============================================================================
END MODULE

!===============================================================================
PROGRAM test
  USE scalar_mod
  USE vector_mod
  USE mesh_mod
  USE io_mod
  USE omp_timer_mod
  USE operators_mod
  USE stagger_mod
  implicit none
  type(scalar_t):: a, b
  type(vector_t):: u, v
  type(mesh_t):: m
  real, dimension(:,:,:,:), pointer:: f, g, h, h0
  type(vector_t):: v1, v2, v3, v4, v5, v6, v7, v8, v9
  type(scalar_t):: s1, s2
  type(stagger_t):: d
  integer:: i, n, dims(3)
  real:: ds(3)
  namelist/test_params/ n
  !-----------------------------------------------------------------------------
  call io%init
  n = 20
  ds = 1.0
  read (io%input, test_params)
  write (*,test_params)
  dims = n
  call tic
  call m%init (n=dims,ng=[2,2,2])
  call timer%toc("init",n**3)
  allocate (f(n,n,n,3), g(n,n,n,3), h(n,n,n,3), h0(n,n,n,3))
  call timer%toc("allocate",n**3)
  do i=1,3
    f(:,:,:,i) = 1.0+i
    g(:,:,:,i) = 2.0-i
  end do
  call timer%toc("set",n**3)
  call cross1 (f, g, h)
  call timer%toc("cross1",n**3)
  call tic
  h0 = h
  call timer%toc("copy",n**3)
  call tic
  call cross2 (f, g, h)
  call timer%toc("cross2",n**3)
  call tic
  call cross1 (f, g, h)
  call timer%toc("cross1",n**3)
  call tic
  h =  crossf (f, g)
  call timer%toc("crossf",n**3)
  print*,minval(h-h0),maxval(h-h0)
  call tic
  h =  cross (f, g)
  call timer%toc("cross",n**3)
  print*,minval(h-h0),maxval(h-h0)
  call tic
  h =  crossf (f, g)
  call timer%toc("crossf",n**3)
  print*,minval(h-h0),maxval(h-h0)
  !-----------------------------------------------------------------------------
  call v1%init(dims,ds,'v1')
  call v2%init(dims,ds,'v2')
  print*,'mk1',v2%allocated,' ',trim(v2%name),v2%x%n,v2%x%l
  v1 = 1.0
  v2 = 1.0
  call tic
  v2 = curl(v1)
  call timer%toc("curl",n**3)
  print*,'mk2',v2%allocated,' ',trim(v2%name),v2%x%n,v2%x%l
  !-----------------------------------------------------------------------------
  !call init_vectors (dims, ds, v1, v2, v3, v4, v5, v6, v7, v8, v9)
  call init_vectors (dims, ds, v1, v2, v3)
  call s1%init(dims,ds,'s1')
  call tic
  s1 = div(v3)
  print*,'mk3',v2%allocated,' ',trim(v2%name),v2%x%n,v2%x%l
  print*,'mk3',v3%allocated,' ',trim(v3%name),v3%x%n,v3%x%l
  print*,'mk3',s1%allocated,' ',trim(s1%name),s1%n,s1%l
  call timer%toc("div",n**3)
  v3 = v2
  !-----------------------------------------------------------------------------
  call tic
  v1 = grad(s1)
  call timer%toc("grad",n**3)
  !-----------------------------------------------------------------------------
  call tic
  v1 = down(s1)
  call timer%toc("down",n**3)
  !-----------------------------------------------------------------------------
  v1 = 1.0
  call tic
  v1 = exp(v1)
  call timer%toc("exp",n**3)
  !-----------------------------------------------------------------------------
  call tic
  v1 = log(v1)
  call timer%toc("log",n**3)
  !-----------------------------------------------------------------------------
  call tic
  v1 = edge(s1)
  call timer%toc("edge",n**3)
  !-----------------------------------------------------------------------------
  block
  type(scalar_t):: d, s, dddt, dsdt
  type(vector_t):: p, B, dpdt, dBdt
  call init_scalars (dims, ds, d, s, dddt, dsdt)
  call init_vectors (dims, ds, p, B, dpdt, dBdt)
  call tic
  do i=1,10
    call pde (d, dddt, s, dsdt, p, dpdt, B, dBdt)
  end do
  call timer%toc("mhd",10*n**3)
  call tic
  do i=1,10
    call pde (d, dddt, s, dsdt, p, dpdt)
  end do
  call timer%toc("hd",10*n**3)
  do i=1,10
    call pde (d, dddt, s, dsdt, p, dpdt, B, dBdt)
  end do
  call timer%toc("mhd",10*n**3)
  end block
  !-----------------------------------------------------------------------------
  call cleanup_temp
END PROGRAM
