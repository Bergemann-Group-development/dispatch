!> =============================================================================
MODULE operators_mod
  implicit none
CONTAINS
!> =============================================================================
FUNCTION crossf (f, g) result (h)
  real, dimension(:,:,:,:), intent(in):: f, g
  real, dimension(size(f,1),size(f,2),size(f,3),size(f,4)):: h
  h(:,:,:,1) = f(:,:,:,2)*g(:,:,:,3)-f(:,:,:,3)*g(:,:,:,2)
  h(:,:,:,2) = f(:,:,:,3)*g(:,:,:,1)-f(:,:,:,1)*g(:,:,:,3)
  h(:,:,:,3) = f(:,:,:,1)*g(:,:,:,2)-f(:,:,:,2)*g(:,:,:,1)
END FUNCTION
!> =============================================================================
SUBROUTINE cross1 (f, g, h)
  real, dimension(:,:,:,:), pointer:: f, g, h
  h(:,:,:,1) = f(:,:,:,2)*g(:,:,:,3)-f(:,:,:,3)*g(:,:,:,2)
  h(:,:,:,2) = f(:,:,:,3)*g(:,:,:,1)-f(:,:,:,1)*g(:,:,:,3)
  h(:,:,:,3) = f(:,:,:,1)*g(:,:,:,2)-f(:,:,:,2)*g(:,:,:,1)
END SUBROUTINE cross1
!> =============================================================================
SUBROUTINE cross2 (f, g, h)
  real, dimension(:,:,:,:), pointer:: f, g, h
  integer:: i1,i2,i3
  do i1=1,3
    i2=mod(i1,3)+1
    i3=mod(i2,3)+1
    h(:,:,:,i1) = f(:,:,:,i2)*g(:,:,:,i3)-f(:,:,:,i3)*g(:,:,:,i2)
  end do
END SUBROUTINE cross2
!> =============================================================================
END MODULE

!> =============================================================================
PROGRAM test
  USE scalar_mod
  USE vector_mod
  USE mesh_mod
  USE io_mod
  USE omp_timer_mod
  USE operators_mod
  implicit none
  type(scalar_t):: a, b
  type(vector_t):: u, v
  type(mesh_t):: m
  real, dimension(:,:,:,:), pointer:: f, g, h
  integer:: n
  namelist/test_params/ n
  !-----------------------------------------------------------------------------
  call io%init
  n = 20
  read (io%input, test_params)
  write (*,test_params)
  call tic
  call m%init (n=[n,n,n],ng=[2,2,2])
  call toc("init",n**3)
  allocate (f(n,n,n,3), g(n,n,n,3), h(n,n,n,3))
  call toc("allocate",3*n**3)
  f = 1.0
  g = 2.0
  call toc("set",3*n**3)
  call cross1 (f, g, h)
  call toc("cross1",3*n**3)
  call cross2 (f, g, h)
  call toc("cross2",3*n**3)
  call cross1 (f, g, h)
  f = h
  call toc("cross1",3*n**3)
  h =  crossf (f, g)
  call toc("crossf",3*n**3)
END PROGRAM
