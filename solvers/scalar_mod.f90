!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE scalar_mod
  USE io_mod
  USE trace_mod
  implicit none
  public
  !
  integer, private:: verbose=0
  type, private:: node_t
    type(node_t), pointer:: next=>null()
  end type
  type(node_t), private, pointer:: head=>null(), tail=>null()
  integer, private, save:: id=0
  !
  integer, private:: nalloc=0
  interface log
    module procedure scalar_log_scalar
  end interface
  interface exp
    module procedure scalar_exp_scalar
  end interface
CONTAINS

FUNCTION scalar_log_scalar(in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  integer:: i
  do i=1,size(in,3)
    out(:,:,i) = log(in(:,:,i))
  end do
END FUNCTION

FUNCTION scalar_exp_scalar(in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  integer:: i
  do i=1,size(in,3)
    out(:,:,i) = exp(in(:,:,i))
  end do
END FUNCTION

!===============================================================================
SUBROUTINE allocate_scalars_a (n, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)
  integer:: n(3)
  real, dimension(:,:,:), allocatable, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
  !.............................................................................
  call trace_begin('allocate_scalars_a')
  if (present(v1)) allocate(v1(n(1),n(2),n(3)))
  if (present(v2)) allocate(v2(n(1),n(2),n(3)))
  if (present(v3)) allocate(v3(n(1),n(2),n(3)))
  if (present(v4)) allocate(v4(n(1),n(2),n(3)))
  if (present(v5)) allocate(v5(n(1),n(2),n(3)))
  if (present(v6)) allocate(v6(n(1),n(2),n(3)))
  if (present(v7)) allocate(v7(n(1),n(2),n(3)))
  if (present(v8)) allocate(v8(n(1),n(2),n(3)))
  if (present(v9)) allocate(v9(n(1),n(2),n(3)))
  if (present(v10)) allocate(v10(n(1),n(2),n(3)))
  if (present(v1)) call io%bits_mem(storage_size(v1),product(shape(v1)),'v1')
  if (present(v2)) call io%bits_mem(storage_size(v2),product(shape(v2)),'v2')
  if (present(v3)) call io%bits_mem(storage_size(v3),product(shape(v3)),'v3')
  if (present(v4)) call io%bits_mem(storage_size(v4),product(shape(v4)),'v4')
  if (present(v5)) call io%bits_mem(storage_size(v5),product(shape(v5)),'v5')
  if (present(v6)) call io%bits_mem(storage_size(v6),product(shape(v6)),'v6')
  if (present(v7)) call io%bits_mem(storage_size(v7),product(shape(v7)),'v7')
  if (present(v8)) call io%bits_mem(storage_size(v8),product(shape(v8)),'v8')
  if (present(v9)) call io%bits_mem(storage_size(v9),product(shape(v9)),'v9')
  if (present(v10)) call io%bits_mem(storage_size(v10),product(shape(v10)),'v10')
  call trace_end
END SUBROUTINE

!===============================================================================
SUBROUTINE allocate_scalars (n, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)
  integer:: n(3)
  real, dimension(:,:,:), pointer, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
  !.............................................................................
  call trace_begin('allocate_scalars')
  if (present(v1)) allocate(v1(n(1),n(2),n(3)))
  if (present(v2)) allocate(v2(n(1),n(2),n(3)))
  if (present(v3)) allocate(v3(n(1),n(2),n(3)))
  if (present(v4)) allocate(v4(n(1),n(2),n(3)))
  if (present(v5)) allocate(v5(n(1),n(2),n(3)))
  if (present(v6)) allocate(v6(n(1),n(2),n(3)))
  if (present(v7)) allocate(v7(n(1),n(2),n(3)))
  if (present(v8)) allocate(v8(n(1),n(2),n(3)))
  if (present(v9)) allocate(v9(n(1),n(2),n(3)))
  if (present(v10)) allocate(v10(n(1),n(2),n(3)))
  if (present(v1)) call io%bits_mem(storage_size(v1),product(shape(v1)),'v1')
  if (present(v2)) call io%bits_mem(storage_size(v2),product(shape(v2)),'v2')
  if (present(v3)) call io%bits_mem(storage_size(v3),product(shape(v3)),'v3')
  if (present(v4)) call io%bits_mem(storage_size(v4),product(shape(v4)),'v4')
  if (present(v5)) call io%bits_mem(storage_size(v5),product(shape(v5)),'v5')
  if (present(v6)) call io%bits_mem(storage_size(v6),product(shape(v6)),'v6')
  if (present(v7)) call io%bits_mem(storage_size(v7),product(shape(v7)),'v7')
  if (present(v8)) call io%bits_mem(storage_size(v8),product(shape(v8)),'v8')
  if (present(v9)) call io%bits_mem(storage_size(v9),product(shape(v9)),'v9')
  if (present(v10)) call io%bits_mem(storage_size(v10),product(shape(v10)),'v10')
  call trace_end
END SUBROUTINE

!===============================================================================
SUBROUTINE deallocate_scalars_a (v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)
  integer:: n(3)
  real, dimension(:,:,:), allocatable, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
  !.............................................................................
  call trace_begin('deallocate_scalars_a')
  if (present(v1)) call io%bits_mem(storage_size(v1),-product(shape(v1)),'v1')
  if (present(v2)) call io%bits_mem(storage_size(v2),-product(shape(v2)),'v2')
  if (present(v3)) call io%bits_mem(storage_size(v3),-product(shape(v3)),'v3')
  if (present(v4)) call io%bits_mem(storage_size(v4),-product(shape(v4)),'v4')
  if (present(v5)) call io%bits_mem(storage_size(v5),-product(shape(v5)),'v5')
  if (present(v6)) call io%bits_mem(storage_size(v6),-product(shape(v6)),'v6')
  if (present(v7)) call io%bits_mem(storage_size(v7),-product(shape(v7)),'v7')
  if (present(v8)) call io%bits_mem(storage_size(v8),-product(shape(v8)),'v8')
  if (present(v9)) call io%bits_mem(storage_size(v9),-product(shape(v9)),'v9')
  if (present(v10)) call io%bits_mem(storage_size(v10),-product(shape(v10)),'v10')
  if (present(v1)) deallocate(v1)
  if (present(v2)) deallocate(v2)
  if (present(v3)) deallocate(v3)
  if (present(v4)) deallocate(v4)
  if (present(v5)) deallocate(v5)
  if (present(v6)) deallocate(v6)
  if (present(v7)) deallocate(v7)
  if (present(v8)) deallocate(v8)
  if (present(v9)) deallocate(v9)
  if (present(v10)) deallocate(v10)
  call trace_end
END SUBROUTINE

!===============================================================================
SUBROUTINE deallocate_scalars (v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)
  integer:: n(3)
  real, dimension(:,:,:), pointer, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
  !.............................................................................
  call trace_begin('deallocate_scalars')
  if (present(v1)) call io%bits_mem(storage_size(v1),-product(shape(v1)),'v1')
  if (present(v2)) call io%bits_mem(storage_size(v2),-product(shape(v2)),'v2')
  if (present(v3)) call io%bits_mem(storage_size(v3),-product(shape(v3)),'v3')
  if (present(v4)) call io%bits_mem(storage_size(v4),-product(shape(v4)),'v4')
  if (present(v5)) call io%bits_mem(storage_size(v5),-product(shape(v5)),'v5')
  if (present(v6)) call io%bits_mem(storage_size(v6),-product(shape(v6)),'v6')
  if (present(v7)) call io%bits_mem(storage_size(v7),-product(shape(v7)),'v7')
  if (present(v8)) call io%bits_mem(storage_size(v8),-product(shape(v8)),'v8')
  if (present(v9)) call io%bits_mem(storage_size(v9),-product(shape(v9)),'v9')
  if (present(v10)) call io%bits_mem(storage_size(v10),-product(shape(v10)),'v10')
  if (present(v1)) deallocate(v1)
  if (present(v2)) deallocate(v2)
  if (present(v3)) deallocate(v3)
  if (present(v4)) deallocate(v4)
  if (present(v5)) deallocate(v5)
  if (present(v6)) deallocate(v6)
  if (present(v7)) deallocate(v7)
  if (present(v8)) deallocate(v8)
  if (present(v9)) deallocate(v9)
  if (present(v10)) deallocate(v10)
  call trace_end
END SUBROUTINE

!===============================================================================
SUBROUTINE scalar_stats (f, label)
  real, dimension(:,:,:):: f
  character(len=*):: label
  integer:: ix, iy, iz
  real(8):: s0, s1, s2, s3
  s1 = 0d0
  s2 = 0d0
  s0 = f(1,1,1)
  s3 = f(1,1,1)
  do iz=1,size(f,3)
  do iy=1,size(f,2)
  do ix=1,size(f,1)
    s1 = s1 + f(ix,iy,iz)
    s2 = s2 + f(ix,iy,iz)**2
    s0 = min(s0,f(ix,iy,iz))
    s3 = max(s3,f(ix,iy,iz))
  end do
  end do
  end do
  s1 = s1/size(f)
  s2 = sqrt(s2/size(f)-s1**2)
  print*,'scalar_stats:', trim(label), ': min, aver, rms, max =', s0, s1, s2, s3
END SUBROUTINE

!===============================================================================
END MODULE scalar_mod
