MODULE stagger_2nd
CONTAINS
!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  if (size(a,1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=2,size(a,1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
    end do
    end do
    end do
  else
    b = a
  end if
END FUNCTION xdn
!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  if (size(a,2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=1,size(a,3)
    do iy=2,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
    end do
    end do
    end do
  else
    b = a
  end if
END FUNCTION ydn
!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  if (size(a,3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=2,size(a,3)
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
    end do
    end do
    end do
  else
    b = a
  end if
END FUNCTION zdn
END MODULE stagger_2nd

MODULE stagger_4th
CONTAINS
!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz
!...............................................................................
  c1 = 0.5
  c2 = 0.5
  if (size(a,1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=3,size(a,1)-1
      b(ix,iy,iz) = c1*(a(ix  ,iy,iz) + a(ix-1,iy,iz)) &
                  + c2*(a(ix+1,iy,iz) + a(ix-2,iy,iz))
    end do
    end do
    end do
  else
    b = a
  end if
END FUNCTION xdn
!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz
!...............................................................................
  c1 = 0.5
  c2 = 0.5
  if (size(a,2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=1,size(a,3)
    do iy=3,size(a,2)-1
    do ix=1,size(a,1)
      b(ix,iy,iz) = c1*(a(ix,iy  ,iz) + a(ix,iy-1,iz)) &
                  + c2*(a(ix,iy+1,iz) + a(ix,iy-2,iz))
    end do
    end do
    end do
  else
    b = a
  end if
END FUNCTION ydn
!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz
!...............................................................................
  c1 = 0.5
  c2 = 0.5
  if (size(a,3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=3,size(a,3)-1
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c1*(a(ix,iy,iz  ) + a(ix,iy,iz-1)) &
                  + c2*(a(ix,iy,iz+1) + a(ix,iy,iz-2))
    end do
    end do
    end do
  else
    b = a
  end if
END FUNCTION zdn
END MODULE stagger_4th

MODULE stagger_6th
CONTAINS
!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2, c3
  integer :: ix, iy, iz
!...............................................................................
  c1 = 0.5
  c2 = 0.2
  c3 = 0.1
  if (size(a,1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=4,size(a,1)-2
      b(ix,iy,iz) = c1*(a(ix  ,iy,iz) + a(ix-1,iy,iz)) &
                  + c2*(a(ix+1,iy,iz) + a(ix-2,iy,iz)) &
                  + c3*(a(ix+2,iy,iz) + a(ix-3,iy,iz))
    end do
    end do
    end do
  else
    b = a
  end if
END FUNCTION xdn
!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2, c3
  integer :: ix, iy, iz
!...............................................................................
  c1 = 0.5
  c2 = 0.2
  c3 = 0.1
  if (size(a,2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iz=1,size(a,3)
    do iy=4,size(a,2)-2
    do ix=1,size(a,1)
      b(ix,iy,iz) = c1*(a(ix,iy  ,iz) + a(ix,iy-1,iz)) &
                  + c2*(a(ix,iy+1,iz) + a(ix,iy-2,iz)) &
                  + c3*(a(ix,iy+2,iz) + a(ix,iy-3,iz))
    end do
    end do
    end do
  else
    b = a
  end if
END FUNCTION ydn
!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2, c3
  integer :: ix, iy, iz
!...............................................................................
  c1 = 0.5
  c2 = 0.2
  c3 = 0.1
  if (size(a,3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    do iy=1,size(a,2)
    do iz=4,size(a,3)-2
    do ix=1,size(a,1)
      b(ix,iy,iz) = c1*(a(ix,iy,iz  ) + a(ix,iy,iz-1)) &
                  + c2*(a(ix,iy,iz+1) + a(ix,iy,iz-2)) &
                  + c3*(a(ix,iy,iz+2) + a(ix,iy,iz-3))
    end do
    end do
    end do
  else
    b = a
  end if
END FUNCTION zdn
END MODULE stagger_6th

!===============================================================================
!> Test program for MPI parallel I/O.  gsize3 is the global size of the array,
!> lsize3 is the rank-thread-local size of the array, and start3 is the offset
!> into the global array.  Each thread handles a separate array f, with an 
!> offset (in this test) that increases with MPI rank in the y-direction and 
!> with OMP thread in the z-direction
!===============================================================================
PROGRAM test_stagger
  implicit none
  integer, parameter:: n=36
  real, pointer, dimension(:,:,:):: f, d
  !-----------------------------------------------------------------------------
  allocate (f(n,n,n), d(n,n,n))
  f=1.
  print*,'n =',n
  print *,'2nd order:'
  call test_2nd (n, f, d)
  print *,'4th order:'
  call test_4th (n, f, d)
  print *,'6th order:'
  call test_6th (n, f, d)
END PROGRAM

SUBROUTINE test_2nd (nn, f, d)
  USE stagger_2nd
  USE timer_mod
  USE omp_timer_mod
  implicit none
  integer(8), parameter:: n=36
  real, dimension(n,n,n):: f, d
  integer:: nn, m
  real(8):: wt
  m=0
  wt = wallclock()
  call tic
  do while ((wallclock()-wt) < 1.0)
    d=xdn(f)
    m=m+1
  end do
  call timer%toc('xdn',m*n**3)
  m=0
  wt = wallclock()
  call tic
  do while ((wallclock()-wt) < 1.0)
    d=ydn(f)
    m=m+1
  end do
  call timer%toc('ydn',m*n**3)
  m=0
  wt = wallclock()
  call tic
  do while ((wallclock()-wt) < 1.0)
    d=zdn(f)
    m=m+1
  end do
  call timer%toc('zdn',m*n**3)
END SUBROUTINE

SUBROUTINE test_4th (nn, f, d)
  USE stagger_4th
  USE timer_mod
  USE omp_timer_mod
  implicit none
  integer(8), parameter:: n=36
  real, dimension(n,n,n):: f, d
  integer:: nn, m
  real(8):: wt
  m=0
  wt = wallclock()
  call tic
  do while ((wallclock()-wt) < 1.0)
    d=xdn(f)
    m=m+1
  end do
  call timer%toc('xdn',m*n**3)
  m=0
  wt = wallclock()
  call tic
  do while ((wallclock()-wt) < 1.0)
    d=ydn(f)
    m=m+1
  end do
  call timer%toc('ydn',m*n**3)
  m=0
  wt = wallclock()
  call tic
  do while ((wallclock()-wt) < 1.0)
    d=zdn(f)
    m=m+1
  end do
  call timer%toc('zdn',m*n**3)
END SUBROUTINE

SUBROUTINE test_6th (nn, f, d)
  USE stagger_6th
  USE timer_mod
  USE omp_timer_mod
  implicit none
  integer(8), parameter:: n=36
  real, dimension(n,n,n):: f, d
  integer:: nn, m
  real(8):: wt
  m=0
  wt = wallclock()
  d=0.0
  call tic
  do while ((wallclock()-wt) < 1.0)
    d=d+xdn(f+m)
    m=m+1
  end do
  call timer%toc('xdn',m*n**3)
  m=0
  wt = wallclock()
  d=0.0
  call tic
  do while ((wallclock()-wt) < 1.0)
    d=d+ydn(f+m)
    m=m+1
  end do
  call timer%toc('ydn',m*n**3)
  m=0
  wt = wallclock()
  d=0.0
  call tic
  do while ((wallclock()-wt) < 1.0)
    d=d+zdn(f+m)
    m=m+1
  end do
  call timer%toc('zdn',m*n**3)
END SUBROUTINE

