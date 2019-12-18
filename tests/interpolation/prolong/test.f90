!===============================================================================
!===============================================================================
MODULE prolong_mod
  USE omp_lib
  implicit none
  private
  type, public:: prolong_t
    contains
    procedure, nopass:: x_interp1
    procedure, nopass:: x_interp3
    generic:: x_interp => x_interp1, x_interp3
    procedure, nopass:: y_interp1
    procedure, nopass:: y_interp3
    generic:: y_interp => y_interp1, y_interp3
    procedure, nopass:: z_interp1
    procedure, nopass:: z_interp3
    generic:: z_interp => z_interp1, z_interp3
    procedure:: interp1
    procedure:: interp3
    generic:: interp => interp1, interp3
    procedure, nopass:: test
  end type
  type(prolong_t), public:: prolong
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE x_interp1 (in, w, l, u, out, lo)
  real, intent(in):: in(:,:,:), w(:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), lo
  !............................................................................
  integer:: i1, i2, i3, j1, n(3)
  real:: p, q
  !----------------------------------------------------------------------------
  n = shape(out)
  out = 0.0
  j1 = lo
  do i1=l(1),u(1)
    do i3=l(3),u(3)
      if (j1<1 .or. j1>n(1)) then
        cycle
      else if (j1==n(1)) then
        do i2=l(2),u(2)
          out(j1  ,i2,i3) = out(j1  ,i2,i3) + in(i1-1,i2,i3)*    w(i1-1) &
                                            + in(i1  ,i2,i3)*(1.-w(i1  ))
        end do
      else if (i1==1) then
        do i2=l(2),u(2)
          out(j1  ,i2,i3) = out(j1  ,i2,i3) + in(i1  ,i2,i3)
          out(j1+1,i2,i3) = out(j1+1,i2,i3) + in(i1  ,i2,i3)*    w(i1  ) &
                                            + in(i1+1,i2,i3)*(1.-w(i1+1))
        end do
      else if (i1==n(1)) then
        do i2=l(2),u(2)
          out(j1  ,i2,i3) = out(j1  ,i2,i3) + in(i1-1,i2,i3)*    w(i1-1) &
                                            + in(i1  ,i2,i3)*(1.-w(i1  ))
          out(j1+1,i2,i3) = out(j1+1,i2,i3) + in(i1  ,i2,i3)
        end do
      else
        do i2=l(2),u(2)
          out(j1  ,i2,i3) = out(j1  ,i2,i3) + in(i1-1,i2,i3)*    w(i1-1) &
                                            + in(i1  ,i2,i3)*(1.-w(i1  ))
          out(j1+1,i2,i3) = out(j1+1,i2,i3) + in(i1  ,i2,i3)*    w(i1  ) &
                                            + in(i1+1,i2,i3)*(1.-w(i1+1))
        end do
      end if
    end do
    j1 = j1 + 2
  end do
  !print '("X",5x,3(a,":",2i4,5x),f10.6)', &
    !'x', lo,j1-1, 'y', l(2),u(2), 'z',l(3),u(3), sum(out)/product(shape(out))
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE x_interp3 (in, w, l, u, out, lo)
  real, intent(in):: in(:,:,:), w(:,:,:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), lo
  !............................................................................
  integer:: i1, i2, i3, j1, n(3)
  real:: p, q
  !----------------------------------------------------------------------------
  n = shape(out)
  out = 0.0
  j1 = lo
  !print '(3(4x,a,":",3i4))', 'x_interp3 lt', l, 'ut', u, 'ls', lo
  do i1=l(1),u(1)
    do i3=l(3),u(3)
      if (j1<1 .or. j1>n(1)) then
        cycle
      else if (j1==n(1)) then
        do i2=l(2),u(2)
          out(j1  ,i2,i3) = out(j1  ,i2,i3) + in(i1-1,i2,i3)*    w(i1-1,i2,i3) &
                                            + in(i1  ,i2,i3)*(1.-w(i1  ,i2,i3))
        end do
      else if (i1==1) then
        do i2=l(2),u(2)
          out(j1  ,i2,i3) = out(j1  ,i2,i3) + in(i1  ,i2,i3)
          out(j1+1,i2,i3) = out(j1+1,i2,i3) + in(i1  ,i2,i3)*    w(i1  ,i2,i3) &
                                            + in(i1+1,i2,i3)*(1.-w(i1+1,i2,i3))
        end do
      else if (i1==n(1)) then
        do i2=l(2),u(2)
          out(j1  ,i2,i3) = out(j1  ,i2,i3) + in(i1-1,i2,i3)*    w(i1-1,i2,i3) &
                                            + in(i1  ,i2,i3)*(1.-w(i1  ,i2,i3))
          out(j1+1,i2,i3) = out(j1+1,i2,i3) + in(i1  ,i2,i3)
        end do
      else
        do i2=l(2),u(2)
          out(j1  ,i2,i3) = out(j1  ,i2,i3) + in(i1-1,i2,i3)*    w(i1-1,i2,i3) &
                                            + in(i1  ,i2,i3)*(1.-w(i1  ,i2,i3))
          out(j1+1,i2,i3) = out(j1+1,i2,i3) + in(i1  ,i2,i3)*    w(i1  ,i2,i3) &
                                            + in(i1+1,i2,i3)*(1.-w(i1+1,i2,i3))
        end do
      end if
    end do
    j1 = j1 + 2
  end do
  !print '("X",5x,3(a,":",2i4,5x),f10.6)', &
    !'x', lo,j1-1, 'y', l(2),u(2), 'z',l(3),u(3), sum(out)/product(shape(out))
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE y_interp1 (in, w, l, u, out, lo)
  real, intent(in):: in(:,:,:), w(:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), lo
  !............................................................................
  integer:: i1, i2, i3, j2, n(3)
  !----------------------------------------------------------------------------
  n = shape(out)
  out = 0.0
  j2 = lo
  do i2=l(2),u(2)
    do i3=l(3),u(3)
      if (j2<1 .or. j2>n(2)) then
        cycle
      else if (j2==n(2)) then
        do i1=l(1),u(1)
          out(i1,j2  ,i3) = out(i1,j2  ,i3) + in(i1,i2-1,i3)*    w(i2-1) &
                                            + in(i1  ,i2,i3)*(1.-w(i2  ))
        end do
      else if (i2==1) then
        do i1=l(1),u(1)
          out(i1,j2  ,i3) = out(i1,j2  ,i3) + in(i1  ,i2,i3)
          out(i1,j2+1,i3) = out(i1,j2+1,i3) + in(i1  ,i2,i3)*    w(i2  ) &
                                            + in(i1,i2+1,i3)*(1.-w(i2+1))
        end do
      else if (i2==n(2)) then
        do i1=l(1),u(1)
          out(i1,j2  ,i3) = out(i1,j2  ,i3) + in(i1,i2-1,i3)*    w(i2-1) &
                                            + in(i1  ,i2,i3)*(1.-w(i2  ))
          out(i1,j2+1,i3) = out(i1,j2+1,i3) + in(i1  ,i2,i3)
        end do
      else
        do i1=l(1),u(1)
          out(i1,j2  ,i3) = out(i1,j2  ,i3) + in(i1,i2-1,i3)*    w(i2-1) &
                                            + in(i1  ,i2,i3)*(1.-w(i2  ))
          out(i1,j2+1,i3) = out(i1,j2+1,i3) + in(i1  ,i2,i3)*    w(i2  ) &
                                            + in(i1,i2+1,i3)*(1.-w(i2+1))
        end do
      end if
    end do
    j2 = j2 + 2
  end do
  !print '("Y",5x,3(a,":",2i4,5x),f10.6)', &
    !'x', l(1),u(1), 'y', lo,j2-1, 'z',l(3),u(3), sum(out)/product(shape(out))
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE y_interp3 (in, w, l, u, out, lo)
  real, intent(in):: in(:,:,:), w(:,:,:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), lo
  !............................................................................
  integer:: i1, i2, i3, j2, n(3)
  !----------------------------------------------------------------------------
  n = shape(out)
  out = 0.0
  j2 = lo
  !print '(3(4x,a,":",3i4))', 'y_interp3 lt', l, 'ut', u, 'ls', lo
  do i2=l(2),u(2)
    do i3=l(3),u(3)
      if (j2<1 .or. j2>n(2)) then
        cycle
      else if (j2==n(2)) then
        do i1=l(1),u(1)
          out(i1,j2  ,i3) = out(i1,j2  ,i3) + in(i1,i2-1,i3)*    w(i1,i2-1,i3) &
                                            + in(i1  ,i2,i3)*(1.-w(i1  ,i2,i3))
        end do
      else if (i2==1) then
        do i1=l(1),u(1)
          out(i1,j2  ,i3) = out(i1,j2  ,i3) + in(i1  ,i2,i3)
          out(i1,j2+1,i3) = out(i1,j2+1,i3) + in(i1  ,i2,i3)*    w(i1  ,i2,i3) &
                                            + in(i1,i2+1,i3)*(1.-w(i1,i2+1,i3))
        end do
      else if (i2==n(2)) then
        do i1=l(1),u(1)
          out(i1,j2  ,i3) = out(i1,j2  ,i3) + in(i1,i2-1,i3)*    w(i1,i2-1,i3) &
                                            + in(i1  ,i2,i3)*(1.-w(i1  ,i2,i3))
          out(i1,j2+1,i3) = out(i1,j2+1,i3) + in(i1  ,i2,i3)
        end do
      else
        do i1=l(1),u(1)
          out(i1,j2  ,i3) = out(i1,j2  ,i3) + in(i1,i2-1,i3)*    w(i1,i2-1,i3) &
                                            + in(i1  ,i2,i3)*(1.-w(i1  ,i2,i3))
          out(i1,j2+1,i3) = out(i1,j2+1,i3) + in(i1  ,i2,i3)*    w(i1  ,i2,i3) &
                                            + in(i1,i2+1,i3)*(1.-w(i1,i2+1,i3))
        end do
      end if
    end do
    j2 = j2 + 2
  end do
  !print '("Y",5x,3(a,":",2i4,5x),f10.6)', &
    !'x', l(1),u(1), 'y', lo,j2-1, 'z',l(3),u(3), sum(out)/product(shape(out))
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE z_interp1 (in, w, l, u, out, lo)
  real, intent(in):: in(:,:,:), w(:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), lo
  !............................................................................
  integer:: i1, i2, i3, j3, n(3)
  !----------------------------------------------------------------------------
  n = shape(out)
  out = 0.0
  j3 = lo
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      if (j3<1 .or. j3>n(3)) then
        cycle
      else if (j3==n(3)) then
        do i1=l(1),u(1)
          out(i1,i2,j3  ) = out(i1,i2,j3  ) + in(i1,i2,i3-1)*    w(i3-1) &
                                            + in(i1,i2,i3  )*(1.-w(i3  ))
        end do
      else if (i3==1) then
        do i1=l(1),u(1)
          out(i1,i2,j3  ) = out(i1,i2,j3  ) + in(i1,i2,i3  )
          out(i1,i2,j3+1) = out(i1,i2,j3+1) + in(i1,i2,i3  )*    w(i3  ) &
                                            + in(i1,i2,i3+1)*(1.-w(i3+1))
        end do
      else if (i3==n(3)) then
        do i1=l(1),u(1)
          out(i1,i2,j3  ) = out(i1,i2,j3  ) + in(i1,i2,i3-1)*    w(i3-1) &
                                            + in(i1,i2,i3  )*(1.-w(i3  ))
          out(i1,i2,j3+1) = out(i1,i2,j3+1) + in(i1,i2,i3  )
        end do
      else
        do i1=l(1),u(1)
          out(i1,i2,j3  ) = out(i1,i2,j3  ) + in(i1,i2,i3-1)*    w(i3-1) &
                                            + in(i1,i2,i3  )*(1.-w(i3  ))
          out(i1,i2,j3+1) = out(i1,i2,j3+1) + in(i1,i2,i3  )*    w(i3  ) &
                                            + in(i1,i2,i3+1)*(1.-w(i3+1))
        end do
      end if
    end do
    j3 = j3 + 2
  end do
  !print '("Z",5x,3(a,":",2i4,5x),3f10.6)', &
    !'x', l(1),u(1), 'y', l(2),u(2), 'z',lo,j3-1, sum(out)/product(shape(out)), &
    !minval(out), maxval(out)
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE z_interp3 (in, w, l, u, out, lo)
  real, intent(in):: in(:,:,:), w(:,:,:)
  real, intent(out):: out(:,:,:)
  integer, intent(in):: l(3), u(3), lo
  !............................................................................
  integer:: i1, i2, i3, j3, n(3)
  !----------------------------------------------------------------------------
  n = shape(out)
  out = 0.0
  j3 = lo
  !print '(3(4x,a,":",3i4))', 'z_interp3 lt', l, 'ut', u, 'ls', lo
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      if (j3<1 .or. j3>n(3)) then
        cycle
      else if (j3==n(3)) then
        do i1=l(1),u(1)
          out(i1,i2,j3  ) = out(i1,i2,j3  ) + in(i1,i2,i3-1)*    w(i1,i2,i3-1) &
                                            + in(i1,i2,i3  )*(1.-w(i1,i2,i3  ))
        end do
      else if (i3==1) then
        do i1=l(1),u(1)
          out(i1,i2,j3  ) = out(i1,i2,j3  ) + in(i1,i2,i3  )
          out(i1,i2,j3+1) = out(i1,i2,j3+1) + in(i1,i2,i3  )*    w(i1,i2,i3  ) &
                                            + in(i1,i2,i3+1)*(1.-w(i1,i2,i3+1))
        end do
      else if (i3==n(3)) then
        do i1=l(1),u(1)
          out(i1,i2,j3  ) = out(i1,i2,j3  ) + in(i1,i2,i3-1)*    w(i1,i2,i3-1) &
                                            + in(i1,i2,i3  )*(1.-w(i1,i2,i3  ))
          out(i1,i2,j3+1) = out(i1,i2,j3+1) + in(i1,i2,i3  )
        end do
      else
        do i1=l(1),u(1)
          out(i1,i2,j3  ) = out(i1,i2,j3  ) + in(i1,i2,i3-1)*    w(i1,i2,i3-1) &
                                            + in(i1,i2,i3  )*(1.-w(i1,i2,i3  ))
          out(i1,i2,j3+1) = out(i1,i2,j3+1) + in(i1,i2,i3  )*    w(i1,i2,i3  ) &
                                            + in(i1,i2,i3+1)*(1.-w(i1,i2,i3+1))
        end do
      end if
    end do
    j3 = j3 + 2
  end do
  !print '("Z",5x,3(a,":",2i4,5x),3f10.6)', &
    !'x', l(1),u(1), 'y', l(2),u(2), 'z',lo,j3-1, sum(out)/product(shape(out)), &
    !minval(out), maxval(out)
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE interp1 (self, source, w, ls, lt, ut, target)
  class(prolong_t):: self
  real, intent(in):: source(:,:,:), w(:)
  real, intent(out):: target(:,:,:)
  integer, intent(in):: lt(3), ut(3), ls(3)
  !............................................................................
  integer:: m(3), l(3), u(3)
  real, allocatable:: target1(:,:,:), target2(:,:,:)
  !----------------------------------------------------------------------------
  m = shape(source)
  allocate (target1(m(1),m(2),m(3)), target2(m(1),m(2),m(3)))
  l = lt
  u = [ut(1),ut(2),ls(3)+(ut(3)-lt(3)+1)/2-1]
  call self%z_interp (source , w, l, u, target1, lt(1))
  u = [ut(1),ls(2)+(ut(2)-lt(2)+1)/2-1,ut(3)]
  call self%y_interp (target1, w, l, u, target2, lt(2))
  u = [ls(1)+(ut(1)-lt(1)+1)/2-1,ut(2),ut(3)]
  call self%x_interp (target2, w, l, u, target , lt(3))
  deallocate (target1, target2)
  u = ut
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE interp3 (self, source, w, ls, lt, ut, target)
  class(prolong_t):: self
  real, intent(in):: source(:,:,:), w(:,:,:)
  real, intent(out):: target(:,:,:)
  integer, intent(in):: lt(3), ut(3), ls(3)
  !............................................................................
  integer:: m(3), l(3), u(3)
  real, allocatable:: target1(:,:,:), target2(:,:,:)
  !----------------------------------------------------------------------------
  m = shape(source)
  allocate (target1(m(1),m(2),m(3)), target2(m(1),m(2),m(3)))
  l = lt
  u = [ut(1),ut(2),ls(3)+(ut(3)-lt(3)+1)/2-1]
  call self%z_interp (source , w, l, u, target1, lt(1))
  u = [ut(1),ls(2)+(ut(2)-lt(2)+1)/2-1,ut(3)]
  call self%y_interp (target1, w, l, u, target2, lt(2))
  u = [ls(1)+(ut(1)-lt(1)+1)/2-1,ut(2),ut(3)]
  call self%x_interp (target2, w, l, u, target , lt(3))
  deallocate (target1, target2)
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE test (n)
  integer, intent(in):: n
  real, allocatable, dimension(:,:,:):: in, out, w3
  real, allocatable, dimension(:):: w1
  integer:: l(3), u(3), m(3), lo, lc(3), i
  real(8):: start, ns
  !-----------------------------------------------------------------------------
  l = 1
  u = n/2
  lo = 1
  !-----------------------------------------------------------------------------
  call alloc
  start = omp_get_wtime()
  call prolong%x_interp (in, w1, l, u, out, lo)
  call results ('x_interp1')
  !-----------------------------------------------------------------------------
  call alloc
  start = omp_get_wtime()
  call prolong%y_interp (in, w1, l, u, out, lo)
  call results ('y_interp1')
  !-----------------------------------------------------------------------------
  call alloc
  start = omp_get_wtime()
  call prolong%z_interp (in, w1, l, u, out, lo)
  call results ('z_interp1')
  !-----------------------------------------------------------------------------
  call alloc
  start = omp_get_wtime()
  call prolong%x_interp (in, w3, l, u, out, lo)
  call results ('x_interp3')
  !-----------------------------------------------------------------------------
  call alloc
  start = omp_get_wtime()
  call prolong%y_interp (in, w3, l, u, out, lo)
  call results ('y_interp3')
  !-----------------------------------------------------------------------------
  call alloc
  start = omp_get_wtime()
  call prolong%z_interp (in, w3, l, u, out, lo)
  call results ('z_interp3')
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call alloc
  start = omp_get_wtime()
  call prolong%interp (in, w3, [1,1,1], [1  ,1  ,1  ], [2  ,2  ,2  ], out) ! LLL
  call prolong%interp (in, w3, [1,1,1], [1  ,3  ,1  ], [2  ,n-2,2  ], out) ! LML
  call prolong%interp (in, w3, [1,1,1], [1  ,n-1,1  ], [2  ,n  ,2  ], out) ! LUL
  call prolong%interp (in, w3, [1,1,1], [3  ,n-1,1  ], [n-2,n  ,2  ], out) ! MUL
  call prolong%interp (in, w3, [1,1,1], [n-1,n-1,1  ], [n  ,n  ,2  ], out) ! UUL
  call prolong%interp (in, w3, [1,1,1], [n-1,3  ,1  ], [n  ,n-2,2  ], out) ! UML
  call prolong%interp (in, w3, [1,1,1], [n-1,1  ,1  ], [n  ,2  ,2  ], out) ! ULL
  call prolong%interp (in, w3, [1,1,1], [3  ,1  ,1  ], [n-2,2  ,2  ], out) ! MLL
  call prolong%interp (in, w3, [1,1,1], [3  ,3  ,1  ], [n-2,n-2,2  ], out) ! MML
  !
  call prolong%interp (in, w3, [1,1,1], [1  ,1  ,n-1], [2  ,2  ,n  ], out) ! LLU
  call prolong%interp (in, w3, [1,1,1], [1  ,3  ,n-1], [2  ,n-2,n  ], out) ! LMU
  call prolong%interp (in, w3, [1,1,1], [1  ,n-1,n-1], [2  ,n  ,n  ], out) ! LUU
  call prolong%interp (in, w3, [1,1,1], [3  ,n-1,n-1], [n-2,n  ,n  ], out) ! MUU
  call prolong%interp (in, w3, [1,1,1], [n-1,n-1,n-1], [n  ,n  ,n  ], out) ! UUU
  call prolong%interp (in, w3, [1,1,1], [n-1,3  ,n-1], [n  ,n-2,n  ], out) ! UMU
  call prolong%interp (in, w3, [1,1,1], [n-1,1  ,n-1], [n  ,2  ,n  ], out) ! ULU
  call prolong%interp (in, w3, [1,1,1], [3  ,1  ,n-1], [n-2,2  ,n  ], out) ! MLU
  call prolong%interp (in, w3, [1,1,1], [3  ,3  ,n-1], [n-2,n-2,n  ], out) ! MMU
  !
  call prolong%interp (in, w3, [1,1,1], [1  ,1  ,3  ], [2  ,2  ,n-2], out) ! LLM
  call prolong%interp (in, w3, [1,1,1], [3  ,1  ,3  ], [n-2,2  ,n-2], out) ! MLM
  call prolong%interp (in, w3, [1,1,1], [n-1,1  ,3  ], [n  ,2  ,n-2], out) ! ULM
  call prolong%interp (in, w3, [1,1,1], [1  ,n-1,3  ], [2  ,n  ,n-2], out) ! LUM
  call prolong%interp (in, w3, [1,1,1], [3  ,n-1,3  ], [n-2,n  ,n-2], out) ! MUM
  call prolong%interp (in, w3, [1,1,1], [n-1,n-1,3  ], [n  ,n  ,n-2], out) ! UUM
  call prolong%interp (in, w3, [1,1,1], [1  ,3  ,3  ], [2  ,n-2,n-2], out) ! LMM
  call prolong%interp (in, w3, [1,1,1], [n-1,3  ,3  ], [n  ,n-2,n-2], out) ! UMM
  !
  call results ('guard zones')
  !-----------------------------------------------------------------------------
  call alloc
  start = omp_get_wtime()
  do i=1,26
    call prolong%interp (in, w1, [3,3,3], [1,1,1], [n,n,n], out)
  end do
  call results ('all 26 w1  ')
  print *, n, sum(in)/product(m)
  print *, n, sum(out)/product(m)
  !-----------------------------------------------------------------------------
  call alloc
  start = omp_get_wtime()
  do i=1,26
    call prolong%interp (in, w3, [3,3,3], [1,1,1], [n,n,n], out)
  end do
  call results ('all 26 w3  ')
  m = n
  print *, n, sum(in)/product(m)
  print *, n, sum(out)/product(m)
  !-----------------------------------------------------------------------------
contains
  subroutine results (label)
    real(8):: mus, ms, ns
    character(len=*):: label
    mus = (omp_get_wtime()-start)*1e6
    ns = mus*1e3/product(m)
    ms = 1e-3*mus
    print '(a," with",i4," point:",f7.1," ns/pt",f7.1," mus,",f7.2," ms")', &
      label, n, ns, mus, ms
  end subroutine
  subroutine alloc
    if (allocated(in)) then
      !deallocate (in, out, w1, w3)
    end if
    m = n
    if (.not.allocated(in)) then
      allocate (in (m(1),m(2),m(3)))
      allocate (out(m(1),m(2),m(3)))
      allocate (w3 (m(1),m(2),m(3)))
      allocate (w1 (n))
    end if
    in = 1.0
    w1 = 0.1
    w3 = 0.1
  end subroutine
END SUBROUTINE test

END MODULE prolong_mod

PROGRAM test
  USE prolong_mod
  call prolong%test (16)
  call prolong%test (32)
END PROGRAM
