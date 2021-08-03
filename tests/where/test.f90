PROGRAM test
  USE timer_mod
  implicit none
  real, allocatable:: a(:,:,:)
  integer, allocatable:: i(:,:,:)
  logical, allocatable:: l(:,:,:)
  integer:: m, i1, i2, i3
  !-----------------------------------------------------------------------------
  m = 33
  allocate (a(m,m,m))
  allocate (i(m,m,m))
  allocate (l(m,m,m))
  call timer%tic
  do i3=1,m
  do i2=1,m
  do i1=1,m
    a(i1,i2,i3) = 1.0
    i(i1,i2,i3) = mod(i1+i2+i3,2)
    l(i1,i2,i3) = i(i1,i2,i3)==1
  end do
  end do
  end do
  call timer%toc ('alloc')
  !-----------------------------------------------------------------------------
  call timer%tic
  do i3=1,m
  do i2=1,m
  do i1=1,m
    if (l(i1,i2,i3)) then
      a(i1,i2,i3) = 0.0
    end if
  end do
  end do
  end do
  call timer%toc ('if')
  !-----------------------------------------------------------------------------
  call timer%tic
  do i3=1,m
  do i2=1,m
  do i1=1,m
    a(i1,i2,i3) = merge(a(i1,i2,i3),0.0,l(i1,i2,i3))
  end do
  end do
  end do
  call timer%toc ('merge')
  !-----------------------------------------------------------------------------
  call timer%tic
  where (i==1) a = 0.0
  call timer%toc ('where(i)')
  !-----------------------------------------------------------------------------
  call timer%tic
  where (l) a = 0.0
  call timer%toc ('where(l)', m**3)
END PROGRAM
