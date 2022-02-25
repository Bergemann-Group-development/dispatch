!===============================================================================
PROGRAM test
  USE omp_mod
  USE mpi_mod
  USE io_mod
  USE timer_mod
  USE patch_mod
  USE point_sources_mod
  !.............................................................................
  implicit none
  type(patch_t), pointer :: patch
  type(point_source_t)   :: point_source
  integer                :: n, m(3), i1, i2, i3, ax(1)
  logical, allocatable   :: mask(:,:,:), axis1(:,:,:), axis2(:,:,:)
  real(8), allocatable   :: distance(:,:,:,:)
  integer, allocatable   :: axis(:,:,:), ii(:,:,:)
  !-----------------------------------------------------------------------------
  call omp%init
  call mpi%init
  call io%init
  !.............................................................................
  allocate (patch)
  call patch%pre_init
  patch%size = 1d0
  call patch%init
  print *, 'position:', patch%position
  print *, '    size:', patch%size
  !-----------------------------------------------------------------------------
  ! Measure the cost of setting memory to something
  !-----------------------------------------------------------------------------
  m = patch%gn
  allocate (mask(m(1),m(2),m(3)))
  allocate (ii(m(1),m(2),m(3)))
  n = product(shape(mask))
  call timer%tic
  patch%mem(:,:,:,1,1,1) = exp(1.0+2.0+3.0)
  call timer%toc ('  set 1', n)
  !-----------------------------------------------------------------------------
  ! Measure the cost of setting memory to something a 2nd time
  !-----------------------------------------------------------------------------
  call timer%tic
  patch%mem(:,:,:,1,1,1) = 2.0
  call timer%toc ('  set 2', n)
  !-----------------------------------------------------------------------------
  ! Measure the cost of testing & setting
  !-----------------------------------------------------------------------------
  call timer%tic
  do i3=1,m(3)
    do i2=1,m(2)
      do i1=1,m(1)
        ii(i1,i2,i3) = i1+i2+i3
        if (i1 <= i2) then
          patch%mem(i1,i2,i3,1,1,1) = exp((i1+i2+i3)*0.01)
        else
          patch%mem(i1,i2,i3,1,1,1) = exp((i1+i2+i3)*0.02)
        end if
      end do
    end do
  end do
  call timer%toc ('if-then', n)
  !-----------------------------------------------------------------------------
  ! Measure the cost of making a mask
  !-----------------------------------------------------------------------------
  call timer%tic
  do i3=1,m(3)
    do i2=1,m(2)
      do i1=1,m(1)
        if (i1 <= i2) then
          mask(i1,i2,i3) = .true.
        else
          mask(i1,i2,i3) = .false.
        end if
      end do
    end do
  end do
  call timer%toc ('   mask', n)
  !-----------------------------------------------------------------------------
  ! Measure the cost of setting memory to 2, with a mask
  !-----------------------------------------------------------------------------
  call timer%tic
  where (mask)
    patch%mem(:,:,:,1,1,1) = exp(ii*0.01)
  else where
    patch%mem(:,:,:,1,1,1) = exp(ii*0.02)
  end where
  call timer%toc ('  where', n)
  !-----------------------------------------------------------------------------
  ! Measure the cost of a merge
  !-----------------------------------------------------------------------------
  call timer%tic
  do i3=1,m(3)
    do i2=1,m(2)
      do i1=1,m(1)
        patch%mem(i1,i2,i3,1,1,1) = merge (&
                                    exp((i1+i2+i3)*0.01), &
                                    exp((i1+i2+i3)*0.02), &
                                    mask(i1,i2,i3))
      end do
    end do
  end do
  call timer%toc ('  merge', n)
  !-----------------------------------------------------------------------------
  ! Set up a point source outside the patch, along space diagonal
  !-----------------------------------------------------------------------------
  allocate (distance(3,m(1),m(2),m(3)))
  point_source%position = [-2,-2,-2]
  call timer%tic
  do i3=1,m(3)
    do i2=1,m(2)
      do i1=1,m(1)
        distance(:,i1,i2,i3) = ([i1,i2,i3]-patch%offset)*patch%ds + &
          patch%position - point_source%position
      end do
    end do
  end do
  call timer%toc ('   dist', n)
  !-----------------------------------------------------------------------------
  ! Measure the nearest axis
  !-----------------------------------------------------------------------------
  allocate (axis(m(1),m(2),m(3)))
  call timer%tic
  do i3=1,m(3)
    do i2=1,m(2)
      do i1=1,m(1)
        ax = maxloc(abs(distance(:,i1,i2,i3)))
        axis(i1,i2,i3) = ax(1)
      end do
    end do
  end do
  call timer%toc ('   axis', n)
  print *, 'axis counts:', count(axis==1), count(axis==2), count(axis==3)
  !-----------------------------------------------------------------------------
  ! Measure using the nearest axis
  !-----------------------------------------------------------------------------
  call timer%tic
  where (axis==1)
    patch%mem(:,:,:,1,1,1) = exp(ii*0.01)
  else where (axis==2)
    patch%mem(:,:,:,1,1,1) = exp(ii*0.02)
  else where
    patch%mem(:,:,:,1,1,1) = exp(ii*0.03)
  end where
  call timer%toc ('where 1', n)
  !-----------------------------------------------------------------------------
  ! Measure using the nearest axis
  !-----------------------------------------------------------------------------
  allocate (axis1(m(1),m(2),m(3)))
  allocate (axis2(m(1),m(2),m(3)))
  axis1 = axis==1
  axis2 = axis==2
  call timer%tic
  where (axis1)
    patch%mem(:,:,:,1,1,1) = exp(ii*0.01)
  else where (axis2)
    patch%mem(:,:,:,1,1,1) = exp(ii*0.02)
  else where
    patch%mem(:,:,:,1,1,1) = exp(ii*0.03)
  end where
  call timer%toc ('where 2', n)
  !-----------------------------------------------------------------------------
  ! Measure the nearest axis
  !-----------------------------------------------------------------------------
  call timer%tic
  do i3=1,m(3)
    do i2=1,m(2)
      do i1=1,m(1)
        ax = maxloc(abs(distance(:,i1,i2,i3)))
        select case (ax(1))
        case(1)
          patch%mem(i1,i2,i3,1,1,1) = exp((i1+i2+i3)*0.01)
        case(2)
          patch%mem(i1,i2,i3,1,1,1) = exp((i1+i2+i3)*0.02)
        case(3)
          patch%mem(i1,i2,i3,1,1,1) = exp((i1+i2+i3)*0.03)
        end select
      end do
    end do
  end do
  call timer%toc (' select', n)
  !-----------------------------------------------------------------------------
  ! Measure the nearest axis
  !-----------------------------------------------------------------------------
  call timer%tic
  do i3=1,m(3)
    do i2=1,m(2)
      do i1=1,m(1)
        ax = maxloc(abs(distance(:,i1,i2,i3)))
        if (ax(1)==1) then
          patch%mem(i1,i2,i3,1,1,1) = exp((i1+i2+i3)*0.01)
        else if (ax(1)==2) then
          patch%mem(i1,i2,i3,1,1,1) = exp((i1+i2+i3)*0.02)
        else
          patch%mem(i1,i2,i3,1,1,1) = exp((i1+i2+i3)*0.02)
        end if
      end do
    end do
  end do
  call timer%toc ('if/else', n)
END PROGRAM
