!===============================================================================
!===============================================================================
PROGRAM test
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer:: x, y, z, n, dstep(3), ray, pt, p(3), r1,r2,r3, iy
  integer:: ray1, ray2, i1,i2
  integer,dimension(:), allocatable :: npt
  real(8), dimension(:,:,:), pointer:: mem, cb2
  integer, dimension(:,:,:), pointer :: cell
  integer, dimension(:,:,:,:), pointer :: cell3d
  real(8), dimension(:,:), pointer:: cb, new
  integer, dimension(:,:), pointer:: npt2
  logical :: cont
  !-----------------------------------------------------------------------------
 
  call omp%init
  n = 51
  allocate(mem(n,n,n))
  allocate(new(n,n*n))
  allocate(cb(n,n*n))
  allocate(cell(3,n,n*n))
  allocate(cell3d(3,n,n,n))
  allocate(cb2(n,n,n))
  allocate(npt(n*n))
  allocate(npt2(n,n))
  
  dstep = 1
  dstep(3) = 0
  do x=1,n
    do y=1,n
      do z = 1,n
      mem(x,y,z) = 100*x+y+real(z)*0.01
      end do 
    end do
  end do
  
  
  ray = 1
  do y= 1, n
    do z = 1,n
      cont = .true.
      p(1) = 1
      p(2) = y
      p(3) = z
      pt = 1
      do while (cont)
        cell(:,pt,ray) = p
        if (any(p>n)) then
          cont = .false.
          npt(ray) = pt - 1
        else
          p = p + dstep
          pt = pt+1
        end if      
      end do
      ray = ray+1
    end do
  end do
  
  do y=1,n
    do z=1,n
      cont = .true.
      p(1) = 1
      p(2) = y
      p(3) = z
      pt = 1
      ray1 = y; ray2 = z
      do while (cont)
        cell3d(:,ray1,ray2,pt) = p
        if (any(p>n)) then
          cont = .false.
          npt2(y,z) = pt-1
        else
          p = p + dstep
          pt = pt+1
        end if      
      end do  
    end do
  end do
  
  i1=size(cb,1)
  i2=size(cb,2)
  i1=(i1*i2)/2
  call tic
  do x=1,n
    do y=1,npt(x)
      p = cell(:,x,y)
      cb(x,y) = mem(p(1),p(2),p(3))   
    end do
  end do
  call timer%toc('cell based lookup in 2D', i1)
       do x=1,n
    do y = 1, npt(x)
      print *, cb(x,y), 'cell', cell(:,x,y)
    end do
  end do

  
  call tic
  do x=1,n
    y=npt(x)
    p = cell(:,x,y)
    cb(x,1:y) = mem(p(1),p(2),:)   
  end do
  call timer%toc('cell based lookup in 2D, striding over z', i1)
  
   
  
  call tic
  do x=1,n, 5
    do y=1,npt(x)
      p = cell(:,x,y)
      cb(x,y) = mem(p(1),p(2),p(3))   
    end do
    iy=x+1
    do y=1,npt(iy)
      p = cell(:,iy,y)
      cb(iy,y) = mem(p(1),p(2),p(3))   
    end do
    iy=x+2
    do y=1,npt(iy)
      p = cell(:,iy,y)
      cb(iy,y) = mem(p(1),p(2),p(3))   
    end do
    iy=x+3
    do y=1,npt(iy)
      p = cell(:,iy,y)
      cb(iy,y) = mem(p(1),p(2),p(3))   
    end do
    iy=x+4
    do y=1,npt(iy)
      p = cell(:,iy,y)
      cb(iy,y) = mem(p(1),p(2),p(3))   
    end do
  end do
  call timer%toc('cell based lookup in 2D, 1st loop unrolled', i1)
  
  
    call tic
  do x=1,n, 5
    y=npt(x)
    p = cell(:,y,x)
    cb(x,1:y) = mem(p(1),p(2),:)   
    iy=x+1
    y=npt(iy)
    p = cell(:,y,iy)
    cb(iy,1:y) = mem(p(1),p(2),:) 
    iy=x+2
    y=npt(iy)
    p = cell(:,y,iy)
    cb(iy,1:y) = mem(p(1),p(2),:) 
    iy=x+3
    y=npt(iy)
    p = cell(:,y,iy)
    cb(iy,1:y) = mem(p(1),p(2),:) 
    iy=x+4
    y=npt(iy)
    p = cell(:,y,iy)
    cb(iy,1:y) = mem(p(1),p(2),:) 
  end do
  
  call timer%toc('cell based lookup in 2D, striding over z, 1st loop unrolled', i1)
  
  call tic
  do y=1,n
    do z=1,n
      do x=1,npt2(y,z)
      p = cell3d(:,y,z,x)
      cb2(x,y,z) = mem(p(1),p(2),p(3))
      end do   
    end do
  end do
  call timer%toc('cell based lookup in 3D',i1)
  
  ! first test assumes dstep(1,1,0) - two positives and 0
   iy=0; r3 = dstep(2)
  call tic
  do x=1,n   ! this test is along X coordinate, which represents points in ray domain
    iy =iy + r3
    do y = iy,n
      r1=y*n
      do z = 1,n
        new(x,r1+z)=mem(x,y,z)
      end do
    end do

  end do
  call timer%toc('plane based lookup',i1)
  
     iy=0; r3 = n
  call tic
  do x=1,n   ! this test is along X coordinate, which represents points in ray domain
    iy =iy + dstep(2)
    do y = iy,n
      r1=y*n
      new(x,r1:r1+z)=mem(x,y,:)

    end do
  end do
  call timer%toc('plane based lookup2: removing inner loop',i1)
  
END
