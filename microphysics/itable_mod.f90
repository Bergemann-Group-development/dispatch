!===============================================================================
!> Module handling irregular tables effectively, by treating them as square
!> tables, while only allocating arrays in limited index ranges.
!===============================================================================
MODULE itable_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Data type for linear scales, used for axis values
  !-----------------------------------------------------------------------------
  type, public:: independent_t
   integer:: n
   character(len=32):: name
   real:: min=huge(1.), step, max=-huge(1.)
  contains
    procedure:: read => read_independent
    procedure:: test => test_independent
  end type
  !-----------------------------------------------------------------------------
  ! Data type for table values
  !-----------------------------------------------------------------------------
  type, public:: dependent_t
    integer:: n, nv
    character(32):: name
    real, pointer, contiguous:: v(:,:)
  contains
    procedure:: read => read_values
    procedure:: test => test_values
  end type
  !-----------------------------------------------------------------------------
  ! Data type for table, with index offsets
  !-----------------------------------------------------------------------------
  type, public:: itable_t
    integer:: n1, n2, nv
    type(independent_t):: scale1, scale2
    type(independent_t), pointer:: scales(:)
    type(dependent_t), pointer:: values(:)
  contains
    procedure:: init
    procedure:: lookup1_1d
    procedure:: lookup1_2d
    procedure:: lookup1_3d
    procedure:: lookup_3d
    generic:: lookup => lookup1_1d, lookup1_2d, lookup1_3d, lookup_3d
  end type
  integer:: verbose = 0
CONTAINS

!===============================================================================
!> Read a table with a regularly spaced 1st access and individual segments of
!> regularly spaced values on the 2nd axis
!===============================================================================
SUBROUTINE init (self, file)
  class(itable_t):: self
  character(len=*):: file
  !.............................................................................
  integer:: i, iv, n, u, o
  type(dependent_t), pointer:: values
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%init')
  verbose = 1
  !-----------------------------------------------------------------------------
  ! Generate test values
  !-----------------------------------------------------------------------------
  if (file == 'test.eos') then
    self%nv = 16
    call self%scale1%test ()
    call self%scale2%test ()
    allocate (self%values(self%scale1%n))
    do i=1,self%scale1%n
      allocate (self%scales(i))
      call self%scales(i)%test (i-1)
      o = floor((self%scales(i)%min-self%scale2%min)/self%scale2%step)
      values => self%values(i)
      allocate (self%values(i)%v(self%nv,o:o+self%scales(i)%n-1))
      do iv=1,self%nv
        call self%values(i)%test (i, iv)
      end do
    end do
  !-----------------------------------------------------------------------------
  ! Read EOS file created by Python package dispatch.EOS
  !-----------------------------------------------------------------------------
  else
    u = io_unit%tmp
    open (u,file=file,form='unformatted',status='old')
    read (u) self%nv
    if (verbose>0) print *,'nv',self%nv
    write (stdlog,'(a,i2)') 'itable_t%read: nv =', self%nv
    !---------------------------------------------------------------------------
    ! Read linear scales specifying the full extend of axis 1 and axis 2
    !---------------------------------------------------------------------------
    call self%scale1%read (u)
    call self%scale2%read (u)
    !---------------------------------------------------------------------------
    ! Read linear scales for each value along axis 1, and allocate arrays for
    ! values, with index offsets
    !---------------------------------------------------------------------------
    allocate (self%values(self%scale1%n))
    do i=1,self%scale1%n
      allocate (self%scales(i))
      call self%scales(i)%read (u)
      o = floor((self%scales(i)%min-self%scale2%min)/self%scale2%step)
      allocate (self%values(i)%v(self%nv,o:o+self%scales(i)%n-1))
      do iv=1,self%nv
        call self%values(i)%read (u, iv)
      end do
    end do
    close(u)
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Read a byte-by-byte representation of a string
!===============================================================================
SUBROUTINE read_str (u, str)
  integer:: u, n, m, i
  character(len=1):: c(32)
  character(len=*) str
  !-----------------------------------------------------------------------------
  call trace%begin ('dependent_t%read')
  read (u) m
  read (u) (c(i),i=1,m)
  do i=1,m
    str(i:i)=c(i)
  end do
END SUBROUTINE read_str

!===============================================================================
!> Read scale parameters
!===============================================================================
SUBROUTINE read_independent (self, u)
  class(independent_t):: self
  integer:: u
  !.............................................................................
  integer:: n, m
  real:: amin, amax, step
  !-----------------------------------------------------------------------------
  call trace%begin ('independent_t%read')
  call read_str (u, self%name)
  read (u) self%n
  read (u) self%min, self%step, self%max
  if (verbose > 0) then
    write (stdlog,'(a,i4,1p,3e12.3)') &
      'independent_t%read:', self%n, self%min, self%step, self%max
  end if
  call trace%end ()
END SUBROUTINE read_independent

!===============================================================================
!> Generaate test scale parameters
!===============================================================================
SUBROUTINE test_independent (self, i)
  class(independent_t):: self
  integer, optional:: i
  !.............................................................................
  integer:: n, m
  real:: amin, amax, step
  !-----------------------------------------------------------------------------
  call trace%begin ('independent_t%test')
  self%name = 'indep'
  self%n = 32
  if (present (i)) then
    self%min = i
    self%max = i + self%n
  else
    self%min = 0.0
    self%max = self%n
  end if
  self%step = 1.0
  call trace%end ()
END SUBROUTINE test_independent

!===============================================================================
!> Read table values
!===============================================================================
SUBROUTINE read_values (self, u, iv)
  class(dependent_t):: self
  integer:: u, iv
  !.............................................................................
  integer:: n, m, i
  character(len=1):: c(32)
  character(len=32):: name
  !-----------------------------------------------------------------------------
  call trace%begin ('dependent_t%read')
  read (u) m
  read (u) (c(i),i=1,m)
  do i=1,m
    self%name(i:i)=c(i)
  end do
  read (u) self%n
  read (u) self%v(iv,:)
  if (verbose > 1) then
    write (stdlog,'(a6,2x,3f11.2)') name(1:m), self%v(1,1:3)
  end if
  call trace%end ()
END SUBROUTINE read_values

!===============================================================================
!> Test table values
!===============================================================================
SUBROUTINE test_values (self, i1, iv)
  class(dependent_t):: self
  integer:: i1, iv
  integer:: i2
  !.............................................................................
  integer:: n, m, i
  character(len=1):: c(32)
  character(len=32):: name
  !-----------------------------------------------------------------------------
  call trace%begin ('dependent_t%test')
#ifdef __GFORTRAN__
  write (self%name,'(a,i1)') 'test', iv
#endif
  self%n = 30
  do i2=lbound(self%v,2),ubound(self%v,2)
    self%v(iv,i2) = i1*100 + i2
  end do
  call trace%end ()
END SUBROUTINE test_values

!===============================================================================
!> 1-D lookup
!===============================================================================
SUBROUTINE lookup1_1d (self, iv, x1, x2, v)
  class(itable_t):: self
  integer:: iv
  real:: x1(:), x2(:), v(:)
  !.............................................................................
  real   , pointer, contiguous:: p1(:), p2(:), a(:,:)
  integer, pointer, contiguous:: q1(:), q2(:)
  integer:: i, n
  integer, save:: itimer=0
  real, pointer:: values1(:), values2(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%lookup1_1d', itimer=itimer)
  n = size(x1)
  allocate (p1(n), p2(n), q1(n), q2(n), a(n,4))
  !-----------------------------------------------------------------------------
  ! Compute interpolation indices and weights in a vectorized loop
  !-----------------------------------------------------------------------------
  do i=1,n
    p1(i) = (x1(i)-self%scale1%min)/self%scale1%step
    p2(i) = (x2(i)-self%scale2%min)/self%scale2%step
    q1(i) = p1(i)
    p1(i) = p1(i) - q1(i)
    q2(i) = p2(i)
    p2(i) = p2(i) - q2(i)
  end do
  !-----------------------------------------------------------------------------
  ! Pick up corner values in a loop with indirect addressing (slow)
  !-----------------------------------------------------------------------------
  do iv=1,self%nv
    do i=1,n
      values1 => self%values(q1(i)  )%v(iv,:)
      values2 => self%values(q1(i)+1)%v(iv,:)
      a(i,1) = values1(q2(i)  )
      a(i,2) = values1(q2(i)+1)
      a(i,3) = values2(q2(i)  )
      a(i,4) = values2(q2(i)+1)
    end do
    !---------------------------------------------------------------------------
    ! Bilinear interpolation
    !---------------------------------------------------------------------------
    do i=1,n
      v(i) = (1.-p1(i))*(a(i,1)*(1.-p2(i))+a(i,2)*p2(i)) &
               + p1(i) *(a(i,3)*(1.-p2(i))+a(i,4)*p2(i))
    end do
  end do
  deallocate (p1, p2, q1, q2, a)
  call trace%end (itimer)
END SUBROUTINE lookup1_1d

!===============================================================================
!> 2-D lookup
!===============================================================================
SUBROUTINE lookup1_2d (self, iv, x1, x2, v)
  class(itable_t):: self
  integer:: iv
  real:: x1(:,:), x2(:,:), v(:,:)
  !.............................................................................
  real   , pointer, contiguous:: p1(:), p2(:), a(:,:)
  integer, pointer, contiguous:: q1(:), q2(:)
  integer:: i, n, i2, n2
  integer, save:: itimer=0
  real, pointer:: values1(:), values2(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%lookup1_2d', itimer=itimer)
  n  = size(x1,1)
  n2 = size(x1,2)
  allocate (p1(n), p2(n), q1(n), q2(n), a(n,4))
  !-----------------------------------------------------------------------------
  ! Compute interpolation indices and weights in a vectorized loop
  !-----------------------------------------------------------------------------
  do i2=1,n2
  do i=1,n
    p1(i) = (x1(i,i2)-self%scale1%min)/self%scale1%step
    p2(i) = (x2(i,i2)-self%scale2%min)/self%scale2%step
    q1(i) = p1(i)
    p1(i) = p1(i) - q1(i)
    q2(i) = p2(i)
    p2(i) = p2(i) - q2(i)
  end do
  !-----------------------------------------------------------------------------
  ! Pick up corner values in a loop with indirect addressing (slow)
  !-----------------------------------------------------------------------------
  do i=1,n
    values1 => self%values(q1(i)  )%v(iv,:)
    values2 => self%values(q1(i)+1)%v(iv,:)
    a(i,1) = values1(q2(i)  )
    a(i,2) = values1(q2(i)+1)
    a(i,3) = values2(q2(i)  )
    a(i,4) = values2(q2(i)+1)
  end do
  !-----------------------------------------------------------------------------
  ! Bilinear interpolation
  !-----------------------------------------------------------------------------
  do i=1,n
    v(i,i2) = (1.-p1(i))*(a(i,1)*(1.-p2(i))+a(i,2)*p2(i)) &
                + p1(i) *(a(i,3)*(1.-p2(i))+a(i,4)*p2(i))
  end do
  end do
  deallocate (p1, p2, q1, q2, a)
  call trace%end (itimer)
END SUBROUTINE lookup1_2d

!===============================================================================
!> 3-D lookup
!===============================================================================
SUBROUTINE lookup1_3d (self, iv, x1, x2, v)
  class(itable_t):: self
  real:: x1(:,:,:), x2(:,:,:), v(:,:,:)
  !.............................................................................
  real   , pointer, contiguous:: p1(:), p2(:), a(:,:)
  integer, pointer, contiguous:: q1(:), q2(:)
  integer:: i, n, iv, i2, n2, i3, n3, nv
  integer, save:: itimer=0
  real, pointer:: values1(:,:), values2(:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%lookup1_3d', itimer=itimer)
  n  = size(v,1)
  n2 = size(v,2)
  n3 = size(v,3)
  allocate (p1(n), p2(n), q1(n), q2(n), a(n,4))
  !-----------------------------------------------------------------------------
  ! Compute interpolation indices and weights in a vectorized loop
  !-----------------------------------------------------------------------------
  do i3=1,n3
  do i2=1,n2
#ifdef __GFORTRAN__
    !$omp simd
#endif
    do i=1,n
      p1(i) = (x1(i,i2,i3)-self%scale1%min)/self%scale1%step
      p2(i) = (x2(i,i2,i3)-self%scale2%min)/self%scale2%step
      q1(i) = p1(i)
      p1(i) = p1(i) - q1(i)
      q2(i) = p2(i)
      p2(i) = p2(i) - q2(i)
    end do
    !---------------------------------------------------------------------------
    ! Pick up corner values in a loop with indirect addressing (slow)
    !---------------------------------------------------------------------------
#ifdef __GFORTRAN__
    !$omp simd
#endif
    do i=1,n
      values1 => self%values(q1(i)  )%v
      values2 => self%values(q1(i)+1)%v
      a(i,1) = values1(iv,q2(i)  )
      a(i,2) = values1(iv,q2(i)+1)
      a(i,3) = values2(iv,q2(i)  )
      a(i,4) = values2(iv,q2(i)+1)
    end do
    !---------------------------------------------------------------------------
    ! Bilinear interpolation
    !---------------------------------------------------------------------------
#ifdef __GFORTRAN__
    !$omp simd
#endif
    do i=1,n
      v(i,i2,i3) = (1.-p1(i))*(a(i,1)*(1.-p2(i))+a(i,2)*p2(i)) &
                     + p1(i) *(a(i,3)*(1.-p2(i))+a(i,4)*p2(i))
    end do
  end do
  end do
  deallocate (p1, p2, q1, q2, a)
  call trace%end (itimer)
END SUBROUTINE lookup1_3d

!===============================================================================
!> 3-D lookup; choose btw v(nv,n,n,n) and v(n,n,n,nv)
!===============================================================================
SUBROUTINE lookup_3d (self, x1, x2, v)
  class(itable_t):: self
  real:: x1(:,:,:), x2(:,:,:), v(:,:,:,:)
  integer:: m(4)
  !-----------------------------------------------------------------------------
  m = shape(v)
  if (m(1) == m(2)) then
    call lookup_3d_nnv (self, x1, x2, v)
  else
    call lookup_3d_nvn (self, x1, x2, v)
  end if
END SUBROUTINE lookup_3d

!===============================================================================
!> 3-D lookup, order (n,n,n,nv)
!===============================================================================
SUBROUTINE lookup_3d_nnv (self, x1, x2, v)
  class(itable_t):: self
  real:: x1(:,:,:), x2(:,:,:), v(:,:,:,:)
  !.............................................................................
  real   , pointer, contiguous:: p1(:), p2(:), a(:,:,:)
  integer, pointer, contiguous:: q1(:), q2(:)
  integer:: i, n, iv, i2, n2, i3, n3, nv
  integer, save:: itimer=0
  real, pointer:: values1(:,:), values2(:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%lookup_3d_nnv', itimer=itimer)
  !print *,'lookup_3d_nnv'
  n  = size(v,1)
  n2 = size(v,2)
  n3 = size(v,3)
  nv = size(v,4)
  allocate (p1(n), p2(n), q1(n), q2(n), a(n,nv,4))
  !-----------------------------------------------------------------------------
  ! Compute interpolation indices and weights in a vectorized loop
  !-----------------------------------------------------------------------------
  do i3=1,n3
  do i2=1,n2
#ifdef __GFORTRAN__
    !$omp simd
#endif
    do i=1,n
      p1(i) = (x1(i,i2,i3)-self%scale1%min)/self%scale1%step
      p2(i) = (x2(i,i2,i3)-self%scale2%min)/self%scale2%step
      q1(i) = p1(i)
      p1(i) = p1(i) - q1(i)
      q2(i) = p2(i)
      p2(i) = p2(i) - q2(i)
    end do
    !---------------------------------------------------------------------------
    ! Pick up corner values in a loop with indirect addressing (slow)
    !---------------------------------------------------------------------------
    do i=1,n
      values1 => self%values(q1(i)  )%v
      values2 => self%values(q1(i)+1)%v
#ifdef __GFORTRAN__
      !$omp simd
#endif
      do iv=1,nv
        a(i,iv,1) = values1(iv,q2(i)  )
        a(i,iv,2) = values1(iv,q2(i)+1)
        a(i,iv,3) = values2(iv,q2(i)  )
        a(i,iv,4) = values2(iv,q2(i)+1)
      end do
    end do
    !---------------------------------------------------------------------------
    ! Bilinear interpolation
    !---------------------------------------------------------------------------
    do iv=1,nv
#ifdef __GFORTRAN__
      !$omp simd
#endif
      do i=1,n
        v(i,i2,i3,iv) = (1.-p1(i))*(a(i,iv,1)*(1.-p2(i))+a(i,iv,2)*p2(i)) &
                          + p1(i) *(a(i,iv,3)*(1.-p2(i))+a(i,iv,4)*p2(i))
      end do
    end do
  end do
  end do
  deallocate (p1, p2, q1, q2, a)
  call trace%end (itimer)
END SUBROUTINE lookup_3d_nnv

!===============================================================================
!> 3-D lookup, order (nv,n,n,n)
!===============================================================================
SUBROUTINE lookup_3d_nvn (self, x1, x2, v)
  class(itable_t):: self
  real:: x1(:,:,:), x2(:,:,:), v(:,:,:,:)
  !.............................................................................
  real   , pointer, contiguous:: p1(:), p2(:), a(:,:,:)
  integer, pointer, contiguous:: q1(:), q2(:)
  integer:: i, n, iv, i2, n2, i3, n3, nv
  integer, save:: itimer=0
  real, pointer:: values1(:,:), values2(:,:)
  real:: c1, c2
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%lookup_3d', itimer=itimer)
  nv = size(v,1)
  n  = size(v,2)
  n2 = size(v,3)
  n3 = size(v,4)
  allocate (p1(n), p2(n), q1(n), q2(n), a(nv,n,4))
  !-----------------------------------------------------------------------------
  ! Compute interpolation indices and weights in a vectorized loop
  !-----------------------------------------------------------------------------
  c1 = 1.0/self%scale1%step
  c2 = 1.0/self%scale2%step
  do i3=1,n3
  do i2=1,n2
#ifdef __GFORTRAN__
    !$omp simd
#endif
    do i=1,n
      p1(i) = (x1(i,i2,i3)-self%scale1%min)*c1
      p2(i) = (x2(i,i2,i3)-self%scale2%min)*c2
      q1(i) = p1(i)
      p1(i) = p1(i) - q1(i)
      q2(i) = p2(i)
      p2(i) = p2(i) - q2(i)
    end do
    !---------------------------------------------------------------------------
    ! Pick up corner values in a loop with indirect addressing (slow)
    !---------------------------------------------------------------------------
    do i=1,n
      values1 => self%values(q1(i)  )%v
      values2 => self%values(q1(i)+1)%v
#ifdef __GFORTRAN__
      !$omp simd
#endif
      do iv=1,nv
        a(iv,i,1) = values1(iv,q2(i)  )
        a(iv,i,2) = values1(iv,q2(i)+1)
        a(iv,i,3) = values2(iv,q2(i)  )
        a(iv,i,4) = values2(iv,q2(i)+1)
      end do
    end do
    !---------------------------------------------------------------------------
    ! Bilinear interpolation
    !---------------------------------------------------------------------------
    do i=1,n
#ifdef __GFORTRAN__
      !$omp simd
#endif
      do iv=1,nv
        v(iv,i,i2,i3) = (1.-p1(i))*(a(iv,i,1)*(1.-p2(i))+a(iv,i,2)*p2(i)) &
                          + p1(i) *(a(iv,i,3)*(1.-p2(i))+a(iv,i,4)*p2(i))
      end do
    end do
  end do
  end do
  deallocate (p1, p2, q1, q2, a)
  call trace%end (itimer)
END SUBROUTINE lookup_3d_nvn

END MODULE itable_mod
