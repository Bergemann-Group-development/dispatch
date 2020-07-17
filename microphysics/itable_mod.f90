!===============================================================================
!> Module handling irregular tables effectively, by treating them as square
!> tables, while only allocating arrays in limited index ranges.
!===============================================================================
MODULE itable_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: scale_t
   integer:: n
   real:: min=huge(1.), step, max=-huge(1.)
  contains
    procedure:: read => read_scale
  end type
  type, public:: values_t
    integer:: n, nv
    character(32):: name
    real, pointer, contiguous:: v(:,:)
  contains
    procedure:: read => read_values
  end type
  type, public:: itable_t
    integer:: n1, n2, nv
    type(scale_t):: scale1, scale2
    type(scale_t), pointer:: scales(:)
    type(values_t), pointer:: values(:)
  contains
    procedure:: init
    procedure:: lookup1
    procedure:: lookup2
    procedure:: lookup3
    generic:: lookup => lookup1, lookup2, lookup3
  end type
  integer:: verbose = 0
CONTAINS

!===============================================================================
!> Read scale parameters
!===============================================================================
SUBROUTINE read_scale (self, u)
  class(scale_t):: self
  integer:: u
  !.............................................................................
  integer:: n
  real:: amin, amax, step
  !-----------------------------------------------------------------------------
  call trace%begin ('scale_t%read')
  read (u) self%n
  read (u) self%min, self%step, self%max
  if (verbose > 0) then
    write (stdlog,'(a,i4,1p,3e12.3)') &
      'scale_t%read:', self%n, self%min, self%step, self%max
  end if
  call trace%end ()
END SUBROUTINE read_scale

!===============================================================================
!> Read table values
!===============================================================================
SUBROUTINE read_values (self, u, iv, o)
  class(values_t):: self
  integer:: u, iv, o
  !.............................................................................
  integer:: n, m, i
  character(len=4):: c(32)
  character(len=32):: name
  !-----------------------------------------------------------------------------
  call trace%begin ('values_t%read')
  read (u) m
  read (u) (c(i),i=1,m)
  do i=1,m
    self%name(i:i)=c(i)
  end do
  read (u) self%n
  read (u) self%v(iv,:)
  if (verbose > 1) then
    write (stdlog,'(a6,2x,i6,3f11.2)') name(1:m), o, self%v(1,o:o+2)
  end if
  call trace%end ()
END SUBROUTINE read_values

!===============================================================================
!> Read a complete table
!===============================================================================
SUBROUTINE init (self, file)
  class(itable_t):: self
  character(len=*):: file
  !.............................................................................
  integer:: i, iv, n, u, o
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%init')
  verbose = 1
  u = io_unit%tmp
  open (u,file=file,form='unformatted',status='old')
  read (u) self%nv
  write (stdlog,'(a,i2)') 'itable_t%read: nv =', self%nv
  call self%scale1%read (u)
  call self%scale2%read (u)
  allocate (self%values(self%scale1%n))
  if (verbose == 1) verbose = 0
  do i=1,self%scale1%n
    allocate (self%scales(i))
    call self%scales(i)%read (u)
    o = (self%scales(i)%min-self%scale2%min)/self%scale2%step + 0.5
    allocate (self%values(i)%v(self%nv,o:o+self%scales(i)%n-1))
    do iv=1,self%nv
      call self%values(i)%read (u, iv, o)
    end do
  end do
  close(u)
  call trace%end ()
END SUBROUTINE init

SUBROUTINE lookup1 (self, iv, x1, x2, v)
  class(itable_t):: self
  integer:: iv
  real:: x1(:), x2(:), v(:)
  !.............................................................................
  real   , pointer, contiguous:: p1(:), p2(:), a(:,:)
  integer, pointer, contiguous:: q1(:), q2(:)
  integer:: i, n
  integer, save:: itimer=0
  real, pointer, contiguous:: values1(:), values2(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%lookup1', itimer=itimer)
  n = size(x1)
  allocate (p1(n), p2(n), q1(n), q2(n), a(n,4))
  !-----------------------------------------------------------------------------
  ! Compute interpolation indices and weights in a vectorized loop
  !-----------------------------------------------------------------------------
  do i=1,n
    p1(i) = (x1(i)-self%scale1%min)/self%scale2%step
    p2(i) = (x2(i)-self%scale2%min)/self%scale2%step
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
      a(i,2) = values2(q2(i)  )
      a(i,4) = values2(q2(i)+1)
    end do
    !---------------------------------------------------------------------------
    ! Bilinear interpolation
    !---------------------------------------------------------------------------
    do i=1,n
      v(i) = (1.-p1(i))*(a(i,1)*(1.-p2(i))+a(i,2)*p2(i)) &
                  + p1(i) *(a(i,3)*(1.-p2(i))+a(i,4)*p2(i))
    end do
  deallocate (p1, p2, q1, q2, a)
  call trace%end (itimer)
END SUBROUTINE lookup1

SUBROUTINE lookup2 (self, iv, x1, x2, v)
  class(itable_t):: self
  integer:: iv
  real:: x1(:,:), x2(:,:), v(:,:)
  !.............................................................................
  real   , pointer, contiguous:: p1(:), p2(:), a(:,:)
  integer, pointer, contiguous:: q1(:), q2(:)
  integer:: i, n, i2, n2
  integer, save:: itimer=0
  real, pointer, contiguous:: values1(:), values2(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%lookup2', itimer=itimer)
  n  = size(x1,1)
  n2 = size(x1,2)
  allocate (p1(n), p2(n), q1(n), q2(n), a(n,4))
  !-----------------------------------------------------------------------------
  ! Compute interpolation indices and weights in a vectorized loop
  !-----------------------------------------------------------------------------
  do i2=1,n2
  do i=1,n
    p1(i) = (x1(i,i2)-self%scale1%min)/self%scale2%step
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
    a(i,2) = values2(q2(i)  )
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
END SUBROUTINE lookup2

SUBROUTINE lookup3 (self, iv, x1, x2, v)
  class(itable_t):: self
  integer:: iv
  real:: x1(:,:,:), x2(:,:,:), v(:,:,:)
  !.............................................................................
  real   , pointer, contiguous:: p1(:), p2(:), a(:,:)
  integer, pointer, contiguous:: q1(:), q2(:)
  integer:: i, n, i2, n2, i3, n3, nv
  integer, save:: itimer=0
  real, pointer, contiguous:: values1(:,:), values2(:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('itable_t%lookup3', itimer=itimer)
  n  = size(x1,1)
  n2 = size(x1,2)
  n3 = size(x1,3)
  nv = self%nv
  allocate (p1(n), p2(n), q1(n), q2(n), a(n,4))
  !-----------------------------------------------------------------------------
  ! Compute interpolation indices and weights in a vectorized loop
  !-----------------------------------------------------------------------------
  do i3=1,n3
  do i2=1,n2
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
      a(i,1) = values1(iv,q2(i)  )
      a(i,2) = values1(iv,q2(i)+1)
      a(i,2) = values2(iv,q2(i)  )
      a(i,4) = values2(iv,q2(i)+1)
    end do
    !---------------------------------------------------------------------------
    ! Bilinear interpolation
    !---------------------------------------------------------------------------
    do i=1,n
       v(i,i2,i3) = (1.-p1(i))*(a(i,1)*(1.-p2(i))+a(i,2)*p2(i)) &
                      + p1(i) *(a(i,3)*(1.-p2(i))+a(i,4)*p2(i))
    end do
  end do
  end do
  deallocate (p1, p2, q1, q2, a)
  call trace%end (itimer)
END SUBROUTINE lookup3

END MODULE itable_mod
