!> =============================================================================
!> $Id$
!> =============================================================================
MODULE direct_output_mod
  USE io_mod
  USE trace_mod
  USE omp_timer_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! The output object contains procedures for opening, buffering, writing, and
  ! closing buffered output to a direct access file.
  !-----------------------------------------------------------------------------
  type, public:: direct_output_t
    integer:: dims(4)
    integer:: count
    integer:: nbuf
    integer:: rec=0
    real, dimension(:,:,:,:), pointer :: buffer=>null()
    character(len=128):: filename
    logical:: locked=.false.
    real(8):: position(3)=0d0
    real:: gb=0.0, gbs=0.0
  contains
    procedure:: init
    procedure:: patch
    procedure:: close
    procedure:: test
  end type
  integer(8):: revision=1
  type(direct_output_t), public:: direct_output
CONTAINS

!> =============================================================================
!> Register a filename, total data size, and number of output buffers
!> =============================================================================
SUBROUTINE init (self, filename, dims, nbuf)
  class(direct_output_t)     :: self
  character(len=*), optional :: filename
  integer, optional          :: dims(4), nbuf
  !-----------------------------------------------------------------------------
  call trace_begin('direct_output_mod::init')
  if (present(filename)) then
    print *, 'direct_output_mod: filename = ', trim(filename)
    self%filename = filename
    self%count = 0
  end if
  if (present(dims)) then
    print *, 'direct_output_mod: dims     = ', dims
    self%dims = dims
  end if
  if (present(nbuf)) then
    print *, 'direct_output_mod: nbuf     = ', nbuf
    self%nbuf = nbuf
    self%count = 0
  end if
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE patch (self, data, offset)
  class (direct_output_t):: self
  real, dimension(:,:,:,:) :: data
  integer, dimension(4) :: offset
  !-----------------------------------------------------------------------------
  call trace_begin('direct_output_mod::patch')
  if (self%locked) then
    !$omp critical (direct_cr)
    call output (self, data, offset)
    !$omp end critical (direct_cr)
  else
    call output (self, data, offset)
  end if
  call trace_end
END SUBROUTINE patch

!> =============================================================================
!> =============================================================================
SUBROUTINE output (self, data, offset)
  class (direct_output_t):: self
  integer, dimension(4) :: offset, n, l, u 
  real, dimension(:,:,:,:) :: data
  !-----------------------------------------------------------------------------
  call trace_begin('direct_output_mod::output')
  n = shape(data)
  l = offset
  u = l+n-1
  if (io%debug(3)) &
    print '("direct_output_mod::output",4(4i4,2x))', l, u, self%dims, shape(data)
  if (.not.associated(self%buffer)) then
    allocate (self%buffer(self%dims(1),self%dims(2),self%dims(3),self%dims(4)))
  end if
  self%buffer(l(1):u(1),l(2):u(2),l(3):u(3),l(4):u(4)) = data
  self%count = self%count+1
  if (self%count==self%nbuf) call write_out (self)
  call trace_end
END SUBROUTINE output

!> =============================================================================
!> =============================================================================
SUBROUTINE close (self)
  class(direct_output_t):: self
  !-----------------------------------------------------------------------------
  call trace_begin('direct_output_mod::close')
  deallocate (self%buffer)
  nullify (self%buffer)
  call trace_end
END SUBROUTINE close

!> =============================================================================
!> Write out the buffer and reset the counter
!> =============================================================================
SUBROUTINE write_out (self, rec)
  class (direct_output_t) :: self
  integer, optional       :: rec
  integer                 :: l_rec, recl, i, j, i_rec
  integer(8)              :: recl8
  real(8)                 :: wt
  !-----------------------------------------------------------------------------
  call trace_begin('direct_output_mod::write_out')
  !$omp critical (direct_cr)
  self%locked = .true.
  if (present(rec)) then
    l_rec = rec
  else
    l_rec = 1
  end if
  self%rec = l_rec
  wt = wallclock()
  if (sizeof(self%buffer) < 2e9) then
    recl = 4*product(self%dims)
    open (io_unit%direct, file=trim(self%filename), access='direct', &
          status='unknown', recl=recl)
    write (io_unit%direct, rec=l_rec) self%buffer
  else
    recl8 = 4_8*product(self%dims(1:3))
    if (recl8 < 2e9) then
      recl = recl8
      open (io_unit%direct, file=trim(self%filename), access='direct', &
            status='unknown', recl=recl)
      do i=1,self%dims(4)
        write (io_unit%direct, rec=i+(l_rec-1)*self%dims(4)) self%buffer(:,:,:,i)
      end do
    else
      print *, 'DOUBLE LOOP'
      recl = 4*product(self%dims(1:2))
      open (io_unit%direct, file=trim(self%filename), access='direct', &
            status='unknown', recl=recl)
      do j=1,self%dims(4)
      do i=1,self%dims(3)
        i_rec = i + (j-1)*self%dims(3) + (l_rec-1)*self%dims(4)
        write (io_unit%direct, rec=i_rec) self%buffer(:,:,i,j)
      end do
      end do
    end if
  end if
  wt = max(wallclock()-wt,1d-9)
  self%gb = 4.0*product(self%dims)/1024.**3
  self%gbs = self%gb/wt
  self%count = 0
  close (io_unit%direct)
  self%locked = .false.
 !$omp end critical (direct_cr)
  call trace_end
END SUBROUTINE write_out

!> =============================================================================
!> Test output of 1000 chunks from a single thread
!> =============================================================================
SUBROUTINE test (self)
  class (direct_output_t) :: self
  integer, dimension(4)   :: sz, n, offset
  integer                 :: nbuf, ix, iy, iz, it, n1, n2, n3, n4
  real, pointer           :: d(:,:,:,:)
  !-----------------------------------------------------------------------------
  print *, 'direct_output_mod::test'
  sz = [100,100,100,5]
  n  = [10,10,10,5]
  allocate (d(n(1),n(2),n(3),n(4)))
  nbuf = product(sz)/product(n)
  call self%init ('data/test.dat', sz, nbuf)
  do iz=1,n(3)
    do iy=1,n(2)
      do ix=1,n(1)
        d(:,:,:,:) = (ix-1)+(iy-1)*10+(iz-1)*100
        offset = 1 + [(ix-1)*n(1),(iy-1)*n(2),(iz-1)*n(3),0]
        call self%patch (d, offset)
      end do
    end do
  end do
  call self%close
  close (io_unit%direct, status='delete')
END SUBROUTINE

END MODULE direct_output_mod

