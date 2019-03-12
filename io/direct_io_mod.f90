!===============================================================================
!> $Id$
!===============================================================================
MODULE direct_io_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE patch_mod
  USE omp_timer_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! The output object contains procedures for opening, buffering, writing, and
  ! closing buffered output to a direct access file.
  !-----------------------------------------------------------------------------
  type, public:: direct_io_t
    integer:: id
    integer:: dims(4)
    integer:: count=0
    integer:: nbuf
    integer:: rec=0
    real, dimension(:,:,:,:), allocatable:: buffer
    character(len=128):: filename
    logical:: locked=.false.
    real(8):: position(3)=0d0
    real:: gb=0.0, gbs=0.0
  contains
    procedure:: init
    procedure:: output
    procedure:: input
    procedure:: out
    procedure:: in
    procedure:: close
  end type
  integer(8):: revision=1
  integer:: verbose=1
  type(direct_io_t), public:: direct_io
CONTAINS

!===============================================================================
!> The unigrid size is available from the mesh box and patch sizes
!===============================================================================
SUBROUTINE init (self, patch, filename)
  class(direct_io_t)         :: self
  class(patch_t):: patch
  character(len=*), optional :: filename
  integer:: mpatch(3)
  !-----------------------------------------------------------------------------
  call trace_begin('direct_io_t%init')
  !$omp critical (direct_cr)
  !
  if (present(filename)) then
    if (verbose>1) print *, 'direct_io_t: filename = ', trim(filename)
    self%filename = filename
  end if
  !----------------------------------------------------------------------------
  ! Cartesian dimensions of patch arrangement and unigrid dimensions self%dims
  !----------------------------------------------------------------------------
  mpatch = nint(patch%mesh%b/patch%mesh%s)
  if (patch%no_mans_land) then
    self%dims(1:3) = mpatch*patch%mesh%n
  else
    self%dims(1:3) = mpatch*(patch%mesh%n-1)+1
  end if
  self%id = patch%id
  self%dims(4) = patch%nv
  self%nbuf = product(mpatch)
  !
  self%count = 0
  !$omp end critical (direct_cr)
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Output a patch to a direct access unigrid data file
!===============================================================================
SUBROUTINE output (self, patch)
  class(direct_io_t):: self
  class(patch_t):: patch
  !-----------------------------------------------------------------------------
  character(len=64)    :: filename
  real, pointer        :: buf(:,:,:,:)
  integer              :: n(3), l(3), u(3), ng(3)=0, jt(2), id=0, offset(4)=1
  real                 :: pt(2)
  !-----------------------------------------------------------------------------
  call trace_begin('patch_t%output_direct')
  !$omp critical (direct_cr)
  !-----------------------------------------------------------------------------
  ! Write a small text file with info when the last chunk is being entered
  !-----------------------------------------------------------------------------
  if (self%count==self%nbuf-1) then
    write (filename,'(a,i5.5,"/",i5.5)') trim(io%outputname), patch%iout, id
    if (patch%no_mans_land) then
      io%format = 1
    else
      io%format = 2
    end if
    open (io%data_unit, file=trim(filename)//'.txt', form='formatted', status='unknown')
    write (io%data_unit,'(i2.2)') io%format
    write (io%data_unit,'(a)') trim(patch%kind)
    write (io%data_unit,'(a)') trim(patch%eos)
    write (io%data_unit,'(a)') trim(patch%opacity)
    write (io%data_unit,'(1p,2e18.10,6(2x,3e18.10)1x,1x,2e10.3)') patch%out_next, patch%dtime, &
          self%position, patch%ds, patch%velocity, patch%llc_nat, patch%llc_cart, &
          patch%mesh%centre_nat, patch%quality, patch%gamma
    l = 1
    u = self%dims(1:3)
    n = u-l+1
    id = mpi%rank
    write (io%data_unit,'(6(3i10.1,2x),3i2.1,1x,1i2.1,i3)') int(patch%box/patch%ds+0.5), &
          n, l, u, n, n, ng, patch%nv, patch%level
    write (io%data_unit,'(3(1i8.1,1x),i2.1)') id, patch%iout, patch%istep, patch%mesh_type
    close (io%data_unit)
    self%filename = trim(filename)//'.dat'
  end if
  !-----------------------------------------------------------------------------
  ! Pass a chunk to the procedure that eventually writes a binary file
  !-----------------------------------------------------------------------------
  l = patch%mesh%li
  u = patch%mesh%ui
  n = u-l+1
  if (verbose>1) print *, 'output:', patch%id, l, u
  allocate (buf(n(1),n(2),n(3),patch%nv))
  call patch%time_interval (patch%out_next, jt, pt)
  buf = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,jt(1),1)*pt(1) &
      + patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,jt(2),1)*pt(2)
  offset(1:3) = 1 + patch%ipos*patch%ncell
  call self%out (buf, offset)
  deallocate (buf)
  !$omp end critical (direct_cr)
  call trace_end
END SUBROUTINE output

!===============================================================================
!> Output a patch to a direct access unigrid data file
!===============================================================================
SUBROUTINE input (self, patch, ok)
  class(direct_io_t):: self
  class(patch_t):: patch
  logical:: ok
  !-----------------------------------------------------------------------------
  character(len=64)    :: filename
  real, pointer        :: buf(:,:,:,:)
  integer              :: n(3), l(3), u(3), ng(3)=0, jt(2), id=0, offset(4)=1, iof
  real                 :: pt(2)
  logical              :: exist
  real(8), save        :: time
  !-----------------------------------------------------------------------------
  call trace_begin('patch_t%input_direct')
  !$omp critical (direct_in_cr)
  !-----------------------------------------------------------------------------
  ! Read a small text file with info, to get the snapshot time
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/",i5.5)') trim(io%outputname), patch%restart, id
  inquire (file=trim(filename)//'.dat', exist=exist)
  if (exist) then
    if (self%count==0) then
      open (io%data_unit, file=trim(filename)//'.txt', form='formatted', status='old')
      read (io%data_unit,'(i2.2)') iof
      read (io%data_unit,'(a)')
      read (io%data_unit,'(a)')
      read (io%data_unit,'(a)')
      read (io%data_unit,'(1p,2e18.10,6(2x,3e18.10)1x,1x,2e10.3)') time
      close (io%data_unit)
      self%filename = trim(filename)//'.dat'
    end if
    patch%time = time
    !-----------------------------------------------------------------------------
    ! Get a chunk from the procedure that reads a whole piece
    !-----------------------------------------------------------------------------
    l = patch%mesh%li
    u = patch%mesh%ui
    n = u-l+1
    allocate (buf(n(1),n(2),n(3),patch%nv))
    !call patch%time_interval (patch%out_next, jt, pt)
    offset(1:3) = 1 + patch%ipos*patch%ncell
    call self%in (buf, offset)
    patch%mem(    :    ,    :    ,    :    ,:,patch%it,1) = 0.0
    patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,patch%it,1) = buf
    if (verbose>1) print *, 'id, dmin, dmax', patch%id, minval(buf(:,:,:,1)), maxval(buf(:,:,:,1)), l, u
    if (verbose>1) print *, 'id, dmin, dmax', patch%id, patch%fminval(patch%idx%d), patch%fmaxval(patch%idx%d), shape(buf)
    deallocate (buf)
    ok = .true.
    if (self%count==self%nbuf) then
      call self%close
    end if
  else
    ok = .false.
  end if
  if (io%verbose > 0) &
    print *, 'input_direct:', trim(filename)//'.dat', exist, ok, self%count
  !$omp end critical (direct_in_cr)
  call trace_end
END SUBROUTINE input

!===============================================================================
!> Interface to optionally enforce a critical region
!===============================================================================
SUBROUTINE out (self, data, offset)
  class (direct_io_t):: self
  real, dimension(:,:,:,:) :: data
  integer, dimension(4) :: offset
  !-----------------------------------------------------------------------------
  call trace_begin('direct_io_t%out')
  if (self%locked) then
    call out_real (self, data, offset)
  else
    call out_real (self, data, offset)
  end if
  call trace_end
END SUBROUTINE out

!===============================================================================
!> Buffer one patch, accumulating until the buffer is complete
!===============================================================================
SUBROUTINE out_real (self, data, offset)
  class (direct_io_t):: self
  integer, dimension(4) :: offset, n, l, u 
  real, dimension(:,:,:,:) :: data
  !-----------------------------------------------------------------------------
  call trace_begin('direct_io_t%out_real')
  n = shape(data)
  l = offset
  u = l+n-1
  if (io%debug(3)) &
    print '("direct_io_t%out",4(4i4,2x))', l, u, self%dims, shape(data)
  if (.not.allocated(self%buffer)) then
    allocate (self%buffer(self%dims(1),self%dims(2),self%dims(3),self%dims(4)))
  end if
  self%buffer(l(1):u(1),l(2):u(2),l(3):u(3),l(4):u(4)) = data
  self%count = self%count+1
  if (self%count==self%nbuf) call write_out (self)
  call trace_end
END SUBROUTINE out_real

!===============================================================================
!> Interface to optionally enforce a critical region
!===============================================================================
SUBROUTINE in (self, data, offset)
  class (direct_io_t):: self
  real, dimension(:,:,:,:) :: data
  integer, dimension(4) :: offset
  !-----------------------------------------------------------------------------
  call trace_begin('direct_io_t%in')
  if (self%locked) then
    call in_real (self, data, offset)
  else
    call in_real (self, data, offset)
  end if
  call trace_end
END SUBROUTINE in

!===============================================================================
!> Buffer one patch, accumulating until the buffer is complete
!===============================================================================
SUBROUTINE in_real (self, data, offset)
  class (direct_io_t):: self
  integer, dimension(4) :: offset, n, l, u 
  real, dimension(:,:,:,:) :: data
  !-----------------------------------------------------------------------------
  call trace_begin('direct_io_t%in_real')
  n = shape(data)
  l = offset
  u = l+n-1
  if (io%debug(3)) &
    print '("direct_io_t%in",4(4i4,2x))', l, u, self%dims, shape(data)
  if (.not.allocated(self%buffer)) then
    allocate (self%buffer(self%dims(1),self%dims(2),self%dims(3),self%dims(4)))
  end if
  if (self%count==0) call read_in (self)
  data = self%buffer(l(1):u(1),l(2):u(2),l(3):u(3),l(4):u(4))
  self%count = self%count+1
  call trace_end
END SUBROUTINE in_real

!===============================================================================
!> Close an output = deallocate buffer
!===============================================================================
SUBROUTINE close (self)
  class(direct_io_t):: self
  !-----------------------------------------------------------------------------
  call trace_begin('direct_io_t%close')
  if (allocated(self%buffer)) then
    if (verbose>1) print *, 'deallocating buffer'
    deallocate (self%buffer)
  end if
  call trace_end
END SUBROUTINE close

!===============================================================================
!> Write out the buffer and reset the counter
!===============================================================================
SUBROUTINE write_out (self, rec)
  class (direct_io_t) :: self
  integer, optional       :: rec
  integer                 :: l_rec, recl, i, j, i_rec
  integer(8)              :: recl8
  real(8)                 :: wt
  !-----------------------------------------------------------------------------
  call trace_begin('direct_io_t%write_out')
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
      if (verbose>1) print *, 'DOUBLE LOOP'
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
  call trace_end
END SUBROUTINE write_out

!===============================================================================
!> Read in the buffer and reset the counter
!===============================================================================
SUBROUTINE read_in (self, rec)
  class (direct_io_t) :: self
  integer, optional       :: rec
  integer                 :: l_rec, recl, i, j, i_rec
  integer(8)              :: recl8
  real(8)                 :: wt
  !-----------------------------------------------------------------------------
  call trace_begin('direct_io_t%read_in')
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
          status='old', recl=recl)
    read (io_unit%direct, rec=l_rec) self%buffer
  else
    recl8 = 4_8*product(self%dims(1:3))
    if (recl8 < 2e9) then
      recl = recl8
      open (io_unit%direct, file=trim(self%filename), access='direct', &
            status='old', recl=recl)
      do i=1,self%dims(4)
        read (io_unit%direct, rec=i+(l_rec-1)*self%dims(4)) self%buffer(:,:,:,i)
      end do
    else
      if (verbose>1) print *, 'DOUBLE LOOP'
      recl = 4*product(self%dims(1:2))
      open (io_unit%direct, file=trim(self%filename), access='direct', &
            status='old', recl=recl)
      do j=1,self%dims(4)
      do i=1,self%dims(3)
        i_rec = i + (j-1)*self%dims(3) + (l_rec-1)*self%dims(4)
        read (io_unit%direct, rec=i_rec) self%buffer(:,:,i,j)
      end do
      end do
    end if
  end if
  print *, 'read_in:', shape(self%buffer), minval(self%buffer), maxval(self%buffer)
  wt = max(wallclock()-wt,1d-9)
  self%gb = 4.0*product(self%dims)/1024.**3
  self%gbs = self%gb/wt
  self%count = 0
  close (io_unit%direct)
  self%locked = .false.
  call trace_end
END SUBROUTINE read_in

END MODULE direct_io_mod

