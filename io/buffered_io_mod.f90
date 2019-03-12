!===============================================================================
!> $Id$
!> This module handles buffering and output of data to disk.  It receives
!> call of the type "call output%buffer (id, it, data)", where id is the task
!> id, iout is the output number, and data is a pointer to a one-dimensional
!> buffer. The actual content of the buffers has no relevance for this module,
!> and is left to the procedure that packs the information into the buffer.
!> For simplicity, the buffer size (ndata) is assumed to remain constant.
!===============================================================================
MODULE buffered_io_mod
  USE io_mod
  USE patch_mod
  USE trace_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! The dump data type holds a pointer to the data to be dumped to disk, with id
  ! the task id, it the time slot, and buffer the packed task info.  It is
  ! assumed that the caller packed data into a temporarily allocated buffer,
  ! which is deallocated once the data has been written to disk.
  !-----------------------------------------------------------------------------
  type, public:: dump_t
    class(dump_t), pointer:: next => null()
    integer:: id, iout
    real(8):: time
    integer:: ndata
    real(4), dimension(:), pointer:: data
  end type
  !-----------------------------------------------------------------------------
  ! The output object contains procedures for opening, buffering, writing, and
  ! closing buffered output to a direct access file dumps.dat, with index data
  ! saved to a separate index.dat file.
  !-----------------------------------------------------------------------------
  type, public:: buffered_io_t
    integer:: ndata
    integer:: iout = -1
    integer(8):: record  = 0
    integer:: n_buffer   = 0
    integer:: max_buffer = 10
    class(dump_t), pointer:: head  => null()
    class(dump_t), pointer:: tail  => null()
  contains
    procedure:: open
    procedure, nopass:: output
    procedure:: buffer
    procedure:: check
    procedure:: write
    procedure:: close
    procedure:: find
    procedure:: test
  end type
  integer(8):: revision=1
  type(buffered_io_t), public:: buffered_io
CONTAINS

!===============================================================================
!> Write results to disk, using the buffered I/O in output_mod
!===============================================================================
SUBROUTINE output (patch)
  class(patch_t):: patch
  !.............................................................................
  type(header_t):: header
  integer:: n_buf, n_data, ibuf, it1, it2, iv
  real(4), dimension(:), pointer:: buffer
  real, dimension(:,:,:), pointer:: var
  real:: pt
  !-----------------------------------------------------------------------------
  ! Compute size of buffer and allocate -- it will be freed by the writer
  !-----------------------------------------------------------------------------
  call trace_begin ('patch_t%output_new')
  !$omp critical (buffered_cr)
  n_buf = product(patch%gn)
  n_data = n_header + patch%nv*n_buf
  allocate (buffer(n_data), var(patch%gn(1),patch%gn(2),patch%gn(3)))
  !-----------------------------------------------------------------------------
  ! Copy relevant patch info to sequenced header, and copy that to the buffer
  !-----------------------------------------------------------------------------
  call patch%patch_to_header (header)
  ibuf = 1
  call anonymous_copy (n_header, header, buffer(ibuf))
  ibuf = ibuf + n_header
  !-----------------------------------------------------------------------------
  ! Copy the variables to the output buffer; these are first interpolated in time
  !-----------------------------------------------------------------------------
  it1 = patch%iit(patch%nt-2)
  it2 = patch%iit(patch%nt-1)
  pt = (patch%out_next-patch%t(it1))/max(patch%t(it2)-patch%t(it1),1d-30)
  do iv=1,patch%nv
    var = patch%mem(:,:,:,iv,it1,1)*(1.-pt) + patch%mem(:,:,:,iv,it2,1)*pt
    call anonymous_copy (n_buf, var, buffer(ibuf))
    ibuf = ibuf + n_buf
  end do
  deallocate (var)
  call buffered_io%buffer (patch%id, patch%iout, patch%out_next, buffer)
  call buffered_io%check
  !$omp end critical (buffered_cr)
  call trace_end
END SUBROUTINE output

!===============================================================================
!> Open an index file and a data file.  The first record of the index file starts
!> with the data format revision number and the buffer size.
!===============================================================================
SUBROUTINE open (self, nwords, iout)
  class(buffered_io_t):: self
  integer:: nwords, iout
  logical, save:: first_time=.true., exist
  character(len=64):: filename
  !-----------------------------------------------------------------------------
  call trace_begin('buffered_io_mod::open')
  self%iout= iout
  self%ndata = nwords
  write (filename,'(a,i5.5,"/buffered")') trim(io%outputname), iout
  if (io%verbose>0) &
    print'(a)', 'OPEN: '//trim(filename)
  !---------------------------------------------------------------------------
  ! Check if the index file exists.  If it does, read the record counter, and
  ! if not, reset the counter and write it out.
  !---------------------------------------------------------------------------
  if (io%verbose>1) &
    print*, 'index file: ', trim(filename)//'.idx'
  inquire (file=trim(filename)//'.idx', exist=exist)
  if (exist) then
    open (io_unit%index, file=trim(filename)//'.idx', &
      form='unformatted', access='direct', status='unknown', recl=2*io%word_size)
    read (io_unit%index, rec=3) self%record
  else
    open (io_unit%index, file=trim(filename)//'.idx', &
      form='unformatted', access='direct', status='unknown', recl=2*io%word_size)
    self%record = 0
    write (io_unit%index, rec=3) self%record
  end if
  write (io_unit%index, rec=1) revision
  write (io_unit%index, rec=2) int(nwords,kind=8)
  !---------------------------------------------------------------------------
  ! Open the data file for writing
  !---------------------------------------------------------------------------
  if (io%verbose>1) &
    print*, ' dump file: ', trim(filename)//'.dat'
  open (io_unit%dump,  file=trim(filename)//'.dat', &
    form='unformatted', access='direct', status='unknown', recl=nwords*io%word_size)
  call trace_end
END SUBROUTINE open

!===============================================================================
!> Append a dump instance to a linked list.
!===============================================================================
SUBROUTINE buffer (self, id, iout, time, data)
  class(buffered_io_t):: self
  integer:: id, iout
  real(8):: time
  real(4), dimension(:), pointer:: data
  class(dump_t), pointer:: dump
  !-----------------------------------------------------------------------------
  call trace_begin('buffered_io_mod::buffer')
  allocate (dump)
  dump%id = id
  dump%iout = iout
  dump%time = time
  dump%ndata = size(data)
  dump%data => data
  if (associated(self%tail)) then
    self%tail%next => dump
  else
    self%head => dump
  end if
  self%tail => dump
  self%n_buffer = self%n_buffer+1
  if (io%verbose>1) &
    print*,'buffered_io_mod::buffer: id, iout, n_buffer=', id, iout, self%n_buffer
  call trace_end
  call self%check
END SUBROUTINE buffer

!===============================================================================
!> Check if the number of dumps has reached the maximum => write dumps to disk.
!> This is done as a task-list independent step by all threads.
!===============================================================================
SUBROUTINE check (self)
  class(buffered_io_t):: self
  class(dump_t), pointer:: dump
  !-----------------------------------------------------------------------------
  if (self%n_buffer >= self%max_buffer) then
    call self%write
  end if
END SUBROUTINE check

!===============================================================================
!> Write the buffers in the linked list to disk and deallocate the buffers.
!> The index file records, sequentially, the task id and dump number of each
!> dump written.  By reading the entire index file one can search very rapidly
!> for a specific dump, by treating the id+iout combination as an 8 byte integer
!===============================================================================
SUBROUTINE write (self)
  class (buffered_io_t):: self
  class (dump_t), pointer:: dump, old
  !-----------------------------------------------------------------------------
  call trace_begin('buffered_io_mod::write')
  dump => self%head
  if (io%verbose>0) &
    print'(a,i7,i6,i9,g15.6)', 'output: it, nbuf, ndata, time', &
        dump%id, self%n_buffer, dump%ndata, dump%time
  do while (associated(dump))
    !---------------------------------------------------------------------------
    ! If time slot differs, close the previous one and open the current one
    !---------------------------------------------------------------------------
    if (dump%iout /= self%iout) then
      if (self%iout >= 0) call self%close
      call self%open (dump%ndata, dump%iout)
    end if
    !---------------------------------------------------------------------------
    ! Increment the record counter, record its value, and write the data
    !---------------------------------------------------------------------------
    !$omp atomic
    self%record = self%record+1
    write (io_unit%index, rec=3) self%record
    write (io_unit%index, rec=self%record+3) dump%id, dump%iout
    flush (io_unit%index)
    write (io_unit%dump, rec=self%record) dump%data
    flush (io_unit%dump)
    !$omp atomic
    self%n_buffer = self%n_buffer-1
    if (io%verbose>1) &
      print'(a,i7,i6,i9,g15.6)', 'output: id, it, ndata, time', &
        dump%id, dump%iout, dump%ndata, dump%time
    !---------------------------------------------------------------------------
    ! Delete the buffer data, and the buffer
    !---------------------------------------------------------------------------
    old => dump
    dump => dump%next
    deallocate (old%data)
    deallocate (old)
  end do
  nullify (self%head)
  nullify (self%tail)
  call trace_end
END SUBROUTINE write

!===============================================================================
!> Close the files, writing out the total number of records in slot 3 of the
!> index file.
!===============================================================================
SUBROUTINE close (self)
  class (buffered_io_t):: self
  !-----------------------------------------------------------------------------
  if (io%do_legacy) return
  call trace_begin('buffered_io_mod::close')
  close (io_unit%index)
  close (io_unit%dump)
  call trace_end
END SUBROUTINE close

!===============================================================================
!> Find the record that contains the combination of task id and dump number iout.
!> This will typically be done in Python or IDL; this is just a demo.
!===============================================================================
FUNCTION find (self, id, iout)
  class (buffered_io_t):: self
  integer, optional:: id, iout
  integer:: find, loc(1)
  integer(8), dimension(:), pointer:: index
  integer(8):: key, ndata, nrecord, revision, i
  !-----------------------------------------------------------------------------
  call trace_begin('buffered_io_mod::find')
  !-----------------------------------------------------------------------------
  ! Read the number of records and the buffer size
  !-----------------------------------------------------------------------------
  open (io_unit%index, file=trim(io%outputname)//'index.dat', &
    form='unformatted', access='direct', status='old', recl=8)
  read (io_unit%index,rec=1) revision
  read (io_unit%index,rec=2) ndata
  read (io_unit%index,rec=3) nrecord
  print *, revision, ndata, nrecord
  close (io_unit%index)
  !-----------------------------------------------------------------------------
  ! Get the entire index in one read
  !-----------------------------------------------------------------------------
  allocate (index(nrecord))
  open (io_unit%index, file=trim(io%outputname)//'index.dat', &
    form='unformatted', access='direct', status='old', recl=8*(nrecord+3))
  read (io_unit%index,rec=1) revision, ndata, nrecord, index
  close (io_unit%index)
  !-----------------------------------------------------------------------------
  ! Find the id+iout combo
  !-----------------------------------------------------------------------------
  key = int(id,kind=8) + int(iout,kind=8)*2_8**32
  if (io%verbose>1) &
    print*, 'searching for id, iout, key =', id, iout, key
  if (io%verbose>2) then
    do i=1,nrecord
      print *, i, index(i)
    end do
  end if
  loc = minloc(abs(index-key))
  find = loc(1)
  if (io%verbose>1) &
    print*, 'record =', find
  deallocate (index)
  call trace_end
END FUNCTION find

!===============================================================================
!> Test the output module, by buffering enough dumps for a self-triggered write,
!> followed by a write triggered by a close.
!===============================================================================
SUBROUTINE test (self)
  class (buffered_io_t):: self
  real(4), dimension(:), pointer:: data
  integer, dimension(:), pointer:: p
  integer:: n=20, nv=8, id, it, ndata, record
  real(8):: time
  !-----------------------------------------------------------------------------
  call trace_begin ('buffered_io_mod::test')
  ndata = n**3*nv+n**2
  call self%open (ndata,0)
  do it=1,3
    time = it-1d0
    do id=1,9
      allocate (data(ndata))
      call self%buffer (id, it, time, data)
    end do
  end do
  call self%close
  record = self%find (id=2,iout=3)
  print *,'record =', record
  call trace_end
END SUBROUTINE test

END MODULE buffered_io_mod
