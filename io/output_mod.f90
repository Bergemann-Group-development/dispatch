!> =============================================================================
!> $Id$
!> This module handles buffering and output of data to disk.  It receives
!> call of the type "call output%buffer (id, it, data)", where id is the task
!> id, iout is the output number, and data is a pointer to a one-dimensional
!> buffer. The actual content of the buffers has no relevance for this module,
!> and is left to the procedure that packs the information into the buffer.
!> For simplicity, the buffer size (ndata) is assumed to remain constant.
!> =============================================================================
MODULE output_mod
  USE io_mod
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
    type(dump_t), pointer:: next => null()
    integer:: id, iout
    integer:: ndata
    integer, dimension(:), pointer:: data
  end type
  !-----------------------------------------------------------------------------
  ! The output object contains procedures for opening, buffering, writing, and
  ! closing buffered output to a direct access file dumps.dat, with index data
  ! saved to a separate index.dat file.
  !-----------------------------------------------------------------------------
  type, public:: output_t
    integer:: ndata
    integer(8):: record  = 0
    integer:: n_buffer   = 0
    integer:: max_buffer = 5
    type(dump_t), pointer:: head  => null()
    type(dump_t), pointer:: tail  => null()
  contains
    procedure:: open
    procedure:: buffer
    procedure:: check
    procedure:: write
    procedure:: close
    procedure:: find
    procedure:: test
  end type
  integer(8):: revision=1
  type(output_t), public:: output
CONTAINS

!> =============================================================================
!> Open an index file and a data file.  The first record of the index file starts
!> with the data format revision number and the buffer size.
!> =============================================================================
SUBROUTINE open (self, nwords)
  class(output_t):: self
  integer:: nwords
  !-----------------------------------------------------------------------------
  call trace_begin('output_mod::open')
  self%ndata = nwords
  if (io%verbose>0) &
    print*,'index file: ', trim(io%datadir)//io%sep//'index.dat'
  open (io%index_unit, file=trim(io%datadir)//io%sep//'index.dat', &
    form='unformatted', access='direct', status='unknown', recl=2*io%word_size)
  if (io%verbose>0) &
    print*,' dump file: ', trim(io%datadir)//io%sep//'dumps.dat'
  open (io%dump_unit, file=trim(io%datadir)//io%sep//'dumps.dat', &
    form='unformatted', access='direct', status='unknown', recl=nwords*io%word_size)
  write (io%index_unit, rec=1) revision
  write (io%index_unit, rec=2) int(nwords,kind=8)
  self%record = 0
  call trace_end
END SUBROUTINE open

!> =============================================================================
!> Append a dump instance to a linked list.
!> =============================================================================
SUBROUTINE buffer (self, id, iout, data)
  class(output_t):: self
  integer:: id, iout
  integer, dimension(:), pointer:: data
  class(dump_t), pointer:: dump
  !-----------------------------------------------------------------------------
  call trace_begin('output_mod::buffer')
  allocate (dump)
  dump%id = id
  dump%iout = iout
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
    print*,'output_mod::buffer: id, iout, n_buffer=', id, iout, self%n_buffer
  call trace_end
END SUBROUTINE buffer

!> =============================================================================
!> Check if the number of dumps has reached the maximum => write dumps to disk.
!> This is done as a task-list independent step by all threads.
!> =============================================================================
SUBROUTINE check (self)
  class(output_t):: self
  class(dump_t), pointer:: dump
  !-----------------------------------------------------------------------------
  if (self%n_buffer >= self%max_buffer) then
    call self%write
  end if
END SUBROUTINE check

!> =============================================================================
!> Write the buffers in the linked list to disk and deallocate the buffers.
!> The index file records, sequentially, the task id and dump number of each
!> dump written.  By reading the entire index file one can search very rapidly
!> for a specific dump, by treating the id+iout combination as an 8 byte integer
!> =============================================================================
SUBROUTINE write (self)
  class (output_t):: self
  class (dump_t), pointer:: dump, old
  !-----------------------------------------------------------------------------
  call trace_begin('output_mod::write')
  dump => self%head
  do while (associated(dump))
    self%record = self%record+1
    write (io%index_unit,  rec=self%record+3) dump%id, dump%iout
    flush (io%index_unit)
    write (io%dump_unit, rec=self%record) dump%data
    flush (io%dump_unit)
    !$omp atomic
    self%n_buffer = self%n_buffer-1
    if (io%verbose>0) print*,'wrote id, it, ndata', dump%id, dump%iout, dump%ndata
    old => dump
    dump => dump%next
    deallocate (old%data)
    deallocate (old)
  end do
  nullify (self%head)
  nullify (self%tail)
  call trace_end
END SUBROUTINE write

!> =============================================================================
!> Close the files, writing out the total number of records in slot 3 of the
!> index file.
!> =============================================================================
SUBROUTINE close (self)
  class (output_t):: self
  !-----------------------------------------------------------------------------
  call trace_begin('output_mod::close')
  call self%write
  write (io%index_unit,rec=3) self%record
  close (io%index_unit)
  close (io%dump_unit)
  call trace_end
END SUBROUTINE close

!> =============================================================================
!> Find the record that contains the combination of task id and dump number iout.
!> This will typically be done in Python or IDL; this is just a demo.
!> =============================================================================
FUNCTION find (self, id, iout)
  class (output_t):: self
  integer, optional:: id, iout
  integer:: find, loc(1)
  integer(8), dimension(:), pointer:: index
  integer(8):: key, ndata, nrecord, revision, i
  !-----------------------------------------------------------------------------
  call trace_begin('output_mod::find')
  !-----------------------------------------------------------------------------
  ! Read the number of records and the buffer size
  !-----------------------------------------------------------------------------
  open (io%index_unit, file=trim(io%datadir)//io%sep//'index.dat', &
    form='unformatted', access='direct', status='old', recl=8)
  read (io%index_unit,rec=1) revision
  read (io%index_unit,rec=2) ndata
  read (io%index_unit,rec=3) nrecord
  print *, revision, ndata, nrecord
  close (io%index_unit)
  !-----------------------------------------------------------------------------
  ! Get the entire index in one read
  !-----------------------------------------------------------------------------
  allocate (index(nrecord))
  open (io%index_unit, file=trim(io%datadir)//io%sep//'index.dat', &
    form='unformatted', access='direct', status='old', recl=8*(nrecord+3))
  read (io%index_unit,rec=1) revision, ndata, nrecord, index
  close (io%index_unit)
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

!> =============================================================================
!> Test the output module, by buffering enough dumps for a self-triggered write,
!> followed by a write triggered by a close.
!> =============================================================================
SUBROUTINE test (self)
  class (output_t):: self
  integer, dimension(:), pointer:: data
  integer, dimension(:), pointer:: p
  integer:: n=20, nv=8, id, it, ndata, record
  !-----------------------------------------------------------------------------
  call trace_begin ('output_mod::test')
  ndata = n**3*nv+n**2
  call self%open (ndata)
  do it=1,3
    do id=1,9
      allocate (data(ndata))
      call self%buffer (id, it, data)
    end do
  end do
  call self%close
  record = self%find (id=2,iout=3)
  print *,'record =', record
  call trace_end
END SUBROUTINE test

END MODULE output_mod
