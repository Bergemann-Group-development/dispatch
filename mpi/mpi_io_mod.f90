!*******************************************************************************
!> Module for handling blocking and non-blocking MPI parallel I/O to a single file.
!> The module is initialized by specifying the chunk size (which can be 1- or
!> 3-dimensional). A chunk could be just a patch, or a single chunk of data per
!> rank.
!>
!> The module can be used to handle any collection of patches, by using the
!> patch%id number to map each patch into the file.  If the patch sizes are fixed,
!> this is a unique mapping, and if each snapshot is written into a separate file
!> this is also sufficiently general.  The header part of the file may be written
!> to a separate file.
!>
!> Typical multi-threaded use:
!>
!>   type(mpi_file_t):: file                         ! holds file info
!>   type(mpi_io_t):: thread_io                         ! holds thread info
!>   !$omp threadprivate (thread_io)                    ! threadprivate data type
!>   ...
!>   call file%openw ('some_file.dat')                  ! implicit MPI barrier
!>   call thread_io%use (file)                          ! all threads use the same file
!>   ...
!>   call thread_io%init (snapshot_size, patch_size)    ! define offset factors
!>   call thread_io%write (patch_data, iout1, id1)      ! write snapshot iout1, patch id1
!>   ...
!>   call thread_io%read (patch_data, iout2, id2)       ! read snapshot iout2, patch id2
!>   ...
!>   !$omp barrier                                      ! all threads must agree
!>   call file%close                                    ! implicit MPI barrier
!>
!*******************************************************************************
MODULE mpi_io_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE trace_mod
  USE mpi_mod, only: mp=>mpi
  USE mpi_file_mod
  USE dll_mod
#ifdef MPI
  USE mpi, only: MPI_OFFSET_KIND, MPI_REAL, MPI_STATUS_SIZE, MPI_THREAD_MULTIPLE
#endif MPI
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! dll_node extension with a request and buffer pointer
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: iwrite_t
    integer:: req, handle, words
#ifdef MPI
    integer(kind=MPI_OFFSET_KIND):: pos
#else
    integer(8):: pos
#endif
    real, dimension(:,:,:), pointer:: buffer3 => null()
    real, dimension(:,:,:,:), pointer:: buffer => null()
  end type
  !-----------------------------------------------------------------------------
  ! dll_t extension with a procedure for testing request completion
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: iwrite_list_t
    logical:: first_time=.true.
    logical:: active=.false.
    integer:: thread=-1
  contains
    procedure:: check
  end type
  !-----------------------------------------------------------------------------
  ! MPI I/O data type
  !-----------------------------------------------------------------------------
  type, public:: mpi_io_t
    type(mpi_file_t):: file
    integer:: req, nwrite=-1
    integer(8):: rec_words
    integer:: chunk_words
    integer:: err
    real, dimension(:,:,:), allocatable:: buf
    type(iwrite_list_t), public:: iwrite_list
#ifdef MPI
    integer:: status(MPI_STATUS_SIZE)
#else
    integer:: status
#endif
  contains
    procedure:: use
    procedure:: init
    procedure:: set
    procedure:: check_init
    procedure, private:: write3
    procedure, private:: write4
    generic, public :: write => write3, write4
    procedure, private:: iwrite3
    procedure, private:: iwrite4
    generic, public :: iwrite => iwrite3, iwrite4
    procedure, private:: read3
    procedure, private:: read4
    generic, public :: read => read3, read4
    procedure:: assert
  end type
  character(len=32), save:: fmt='(1x,a,4i8,1p,5e12.3)'
  integer, save:: verbose=0
  logical, save:: direct=.false.
  type(mpi_io_t), public:: mpi_io
CONTAINS

!===============================================================================
!> Initialize the local size self%lsize, global size self%gsize
!===============================================================================
SUBROUTINE use (self, file)
  class(mpi_io_t):: self
  type(mpi_file_t):: file
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%use')
  self%file = file
  call trace%end()
END SUBROUTINE use

!===============================================================================
!> Initialize the local size self%lsize, global size self%gsize
!===============================================================================
SUBROUTINE init (self)
  class(mpi_io_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /mpi_io_params/ verbose, direct
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%init')
  !-----------------------------------------------------------------------------
  ! Note that when this is called from mpi_io_t%init(), it is already in a
  ! critical (input_cr) region, so must either use a different name, or no
  ! critical region
  !-----------------------------------------------------------------------------
  !$omp critical (mpi_io_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, mpi_io_params, iostat=iostat)
    if (io%master) write (io%output, mpi_io_params)
  end if
  !$omp end critical (mpi_io_cr)
  !$omp critical (iwrite_cr)
  if (self%iwrite_list%first_time) then
    self%iwrite_list%first_time = .false.
    call self%iwrite_list%init ('iwrite_list')
    if (io%master) &
      print *, 'mpi_io_t%init done'
  end if
  !$omp end critical (iwrite_cr)
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize the local size self%lsize, global size self%gsize
!===============================================================================
SUBROUTINE set (self, rec_words, chunk_words, nwrite)
  class(mpi_io_t):: self
  integer(8), optional:: rec_words, chunk_words
  integer, optional:: nwrite
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%set')
  if (present(rec_words)) then
    self%rec_words = rec_words
    if (verbose > 2) print *, 'mpi_io_t%init: setting rec_words =', rec_words
  end if
  if (present(chunk_words)) then
    self%chunk_words = chunk_words
    if (verbose > 2) print *, 'mpi_io_t%init: setting chunk_words =', chunk_words
  end if
  if (present(nwrite)) then
    mpi_io%nwrite = nwrite
    if (verbose > 2) print *, 'mpi_io_t%init: setting nwrite =', nwrite
  end if
  call trace%end()
END SUBROUTINE set

!===============================================================================
!> Initialize the local size self%lsize, global size self%gsize
!===============================================================================
SUBROUTINE check_init (self)
  class(mpi_io_t):: self
  !-----------------------------------------------------------------------------
  call self%assert ('mpi_io_t:rec_words is not positive', err=self%rec_words <= 0)
  call self%assert ('mpi_io_t:chunk_words is not positive', err=self%chunk_words <= 0)
END SUBROUTINE check_init

!===============================================================================
!> Write a chunk (self%lsize) out to the open file
!===============================================================================
SUBROUTINE write3 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer:: rec, id
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%write')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (verbose>1) &
    write (io_unit%output, '(a,3i8,i15,1p,2e14.6)') ' mpi_io_t%write3: rec, id, words, pos =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_file_write_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                            self%status, self%err)
  else
    !$omp critical (mpi_cr)
    call MPI_file_write_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                            self%status, self%err)
    !$omp end critical (mpi_cr)
  end if
  call self%assert ('mpi_io_mod::write MPI_File_write_at '//trim(self%file%filename))
  call trace%end()
#endif
END SUBROUTINE write3

!===============================================================================
!> Write a chunk (self%lsize) out to the open file
!===============================================================================
SUBROUTINE write4 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:,:):: f
  integer:: rec, id
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%write')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (verbose > 1) &
    write (io_unit%output, '(a,2i6,i12,i15,1p,2e14.6)') &
      ' mpi_io_t%write4: rec, id, words, pos, min/max =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  if (verbose > 0) then
    write(io_unit%mpi,'(g12.5,i4,2x,a)') wallclock(), omp%thread, &
      'mpi_io_t%write4: writing to '//trim(self%file%filename)
    flush (io_unit%mpi)
  end if
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_file_write_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                            self%status, self%err)
  else
    !$omp critical (mpi_cr)
    call MPI_file_write_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                            self%status, self%err)
    !$omp end critical (mpi_cr)
  end if
  call self%assert ('mpi_io_mod::write MPI_File_write_at '//trim(self%file%filename))
  call trace%end()
#endif
END SUBROUTINE write4

!===============================================================================
!> Write a chunk (self%lsize) out to the open file
!===============================================================================
SUBROUTINE iwrite3 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:), pointer:: f
  integer:: rec, id, unit
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  class(iwrite_t), pointer:: item
  class(dll_node_t), pointer:: node
  integer:: n
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%iwrite', itimer=itimer)
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (verbose > 1) then
    if (verbose > 2) then
      unit = io_unit%output
    else
      unit = io_unit%mpi
    end if
    write (unit,'(a,3i8,i15,1p,2e14.6)') &
      ' mpi_io_t%read3: rec, id, words, pos =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  end if
  !-----------------------------------------------------------------------------
  ! Add details we need in order to check for buffers to deallocate
  !-----------------------------------------------------------------------------
  allocate (item)
  item%handle = self%file%handle
  item%words  = self%chunk_words
  item%pos    = pos
  item%req    = id
  item%buffer3 => f
  !$omp critical (iwrite_cr)
  node => item
  call self%iwrite_list%append (node)
  self%iwrite_list%active = .true.
  n = self%iwrite_list%n
  !$omp end critical (iwrite_cr)
  if (verbose>1) &
    write(io%output,'(f12.6,i5,2x,a,i7,i4)') &
      wallclock(), omp%thread, 'appended req, n =', item%req, n
  call trace%end (itimer)
#endif
END SUBROUTINE iwrite3

!===============================================================================
!> Write a chunk (self%lsize) out to the open file
!===============================================================================
SUBROUTINE iwrite4 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:,:), pointer:: f
  integer:: rec, id, unit
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  class(iwrite_t), pointer:: item
  class(dll_node_t), pointer:: node
  integer:: n
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%iwrite', itimer=itimer)
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (verbose > 1) then
    if (verbose > 2) then
      unit = io_unit%output
    else
      unit = io_unit%mpi
    end if
    write (unit,'(a,3i8,i15,1p,2e14.6)') &
      ' mpi_io_t%read3: rec, id, words, pos =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  end if
  !-----------------------------------------------------------------------------
  ! Add details we need in order to check for buffers to deallocate
  !-----------------------------------------------------------------------------
  allocate (item)
  item%handle = self%file%handle
  item%words  = self%chunk_words
  item%pos    = pos
  item%req    = id
  item%buffer => f
  !$omp critical (iwrite_cr)
  node => item
  call self%iwrite_list%append (node)
  self%iwrite_list%active = .true.
  n = self%iwrite_list%n
  !$omp end critical (iwrite_cr)
  if (verbose>1) &
    write(io%output,'(f12.6,i5,2x,a,i7,i4)') &
      wallclock(), omp%thread, 'appended req, n =', item%req, n
  call trace%end (itimer)
#endif
END SUBROUTINE iwrite4

!===============================================================================
!> Read a chunk (self%lsize) from the open file
!===============================================================================
SUBROUTINE read3 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer rec, id, unit
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  if (rec < 1) &
    call io%abort ('mpi_io_t%read4: illegal record number')
  call trace%begin ('mpi_io_t%read3')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_File_read_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                           self%status, self%err)
  else
    !$omp critical (mpi_cr)
    call MPI_File_read_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                           self%status, self%err)
    !$omp end critical (mpi_cr)
  end if
  if (verbose > 1) then
    if (verbose > 2) then
      unit = io_unit%output
    else
      unit = io_unit%mpi
    end if
    write (unit,'(a,3i8,i15,1p,2e14.6)') &
      ' mpi_io_t%read3: rec, id, words, pos =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  end if
  call trace%end()
#endif
END SUBROUTINE read3

!===============================================================================
!> Read a chunk (self%lsize) from the open file
!===============================================================================
SUBROUTINE read4 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:,:):: f
  integer rec, id, unit
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  if (rec < 1) &
    call io%abort ('mpi_io_t%read4: illegal record number')
  call trace%begin ('mpi_io_t%read4')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_File_read_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                         self%status, self%err)
  else
    !$omp critical (mpi_cr)
    call MPI_File_read_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                         self%status, self%err)
    !$omp end critical (mpi_cr)
  end if
  if (verbose > 1) then
    if (verbose > 2) then
      unit = io_unit%output
    else
      unit = io_unit%mpi
    end if
    write (unit,'(a,3i8,i15,1p,2e14.6)') &
      ' mpi_io_t%read4: rec, id, words, pos =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  end if
  call trace%end()
#endif
END SUBROUTINE read4

!===============================================================================
!> Assert no error
!===============================================================================
SUBROUTINE assert (self, label, err)
  class(mpi_io_t):: self
  character(len=*):: label
  logical, optional:: err
  !-----------------------------------------------------------------------------
  if (present(err)) then
    if (err) then
      call mp%abort (label)
    end if
  else if (self%err /= 0) then
    print*,mp%rank,omp%thread,' error =',self%err
    call mp%abort (label)
  end if
END SUBROUTINE assert

!===============================================================================
!> Check for queued I/O, and do it all with one single thread, to avoid
!> slowing down the other threads.  Other threads are free to append more
!> items at the end of the iwrite_list in the mean time -- this does not
!> disturb the thread doing I/O.
!===============================================================================
SUBROUTINE check (self)
  class(iwrite_list_t):: self
  class(dll_node_t), pointer:: node, next
  class(iwrite_t), pointer:: item
#ifdef MPI
  integer:: status(MPI_STATUS_SIZE), err
#else
  integer:: status(8), err
  integer:: MPI_REAL=1
#endif
  integer, save:: itimer=0, nprint=3
  integer:: done, omp_get_thread_num, rec
  logical:: ok
  real(8):: start, wc, dwc
  !-----------------------------------------------------------------------------
  ! Most threads will exit immediately here
  !-----------------------------------------------------------------------------
  if (.not.io%do_output .or. .not.self%active) &
    return
  if (self%n < mpi_io%nwrite .or. self%thread >=0) then
    if (verbose > 1 .and. self%thread >= 0) &
      write(io%output,'(f12.6,i5,2x,a,i7,i4)') &
        wallclock(), omp%thread, ' returning immediately', self%n, self%thread
    return
  end if
  call trace%begin ('iwrite_list_t%check', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! More than one thread may succeed getting to here, so we pick the first of
  ! them in a critical region -- the rest will find ok to be false and exit
  !-----------------------------------------------------------------------------
  !$omp critical (check_cr)
  ok = (self%thread == -1)
  if (ok) &
    self%thread = omp_get_thread_num()
  !$omp end critical (check_cr)
  !-----------------------------------------------------------------------------
  ! The thread that wins the race starts writing out data, while also allowing
  ! more buffers to be appended simultaneously
  !-----------------------------------------------------------------------------
  if (ok) then
    if (verbose > 0 .and. .not.io_unit%do_validate) &
      write(io%output,'(1p,g14.6,i6,3x,a,i6)') &
        wallclock(), self%thread, 'I/O start, n =', self%n
    start = wallclock()
    done = 0
    !!omp critical (iwrite_cr)
    node => self%head
    !!omp end critical (iwrite_cr)
    do while (self%n > 0 .and. associated(node))
      err = 0
      select type (node)
      class is (iwrite_t)
        item => node
        if (associated(item%buffer) .or. associated(item%buffer3)) then
          wc = wallclock()
          if (direct) then
            rec = 1 + item%pos/item%words/4
            write (io_unit%direct, rec=rec) item%buffer
            deallocate (item%buffer)
          else if (associated(item%buffer3)) then
#ifdef MPI
            if (mp%mode == MPI_THREAD_MULTIPLE) then
              call MPI_File_write_at (item%handle, item%pos, item%buffer3, &
                                      item%words, MPI_REAL, status, err)
            else
              !$omp critical (mpi_cr)
              call MPI_File_write_at (item%handle, item%pos, item%buffer3, &
                                      item%words, MPI_REAL, status, err)
              !$omp end critical (mpi_cr)
            end if
#endif
            call io%assert (err==0, 'MPI_File_write_at: error code non-zero')
            deallocate (item%buffer3)
            nullify (item%buffer3)
          else
#ifdef MPI
            if (mp%mode == MPI_THREAD_MULTIPLE) then
              call MPI_File_write_at (item%handle, item%pos, item%buffer, &
                                      item%words, MPI_REAL, status, err)
            else
              !$omp critical (mpi_cr)
              call MPI_File_write_at (item%handle, item%pos, item%buffer, &
                                      item%words, MPI_REAL, status, err)
              !$omp end critical (mpi_cr)
            end if
#endif
            call io%assert (err==0, 'MPI_File_write_at: error code non-zero')
            deallocate (item%buffer)
          end if
          nullify (item%buffer)
          !---------------------------------------------------------------------
          if (verbose > 0 .and. nprint > 0 .and. .not.io_unit%do_validate) then
            nprint = nprint-1
            dwc = wallclock()-wc
            if (direct) then
              write(io%output,'(10x,a,f12.6)') 'direct access write time:', dwc
            else
              write(io%output,'(10x,a,f12.6)') 'MPI_File_write_at time:', dwc
            end if
          end if
          !---------------------------------------------------------------------
          ! Spin for similar time, to yield to other MPI processes
          !---------------------------------------------------------------------
          !wc = wallclock()
          !do while (wallclock()-wc < dwc)
          !end do
        else
          write(io%output,*) &
            self%thread, 'buffer not associated on req', item%req
        end if
      end select
      done = done+1
      if (.not.io_unit%do_validate) then
        if (verbose > 1) &
          write(io%output,'(1p,g14.6,i6,3x,a,2i7)') &
            wallclock(), self%thread, 'I/O req =', item%req, done
        if (verbose > 0 .and. mod(done,1000) == 0) &
          write(io%output,'(1p,g14.6,i6,3x,a,i7)') &
            wallclock(), self%thread, 'I/O done =', done
      end if
      !!omp critical (iwrite_cr)
      next => node%next
      call self%remove (node)
      deallocate (node)
      node => next
      !!omp end critical (iwrite_cr)
    end do
    if (verbose > 0 .and. .not.io_unit%do_validate) &
      write(io%output,'(1p,g14.6,i6,3x,a,i6,0p,f12.3," s")') &
        wallclock(), self%thread, 'I/O final, n =', done, wallclock()-start
    !!omp critical (iwrite_cr)
    self%thread = -1
    !!omp end critical (iwrite_cr)
  else
    if (verbose > 0) &
      write(io%output,'(1p,g14.6,i6,3x,a,i7)') wallclock(), omp_get_thread_num(), 'I/O skip'
  end if
  call trace%end (itimer)
END SUBROUTINE check

END MODULE mpi_io_mod
