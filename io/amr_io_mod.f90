!===============================================================================
!> Module handling convenient AMR I/O.
!>
!> Each rank writes metadata to the file data/run/sssss/rank_rrrrr_tasks.nml
!> which contains namelists for each patch, including the offset into the
!> file run/SSSSS/rank_RRRRR_amr.dat, with the task data. These files are¨
!> consilidated into a single data/run/SSSSS/amr.nml in connection with
!> writing the data/run/SSSSS/amr.dat binary file.
!>
!> The metadata is sufficient to pre_init the tasks, and add them to a rank
!> tasks list.  In a loop through the task list the tasks are then allocated,
!> and the task offset is used to direct the rank to read in the task data.
!>
!> This can all be done "per rank", independent of other ranks, and a barrier
!> is only needed to ensure all ranks are done, before starting to execute.
!>
!> In case there are fewer ranks than in the saved snapshot (as given my the
!> data/run/parmeters.nml namelist file), each rank also reads in data from
!> the rank with number task%rank + i*mpi%size, until all saved rank data has
!> been read.
!>
!> Calling sequence:
!>
!> restart_amr_%init
!>   tasks_t%pre_init (reads task_nml, with type, rank, offset)
!>   amr_io%input
!>
!===============================================================================
MODULE amr_io_mod
  USE trace_mod
  USE io_unit_mod
  USE io_mod
  USE kinds_mod
  USE task_mod
  USE patch_mod
  USE mpi_file_mod
  USE mpi_mesg_mod
  USE mpi_mod
  USE dll_mod
  implicit none
  private
  !.............................................................................
  type, public:: amr_io_t
    integer:: iout                           ! snapshot number
    real(8):: time=0d0, out_next=0d0         ! time and output time
    integer(8):: size=0_8                    ! output size in bytes
    integer(8):: offset=0_8                  ! offset size in bytes
    type(mpi_counter_t):: file_offset        ! rank offset in variable
    type(mesg_list_t):: mesg_list
    type(mpi_file_t):: file                  ! output file
  contains
    procedure:: init
    procedure:: output
    procedure, private:: output_data
    procedure, private:: output_metadata
    procedure, private:: consolidate_metadata
    procedure:: input
    procedure:: open_file
    procedure:: close_file
  end type
  !.............................................................................
  integer:: verbose=0                        ! verbosity of log info
  type(amr_io_t), public:: amr_io            ! global instance
CONTAINS

!===============================================================================
!> Initialize AMR I/O, by initializing the rank counter; the counter resets to
!> the initial value after hitting zero.
!===============================================================================
SUBROUTINE init (self, task)
  class(amr_io_t):: self
  class(task_t):: task
  namelist /amr_io_params/ verbose
  integer:: iostat
  logical, save:: first_time = .true.
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('amr_io_t%init')
  call io%lock%set ('amr_io_t')
  if (first_time) then
    rewind (io%input)
    read (io%input, amr_io_params, iostat=iostat)
    call io%namelist_warning ('amr_io_params', iostat)
    write (io%namelist, amr_io_params)
    !...........................................................................
    call io%counter%init (mpi%size)
    call self%file_offset%init (0_8)
    first_time = .false.
  end if
  call io%lock%unset ('amr_io_t')
  call self%mesg_list%init
  self%size = 0
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  if (.not.associated(task%mesg)) then
    call task%allocate_mesg
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> This procedure is called each time a task is at its task%out_next time.
!> Each task has a task%pack() procedure, which is also used for communication.
!> In addition each task should have a task%output_metadata (unit) procedure,¨
!> which writes out whatever metadata that is needed to construct it.
!===============================================================================
SUBROUTINE output (self, task, count)
  class(amr_io_t):: self
  class(patch_t):: task
  integer:: count, remains
  type(mesg_t), pointer:: mesg
  character(120):: filename
  integer, save:: unit=0
  !-----------------------------------------------------------------------------
  ! Check if we still have unfinished output on this or other ranks
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('amr_io_t%output')
  !-----------------------------------------------------------------------------
  ! If the self time is at its output time we may have pending bussines
  !-----------------------------------------------------------------------------
  if (self%time >= self%out_next) then
    !-----------------------------------------------------------------------------
    ! Open a file for metadata output, if not already open
    !-----------------------------------------------------------------------------
    if (unit==0) then
      write (filename,'(a,i5.5,"/rank_",i5.5,"_amr.nml")') &
        trim(io%run_dir), task%iout, mpi%rank
      print *, trim(filename)
      open (newunit=unit,file=filename,form='formatted',status='unknown')
      call self%output_metadata (task, unit, snapshot=.true.)
      self%offset = 0
    end if
    !---------------------------------------------------------------------------
    ! If the task time is at its output time, write out meta-data, and save
    ! its output data in a message list, while reserving an offset slot
    !---------------------------------------------------------------------------
    if (task%time >= task%out_next) then
      call self%output_metadata (task, unit)
      call task%pack
      call self%mesg_list%add (task%mesg)
      task%mesg%offset = self%offset            ! FIXME: -> MPI
      self%offset = self%offset + task%mesg%nbuf*4_8
      self%size = self%size + task%mesg%nbuf
      task%out_next = task%out_next + io%out_time
    end if
    !---------------------------------------------------------------------------
    ! When all tasks on this rank are done, decrement the rank count
    !---------------------------------------------------------------------------
    if (count == 0) then
      remains = io%counter%update(-1)
      !-------------------------------------------------------------------------
      ! If all ranks are done do collective output, and update task%next_out
      !-------------------------------------------------------------------------
      if (remains <= 0) then
        close (unit)
        unit = 0
        call self%consolidate_metadata
        call self%output_data
        self%out_next = self%out_next + io%out_time
        !-----------------------------------------------------------------------
        if (mpi%master) then
          call self%file_offset%reset (0_8)
        end if
      end if
    end if
  end if
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!> This procedure is called when all ranks are ready to do I/O.  For now, we do
!> it in a blocking fashion
!===============================================================================
SUBROUTINE output_data (self)
  class(amr_io_t):: self
  type(mesg_t), pointer:: mesg, next
  real(8):: size
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('amr_io_t%output_metadata')
  !-----------------------------------------------------------------------------
  ! Update the file offset
  !-----------------------------------------------------------------------------
  self%offset = self%file_offset%update (self%size*4_8)
  !-----------------------------------------------------------------------------
  ! Write out the stored task buffers for each task (FIXME: buffer output)
  !-----------------------------------------------------------------------------
  call self%open_file (io%data_dir, self%iout)
  mesg => self%mesg_list%head
  do while (associated(mesg))
    call self%file%write (mesg%offset, mesg%nbuf, mesg%buffer)
    mesg => mesg%next
  end do
  self%offset = 0_8
  self%size = 0
  !-----------------------------------------------------------------------------
  ! Close the file, which ends the critical collective call section
  !-----------------------------------------------------------------------------
  call self%close_file
  call trace%end ()
END SUBROUTINE output_data

!===============================================================================
!> Consolidates all the rank metadata files into one file
!===============================================================================
SUBROUTINE consolidate_metadata (self)
  class(amr_io_t):: self
  !-----------------------------------------------------------------------------
  integer:: unit_in, unit_out, rank, iostat
  character(120):: line, filename
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%consolidate__metadata')
  if (mpi%master) then
    call io%lock%set ('amr_io_t%consolidate_metadata')
    !---------------------------------------------------------------------------
    call io%check_trailing_sep (io%run_dir)
    write (filename,'(a,i5.5,"/amr.nml")') trim(io%run_dir), self%iout
print *,filename
    open (newunit=unit_out, file=filename, form='formatted', status='unknown')
    if (verbose > 0) &
      write (stdlog,'(a,a)') 'amr_io_t%output: writing ', trim(filename)
    !---------------------------------------------------------------------------
    ! Loop over rank files
    !---------------------------------------------------------------------------
    do rank=0,mpi%size-1
      write (filename,'(a,i5.5,"/rank_",i5.5,"_amr.nml")') &
        trim(io%run_dir), self%iout, mpi%rank
print *,filename
      if (verbose > 1) &
        write (stdlog,'(a,a)') 'amr_io_t%output: reading ', trim(filename)
      open (newunit=unit_in, file=filename, form='formatted', status='old')
      do while (.true.)
        read (unit_in,'(a120)', iostat=iostat) line
        if (iostat /= 0) then
          exit
        end if
        write (unit_out,'(a)') trim(line)
      end do
      close (unit_in, status='delete')
    end do
    close (unit_out)
    call io%lock%unset ('amr_io_t%oconsolidate_metadatautput')
  end if
  call trace%end ()
END SUBROUTINE consolidate_metadata

!===============================================================================
SUBROUTINE open_file (self, dir, snap)
  class(amr_io_t):: self
  character(120):: dir, filename
  integer:: snap
  !-----------------------------------------------------------------------------
  call io%check_trailing_sep (dir)
  write (filename,'(a,i5.5,"/amr.dat")') trim(dir), snap
print *, filename  ! DEBUG
  call self%file%open (filename)
END SUBROUTINE open_file

!===============================================================================
SUBROUTINE close_file (self)
  class(amr_io_t):: self
  !-----------------------------------------------------------------------------
  call self%file%close
END SUBROUTINE close_file

!===============================================================================
!> AMR I/O input.  The patch is assumed to exist, with correct dimensions of
!> %mem() and patch%mesg%buffer, based on the reading of metadata.  All the is
!> required of the calling procedure is to call this procedure as many times
!> as there are patches to read.
!===============================================================================
SUBROUTINE input (self, task)
  class(amr_io_t):: self
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%input')
  !-----------------------------------------------------------------------------
  ! Read the stored task buffer for a task and unpack
  !-----------------------------------------------------------------------------
  call self%file%read (task%amr_offset, task%mesg%nbuf, task%mesg%buffer)
  call task%unpack (task%mesg)
  call trace%end ()
END SUBROUTINE input

!===============================================================================
!> write out metadata for a task to the tiven unit, with a short generic
!> amr_nml namelist, followed by a task-specific namelist with arbitrary
!> content.  This is called from amr_io_t%output(), once for each task.
!===============================================================================
SUBROUTINE output_metadata (self, task, unit, snapshot)
  class(amr_io_t):: self
  class(task_t):: task
  integer:: id, unit
  logical, optional:: snapshot
  !.............................................................................
  character(16):: type
  integer:: rank, nml_version
  integer(8):: offset
  namelist /amr_nml/ id, type, rank, offset, nml_version
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%output_metadata')
  call io%lock%set ('output_metadata')
  id = task%id
  type = task%type
  rank = mpi%rank
  offset = self%offset
  nml_version = io%nml_version
  write (unit, amr_nml)
  call task%output_metadata (unit, snapshot)
  call io%lock%unset ('output_metadata')
  call trace%end ()
END SUBROUTINE output_metadata

END MODULE amr_io_mod
