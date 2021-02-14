!===============================================================================
!> Module handling convenient AMR I/O.  We wish to have a file format that makes
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
    procedure:: update
    procedure:: output
    procedure:: output_metadata
    procedure:: input
  end type
  !.............................................................................
  integer:: verbose=0                        ! verbosity of log info
  type(amr_io_t), public:: amr_io            ! global instance
CONTAINS

!===============================================================================
!> Initialize AMR I/O, by initializing the rank counter; the counter resets to
!> the initial value after hitting zero.
!===============================================================================
SUBROUTINE init (self, patch)
  class(amr_io_t):: self
  class(patch_t):: patch
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
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> This procedure is called each time a task is at its task%out_next time.
!> Each task has a task%output_pack() procedure, which sets an anonymous
!> pointer to a buffer array of arbitrary shape -- for patches it is simply
!> a copy of the mem array time slot.
!===============================================================================
SUBROUTINE update (self, task, count)
  class(amr_io_t):: self
  class(patch_t):: task
  integer:: count, remains
  type(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  ! Check if we still have unfinished output on this or other ranks
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('amr_io_t%update')
  !-----------------------------------------------------------------------------
  ! If the self time is at its output time we have bussines
  !-----------------------------------------------------------------------------
  if (self%time >= self%out_next) then
    !---------------------------------------------------------------------------
    ! If the task time is at its output time, write out meta-data, and save 
    ! its output data in the output list
    !---------------------------------------------------------------------------
    if (task%time >= task%out_next) then
      call task%pack
      call self%mesg_list%add (task%mesg)
      self%size = self%size + task%mesg%nbuf
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
        call self%output
        self%out_next = self%out_next + io%out_time
        !-----------------------------------------------------------------------
        if (mpi%master) then
          call self%file_offset%reset (0_8)
        end if
      end if
    end if
  end if
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> This procedure is called when all ranks are ready to do I/O.  For now, we do
!> it in a blocking fashion
!===============================================================================
SUBROUTINE output (self)
  class(amr_io_t):: self
  type(mesg_t), pointer:: mesg, next
  real(8):: size
  !-----------------------------------------------------------------------------
  if (.not.io%do_output) return
  call trace%begin ('amr_io_t%fill_buffer')
  call self%output_metadata
  !-----------------------------------------------------------------------------
  ! Update the file offset
  !-----------------------------------------------------------------------------
  self%size = self%size*4_8
  self%offset = self%file_offset%update (self%size)
  !-----------------------------------------------------------------------------
  ! Write out the stored task buffers for each task (FIXME: buffer output)
  !-----------------------------------------------------------------------------
  write (self%file%name,'(a,i5.5,"/amr.dat")') trim(io%rundir), self%iout
  call self%file%open
  mesg => self%mesg_list%head
  do while (associated(mesg))
    next => mesg%next
    call self%file%write (self%offset, mesg%nbuf, mesg%buffer)
    call mesg%dealloc
    deallocate (mesg)
    mesg => next
  end do
  self%offset = 0_8
  self%size = 0
  !-----------------------------------------------------------------------------
  ! Close the file, which ends the critical collective call section
  !-----------------------------------------------------------------------------
  call self%file%close
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!> Write out metadata; a file with one namelist for each task, with the task id
!> and the task buffer size (in integer(4) units)
!===============================================================================
SUBROUTINE output_metadata (self)
  class(amr_io_t):: self
  integer:: id, size
  character(len=80):: file
  type(mesg_t), pointer:: mesg
  namelist /amr_nml/ id, size
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%output_metadata')
  write (file,'(a,i5.5,"/rank_",i5.5,"_amr.mml")') &
    trim(io%rundir), self%iout, mpi%rank
  open (unit=io_unit%data, file=file, form='formatted', status='unknown')
  self%size = 0
  mesg => self%mesg_list%head
  do while (associated(mesg))
    id = mesg%id
    size = mesg%nbuf
    self%size = self%size + size
    write (io_unit%data, amr_nml)
    mesg => mesg%next
  end do
  self%offset = 0_8
  !-----------------------------------------------------------------------------
  ! Close the file, which ends the critical collective call section
  !-----------------------------------------------------------------------------
  call self%file%close
  call trace%end ()
END SUBROUTINE output_metadata

!===============================================================================
!> AMR I/O input.  The patch is assumed to exist, with correct dimensions of
!> %mem(), based on the reading of metadata.
!===============================================================================
SUBROUTINE input (self, patch)
  class(amr_io_t):: self
  class(patch_t):: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%input')
  call trace%end ()
END SUBROUTINE input

END MODULE amr_io_mod
