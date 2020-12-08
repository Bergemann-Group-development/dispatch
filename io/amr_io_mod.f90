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
    type(mpi_counter_t):: ranks_counter      ! number of ranks not ready
    type(mpi_counter_t):: file_offset        ! rank offset in variable
    type(dll_t):: list                       ! output list
    type(mpi_file_t):: file                  ! output file
  contains
    procedure:: init
    procedure:: update
    procedure:: output
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
  call trace%begin ('amr_io_t%init')
  call io%lock%set ('amr_io_t')
  if (first_time) then
    rewind (io%input)
    read (io%input, amr_io_params, iostat=iostat)
    call io%namelist_warning ('amr_io_params', iostat)
    write (io%namelist, amr_io_params)
    !...........................................................................
    call self%ranks_counter%init (mpi%size)
    call self%file_offset%init (0_8)
    call self%list%init ('amr_io')
    first_time = .false.
  end if
  call io%lock%unset ('amr_io_t')
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
  class(dll_node_t), pointer:: task_node
  !-----------------------------------------------------------------------------
  ! Check if we still have unfinished output on this or other ranks
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%update')
  !---------------------------------------------------------------------------
  ! If the self time is at its output time we have bussines
  !---------------------------------------------------------------------------
  if (self%time >= self%out_next) then
    !---------------------------------------------------------------------------
    ! If the task time is at its output time, write out meta-data, and save 
    ! its output data in the output list
    !---------------------------------------------------------------------------
    if (task%time >= task%out_next) then
      call task%output_metadata ()
      allocate (task_node)
      call task%pack_output (task_node%car)
      call self%list%append (task_node)
    end if
    !---------------------------------------------------------------------------
    ! When all tasks on this rank are done, decrement the rank count
    !---------------------------------------------------------------------------
    if (count == 0) then
      remains = self%ranks_counter%update(-1)
      !-------------------------------------------------------------------------
      ! If all ranks are done do collective output, and update task%next_out
      !-------------------------------------------------------------------------
      if (remains == 0) then
        call self%output
        self%out_next = self%out_next + io%out_time
        !-----------------------------------------------------------------------
        if (mpi%master) then
          call self%ranks_counter%reset (mpi%size)
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
  class(dll_node_t), pointer:: node, next
  real(8):: size
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%fill_buffer')
  !-----------------------------------------------------------------------------
  ! Compute the total size on the rank, and update global offset with it
  !-----------------------------------------------------------------------------
  node => self%list%head
  self%size = 0
  do while (associated(node))
    self%size = self%size + product(shape(node%car))*storage_size(node%car)
    node => node%next
  end do
  !-----------------------------------------------------------------------------
  ! Convert to bytes and update the file offset
  !-----------------------------------------------------------------------------
  self%size = self%size*4_8
  self%offset = self%file_offset%update (self%size)
  !-----------------------------------------------------------------------------
  ! Write out the stored task buffers for each task (FIXME: buffer output)
  !-----------------------------------------------------------------------------
  write (self%file%name,'(a,i5.5,"/amr.dat")') trim(io%rundir), self%iout
  call self%file%open
  node => self%list%head
  do while (associated(node))
    next => node%next
    call self%file%write (node%car)
    self%offset = self%offset + self%size
    deallocate (node%car)
    deallocate (node)
    node => next
  end do
  self%list%n = 0  
  self%offset = 0_8
  !-----------------------------------------------------------------------------
  ! Close the file, which ends the critical collective call section
  !-----------------------------------------------------------------------------
  call self%file%close
  call trace%end ()
END SUBROUTINE output

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