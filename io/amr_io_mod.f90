!===============================================================================
!> Module handling convenient AMR I/O.  We wish to have a file format that makes
!> it easy to read when restarting jobs, and which makes it easy to read with
!> Python.  We can rely on the existing mechanism to write patch metadata in the
!> form of file run/snapno/rank_rankno_patches.nml with namelist data for each
!> patch.  This may be replaced later with a file written with mpi_buffer_mod.
!>
!> We can reuse the functionality in mpi_file_mod, and build a task list with
!> patch data, ready for output.  We also reuse the mechanism which counts down 
!> the number of tasks still remaining to be added to the list -- this saves
!> debugging efforts, since it is known to work and is rather delicate, since
!> new AMR patches may be added after the I/O process has started, until it is
!> complete, with all existing tasks having passed the current out_next value.
!>
!> Note that the output procedure is called from inside a critical region in
!> data_io_mod.f90, and hence everything in output and the procedures it calls
!> is trivially thread safe; only the last thread to reach next_out enters here,
!> and while it is active other threads are waiting at the start of the critical
!> region.  When they enter, they find that the output has been completed.
!>
!> Call hierarchy:
!>
!> experiment_t%output                  ! request patch output
!>   gpatch_t%output                    ! catch output request
!>     data_io_t%outpout                ! select output method
!>       amr_io_t%output                ! check if final thread & rank
!>         amr_io_t%output_list         ! run through list
!>           mpi_file_t%open            ! open data/run/IOUT/snapshot.dat
!>           amr_io%t%output_buffer     ! fill buffer for variable iv
!>             mpi_file_t%write         ! compute offset and size
!>               MPI_File_write_at      ! write out
!>           mpi_file_t%close           ! close file
!===============================================================================
MODULE amr_io_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE list_mod
  USE mpi_file_mod
  USE mpi_mod
  implicit none
  private
  !.............................................................................
  type, public:: amr_io_t
    integer:: n(3)
    integer(8):: var_size                    ! words per variable in the file
    type(mpi_counter_t):: ranks
    type(list_t), pointer:: list => null()
    real, allocatable:: buffer(:,:,:,:,:)
    type(mpi_file_t):: file
  contains
    procedure:: init
    procedure:: output
    procedure:: check
    procedure:: output_list
    procedure:: output_buffer
    procedure:: input
  end type
  !.............................................................................
  integer:: verbose=0
  type(amr_io_t), public:: amr_io
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
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, amr_io_params, iostat=iostat)
    write (io%output, amr_io_params)
  end if
  call self%ranks%init (mpi%size)
  if (verbose > 0) &
    write (stdout,"('amr_io_t%init: ',a,i6)") &
     'ranks%i =', self%ranks%i
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> AMR I/O output. For the most exact restart results, the two time slices
!> bracketing the out_next time should both be saved -- interpolation in time
!> would not conserve mass, momentum, energu, and magnetic flux divergence.
!===============================================================================
SUBROUTINE output (self, patch, count)
  class(amr_io_t):: self
  class(patch_t):: patch
  class(task_t), pointer:: task
  class(patch_t), pointer:: out
  integer:: count, l(3), u(3), prv, now, remains
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%output')
  !-----------------------------------------------------------------------------
  ! Allocate a task with two time slots, for the two bracketing time slices,
  ! and with a size that optionally include the guard zones.
  !-----------------------------------------------------------------------------
  allocate (out, source=patch)
  out%nw = 1
  out%nt = 2
  if (io%guard_zones) then
    l = patch%mesh%lb
    u = patch%mesh%ub
  else
    l = patch%mesh%li
    u = patch%mesh%ui
  end if
  self%n(1:3) = u-l+1
  allocate (out%mem(self%n(1),self%n(2),self%n(3),out%nv,out%nt,out%nw))
  now = patch%it
  if (patch%istep > 1) then
    prv = 1 + modulo(now-2,patch%nt)
  else
    prv = now
  end if
  !-----------------------------------------------------------------------------
  ! Copy the two bracketing time slots to the output task
  !-----------------------------------------------------------------------------
  out%mem(:,:,:,:,1,1) = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,prv,1)
  out%mem(:,:,:,:,2,1) = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,now,1)
  allocate (out%t(2), out%dt(2))
  out%t(1)  = patch%t(prv)
  out%t(2)  = patch%t(now)
  out%dt(1) = patch%dt(prv)
  out%dt(2) = patch%dt(now)
  !-----------------------------------------------------------------------------
  ! Append the task copy to a list_t instance
  !-----------------------------------------------------------------------------
  if (.not.associated(self%list)) then
    allocate (self%list)
    call self%list%init
  end if
  task => out
  call self%list%append_task (task)
  !-----------------------------------------------------------------------------
  ! When the last task has been appended, decrement the counter that counts
  ! the number of ranks that are not yet ready
  !-----------------------------------------------------------------------------
  if (verbose > 1) &
    write (stdout,"('amr_io_t%output: ',a,i6)") &
     'count =', count
  if (count == 0) then
    remains = self%ranks%update(-1) - 1
    if (verbose > 0) &
      write (stdout,"('amr_io_t%output: ',a,2i6)") &
       'remains =', remains, self%ranks%i
    !---------------------------------------------------------------------------
    ! If other ranks remain, set a flag forcing data_io_t%output() to make
    ! repeated calls to self%check
    !---------------------------------------------------------------------------
    if (remains > 0) then
      io%needs_check = .true.
    else
      call self%check
    end if
  end if
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!> Check if all ranks are ready to do I/O
!===============================================================================
SUBROUTINE check (self)
  class(amr_io_t):: self
  integer:: remains
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%check')
  remains = self%ranks%update(0)
  if (remains == 0) then
    call self%output_list
    call self%list%reset
    io%needs_check = .false.
    if (mpi%master) then
      remains = self%ranks%update (mpi%size)
    end if
  end if
  call trace%end ()
END SUBROUTINE check

!===============================================================================
!> Actual AMR I/O output.  This procedure is called when all ranks are ready to
!> do I/O.  For now, we do it in a blocking fashion
!===============================================================================
SUBROUTINE output_list (self)
  class(amr_io_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  integer:: i, iv, n(5)
  character(len=64):: filename
  !-----------------------------------------------------------------------------
  ! For each variable, allocate output buffer, fill it, and write it out
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%fill_buffer')
  n(1:3) = self%n
  n(4) = 2
  n(5) = self%list%n
  allocate (self%buffer(n(1),n(2),n(3),n(4),n(5)))
  !-----------------------------------------------------------------------------
  ! Open the output file and loop over variables
  !-----------------------------------------------------------------------------
  link => self%list%head
  task => link%task
  write (filename,'(a,i5.5,"/snapshot.dat")') &
    trim(io%outputname), task%iout
  if (verbose > 0) &
    write (stdout,"('amr_io_t%fill_buffer: ',a)") &
      'filename = '//trim(filename)
  call self%file%openw (filename)
  do iv=1,n(4)
    i = 1
    link => self%list%head
    do while (associated(link))
      task => link%task
      select type (task)
      class is (patch_t)
        self%buffer(:,:,:,:,i) = task%mem(:,:,:,iv,:,1)
        i = i + task%nt
      end select
      link => link%next
    end do
    call self%output_buffer (iv)
  end do
  deallocate (self%buffer)
  call self%file%close
  call trace%end ()
END SUBROUTINE output_list

!===============================================================================
!> Actual AMR I/O output.  Compute offset into file, open the file, output the
!> buffer, and close the file.   The file operations are collective, and hence
!> this needs to be done by just one thread, while all other threads either wait
!> (conservative choice), or continue.
!===============================================================================
SUBROUTINE output_buffer (self, iv)
  class(amr_io_t):: self
  integer:: iv, rank_size, task_size
  integer(8):: offset
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%output_buffer')
  task_size = product(shape(self%buffer))
  rank_size = task_size*self%list%n
  self%var_size = rank_size*mpi%size
  offset = iv*self%var_size + mpi%rank*rank_size
  if (verbose > 0) &
    write (stdout,"('amr_io_t%output_buffer: ',a,5i12)") &
     'task_size, rank_size, self%var_size, iv, offset =', &
      task_size, rank_size, self%var_size, iv, offset
  !call self%file%write (4_8*offset, rank_size, self%buffer)
  call trace%end ()
END SUBROUTINE output_buffer 

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
