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
!>     data_io_t%output                 ! select output method
!>       amr_io_t%output                ! check if final thread & rank
!>         amr_io_t%output_list         ! run through list
!>           mpi_file_t%open            ! open data/run/IOUT/snapshot.dat
!>           amr_io%t%output_buffer     ! fill buffer for variable iv
!>             mpi_file_t%write         ! compute offset and size
!>               MPI_File_write_at      ! write out
!>           mpi_file_t%close           ! close file
!>
!> OpenMPI aspects:  All threads call experiment_t%output repeatedly, but most
!> of the time they return from data_io%output, after finding that the time has
!> not advanced enough.  When reaching a time where they should add their data
!> to the output list, they must operate one-at-a-time, so that part should be
!> inside a unique critical region (or should use a lock).
!>
!> data_io%output indeed has name critical regions around all calls that end up
!> inside amr_io_mod, so only one thread can ever be inside the procedures below.
!> But we must also make sure that when that thread decides to do collective MPI
!> calls, which may only be done by one thread on some systems (laptops), all
!> other threads must be held up at the very same critical region (data_io_cr).
!> This is arranged by using an OMP lock in data_io_mod -- cf. comments there.
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
  USE time_slices_mod
  implicit none
  private
  !.............................................................................
  type, public:: amr_io_t
    integer:: iout                           ! snapshot number
    integer:: n(3)                           ! 3-D array dims in output
    integer:: nv                             ! number of variables in output
    integer:: task_size                      ! words per variable in a task
    integer:: var_size                       ! words per variable in the run
    integer:: rank_size                      ! words per variable in the rank
    integer(8):: rank_offset                 ! offset from variable start
    integer(8):: task_offset=0_8             ! offset in buffer
    type(mpi_counter_t):: ranks_counter      ! number of ranks not ready
    type(mpi_counter_t):: offset_counter     ! rank offset in variable
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
  integer:: verbose=0                        ! verbosity of log info
  integer:: nt=2                             ! number of time slots in output
  type(amr_io_t), public:: amr_io            ! global instance
CONTAINS

!===============================================================================
!> Initialize AMR I/O, by initializing the rank counter; the counter resets to
!> the initial value after hitting zero.
!===============================================================================
SUBROUTINE init (self, patch)
  class(amr_io_t):: self
  class(patch_t):: patch
  namelist /amr_io_params/ verbose, nt
  integer:: iostat
  logical, save:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, amr_io_params, iostat=iostat)
    write (io%output, amr_io_params)
    call self%ranks_counter%init (mpi%size)
    call self%offset_counter%init (0_8)
    if (nt == 2) then
      time_slices%order = 2
    end if
  end if
  if (verbose > 0) &
    write (stdout,"('amr_io_t%init: ',a,i6)") &
     'ranks%i =', self%ranks_counter%i
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
  integer:: count, l(3), u(3), n(3), prv, now, remains, it, iv
  real, allocatable:: tmp(:,:,:,:)
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%output')
  !-----------------------------------------------------------------------------
  ! Allocate a task with two time slots, for the two bracketing time slices,
  ! and with a size that optionally include the guard zones.
  !-----------------------------------------------------------------------------
  allocate (out)
  out%id = patch%id
  out%nt = nt
  out%nv = patch%nv
  out%nw = 1
  if (io%guard_zones) then
    l = patch%mesh%lb
    u = patch%mesh%ub
  else
    l = patch%mesh%li
    u = patch%mesh%ui
  end if
  n(1:3) = u-l+1
  allocate (out%mem(n(1),n(2),n(3),out%nv,out%nt,out%nw))
  if (first_time) then
  end if
  now = patch%it
  if (patch%istep > 1) then
    prv = 1 + modulo(now-2,patch%nt)
  else
    prv = now
  end if
  allocate (out%t(2), out%dt(2))
  if (nt == 2) then
    io%nml_version = 3
    !---------------------------------------------------------------------------
    ! Copy the two bracketing time slots to the output task
    !---------------------------------------------------------------------------
    do iv=1,patch%nv
      call convert_variables (patch, out%mem(:,:,:,iv,1,1), iv, prv)
      call convert_variables (patch, out%mem(:,:,:,iv,2,1), iv, now)
    end do
    out%t(1)  = patch%t(prv)
    out%t(2)  = patch%t(now)
    out%dt(1) = patch%dt(prv)
    out%dt(2) = patch%dt(now)
  else
    !---------------------------------------------------------------------------
    ! Alternatively, use time-slice interpolation to a single output time
    !---------------------------------------------------------------------------
    allocate (tmp(n(1),n(2),n(3),patch%nt-1))
    do iv=1,patch%nv
      do it=1,patch%nt-1
        call convert_variables (patch, tmp(:,:,:,it), iv, it)
      end do
      call time_slices%interpolate (patch, tmp, out%mem(:,:,:,iv,1,1))
    end do
    deallocate (tmp)
  end if
  self%task_size = product(n)*nt
  self%task_offset = self%task_offset + self%task_size
  patch%amr_offset = self%task_offset
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
    !---------------------------------------------------------------------------
    ! Increase the offset into the variable block by the size of the output
    ! buffer (computed when allocated) in this rank
    !---------------------------------------------------------------------------
    self%n = n
    self%nv = patch%nv
    remains = self%ranks_counter%update(-1) - 1
    !---------------------------------------------------------------------------
    ! If other ranks remain, set a flag forcing data_io_t%output() to make
    ! repeated calls to self%check
    !---------------------------------------------------------------------------
    if (verbose > 1) &
      write (stdout,"('amr_io_t%output: ',a,i6)") &
       'remains =', remains
    self%iout = patch%iout
    self%rank_size = self%task_size*self%list%n
    if (remains == 0) then
      call self%check
    else
      io%needs_check = .true.
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
  remains = self%ranks_counter%update(0)
  if (verbose > 1) &
    write (stdout,"('amr_io_t%check: ',a,i6)") &
     'remains =', remains
  if (remains == 0) then
    !---------------------------------------------------------------------------
    ! Now that all ranks are ready, the total size per variable may be computed,
    ! and the contents of the task list may be written out
    !---------------------------------------------------------------------------
    if (verbose > 0) &
      write (stdout,"('amr_io_t%check: ',a,i10)") &
        'var_size =', self%var_size
    call self%output_list
    !---------------------------------------------------------------------------
    ! As output is now complete, the output task list should be reset, and the
    ! MPI master should reset the rank count
    !---------------------------------------------------------------------------
    call self%list%reset
    io%needs_check = .false.
    self%task_offset = 0_8
    if (mpi%master) then
      call trace%begin ('counter_t%begin')
      call self%ranks_counter%reset (mpi%size)
      call self%offset_counter%reset (0_8)
      call trace%end ()
    end if
  end if
  call trace%end ()
END SUBROUTINE check

!===============================================================================
!> This procedure is called when all ranks are ready to do I/O.  For now, we do
!> it in a blocking fashion
!===============================================================================
SUBROUTINE output_list (self)
  class(amr_io_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  integer:: i, iv, n(5)
  integer(8):: rank_offset, var_size
  character(len=64):: filename
  namelist /offset_nml/ rank_offset, var_size
  !-----------------------------------------------------------------------------
  ! For each variable, allocate output buffer, fill it, and write it out
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%fill_buffer')
  n(1:3) = self%n
  n(4) = nt
  n(5) = self%list%n
  allocate (self%buffer(n(1),n(2),n(3),n(4),n(5)))
  !-----------------------------------------------------------------------------
  ! Compute offset and size
  !-----------------------------------------------------------------------------
  self%rank_offset = self%offset_counter%update (self%rank_size)
  self%var_size = self%offset_counter%update (0_8)
  io%gb_out = self%var_size*self%nv*4.0/1024.**3
   if (verbose > 0) &
    write (stdout,"('amr_io_t%output: ',a,5i10)") &
     'task_size, rank_size, var_size, rank_offset =', &
      self%task_size, self%rank_size, self%var_size, self%rank_offset
  !-----------------------------------------------------------------------------
  ! Open the output file and loop over variables
  !-----------------------------------------------------------------------------
  link => self%list%head
  task => link%task
  write (filename,'(a,i5.5,"/snapshot.dat")') &
    trim(io%outputname), self%iout
  if (verbose > 0) &
    write (stdout,"('amr_io_t%fill_buffer: ',a)") &
      'filename = '//trim(filename)
  !$omp atomic write
  io%halt = .true.
  !$omp end atomic
  call self%file%openw (filename)
  do iv=1,self%nv
    i = 1
    link => self%list%head
    do while (associated(link))
      task => link%task
      select type (task)
      class is (patch_t)
!print *, shape(self%buffer), 'buffer'
!print *, shape(task%mem), 'mem'
        self%buffer(:,:,:,:,i) = task%mem(:,:,:,iv,:,1)
        i = i + 1
      end select
      link => link%next
    end do
    call self%output_buffer (iv)
  end do
  deallocate (self%buffer)
  !-----------------------------------------------------------------------------
  ! Close the file, which ends the critical collective call section, so we can
  ! now remove the io%halt trap, which prevents other threads from getting into
  ! trouble in the mean time.
  !-----------------------------------------------------------------------------
  call self%file%close
  !$omp atomic write
  io%halt = .false.
  !$omp end atomic
  !-----------------------------------------------------------------------------
  ! Open the file "data/run/rank_rrrrr_patches.nml", appending offset
  !-----------------------------------------------------------------------------
  rank_offset = self%rank_offset
  var_size = self%var_size
  flush (io_unit%nml2)
  write (io_unit%nml2, offset_nml)
  flush (io_unit%nml2)
  call trace%end ()
END SUBROUTINE output_list

!===============================================================================
!> Actual AMR I/O output.  Compute offset into file, open the file, output the
!> buffer, and close the file.   The file operations are collective, and hence
!> this needs to be done by all ranks together, and by just one thread in each, 
!> while all other threads either wait (conservative choice), or continue.
!===============================================================================
SUBROUTINE output_buffer (self, iv)
  class(amr_io_t):: self
  integer:: iv
  integer(8):: offset
  !-----------------------------------------------------------------------------
  call trace%begin ('amr_io_t%output_buffer')
  !-----------------------------------------------------------------------------
  ! The offset in the snapshot file to this buffer is the offset over total
  ! variable size for this snapshot, plus the rank offset in a variable
  !-----------------------------------------------------------------------------
  offset = (iv-1)*self%var_size + self%rank_offset
  if (verbose > 0) &
    write (stdout,"('amr_io_t%output_buffer: ',a,2i12)") &
     'iv, offset =', iv, offset
  call self%file%write (4_8*offset, self%rank_size, self%buffer)
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

!===============================================================================
! Copy variables from `patch%mem` to `mem` array, converting variables based
! on the I/O format parameter.  In the process, the time slot ordering is
! rearranged into increasing time order.
!===============================================================================
SUBROUTINE convert_variables (patch, mem, iv, it)
  class(patch_t):: patch
  real:: mem(:,:,:)
  integer:: iv, it, jt
  integer:: io_format, l(3), u(3), i, j, k
  !...........................................................................
  call trace%begin ('amr_io-t%convert_variables')
  if (io%guard_zones) then
    l = patch%mesh%lb
    u = max(patch%mesh%ub,patch%mesh%gn)
  else
    l = patch%mesh%li
    u = patch%mesh%ui
  end if
  io_format = io%format
  if (io_format >= 10) then
    io_format = io_format-10
  else if (io_format >= 6) then
    io_format = io_format-6
  end if
  jt = patch%iit(it)
  if (io_format >= 2 .and. iv==patch%idx%d) then
    ! --- Convert density to log density ---
    mem(:,:,:) = log( &
        patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv         ,jt,1))
  else if ((io%format == 0 .or. io%format==2) .and. patch%solver_is('ramses') &
      .and. iv >= patch%idx%px .and. iv <= patch%idx%pz) then
    ! --- Convert momentum to velocity ---
    mem(:,:,:) &
      = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv         ,jt,1) &
      / patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,jt,1)
  else
    mem(:,:,:) = &
        patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv         ,jt,1)
  end if
  call trace%end()
END SUBROUTINE convert_variables

END MODULE amr_io_mod
