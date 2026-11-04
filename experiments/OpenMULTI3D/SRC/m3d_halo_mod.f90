!===============================================================================
!> Task links of OpenMULTI3D: neighbour lists, the time-ordered ready queue and the
!> copy of neighbour data into a task.
!===============================================================================
MODULE m3d_halo_mod
  USE os_mod
  USE list_mod
  USE halo_mod
  USE link_mod
  USE task_mod
  USE patch_mod, only:patch_t
  USE gpatch_mod, only:gpatch_t
  USE mesh_mod, only:mesh_t
  ! USE interp_mod
  USE io_mod
  USE timer_mod
  USE trace_mod
  USE mpi_mod
  USE shared_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE bits_mod
  USE console_mod
  USE handyfuncs_mod, only:option
  implicit none
  private

  !===============================================================================
  !> A halo_t link whose download copies whole neighbour patches (patch_t%copy_to).
  !===============================================================================
  type, public, extends(halo_t):: m3d_halo_t
  contains
    procedure:: info
    procedure:: append_to
    procedure:: reset_status
    procedure:: update
    procedure:: download
    procedure:: check_status
  end type

  logical, public, save :: print_finished=.false. !< report each task when it finishes
  logical, public, save :: print_lock=.false.     !< report every lock and unlock of a patch's memory

  public :: lock_mem, unlock_mem, add_nbor_by_key

CONTAINS

!===============================================================================
!> Add the task stored under key in the task list's hash table, if any, to the
!> nbor list nbors of task, and count it
!===============================================================================
SUBROUTINE add_nbor_by_key (task, key, nbors)
  class(gpatch_t):: task         !< task whose nbor list grows
  integer:: key(4)               !< hash key of the neighbour
  class(nbor_t), pointer:: nbors !< nbor list
  class(link_t), pointer:: nblink
  logical:: added
  !-----------------------------------------------------------------------------
  call task%task_list%hash_table%get (key, nblink)
  if (associated(nblink)) then
    select type (link => task%link)
    class is (halo_t)
      call link%conditionally_add_nbor (nblink, nbors, added)
      if (added) then
        !$omp atomic update
        task%n_nbors = task%n_nbors + 1
      end if
    end select
  end if
END SUBROUTINE add_nbor_by_key

  !===============================================================================
  !> Inserts the link into list, which is ordered by task time, after the links of
  !> equal time.
  !===============================================================================
  SUBROUTINE append_to(self, list)
    class (m3d_halo_t), target :: self
    class (list_t), pointer :: list !< list to insert into
    class (link_t), pointer :: next, prev
    class (task_t), pointer :: task
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('m3d_halo_t%append_to', itimer=itimer)
    nullify(self%next)
    nullify(self%prev)
    call                                      list%lock%set ('list_t%append_link')
    nullify (prev)
    next => list%head
    do while (associated(next))
      if (next%task%time > self%task%time) then
        self%next => next
        self%prev => prev
        next%prev => self
        exit
      endif
      prev => next
      next => next%next
    end do
    if (associated(prev)) then
      prev%next => self
      self%prev => prev
    else
      list%head => self
    endif
    if (.not. associated(next)) then
      list%tail => self
    endif
    list%n = list%n + 1
    call                                    list%lock%unset ('list_t%append_link')
    call trace%end(itimer)
  END SUBROUTINE append_to

  !===============================================================================
  !> Calls the task's pre_update and then halo_t's update.
  !===============================================================================
  SUBROUTINE update (self, all_tgt, all_src, only)
    class(m3d_halo_t)       :: self
    class(nbor_t), pointer  :: nbor
    logical, optional       :: all_tgt !< passed to halo_t's update
    logical, optional       :: all_src !< passed to halo_t's update
    integer, optional       :: only    !< passed to halo_t's update
    !-----------------------------------------------------------------------------
    call self%task%pre_update()
    !call self%halo_t%update(all_tgt, all_src, only)
    call halo_update(self, all_tgt, all_src, only)
  END SUBROUTINE update

  !===============================================================================
  !> Locks the patch's memory, for reading only with read.
  !===============================================================================
  SUBROUTINE lock_mem (self, label, read)
    class(patch_t) :: self
    character(len=*) :: label !< label of the lock, for diagnostics
    logical, optional :: read !< read lock (default .false.)
    !-----------------------------------------------------------------------------
    if (option(read, .false.)) then
      call self%mem_lock(1)%set_read (label)
    else
      call self%mem_lock(1)%set (label)
    endif
    if (print_lock) print '(A20, I6, A10)', trim(self%kind), self%id, 'locked'
  END SUBROUTINE lock_mem

  !===============================================================================
  !> Releases a lock_mem lock of the patch's memory.
  !===============================================================================
  SUBROUTINE unlock_mem (self, label, read)
    class(patch_t) :: self
    character(len=*) :: label !< label of the lock, for diagnostics
    logical, optional :: read !< read lock (default .false.)
    !-----------------------------------------------------------------------------
    if (option(read, .false.)) then
      call self%mem_lock(1)%unset_read (label)
    else
      call self%mem_lock(1)%unset (label)
    endif
    if (print_lock) print '(A20, I6, A10)', trim(self%kind), self%id, 'unlocked'
  END SUBROUTINE unlock_mem

  !===============================================================================
  !> Copies the data of every neighbour in nbors that is marked for download, under a
  !> read lock of its memory, into the task.
  !===============================================================================
  SUBROUTINE download (self, nbors, only, all_tgt, all_src)
    class(m3d_halo_t)      :: self
    class(nbor_t), pointer :: nbors   !< nbor list
    integer,      optional :: only    !< unused, part of the halo_t interface
    logical,      optional :: all_tgt !< unused, part of the halo_t interface
    logical,      optional :: all_src !< unused, part of the halo_t interface
    !.............................................................................
    class(nbor_t), pointer :: tnbors
    integer                :: unit, tid
    logical                :: all_tgt_L, all_src_L
    class(nbor_t), pointer :: nbor
    class(task_t), pointer :: target, source, patch
    integer, save          :: itimer=0
    !-----------------------------------------------------------------------------
    ! Decide whether to use a coarse clone or not
    !-----------------------------------------------------------------------------
    call trace%begin ('halo_t%download', itimer=itimer)
    all_tgt_L = shared%option (all_tgt)
    all_src_L = shared%option (all_src)
    unit = merge (stdout, stdlog, self%verbose > 2)
    target => self%task
    tid = target%id
    !-----------------------------------------------------------------------------
    ! If the task does not yet have support it will soon have it, and until then
    ! we make sure to fill all guard zones with values from the parent task, which
    ! has at least time-extrapolated guard zone values, and possibly already
    ! authoritative guard zone values
    !-----------------------------------------------------------------------------
    call self%set_filled (target)
    !-----------------------------------------------------------------------------
    ! Pass on, or loop over nbors
    !-----------------------------------------------------------------------------
    nbor => nbors
    do while (associated(nbor))
      if (nbor%download) then
        if (nbor%task%is_clear (bits%remove)) then
          source => nbor%task
          if (self%verbose > 2) then
            write (unit,'(a,2(i6,":",i2.2,3i3,2x))') &
              'halo_t%download: tid,sid=', &
              target%id, target%key, source%id, source%key
            flush (unit)
          end if
          select type (target)
           class is (patch_t)
            select type (source)
             class is (patch_t)
              call lock_mem(source, 'download', read=.true.)
              !call lock_mem(target, 'download')
              !target%mem = source%mem
              call source%copy_to(target)
              !call unlock_mem(target, 'download')
              call unlock_mem(source, 'download', read=.true.)
            end select
          end select
        end if
      end if
      nbor => nbor%next
    end do
    call self%get_filled (target)
    !.............................................................................
    if (self%verbose > 2) then
      write (unit,'(f12.6,a,3i5)') &
        wallclock(), pad_trace(1:indent)// &
        'halo_t%download: n_same,prolong,restrict=', self%n_op
      self%n_op = 0
    end if
    call trace%end (itimer)
  END SUBROUTINE download

  !===============================================================================
  !> Writes the task's info and one line per neighbour to unit.
  !===============================================================================
  SUBROUTINE info (self, label, unit)
    class(m3d_halo_t):: self
    character(*), optional:: label !< unused, part of the halo_t interface
    integer, optional:: unit       !< unit to write to (default io_unit%mpi)
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nt
    real(8):: distance(3)
    integer:: unit_L
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%info')
    unit_L = shared%option (unit, io_unit%mpi)
    call self%task%task_info(unit_L)
    nbor => self%nbor
    do while (associated(nbor))
      nt => nbor%task
      distance = nt%distance(self%task)
      write(unit_L,'(3x,a16,a,i6,":",i2.2,i4,2x,3l1,2x,2l1,3x,3f9.3)') &
        nt%kind, ' nbor_info: id, rank, level, needed, needs_me, download, BV, pos =', &
        nt%id, nt%level, nt%rank, nbor%needed, nbor%needs_me, &
        nbor%download, nt%is_set(bits%boundary), &
        nt%is_set(bits%virtual), distance/ &
        (0.5d0*(self%task%size+nt%size))
      nbor => nbor%next
    end do
    call trace%end ()
  END SUBROUTINE info

  !===============================================================================
  !> halo_t's check_status, which also reports a finished task with print_finished.
  !===============================================================================
  SUBROUTINE check_status (self, label)
    class(m3d_halo_t):: self
    character(*):: label !< passed to halo_t's check_status
    !.............................................................................
    class(task_t), pointer:: task
    !-----------------------------------------------------------------------------
    call self%halo_t%check_status(label)
    if (print_finished .and. self%task%finished) then
      task => self%task
      select type (task)
       class is (gpatch_t)
        if (console_on(ch_m3d, 2, all_ranks=.true.)) then
          write(console_unit,*) task%id, ' has finished!', task%task_list%nfinished, '/', task%task_list%na, timer%done
        endif
      end select
    endif
  END SUBROUTINE check_status

  !===============================================================================
  !> Reset the status bits and flags on a linked task and an nbor list
  !===============================================================================
  SUBROUTINE reset_status (self, nbor_list, check)
    class(m3d_halo_t):: self
    class(link_t), pointer:: nblink
    class(nbor_t), pointer:: nbor_list !< nbor list of the task
    logical, optional:: check          !< unused, part of the halo_t interface
    !.............................................................................
    class(nbor_t), pointer:: nb
    class(task_t), pointer:: task
    type(task_t):: new
    !-----------------------------------------------------------------------------
    call trace%begin ('link_t%reset_status')
    !-----------------------------------------------------------------------------
    ! Set new status, as a modified copy of self%status
    !-----------------------------------------------------------------------------
    task => self%task
    new%status = task%status
    call new%clear (bits%boundary+bits%virtual)
    if (task%rank == mpi%rank) then
      call new%set (bits%internal)
    else
      call new%set (bits%external)
    end if
    !-----------------------------------------------------------------------------
    ! Loop over nbor list and modify bits; B and V tasks have only that bit
    !-----------------------------------------------------------------------------
    nb => nbor_list
    do while (associated(nb))
      if (self%verbose > 1) then
        write (stdlog,'(a,2(i7,i5,2x))') &
          'reset_status: nbid,rank=', task%id, task%rank, nb%task%id, nb%task%rank
      end if
      if (nb%needs_me .and. nb%task%is_clear(bits%no_download)) then
        if (task%rank == mpi%rank) then
          if (nb%task%rank /= mpi%rank) then
            call new%set (bits%boundary)
            call new%clear (bits%internal+bits%external)
          end if
        else
          if (nb%task%rank == mpi%rank) then
            call new%set (bits%virtual)
            call new%clear (bits%internal+bits%external)
          end if
        end if
      endif
      nb => nb%next
    end do
    !-----------------------------------------------------------------------------
    ! If any change, call check_status() to take proper action
    !-----------------------------------------------------------------------------
    if (new%status /= task%status) then
      !$omp atomic write
      task%status = new%status
      task%previous_status = task%status
    end if
    !-----------------------------------------------------------------------------
    !if (shared%option(check)) then
    !  call self%check_nbor_relations
    !end if
    !.............................................................................
    if (self%verbose > 0) then
      write(stdlog,'(a,i6,i4,3x,4l1)') &
        'reset_status: id,rank,IBVE=', self%task%id, self%task%rank, &
        self%task%is_set(bits%internal), &
        self%task%is_set(bits%boundary), &
        self%task%is_set(bits%virtual ), &
        self%task%is_set(bits%external)
    end if
    call trace%end ()
  END SUBROUTINE reset_status


END MODULE
