MODULE m3d_halo_mod
    USE os_mod
    USE list_mod
    USE halo_mod
    USE link_mod
    USE task_mod
    USE gpatch_mod, only:gpatch_t
    USE io_mod
    USE timer_mod
    USE trace_mod
    USE mpi_mod
    USE shared_mod
    USE omp_mod
    USE omp_lock_mod
    USE omp_timer_mod
    USE bits_mod
    USE guards_mod
    implicit none
    private
  
    type, public, extends(halo_t):: m3d_halo_t
    contains
      procedure:: info
      procedure:: append_to
      procedure:: reset_status
      procedure:: update
      procedure:: check_status
    end type

    logical, public, save :: print_finished=.false.

CONTAINS

  !===============================================================================
  SUBROUTINE append_to(self, list)
    class (m3d_halo_t), target :: self
    class (list_t), pointer :: list
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
  !> Reset the status bits and flags on a linked task and an nbor list
  !===============================================================================
  SUBROUTINE reset_status (self, nbor_list, check)
    class(m3d_halo_t):: self
    class(link_t), pointer:: nblink
    class(nbor_t), pointer:: nbor_list
    logical, optional:: check
    !.............................................................................
    class(nbor_t), pointer:: nb
    class(task_t), pointer:: task
    type(task_t):: new
    logical:: flags(3,2)
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
      if (nb%needs_me) then
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

  !===============================================================================
  SUBROUTINE update (self, all_tgt, all_src, only)
    class(m3d_halo_t)           :: self
    class(nbor_t), pointer  :: nbor
    logical, optional       :: all_tgt, all_src
    integer, optional       :: only
    !-----------------------------------------------------------------------------
    call self%task%pre_update()
    call self%halo_t%update(all_tgt, all_src, only)
  END SUBROUTINE update

  !===============================================================================
  SUBROUTINE info (self, label, unit)
    class(m3d_halo_t):: self
    character(*), optional:: label
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nt
    real(8):: distance(3)
    integer, optional:: unit
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
  SUBROUTINE check_status (self, label)
    class(m3d_halo_t):: self
    character(*):: label
    !.............................................................................
    class(task_t), pointer:: task
    type(task_t):: previous
    !-----------------------------------------------------------------------------
    call self%halo_t%check_status(label)
    if (print_finished .and. self%task%finished) then
      task => self%task
      select type (task)
      class is (gpatch_t)
        if (task%verbose > 4) then
          write(stdout,*) task%id, ' has finished!', task%task_list%nfinished, '/', task%task_list%na, timer%done
        endif
      end select
    endif
  END SUBROUTINE check_status

END MODULE