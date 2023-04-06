MODULE m3d_halo_mod
    USE os_mod
    USE list_mod
    USE halo_mod
    USE link_mod
    USE task_mod
    USE patch_mod, only:patch_t
    USE gpatch_mod, only:gpatch_t
    USE mesh_mod, only:mesh_t
    USE interp_mod
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
      procedure:: download
      procedure:: source_to_target
      procedure:: same
      procedure:: no_t
      procedure:: interp_to_higher_level
      procedure:: binning_to_lower_level
      procedure:: check_status
    end type

    logical, public, save :: print_finished=.false., print_download=.false.

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
  SUBROUTINE update (self, all_tgt, all_src, only)
    class(m3d_halo_t)       :: self
    class(nbor_t), pointer  :: nbor
    logical, optional       :: all_tgt, all_src
    integer, optional       :: only
    !-----------------------------------------------------------------------------
    call self%task%pre_update()
    !call self%halo_t%update(all_tgt, all_src, only)
    call halo_update(self, all_tgt, all_src, only)
  END SUBROUTINE update

  !===============================================================================
  !> Comprehensive entry point
  !===============================================================================
  SUBROUTINE download (self, nbors, only, all_tgt, all_src)
    class(m3d_halo_t)      :: self
    class(nbor_t), pointer :: nbors
    integer,      optional :: only
    logical,      optional :: all_tgt, all_src
    !.............................................................................
    class(nbor_t), pointer :: tnbors
    integer                :: unit, n_same, tid, i
    logical                :: use_coarse, all_tgt_L, all_src_L
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
          !call guards%source_to_target (source, target, n_op=self%n_op, &
          !  only=only, all_src=all_src, all_tgt=all_tgt, &
          !  nbstep=nbor%check_step)
          call self%source_to_target (source, target, only=only)
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
  SUBROUTINE source_to_target(self, source, target, only)
    class(m3d_halo_t)      :: self
    class(task_t)          :: source, target
    integer, optional      :: only
    integer, save          :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_halo_t%source_to_target', itimer=itimer)
    !.............................................................................
    select type (target)
    class is (patch_t)
      select type (source)
      class is (patch_t)
        call self%same(target, source, only=only)
      end select
    end select
    call trace%end (itimer=itimer)
  END SUBROUTINE source_to_target

  !===============================================================================
  SUBROUTINE same (self, tgt, src, only)
    class(m3d_halo_t)      :: self
    class(patch_t)         :: tgt, src
    integer, optional      :: only
    !.............................................................................
    class(mesh_t), pointer, dimension(:) :: tm, sm
    real(8), dimension(3)                :: sm_p, tl, tu, sl, su, fl, fu, pl, pu,d
    integer, dimension(3)                :: ls, us, lt, ut
    real(8)                              :: off
    integer, save          :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_halo_t%same', itimer=itimer)

    if (print_download) then
      write(stdout, '(I4, A, I4)') tgt%id, ' downloading from ', src%id
    endif

    tm => tgt%mesh                             ! tgt mesh
    sm => src%mesh                             ! src mesh
    !off = merge(0.5, -0.5, src%level > tgt%level)
    d  = tgt%distance_to (src)                 ! wrapped vector (src-tgt)
    tl =   - tm%s*0.5d0 - tm%ng*tm%d           ! tgt lower corner
    tu =   + tm%s*0.5d0 + tm%ng*tm%d           ! tgt upper corner
    sl = d - sm%s*0.5d0                        ! src lower corner
    su = d + sm%s*0.5d0                        ! src upper corner
    fl = max(tl,sl)                            ! overlap lower corner
    fu = min(tu,su)                            ! overlap upper corner
    call io%assert (all(fu > fl), 'guards_t%differ: corner reversal')

    lt = ceiling((fl  )/tm%d + tm%o)      ! tgt lowwe index
    ut = floor  ((fu  )/tm%d + tm%o)      ! tgt upper index

    if (tgt%level == src%level) then
      ls = ceiling((fl-d)/sm%d + sm%o)      ! src lower index
      us = floor  ((fu-d)/sm%d + sm%o)      ! src upper index
    else
      !---------------------------------------------------------------------------
      ! Restrict operation:  The source cells are smaller, and are hence closer
      ! to the cell corners.
      !---------------------------------------------------------------------------
      fl = (lt-tm%o)*tm%d                           ! floating lower
      fu = (ut-tm%o)*tm%d                           ! floating upper
      ls = floor  ((fl-d)/sm%d + sm%o)              ! src lower index
      us = ceiling((fu-d)/sm%d + sm%o)              ! src upper index
    endif

    ls = max(ls,sm%li)                            ! bounds
    us = min(us,sm%ui)                            ! bounds
    lt = max(lt,tm%lb)                            ! bounds
    ut = min(ut,tm%ub)                            ! bounds
    if (self%verbose > 2) then
      write (stdlog,'(a,2i6,4(2x,3i7))') 'guards_t%same: tgt,src,ls,us,lt,ut=', &
        tgt%id, src%id, ls, us, lt, ut
      call tgt%stats ('guards_t%same tgt')
      call src%stats ('guards_t%same src')
    end if
    call io%assert (all(lt >= tm%lb), 'guards_t%same: lt breach')
    call io%assert (all(ut <= tm%ub), 'guards_t%same: ut breach')
    call io%assert (all(ls >= sm%lb), 'guards_t%same: ls breach')
    call io%assert (all(us <= sm%ub), 'guards_t%same: us breach')

    if (src%n(3)==1) then
      ls(3) = 1
      us(3) = 1
    endif

    if (src%level == tgt%level) then
      call self%no_t (tgt, src, ls, us, src%mem, tgt%mem, lt, ut, only=only)
    elseif (src%level == tgt%level-1) then
      call self%interp_to_higher_level (tgt, src, ls, us, src%mem, tgt%mem, lt, ut, only=only)
    elseif (src%level == tgt%level+1) then
      call self%binning_to_lower_level (tgt, src, ls, us, src%mem, tgt%mem, lt, ut, only=only)
    endif

    if (associated(tgt%filled)) &
    tgt%filled (lt(1):ut(1),lt(2):ut(2),lt(3):ut(3)) = 10*src%level
    call trace%end (itimer=itimer)
  END SUBROUTINE same

  !===============================================================================
  !> Not interpolating in time, but simply copying over the relevant source region
  !===============================================================================
  SUBROUTINE no_t (self, target, source, ls, us, in, out, lt, ut, only, lock)
    class(m3d_halo_t)      :: self
    class(patch_t)                  :: source, target
    integer, dimension(3)           :: ls, us
    real(4), dimension(:,:,:,:,:,:) :: in 
    real(4), dimension(:,:,:,:,:,:) :: out
    integer, dimension(3)           :: lt, ut
    integer, optional               :: only
    type(lock_t), optional          :: lock
    !.............................................................................
    integer, dimension(3)           :: l, u, o
    integer                         :: i1, i2, i3, i4, i5, it, id
    integer                         :: iv0, iv1
    type(lock_t), pointer           :: slock
    integer, save                   :: itimer=0, ntimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('interp_t%no_t', itimer=itimer)
    if (present(only)) then
      iv0 = only
      iv1 = only
    else
      iv0 = 1
      iv1 = target%nv
    end if
    !---------------------------------------------------------------------------
    ! Shifted from source index space
    !---------------------------------------------------------------------------
    if (present(lock))                                      call lock%set ('no_t')
    !.........................................................................
    id = target%idx%d
    o = lt-ls
    do i5=1,source%nw
    do i4=iv0,iv1                                         ! loop over variables
    do i3=ls(3),us(3)
    do i2=ls(2),us(2)
    !$omp simd
    do i1=ls(1),us(1)
      out(i1+o(1),i2+o(2),i3+o(3),i4,1,i5) = in(i1,i2,i3,i4,1,i5)
    end do
    end do
    end do
    end do
    end do
    if (present(lock))                                    call lock%unset ('no_t')
    !.............................................................................
    call trace%end (itimer)
  END SUBROUTINE no_t

  !===============================================================================
  SUBROUTINE interp_to_higher_level (self, target, source, ls, us, in, out, lt, ut, only, lock)
    class(m3d_halo_t)      :: self
    class(patch_t)                  :: source, target
    integer, dimension(3)           :: ls, us
    real(4), dimension(:,:,:,:,:,:) :: in 
    real(4), dimension(:,:,:,:,:,:) :: out
    integer, dimension(3)           :: lt, ut
    integer, optional               :: only
    type(lock_t), optional          :: lock
    !.............................................................................
    integer, dimension(3)           :: l, u, o, i, sn, tn
    real(4), dimension(:,:), allocatable :: tmp_in
    real(4), dimension(:,:), allocatable :: tmp_out
    real(4), dimension(3,3)         :: kernel
    integer                         :: i1, i2, i3, i4, i5, it
    integer                         :: iv0, iv1, o1, o2, o3, k1, k2
    type(lock_t), pointer           :: slock
    integer, save                   :: itimer=0, ntimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_halo_t%interp_to_higher_level', itimer=itimer)
    if (present(only)) then
      iv0 = only
      iv1 = only
    else
      iv0 = 1
      iv1 = target%nv
    end if
    sn = source%n
    tn = target%n
    allocate(tmp_in (sn(1), sn(2)))
    allocate(tmp_out(tn(1), tn(2)))
    kernel(:,1) = [0.25,0.50,0.25]
    kernel(:,2) = [0.50,1.00,0.50]
    kernel(:,3) = [0.25,0.50,0.25]
    !---------------------------------------------------------------------------
    ! Shifted from source index space
    !---------------------------------------------------------------------------
    if (present(lock))                                      call lock%set ('no_t')
    !.........................................................................
    o = lt-ls
    do i5=1,source%nw
    do i4=iv0,iv1                                         ! loop over variables
    do i3=ls(3),us(3)
      o3 = i3+o(3)
      tmp_in  = in(ls(1):us(1),ls(2):us(2),i3,i4,1,i5)
      tmp_out = 0
      do i2=1,sn(2)
      do i1=1,sn(1)
        do k2=1,3
          o2 = i2*2 + k2 - 3
          if (o2 < 1) o2 = tn(2) + o2
          do k1=1,3
            o1 = i1*2 + k1 - 3
            if (o1 < 1) o1 = tn(1) + o1
            tmp_out(o1,o2) = tmp_out(o1,o2) + kernel(k1,k2)*tmp_in(i1,i2)
          end do
        end do
      end do
      end do
      out(lt(1):ut(1),lt(2):ut(2),o3,i4,1,i5) = tmp_out
    end do
    end do
    end do
    if (present(lock))                                    call lock%unset ('no_t')
    !.............................................................................
    call trace%end (itimer)
  END SUBROUTINE interp_to_higher_level

  !===============================================================================
  SUBROUTINE binning_to_lower_level (self, target, source, ls, us, in, out, lt, ut, only, lock)
    class(m3d_halo_t)      :: self
    class(patch_t)                  :: source, target
    integer, dimension(3)           :: ls, us
    real(4), dimension(:,:,:,:,:,:) :: in 
    real(4), dimension(:,:,:,:,:,:) :: out
    integer, dimension(3)           :: lt, ut
    integer, optional               :: only
    type(lock_t), optional          :: lock
    !.............................................................................
    integer, dimension(3)           :: l, u, i, n
    real(4), dimension(:,:), allocatable :: tmp_in
    real(4), dimension(:,:), allocatable :: tmp_out
    real(4), dimension(3,3)         :: kernel
    integer                         :: i1, i2, i3, i4, i5, it
    integer                         :: iv0, iv1, k1, k2
    type(lock_t), pointer           :: slock
    integer, save                   :: itimer=0, ntimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_halo_t%binning_to_lower_level', itimer=itimer)
    if (present(only)) then
      iv0 = only
      iv1 = only
    else
      iv0 = 1
      iv1 = target%nv
    end if
    n = us-ls+1
    allocate(tmp_in (n(1), n(2)))
    allocate(tmp_out(ut(1)-lt(1)+1, ut(2)-lt(2)+1))
    !---------------------------------------------------------------------------
    ! Shifted from source index space
    !---------------------------------------------------------------------------
    if (present(lock))                                      call lock%set ('no_t')
    !.........................................................................
    do i5=1,source%nw
    do i4=iv0,iv1                                         ! loop over variables
      tmp_out = 0
      do i3=ls(3),us(3)
        tmp_in  = in(ls(1):us(1),ls(2):us(2),i3,i4,1,i5)
        tmp_out = tmp_out + tmp_in(1:n(1):2, 1:n(2):2)
        tmp_out = tmp_out + tmp_in(1:n(1):2, 2:n(2):2)
        tmp_out = tmp_out + tmp_in(2:n(1):2, 1:n(2):2)
        tmp_out = tmp_out + tmp_in(2:n(1):2, 2:n(2):2)
      end do
      out(lt(1):ut(1),lt(2):ut(2),lt(3),i4,1,i5) = tmp_out / (4*n(3))
    end do
    end do
    if (present(lock))                                    call lock%unset ('no_t')
    !.............................................................................
    call trace%end (itimer)
  END SUBROUTINE binning_to_lower_level

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
        if (task%verbose>1) then
          write(stdout,*) task%id, ' has finished!', task%task_list%nfinished, '/', task%task_list%na, timer%done
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