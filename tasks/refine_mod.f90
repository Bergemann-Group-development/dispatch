!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> This module handles checking max change between neighboring points.  Each
!> instance of it needs an index to the variable, whether to take the log, and
!> a value for the max allowed change.  It might be applied to log density, log
!> pressure, log opacity, velocity amplitude, magnetic field amplitude, etc.
!>
!> task_list_t%update
!>   refine_t%check_current
!>     if (associated(task%parent)) &
!>       task%parent%lock%set ('parent')               + parent locked
!>     task%lock%set ('task')                          | + child locked
!>     refine_t%refine_needed                          | !
!>       refine_t%check_Jeans                          | |
!>     ...                                             | |
!>     tasklist%lock%set ('split_patch')               | | + tasklist locked
!>     split_patch                                     | | |
!>       ...                                           | | |
!>       tasklist%init_all_nbors                       | | |
!>       tasklist%reset_status                         | | |
!>       child%print_nbors                             | | |
!>       ...                                           | | |
!>     check_level_constraints                         | | |
!>     tasklist%lock%unset ('split_patch')             | | - tasklist unlocked
!>     ...                                             | |
!>     tasklist%lock%set ('remove_patch')              | | + tasklist locked
!>     remove_patch                                    | | |
!>       tasklist%remove_link (this)                   | | |
!>       tasklist%init_nbor_nbors (this)               | | |
!>       tasklist%check_nbor_nbors (this)              | | |
!>     check_level_constraints                         | | !
!>     tasklist%lock%unset ('remove_patch')            | | - tasklist unloccked
!>     ...                                             | |
!>     if (associated(task%parent)) &                  | |
!>     task%lock%unset ('task')                        ! - child locked
!>       task%parent%lock%unset ('parent')             - parent locked
!>   ...
!>
!> How can locking prevent a download from finding that a patch has disappeared?
!> If download causes a lock on the target link, that is enough to prevent any
!> nbor from refining, if and only if refining and derefining requires a lock on
!> all nbors of the affected patch.
!>
!> Special attention is required for the parent task, which must not be allowed
!> to download from interior patches that have disappeared in the mean time.
!> Hence a download to the parent patch must be prevented from even starting,
!> while some child patch is being refined or derefined
!>
!> If a download to the parent patch has started, then it is locked by that
!> thread, and therefore nothing should be able to change.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE refine_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE mpi_mod
  USE bits_mod
  USE kinds_mod
  USE link_mod
  USE list_mod
  USE task_mod
  USE gpatch_mod
  USE solver_mod
  USE extras_mod
  USE scene_mod
  USE download_mod
  USE counters_mod
  implicit none
  private
  !
  type:: refine_t
    integer:: levelmin=1                      ! Min level of refinment
    integer:: levelmax=1                      ! Max level of refinment
    integer:: ratio
  contains
    procedure:: init
    procedure:: refine_needed
    procedure:: check_current
    procedure:: setup_mask
    procedure:: check_cubical
    procedure:: check_jeans
    procedure:: shock_and_cd_detector
    procedure:: refine_region
    procedure:: split_patch
    procedure:: remove_patch
    procedure:: gradient_detector
    procedure:: check_nbor_support
  end type
  !
  integer:: check_interval=10                 ! Number of steps betwen checks
  real:: min_dx=0.0, max_dx=1e32              ! Min and max grid size
  real:: refine_from_t=0.0                    ! Time limits for occurence of refinement
  real:: refine_until_t=1e30                  ! Time limits for occurence of refinement
  real:: max_jeans=0.0                        ! Jeans number below which a patch is flagged for refinement
  real:: min_jeans=0.0                        ! Jeans number above which a patch is flagged for derefinement
  real:: max_shock=0.0, max_contact=0.0       ! Tolerance values for shock/discontinuity detector refinement
  real:: min_shock=0.0, min_contact=0.0       ! Tolerance values for shock/discontinuity detector derefinement
  integer, parameter:: ngradvar=10            ! Number of gradient variables
  real, dimension(3):: region_llc=0.0         ! Sometimes you want to refine a region (manually)
  real, dimension(3):: region_urc=0.0
  real, dimension(ngradvar):: max_grad=0.0    ! Tolerance values for gradient detector refinement
  real, dimension(ngradvar):: min_grad=0.0    ! Tolerance values for gradient detector derefinement
  character(len=32), dimension(ngradvar):: grad_var=""
  logical:: on=.false.                        ! Flag to turn the entire process on or off 
  logical:: force_cubical=.false.             ! Force cubical patches
  integer:: verbose=0                         ! Noise level
  type(refine_t), public:: refine
CONTAINS

!===============================================================================
!> Initialize refinement parameters
!===============================================================================
SUBROUTINE init (self)
  class(refine_t):: self
  integer:: iostat
  integer, save:: levelmin, levelmax
  integer, save:: ratio=2                     ! How does a patch split when it is refined
  logical, save:: first_time = .true.
  namelist /refine_params/ on, levelmin, levelmax, check_interval, ratio, &
    refine_from_t, refine_until_t, min_dx, max_dx, max_jeans, min_jeans, &
    max_shock, min_shock, max_contact, min_contact, max_grad, min_grad, &
    grad_var, region_llc, region_urc, force_cubical, verbose
  character(len=120):: id = &
    '$Id$ tasks/refine_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%init')
  call trace%print_id (id)
  !$omp critical (refine_cr)
  if (first_time) then
    first_time=.false.
    levelmin = self%levelmin
    levelmax = self%levelmax
    rewind (io%input)
    read (io%input, refine_params, iostat=iostat)
    if (io%master) then
      if (iostat > 0) call io%namelist_warning ('refine_params')
      write (*, refine_params)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! If AMR is active, we must use check_filled, since nbor levels are in
  ! decreasing order, to save computing time
  !-----------------------------------------------------------------------------
  if (levelmax > levelmin) &
    download%check_filled = .true.
  !$omp end critical (refine_cr)
  self%levelmin = levelmin
  self%levelmax = levelmax
  self%ratio = ratio
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Check if refinement is desired and needed for the current task
!===============================================================================
FUNCTION refine_needed (self, tasklist, link) RESULT (refine)
  class(refine_t):: self
  class(list_t):: tasklist
  class(link_t), pointer:: link
  integer:: refine
  !.............................................................................
  class(gpatch_t), pointer:: patch
  real:: dxmin, dxmax
  integer, save:: itimer=0
  !logical:: new_sink
  !----------------------------------------------------------------------------
  refine = 0
  if (.not. on) return
  call trace%begin('refine_t%needed', itimer=itimer)
  !----------------------------------------------------------------------------
  patch => gpatch%cast2gpatch (link%task)
  dxmin = minval(patch%ds)
  dxmax = maxval(patch%ds)
  if (verbose > 2) then
    write (io_unit%output,'(a,i6,2i4,l4,1p,4e11.2)') &
      'refine_t%refine_needed: id, istep, refine_next =', &
      patch%id, patch%istep, patch%check_refine_next, &
      patch%is_set(bits%virtual), min_dx, dxmin, max_dx, dxmax
  end if
  !----------------------------------------------------------------------------
  ! Filter out patches that do not match the necessary criteria
  !----------------------------------------------------------------------------
  if (patch%is_set(bits%virtual)) then
    call trace%end (itimer)
    return
  end if
  if (patch%istep < patch%check_refine_next) then
    call trace%end (itimer)
    return
  end if
  if (dxmin < min_dx .or. dxmax >= max_dx) then
    call trace%end (itimer)
    return
  end if
  if (patch%time < refine_from_t .or. patch%time > refine_until_t) then
    call trace%end (itimer)
    return
  end if
  !------------------------------------------------------------------------
  ! Refinement criteria, starting with assumed derefinement
  !------------------------------------------------------------------------
  refine = -1
  !
  ! Extras refinement
  refine = max(refine, extras%check_refine (patch))
  !
  ! Ensure consistent, cubical patch size (e.g. after RAMSES startup)
  if (force_cubical) refine = max (refine, self%check_cubical(patch))
  !
  ! Jeans criterion
  if (max_jeans > 0.0) refine = max (refine, self%check_Jeans(patch))
  !
  ! shock and/or contact discontinuity detector
  if (max_shock > 0.0 .or. max_contact > 0.0) &
    refine = max (refine, self%shock_and_cd_detector(patch))
  !
  ! gradient detector
  if (any(max_grad > 0.0)) refine = max (refine, self%gradient_detector(patch))
  !
  ! refine on (a) user-specified region(s)
  if (any(region_llc /= 0.0) .and. any(region_urc /= 0.0)) &
    refine = max (refine, self%refine_region(patch))
  !
  ! Check that neighbors of `self` are not more than one level different.
  ! If they are, mark this patch for splitting so that support can be
  ! provided to the neighbors of `patch`.
  refine = max (refine ,self%check_nbor_support (tasklist, link))
  !
  patch%check_refine_next = patch%istep + check_interval
  call trace%end (itimer)
END FUNCTION refine_needed

!===============================================================================
!> Check if refinement is desired and needed for the current task; if so, push
!> new tasks onto the (ready) queue with the same time as the parent, i.e., to
!> the head of the queue.  If the task is a leaf task, and the level of refinement
!> does not appear to be needed, request derefinement.
!===============================================================================
SUBROUTINE check_current (self, tasklist, link, was_refined, was_derefined)
  class(refine_t):: self
  class(list_t):: tasklist
  class(link_t), pointer:: link, nbor, first_nbor
  class(gpatch_t), pointer:: patch
  logical, intent(out):: was_refined, was_derefined
  real:: dxmin, dxmax
  integer:: refine, id
  integer, save:: itimer=0
  !.............................................................................
  was_refined = .false.
  was_derefined = .false.
  if (.not. on) return
  if (verbose > 2) &
    print *, link%task%id, 'refine_t%check_current, at t =', &
      link%task%time, link%task%level, self%levelmax
  !---------------------------------------------------------------------------
  call trace%begin('refine_t%check_current', itimer=itimer)
  patch => gpatch%cast2gpatch (link%task)
  id = patch%id
  !---------------------------------------------------------------------------
  ! Lock the task, and its parent task, if any, for check_refine. This prevents
  ! the parent and child from considering refine/derefine simultaneously.  If
  ! the thread working on the parent gets a lock on it, then the child will
  ! have to wait, and in the mean time it should not be locked, to allow for
  ! example download access to the child while it is waiting. 
  !
  ! If the thread working on the parent gets a lock on it, this forces any
  ! thread that works on a child patch to wait here, until it can get a lock
  ! on the parent link.
  !
  ! Note that the thread that comes here does NOT (yet) have a lock on the task
  ! list, so it should safely be able to put a lock on its own task.
  !
  ! To prevent confusion, a thread calling download_t%download_link should
  ! check if any source task has a parent, and if so acquire a lock on it, to
  ! before proceeding to go through the nbor list.  This should prevent the
  ! confusion that would arise if the nbor list was modified during the
  ! download.
  !---------------------------------------------------------------------------
  !if (associated(patch%parent)) then
  !  associate (parent => patch%parent)
  !  select type (parent)
  !  class is (gpatch_t)
  !  call parent%link%lock%set ('parent')
  !  end select
  !  end associate
  !end if
  !call patch%link%lock%set ('check_current')
  !-----------------------------------------------------------------------------
  if (patch%is_clear(bits%not_leaf)) then
    if (verbose > 2) &
      tasklist%verbose = verbose
    !---------------------------------------------------------------------------
    ! Setup mask, and check if refinement or derefinement is needed
    !---------------------------------------------------------------------------
    if (patch%level < self%levelmax) & 
      call self%setup_mask (patch)
    refine = self%refine_needed (tasklist, link)
    if (associated(patch%mask)) &
      deallocate (patch%mask)
    !---------------------------------------------------------------------------
    ! Optional custom refine
    !---------------------------------------------------------------------------
    if (patch%manual_refine) then
      call patch%custom_refine
      call check_level_constraint
    end if
    !---------------------------------------------------------------------------
    ! A negative refine value indicates that the patch is not needed
    !---------------------------------------------------------------------------
    if (refine < 0 .and. patch%level > self%levelmin) then
      if (tasklist%verbose > 2) &
        write (io%output,*) wallclock(),' thread',omp%thread,' waitfor tasklist(d)'
      if (tasklist%verbose > 2) &
        write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(d)'
      if (verbose > 0) &
        write (io%output,'(a,i6,i3,2x,a)') &
          'id, level =', patch%id, patch%level, 'derefine'
      !-------------------------------------------------------------------------
      ! The task_list must be locked while the patch is removed, and while the
      ! question of consistent support is considered.  The nbors of the task
      ! must also be locked, to prevent downloads from happening, until the
      ! new nbor lists are in place.  There is no need to lock the link itself,
      ! since it, and the parent link, are already locked.
      !-------------------------------------------------------------------------
      !$omp critical (refine_nbor_cr)
      first_nbor => link%nbor
      nbor => first_nbor
      call link%lock%set ('refine_t%split_patch link tasks')
      do while (associated(nbor))
        call nbor%link%lock%set ('refine_t remove_patch parent nbor tasks')
        nbor => nbor%next
      end do
      !-------------------------------------------------------------------------
      call tasklist%lock%set ('refine_t remove_patch tasklist')
      call self%remove_patch (tasklist, link)
      call check_level_constraint
      call tasklist%lock%unset ('refine_t remove_patch tasklist')
      !-------------------------------------------------------------------------
      ! The nbor lists are now different, in that they lack the removed task
      ! (their locks are unlocked before removal), but the ones that remain
      ! need to be unlocked
      !-------------------------------------------------------------------------
      nbor => first_nbor
      do while (associated(nbor))
        !write (io_unit%output,'(a,2i6,i3,i4)') &
        !  'unlock link: task, lock, level, thread =', &
        !  nbor%task%id, nbor%link%lock%id, nbor%link%lock%level, omp%thread
        !flush (io_unit%output)
        call nbor%link%lock%unset ('refine_t remove_patch parent nbor tasks')
        nbor => nbor%next
      end do
      if (associated(link)) &
        call link%lock%unset ('refine_t%split_patch link tasks')
      !$omp end critical (refine_nbor_cr)
      !-------------------------------------------------------------------------
      if (tasklist%verbose > 2) &
        write (io%output,*) wallclock(),' thread',omp%thread,' unlocked tasklist(d)'
      was_derefined = .true.
    !---------------------------------------------------------------------------
    ! A positive refine value indicates that the patch needs to be refined
    !---------------------------------------------------------------------------
    else if (refine > 0 .and. patch%level < self%levelmax) then
      if (tasklist%verbose > 2) &
        write (io%output,*) wallclock(),' thread',omp%thread,' waitfor tasklist(r)'
      if (tasklist%verbose > 2) &
        write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(r)'
      if (verbose > 0) &
        write (io%output,'(a,i6,i3,2x,a)') &
          'id, level =', patch%id, patch%level, 'refine'
      !-------------------------------------------------------------------------
      ! The task_list must be locked while the patch is removed, and while the
      ! question of consistent support is considered.  Also, the current links
      ! all nbor links must be locked during refine (the current link itself is
      ! already locked). Such a lock of all nbors must be completed before the
      ! next thread is allowed to try the same, in case they have competing
      ! interests.
      !-------------------------------------------------------------------------
      !$omp critical (refine_nbor_cr)
      nbor => link%nbor
      call link%lock%set ('refine_t%split_patch link tasks')
      do while (associated(nbor))
        call nbor%link%lock%set ('refine_t%split_patch nbor tasks')
        nbor => nbor%next
      end do
      !-------------------------------------------------------------------------
      call tasklist%lock%set ('refine_t split_patch tasklist')
      call self%split_patch (tasklist, link)
      call check_level_constraint
      call tasklist%lock%unset ('refine_t split_patch tasklist')
      !-------------------------------------------------------------------------
      ! The nbor list now also contains the new child links, which should be
      ! locked when created.  All nbor links are unlocked here.
      !-------------------------------------------------------------------------
      nbor => link%nbor
      do while (associated(nbor))
        call nbor%link%lock%unset ('refine_t%split_patch nbor tasks')
        nbor => nbor%next
      end do
      call link%lock%unset ('refine_t%split_patch link tasks')
      !$omp end critical (refine_nbor_cr)
      !-------------------------------------------------------------------------
      if (tasklist%verbose > 2) &
        write (io%output,*) wallclock(),' thread',omp%thread,' unlocked tasklist(r)'
      was_refined = .true.
    end if
    if (verbose > 2) &
      tasklist%verbose = 0
  end if
  !
  !---------------------------------------------------------------------------
  ! Unlock first the parent task, and then the task, if it still exists
  !---------------------------------------------------------------------------
  if (associated(link)) then
  !  call link%lock%unset ('check_current')
  !  if (associated(link%task%parent)) then
  !    associate (parent => link%task%parent)
  !    select type (parent)
  !    class is (gpatch_t)
  !    call parent%link%lock%unset ('parent')
  !    end select
  !    end associate
  !  end if
  end if
  call trace%end (itimer)
contains
  !=============================================================================
  !> Iteratively check the nbor level constraint, splitting patches until the
  !> there are no more violations.  Since this inserts new patches in the patch
  !> list, and modifies the queue, it must be done in a critical region, blocking
  !> other updates to the patch list.
  !=============================================================================
  subroutine check_level_constraint
    class(link_t), pointer:: link
    integer:: iteration
    logical:: still_needed
    !---------------------------------------------------------------------------
    iteration = 0
    still_needed = .true.
    do while (still_needed)
      still_needed = .false.
      link => tasklist%head
      do while (associated(link))
        if (link%task%is_clear (bits%not_leaf)) then
          if (self%check_nbor_support (tasklist, link) > 0) then
            if (verbose > 3) &
              print *, link%task%id, &
              'check_level_constraint: patch needs splitting at level', &
              link%task%level
            call self%split_patch (tasklist, link)
            still_needed = .true.
            exit
          end if
        end if
        link => link%next
      end do
      iteration = iteration+1
      if (verbose > 3) &
        print *, 'check_level_constraint, iteration =', &
        iteration, still_needed
    end do
  end subroutine check_level_constraint
END SUBROUTINE check_current

!===============================================================================
!> Setup a mask with false values where a higher level patch overlaps this patch.
!> FIXME: Currently ignores overlapping patches with centers outside this patch.
!===============================================================================
SUBROUTINE setup_mask (self, patch)
  class(refine_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  class(link_t), pointer:: nbor
  integer, dimension(3):: m, l, u
  !-----------------------------------------------------------------------------
  ! Check for higher level patches inside this
  !-----------------------------------------------------------------------------
  nbor => patch%link%nbor
  do while (associated(nbor))
    if (patch%contains(nbor%task%position)) then
      if (.not.associated(patch%mask)) then
        l = patch%mesh%li
        u = patch%mesh%ui
        allocate (patch%mask(l(1):u(1),l(2):u(2),l(3):u(3)))
        patch%mask = .true.
      end if
      associate (nbpatch => nbor%task)
      select type (nbpatch)
      class is (gpatch_t)
      l = patch%index_only (nbpatch%position + nbpatch%mesh%lf)
      u = patch%index_only (nbpatch%position + nbpatch%mesh%uf)
      l = max(l,patch%mesh%li)
      u = min(u,patch%mesh%ui)
      end select
      end associate
      patch%mask(l(1):u(1),l(2):u(2),l(3):u(3)) = .false.
    end if
    nbor => nbor%next
  end do
END SUBROUTINE setup_mask

!===============================================================================
!> Create a new patch using an existing "parent" patch.
!> This procedure is used frequently during refinement.
!===============================================================================
SUBROUTINE make_child_patch (self, tlist, parent_link, parent, pos, child_size, i)
  class(refine_t):: self
  class(list_t) :: tlist
  class(link_t), pointer :: parent_link, child_link, child_nbor
  class(gpatch_t), pointer :: parent, child
  real(8), dimension(3), intent(in) :: pos, child_size
  integer:: i, j, count, l(3), u(3)
  !------------------------------------------------------------------------
  call trace%begin('refine_mod::make_child_patch')
  !------------------------------------------------------------------------
  ! Allocate a new task, 1/nsplit the size of original one, and adjust its
  ! position in each direction with "pos", in units of the new size.
  !------------------------------------------------------------------------
  ! Using sourced allocation below effectively clones the parent.  This is
  ! necessary, to pick up properties that may have been set at the experiment
  ! level, which would otherwise be left undefined, or in inconsistent state
  !------------------------------------------------------------------------
  allocate (child, source=parent)
  child%mem_allocated = .false.
  nullify (child%mem)
  nullify (child%mesh)
  child%id = 0                                       ! to ensure it gets its own ID
  child%size = child_size                            ! new size
  child%position = pos                               ! patch position
  child%llc_cart = child%position - merge(0.5,0.0,child%n>1) * child%size  ! lower left corner
  child%llc_nat = child%llc_cart                     ! FIXME!: this form only applies for Cartesian coords.
  !------------------------------------------------------------------------
  ! Force consistent resolution in new patches, overriding the resolution
  ! inherited from the parent patch.
  !------------------------------------------------------------------------
  child%n = scene%patch_n
  !------------------------------------------------------------------------
  ! Reset flags (FIXME: Would it be better to just reset `status` to zero?)
  !------------------------------------------------------------------------
  call child%clear(bits%root_grid)
  call child%clear(bits%ready)
  call child%clear(bits%busy)
  child%boundaries%bits = 0
  !------------------------------------------------------------------------
  ! Initialise the new patch
  !------------------------------------------------------------------------
  child%level = parent%level + 1                     ! must be set before init!
  call child%init()                                  ! patch, for a new id
  child%level = parent%level + 1                     ! override possible calculation in `patch_t%init`.
  child%time  = parent%time
  child%it    = 1
  child%iit   = 1
  child%dt(:) = parent%dt(parent%it) / self%ratio    ! reset the time-step
  child%t(:) = child%time
  child%parentid = parent%id                         ! for debugging
  child%check_refine_next = 0                        ! allow recursive refine
  !------------------------------------------------------------------------
  ! Immediately add output of the new patch to the previous snapshot, for
  ! debugging purposes.
  !------------------------------------------------------------------------
  !if (child%restart <= 0) then
  !  child%out_next = child%time                      ! force output now
  !  child%iout = max(0,child%iout - 1)               ! correct output counter
  !end if
  !------------------------------------------------------------------------
  ! Increment I/O write counter, to avoid repeated writing of metadata
  !------------------------------------------------------------------------
  io%nwrite = io%nwrite + 1
  count = counters%increment (parent%iout+1, io%nwrite)
  !------------------------------------------------------------------------
  ! Append links to the new patch in the patch list, sorted by quality
  ! Neighbour lists are constructed at the end of `split_patch`.
  !------------------------------------------------------------------------
  allocate (child_link)
  call child_link%init                               ! initialise the lock
  call child_link%lock%set ('child init')            ! immediately lock it
  child_link%task => child                           ! associate task with link
  child_link%parent => parent_link                   ! used in remove_patch
  child%link => child_link                           ! associate link with task
  child%parent => parent                             ! task relation
  call tlist%append_link (child_link)                ! which is added to the patch list
  io%ntask = tlist%na                                ! update the task count
  !------------------------------------------------------------------------
  ! Make parent a neighbour of new child patch.
  ! This will be overwritten in `init_nbor_nbors` but is necessary here to permit
  ! interpolation of values from parent to child.
  !------------------------------------------------------------------------
  allocate (child_nbor)
  child_nbor%link => parent_link
  child_nbor%task => parent
  child_link%nbor => child_nbor
  !------------------------------------------------------------------------
  ! Download values from parent and neighbors, except at time = 0, to allow
  ! recursive refinement there.
  ! If time = 0, patch should have been filled with values during `child%init()`
  ! call.
  !------------------------------------------------------------------------
  if (parent%time > 0.0 .and. child%restart <= 0) &
    call download%download_link (child_link, all_cells=.true.)
  if (verbose > 1) then
    l(:) = parent%mesh(:)%li
    u(:) = parent%mesh(:)%ui
    if (i==0) &
    print 1,'refine_t%make_child_patch: parent, level, iout, time, size(3), position(3) =', &
      parent%id, parent%level, parent%iout, parent%time, parent%size, parent%position , &
      maxval(parent%mem(l(1):u(1),l(2):u(2),l(3):u(3),parent%idx%d,parent%it,1))
    1 format(a,3i6,f10.6,2(2x,3f8.4),1pe12.2)
    print 1,'                            child, level, iout, time, size(3), position(3) =', &
      child%id, child%level, child%iout, child%time, child%size, child%position , &
      maxval(child%mem(l(1):u(1),l(2):u(2),l(3):u(3),child%idx%d,child%it,1))
  end if
  call trace%end()
END SUBROUTINE make_child_patch

!===============================================================================
!> Split a patch into two (parameter) in each direction, merging the new patches
!> into the list.  Note that the link to the original patch is kept, with a "leaf"
!> flag being cleared, indicating it is not longer a "leaf" patch.
!===============================================================================
SUBROUTINE split_patch (self, tlist, this)
  class(refine_t):: self
  class(list_t):: tlist
  class(link_t), pointer:: this, link, tail
  class(gpatch_t), pointer:: patch
  integer:: ix, iy, iz
  integer:: nsplit(3), i, j(3), j1(3), j2(3)
  real(8), dimension(3):: pos, rpos, split_size
  logical:: do_trace
  !.............................................................................
  if (this%task%is_set(bits%not_leaf)) return
  do_trace = io%do_trace
  if (verbose > 3) &
    io%do_trace=.true.
  call trace%begin('refine_t%split_patch', 1)
  !-----------------------------------------------------------------------------
  ! Save the current task list tail, which is where the child tasks will be 
  ! appended
  !-----------------------------------------------------------------------------
  tail => tlist%tail
  !-----------------------------------------------------------------------------
  ! Decide how to split
  !-----------------------------------------------------------------------------
  patch => gpatch%cast2gpatch (this%task)
  if (patch%msplit>0) then
    nsplit = patch%msplit
  else
    nsplit = self%ratio
  end if
  where (patch%size==0 .or. patch%n==1)
    nsplit = 1
  end where
  !-----------------------------------------------------------------------------
  ! Perform the split
  !-----------------------------------------------------------------------------
  split_size = patch%size / nsplit
  i = 0
  do ix=1,nsplit(1)
    rpos(1) = (ix-(nsplit(1)+1)*0.5) * split_size(1)
    do iy=1,nsplit(2)
      rpos(2) = (iy-(nsplit(2)+1)*0.5) * split_size(2)
      do iz=1,nsplit(3)
        rpos(3) = (iz-(nsplit(3)+1)*0.5) * split_size(3)
        pos = rpos + patch%position
        j = patch%index_only(pos)
        j1 = j - patch%mesh%n/4 + 1
        j2 = j1 + patch%mesh%n/2 - 1
        !write (io_unit%output,'(a,i6,1p,2g15.4,0p,3(2x,3f9.4),2x,3(i4,":",i2))') 'xxx', &
        !  patch%id, patch%fmaxval(patch%mem(:,:,:,patch%idx%d,patch%it,1)), &
        !  maxval(patch%mem(j1(1):j2(1),j1(2):j2(2),j1(3):j2(3),patch%idx%d,patch%it,1)), &
        !  rpos, pos, patch%position, &
        !  j1(1), j2(1), j1(2), j2(2), j1(3), j2(3)
        write (io_unit%output,'(i4,1p,2g14.4,0p,2x,3(i4,":",i2),2x,a)') &
          patch%id, patch%fmaxval(patch%mem(:,:,:,patch%idx%d,patch%it,1)), &
          maxval(patch%mem(j1(1):j2(1),j1(2):j2(2),j1(3):j2(3),patch%idx%d,patch%it,1)), &
          j1(1), j2(1), j1(2), j2(2), j1(3), j2(3), 'xxx'
        if (associated(patch%mask) .and. &
         all(j >= patch%mesh%li .and. j <= patch%mesh%ui)) then
          if (patch%mask(j(1),j(2),j(3))) then
            call make_child_patch(self, tlist, this, patch, pos, split_size, i)
            i = i+1
          end if
        else
          call make_child_patch(self, tlist, this, patch, pos, split_size, i)
          i = i+1
        end if
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! At this point patch is (possibly again) not a leaf patch, although that
  ! may change after the next derefine opportunity
  !-----------------------------------------------------------------------------
  call patch%set(bits%not_leaf)
  !-----------------------------------------------------------------------------
  ! Initialize the nbor lists of the parent patch and all its nbors, based on
  ! positions, then reset their status based on the nbor lists, and count
  !-----------------------------------------------------------------------------
  call tlist%init_nbor_nbors (this)
  call tlist%reset_status
  !-----------------------------------------------------------------------------
  ! Also, expensive way to make sure the state of readiness is consistent (FIXME)
  !-----------------------------------------------------------------------------
  link => tlist%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call tlist%check_ready (link)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Optionally, print children task list info
  !-----------------------------------------------------------------------------
  if (verbose > 2) then
    write (io%output,*) 'child tasks:'
    link => tail
    do while (associated(link))
      call link%print_nbors
      link => link%next
    end do
  end if
  !-----------------------------------------------------------------------------
  ! FIXME: special care may be needed in multi-rank context
  !-----------------------------------------------------------------------------
  call tlist%init_bdries
  !
  if (verbose > 3) call tlist%print_tasks()
  !
  call trace%end()
  io%do_trace = do_trace
END SUBROUTINE split_patch

!===============================================================================
!> Remove a patch that has been deemed not necessary.  At this point, the task
!> is not in the ready queue, but is still in the task list.  The proper order
!> of removal is then:
!>
!>   * remove the task link from the task list
!>   * update the nbor lists of all its nbors
!>   * deallocate the task memory
!>   * deallocate the task
!>   * deallocate the nbor list
!>   * deallocate the task link
!===============================================================================
SUBROUTINE remove_patch (self, tlist, this)
  class(refine_t):: self
  class(list_t):: tlist
  class(link_t), pointer:: this, link, nbor, parent
  class(task_t), pointer:: task
  integer:: count, io_verbose
  !-----------------------------------------------------------------------------
  ! Start processing, blocking (at least for now) all queue-related activity
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%remove_patch')
  !-----------------------------------------------------------------------------
  ! Send a suicide note to virtual nbors, if this is a boundary patch
  !-----------------------------------------------------------------------------
  if (this%task%is_set (bits%boundary)) then
    call this%task%set (bits%remove)
    call this%send_to_vnbors
  end if
  !-----------------------------------------------------------------------------
  ! Change the status of the parent patch to leaf
  !-----------------------------------------------------------------------------
  parent => this%parent
  call parent%task%clear (bits%not_leaf)
  !-----------------------------------------------------------------------------
  ! Decrement the I/O write counter
  !-----------------------------------------------------------------------------
  io%nwrite = io%nwrite - 1
  count = counters%decrement (this%task%iout+1, io%nwrite)
  !-----------------------------------------------------------------------------
  ! Remove the link from the task list, update the nbor lists of nbors, then
  ! deallocate this link (and the associated task memory and task)
  !-----------------------------------------------------------------------------
  call tlist%remove_link (this)
  call tlist%init_nbor_nbors (this)
  call tlist%check_nbor_nbors (this)
  call this%lock%unset ('delete')
  call this%delete (this)
  nullify (this)
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE remove_patch

!===============================================================================
!> Check the neighbors of the task associated with `link` to see if they can
!> provide support to newly refined patches made from `link`.
!> The current definition of "support" is a difference in levels of one or less.
!>
!> In the absence of the optional argument, `check_only`, neighbour patches will
!> be split as needed.
!===============================================================================
FUNCTION check_nbor_support (self, tlist, link, check_only) RESULT (refine_it)
  USE task_mod
  class(refine_t):: self
  class(list_t):: tlist
  class(link_t), pointer:: link
  integer:: refine_it
  !.............................................................................
  class(link_t), pointer:: nbor
  class(task_t), pointer:: this, npatch
  logical, optional:: check_only
  logical:: split_if_necessary
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%check_nbor_support')
  if (present(check_only)) then
    split_if_necessary = .not. check_only
  else
    split_if_necessary = .true.
  end if
  refine_it = -1
  this => link%task
  nbor => link%nbor
  do while (associated(nbor))
    npatch => nbor%task
    ! support is defined as having a level difference no greater than one
    if ((npatch%level - this%level) > 1) then
      if (verbose > 0) print &
        "('Patch',i5.1,' at level ',i2.1,' cannot provide support to patch',i5.1" &
        //",' at level',i2.1)", this%id, this%level, npatch%id, npatch%level
      if (split_if_necessary) then
        ! mark this patch for splitting to ensure support for neighbor patches
        if (verbose > 0) &
          print "('Splitting',i6.1,' to provide support')",this%id
        refine_it = 1
        exit
      end if
    end if
    nbor => nbor%next
  end do
  call trace%end()
END FUNCTION check_nbor_support

!===============================================================================
!> Check that the patch dimensions are the desired cubical ones
!===============================================================================
FUNCTION check_cubical (self, patch) RESULT (refine_it)
  class(refine_t):: self
  class(gpatch_t), pointer:: patch
  integer:: refine_it
  !.............................................................................
  if (all(patch%n==scene%patch_n)) then
    refine_it = -1
  else
    refine_it = 1
    patch%msplit = 1
  end if
END FUNCTION check_cubical

!===============================================================================
!> Calculate the Jeans length and compare it to the patch resolution.
!===============================================================================
FUNCTION check_Jeans (self, patch) RESULT(refine_it)
  USE math_mod
  USE scaling_mod
  class(refine_t):: self
  class(gpatch_t), pointer:: patch
  integer:: refine_it, imin(3), l(3), u(3)
  real:: Jeans_length, Jeans_num, dmin
  real, dimension(:,:,:), allocatable:: c_sound2, pg
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e
  !.............................................................................
  call trace%begin ('refine_t%check_check_jeans')
  l = patch%mesh%lb
  u = patch%mesh%ub
  d => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%it,1)
  allocate(c_sound2(size(d,1),size(d,2),size(d,3)), &
           pg(patch%gn(1),patch%gn(2),patch%gn(3)))
  select type (patch)
  class is (solver_t)
    pg = patch%gas_pressure()
  class default
    call mpi%abort("Unrecognised patch type. Abort!")
  end select
  dmin = minval(d)
  if (dmin <= 0.0) then
    imin = minloc(d)
    write (io%output,*) &
      'check_Jeans ERROR: id, dmin, imin =', patch%id, dmin, imin
    Jeans_num = min_jeans
  else
    c_sound2 = patch%gamma * pg(l(1):u(1),l(2):u(2),l(3):u(3)) / d
    if (associated(patch%mask)) then
      Jeans_length = minval(sqrt(math%pi * c_sound2 / (d*scaling%grav)), &
        mask=patch%mask)
    else
      Jeans_length = patch%fminval(sqrt(math%pi * c_sound2 / (d*scaling%grav)))
    end if
    Jeans_num = minval(patch%mesh%d) / Jeans_length
  end if
  refine_it = 0
  if (Jeans_num > max_jeans) refine_it = +1
  if (Jeans_num < min_jeans) refine_it = -1
  if (verbose > 2) &
    write (io_unit%output,'(a,i6,1p,3e9.2,i4)') &
      'Jeans criterion: id, max, num, min, refine', &
      patch%id, max_jeans, Jeans_num, min_jeans, refine_it
  deallocate (c_sound2, pg)
  nullify (d)
  call trace%end()
END FUNCTION check_Jeans

!===============================================================================
!> Refine on shocks, with threshold set by tolerance parameter `tol_shock`.
!===============================================================================
FUNCTION shock_and_cd_detector (self, patch) RESULT (refine)
  USE mesh_mod
  class(refine_t):: self
  class(gpatch_t), pointer:: patch
  integer:: refine
  !.............................................................................
  class(mesh_t), pointer:: mx, my, mz
  integer:: ix, iy, iz, ione=0, jone=0, kone=0
  real(kind=KindScalarVar):: dp1, dp2, dp3, dv1, dv2, dv3, dd1, dd2, dd3
  real, parameter:: small = 1.0e-10
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: pg
  real(kind=KindScalarVar), dimension(:,:,:,:), allocatable:: v
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e
  logical:: derefine(3)
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%shocks_and_cd_detector')

  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)

  if (mx%n > 1) ione = 1
  if (my%n > 1) jone = 1
  if (mz%n > 1) kone = 1

  allocate(pg(patch%mesh(1)%gn,patch%mesh(2)%gn,patch%mesh(3)%gn), &
           v(patch%mesh(1)%gn,patch%mesh(2)%gn,patch%mesh(3)%gn,3))
  ! we need the gas pressure, density and velocity
  select type (patch)
  class is (solver_t)
    d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
    e => patch%mem(:,:,:,patch%idx%e,patch%it,1)
    pg = patch%gas_pressure()
    v  = patch%gas_velocity()
  class default
    call mpi%abort("Unrecognised patch type. Abort!")
  end select
  !
  refine = 0
  do iz=mz%li,mz%ui
    do iy=my%li,my%ui
      do ix=mx%li,mx%ui
        dp1 = abs(pg(ix+ione,iy,iz) - pg(ix-ione,iy,iz)) &
            / (min(pg(ix+ione,iy,iz), pg(ix-ione,iy,iz)) + small )
        dp2 = abs(pg(ix,iy+jone,iz) - pg(ix,iy-jone,iz)) &
            / (min(pg(ix,iy+jone,iz), pg(ix,iy-jone,iz)) + small )
        dp3 = abs(pg(ix,iy,iz+kone) - pg(ix,iy,iz-kone)) &
            / (min(pg(ix,iy,iz+kone), pg(ix,iy,iz-kone)) + small )
        !
        if (max_shock > 0.0) then
          dv1 = v(ix  ,iy,iz,1) - v(ix+ione,iy,iz,1)
          dv2 = v(ix,iy  ,iz,2) - v(ix,iy+jone,iz,2)
          dv3 = v(ix,iy,iz  ,3) - v(ix,iy,iz+kone,3)
          !
          if ( dp1 > max_shock .and. dv1 > small ) refine = 1
          if ( dp2 > max_shock .and. dv2 > small ) refine = 1
          if ( dp3 > max_shock .and. dv3 > small ) refine = 1
          !
          if ( dp1 < min_shock .and. dv1 > small ) derefine(1) = .true.
          if ( dp2 < min_shock .and. dv2 > small ) derefine(2) = .true.
          if ( dp3 < min_shock .and. dv3 > small ) derefine(3) = .true.
          if (all(derefine)) refine = -1
        end if
        !
        if (max_contact > 0.0) then
          dd1 = abs(d(ix+ione,iy,iz) - d(ix-ione,iy,iz)) &
              / (min(d(ix+ione,iy,iz), d(ix-ione,iy,iz)) + small )
          dd2 = abs(d(ix,iy+jone,iz) - d(ix,iy-jone,iz)) &
              / (min(d(ix,iy+jone,iz), d(ix,iy-jone,iz)) + small )
          dd3 = abs(d(ix,iy,iz+kone) - d(ix,iy,iz-kone)) &
              / (min(d(ix,iy,iz+kone), d(ix,iy,iz-kone)) + small )
          !
          if ( dd1 > max_contact .and. dp1 < 0.01*max_contact ) refine = 1
          if ( dd2 > max_contact .and. dp2 < 0.01*max_contact ) refine = 1
          if ( dd3 > max_contact .and. dp3 < 0.01*max_contact ) refine = 1
          !
          if ( dd1 < min_contact .and. dp1 < 0.01*max_contact ) derefine(1) = .true.
          if ( dd2 < min_contact .and. dp2 < 0.01*max_contact ) derefine(2) = .true.
          if ( dd3 < min_contact .and. dp3 < 0.01*max_contact ) derefine(3) = .true.
          if (all(derefine)) refine = -1
        end if
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Cancel actions that would go outside the interval [levelmin...levelmax]
  !-----------------------------------------------------------------------------
  if (patch%level == self%levelmax) refine = min(refine,0)
  if (patch%level == self%levelmin) refine = max(refine,0)
  if (verbose > 1) then
    if (refine > 0) then
      write (io%output,*) 'shock_and_cd_detector: refine', patch%id, patch%level, self%levelmax
    else if (refine < 0) then
      write (io%output,*) 'shock_and_cd_detector: derefine', patch%id, patch%level, self%levelmin
    end if
  end if
  !
  nullify (mx, my, mz, d)
  deallocate (pg, v)
  !
  call trace%end()
END FUNCTION shock_and_cd_detector

!===============================================================================
!> Refine on shocks, with threshold set by tolerance parameter `max_shock`.
!===============================================================================
FUNCTION gradient_detector (self, patch) RESULT (refine)
  USE mesh_mod
  USE solver_mod
  class(refine_t):: self
  class(gpatch_t), pointer:: patch
  integer:: refine
  !.............................................................................
  class(mesh_t), pointer:: mx, my, mz
  integer:: i, j, k, ione=0, jone=0, kone=0, im1, jm1, km1, ip1, jp1, kp1, n
  real, parameter:: small=1e-10
  real(kind=KindScalarVar):: dqmax
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: q, dq1, dq2, dq3
  real, parameter:: rsmall = 1.0e-10
  logical:: is_unsigned
  character(len=1):: unsigned(3) = ['d', 'p', 'e']
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%gradient_detector')
  allocate (dq1(patch%gn(1),patch%gn(2),patch%gn(3)), &
            dq2(patch%gn(1),patch%gn(2),patch%gn(3)), &
            dq3(patch%gn(1),patch%gn(2),patch%gn(3)), &
            q  (patch%gn(1),patch%gn(2),patch%gn(3)))
  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)
  if (mx%n > 1) ione = 1
  if (my%n > 1) jone = 1
  if (mz%n > 1) kone = 1
  refine = 0
  do n=1,size(grad_var)
    if (trim(grad_var(n)) == '') exit
    if (trim(grad_var(n)) == 'd')  q = patch%mem(:,:,:,patch%idx%d,patch%it,1)
    if (trim(grad_var(n)) == 'e')  q = patch%mem(:,:,:,patch%idx%e,patch%it,1)
    if (trim(grad_var(n)) == 'p1') q = patch%mem(:,:,:,patch%idx%px,patch%it,1)
    if (trim(grad_var(n)) == 'p2') q = patch%mem(:,:,:,patch%idx%py,patch%it,1)
    if (trim(grad_var(n)) == 'p3') q = patch%mem(:,:,:,patch%idx%pz,patch%it,1)
    if (trim(grad_var(n)) == 'b1') q = patch%mem(:,:,:,patch%idx%bx,patch%it,1)
    if (trim(grad_var(n)) == 'b2') q = patch%mem(:,:,:,patch%idx%by,patch%it,1)
    if (trim(grad_var(n)) == 'b3') q = patch%mem(:,:,:,patch%idx%bz,patch%it,1)
    select type(patch)
    class is (solver_t)
      if (trim(grad_var(n)) == 'p')  q = patch%gas_pressure()
      if (trim(grad_var(n)) == 'v1') q = patch%gas_velocity(1)
      if (trim(grad_var(n)) == 'v2') q = patch%gas_velocity(2)
      if (trim(grad_var(n)) == 'v3') q = patch%gas_velocity(3)
    class default
      call mpi%abort("Unrecognised patch type. Abort!")
    end select
    is_unsigned = .false.
    do i=1,size(unsigned)
      if (grad_var(n) == unsigned(i)) then
        is_unsigned = .true.
        exit
      end if
    end do
    ! calculate the gradient in each direction
    do k=mz%lo,mz%uo
      km1 = k - kone
      kp1 = k + kone
      do j=my%lo,my%uo
        jm1 = j - jone
        jp1 = j + jone
        if (is_unsigned) then
          do i=mx%lo,mx%uo
            im1 = i - ione
            ip1 = i + ione
            dq1(i,j,k) = abs ( q(ip1,j,k) - q(im1,j,k) ) &
                             / ( max ( q(ip1,j,k), q(im1,j,k) ) &
                               + small )
            dq2(i,j,k) = abs ( q(i,jp1,k) - q(i,jm1,k) ) &
                             / ( max ( q(i,jp1,k), q(i,jm1,k) ) &
                               + small )
            dq3(i,j,k) = abs ( q(i,j,kp1) - q(i,j,km1) ) &
                             / ( max ( q(i,j,kp1), q(i,j,km1) ) &
                               + small )
          end do
        else
          do i=mx%lo,mx%uo
            im1 = i - ione
            ip1 = i + ione
            dq1(i,j,k) = abs ( q(ip1,j,k) - q(im1,j,k) )
            dq2(i,j,k) = abs ( q(i,jp1,k) - q(i,jm1,k) )
            dq3(i,j,k) = abs ( q(i,j,kp1) - q(i,j,km1) )
          end do
        end if
      end do
    end do
    ! choose max. gradient; if the max. is greater than the threshold, flag
    ! for refinement
    do k=mz%li+kone,mz%ui-kone
      km1 = k - kone
      kp1 = k + kone
      do j=my%li+jone,my%ui-jone
        jm1 = j - jone
        jp1 = j + jone
        do i=mx%li+ione,mx%ui-ione
          im1 = i - ione
          ip1 = i + ione
          dqmax = max ( dq1(i,j,k), dq1(ip1,j,k), dq1(im1,j,k) &
                      , dq2(i,j,k), dq2(i,jp1,k), dq2(i,jm1,k) &
                      , dq3(i,j,k), dq3(i,j,kp1), dq3(i,j,km1) )
          if (patch%level < self%levelmax .and. dqmax > max_grad(n)) refine = 1
          if (patch%level > self%levelmin .and. dqmax < min_grad(n)) refine = -1
        end do
      end do
    end do
    if (verbose > 1) then
      if (refine > 0) then
        write (io%output,*) 'gradient_detector: refine ', grad_var(n), patch%id
      else if (refine < 0) then
        write (io%output,*) 'gradient_detector: derefine ', grad_var(n), patch%id
      end if
    end if
  end do
  deallocate (dq1, dq2, dq3, q)
  call trace%end()
END FUNCTION gradient_detector

!===============================================================================
!> Given a box defined by its lower and upper corners, flag a patch for refinement
!> if it lies inside this box, up to a maximum level.
!> This particular criteria is primarily used for testing and debugging.
!===============================================================================
FUNCTION refine_region (self, patch) RESULT (refine)
  USE mesh_mod
  class(refine_t) :: self
  class(gpatch_t), pointer :: patch
  integer:: refine
  class(mesh_t), pointer:: m1, m2, m3, m(:)
  logical :: intersect
  logical :: refine_it
  integer :: ix, iy, iz
  real :: x, y, z
  real, dimension(3) :: p_llc, p_urc, b_llc, b_urc
  !.............................................................................
  call trace%begin('refine_t%refine_on_box')

  refine_it = .false.
  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)
  m => patch%mesh(:)

  ! check if the patch (as a whole) intersects the box-shaped region of interest
  do ix=1,3
    associate (mr => patch%mesh(ix)%r)
    p_llc(ix) = m(ix)%p + mr(m(ix)%li) - 0.5 * m(ix)%d
    p_urc(ix) = m(ix)%p + mr(m(ix)%ui) + 0.5 * m(ix)%d
    end associate
  end do
  b_llc = region_llc
  b_urc = region_urc

  intersect = all((p_urc - b_llc) > 0) .and. all((b_urc - p_llc) > 0)
  if (.not. intersect) then
    call trace_end
    return
  endif

  do iz=m3%li,m3%ui
    z = m3%p + m3%r(iz)
    do iy=m2%li,m2%ui
      y = m2%p + m2%r(iy)
      do ix=m1%li,m1%ui
        x = m1%p + m1%r(ix)
        refine_it = refine_it .or. (x >= b_llc(1) .and. x <= b_urc(1) .and. &
                                    y >= b_llc(2) .and. y <= b_urc(2) .and. &
                                    z >= b_llc(3) .and. z <= b_urc(3) )
      end do
    end do
  end do
  refine = merge(1,0,refine_it)
  call trace%end()
END FUNCTION refine_region

END MODULE refine_mod
