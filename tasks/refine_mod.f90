!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> This module handles checking max change between neighboring points.  Each
!> instance of it needs an index to the variable, whether to take the log, and
!> a value for the max allowed change.  It might be applied to log density, log
!> pressure, log opacity, velocity amplitude, magnetic field amplitude, etc.
!>
!> task_list_t%update
!>   check_current
!>     refine_needed
!>     selective_refine
!>       if ... refine
!>         make_child_patch
!>           tasklist%lock%set
!>           init_nbors ...
!>           check_nbors ...
!>           tasklist%lock%unset
!>       if ... derefine
!>         remove_and_reset
!>           tasklist%lock%set
!>           init_nbors ...
!>           check_nbors ...
!>           tasklist%lock%unset
!>   download_link
!>   update
!>   ...
!>
!> Support criterion
!> -----------------
!> A patch should be able to load guard zone values from nbors differing at
!> most by one AMR-level, to avoid too abrupt changes of resolution at patch
!> boundaries.  The criterion to test must be formulated as a criterion on
!> the patch that should potentially be created, and not as a criterion on
!> the that needs the support (to avoid having to find a sub-region of the 
!> parent of the parent -- and thus have to maintain parent relationships)
!>
!> On the other hand, the criterion needs to be checked at the point in time
!> when a new, higher level patch is to be created, to check if this requires
!> that an nbor patch should also be refined.  The typical situation to look
!> for is when creating a patch at level L+1 requires of the patch whose sub-
!> region is being investigated requires a new patch at level L -- a so far
!> non-existing nbor of the patch at level L under investigation.  The logic
!> is, in brief:
!>
!>   for candidate in sub_regions_of (patch):
!>     if needed(candidate):
!>       attach_preliminary_nbor_list_to (candidate)
!>       for nbor in candidate.nbors:
!>         if nbor.patch.level==patch_candidate+2:
!>           for nbor_candidate in sub_region_of (nbor.patch):
!>             if overlap_between (candidate,nbor_candidate):
!>               create(nbor_candidate)
!>               candidate.add_nbor (nbor_candidate)
!>               update_nbor_lists
!>
!> Note that this procedure does not assume a specific (e.g. side-by-side)
!> arrangement of patches, and hence allows for moving patches.  Thus, it also
!> does not rely on the existence of a parent-child relation between patches.
!>
!> OpenMP synchronization is implemented with two (nested) OMP locks; one
!> (check_lock) responsible for synchronizing calls to list_t%check_nbors,
!> and one (download_lock) responsible for synchronizing filling of guard
!> zones (as well as restrict and prolong operations).
!>
!> When this module decides to refine or derefine a task, it acquires both
!> locks, and hence blocks check_nbors and download_link() calls while the
!> operation is on-going.  Since no time advance is involved in the refine
!> / derefine, the blocking cannot cause a hang.
!>
!> Since refine/derefine is not done every timestep, the impact on running
!> tasks from the locking is not severe, and since two different locks are
!> used, one avoids the significant impact that would be the result of
!> locking out check_nbors and download_link from occuring at the same time.
!>
!> The impact on download_link could be entirely removed by always using
!> the sorted nbor list, and by collecting patches to be removed on a GC
!> list, removing them a few AMR cycles later.
!>
!> The absolute minimum amount of synchronization would be to just make sure
!> that an nbor list 1) can never change exactly when it is accessed, 2) to
!> immediately make a copy of it, 3) to make sure no to remove any task on the
!> nbor list while it is being used, and 4) to make sure that rotating the
!> memory slots will not invalidate using the nbor list, as is.
!>
!> An nbor list relevant for a particular task update is accessed in the
!> following order:
!>
!> 1) During a link_t%check_nbors() call, made after an nbor of the task in
!>    question has been updated, link_t%check_ready() is called, and runs
!>    through the tasks nbor list, comparing task%time values, finding that
!>    all nbors of the task are advanced enough to provide boundary zone data,
!>    or other data that the task needs.  The task is then added to the ready
!>    queue.  The comparisons of times requires only that task_t%time is updated
!>    atomically. These updates happen in task_t%rotate, which also manipulates
!>    the task_t%iit(:) array, which holds the indices to the memory slots.
!>
!> 2) After a thread picks that task off the head of the ready queue, the first
!>    action (in task_list_t%update) is to check if the task should be refined
!>    or derefined.  The test itself does not use the nbor lists, but if either
!>    the current thread, or a thread working on an nbor decides to refine or
!>    derefine, then the nbor list of the task may become modified.
!>
!>    a) If it is refined, the data upon which the decision to use it as a
!>    source of guard zone data still exists (in the parent task), and the
!>    refined data in the child task(s) is initially just a prolongation of the
!>    parent data, so does not yet provide any improvement.  Therefore, there
!>    is no problem, and the insertion of the new nbor list should be delayed
!>    until it is entirely complete, at which point the insertion can be done
!>    using a very brief locking of the link, while changeing task%link%nbor to
!>    point to the new first nbor.
!>
!>    b) If it is derefined, one just has to leave the nbor list and the task
!>    there until after they have been used, in guard zone download, or in
!>    check_ready().
!>
!> 3) The next step is download of data, which relies on data from all source
!>    tasks, one of which could possible be in the process of being derefined.
!>    To minimize the acess time needed, the download procedure should copy the
!>    nbor list, a step that is anyway needed to sort the list into level order.
!>    Then, to guarantee existence of the derefined patch, it has to either be
!>    prevented from going away by locking, or else it must be put in a garbage
!>    collector, which only removes patches when they are no longer needed.
!>
!>    The latter option is the easier one, since it avoid  messy synchronization.
!>    A sync would require that derefinement is prevented in the interval btw
!>    the task has been put in the ready queue, and until all nbors that might
!>    possibly need it has done their downloads; a tall order.
!>
!> As part of putting the task in the ready queue, a copy of the nbor list
!> should be made -- this could be exactly the sorted copy that will be needed
!> in the download step in any case.
!>
!> Keeping track of "no longer needed" can be done with an atomically updated
!> counter in the task, which starts out with a value of one, and is incremented
!> when a task is put on a sorted, temporary nbor list, pending use in download.
!> When used as a source in download, the counter is decremented, endingg up on
!> one again, unless the task has been derefined away in the meantime, in which
!> case the counter ends on zero, and the task is then deallocated / deleted.
!>
!> Updating an nbor list with list_t%init_nbors() should be done in two steps:
!> first the head of the list is allocated, without linking it to link%nbor,
!> and the rest of the tasks are appended (or prepended).  Then, in an lock
!> protected assignment, link%nbor is pointed to the head of the new nbor list.
!> Before that, a link is saved to the head of the existing nbor list, and after
!> the atomic reassignment, the old list is deleted.
!>
!> A derefined task is immediately removed from the task list, and hence can no
!> longer be added to re-initialized nbor lists.  It thus ceases to be needed
!> when the last task that still has it on its nbor list copy deallocates its
!> nbor list, after using it in dowload_link().
!>
!> Summary:  A sorted copy of the nbor list is created in connection with the
!> queue_by_time() call in check_ready().  At the same time a counter task_t%
!> n_needed is incremented for all tasks in the temporary nbor list.  When the
!> task is used as a source in download_link() the counter is decremented
!> (atomically), and if the task is derefined the counter is decremented there
!> as well, leading to a zero cound when the task is no longer on an nbor list
!> copy.  Then, and only then, the task may be deleted.
!>
!> FIXME:  Care should be extended also to updates over MPI.
!> FIXME:  Use of the deprecated shared_mod should be removed
!>
!> MPI specific handling:  Boundary patches are always sent to vnbors, and if
!> they are new on the recieving rank, new nbor lists are generated, so no new
!> actions should be needed there.   If a boundary task is removed, the vnbors
!> are send a copy with the "remove" bit set, and they should react correspond-
!> ingly, by doing essentially what remove_and_reset does.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE refine_mod
  USE iso_fortran_env, only: int8
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
  USE patch_mod
  USE solver_mod
  USE extras_mod
  USE shared_mod
  USE download_mod
  USE timer_mod
  USE data_io_mod
  implicit none
  private
  !
  type:: refine_t
    logical:: on                              ! allow external access
    integer:: levelmin=0                      ! Min level of refinment
    integer:: levelmax=0                      ! Max level of refinment
    integer:: ratio=2                         ! default
    real:: safe_ratio=2.1                     ! safe derefine ratio
  contains
    procedure:: init
    procedure:: refine_needed
    procedure:: check_current
    procedure:: check_support
    procedure:: need_to_support
    procedure:: check_cubical
    procedure:: check_jeans
    procedure:: shock_and_cd_detector
    procedure:: refine_region
    procedure:: gradient_detector
    procedure:: sanity_check
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
  real:: min_compress=0.0, max_compress=0.0   ! Tolerance values for compression (negative divergence)
  real:: min_vorticity=0.0, max_vorticity=0.0 ! Tolerance values for vorticity
  integer, parameter:: ngradvar=10            ! Number of gradient variables
  real, dimension(3):: region_llc=0.0         ! Sometimes you want to refine a region (manually)
  real, dimension(3):: region_urc=0.0
  real, dimension(ngradvar):: max_grad=0.0    ! Tolerance values for gradient detector refinement
  real, dimension(ngradvar):: min_grad=0.0    ! Tolerance values for gradient detector derefinement
  character(len=32), dimension(ngradvar):: grad_var=""
  logical:: on=.false.                        ! Flag to turn the entire process on or off
  logical:: force_cubical=.false.             ! Force cubical patches
  logical:: detailed_timer=.false.            ! Lump all AMR together when false
  integer:: n_locks=1                         ! number of locks to use
  integer:: verbose=0                         ! Noise level
  type(refine_t), public:: refine
CONTAINS

!===============================================================================
!> Initialize refinement parameters (called from task_list_t%init_levels)
!===============================================================================
SUBROUTINE init (self)
  class(refine_t):: self
  integer:: iostat
  integer, save:: levelmin=0, levelmax=0
  integer, save:: ratio=2                     ! How does a patch split when it is refined
  logical, save:: first_time = .true.
  namelist /refine_params/ verbose, on, levelmin, levelmax, check_interval, &
    ratio, refine_from_t, refine_until_t, min_dx, max_dx, max_jeans, min_jeans, &
    max_shock, min_shock, max_contact, min_contact, max_grad, min_grad, &
    max_vorticity, min_vorticity, max_compress, min_compress, &
    grad_var, region_llc, region_urc, force_cubical, n_locks, detailed_timer
  character(len=120):: id = &
    '$Id$ tasks/refine_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%init')
  call trace%print_id (id)
  !$omp critical (input_cr)
  if (first_time) then
    first_time=.false.
    levelmin = shared%levelmin
    levelmax = shared%levelmax
    rewind (io%input)
    read (io%input, refine_params, iostat=iostat)
    write (io%output, refine_params)
    if (on) then
      omp_lock%links = .true.
    else
      omp_lock%links = .false.
    end if
  end if
  !-----------------------------------------------------------------------------
  ! If AMR is active, we must use check_filled, since nbor levels are in
  ! decreasing order, to save computing time
  !-----------------------------------------------------------------------------
  !$omp end critical (input_cr)
  self%levelmin = levelmin
  self%levelmax = levelmax
  shared%levelmin = levelmin
  shared%levelmax = levelmax
  self%ratio = ratio
  self%on = on
  !-----------------------------------------------------------------------------
  ! Initialize the lock in the garbage collector (shared by link_mod)
  !-----------------------------------------------------------------------------
  call garbage%init
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
  class(solver_t), pointer:: patch
  class(extras_t), pointer:: extras
  real:: dxmin, dxmax
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  refine = 0
  if (.not. on) return
  call trace%begin('refine_t%needed', itimer=itimer, detailed_timer=detailed_timer)
  !----------------------------------------------------------------------------
  ! Refinement criteria, starting with assumed derefinement
  !------------------------------------------------------------------------
  patch => solver%cast2solver (link%task)
  refine = -1
  !
  ! Extras refinement
  extras => patch
  refine = max(refine, extras%check_refine (extras))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'extras'
  !
  ! Ensure consistent, cubical patch size (e.g. after RAMSES startup)
  if (force_cubical) refine = max (refine, self%check_cubical(patch))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'cubical'
  !
  ! Jeans criterion
  if (max_jeans > 0.0) refine = max (refine, self%check_Jeans(patch))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'Jeans'
  !
  ! shock and/or contact discontinuity detector
  if (max_shock > 0.0 .or. max_contact > 0.0) &
    refine = max (refine, self%shock_and_cd_detector(patch))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'shock'
  !
  ! compression detector
  if (max_compress > 0.0) &
    refine = max (refine, check_compress (self, patch))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'compress'
  !
  ! vorticity detector
  if (max_vorticity > 0.0) &
    refine = max (refine, check_vorticity (self, patch))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'vorticity'
  !
  ! gradient detector
  if (any(max_grad > 0.0)) refine = max (refine, self%gradient_detector(patch))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'gradient'
  !
  ! refine on (a) user-specified region(s)
  if (any(region_llc /= 0.0) .and. any(region_urc /= 0.0)) &
    refine = max (refine, self%refine_region(patch))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'user'
  !
  ! Check if support constraints require refinement
  refine = max (refine, need_to_support (self, tasklist, link))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'need to support'
  !
  ! Check if support constraints allow derefinement
  refine = max (refine, check_derefine_support (self, tasklist, link))
  if (verbose > 3) write (io_unit%output,*) &
    patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'derefine support'
  !
  call trace%end (itimer, detailed_timer=detailed_timer)
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
  class(link_t), pointer:: link, nbor
  class(solver_t), pointer:: patch
  logical, intent(out):: was_refined, was_derefined
  real:: dxmin, dxmax
  integer:: refine, id, n_added, l(3), u(3)
  integer, save:: itimer=0
  !.............................................................................
  was_refined = .false.
  was_derefined = .false.
  !-----------------------------------------------------------------------------
  ! Punt on refinement if the task is very new, or has the init_nbors bit set
  !-----------------------------------------------------------------------------
  if (.not.on .or. &
      link%task%istep < 3 .or. &
      link%task%is_set (bits%init_nbors)) &
        return
  !---------------------------------------------------------------------------
  call trace%begin('refine_t%check_current', itimer=itimer)
  patch => solver%cast2solver (link%task)
  id = patch%id
  if (verbose > 3) &
    write (io_unit%output,'(i6,2x,a,f12.6,2i4)') &
      link%task%id, 'refine_t%check_current, at t =', &
      link%task%time, link%task%level, self%levelmax
  !----------------------------------------------------------------------------
  dxmin = minval(patch%ds)
  dxmax = maxval(patch%ds)
  if (verbose > 3) then
    write (io_unit%output,'(a,i6,2i4,l4,1p,4e11.2)') &
      'refine_t%refine_needed: id, istep, refine_next =', &
      patch%id, patch%istep, patch%check_refine_next, &
      patch%is_set(bits%virtual), min_dx, dxmin, max_dx, dxmax
  end if
  !----------------------------------------------------------------------------
  ! Filter out patches that do not match the necessary criteria
  !----------------------------------------------------------------------------
  if (patch%is_clear(bits%virtual)  .and. &
      dxmin >  min_dx               .and. &
      dxmax <= max_dx               .and. &
      patch%time >= refine_from_t   .and. &
      patch%time <  refine_until_t  .and. &
      patch%istep >= patch%check_refine_next) then
    if (verbose > 3) &
      tasklist%verbose = verbose
    !---------------------------------------------------------------------------
    ! Before testing any criteria, check if patch%irefine is set
    !---------------------------------------------------------------------------
    if (.not.allocated(patch%irefine)) &
      allocate (patch%irefine(patch%gn(1),patch%gn(2),patch%gn(3)))
    patch%irefine = -1
    l = patch%mesh%li
    u = patch%mesh%ui
    refine = maxval(patch%irefine(l(1):u(1),l(2):u(2),l(3):u(3)))
    if (verbose > 1 .and. &
        any(patch%irefine(l(1):u(1),l(2):u(2),l(3):u(3)) > 0)) then
      write (stdout,'(3x,a,i6,":",i2,2x,a)')  &
        'refine'    ,  patch%id, patch%level, &
        'requested for support'
    end if
    !---------------------------------------------------------------------------
    ! Check all refinement criteria.  If patch%irefine is already marked for
    ! refinement, then it has been done to request support of new level+2 patches
    !---------------------------------------------------------------------------
    refine = max(refine,self%refine_needed (tasklist, link))
    !---------------------------------------------------------------------------
    if (verbose > 3) write (io_unit%output,*) &
      patch%id, refine, minval(patch%irefine), maxval(patch%irefine), 'mask'
    !---------------------------------------------------------------------------
    ! Apply refinement criteria on potential child patch parts, but only if
    ! either refinement has been indicated, or if derfinement is allowed
    !---------------------------------------------------------------------------
    if (refine /= 0) then
      call selective_refine (self, tasklist, patch, refine, was_refined, was_derefined)
      if (was_derefined) then
        call trace%end (itimer)
        return
      end if
    end if
    !---------------------------------------------------------------------------
    ! Update next refinement step counter and deallocate local irefine
    !---------------------------------------------------------------------------
    patch%check_refine_next = patch%istep + check_interval
    deallocate (patch%irefine)
  end if
  call trace%end (itimer)
END SUBROUTINE check_current

!===============================================================================
!> Loop over the ratio x ratio x ratio positions and determine whether to refine
!> or derefine.
!===============================================================================
SUBROUTINE selective_refine (self, tasklist, patch, refine, was_refined, was_derefined)
  class(refine_t):: self
  class(list_t):: tasklist
  class(solver_t), pointer:: patch
  integer:: refine
  logical:: was_refined, was_derefined
  !.............................................................................
  class(solver_t), pointer:: child
  class(link_t), pointer:: nbor, link, nbors, tail
  class(patch_t), pointer:: nbpatch
  integer:: i, j, k
  integer, dimension(3):: l, u, n, ii
  real(8), dimension(3):: size, pos
  logical:: not_already_refined
  integer:: refine_it, n_added
  integer(int8), allocatable:: refined(:,:,:)
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%selective_refine', itimer=itimer, detailed_timer=detailed_timer)
  was_refined = .false.
  was_derefined = .false.
  link => patch%link
  call link%lock%set ('selective_refine')
  nbors => link%nbor
  tail => tasklist%tail
  nbor => patch%link%nbor
  !-----------------------------------------------------------------------------
  ! As preparation, mark all points that already are refined
  !-----------------------------------------------------------------------------
  allocate (refined(patch%gn(1),patch%gn(2),patch%gn(3)))
  refined = 0
  do while (associated(nbor))
    nbpatch => solver%cast2solver(nbor%task)
    if (patch%contains(nbpatch%position) .and. &
        patch%level == nbpatch%level-1) then
      l = patch%index_only (nbpatch%position + nbpatch%mesh%lf)
      u = patch%index_only (nbpatch%position + nbpatch%mesh%uf)
      !-------------------------------------------------------------------------
      ! A patch with existing child patches should not be derefined.  To prevent
      ! this, refined regions inside are set to at least 0 in irefine, and to
      ! prevent creating doubles, the refined array is set to 1
      !-------------------------------------------------------------------------
      refined(l(1):u(1),l(2):u(2),l(3):u(3)) = 1
      patch%irefine(l(1):u(1),l(2):u(2),l(3):u(3)) = &
        max(0,patch%irefine(l(1):u(1),l(2):u(2),l(3):u(3)))
      !-------------------------------------------------------------------------
      if (verbose > 2) then
        write (io_unit%output,'(a,i6,a,i6,2(3x,3i4))') &
          'derefine: patch', patch%id, ' contains', nbpatch%id, l, u
      end if
      !-------------------------------------------------------------------------
    end if
    nbor => nbor%next
  end do
  call link%lock%unset ('selective_refine')
  !-----------------------------------------------------------------------------
  ! First, consider if any part should be refined
  !-----------------------------------------------------------------------------
  if (patch%level < self%levelmax .and. refine > 0) then
    if (any(patch%irefine > 0)) then
      do k=1,self%ratio
        do j=1,self%ratio
          do i=1,self%ratio
            !-------------------------------------------------------------------
            ! Size and position
            !-------------------------------------------------------------------
            size = patch%size/self%ratio
            pos = patch%position + ([i,j,k]-(self%ratio+1)*0.5)*size
            !-------------------------------------------------------------------
            ! Lower and upper index in parent patch
            !-------------------------------------------------------------------
            n = patch%mesh%n/self%ratio
            l = patch%mesh%li + ([i,j,k]-1)*n
            u = l + n - 1
            !-------------------------------------------------------------------
            ! Should this part be refined or not?
            !-------------------------------------------------------------------
            refine_it = maxval (patch%irefine(l(1):u(1),l(2):u(2),l(3):u(3)))
            !-------------------------------------------------------------------
            ! Is this part refined before, or should it be now?
            !-------------------------------------------------------------------
            ii = (l+u)/2
            not_already_refined = refined(ii(1),ii(2),ii(3)) == 0
            if (not_already_refined .and. refine_it > 0) then
              !io%do_trace = .true.
              !-----------------------------------------------------------------
              ! Create a refined patch here
              !-----------------------------------------------------------------
              if (verbose > 1) &
                write (io_unit%output,'(a,i6,1p2g16.5,2x,3(i4,":",i2))') 'refine', &
                  patch%id, patch%fmaxval(patch%mem(:,:,:,patch%idx%d,patch%it,1)), &
                  maxval(patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%it,1)), &
                  l(1),u(1),l(2),u(2),l(3),u(3)
              call make_child_patch (self, tasklist, patch, pos, size, child)
              !-----------------------------------------------------------------
              if (verbose > 0) then
                associate (unit => merge (stdout, io_unit%output, verbose>1))
                write (unit,'(2(3x,a,i6,":",i2),2x,a,2l1)')  &
                  'refine'    ,  patch%id, patch%level, &
                  'to'        ,  child%id, child%level, &
                  'BV:', child%is_set(bits%boundary), child%is_set(bits%virtual)
                flush (unit)
                end associate
              end if
              !----------------------------------------------------------------- 
              was_refined= .true.
            end if
          end do
        end do
      end do
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Then, consider if the whole patch should be derefined.  This requires that
  ! it has no higher level patch inside (would trigger patch%refine > 0 there),
  ! and secondly that irefine(:,:,:) indicates that derefine would be OK
  !-----------------------------------------------------------------------------
  if (refine == -1                .and. &               ! derefine possible
      patch%level > self%levelmin .and. &               ! relevant?
      patch%check_refine_next > 0 .and. &               ! not immediately!
      patch%istep >= 4            .and. &               ! must take 4 steps
      .not.was_refined) then                            ! not after a refine!
    l = patch%mesh%li
    u = patch%mesh%ui
    if (all(patch%irefine(l(1):u(1),l(2):u(2),l(3):u(3)) == -1)) then
      !-----------------------------------------------------------------------
      if (verbose > 0) then
        associate (unit => merge (stdout, io_unit%output, verbose > 1))
        write (io_unit%output,'(1x,a,i6,":",i2,2x,"B:",l1,3x,a,i6,3x,a,f12.6)') &
          'derefine', patch%id, patch%level, patch%is_set (bits%boundary), &
          'istep:', patch%istep, 'time:', patch%time
        flush (unit)
        end associate
      end if
      !-------------------------------------------------------------------------
      ! If task is to be removed, call the comprehensive remove procedure, which
      ! includes calling set_init_nbors(), queueing all nbors for immediate
      ! update of their nbor lists.
      !-------------------------------------------------------------------------
      call data_io%update_counters (patch, -1)
      call patch%log ('remove')
      call tasklist%remove_and_reset (link)
      was_derefined = .true.
    end if
  end if
  deallocate (refined)
  call trace%end (itimer, detailed_timer=detailed_timer)
END SUBROUTINE selective_refine

!===============================================================================
!> Check if new%task is a patch, and if it overlaps with other level L-2 patches,
!> with points not covered by L-1 patches, indicating that a new L-1 patch should
!> be created from the L-2 patch.
!>
!> To implement this, we need to look at the regions that overlap with an L-2
!> patch, and then run through the nbor list L-1 patches, to see if they cover
!> the area.  If not, "make an impression" on the L-2 patch, by setting the
!> corresponding cells in its %xrefine map.
!>
!> NOTE: This mechanism is currently NOT used; it is left for now, to be removed
!> when need_to_support() has been validated as being sufficient -- it looks
!> at the same issue from the point of view of the L-2 patch, which is better,
!> since it eliminates the need for synchronizing access to %xrefine.
!===============================================================================
SUBROUTINE check_support (self, tlist, new, n_added)
  class(refine_t):: self
  class(list_t):: tlist
  class(link_t), pointer:: new
  integer:: n_added
  !.............................................................................
  class(link_t), pointer:: old, nbor
  class(patch_t), pointer:: old_patch, new_patch, nbor_patch
  type(patch_t):: patch
  integer:: l(3), u(3)
  logical:: overlap, guards
  integer, allocatable, dimension(:,:,:):: trefine
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  ! Return if there are less than two levels down to levelmin
  !-----------------------------------------------------------------------------
  n_added = 0
  if (new%task%level <= self%levelmin+1) return
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_support', itimer=itimer, detailed_timer=detailed_timer)
!write (io_unit%log,*) 'unexpected call to check_support for', new%task%id, &
!new%task%is_set(bits%boundary),  new%task%is_set(bits%virtual)
  !-----------------------------------------------------------------------------
  ! Select patch-based tasks, and assume for now it is OK for updates
  !-----------------------------------------------------------------------------
  new_patch => patch%cast2patch (new%task)
  if (associated(new_patch)) then
    if (verbose > 2) &
      write (stdout,*) 'check_support: id =', &
      new_patch%id, new_patch%is_set (bits%support)
    !---------------------------------------------------------------------------
    ! Search task list for patch-based overlapping L-2 patches
    !---------------------------------------------------------------------------
    call new_patch%clear (bits%support)
    call tlist%lock%set ('check_support')
    old => tlist%head
    do while (associated(old))
      if (old%task%level <= new%task%level-2) then
        old_patch => patch%cast2patch (old%task)
        guards = .true.
        if (old_patch%get_overlap (new_patch, guards, l, u)) then
          if (verbose > 2) &
            write (stdout,'(a,i6,3(2x,i2,":",i2))') &
            '  level-2 overlap: id =', old_patch%id, &
            l(1),u(1), l(2),u(2), l(3),u(3)
          !---------------------------------------------------------------------
          ! Allocate a temporary integer mask, to assemble evidence in
          !---------------------------------------------------------------------
          allocate (trefine(old_patch%gn(1),old_patch%gn(2),old_patch%gn(3)))
          trefine = 0
          trefine (l(1):u(1), l(2):u(2), l(3):u(3)) = 1        
          !---------------------------------------------------------------------
          ! Find all L-1 nbors of the new patch, and cancel the refinement
          !---------------------------------------------------------------------
          nbor => new%nbor
          do while (associated(nbor))
            if (nbor%task%level == new%task%level-1) then
              nbor_patch => patch%cast2patch (nbor%task)
              guards = .false.
              if (old_patch%get_overlap (nbor_patch, guards, l, u)) then
                if (verbose > 2) &
                  write (stdout,'(a,i6,3(2x,i2,":",i2))') &
                    '    level-1 overlap: id =', nbor_patch%id, &
                         l(1),u(1), l(2),u(2), l(3),u(3)
                trefine (l(1):u(1), l(2):u(2), l(3):u(3)) = 0    
              end if
            end if
            nbor => nbor%next
          end do
          !---------------------------------------------------------------------
          ! If any refinement marking remains, emit messages, set bits%support
          !---------------------------------------------------------------------
          l = old_patch%li
          u = old_patch%ui
          if (any(trefine(l(1):u(1), l(2):u(2), l(3):u(3)) == 1)) then
            call new_patch%set (bits%support)
            if (verbose > 1) then
              write (stdout,'(2(a,i6,":",i2,2x))') &
                'id =', old_patch%id, old_patch%level, 'support request from', &
                new_patch%id, new_patch%level
              write (io_unit%log,'(2(a,i6,":",i2,2x))') &
                'id =', old_patch%id, old_patch%level, 'support request from', &
                new_patch%id, new_patch%level
            end if
          else if (verbose > 3) then
            write (stdout,'(2(a,i6,":",i2,2x))') &
              'id =', old_patch%id, old_patch%level, 'support not requested from', &
              new_patch%id, new_patch%level
          end if
          deallocate (trefine)
        end if
      end if
      old => old%next
    end do
    call tlist%lock%unset ('check_support')
  end if
  !-----------------------------------------------------------------------------
  call trace%end (itimer, detailed_timer=detailed_timer)
END SUBROUTINE check_support

!===============================================================================
!> Decide if there are level L+2 tasks we need to support be refining
!===============================================================================
FUNCTION need_to_support (self, tlist, link) RESULT (refine)
  class(refine_t):: self
  class(list_t):: tlist
  class(link_t), pointer:: link
  integer:: refine
  !-----------------------------------------------------------------------------
  class(link_t), pointer:: nbor1, nbor2, nbor3
  class(patch_t), pointer:: task, nbtask
  type(patch_t):: patch
  logical:: guards
  integer:: l(3), u(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%need_to_support')
  task => patch%cast2patch(link%task)
  !-----------------------------------------------------------------------------
  ! Find L+2 and higher tasks, and check for overlap of their guard zones
  !-----------------------------------------------------------------------------
  call tlist%lock%set ('need_to_support')
  nbor1 => tlist%head
  do while (associated(nbor1))
    if (nbor1%task%level >= link%task%level+2) then
      nbtask => patch%cast2patch (nbor1%task)
      guards = .true.
      if (task%get_overlap (nbtask, guards, l, u)) then
        task%irefine(l(1):u(1),l(2):u(2),l(3):u(3)) = 1
      end if
    end if
    nbor1 => nbor1%next
  end do
  !-----------------------------------------------------------------------------
  ! Find L+1 nbor-nbors, check for overlap of their interiors
  !-----------------------------------------------------------------------------
  call link%lock%set ('need_to_support')
  nbor1 => link%nbor
  do while (associated(nbor1))
    if (nbor1%task%level == link%task%level+1) then
      nbtask => patch%cast2patch (nbor1%task)
      guards = .false.
      if (task%get_overlap (nbtask, guards, l, u)) then
        task%irefine(l(1):u(1),l(2):u(2),l(3):u(3)) = 0
      end if
    end if
    nbor1 => nbor1%next
  end do
  call link%lock%unset ('need_to_support')
  call tlist%lock%unset ('need_to_support')
  l = task%mesh%li
  u = task%mesh%ui
  refine = maxval(task%irefine(l(1):u(1),l(2):u(2),l(3):u(3)))
  call trace%end ()
END FUNCTION need_to_support

!===============================================================================
!> To check if a patch at level L may be removed w/o patches at level L+1 loosing
!> support it is enough to check if there are nbors at level L+1, because if
!> there are they would have, at the place where it had level L support, instead
!> only a level L-1 patch.  
!===============================================================================
FUNCTION check_derefine_support (self, tlist, link, check_only) RESULT (refine)
  USE task_mod
  class(refine_t):: self
  class(list_t):: tlist
  class(link_t), pointer:: link
  logical, optional:: check_only
  integer:: refine
  !.............................................................................
  class(link_t), pointer:: nbor, nbors
  integer:: n_nbor1, n_nbor2, n_added
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  refine = -1
  if (link%task%level == self%levelmin) &
    return
  call trace%begin('refine_t%check_derefine_support', itimer=itimer, detailed_timer=detailed_timer)
  !-----------------------------------------------------------------------------
  call tlist%lock%set ('check_derefine_support')
  nbor => link%nbor
  do while (associated(nbor))
    if (nbor%task%level == link%task%level+1) then
      refine = 0
      exit
    end if
    nbor => nbor%next
  end do
  if (verbose > 2) &
    write (stdout,*) link%task%id, 'check_derefine_support =', refine
  !-----------------------------------------------------------------------------
  ! Release lock
  !-----------------------------------------------------------------------------
  call tlist%lock%unset ('check_derefine_support')
  !-----------------------------------------------------------------------------
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION check_derefine_support

!===============================================================================
!> Create a new patch using an existing "parent" patch.
!> This procedure is used frequently during refinement.
!===============================================================================
SUBROUTINE make_child_patch (self, tlist, parent, pos, child_size, &
                             child)
  class(refine_t):: self
  class(list_t) :: tlist
  class(solver_t), pointer :: parent, child
  real(8), dimension(3), intent(in) :: pos, child_size
  !.............................................................................
  class(link_t), pointer :: child_link, child_nbor
  integer:: j, count, l(3), u(3), n_added
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%make_child_patch', itimer=itimer, detailed_timer=detailed_timer)
  if (parent%rank /= mpi%rank) then
    print *, 'make_child_patch: mpi%rank, parent%rank =', mpi%rank, parent%rank
    call io%abort ('make_child_patch: wrong MPI rank')
  end if
  !-----------------------------------------------------------------------------
  ! Allocate a new task, 1/nsplit the size of original one, and adjust its
  ! position in each direction with "pos", in units of the new size.
  !-----------------------------------------------------------------------------
  ! Using sourced allocation below effectively clones the parent.  This is
  ! necessary, to pick up properties that may have been set at the experiment
  ! level, which would otherwise be left undefined, or in inconsistent state
  !-----------------------------------------------------------------------------
  call parent%lock%set ('parent')
  allocate (child, source=parent)
  call parent%lock%unset ('parent')
  call child%clone_mem_accounting
 !-----------------------------------------------------------------------------
  child%mem_allocated = .false.                      ! make sure to get ..
  nullify (child%mem)                                ! .. new mem
  nullify (child%mesh)                               ! .. new mesh
  nullify (child%lock)                               ! .. new lock
  nullify (child%mem_lock)                           ! .. new mem_lock
  child%id = 0                                       ! ensure a new ID
  child%size = child_size                            ! new size
  child%position = pos                               ! patch position
  child%llc_cart = child%position - &                ! mid point position
    merge(0.5,0.0,child%n>1) * child%size            ! lower left corner
  child%llc_nat = child%llc_cart                     ! FIXME!: this form only applies for Cartesian coords.
  child%n_needed = 1                                 ! reset to initial
  child%istep = 0                                    ! for diagnostics
  !-----------------------------------------------------------------------------
  ! Force consistent resolution in new patches, overriding the resolution
  ! inherited from the parent patch.
  !-----------------------------------------------------------------------------
  child%n = shared%patch_n
  !-----------------------------------------------------------------------------
  ! Clear flags -- explictly to show which ones are relevant
  !-----------------------------------------------------------------------------
  !child%status = 0
  call child%clear(bits%root_grid)
  call child%clear(bits%ready)
  call child%clear(bits%busy)
  call child%clear(bits%frozen)
  child%boundaries%bits = 0
  !-----------------------------------------------------------------------------
  ! Append links to the new patch in the patch list, sorted by quality
  ! Neighbour lists are constructed at the end of `split_patch`.
  !-----------------------------------------------------------------------------
  allocate (child_link)
  call child_link%init                               ! initialise the lock
  child_link%task => child                           ! associate task with link
  child_link%parent => parent%link                   ! used in remove_patch
  child%link => child_link                           ! associate link with task
  child%parent => parent                             ! task relation
  !-----------------------------------------------------------------------------
  ! Initialise the new patch
  !-----------------------------------------------------------------------------
  child%level = parent%level + 1                     ! must be set before init!
  call child%set (bits%frozen)
  call child%init()                                  ! patch, for a new id
  call child%clear (bits%frozen)
  child%level = parent%level + 1                     ! override patch_t%init
  child%time  = parent%time                          ! same time
  child%it    = 1                                    ! time index
  child%new   = 2                                    ! new time slot
  child%iit   = 1                                    ! single time slice
  child%iit(child%nt) = child%new                    ! new slot
  child%dt(:) = parent%dt(parent%it) / self%ratio    ! reset the time-step
  child%t(:) = child%time                            ! sinle time
  child%parentid = parent%id                         ! for debugging
  child%check_refine_next = 0                        ! allow recursive refine
  child%iout = child%time/io%out_time+1              ! previous output index
  child%out_next = child%iout*io%out_time            ! next output time
  child%dnload_time = -1.0                           ! for duplicate call test
  !-----------------------------------------------------------------------------
  ! Increment I/O write counter, to avoid repeated writing of metadata
  !-----------------------------------------------------------------------------
  call data_io%update_counters (child, +1)
  call child%log ('created')
  !-----------------------------------------------------------------------------
  if (verbose > 0) &
    write (io_unit%log,'(f12.5,2x,a,2i6,1p,g14.5)') &
      wallclock(), 'refine_t%make_child_patch, parent, child, time =', &
      parent%id, child%id, child%time
  !-----------------------------------------------------------------------------
  ! Check if the new patch, which might have been created to support another
  ! patch, possibly should be immediately frozen. If that is the case, it still
  ! needs to be added to the task list, in case it is needed for support, but by
  ! setting bits%frozen it is prevented from being added to the ready queue.
  !-----------------------------------------------------------------------------
  if (child%time > io%end_time) then
    call child%set (bits%frozen)
  end if
  !-----------------------------------------------------------------------------
  ! Set the bits%support here, which causes it to be set on the rank nbors
  ! if this is a boundary task.  There, a check for support is then triggered,
  ! where the nbor lists also in the virtual nbors are refreshed.  On the owner
  ! rank the bits%support will remain until next time the task is updated, when 
  ! code in dispatcher0_t%update() sees the bit, call init_nbors(), and clears
  ! the bit.
  !-----------------------------------------------------------------------------
  call child%set (bits%support)
  !-----------------------------------------------------------------------------
  ! Make parent a neighbor of the new child patch.  This will be overwritten
  ! in `init_nbor_nbors` but is necessary here to permit interpolation of
  ! values from parent to child.
  !-----------------------------------------------------------------------------
  allocate (child_nbor)
  child_nbor%link => parent%link
  child_nbor%task => parent
  nullify (child_link%nbor)
  child_link%nbors_by_level => child_nbor
  !------------------------------------------------------------------------
  ! Download values from parent and neighbors, except at time = 0, to allow
  ! recursive refinement there. If time = 0, patch should have been filled
  ! with values during `child%init()` call.
  !------------------------------------------------------------------------
  !$omp atomic update
  parent%n_needed = parent%n_needed + 2
  if (parent%time > 0.0 .and. child%restart <= 0) then
    call download%download_link (child_link, all_cells=.true.)
  end if
  !-----------------------------------------------------------------------------
  ! Add the child task to the task list, including giving it an nbor list etc,
  ! and add it also to the ready queue
  !-----------------------------------------------------------------------------
  call tlist%add_new_link (child_link)
  call tlist%queue_by_time (child_link)
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    l = parent%index_only (child%position + child%mesh%lf)
    u = parent%index_only (child%position + child%mesh%uf)
    l = max(l,parent%mesh%li)
    u = min(u,parent%mesh%ui)
    print 1,'make_child_patch: parent, level, iout, time, size(3), position(3) =', &
      parent%id, parent%level, parent%time, parent%size, parent%position , &
      maxval(parent%mem(l(1):u(1),l(2):u(2),l(3):u(3),parent%idx%d,parent%it,1)), &
      l, u, parent%contains(child)
    1 format(a,i6,i3,f9.6,2(1x,3f8.4),1pe12.2,2(1x,3i3),l2)
    l(:) = parent%mesh(:)%li
    u(:) = parent%mesh(:)%ui
    print 1,'                   child, level, iout, time, size(3), position(3) =', &
      child%id, child%level, child%time, child%size, child%position , &
      maxval(child%mem(l(1):u(1),l(2):u(2),l(3):u(3),child%idx%d,child%it,1))
    flush (io%output)
  end if
  !-----------------------------------------------------------------------------
  call trace%end (itimer, detailed_timer=detailed_timer)
END SUBROUTINE make_child_patch

!===============================================================================
!> Check that the patch dimensions are the desired cubical ones
!===============================================================================
FUNCTION check_cubical (self, patch) RESULT (refine_it)
  class(refine_t):: self
  class(solver_t), pointer:: patch
  integer:: refine_it
  !.............................................................................
  if (all(patch%n==shared%patch_n)) then
    refine_it = -1
  else
    refine_it = 1
    patch%msplit = 1
  end if
END FUNCTION check_cubical

!===============================================================================
!> Refine / derefine on ratio of Jeans length to cell size (max_jeans/min_jeans)
!===============================================================================
FUNCTION check_Jeans (self, patch) RESULT(refine_it)
  USE math_mod
  USE scaling_mod
  class(refine_t):: self
  class(solver_t), pointer:: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: Jeans_num, dmin, c
  real, dimension(:,:,:), allocatable:: pg
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e
  real, allocatable:: test(:,:,:)
  logical:: masked
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_jeans')
  call self%sanity_check (min_jeans, max_jeans, 'check_Jeans')
  l = 1
  u = patch%gn
  d => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%it,1)
  allocate(pg(patch%gn(1),patch%gn(2),patch%gn(3)))
  select type (patch)
  class is (solver_t)
    pg = patch%gas_pressure() + tiny(1.0)
  class default
    call mpi%abort("Unrecognised patch type. Abort!")
  end select
  dmin = patch%fminval(d)
  if (dmin <= 0.0) then
    imin = minloc(d)
    write (io%output,*) &
      'check_Jeans ERROR: id, dmin, imin =', patch%id, dmin, imin
    Jeans_num = min_jeans
  else
    allocate (test(patch%gn(1),patch%gn(2),patch%gn(3)))
    c =  minval(patch%mesh%d)*sqrt(scaling%grav/(math%pi*patch%gamma))
    !print *, patch%id, scaling%grav, 'grav'
    test = sqrt(pg)/(c*d)
    !print *, patch%id, patch%fminval(d), patch%fmaxval(d), 'min/max d'
    !print *, patch%id, patch%fminval(test), patch%fmaxval(test), 'min/max L_J/ds'
    !---------------------------------------------------------------------------
    ! The Jeans number in the unrefined part
    !---------------------------------------------------------------------------
    Jeans_num = minval(test, mask=patch%irefine <= 0)
    !---------------------------------------------------------------------------
    ! Compute a per-cell criterion, which sets patch%irefine to +1 where it
    ! needs refinement, and sets it to 0 where the current resolution is needed.
    ! By implication, where it does not set patch%irefine, and where it remains
    ! at the intial -1 value, derefinement is allowed
    !---------------------------------------------------------------------------
    where (test < min_jeans)
      patch%irefine = +1
    !---------------------------------------------------------------------------
    ! Jeans number not small enough, or overlying refinement
    !---------------------------------------------------------------------------
    else where (test < max_jeans)
      patch%irefine = max (patch%irefine, 0)
    end where
    deallocate (test)
    refine_it = 0
    if (Jeans_num > max_jeans) refine_it = -1
    if (Jeans_num < min_jeans) refine_it = +1
    !---------------------------------------------------------------------------
    if (verbose > 2 .or. &
       (verbose > 1 .and. patch%level < self%levelmax .and. refine_it > 0) .or. &
       (verbose > 1 .and. patch%level > self%levelmin .and. refine_it < 0)) then
      masked = any(patch%irefine >= 0)
      write (io_unit%output,'(a,i6,1p,3e9.2,i4,e9.2,l3)') &
        'Jeans criterion: id, max, num, min, refine, maxval(d)', &
        patch%id, max_jeans, Jeans_num, min_jeans, refine_it, maxval(d), &
        masked
    end if
    !---------------------------------------------------------------------------
  end if
  deallocate (pg)
  call trace%end()
END FUNCTION check_Jeans

!===============================================================================
!> Refine / derefine on vorticity (velocity difference per cell)
!===============================================================================
FUNCTION check_compress (self, patch) RESULT(refine_it)
  USE math_mod
  USE scaling_mod
  class(refine_t):: self
  class(solver_t), pointer:: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: w_num, wmin
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: test
  logical:: masked
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_compress', itimer=itimer, detailed_timer=detailed_timer)
  l = 1
  u = patch%gn
  allocate (test(patch%gn(1),patch%gn(2),patch%gn(3)))
  select type (patch)
  class is (solver_t)
    call patch%compression_magnitude (test)
  class default
    call mpi%abort("Unrecognised patch type. Abort!")
  end select
  test = minval(patch%mesh%d)*test
  w_num = maxval(test, mask=patch%irefine <= 0)
  !---------------------------------------------------------------------------
  ! Evaluate a per-cell criterion, which sets patch%irefine to +1 where it
  ! needs refinement, and sets it to 0 where the current resolution is needed.
  ! By implication, where it does not set patch%irefine, and where it remains
  ! at the intial -1 value, derefinement is allowed
  !---------------------------------------------------------------------------
  where (test > max_compress)
    patch%irefine = +1
  !---------------------------------------------------------------------------
  ! compression number not small enough, or overlying refinement
  !---------------------------------------------------------------------------
  else where (test > min_compress)
    patch%irefine = max (patch%irefine, 0)
  end where
  refine_it = 0
  if (w_num > max_compress) refine_it = +1
  if (w_num < min_compress) refine_it = -1
  !---------------------------------------------------------------------------
  if (verbose > 2 .or. &
     (verbose > 1 .and. patch%level < self%levelmax .and. refine_it > 0) .or. &
     (verbose > 1 .and. patch%level > self%levelmin .and. refine_it < 0)) then
    masked = any(patch%irefine >= 0)
    write (io_unit%output,'(a,i6,1p,3e9.2,i4,e9.2,l3)') &
      'compression criterion: id, max, num, min, refine, maxval(test)', &
      patch%id, max_compress, w_num, min_compress, refine_it, maxval(test), &
      masked
  end if
  deallocate (test)
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION check_compress

!===============================================================================
!> Refine / derefine on vorticity (velocity difference per cell)
!===============================================================================
FUNCTION check_vorticity (self, patch) RESULT(refine_it)
  USE math_mod
  USE scaling_mod
  class(refine_t):: self
  class(solver_t), pointer:: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: w_num, wmin
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: test
  logical:: masked
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_vorticity', itimer=itimer, detailed_timer=detailed_timer)
  l = 1
  u = patch%gn
  allocate (test(patch%gn(1),patch%gn(2),patch%gn(3)))
  select type (patch)
  class is (solver_t)
    call patch%vorticity_magnitude (test)
  class default
    call mpi%abort("Unrecognised patch type. Abort!")
  end select
  test = minval(patch%mesh%d)*test
  w_num = maxval(test, mask=patch%irefine <= 0)
  !---------------------------------------------------------------------------
  ! Evaluate a per-cell criterion, which sets patch%irefine to +1 where it
  ! needs refinement, and sets it to 0 where the current resolution is needed.
  ! By implication, where it does not set patch%irefine, and where it remains
  ! at the intial -1 value, derefinement is allowed
  !---------------------------------------------------------------------------
  where (test > max_vorticity)
    patch%irefine = +1
  !---------------------------------------------------------------------------
  ! vorticity number not small enough, or overlying refinement
  !---------------------------------------------------------------------------
  else where (test > min_vorticity)
    patch%irefine = max (patch%irefine, 0)
  end where
  refine_it = 0
  if (w_num > max_vorticity) refine_it = +1
  if (w_num < min_vorticity) refine_it = -1
  !---------------------------------------------------------------------------
  if (verbose > 2 .or. &
     (verbose > 1 .and. patch%level < self%levelmax .and. refine_it > 0) .or. &
     (verbose > 1 .and. patch%level > self%levelmin .and. refine_it < 0)) then
    masked = any(patch%irefine >= 0)
    write (io_unit%output,'(a,i6,1p,3e9.2,i4,e9.2,l3)') &
      'vorticity criterion: id, max, num, min, refine, maxval(test)', &
      patch%id, max_vorticity, w_num, min_vorticity, refine_it, maxval(test), &
      masked
  end if
  deallocate (test)
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION check_vorticity

!===============================================================================
!> Refine on shocks, with threshold set by tolerance parameter `tol_shock`.
!===============================================================================
FUNCTION shock_and_cd_detector (self, patch) RESULT (refine)
  USE mesh_mod
  class(refine_t):: self
  class(solver_t), pointer:: patch
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
  class(solver_t), pointer:: patch
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
  class(solver_t), pointer :: patch
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
    call trace%end ()
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

!===============================================================================
!> Check that the derefinement criterion is sane
!===============================================================================
SUBROUTINE sanity_check (self, min, max, label, reverse)
  class(refine_t):: self
  real:: min, max
  character(len=*):: label
  logical, optional:: reverse
  !.............................................................................
  if (max == 0.0) then
    if (present(reverse)) then
      min = max/self%safe_ratio
    else
      max = min*self%safe_ratio
    end if
  else if (max < min*2.0) then
    if (present(reverse)) then
      min = max/2.0
    else
      max = min*2.0
    end if
    if (io%master) then
      write (stdout,*) 'WARNING: '//trim(label)//' max value set to', max
    end if
  end if
END SUBROUTINE sanity_check

END MODULE refine_mod
