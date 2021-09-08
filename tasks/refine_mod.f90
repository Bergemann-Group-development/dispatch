!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> This module handles checking various refinement criteria.  Call hierarchy:
!>
!> task_list_t%update
!>   check
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
!> the patch that should potentially be created, and not as a criterion on the
!> one that needs the support (to avoid having to find a sub-region of the
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
  USE timer_mod
  USE data_io_mod
  USE math_mod
  USE scaling_mod
  USE eos_mod
  USE region_mod
  implicit none
  private
  !
  type:: refine_t
    logical:: on                              ! allow external access
    integer:: levelmin=0                      ! Min level of refinment
    integer:: levelmax=0                      ! Max level of refinment
    integer:: ratio=2                         ! default
    real:: safe_ratio=2.1                     ! safe derefine ratio
    class(list_t), pointer:: task_list => null()
  contains
    procedure:: init
    procedure:: check
    procedure:: refine_needed
    procedure:: selective_refine
    procedure:: post_process
    procedure:: make_child_patch
    procedure:: check_jeans
    procedure:: shock_and_cd_detector
    procedure:: evaluate
    procedure:: refine_region
    procedure:: gradient_detector
    procedure:: check_derefine_support
    procedure:: check_support
  end type
  !
  integer:: check_interval=3                  ! Number of steps betwen checks
  real:: min_dx=0.0, max_dx=1e32              ! Min and max grid size
  real:: refine_from_t=0.0                    ! Time limits for occurence of refinement
  real:: refine_until_t=1e30                  ! Time limits for occurence of refinement
  real:: min_jeans=0.0                        ! Jeans number above which a patch is flagged for derefinement
  real:: max_shock=0.0, max_contact=0.0       ! Tolerance values for shock/discontinuity detector refinement
  real:: min_shock=0.0, min_contact=0.0       ! Tolerance values for shock/discontinuity detector derefinement
  real:: max_compress=0.0                     ! Tolerance value for compression (negative divergence)
  integer:: max_compress_level=-1             ! Max refinement level on compression
  real:: max_vorticity=0.0                    ! Tolerance value for vorticity
  real:: derefine_limit=0.8                   ! Fraction of criterion to allow derefinement
  integer, parameter:: ngradvar=10            ! Number of gradient variables
  real, dimension(3):: region_llc=0.0         ! Sometimes you want to refine a region (manually)
  real, dimension(3):: region_urc=0.0         ! upper-right-corner of region
  integer:: region_maxlevel=5                 ! max region refinement level
  real, dimension(ngradvar):: max_grad=0.0    ! Tolerance values for gradient detector refinement
  real, dimension(ngradvar):: min_grad=0.0    ! Tolerance values for gradient detector derefinement
  character(len=32), dimension(ngradvar):: grad_var=""
  logical:: force_cubical=.false.             ! Force cubical patches
  logical:: detailed_timer=.false.            ! Lump all AMR together when false
  logical:: do_dump=.false.                   ! dump parent & child
  integer:: n_locks=1                         ! number of locks to use
  integer:: verbose=0                         ! Noise level
  integer:: n_criteria=0                      ! Number of active criteria
  integer:: wait_refine_nsteps=0              ! How many time steps to wait before allowing refinement
  logical:: quadrature=.false.                ! Combine criteria in quadrature
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
  logical:: on=.false.                        ! Flag to turn the entire process on or off
  logical, save:: first_time = .true.
  namelist /refine_params/ verbose, on, levelmin, levelmax, check_interval, &
    ratio, refine_from_t, refine_until_t, min_dx, max_dx, derefine_limit, &
    min_jeans, max_shock, min_shock, max_contact, min_contact, max_grad, &
    min_grad, max_vorticity, max_compress, max_compress_level, quadrature, &
    grad_var, region_llc, region_urc, region_maxlevel, force_cubical, n_locks, &
    wait_refine_nsteps, detailed_timer, do_dump
  character(len=120):: id = &
    '$Id$'
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%init')
  call trace%print_id (id)
  call io%lock%set ('refine_params')
  if (first_time) then
    first_time=.false.
    levelmin = shared%levelmin
    levelmax = shared%levelmax
    rewind (io%input)
    read (io%input, refine_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, refine_params)
    end if
    call io%namelist_warning ('refine_params', iostat)
    if (max_compress_level == -1) &
      max_compress_level = levelmax
    write (io%namelist, refine_params)
    !---------------------------------------------------------------------------
    ! If AMR is active, we must use check_filled, since nbor levels are in
    ! decreasing order, to save computing time.
    !---------------------------------------------------------------------------
    if (on) then
      omp_lock%links = .true.
    else
      omp_lock%links = .false.
    end if
    shared%refine_on = on
    shared%f_zoom = 1.0/ratio
    shared%refine_ratio = ratio
    if (min_jeans     > 0.0) n_criteria = n_criteria+1
    if (max_compress  > 0.0) n_criteria = n_criteria+1
    if (max_vorticity > 0.0) n_criteria = n_criteria+1
    if (max_shock     > 0.0) n_criteria = n_criteria+1
    if (max_contact   > 0.0) n_criteria = n_criteria+1
    if (any(region_llc/=0) .or.any(region_urc/=0)) n_criteria = n_criteria+1
    write (io%output,'(a,i4)') ' number of AMR criteria:', n_criteria
    if (n_criteria ==0) self%on=.false.
  end if
  call io%lock%unset ('refine_params')
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
!> Check if refinement is desired and needed for the current task; if so, push
!> new tasks onto the (ready) queue with the same time as the parent, i.e., to
!> the head of the queue.  If the task is a leaf task, and the level of
!> refinement does not appear to be needed, request derefinement.
!===============================================================================
SUBROUTINE check (self, tasklist, link, was_refined, was_derefined)
  class(refine_t):: self
  class(list_t), target:: tasklist
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task
  class(solver_t), pointer:: patch
  class(list_t), pointer:: child_list
  logical, intent(out):: was_refined, was_derefined
  real:: dxmin, dxmax
  integer:: refine, id, level, id0
  integer, save:: itimer=0
  !.............................................................................
  was_refined = .false.
  was_derefined = .false.
  !-----------------------------------------------------------------------------
  ! Punt on refinement if the task is very new, or has the init_nbors bit set
  !-----------------------------------------------------------------------------
  if (.not.self%on .or. &
      link%task%istep < wait_refine_nsteps) then
        return
  end if
  if (verbose > 2) then
    if (link%task%is_set (bits%init_nbors)) &
      write(io_unit%log,'(i6,f12.6,4x,a,3x,2l1)') &
        link%task%id, link%task%time, 'has bits%init_nbors', &
        link%task%is_set (bits%boundary), &
        link%task%is_set (bits%virtual)
  end if
  !---------------------------------------------------------------------------
  call trace%begin('refine_t%check', itimer=itimer)
  self%task_list => tasklist
  id0 = link%task%id
  task => link%task
  select type (task)
  class is (solver_t)
  patch => task
  end select
  id = patch%id
  level = patch%level
  !-----------------------------------------------------------------------------
  dxmin = minval(patch%ds)
  dxmax = maxval(patch%ds)
  if (verbose > 1) then
    write (io_unit%log,'(f12.6,a,i6,":",i2.2,i6,2i4,l3,1p,2(2x,2e11.2))') &
      wallclock(), pad_trace(1:indent)// &
      'refine_t%check: id,id0,istep,next,V,dx_min,dx_max =', &
      id, level, id0, patch%istep, patch%check_refine_next, &
      patch%is_set(bits%virtual), min_dx, dxmin, max_dx, dxmax
  end if
  !----------------------------------------------------------------------------
  ! Filter out patches that do not match the necessary criteria
  !----------------------------------------------------------------------------
  if (patch%is_clear(bits%virtual)  .and. &
      dxmin >  min_dx               .and. &
      dxmax <= max_dx               .and. &
      patch%time >= refine_from_t   .and. &
      patch%time <  refine_until_t) then
    if (verbose > 3) &
      self%task_list%verbose = verbose
    !---------------------------------------------------------------------------
    ! Check all refinement criteria.  If patch%frefine is already marked for
    ! refinement, then it has been done to request support of new level+2 patches
    !---------------------------------------------------------------------------
    call patch%support_arrays ()
    refine = self%refine_needed (link)
    !---------------------------------------------------------------------------
    ! Apply refinement criteria on potential child patch parts, but only if
    ! either refinement has been indicated, or if derefinement is allowed
    !---------------------------------------------------------------------------
    if (refine /= 0) then
      allocate (child_list)
      call child_list%init ()
      call self%selective_refine (patch, child_list, refine, &
        was_refined, was_derefined)
      call self%post_process (patch, child_list)
      call child_list%dealloc ()
      deallocate (child_list)
    end if
    if (verbose > 2) &
      write (io_unit%log,*) 'deallocated irefine in refine', mpi%rank, id
    !---------------------------------------------------------------------------
  end if
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE check

!===============================================================================
!> Check if refinement is desired and needed for the current task
!===============================================================================
FUNCTION refine_needed (self, link) RESULT (refine)
  class(refine_t)          :: self
  class(link_t), pointer   :: link
  integer                  :: refine
  !.............................................................................
  class(task_t), pointer   :: task
  class(solver_t), pointer :: patch
  class(extras_t), pointer :: extras
  integer(kind=int8), pointer, contiguous:: request(:,:,:)
  integer, dimension(3)    :: l, u
  integer                  :: id, level, lrefine
  integer, save            :: itimer=0
  !-----------------------------------------------------------------------------
  refine = 0
  if (.not. self%on) return
  call trace%begin('refine_t%needed', itimer=itimer, detailed_timer=detailed_timer)
  !-----------------------------------------------------------------------------
  ! Refinement criteria, starting with assumed derefinement
  !-----------------------------------------------------------------------------
  task => link%task
  id = task%id
  level = task%level
  select type (task)
  class is (solver_t)
  patch => task
  refine = -1
  !-----------------------------------------------------------------------------
  ! Support refinement
  !-----------------------------------------------------------------------------
  refine = max(refine, self%check_support (link))
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, 0., 0., 'support'
    1 format(f12.6,a,"refine_needed: id,refine,min,max,type=", &
            i6,":",i2.2,i3,1p,2e11.2,2x,a)
  !-----------------------------------------------------------------------------
  ! Check the rest of the criteria less often, but allow derefine only when
  ! all criteria are used
  !-----------------------------------------------------------------------------
  if (patch%istep < patch%check_refine_next) then
    refine = max(0,refine)
    call trace%end (itimer, detailed_timer=detailed_timer)
    return
  end if
  !-----------------------------------------------------------------------------
  patch%check_refine_next = patch%istep + check_interval
  call patch%support_arrays ()
  patch%frefine = -derefine_limit
  !-----------------------------------------------------------------------------
  ! Extras refinement
  !-----------------------------------------------------------------------------
  extras => task
  end select
  refine = max(refine, extras%check_refine ())
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, &
    minval(patch%frefine), maxval(patch%frefine), 'extras'
  !-----------------------------------------------------------------------------
  ! Jeans criterion
  !-----------------------------------------------------------------------------
  if (min_jeans > 0.0) refine = max (refine, self%check_Jeans(patch))
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, &
    minval(patch%frefine), maxval(patch%frefine), 'Jeans'
  !-----------------------------------------------------------------------------
  ! shock and/or contact discontinuity detector
  !-----------------------------------------------------------------------------
  if (max_shock > 0.0 .or. max_contact > 0.0) &
    refine = max (refine, self%shock_and_cd_detector(patch))
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, &
    minval(patch%frefine), maxval(patch%frefine), 'shock'
    !---------------------------------------------------------------------------
  ! compression detector
  !-----------------------------------------------------------------------------
  if (max_compress > 0.0) &
    refine = max (refine, check_compress (self, patch))
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, &
    minval(patch%frefine), maxval(patch%frefine), 'compress'
  !-----------------------------------------------------------------------------
  ! vorticity detector
  !-----------------------------------------------------------------------------
  if (max_vorticity > 0.0) &
    refine = max (refine, check_vorticity (self, patch))
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, &
    minval(patch%frefine), maxval(patch%frefine), 'vorticity'
  !-----------------------------------------------------------------------------
  ! gradient detector
  !-----------------------------------------------------------------------------
  if (any(max_grad > 0.0)) &
    refine = max (refine, self%gradient_detector(patch))
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, &
    minval(patch%frefine), maxval(patch%frefine), 'gradient'
  !-----------------------------------------------------------------------------
  ! refine on (a) user-specified region(s)
  !-----------------------------------------------------------------------------
  if (any(region_llc /= 0.0) .and. any(region_urc /= 0.0)) &
    refine = max (refine, self%refine_region(patch))
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, &
    minval(patch%frefine), maxval(patch%frefine), 'user'
  ! !-----------------------------------------------------------------------------
  ! ! Check if support constraints allow derefinement
  ! !-----------------------------------------------------------------------------
  ! refine = max (refine, self%check_derefine_support (link))
  ! if (verbose > 1) write (stdlog,1) &
  !   wallclock(), pad_trace(1:indent), id, level, refine, &
  !   minval(patch%frefine), maxval(patch%frefine), 'derefine support'
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION refine_needed

!===============================================================================
!> Loop over the ratio x ratio x ratio positions and determine whether to refine
!> or derefine.  Child tasks that are created are added to a temporary task
!> list, and are afterwards given new nbor lists, via which their new nbors
!> in turn are given notice to create new nbor lists.  This includes the current
!> task, which will do so after the task update.  With the new child tasks in
!> its nbor list it will normally be prevented from advancing further, until
!> the child tasks have caught up.
!>
!> As part of their creation, the new child tasks are immediately advanced one
!> step, since they are created with filled guard zones.  As part of that
!> process they are also added to the task list, thus making it possible to
!> give them complete nbor lists at the bottom of this procedure, where they
!> are also scheduled for further updates, as part of the post-update task
!> handling.
!===============================================================================
SUBROUTINE selective_refine (self, patch, child_list, refine, &
                             was_refined, was_derefined)
  class(refine_t)          :: self
  class(solver_t), pointer :: patch
  class(list_t),   pointer :: child_list
  integer                  :: refine
  logical                  :: was_refined, was_derefined
  !.............................................................................
  class(solver_t), pointer :: child
  class(link_t),   pointer :: link, old_link, next
  class(nbor_t),   pointer :: nbor
  class(task_t),   pointer :: task
  integer                  :: i, j, k, id, level, id1
  integer, dimension(3)    :: l, u, n
  real(8), dimension(3)    :: size, pos
  logical                  :: already_refined, refined
  integer                  :: refine_it, ratio_per_dim(3), key(4)
  integer                  :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%selective_refine', itimer=itimer, &
                    detailed_timer=detailed_timer)
  refine = -1
  was_refined = .false.
  was_derefined = .false.
  refined = .false.
  id = patch%id
  level = patch%level
  link => patch%link
  call link%lock%set ('selective_refine')
  nbor => patch%link%nbor
  call link%lock%unset ('selective_refine')
  ratio_per_dim = merge(self%ratio, 1, patch%n>1)
  !-----------------------------------------------------------------------------
  ! First, consider if any part should be refined
  !-----------------------------------------------------------------------------
  if (patch%level < self%levelmax) then
    l = patch%mesh%li
    u = patch%mesh%ui
    if (maxval(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3))) >= 1.0) then
      do k=1,ratio_per_dim(3)
        do j=1,ratio_per_dim(2)
          do i=1,ratio_per_dim(1)
            refine = max (refine, check_support2 (self, patch, [i,j,k]))
            !-------------------------------------------------------------------
            ! Lower and upper index in parent patch
            !-------------------------------------------------------------------
            n = patch%mesh%n/ratio_per_dim
            l = patch%mesh%li + ([i,j,k]-1)*n
            u = max(1, l + n - 1)
            !-------------------------------------------------------------------
            ! Should this part be refined or not, according to patch%frefine
            ! accumulated from all active conditions?
            !-------------------------------------------------------------------
            refine_it = floor(maxval( &
              patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3)) ))
            refine = max(refine, refine_it)
            !-------------------------------------------------------------------
            ! Size and position
            !-------------------------------------------------------------------
            size = patch%size/ratio_per_dim
            pos = patch%position + ([i,j,k]-(ratio_per_dim+1)*0.5)*size
            key = patch%hash_key (level=patch%level+1, size=size, position=pos)
            call self%task_list%hash_table%get (key, old_link)
            !-------------------------------------------------------------------
            ! Is this part refined before, or should it be now?
            !-------------------------------------------------------------------
            if (associated(old_link)) then
              if (verbose > 1) &
                write (stdlog,'(a,i6,3(1x,i2,":",i2.2))') &
                  'selective_refine: existed=', &
                  old_link%task%id, l(1),u(1), l(2),u(2), l(3),u(3)
              refined = .true.
            else
              if (refine_it > 0) then
                was_refined = .true.
                refined = .true.
                !---------------------------------------------------------------
                ! Create a refined patch here
                !---------------------------------------------------------------
                if (verbose > 1) &
                  write (io_unit%log,'(f12.6,a,i6,1p2g16.5,2x,3(i4,":",i2.2))') &
                    wallclock(), pad_trace(1:indent)// &
                    'refine', patch%id, &
                    patch%fmaxval(patch%mem(:,:,:,patch%idx%d,patch%it,1)), &
                    maxval(patch%mem(l(1):u(1),l(2):u(2),l(3):u(3), patch%idx%d, &
                    patch%it,1)), l(1),u(1),l(2),u(2),l(3),u(3)
                call patch%log ('refine', 2)
                call make_child_patch (self, patch, pos, size, child)
                !---------------------------------------------------------------
                ! Allocate a tempory task link and append it to the child list,
                ! making sure not to change the task%link, which is needed below
                !---------------------------------------------------------------
                allocate (link)
                link%task => child
                call child_list%append (link, keep=.true.)
                !---------------------------------------------------------------
                if (verbose > 0) &
                  write (io_unit%log,'(f12.6,2(a,i6,":",i2.2,2x),a,2l1,1p,e14.5)')  &
                    wallclock(), pad_trace(1:indent)// &
                    'refine:'   ,  id, level, &
                    'to'        ,  child%id, child%level, &
                    'BV:', child%is_set(bits%boundary), child%is_set(bits%virtual), &
                    child%time
                flush (io_unit%log)
                !---------------------------------------------------------------
              end if
            end if
          end do
        end do
      end do
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Then, consider if the whole patch should be derefined.  This requires that
  ! it has no higher level patch inside (would trigger patch%refine > 0 there),
  ! and secondly that frefine(:,:,:) indicates that derefine would be OK
  !-----------------------------------------------------------------------------
  if (refine == -1                .and. &               ! derefine possible
      patch%level > self%levelmin .and. &               ! relevant?
      patch%check_refine_next > 0 .and. &               ! not immediately!
      patch%istep >= 4            .and. &               ! must take 4 steps
      .not.refined) then                                ! not if refined!
    l = patch%mesh%li
    u = patch%mesh%ui
    if (maxval(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3))) < 0.0) then
      !-----------------------------------------------------------------------
      if (verbose > 1) then
        write (io_unit%log,'(1x,a,i6,":",i2,2x,"B:",l1,3x,a,i6,3x,a,f12.6)') &
          'derefining!:', id, level, patch%is_set (bits%boundary), &
          'istep:', patch%istep, 'time:', patch%time
        flush (io_unit%log)
      else if (verbose > 0) then
        write (io_unit%log,'(f12.6,a,i6,":",i2.2,16x,a,2l1,1p,e14.5)') &
          wallclock(), pad_trace(1:indent)// &
          'derefine:' ,  id, level, 'BV:', &
          patch%is_set(bits%boundary), patch%is_set(bits%virtual), patch%time
      end if
      !-------------------------------------------------------------------------
      ! If task is to be removed, call the comprehensive remove procedure, which
      ! includes calling set_init_nbors(), queueing all nbors for immediate
      ! update of their nbor lists.
      !-------------------------------------------------------------------------
      write (stderr,'(a,i6,":",i2.2,f12.6)') 'derefine:', &
        patch%id, patch%level, patch%time
      call patch%log ('derefine', 2)
      call data_io%update_counters (patch, -1)
      call self%task_list%remove_and_reset (link)
      was_derefined = .true.
    end if
  end if
  call trace%end (itimer, detailed_timer=detailed_timer)
END SUBROUTINE selective_refine

!============================================ ===================================
!> Post process the list of child tasks created by refine%check()
!===============================================================================
SUBROUTINE post_process (self, patch, child_list)
  class(refine_t):: self
  class(solver_t), pointer:: patch
  class(list_t), pointer:: child_list
  !.............................................................................
  class(link_t), pointer:: link, next
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%post_process')
  !-----------------------------------------------------------------------------
  ! The new tasks in the child_list need to be 
  ! 1) added to the task list (done above)
  ! 2) have new nbor lists, which include the new siblings
  ! 3) have bits%no_support set on or off
  ! 4) be added to the queue if ready to update
  ! Since 4) depends on all new tasks having done 2) and 3) it is not done
  ! here rather than inside add_new_link
  !-----------------------------------------------------------------------------
  if (associated(child_list%head)) then
    !---------------------------------------------------------------------------
    ! Do the standard set of actions after adding new tasks
    !---------------------------------------------------------------------------
    link => child_list%head
    do while (associated(link))
      task => link%task
      select type (task)
      class is (patch_t)
        omp%thread_task = task%link%task%id
        call self%task_list%add_new_link (task%link)
        !.......................................................................
        if (verbose > 1) then
          write (stdout,'(a,i6,a,i6,":",i2.2,3i3,l3)') &
            'refine_mod:',patch%id,' add_new_link', task%id, task%key, &
            task%is_clear (bits%no_support)
        end if
        !.......................................................................
      end select
      link => link%next
    end do
  end if
  omp%thread_task = patch%id
  !-----------------------------------------------------------------------------
  ! Deallocate and delete the child list, after check_ready()
  !-----------------------------------------------------------------------------
  link => child_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (patch_t)
      !-------------------------------------------------------------------------
      ! Clear the busy bit, now that the task is ready to be scheduled.  The
      ! busy bit is set before the task is added to the task list, to avoid
      ! that it is sheduled by check_nbor calls in a multi-threaded setting.
      !-------------------------------------------------------------------------
      call link%task%clear (bits%busy)
      call self%task_list%check_ready (task%link)
    end select
    next => link%next
    deallocate (link)
    link => next
  end do
  call trace%end ()
END SUBROUTINE post_process

!===============================================================================
!> To check if a patch at level L may be removed w/o patches at level L+1
!> loosing support it is enough to check if there are nbors at level L+1,
!> because if there are they would have, at the place where it gave level L
!> support, instead only a level L-1 patch.  No locking is needed, since the
!> link task is the one being updated, and no other thread is allowed to change
!> the nbor list.
!===============================================================================
FUNCTION check_derefine_support (self, link, check_only) RESULT (refine)
  USE task_mod
  class(refine_t):: self
  class(link_t), pointer:: link
  logical, optional:: check_only
  integer:: refine
  !.............................................................................
  class(nbor_t), pointer:: nbor, nbors
  integer:: n_nbor1, n_nbor2, n_added
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  refine = -1
  if (link%task%level == self%levelmin) &
    return
  call trace%begin('refine_t%check_derefine_support', itimer=itimer, &
    detailed_timer=detailed_timer)
  !-----------------------------------------------------------------------------
  ! call link%lock%set ('check_derefine_support')
  nbor => link%nbor
  do while (associated(nbor))
    if (nbor%task%level == link%task%level+1 .and. nbor%needs_me) then
      refine = 0
      exit
    end if
    nbor => nbor%next
  end do
  if (verbose > 4) &
    write (stdout,*) link%task%id, 'check_derefine_support =', refine
  !-----------------------------------------------------------------------------
  ! Release lock
  !-----------------------------------------------------------------------------
  ! call link%lock%unset ('check_derefine_support')
  !-----------------------------------------------------------------------------
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION check_derefine_support

!===============================================================================
!> Create a new patch using an existing "parent" patch.
!> This procedure is used frequently during refinement.
!===============================================================================
SUBROUTINE make_child_patch (self, parent, pos, child_size, &
                             child)
  class(refine_t):: self
  class(task_t), pointer :: task
  class(solver_t), pointer :: parent, child
  real(8), dimension(3), intent(in) :: pos, child_size
  !.............................................................................
  class(nbor_t), pointer :: child_nbor
  integer:: j, count, l(3), u(3), n_added, unit
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%make_child_patch', itimer=itimer, detailed_timer=detailed_timer)
  if (parent%rank /= mpi%rank) then
    print *, 'make_child_patch: mpi%rank, parent%rank =', mpi%rank, parent%rank
    call io%abort ('make_child_patch: wrong MPI rank')
  end if
  call parent%log ('refine', 3)
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
  ! Reset child
  !-----------------------------------------------------------------------------
  call child%reset ()
  child%mv = child%nv                                ! extra output (FIXME?)
  child%size = child_size                            ! new size
  child%position = pos                               ! patch position
  child%centre_nat = pos                             ! patch position
  child%llc_cart = child%position - &                ! mid point position
    merge(0.5,0.0,child%n>1) * child%size            ! lower left corner
  child%llc_nat = child%llc_cart                     ! FIXME!: this form only applies for Cartesian coords.
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
  call child%clear(bits%boundary)
  call child%clear(bits%ready)
  call child%clear(bits%frozen)
  !-----------------------------------------------------------------------------
  ! Set the busy flag, to prevent scheduling until the task is ready to update
  !-----------------------------------------------------------------------------
  call child%set(bits%busy)
  if (verbose > 1) then
    call bits%print (stdlog, pad_trace(1:indent)//'refine', child%status)
  end if
  !-----------------------------------------------------------------------------
  ! Append links to the new patch in the patch list, sorted by quality
  ! Neighbor lists are constructed at the end of `split_patch`.
  !-----------------------------------------------------------------------------
  task => child                                      ! class(task_t) argument
  call self%task_list%prepare (task)                          ! allocate & init task%link
  child%link%parent => parent%link                   ! used in remove_patch
  child%parent => parent                             ! task relation
  !-----------------------------------------------------------------------------
  ! Initialise the new patch
  !-----------------------------------------------------------------------------
  child%level = parent%level + 1                     ! must be set before init!
  child%pre_initialized = .false.                    ! force pre_init
  call child%init()                                  ! patch, for a new id
  call child%init_bdries                             ! BCs
  child%boundaries%bits = 0                          ! reset bdry bits
  child%level = parent%level + 1                     ! override patch_t%init
  child%time  = parent%time                          ! same time
  child%dtime = parent%dtime/self%ratio              ! for info only
  child%dt    = child%dtime                          ! copy to all
  child%t     = child%time                           ! copy to all
!print *, child%id, 'created', omp%thread, child%level, child%time
  child%iit   = 1                                    ! initial slots
  child%it    = 1                                    ! slot for current values
  child%new   = 2                                    ! slot for new values
  child%istep = 0                                    ! for diagnostics
  child%n_needed = 1                                 ! reset to initial
  child%check_refine_next = check_interval           ! next check
  child%iit(child%nt) = child%new                    ! slot for new values
  child%parentid = parent%id                         ! for debugging
  child%iout     = parent%iout                       ! parent's output index
  child%out_next = parent%out_next                   ! next output time
  child%dnload_time = -1.0                           ! for duplicate call test
  child%check_refine_next = 0                        ! allow recursive refine
  !-----------------------------------------------------------------------------
  ! Copy the %it slot to %new
  !-----------------------------------------------------------------------------
  !call child%pre_update ()
  !-----------------------------------------------------------------------------
  ! Check consistency
  !-----------------------------------------------------------------------------
  unit = merge (stdout, stdlog, verbose > 2)
  if (any(pos /= child%position  )) write (unit,*) 'ERROR: position'
  if (any(pos /= child%mesh%p    )) write (unit,*) 'ERROR: mesh%p'
  if (any(pos /= child%centre_nat)) write (unit,*) 'ERROR: centre_nat'
  if (any(pos /= child%mesh%centre_nat)) write (unit,*) &
   'ERROR: mesh%centre_nat'
  child%position        = pos
  child%centre_nat      = pos
  child%mesh%p          = pos
  child%mesh%centre_nat = pos
  !-----------------------------------------------------------------------------
  ! Increment I/O write counter, to avoid repeated writing of metadata
  !-----------------------------------------------------------------------------
  call data_io%update_counters (child, +1)
  !-----------------------------------------------------------------------------
  ! Download values from parent and neighbors, except at time = 0, to allow
  ! recursive refinement there. If time = 0, patch should have been filled
  ! with values during `child%init()` call.
  !-----------------------------------------------------------------------------
  !$omp atomic update
  parent%n_needed = parent%n_needed + 2
  !-----------------------------------------------------------------------------
  ! Make parent a neighbor of the new child patch, to permit interpolation of
  ! values from parent to child.
  !-----------------------------------------------------------------------------
  allocate (child_nbor)
  child_nbor%link => parent%link
  child_nbor%task => parent
  child%link%nbor => child_nbor
  if (parent%time > 0.0) then
    !---------------------------------------------------------------------------
    omp%thread_task = child%id
    call child%link%update (all_tgt=.true., all_src=.true.)
    omp%thread_task = parent%id
    !...........................................................................
    if (do_dump) then
      call child%dump (3, ad_hoc=.true.)
      call parent%dump (4, ad_hoc=.true.)
    end if
  else
    if (do_dump) then
      call child%dump (5, ad_hoc=.true.)
      call parent%dump (6, ad_hoc=.true.)
    end if
  end if
  deallocate (child%link%nbor)
  !-----------------------------------------------------------------------------
  ! Add the task to the task list, locking the queue briefly, to avoid that the
  ! task is discovered and testet, before the no_support and ready bits are set
  !-----------------------------------------------------------------------------
  call self%task_list%qlock%set ('add_by_level')
  call self%task_list%add_by_level (child%link, count=.true.)
  call self%task_list%qlock%unset ('add_by_level')
  !-----------------------------------------------------------------------------
  if (verbose > 1 .or. child%id == io%id_debug) then
    if (child%id == io%id_debug) &
      unit = stdout
    l = parent%index_only (child%position + child%mesh%lf)
    u = parent%index_only (child%position + child%mesh%uf)
    l = max(l,parent%mesh%li)
    u = min(u,parent%mesh%ui)
    write(unit,1) &
     'make_child_patch: parent id,level,time,size,pos,min,max(d)=', &
      parent%id, parent%level, parent%time, parent%size, parent%position, &
      minval(child%mem(l(1):u(1),l(2):u(2),l(3):u(3),child%idx%d,child%it,1)), &
      maxval(child%mem(l(1):u(1),l(2):u(2),l(3):u(3),child%idx%d,child%it,1))
    write(unit,2) &
     '                      l, u, contains(child) =', &
      l, u, parent%contains(child)
    1 format(a,i6,i3,f9.6,2(1x,3f8.4),1p,2e12.2)
    2 format(a,2(1x,3i3),l2)
    l(:) = parent%mesh(:)%li
    u(:) = parent%mesh(:)%ui
    write(unit,1) &
     '                   child id,level,time,size,pos,min,max(d)=', &
      child%id, child%level, child%time, child%size, child%position , &
      minval(child%mem(l(1):u(1),l(2):u(2),l(3):u(3),child%idx%d,child%it,1)), &
      maxval(child%mem(l(1):u(1),l(2):u(2),l(3):u(3),child%idx%d,child%it,1))
    flush (unit)
  end if
  if (self%task_list%verbose > 1) then
    write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
      wallclock(), omp%thread, 'created', child%id, child%time, parent%id
  end if
  !-----------------------------------------------------------------------------
  call trace%end (itimer, detailed_timer=detailed_timer)
END SUBROUTINE make_child_patch

!===============================================================================
!> Refine / derefine on ratio of Jeans length to cell size
!===============================================================================
FUNCTION check_Jeans (self, patch) RESULT(refine_it)
  class(refine_t):: self
  class(solver_t), pointer:: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: Jeans_num, dmin, pgmin, ds, c
  real, dimension(:,:,:), pointer, contiguous:: pg, test
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e
  logical:: masked
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_jeans')
  refine_it = -1
  l = 1
  u = patch%gn
  d => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%it,1)
  dmin = minval(d)
  call io%assert (associated(patch), '%patch not associated in check_Jeans')
  allocate(pg(patch%gn(1),patch%gn(2),patch%gn(3)))
  !call patch%gas_pressure_sub (pg)
  pg = patch%gas_pressure()
  pgmin = minval(pg)
  if (pgmin <= 0.0) then
    imin = minloc(pg)
    write (stderr,*) &
      'check_Jeans ERROR: id, pgmin, imin =', patch%id, pgmin, imin
    flush (stderr)
    refine_it = 0
  else if (dmin <= 0.0) then
    imin = minloc(d)
    write (stderr,*) &
      'check_Jeans ERROR: id, dmin, imin =', patch%id, dmin, imin
    flush (stderr)
    refine_it = 0
  else
    !-------------------------------------------------------------------------
    ! Evaluate a per-cell criterion, which adds to patch%frefine
    !-------------------------------------------------------------------------
    allocate (test(patch%gn(1),patch%gn(2),patch%gn(3)))
    ds = maxval(patch%mesh%d, mask=patch%mesh%n > 1)
    c =  ds*sqrt(scaling%grav/(math%pi*eos%gamma))
    test = min_jeans*c*d/sqrt(pg)
    if (verbose > 1) then
      call patch%stats (d   , '   d', stdlog)
      call patch%stats (pg  , '  pg', stdlog)
      call patch%stats (test, 'test', stdlog)
    end if
    !---------------------------------------------------------------------------
    ! Combine with previous tests and evaluate
    !---------------------------------------------------------------------------
    call self%evaluate ('Jeans', patch, test, refine_it)
    deallocate (test)
  end if
  deallocate (pg)
  call trace%end()
END FUNCTION check_Jeans

!===============================================================================
!> Refine / derefine on vorticity (velocity difference per cell)
!===============================================================================
FUNCTION check_compress (self, patch) RESULT(refine_it)
  class(refine_t):: self
  class(solver_t):: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: c, ds, cut_off
  real, dimension(:,:,:), pointer, contiguous:: test
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_compress', itimer=itimer, detailed_timer=detailed_timer)
  allocate (test(patch%gn(1),patch%gn(2),patch%gn(3)))
  call patch%compression_magnitude (test)
  !---------------------------------------------------------------------------
  ! Evaluate a normalized criterion, which becomes 2.0 at max_compress, and
  ! is reduced with a factor 0.5 for each level at or above max_compress_level,
  ! allowing only maintenance of refinement at max_compress_level, if / when
  ! this is the only criterion
  !---------------------------------------------------------------------------
  ds = maxval(patch%mesh%d, mask=patch%mesh%n > 1)
  cut_off = min(1.0,0.5**(1+patch%level-max_compress_level))
  c = ds/max_compress*cut_off
  test = test*c
  !---------------------------------------------------------------------------
  ! Combine with previous tests and evaluate
  !---------------------------------------------------------------------------
  call self%evaluate ('compress', patch, test, refine_it)
  deallocate (test)
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION check_compress

!===============================================================================
!> Refine on support
!===============================================================================
FUNCTION check_vorticity (self, patch) RESULT(refine_it)
  class(refine_t):: self
  class(solver_t):: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: c, ds
  real, dimension(:,:,:), pointer, contiguous:: test
  logical:: masked
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_vorticity', itimer=itimer, detailed_timer=detailed_timer)
  !---------------------------------------------------------------------------
  ! Compute vorticity
  !---------------------------------------------------------------------------
  allocate (test(patch%gn(1),patch%gn(2),patch%gn(3)))
  call patch%vorticity_magnitude (test)
  !---------------------------------------------------------------------------
  ! Evaluate a per-cell criterion, which adds to patch%frefine
  !---------------------------------------------------------------------------
  ds = maxval(patch%mesh%d, mask=patch%mesh%n > 1)
  c = ds/max_vorticity
  test = test*c
  !---------------------------------------------------------------------------
  ! Combine with previous tests and evaluate
  !---------------------------------------------------------------------------
  call self%evaluate ('vorticity', patch, test, refine_it)
  deallocate (test)
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION check_vorticity

!===============================================================================
!> Refine on shocks, with threshold set by tolerance parameter `tol_shock`.
!===============================================================================
FUNCTION shock_and_cd_detector (self, patch) RESULT (refine)
  USE mesh_mod
  class(refine_t):: self
  class(solver_t):: patch
  integer:: refine
  !.............................................................................
  class(mesh_t), pointer:: mx, my, mz
  integer:: ix, iy, iz, ione=0, jone=0, kone=0
  real(kind=KindScalarVar):: dp1, dp2, dp3, dv1, dv2, dv3, dd1, dd2, dd3
  real, parameter:: small = 1.0e-10
  real(kind=KindScalarVar), dimension(:,:,:), pointer, contiguous:: pg
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer, contiguous:: v
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
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  e => patch%mem(:,:,:,patch%idx%e,patch%it,1)
  pg = patch%gas_pressure()
  v  = patch%gas_velocity()
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
      write (io_unit%log,*) 'shock_and_cd_detector: refine', patch%id, patch%level, self%levelmax
    else if (refine < 0) then
      write (io_unit%log,*) 'shock_and_cd_detector: derefine', patch%id, patch%level, self%levelmin
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
  class(solver_t) :: patch
  integer:: refine
  !.............................................................................
  class(mesh_t), pointer:: mx, my, mz
  integer:: i, j, k, ione=0, jone=0, kone=0, im1, jm1, km1, ip1, jp1, kp1, n
  real, parameter:: small=1e-10
  real(kind=KindScalarVar):: dqmax
  real(kind=KindScalarVar), dimension(:,:,:), pointer, contiguous:: q, dq1, dq2, dq3
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
    if (trim(grad_var(n)) == 'p')  q = patch%gas_pressure()
    if (trim(grad_var(n)) == 'v1') q = patch%gas_velocity(1)
    if (trim(grad_var(n)) == 'v2') q = patch%gas_velocity(2)
    if (trim(grad_var(n)) == 'v3') q = patch%gas_velocity(3)
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
        write (io_unit%log,*) 'gradient_detector: refine ', grad_var(n), patch%id
      else if (refine < 0) then
        write (io_unit%log,*) 'gradient_detector: derefine ', grad_var(n), patch%id
      end if
    end if
  end do
  deallocate (dq1, dq2, dq3, q)
  call trace%end()
END FUNCTION gradient_detector

!===============================================================================
!> Add test and evaluate, renormalizing with (by default) a factor 2.5,
!> allowing derefine below 0.0 and refininig for values above 1.5, on a scale
!> starting at -1.0
!===============================================================================
SUBROUTINE evaluate (self, label, patch, test, refine_it)
  class(refine_t):: self
  character(len=*):: label
  class(solver_t):: patch
  real:: test(:,:,:)
  integer:: refine_it, unit
  !...........................................................................
  real:: max_crit, max_test, norm
  integer:: l(3), u(3)
  !---------------------------------------------------------------------------
  ! Add to previous values of patch%frefine, with a normalization such that
  ! test > 1.0 triggers refinement, independent of if quadrature is used when
  ! combining crieria.
  !---------------------------------------------------------------------------
  call trace%begin ('refine_t%evaluate')
  l = patch%mesh%lo
  u = patch%mesh%uo
  norm = (1.+derefine_limit)/n_criteria
  max_test = maxval(test(l(1):u(1),l(2):u(2),l(3):u(3)))
  call patch%lock%set ('refine_t%evaluate')
  if (quadrature) then
    patch%frefine = patch%frefine + norm*test**2
  else
    patch%frefine = patch%frefine + norm*test
  end if
  call patch%lock%unset ('refine_t%evaluate')
  !---------------------------------------------------------------------------
  ! As refine remains <= 0.0 derefinement is allowed, and if it 
  ! exceeds 1.0, refinement is needed
  !---------------------------------------------------------------------------
  max_crit = maxval(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3)))
  refine_it = 0
  if (max_crit <= 0.0) refine_it = -1
  if (max_crit >= 1.0) refine_it = +1
  !---------------------------------------------------------------------------
  if (verbose > 2 .or. &
     (verbose > 1 .and. patch%level < self%levelmax .and. refine_it > 0) .or. &
     (verbose > 1 .and. patch%level > self%levelmin .and. refine_it < 0)) then
    unit = merge (stdout, stdlog, verbose > 2)
    write (unit,'(a22,3x,a,f12.6,i6,2f7.3,i4,2x,a1)') &
      label//' criterion: ', 'time, id, max_test, max_crit, refine', &
      patch%time, patch%id, max_test, max_crit, refine_it, merge ('R',' ',refine_it==1)
  end if
  if (patch%level <= self%levelmin) refine_it = max(0,refine_it)
  if (patch%level >= self%levelmax) refine_it = min(0,refine_it)
  call trace%end ()
END SUBROUTINE evaluate

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
  logical :: refine_it, refine_this_one
  integer :: ix, iy, iz
  real :: x, y, z
  real, dimension(3) :: p_llc, p_urc, b_llc, b_urc
  !.............................................................................
  call trace%begin('refine_t%refine_on_box')
  refine_it = .false.
  refine=0
  if (patch%level >= region_maxlevel) then
    call trace%end ()
    return
  end if
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
  call patch%lock%set ('refine_t%refine_region')
  do iz=m3%li,m3%ui
    z = m3%p + m3%r(iz)
    do iy=m2%li,m2%ui
      y = m2%p + m2%r(iy)
      do ix=m1%li,m1%ui
        x = m1%p + m1%r(ix)
        refine_this_one = (x >= b_llc(1) .and. x <= b_urc(1) .and. &
                           y >= b_llc(2) .and. y <= b_urc(2) .and. &
                           z >= b_llc(3) .and. z <= b_urc(3) )
        refine_it = refine_it .or. refine_this_one
        if (refine_this_one) patch%frefine(ix,iy,iz) = +1
      end do
    end do
  end do
  call patch%lock%unset ('refine_t%refine_region')
  refine = merge(1,0,refine_it)
  call trace%end()
END FUNCTION refine_region

!===============================================================================
!> Check if refinement may be needed for this level L patch to support 
!> overlapping level L+2 patches.  If there are any level L+2 nbors, check
!> if the level L+1 patches needed for support already exist.
!===============================================================================
FUNCTION check_support (self, link) RESULT (refine)
  class(refine_t):: self
  class(link_t), pointer:: link
  integer:: refine
  !.............................................................................
  class(task_t), pointer:: patch, nbpatch
  class(nbor_t), pointer:: nbor
  integer:: i1, i2, i3, ratio(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_support')
  call link%lock%set ('check_support')
  patch => link%task
  select type (patch)
  class is (solver_t)
    ratio = merge(self%ratio, 1, patch%n>1)
    refine = -1
    !---------------------------------------------------------------------------
    ! Loop over nbor tasks and investigate level L+2 cases
    !---------------------------------------------------------------------------
    nbor => link%nbor
    do while (associated(nbor))
      nbpatch => nbor%task
      if (nbpatch%level == patch%level+2) then
        !-----------------------------------------------------------------------
        ! The mere existence of a level L+2 task prevents derefinement
        !-----------------------------------------------------------------------
        refine = max(0,refine)
        !-----------------------------------------------------------------------
        ! Tests of the possible sub-patch positions return 0 if there already
        ! exists one or more sub-patches, and +1 if a new one is needed
        !-----------------------------------------------------------------------
        do i3=1,ratio(3)
        do i2=1,ratio(2)
        do i1=1,ratio(1)
          refine = max (refine, check_support2 (self, patch, [i1,i2,i3]))
        end do
        end do
        end do
      end if
      nbor => nbor%next
    end do
  end select
  call link%lock%unset ('check_support')
  call trace%end ()
END FUNCTION check_support

!===============================================================================
!> Check if refinement is needed to support some level L+2 task that overlaps
!> with index position ii(3) relative to the level L link task
!===============================================================================
FUNCTION check_support2 (self, patch, ii) RESULT (refine)
  class(refine_t):: self
  class(solver_t), pointer:: patch
  integer:: ii(3), refine
  !.............................................................................
  class(task_t), pointer:: nbpatch
  class(link_t), pointer:: nblink
  class(nbor_t), pointer:: nbor
  integer:: ratio(3), key(4)
  real(8), dimension(3):: size, position, distance
  type(region_t):: roa
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_support2')
  ratio = merge(self%ratio, 1, patch%n>1)
  size = patch%size/ratio
  !-----------------------------------------------------------------------------
  ! Assume until disproven that based on support the patch allows derefinement
  !-----------------------------------------------------------------------------
  refine = -1
  nbor => patch%link%nbor
  do while (associated(nbor))
    nbpatch => nbor%task
    if (nbpatch%level == patch%level+2) then
      select type (nbpatch)
      class is (patch_t)
        !---------------------------------------------------------------------
        ! Compoute region-of-authority for the sub-patch and check overlap
        !---------------------------------------------------------------------
        position = patch%position + (ii-(1+ratio)/2.0)*size
        roa%l = position - size/2d0
        roa%u = position + size/2d0
        distance = nbpatch%distance_to (patch)
        if (nbpatch%interest%overlaps (roa, distance)) then
          !-------------------------------------------------------------------
          ! There is overlap, so check if a patch already exists, and
          ! signal that one is neeeded if not
          !-------------------------------------------------------------------
          key = patch%hash_key (position, size, patch%level+1)
          call self%task_list%hash_table%get (key, nblink)
          if (associated(nblink)) then
            refine = 0
          else
            refine = 1
            exit
          end if
        end if
      end select
    end if
    nbor => nbor%next
  end do
  call trace%end ()
END FUNCTION check_support2

END MODULE refine_mod
