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
  USE os_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE mpi_mod
  USE mpi_global_mod
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
  USE dll_mod
  USE experiment_mod
  USE halo_mod
  USE guards_mod
  USE patch_lock_mod
  implicit none
  private
  intrinsic:: log
  type:: refine_t
    logical:: on                              ! allow external access
    integer:: levelmin=0                      ! Min level of refinment
    integer:: levelmax=0                      ! Max level of refinment
    integer:: ratio=2                         ! default
    real:: safe_ratio=2.1                     ! safe derefine ratio
    class(list_t), pointer:: task_list => null()  ! global list
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
    procedure:: custom_refine
    procedure:: refine_region
    procedure:: gradient_detector
    procedure:: check_derefine_support
    procedure:: create_support
    procedure:: create_patch
  end type
  !
  integer, parameter:: max_levels=40          ! max number of levels
  integer:: check_interval=3                  ! Number of steps betwen checks
  real:: min_dx=0.0, max_dx=1e32              ! Min and max grid size
  real:: refine_from_t=0.0                    ! Time limits for occurence of refinement
  real:: refine_until_t=1e30                  ! Time limits for occurence of refinement
  real:: min_jeans=0.0                        ! Jeans number above which a patch is flagged for derefinement
  real:: max_shock=0.0, max_contact=0.0       ! Tolerance values for shock/discontinuity detector refinement
  real:: min_shock=0.0, min_contact=0.0       ! Tolerance values for shock/discontinuity detector derefinement
  real:: max_compress=0.0                     ! Tolerance value for compression (negative divergence)
  real:: max_custom=0.0                       ! Tolerance value for custom test
  integer:: max_compress_level=-1             ! Max refinement level on compression
  real:: derefine_limit=0.4                   ! Fraction of criterion to allow derefinement
  integer, parameter:: ngradvar=10            ! Number of gradient variables
  real, dimension(3,3):: region_llc=0.0       ! Sometimes you want to refine a region (manually)
  real, dimension(3,3):: region_urc=0.0       ! upper-right-corner of region
  integer:: region_maxlevel(3)=5              ! max region refinement level
  real, dimension(ngradvar):: max_grad=0.0    ! Tolerance values for gradient detector refinement
  real, dimension(ngradvar):: min_grad=0.0    ! Tolerance values for gradient detector derefinement
  real, allocatable:: max_vorticity(:)        ! Tolerance value for vorticity
  character(len=32), dimension(ngradvar):: grad_var=""
  logical:: full_grad=.false.                 ! expensive variant
  logical:: force_cubical=.false.             ! Force cubical patches
  logical:: detailed_timer=.false.            ! Lump all AMR together when false
  logical:: do_dump=.false.                   ! dump parent & child
  integer:: n_locks=1                         ! number of locks to use
  integer:: verbose=0                         ! Noise level
  integer:: n_criteria=0                      ! Number of active criteria
  integer:: wait_steps=0                      ! Steps to wait before refinement
  integer, save:: levelmin=0, levelmax=0      ! criteria
  logical:: quadrature=.false.                ! Combine criteria in quadrature
  type(refine_t), public:: refine
CONTAINS

!===============================================================================
!> Initialize refinement parameters (called from task_list_t%init_levels)
!===============================================================================
SUBROUTINE init (self)
  class(refine_t):: self
  integer:: iostat
  integer, save:: ratio=2           ! How does a patch split when it is refined
  logical:: on=.false.              ! Flag to turn the entire process on or off
  integer:: i
  logical, save:: first_time = .true.
  character(len=120):: id = &
    '$Id$'
  namelist /refine_params/ verbose, on, levelmin, levelmax, check_interval, &
    ratio, refine_from_t, refine_until_t, min_dx, max_dx, derefine_limit, &
    min_jeans, max_shock, min_shock, max_contact, min_contact, max_grad, &
    min_grad, max_vorticity, max_compress, max_compress_level, quadrature, &
    grad_var, region_llc, region_urc, region_maxlevel, force_cubical, n_locks, &
    wait_steps, detailed_timer, do_dump, full_grad, max_custom
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%init')
  call io%header ('refine_t%init:')
  call trace%print_id (id)
  call io%lock%set ('refine_params')
  if (first_time) then
    first_time=.false.
    levelmin = shared%levelmin
    levelmax = shared%levelmax
    allocate (max_vorticity(1:max_levels))
    max_vorticity = 0.0
    rewind (io%input)
    read (io%input, refine_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, refine_params)
    end if
    call io%namelist_warning ('refine_params', iostat)
    if (max_compress_level == -1) &
      max_compress_level = levelmax
    !---------------------------------------------------------------------------
    ! If AMR is active, we must use check_filled, since nbor levels are in
    ! decreasing order, to save computing time.
    !---------------------------------------------------------------------------
    if (on) then
      omp_lock%links = .true.
      if (max_vorticity(2)==0) then
        do i=2,max_levels
          max_vorticity(i) = max_vorticity(1)/2.**(i*0.38)
        end do
      end if
      call reallocate (max_vorticity, levelmin, levelmax)
    else
      omp_lock%links = .false.
    end if
    write (io%namelist, refine_params)
    shared%refine_on = on
    shared%f_zoom = 1.0/ratio
    shared%refine_ratio = ratio
    if (min_jeans         > 0.0)  n_criteria = n_criteria+1
    if (max_compress      > 0.0)  n_criteria = n_criteria+1
    if (any(max_vorticity > 0.0)) n_criteria = n_criteria+1
    if (max_shock         > 0.0)  n_criteria = n_criteria+1
    if (max_contact       > 0.0)  n_criteria = n_criteria+1
    if (maxval(max_grad)  > 0.0)  n_criteria = n_criteria+1
    if (any(region_llc/=0) .or.any(region_urc/=0)) n_criteria = n_criteria+1
    write (io%output,'(a,i4)') ' number of AMR criteria:', n_criteria
    !---------------------------------------------------------------------------
    ! Initialize a global array for patch creation locking, with room for
    ! at most on average 8 pending patch creations per rank
    !---------------------------------------------------------------------------
    call patch_lock%init (8)
    !---------------------------------------------------------------------------
    ! Initialize the lock in the garbage collector (shared by link_mod)
    !---------------------------------------------------------------------------
    call garbage%init
  end if
  call io%lock%unset ('refine_params')
  self%levelmin = levelmin
  self%levelmax = levelmax
  shared%levelmin = levelmin
  shared%levelmax = levelmax
  self%ratio = ratio
  self%on = on
  call trace%end()
contains
  !=============================================================================
  subroutine reallocate (f, min, max)
    real, allocatable:: f(:), tmp(:)
    integer:: min, max
    !---------------------------------------------------------------------------
    allocate (tmp(levelmin:levelmax))
    tmp = max_vorticity(1:levelmax-levelmin+1)
    deallocate (max_vorticity)
    allocate (max_vorticity(levelmin:levelmax))
    max_vorticity = tmp
    deallocate (tmp)
  end subroutine reallocate
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
  class(experiment_t), pointer:: patch
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
  task => link%task
  if (task%time > refine_until_t) then
    self%on = .false.
    shared%refine_on = .false.
  end if
  if (.not.self%on .or. &
      task%istep < wait_steps) then
        return
  end if
  if (verbose > 2) then
    if (task%is_set (bits%init_nbors)) &
      write(stdlog,'(i6,f12.6,4x,a,3x,2l1)') &
        task%id, link%task%time, 'has bits%init_nbors', &
        task%is_set (bits%boundary), &
        task%is_set (bits%virtual)
  end if
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%check', itimer=itimer)
  self%task_list => tasklist
  id0 = task%id
  select type (task)
  class is (experiment_t)
  patch => task
  end select
  id = patch%id
  level = patch%level
  !-----------------------------------------------------------------------------
  dxmin = minval(patch%ds)
  dxmax = maxval(patch%ds)
  if (verbose > 1) then
    write (stdlog,'(a,i6,":",i2.2,3i3,2i7,l3,1p,2(2x,2e11.2))') &
      'refine_t%check: id,istep,next,V,dx_min,dx_max =', &
      id, patch%key, patch%istep, patch%check_refine_next, &
      patch%is_set(bits%virtual), min_dx, dxmin, max_dx, dxmax
  end if
  !-----------------------------------------------------------------------------
  ! Filter out patches that do not match the necessary criteria
  !-----------------------------------------------------------------------------
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
    call os%allocate (patch%test, patch%gn, 'test')
    refine = self%refine_needed (link)
    !---------------------------------------------------------------------------
    ! Force refinment up to levelmin
    !---------------------------------------------------------------------------
    if (level < levelmin) then
      refine = 1
      patch%frefine = 2.0
    end if
    !---------------------------------------------------------------------------
    ! Apply refinement criteria on potential child patch parts, but only if
    ! either refinement has been indicated, or if derefinement is allowed
    !---------------------------------------------------------------------------
    if (refine /= 0) then
      call self%selective_refine (patch, refine, was_refined, was_derefined)
    end if
    call os%deallocate (patch%test, 'test')
    if (verbose > 2) &
      write (stdlog,*) 'deallocated irefine in refine', mpi%rank, id
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
  class(experiment_t), pointer :: patch
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
  class is (experiment_t)
  patch => task
  refine = -1
  call patch%support_arrays ()
  patch%frefine = -derefine_limit
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
  !-----------------------------------------------------------------------------
  ! Custom experiment refinement
  !-----------------------------------------------------------------------------
  if (max_custom > 0) refine = max(refine, self%custom_refine(patch))
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, &
    minval(patch%frefine), maxval(patch%frefine), 'extras'
    1 format(f12.6,a,"refine_needed: id,refine,min,max,type=", & 
             i6,":",i2.2,i3,1p,2e11.2,2x,a) 
  end select
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
  if (any(max_vorticity > 0.0)) &
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
  !---------------------------------------------------------------------------
  ! Check if support constraints allow derefinement
  !---------------------------------------------------------------------------
  refine = max (refine, self%check_derefine_support (link))
  if (verbose > 1) write (stdlog,1) &
    wallclock(), pad_trace(1:indent), id, level, refine, &
    minval(patch%frefine), maxval(patch%frefine), 'derefine support'
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
SUBROUTINE selective_refine (self, patch, refine, &
                             was_refined, was_derefined)
  class(refine_t)          :: self
  class(experiment_t), pointer :: patch
  integer                  :: refine
  logical                  :: was_refined, was_derefined
  !.............................................................................
  class(experiment_t), pointer :: child
  class(link_t),   pointer :: link, old_link, next
  class(nbor_t),   pointer :: nbor
  class(task_t),   pointer :: task
  class(dll_node_t), pointer :: node
    type(dll_t)            :: dll_list
  integer                  :: i, j, k, id, level, id1
  integer, dimension(3)    :: l, u, n
  real(8), dimension(3)    :: size, pos
  character(24)            :: label
  real                     :: frefine_max
  logical                  :: already_refined, refined
  integer                  :: refine_it, ratio_per_dim(3), key(4), n_sub
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
  call dll_list%init ('refine%dll_list')
  if (patch%level < self%levelmax) then
    l = patch%mesh%li
    u = patch%mesh%ui
    if (maxval(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3))) >= 1.0) then
                                        call patch%lock%set ('selective_refine')
      do k=1,ratio_per_dim(3)
        do j=1,ratio_per_dim(2)
          do i=1,ratio_per_dim(1)
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
            key = patch%hash_key (level=patch%level+1, size=size, position=pos) !has_key overloadable until experiment_t
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
              refine = max(0,refine)
            else
              if (refine_it > 0) then
                was_refined = .true.
                refined = .true.
                !---------------------------------------------------------------
                ! Create a refined patch here
                !---------------------------------------------------------------
                if (verbose > 1) &
                  write (stdlog,'(f12.6,a,i6,1p2g16.5,2x,3(i4,":",i2.2))') &
                    wallclock(), pad_trace(1:indent)// &
                    'refine', patch%id, &
                    patch%fmaxval(patch%mem(:,:,:,patch%idx%d,patch%it,1)), &
                    maxval(patch%mem(l(1):u(1),l(2):u(2),l(3):u(3), patch%idx%d, &
                    patch%it,1)), l(1),u(1),l(2),u(2),l(3),u(3)
                call patch%log ('refine')
                call make_child_patch (self, patch, pos, size, child)
                call self%create_support (child, patch, dll_list)
                !---------------------------------------------------------------
                ! Allocate a tempory task link and append it to the child list,
                ! making sure not to change the task%link, which is needed below
                !---------------------------------------------------------------
                allocate (node)
                node%car => child%link
                call dll_list%append (node)
                !---------------------------------------------------------------
                if (verbose > 0) &
                  write (stdlog,'(f12.6,2(a,i6,":",i2.2,2x),a,2l1,1p,e14.5)')  &
                    wallclock(), pad_trace(1:indent)// &
                    'refine:'   ,  id, level, &
                    'to'        ,  child%id, child%level, &
                    'BV:', child%is_set(bits%boundary), child%is_set(bits%virtual), &
                    child%time
                flush (stdlog)
                !---------------------------------------------------------------
              end if
            end if
          end do
        end do
      end do
                                      call patch%lock%unset ('selective_refine')
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Then, consider if the whole patch should be derefined.  This requires that
  ! it has no higher level patch inside (would trigger patch%refine >= 0 there),
  ! and secondly that frefine(:,:,:) indicates that derefine would be OK.
  ! Also check the nbor list for level L+1 nbors
  !-----------------------------------------------------------------------------
  if (refine == -1                .and. &               ! derefine possible
      patch%level > self%levelmin .and. &               ! relevant?
      patch%check_refine_next > 0 .and. &               ! not immediately!
      patch%istep >= 4            .and. &               ! must take 4 steps
      patch%is_clear (bits%init_nbors) .and. &          ! no pending init_nbors
      .not.refined) then                                ! not if refined!
    l = patch%mesh%li
    u = patch%mesh%ui
    frefine_max = maxval(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3)))
    if (verbose == 2) &
      write (stdlog,*) 'derefine: frefine_max=', frefine_max 
    if (frefine_max < 0.0) then
      !-------------------------------------------------------------------------
      ! Check if there are any level L+1 nbors; these would be outside this
      ! patch, and hence not detected by the conditions above
      !-------------------------------------------------------------------------
      n_sub = 0
      nbor => patch%link%nbor
      do while (associated(nbor))
        if (nbor%task%level == patch%level+1) then
          n_sub = n_sub+1
        end if
        nbor => nbor%next
      end do
      if (verbose == 1) &
        write (stdlog,*) 'derefine: frefine_max,n_sub=', frefine_max, n_sub 
      !-------------------------------------------------------------------------
      ! If task is to be removed, call a comprehensive remove procedure, which
      ! includes calling set_init_nbors(), queueing all nbors for immediate
      ! update of their nbor lists.
      !-------------------------------------------------------------------------
      if (n_sub == 0) then
        !.......................................................................
        if (verbose > 0) then
          write (stdlog,'(f12.6,2x,a,i6,":",i2.2,16x,a,2l1,2x,a,1p,e14.5,i3)') &
            wallclock(), &
            'derefine:' ,  id, level, 'BV:', &
            patch%is_set(bits%boundary), patch%is_set(bits%virtual), &
            'time,n_need:', patch%time, patch%n_needed
        end if
        !.......................................................................
        if (io%task_logging > 0) then
          write (label,"('derefine',2x,2a1)") &
            merge('B','b',patch%is_set(bits%boundary)), &
            merge('V','v',patch%is_set(bits%virtual))
          call patch%log (trim(label), 1)
        end if
        call data_io%update_counters (patch, -1)
        call self%task_list%remove_and_reset (link)
        was_derefined = .true.
      end if
    end if
  end if
  call self%post_process (patch, dll_list)
  call trace%end (itimer, detailed_timer=detailed_timer)
END SUBROUTINE selective_refine

!============================================ ==================================
!> Post process the list of new tasks created by refine%check()
!===============================================================================
SUBROUTINE post_process (self, patch, dll_list)
  class(refine_t)              :: self
  class(experiment_t), pointer :: patch
  type(dll_t)                  :: dll_list
  !.............................................................................
  class(link_t),       pointer :: link
  class(task_t),       pointer :: task
  class(dll_node_t),   pointer :: node, next
  class(*),            pointer :: ptr
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%post_process')
  !-----------------------------------------------------------------------------
  ! The new tasks in the dll_list need to be 
  ! 1) added to the task list (done above)
  ! 2) have new nbor lists, which include the new siblings
  ! 3) have bits%no_support set on or off
  ! 4) be added to the queue if ready to update
  ! Since 4) depends on all new tasks having done 2) and 3) it is not done
  ! here rather than inside add_new_link
  !-----------------------------------------------------------------------------
                          !call self%task_list%lock%set ('refine_t%post_process')
  call patch%log ('post_process', 2)
  if (associated(dll_list%head)) then
    !---------------------------------------------------------------------------
    ! Do the standard set of actions after adding new tasks
    !---------------------------------------------------------------------------
    node => dll_list%tail
    do while (associated(node))
      ptr => node%car
      select type (ptr)
      class is (link_t)
        task => ptr%task
      end select
      select type (task)
      class is (patch_t)
        link => task%link
        call self%task_list%add_new_link (link)       ! standard actions
      end select
      node => node%prev
    end do
    omp%thread_task = patch%id
    !---------------------------------------------------------------------------
    ! Deallocate and delete the child list, after check_ready()
    !---------------------------------------------------------------------------
    node => dll_list%tail
    do while (associated(node))
      ptr => node%car
      select type (ptr)
      class is (link_t)
        task => ptr%task
      end select
      select type (task)
      class is (patch_t)
        link => task%link
        !-----------------------------------------------------------------------
        ! Clear the busy bit, now that the task is ready to be scheduled
        !-----------------------------------------------------------------------
        call link%info ('refine%post_process')        ! write to stdlog
        call task%log ('clear ignore', 2)
        if (verbose > 2) &
          call bits%print (stderr, 'before clear', task%status)
        call task%clear (bits%ignore)                 ! clear for action
        if (verbose > 2) &
          call bits%print (stderr, 'after  clear', task%status)
        call task%log ('check_ready', 2)
        call self%task_list%check_ready (link)        ! check if ready for action
      end select
      next => node%prev
      call dll_list%remove (node)
      deallocate (node)
      node => next
    end do
  end if
                        !call self%task_list%lock%unset ('refine_t%post_process')
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
  ! If the task is at levelmax it is free to derefine, while at levemin it is
  ! not.  At other levels derefinement is OK if no higher level nbors exist. 
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%check_derefine_support', itimer=itimer, &
    detailed_timer=detailed_timer)
  if (link%task%level == self%levelmax) then
    refine = -1
  else if (link%task%level == self%levelmin) then
    refine = 0
  else
    refine = -1
                                   call link%lock%set ('check_derefine_support')
    nbor => link%nbor
    do while (associated(nbor))
      if (nbor%task%level == link%task%level+1 .and. nbor%needs_me) then
        refine = 0
        exit
      end if
      nbor => nbor%next
    end do
                                 call link%lock%unset ('check_derefine_support')
    if (verbose > 4) &
      write (stdout,*) link%task%id, 'check_derefine_support =', refine
    !---------------------------------------------------------------------------
  end if
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION check_derefine_support

!===============================================================================
!> Create a new patch using an existing "parent" patch.
!> This procedure is used frequently during refinement.
!===============================================================================
SUBROUTINE make_child_patch (self, parent, pos, child_size, &
                             child)
  class(refine_t):: self
  class(experiment_t), pointer :: parent, child
  real(8), dimension(3), intent(in) :: pos, child_size
  !.............................................................................
  class(task_t), pointer :: task
  class(nbor_t), pointer :: child_nbor
  integer:: j, count, l(3), u(3), n_added, unit, key(4)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%make_child_patch', itimer=itimer, detailed_timer=detailed_timer)
  if (parent%rank /= mpi%rank) then
    print *, 'make_child_patch: mpi%rank, parent%rank =', mpi%rank, parent%rank
    call io%abort ('make_child_patch: wrong MPI rank')
  end if
  call parent%log ('make_child', 3)
  !-----------------------------------------------------------------------------
  ! Allocate a new task, 1/nsplit the size of original one, and adjust its
  ! position in each direction with "pos", in units of the new size.
  !-----------------------------------------------------------------------------
  key = parent%hash_key (position=pos, size=child_size, level=parent%level+1)   ! hash key call overloadable until experiment_t
  call self%create_patch (key, parent, child, pos)
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
  !$omp atomic update
  parent%n_needed = parent%n_needed + 2
  !-----------------------------------------------------------------------------
  if (verbose > 1 .or. child%id == io%id_debug) then
    if (child%id == io%id_debug) &
      unit = stdout
    l = parent%index_only (child%position + child%mesh%lf)
    u = parent%index_only (child%position + child%mesh%uf)
    l = max(l,parent%mesh%li)
    u = min(u,parent%mesh%ui)
    write(unit,1) &
     'parent id,level,time,size,pos,l,u,min,max(d)=', &
      parent%id, parent%level, parent%time, parent%size, parent%position, l,u, &
      minval(child%mem(l(1):u(1),l(2):u(2),l(3):u(3),child%idx%d,child%it,1)), &
      maxval(child%mem(l(1):u(1),l(2):u(2),l(3):u(3),child%idx%d,child%it,1))
    write(unit,2) &
     '                      l, u, contains(child) =', &
      l, u, parent%contains2(child)
    1 format(a,i6,":",i2.2,f10.6,2(1x,3f10.4),2(1x,3i3),1p,2e12.2)
    2 format(a,2(1x,3i3),l2)
    l(:) = parent%mesh(:)%li
    u(:) = parent%mesh(:)%ui
    write(unit,1) &
     ' child id,level,time,size,pos,l,u,min,max(d)=', &
      child%id, child%level, child%time, child%size, child%position , l, u, &
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
  class(experiment_t), pointer:: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: Jeans_num, dmin, pgmin, ds, c
  real, dimension(:,:,:), pointer, contiguous:: pg
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
  call os%allocate (pg, patch%gn, 'pg')
  call patch%gas_pressure_sub (pg)
  pgmin = minval(pg)
  if (pgmin <= 0.0) then
    imin = minloc(pg)
    write (stderr,*) &
      ' check_Jeans ERROR: id, pgmin, imin =', patch%id, pgmin, imin
    flush (stderr)
    refine_it = 0
  else if (dmin <= 0.0) then
    imin = minloc(d)
    write (stderr,*) &
      ' check_Jeans ERROR: id, dmin, imin =', patch%id, dmin, imin
    flush (stderr)
    refine_it = 0
  else
    !---------------------------------------------------------------------------
    ! Evaluate a per-cell criterion, which adds to patch%frefine
    !---------------------------------------------------------------------------
    ds = maxval(patch%mesh%d, mask=patch%mesh%n > 1)
    c =  ds*sqrt(scaling%grav/(math%pi*eos%gamma))
    patch%test = min_jeans*c*d/sqrt(pg)
    if (verbose > 1) then
      call patch%stats (d   , '   d', stdlog)
      call patch%stats (pg  , '  pg', stdlog)
      call patch%stats (patch%test, 'test', stdlog)
    end if
    !---------------------------------------------------------------------------
    ! Combine with previous tests and evaluate
    !---------------------------------------------------------------------------
    call self%evaluate ('Jeans', patch, refine_it)
  end if
  call os%deallocate (pg, 'pg')
  call trace%end()
END FUNCTION check_Jeans

!===============================================================================
!> Refine / derefine on vorticity (velocity difference per cell)
!===============================================================================
FUNCTION check_compress (self, patch) RESULT(refine_it)
  class(refine_t):: self
  class(experiment_t):: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: c, ds, cut_off
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_compress', itimer=itimer, detailed_timer= &
    detailed_timer)
  call patch%compression_magnitude (patch%test)
  !-----------------------------------------------------------------------------
  ! Evaluate a normalized criterion, which becomes 2.0 at max_compress, and
  ! is reduced with a factor 0.5 for each level at or above max_compress_level,
  ! allowing only maintenance of refinement at max_compresspatch_level, if /
  ! when this is the only criterion
  !-----------------------------------------------------------------------------
  ds = maxval(patch%mesh%d, mask=patch%mesh%n > 1)
  cut_off = min(1.0,0.5**(patch%level-max_compress_level))
  c = ds/max_compress*cut_off
  patch%test = patch%test*c
  !-----------------------------------------------------------------------------
  ! Combine with previous tests and evaluate
  !-----------------------------------------------------------------------------
  call self%evaluate ('compress', patch, refine_it)
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION check_compress

!===============================================================================
!> Refine on support
!===============================================================================
FUNCTION check_vorticity (self, patch) RESULT(refine_it)
  class(refine_t):: self
  class(experiment_t):: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: c, ds
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_vorticity', itimer=itimer, detailed_timer=&
    detailed_timer)
  !-----------------------------------------------------------------------------
  ! Compute vorticity
  !-----------------------------------------------------------------------------
  patch%test = 0.0
  call patch%vorticity_magnitude (patch%test)
  !-----------------------------------------------------------------------------
  ! Evaluate a per-cell criterion, which adds to patch%frefine
  !-----------------------------------------------------------------------------
  ds = maxval(patch%mesh%d, mask=patch%mesh%n > 1)
  c = ds/max_vorticity(patch%level)
  patch%test = patch%test*c
  !-----------------------------------------------------------------------------
  ! Combine with previous tests and evaluate
  !-----------------------------------------------------------------------------
  call self%evaluate ('vorticity', patch, refine_it)
  call trace%end (itimer, detailed_timer=detailed_timer)
END FUNCTION check_vorticity

!===============================================================================
!> Refine on shocks, with threshold set by tolerance parameter `tol_shock`.
!===============================================================================
FUNCTION shock_and_cd_detector (self, patch) RESULT (refine)
  USE mesh_mod
  class(refine_t):: self
  class(experiment_t):: patch
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
  call os%allocate (pg, patch%gn, 'pg')
  call os%allocate (v, [patch%gn,3], 'v')
  ! we need the gas pressure, density and velocity
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  e => patch%mem(:,:,:,patch%idx%e,patch%it,1)
  call patch%gas_pressure_sub (pg)
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
      write (stdlog,*) 'shock_and_cd_detector: refine', patch%id, patch%level, &
        self%levelmax
    else if (refine < 0) then
      write (stdlog,*) 'shock_and_cd_detector: derefine', patch%id, patch%level,&
        self%levelmin
    end if
  end if
  !
  nullify (mx, my, mz, d)
  call os%deallocate (pg, 'pg')
  call os%deallocate (v, 'v')
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
  class(experiment_t) :: patch
  integer:: refine
  !.............................................................................
  class(mesh_t), pointer:: mx, my, mz
  integer:: i, j, k, ione=0, jone=0, kone=0, im1, jm1, km1, ip1, jp1, kp1, n, id
  real, parameter:: small=1e-10
  real, dimension(:,:,:), pointer, contiguous:: q, dq
  real, dimension(:,:,:), pointer, contiguous:: qp
  real, parameter:: rsmall = 1.0e-10
  logical:: is_unsigned
  character(len=1):: unsigned(3) = ['d', 'p', 'e']
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%gradient_detector')
  call os%allocate (q , patch%gn, 'q')
  call os%allocate (dq, patch%gn, 'dq')
  patch%test = 0.0
  id = patch%id
  mx => patch%mesh(1)
  my => patch%mesh(2)
  mz => patch%mesh(3)
  if (mx%n > 1) ione = 1
  if (my%n > 1) jone = 1
  if (mz%n > 1) kone = 1
  do n=1,size(grad_var)
    if (trim(grad_var(n)) == '') exit
    if (trim(grad_var(n)) == 'd')  q = patch%mem(:,:,:,patch%idx%d ,patch%it,1)
    if (trim(grad_var(n)) == 'e')  q = patch%mem(:,:,:,patch%idx%e ,patch%it,1)
    if (trim(grad_var(n)) == 'b1') q = patch%mem(:,:,:,patch%idx%bx,patch%it,1)
    if (trim(grad_var(n)) == 'b2') q = patch%mem(:,:,:,patch%idx%by,patch%it,1)
    if (trim(grad_var(n)) == 'b3') q = patch%mem(:,:,:,patch%idx%bz,patch%it,1)
    if (trim(grad_var(n)) == 'E')  q = patch%mem(:,:,:,patch%idx%e ,patch%it,1)&
                                     / patch%mem(:,:,:,patch%idx%d ,patch%it,1)
    if (trim(grad_var(n)) == 'S')  q = patch%mem(:,:,:,patch%idx%s ,patch%it,1)&
                                     / patch%mem(:,:,:,patch%idx%d ,patch%it,1)
    if (trim(grad_var(n)) == 'lnd') &
                               q = log(patch%mem(:,:,:,patch%idx%d ,patch%it,1))
    if (trim(grad_var(n)) == 'lnpg') then
      qp => q
      call patch%gas_pressure_sub (qp)
      q = log(q)
    end if
    is_unsigned = .false.
    do i=1,size(unsigned)
      if (grad_var(n) == unsigned(i)) then
        is_unsigned = .true.
        exit
      end if
    end do
    ! calculate the gradient in each direction
    if (full_grad) then
      do k=mz%li,mz%ui
        km1 = k - kone
        kp1 = k + kone
        do j=my%li,my%ui
          jm1 = j - jone
          jp1 = j + jone
          do i=mx%li,mx%ui
            im1 = i - ione
            ip1 = i + ione
            dq(i,j,k) = max(q(ip1,jm1,km1), q(i  ,jm1,km1), q(ip1,jm1,km1), &
                            q(ip1,j  ,km1), q(i  ,j  ,km1), q(ip1,j  ,km1), &
                            q(ip1,jp1,km1), q(i  ,jp1,km1), q(ip1,jp1,km1), &
                            q(ip1,jm1,k  ), q(i  ,jm1,k  ), q(ip1,jm1,k  ), &
                            q(ip1,j  ,k  ), q(i  ,j  ,k  ), q(ip1,j  ,k  ), &
                            q(ip1,jp1,k  ), q(i  ,jp1,k  ), q(ip1,jp1,k  ), &
                            q(ip1,jm1,kp1), q(i  ,jm1,kp1), q(ip1,jm1,kp1), &
                            q(ip1,j  ,kp1), q(i  ,j  ,kp1), q(ip1,j  ,kp1), &
                            q(ip1,jp1,kp1), q(i  ,jp1,kp1), q(ip1,jp1,kp1)) - &
                        min(q(ip1,jm1,km1), q(i  ,jm1,km1), q(ip1,jm1,km1), &
                            q(ip1,j  ,km1), q(i  ,j  ,km1), q(ip1,j  ,km1), &
                            q(ip1,jp1,km1), q(i  ,jp1,km1), q(ip1,jp1,km1), &
                            q(ip1,jm1,k  ), q(i  ,jm1,k  ), q(ip1,jm1,k  ), &
                            q(ip1,j  ,k  ), q(i  ,j  ,k  ), q(ip1,j  ,k  ), &
                            q(ip1,jp1,k  ), q(i  ,jp1,k  ), q(ip1,jp1,k  ), &
                            q(ip1,jm1,kp1), q(i  ,jm1,kp1), q(ip1,jm1,kp1), &
                            q(ip1,j  ,kp1), q(i  ,j  ,kp1), q(ip1,j  ,kp1), &
                            q(ip1,jp1,kp1), q(i  ,jp1,kp1), q(ip1,jp1,kp1))
          end do
        end do
      end do
    else
      do k=mz%li,mz%ui
        km1 = k - kone
        kp1 = k + kone
        do j=my%li,my%ui
          jm1 = j - jone
          jp1 = j + jone
          do i=mx%li,mx%ui
            im1 = i - ione
            ip1 = i + ione
            dq(i,j,k) = max(q(im1,j  ,k  ),           q(ip1,j  ,k  ), &
                            q(i  ,jm1,k  ), q(i,j,k), q(i  ,jp1,k  ), &
                            q(i  ,j  ,km1),           q(i  ,j  ,kp1)) - &
                        min(q(im1,j  ,k  ),           q(ip1,j  ,k  ), &
                            q(i  ,jm1,k  ), q(i,j,k), q(i  ,jp1,k  ), &
                            q(i  ,j  ,km1),           q(i  ,j  ,kp1))
          end do
        end do
      end do
    end if
    !---------------------------------------------------------------------------
    ! If unsigned, use relative values
    !---------------------------------------------------------------------------
    if (is_unsigned) then
      do k=mz%li,mz%ui
        do j=my%li,my%ui
          do i=mx%li,mx%ui
            dq(i,j,k) = dq(i,j,k)/abs(q(i,j,k))
          end do
        end do
      end do
    end if
    !---------------------------------------------------------------------------
    ! Combine
    !---------------------------------------------------------------------------
    do k=mz%li,mz%ui
      do j=my%li,my%ui
        do i=mx%li,mx%ui
          patch%test(i,j,k) = max(patch%test(i,j,k),dq(i,j,k)/max_grad(n))
        end do
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Combine with previous tests and evaluate
  !-----------------------------------------------------------------------------
  call self%evaluate ('gradient', patch, refine)
  call os%deallocate (q, 'q')
  call os%deallocate (dq, 'dq')
  call trace%end()
END FUNCTION gradient_detector

!===============================================================================
!> Given a box defined by its lower and upper corners, flag a patch for
!> refinement if it lies inside this box, up to a maximum level.
!> This particular criterion is primarily used for testing and debugging,
!> but can also be used to enforce enhanced resolution in a specific region
!> (aka "static refinement").
!===============================================================================
FUNCTION refine_region (self, patch) RESULT (refine)
  USE mesh_mod
  class(refine_t) :: self
  class(experiment_t), pointer :: patch
  integer:: refine
  !.............................................................................
  class(mesh_t), pointer:: m1, m2, m3, m(:)
  logical :: intersect
  logical :: refine_it, refine_this_one
  integer :: ix, iy, iz, k, unit
  real :: x, y, z, vol
  real, dimension(3) :: p_llc, p_urc, b_llc, b_urc
  !.............................................................................
  call trace%begin('refine_t%refine_on_box')
  refine_it = .false.
  refine = -1
  if (patch%level >= maxval(region_maxlevel)) then
    call trace%end ()
    return
  end if
  patch%test = 0.0
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
  do k=1,3
    b_llc = region_llc(k,:)
    b_urc = region_urc(k,:)
    intersect = product(max (0.0, min(b_urc,p_urc)-max(b_llc,p_llc))) > 0.0
    if (verbose > 1) then
      unit = merge (stdout, stdlog, verbose > 2)
      write (unit,1) patch%id, patch%level, k, p_llc(3), b_llc(3), intersect
    1 format("refine_region:",i6,":",i2.2,i3,2f7.1,l3)
    end if
    if (.not. intersect .or. patch%level >= region_maxlevel(k)) then
      cycle
    endif
    do iz=m3%li,m3%ui
      z = m3%p + m3%r(iz)
      do iy=m2%li,m2%ui
        y = m2%p + m2%r(iy)
        do ix=m1%li,m1%ui
          x = m1%p + m1%r(ix)
          refine_this_one = (x >= b_llc(1) .and. x <= b_urc(1) .and. &
                             y >= b_llc(2) .and. y <= b_urc(2) .and. &
                             z >= b_llc(3) .and. z <= b_urc(3) )
          patch%test(ix,iy,iz) = merge (1,0,refine_this_one)
          refine_it = refine_it .or. refine_this_one
        end do
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Combine with previous tests and evaluate
  !-----------------------------------------------------------------------------
  call self%evaluate ('region', patch, refine)
  call trace%end()
END FUNCTION refine_region

!===============================================================================
!> Experiment-specific refinement, defined in e.g. experiment_mod or extras_mod
!===============================================================================
FUNCTION custom_refine (self, patch) RESULT (refine)
  class(refine_t) :: self
  class(experiment_t), pointer :: patch
  integer:: refine
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%custom_refine')
  patch%test = 0.0
  refine = patch%check_refine (patch%test)
  !-----------------------------------------------------------------------------
  ! Combine with previous tests and evaluate
  !-----------------------------------------------------------------------------
  call self%evaluate ('custom', patch, refine)
  call trace%end ()
END FUNCTION custom_refine

!===============================================================================
!> Add test and evaluate, allowing derefine below 0.0 and refininig for values 
!> above 1.0
!===============================================================================
SUBROUTINE evaluate (self, label, patch, refine_it)
  class(refine_t):: self
  character(len=*):: label
  class(experiment_t):: patch
  integer:: refine_it, unit
  !.............................................................................
  real:: max_crit, max_test, norm
  integer:: l(3), u(3)
  !-----------------------------------------------------------------------------
  ! Add to previous values of patch%frefine, with a normalization such that
  ! test > 1.0 triggers refinement, independent of if quadrature is used when
  ! combining criteria
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%evaluate')
  l = patch%mesh%lo
  u = patch%mesh%uo
  norm = (1.+derefine_limit)
  max_test = maxval(patch%test(l(1):u(1),l(2):u(2),l(3):u(3)))
  call patch%lock%set ('refine_t%evaluate')
  if (quadrature) then
    patch%frefine = patch%frefine + norm*patch%test**2
  else
    patch%frefine = patch%frefine + norm*patch%test
  end if
  call patch%lock%unset ('refine_t%evaluate')
  !-----------------------------------------------------------------------------
  ! If refine remains <= 0.0 derefinement is allowed, and if it 
  ! exceeds 1.0, refinement is needed
  !-----------------------------------------------------------------------------
  max_crit = maxval(patch%frefine(l(1):u(1),l(2):u(2),l(3):u(3)))
  refine_it = 0
  if (max_crit <= 0.0) refine_it = -1
  if (max_crit >= 1.0) refine_it = +1
  !-----------------------------------------------------------------------------
  if (verbose > 2 .or. &
     (verbose > 1 .and. patch%level < self%levelmax .and. refine_it > 0) .or. &
     (verbose > 1 .and. patch%level > self%levelmin .and. refine_it < 0)) then
    unit = merge (stdout, stdlog, verbose > 2)
    write (unit,'(a22,3x,a,i6,":",i2.2,3i3,f12.6,2f7.3,i4,2x,a1)') &
      label//' criterion: ', 'id,time,max_test,max_crit,refine:', &
      patch%id, patch%key, patch%time, max_test, max_crit, refine_it, &
      merge ('R',' ', &
      refine_it==1)
  end if
  if (patch%level <= self%levelmin) refine_it = max(0,refine_it)
  if (patch%level >= self%levelmax) refine_it = min(0,refine_it)
  call trace%end ()
END SUBROUTINE evaluate

!===============================================================================
!> Check and if needed create supporting patches for a new child patch.
!> The steps needed are:
!>
!> 1) Loop over the possible positions of level L nbors
!> 2) Add missing tasks to the task list, with status set to avoid activity
!> 3) Construct nbor lists for each of the new tasks
!> 4) Add parent link and task pointers
!> 5) Download values from source tasks on the nbor lists
!> 6) Add the new tasks to the nbor lists of their nbors
!> 7) Schedule the new tasks
!>
!> As long as these tasks are not on the nbor list of any existing task, they
!> will not be reached by any check_nbors() loop over nbors. In the meantime
!> the tasks may be given nbor lists, and the nbor lists may be used to load
!> values from the source tasks, provided that adding the tasks does not imply
!> changes of the nbor list attributes in other nbor lists.
!>
!> To ensure consistent nbor lists and nbor relations, check_nbor_relations()
!> must be executed before releasing the tasks. This is done by calling
!> task_list%update_status, or by calling check_nbor_relations() explicitly
!===============================================================================
RECURSIVE SUBROUTINE create_support (self, child, parent, dll_list)
  class(refine_t)              :: self
  class(experiment_t), pointer :: child, parent
  type(dll_t)                  :: dll_list
  !.............................................................................
  class(experiment_t), pointer :: patch, grand
  class(task_t),       pointer :: task
  class(link_t),       pointer :: link
  class(dll_node_t),   pointer :: node, next
  class(*),            pointer :: ptr
  integer                      :: i, j, k, key(4), key2(4)
  real(8)                      :: corner(3), pos(3)
  !-----------------------------------------------------------------------------
  ! Not relevant for single level refinement
  !-----------------------------------------------------------------------------
  if (parent%level <= shared%levelmin) &
    return
  !-----------------------------------------------------------------------------
  ! Find the relevant "corner", which is child%size/2 in the same direction as
  ! child%position-parent%position
  !-----------------------------------------------------------------------------
  corner = parent%position + 2d0*(child%position - parent%position)
  !-----------------------------------------------------------------------------
  ! Loop over the 8 possible level L patch positions and test for existence
  !-----------------------------------------------------------------------------
                        call self%task_list%lock%set ('refine_t%create_support')
  do k=-1,1,2; do j=-1,1,2; do i=-1,1,2
    pos = corner + [i,j,k]*child%size                 ! support position
    if (outside_box (pos)) &
      cycle
    call parent%periodic_wrap (pos)                   ! wrap position
    key = parent%hash_key (position=pos)              ! support key, overloadable until experiment_t
    call self%task_list%hash_table%get (key, link)    ! exists?
    !---------------------------------------------------------------------------
    ! Create new level L tasks where there are none before
    !---------------------------------------------------------------------------
    if (.not. associated(link)) then                  ! exists?
      key2 = parent%parent%hash_key (position=pos)    ! support position, overloadable until experiment_t
      call self%task_list%hash_table%get (key2, link) ! get grandparent link
      if (.not.associated(link)) &
        call io%abort ('refine_t%create_support: missing grandfather')
      task => link%task                               ! get grandparent task
      select type (task)
      class is (experiment_t)
        grand => task
      end select
      call self%create_patch (key, grand, patch, &    ! support patch
        pos, child=child, label='support')            ! at child time
      allocate (node)                                 ! allocate dll_node
      node%car => patch%link                          ! point to new link
      call dll_list%append (node)                     ! append to the dll
      call self%create_support (patch, grand, dll_list)
    end if
  end do; end do; end do
                      call self%task_list%lock%unset ('refine_t%create_support')
contains
  logical function outside_box (pos)
    real(8):: pos(3)
    outside_box = any(.not.parent%periodic .and. &
       (pos < parent%origin .or. pos > parent%origin + parent%box))
  end function outside_box
END SUBROUTINE create_support

!===============================================================================
!> Create a new experiment_t task with the level and position given by key
!===============================================================================
SUBROUTINE create_patch (self, key, parent, patch, position, child, label)
  class(refine_t):: self
  integer:: key(4)
  class(link_t), pointer:: link
  class(experiment_t), pointer:: parent, patch
  class(experiment_t), pointer, optional:: child
  real(8):: position(3)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  ! Allocate and initialize a patch with the given task key
  !-----------------------------------------------------------------------------
  allocate (patch)
  allocate (halo_t:: link)
  patch%box      = parent%box                   ! total size
  patch%origin   = parent%origin                ! region LLC
  patch%status   = parent%status                ! for now, may change
  patch%periodic = parent%periodic              ! copy to all
  patch%n        = parent%n                     ! same resolution
  patch%rank     = mpi%rank                     ! make sure
  patch%size = parent%size/shared%refine_ratio  ! new size
  patch%ds = patch%size/patch%n                 ! cell size
  link%task => patch                            ! set task pointer
  patch%link => link                            ! set link pointer
  patch%parent => parent                        ! task up-link
  patch%link%parent => parent%link              ! link up-link
  call patch%set_properties (key=key)           ! set level and position
  patch%position = position                     ! make it bitwise exact
  patch%ipos     = parent%ipos
  call link%init ()                             ! initialize link
  call patch%pre_init()                         ! pre-initialize task
  call patch%init ()                            ! initialize
  call patch%support_arrays ()                  ! create arrays for AMR
  call patch%init_bdries ()                     ! BC bits
  if (present(child)) then
    patch%time = child%time                     ! given time
    patch%time  = child%time                    ! child time
    patch%t     = child%time                    ! copy to all slots
    patch%dtime = child%dtime/self%ratio        ! for info only
    patch%dt    = child%dtime                   ! copy to all
    patch%iout  = child%iout                    ! child's output index
    patch%out_next = child%out_next             ! next output time
  else
    patch%time  = parent%time                   ! parent time
    patch%t     = parent%time                   ! copy to all slots
    patch%dtime = parent%dtime/self%ratio       ! for info only
    patch%dt    = parent%dtime                  ! copy to all
    patch%iout  = parent%iout                   ! parent's output index
    patch%out_next = parent%out_next            ! next output time
  end if
  patch%tag   = 0                               ! for initial send
  patch%iit   = 1                               ! initial slots
  patch%it    = 1                               ! slot for current values
  patch%new   = 2                               ! slot for new values
  patch%seq   = -1                              ! triggers tag=0 for new task
  patch%istep = 0                               ! for diagnostics
  patch%n_needed = 1                            ! reset to initial
  patch%check_refine_next = check_interval      ! next check
  patch%iit(patch%nt) = patch%new               ! slot for new values
  patch%dnload_time = -1.0                      ! for duplicate call test
  patch%check_refine_next = 0                   ! allow recursive refine
  patch%task_list => self%task_list             ! pointer to common list
  !-----------------------------------------------------------------------------
  ! Initial values from ICs at t==0, from prolong for t>0 
  !-----------------------------------------------------------------------------
  if (patch%time > 0d0) then                    ! ICs are used at t==0
    !call parent%check_gz_density (verbose, &    ! check parent sanity
    !  label='parent')
    call guards%prolong (patch, parent, &       ! load values from parent
              all_tgt=.true., all_src=.true.)   ! use parent guard zones
    call patch%bndries ()                       ! load BCs
    !call patch%check_gz_density (verbose, &     ! check child sanity
    !  label='prolong')
  end if
  !call patch%set (bits%ignore)                  ! prevent use for now
  call self%task_list%add_by_level (link, count=.true.) ! add to task list
  !$omp atomic update
  timer%refined(patch%level) = &
  timer%refined(patch%level) + 1
  !.............................................................................
  if (verbose > 1) then
    if (present(label)) then
      write (stdout,'(a,i6,":",i2.2,3i3)') &
        'created '//label//' with id,key=', patch%id, patch%key
    else
      write (stdout,'(a,i6,":",i2.2,3i3)') &
        'created patch   with id,key=', patch%id, patch%key
    end if
  end if
END SUBROUTINE create_patch

END MODULE refine_mod
