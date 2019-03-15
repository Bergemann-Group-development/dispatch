!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> This module handles checking max change between neighboring points.  Each
!> instance of it needs an index to the variable, whether to take the log, and
!> a value for the max allowed change.  It might be applied to log density, log
!> pressure, log opacity, velocity amplitude, magnetic field amplitude, etc.
!>
!> task_list_t%update
!>   refine_t%check_current
!>     refine = refine_t%refine_needed
!>       refine_t%check_Jeans
!>       refine_t%check_whatever
!>       ...
!>   ...
!>   if refine > 0:
!>     lock_link_and_nbors
!>     split_patch
!>       tasklist%lock%set
!>       tasklist%init_all_nbors
!>       tasklist%reset_status
!>       tasklist%lock%unset
!>     unlock_link_and_nbors
!>     ...
!>     ...
!>   elif refine < 0:
!>     lock_link_and_nbors
!>     remove_patch
!>       tasklist%lock%set
!>       tasklist%remove_link (this)
!>       tasklist%init_nbor_nbors (this)
!>       tasklist%check_nbor_nbors (this)
!>       tasklist%lock%unset
!>     unlock_link_and_nbors
!>     ...
!>   ...
!>
!> Download causes a lock on the target link, which is enough to prevent any
!> nbor from refining, if and only if refining and derefining requires a lock on
!> all nbors of the affected patch.  This is implemented.
!>
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
    procedure:: check_support
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
  logical:: selective=.true.                  ! Create new patches selectively
  integer:: verbose=0                         ! Noise level
  type(refine_t), public:: refine
CONTAINS

!===============================================================================
!> Initialize refinement parameters
!===============================================================================
SUBROUTINE init (self)
  class(refine_t):: self
  integer:: iostat
  integer, save:: levelmin, levelmax    ! powers of ratio relative to box
  integer, save:: ratio=2               ! split ratio when refined
  logical, save:: first_time = .true.
  namelist /refine_params/ on, levelmin, levelmax, check_interval, ratio, &
    refine_from_t, refine_until_t, min_dx, max_dx, max_jeans, min_jeans, &
    max_shock, min_shock, max_contact, min_contact, max_grad, min_grad, &
    grad_var, region_llc, region_urc, force_cubical, selective, verbose
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
  ! If AMR is active, we must use check_filled, with nbor levels sorted into
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
!> Check if refinement is desired and needed for the current task.  Each
!> criterion is represented by a function, which return +1 to indicate refinement,
!> 0 to indicate keep the patch, and -1 to indicate patch removal.
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
  if (verbose > 3) then
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
  ! Shock and/or contact discontinuity detector
  if (max_shock > 0.0 .or. max_contact > 0.0) &
    refine = max (refine, self%shock_and_cd_detector(patch))
  !
  ! Gradient detector
  if (any(max_grad > 0.0)) refine = max (refine, self%gradient_detector(patch))
  !
  ! Refine on (a) user-specified region(s)
  if (any(region_llc /= 0.0) .and. any(region_urc /= 0.0)) &
    refine = max (refine, self%refine_region(patch))
  !
  ! Do not derefine if the patch is needed for support
  refine = max (refine,self%check_support (link))
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
  if (verbose > 3) &
    write (io_unit%output,'(i6,2x,a,f12.6,2i4)') &
      link%task%id, 'refine_t%check_current, at t =', &
      link%task%time, link%task%level, self%levelmax
  !---------------------------------------------------------------------------
  call trace%begin('refine_t%check_current', itimer=itimer)
  patch => gpatch%cast2gpatch (link%task)
  id = patch%id
  !---------------------------------------------------------------------------
  ! Only check leaf patches
  !-----------------------------------------------------------------------------
  if (patch%is_clear(bits%not_leaf)) then
    if (verbose > 3) &
      tasklist%verbose = verbose
    !---------------------------------------------------------------------------
    ! Setup mask, and check if refinement or derefinement is needed
    !---------------------------------------------------------------------------
    call self%setup_mask (patch)
    !---------------------------------------------------------------------------
    ! Prepare per-cell refinement array
    !---------------------------------------------------------------------------
    if (selective) then
      if (associated(patch%irefine)) then
        print *, 'irefine already allocated', patch%id, shape(patch%irefine)
      else
        allocate (patch%irefine(patch%gn(1),patch%gn(2),patch%gn(3)))
      end if
      patch%irefine = -1
    end if
    refine = self%refine_needed (tasklist, link)
    !---------------------------------------------------------------------------
    ! Where the mask is false, there is overlap with better resolved patches,
    ! so should not be refined or derefined
    !---------------------------------------------------------------------------
    if (selective) then
      where (.not.patch%mask)
        patch%irefine = 0
      end where
      !-------------------------------------------------------------------------
      ! Check criteria, and decide if to use selective refinement
      !-------------------------------------------------------------------------
      if (any(patch%irefine /= 0)) then
        call selective_refine (self, tasklist, patch, was_refined, was_derefined)
          deallocate (patch%irefine)
        call trace%end (itimer)
        return
      end if
    end if
    !---------------------------------------------------------------------------
    ! A positive refine value indicates that the patch needs to be refined
    !---------------------------------------------------------------------------
    if (refine > 0 .and. patch%level < self%levelmax) then
      if (tasklist%verbose > 2) &
        write (io%output,*) wallclock(),' thread',omp%thread,' waitfor tasklist(r)'
      if (tasklist%verbose > 2) &
        write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(r)'
      if (verbose > 0) &
        write (io%output,'(a,i6,i3,2x,a)') &
          'id, level =', patch%id, patch%level, 'refine'
      !-------------------------------------------------------------------------
      ! The current task link, and all its nbor links must be locked during
      ! refine.  All nbors must be locked before the next thread is allowed
      ! to try the same, in case they have competing interests.
      !-------------------------------------------------------------------------
      !$omp critical (refine_nbor_cr)
      nbor => link%nbor
      call link%lock%set ('refine_t%split_patch link tasks')
      do while (associated(nbor))
        call nbor%link%lock%set ('refine_t%split_patch nbor tasks')
        nbor => nbor%next
      end do
      !-------------------------------------------------------------------------
      call self%split_patch (tasklist, link)
      !-------------------------------------------------------------------------
      ! The nbor list now also contains the new child links, which should be
      ! locked when created.  All nbor links are thus unlocked here.
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
    !---------------------------------------------------------------------------
    ! A negative refine value indicates that the patch is not needed
    !---------------------------------------------------------------------------
    else if (refine < 0 .and. patch%level > self%levelmin) then
      if (verbose > 0) &
        write (io%output,'(a,i6,i3,2x,a)') &
          'id, level =', patch%id, patch%level, 'derefine'
      !-------------------------------------------------------------------------
      ! The current task link, and all its nbor links must be locked during
      ! derefine.  All nbors must be locked before the next thread is allowed
      ! to try the same, in case they have competing interests.
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
      call self%remove_patch (tasklist, link)
      !-------------------------------------------------------------------------
      ! The nbor lists are now different, in that they lack the removed task
      ! (their locks are unlocked before removal), but the ones that remain
      ! need to be unlocked.
      !-------------------------------------------------------------------------
      nbor => first_nbor
      do while (associated(nbor))
        call nbor%link%lock%unset ('refine_t remove_patch parent nbor tasks')
        nbor => nbor%next
      end do
      !$omp end critical (refine_nbor_cr)
      !-------------------------------------------------------------------------
      was_derefined = .true.
    end if
    if (associated(patch%mask)) &
      deallocate (patch%mask)
  end if
  call trace%end (itimer)
END SUBROUTINE check_current

!===============================================================================
! Loop over the ratio x ratio x ratio positions and determine wether to refine
! or derefine
!===============================================================================
SUBROUTINE selective_refine (self, tasklist, patch, was_refined, was_derefined)
  class(refine_t):: self
  class(list_t):: tasklist
  logical:: was_refined, was_derefined
  class(gpatch_t), pointer:: patch, nbpatch
  integer:: i, j, k, i1, i2, j1, j2, k1, k2, n1, n2, n3, ii(3), refine
  real:: rpos(3)
  real(8):: pos(3)
  logical:: already_refined
  class(link_t), pointer:: nbor, first_nbor
  class(link_t), pointer:: link
  associate (m1=>patch%mesh(1), m2=>patch%mesh(2), m3=>patch%mesh(3))
  !-----------------------------------------------------------------------------
  ! First, consider if the whole patch should be derefined.  This requires that
  ! it has no higher level patch inside (would trigger patch%mask == false there),
  ! and secondly that irefine(:,:,:) indicates that derefine would be OK
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%selective_refine')
  if (patch%level > self%levelmin) then
    if (all(patch%mask)) then
      if (all(patch%irefine < 0)) then
        !$omp critical (refine_nbor_cr)
        link => patch%link
        first_nbor => link%nbor
        call lock_link_and_nbors (link, 'derefine')
        call remove_patch (self, tasklist, link)
        call unlock_link_and_nbors (link, 'derefine', first_nbor=first_nbor)
        !$omp end critical (refine_nbor_cr)
        was_derefined = .true.
      end if
    end if
  !-----------------------------------------------------------------------------
  ! Now consider if any part should be refined
  !-----------------------------------------------------------------------------
  else if (patch%level < self%levelmax) then
    if (any(patch%irefine > 0)) then
      n1 = m1%n/self%ratio
      n2 = m2%n/self%ratio
      n3 = m3%n/self%ratio
      do k=1,self%ratio
        k1 = m3%li + (k-1)*n3
        k2 = k1 + n3 - 1
        rpos(3) = (k-(self%ratio+1)*0.5)*n3*m3%d
        do j=1,self%ratio
          j1 = m2%li + (j-1)*n2
          j2 = j1 + n2 - 1
          rpos(2) = (j-(self%ratio+1)*0.5)*n2*m2%d
          do i=1,self%ratio
            i1 = m1%li + (i-1)*n1
            i2 = i1 + n1 - 1
            !-----------------------------------------------------------------------
            ! Should this part be refined or not?
            !-----------------------------------------------------------------------
            refine = maxval (patch%irefine(i1:i2,j1:j2,k1:k2))
            rpos(1) = (i-(self%ratio+1)*0.5)*n1*m1%d
            pos = rpos + patch%position
            ii = patch%index_only(pos)
            !print '(3i3,2x,3i3,3(2x,3f9.4))',i,j,k,ii,pos,patch%position,patch%size
            !-----------------------------------------------------------------------
            ! Is this part refined or not?
            !-----------------------------------------------------------------------
            already_refined = .not.patch%mask(ii(1),ii(2),ii(3))
            if (patch%level > self%levelmin) &
              print *, 'selective_refine:', patch%id, i, j, k, refine, already_refined
            if (.not.already_refined .and. refine > 0) then
              !---------------------------------------------------------------------
              ! Create a refined patch here
              !---------------------------------------------------------------------
              write (io_unit%output,'(a,i6,1p2g16.5,2x,3(i4,":",i2))') 'xxx', &
                patch%id, patch%fmaxval(patch%mem(:,:,:,patch%idx%d,patch%it,1)), &
                     maxval(patch%mem(i1:i2,j1:j2,k1:k2,patch%idx%d,patch%it,1)), &
                i1, i2, j1, j2, k1, k2
              !$omp critical (refine_nbor_cr)
              call lock_link_and_nbors (patch%link, 'refine')
              call create_patch (self, tasklist, patch, i, j, k)
              call unlock_link_and_nbors (patch%link, 'refine')
              !$omp end critical (refine_nbor_cr)
              was_refined= .true.
            end if
          end do
        end do
      end do
    end if
  end if
  end associate
  call trace%end
END SUBROUTINE selective_refine

!===============================================================================
!===============================================================================
SUBROUTINE lock_link_and_nbors (link, label)
  class(link_t), pointer:: link, nbor
  character(len=*):: label
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%lock_link_and_nbors')
  nbor => link%nbor
  call link%lock%set ('refine_t%lock_link_and_nbors '//label)
  do while (associated(nbor))
    call nbor%link%lock%set ('refine_t%lock_link_and_nbors '//label)
    nbor => nbor%next
  end do
  call trace%end
END SUBROUTINE lock_link_and_nbors

!===============================================================================
!===============================================================================
SUBROUTINE unlock_link_and_nbors (link, label, first_nbor)
  class(link_t), pointer:: link, nbor
  class(link_t), pointer, optional:: first_nbor
  character(len=*):: label
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%unlock_link_and_nbors')
  if (present(first_nbor)) then
    nbor => first_nbor
  else
    nbor => link%nbor
  end if
  do while (associated(nbor))
    call nbor%link%lock%unset ('refine_t%unlock_link_and_nbors '//label)
    nbor => nbor%next
  end do
  if (associated(link)) &
    call link%lock%unset ('refine_t%unlock_link_and_nbors '//label)
  call trace%end
END SUBROUTINE unlock_link_and_nbors

!===============================================================================
!===============================================================================
SUBROUTINE create_patch (self, tasklist, patch, ix, iy, iz)
  class(refine_t):: self
  class(list_t):: tasklist
  class(gpatch_t), pointer:: patch
  integer:: ix, iy, iz
  !.............................................................................
  class(link_t), pointer:: link
  real(8), dimension(3):: pos, split_size
  !-----------------------------------------------------------------------------
  ! Lock the task list, and compute the position of the new patch
  !-----------------------------------------------------------------------------
  call tasklist%lock%set ('refine_t refine')
  split_size = patch%size/self%ratio
  pos = patch%position + ([ix,iy,iz] - (self%ratio+1)*0.5) * split_size
  call make_child_patch(self, tasklist, patch%link, patch, pos, split_size, 0)
  !-----------------------------------------------------------------------------
  ! At this point the patch is (possibly again) not a leaf patch
  !-----------------------------------------------------------------------------
  call patch%set(bits%not_leaf)
  !-----------------------------------------------------------------------------
  ! Reset the task list and unlock
  !-----------------------------------------------------------------------------
  call reset_task_list (self, tasklist, patch%link)
  call tasklist%lock%unset ('refine_t refine')
END SUBROUTINE create_patch

!===============================================================================
!===============================================================================
SUBROUTINE reset_task_list (self, tasklist, patch_link)
  class(refine_t):: self
  class(list_t):: tasklist
  class(link_t), pointer:: patch_link, link
  !-----------------------------------------------------------------------------
  ! Lock the task list, and compute the position of the new patch
  !-----------------------------------------------------------------------------
  call tasklist%lock%set ('refine_t create_patch')
  ! positions, then reset their status based on the nbor lists, and count
  !-----------------------------------------------------------------------------
  !call tasklist%init_nbor_nbors (patch_link)
  call tasklist%init_all_nbors
  call tasklist%reset_status
  !-----------------------------------------------------------------------------
  ! Expensive way to make sure the state of readiness is consistent (FIXME)
  !-----------------------------------------------------------------------------
  link => tasklist%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call tasklist%check_ready (link)
    link => link%next
  end do
  call tasklist%init_bdries
  !-----------------------------------------------------------------------------
END SUBROUTINE reset_task_list

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
  if (.not.associated(patch%mask)) then
    l = patch%mesh%lb
    u = patch%mesh%ub
    allocate (patch%mask(l(1):u(1),l(2):u(2),l(3):u(3)))
    patch%mask = .true.
  end if
  !-----------------------------------------------------------------------------
  ! Check for higher level patches inside this
  !-----------------------------------------------------------------------------
  nbor => patch%link%nbor
  do while (associated(nbor))
    if (patch%contains(nbor%task%position)) then
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
!> Create a new patch using an existing "parent" patch as the source.
!===============================================================================
SUBROUTINE make_child_patch (self, tasklist, parent_link, parent, pos, child_size, i)
  class(refine_t):: self
  class(list_t) :: tasklist
  class(link_t), pointer :: parent_link, child_link, child_nbor
  class(gpatch_t), pointer :: parent, child
  real(8), dimension(3), intent(in) :: pos, child_size
  integer:: i, count, l(3), u(3)
  !------------------------------------------------------------------------
  call trace%begin('refine_mod::make_child_patch')
  !------------------------------------------------------------------------
  ! Allocate a new task, 1/nsplit the size of original one, and adjust its
  ! position in each direction with "pos", in units of the new size.
  !
  ! Using sourced allocation below effectively clones the parent.  This is
  ! necessary, to pick up properties that may have been set at the experiment
  ! level, which would otherwise be left undefined, or in inconsistent state
  !------------------------------------------------------------------------
  allocate (child, source=parent)
  child%mem_allocated = .false.
  nullify (child%irefine)
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
  ! Increment I/O write counter
  !------------------------------------------------------------------------
  io%nwrite = io%nwrite + 1
  count = counters%increment (parent%iout+1, io%nwrite)
  !------------------------------------------------------------------------
  ! Append links to the new patch in the task list. Neighbor lists are
  ! constructed at the end of `split_patch`.
  !------------------------------------------------------------------------
  allocate (child_link)
  call child_link%init                               ! initialise the lock
  call child_link%lock%set ('child init')            ! immediately lock it
  child_link%task => child                           ! associate task with link
  child_link%parent => parent_link                   ! used in remove_patch
  child%link => child_link                           ! associate link with task
  child%parent => parent                             ! task relation
  call tasklist%append_link (child_link)             ! which is added to the patch list
  io%ntask = tasklist%na                             ! update the task count
  !------------------------------------------------------------------------
  ! Make parent a neighbor of new child patch. This will be overwritten
  ! in `init_nbor_nbors` but is necessary here to permit interpolation
  ! of values from parent to child.
  !------------------------------------------------------------------------
  allocate (child_nbor)
  child_nbor%link => parent_link
  child_nbor%task => parent
  child_link%nbor => child_nbor
  !------------------------------------------------------------------------
  ! Download values from parent and neighbors, except at time = 0, to allow
  ! recursive refinement there.  If time = 0, patch should have been filled
  ! with values during `child%init()` ! call.
  !------------------------------------------------------------------------
  if (parent%time > 0.0 .and. child%restart <= 0) &
    call download%download_link (child_link, all_cells=.true.)
  if (verbose > 0 .and. i==0) &
    write (io%output,'(a,i6,i3,2x,a,2x,i6,i3)') 'id, level =', &
      parent%id, parent%level, 'refine to', child%id, child%level
  if (verbose > 1) then
    if (i==0) &
    l = parent%mesh%li
    u = parent%mesh%ui
    !print 1,'refine_t%make_child_patch: parent, level, iout, time, size(3), position(3) =', &
    !  parent%id, parent%level, parent%iout, parent%time, parent%size, parent%position, &
    !  maxval(parent%mem(l(1):u(1),l(2):u(2),l(3):u(3),parent%idx%d,parent%it,1))
    1 format(a,3i6,f10.6,2(2x,3f8.4),1pe12.2)
    !print 1,'                            child, level, iout, time, size(3), position(3) =', &
    !  child%id, child%level, child%iout, child%time, child%size, child%position, &
    !  maxval(child%mem(l(1):u(1),l(2):u(2),l(3):u(3),child%idx%d,child%it,1))
  end if
  call trace%end()
END SUBROUTINE make_child_patch

!===============================================================================
!> Split a patch into ratio (parameter) in each direction, merging the new patches
!> into the list.  Note that the link to the original patch is kept, with a "leaf"
!> flag being cleared, indicating it is not longer a "leaf" patch.
!===============================================================================
SUBROUTINE split_patch (self, tasklist, this)
  class(refine_t):: self
  class(list_t):: tasklist
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
  call tasklist%lock%set ('refine_t split_patch tasklist')
  tail => tasklist%tail
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
        j1 = j - patch%mesh%n/2
        j2 = j1 + patch%mesh%n - 1
        write (io_unit%output,'(a,i6,1p2g16.5,2x,3(i4,":",i2))') 'xxx', &
          patch%id, patch%fmaxval(patch%mem(:,:,:,patch%idx%d,patch%it,1)), &
          maxval(patch%mem(j1(1):j2(1),j1(2):j2(2),j1(3):j2(3),patch%idx%d,patch%it,1)), &
          j1(1), j2(1), j1(2), j2(2), j1(3), j2(3)
        if (associated(patch%mask)) then
          if (patch%mask(j(1),j(2),j(3))) then
            call make_child_patch(self, tasklist, this, patch, pos, split_size, i)
            i = i+1
          end if
        else
          call make_child_patch(self, tasklist, this, patch, pos, split_size, i)
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
  call tasklist%init_nbor_nbors (this)
  call tasklist%reset_status
  !-----------------------------------------------------------------------------
  ! Also, expensive way to make sure the state of readiness is consistent (FIXME)
  !-----------------------------------------------------------------------------
  link => tasklist%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call tasklist%check_ready (link)
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
  call tasklist%init_bdries
  if (verbose > 3) &
    call tasklist%print_tasks()
  call tasklist%lock%unset ('refine_t split_patch tasklist')
  !-----------------------------------------------------------------------------
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
SUBROUTINE remove_patch (self, tasklist, this)
  class(refine_t):: self
  class(list_t):: tasklist
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
  call tasklist%lock%set ('refine_t%remove_patch tasklist')
  call tasklist%remove_link (this)
  call tasklist%init_nbor_nbors (this)
  call tasklist%check_nbor_nbors (this)
  call this%lock%unset ('delete')
  call this%delete (this)
  nullify (this)
  !-----------------------------------------------------------------------------
  call tasklist%lock%unset ('refine_t%remove_patch tasklist')
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE remove_patch

!===============================================================================
!> Check of link%task is a patch that may be removed (refine = -1) or not (= 0)
!> A patch may not be removed if it contains a patch on the nbor list with a
!> higher level.
!===============================================================================
FUNCTION check_support (self, link) RESULT (refine)
  class(refine_t):: self
  class(link_t), pointer:: link
  integer:: refine
  !.............................................................................
  class(link_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin('refine_t%check_support')
  refine = -1
  associate (task => link%task)
  select type (task)
  class is (gpatch_t)
  nbor => link%nbor
  do while (associated(nbor))
    if (nbor%task%level > task%level .and. task%contains (nbor%task%position)) then
      if (verbose > 1) then
        write (io_unit%output,'(i6,i3,2x,a,i6,i3)') task%id, task%level, &
          'required for support of', nbor%task%id, nbor%task%level
      end if
      refine = 0
      exit
    end if
    nbor => nbor%next
  end do
  end select
  end associate
  if (verbose > 1 .and. link%task%level > self%levelmin) &
    write (io_unit%output,'(a,i6,i4)') &
      'support criterion: id, refine', link%task%id, refine
  call trace%end()
END FUNCTION check_support

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
!> Refine / derefine on ratio of cell size to Jeans lenght (max_jeans/min_jeans)
!===============================================================================
FUNCTION check_Jeans (self, patch) RESULT(refine_it)
  USE math_mod
  USE scaling_mod
  class(refine_t):: self
  class(gpatch_t), pointer:: patch
  integer:: refine_it
  !.............................................................................
  integer:: imin(3), l(3), u(3)
  real:: Jeans_num, dmin, c
  real, dimension(:,:,:), allocatable:: c_sound2, pg
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e
  real, allocatable:: test(:,:,:)
  logical:: masked
  !-----------------------------------------------------------------------------
  call trace%begin ('refine_t%check_check_jeans')
  d => patch%mem(:,:,:,patch%idx%d,patch%it,1)
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
    test = c*d/sqrt(pg)
    if (associated(patch%mask)) then
      Jeans_num = maxval(test, mask=patch%mask)
    else
      Jeans_num = patch%fmaxval(test)
    end if
    !---------------------------------------------------------------------------
    ! Compute a per-cell criterion, which setes patch%irefine to +1 where it
    ! needs refinement, and sets it to 0 where the current resolution is needed.
    ! By implication, where it does not set patch%irefine, and where it remains
    ! at the intial -1 value, derefinement is allowed
    !---------------------------------------------------------------------------
    if (selective) then
      where (test > max_jeans)
        patch%irefine = +1
      !-------------------------------------------------------------------------
      ! Jeans number not small enough, or overlying refinement
      !-------------------------------------------------------------------------
      else where (test > min_jeans)
        patch%irefine = max (patch%irefine, 0)
      end where 
    end if
    deallocate (test)
    refine_it = 0
    if (Jeans_num > max_jeans) refine_it = +1
    if (Jeans_num < min_jeans) refine_it = -1
    !---------------------------------------------------------------------------
    ! A patch with existing child patches has those regions masked out in the
    ! calculation of Jeans_num above, so it might indicate derefinement,
    ! causing the child patch to loose support (or causing support check to
    ! repeatedly call for refinement.  To prevent this, if the mask shows that
    ! there are refine regions inside, derefinement is off table.
    !---------------------------------------------------------------------------
    if (.not.all(patch%mask)) &
      refine_it = max(refine_it,0)
    if (verbose > 2 .or. &
       (verbose > 1 .and. patch%level < self%levelmax .and. refine_it > 0) .or. &
       (verbose > 1 .and. patch%level > self%levelmin .and. refine_it < 0)) then
      if (associated(patch%mask)) then
        masked = .not.all(patch%mask)
      else
        masked = .false.
      end if
      write (io_unit%output,'(a,i6,1p,3e9.2,i4,e9.2,l3)') &
        'Jeans criterion: id, max, num, min, refine, maxval(d)', &
        patch%id, max_jeans, Jeans_num, min_jeans, refine_it, maxval(d), &
        masked
    end if
  end if
  deallocate (pg)
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
!> Refine on gradients
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
