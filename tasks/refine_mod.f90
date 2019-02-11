!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> This module handles checking max change between neighboring points.  Each
!> instance of it needs an index to the variable, whether to take the log, and
!> a value for the max allowed change.  It might be applied to log density, log
!> pressure, log opacity, velocity amplitude, magnetic field amplitude, etc.
!>
!> task_list_t%update
!>   refine_t%check_current
!>     refine_t%refine_needed
!>       refine_t%check_sinkparticles
!>       refine_t%check_Jeans
!>       ...
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE refine_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE bits_mod
  USE kinds_mod
  USE link_mod
  USE list_mod
  USE task_mod
  USE patch_mod
  USE solver_mod
  USE sinkparticles_mod
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
    procedure:: check_cubical
    procedure:: check_jeans
    procedure:: shock_and_cd_detector
    procedure:: refine_region
    procedure:: split_patch
    procedure:: remove_patch
    procedure:: gradient_detector
    procedure:: check_nbor_support
    procedure:: check_sinkparticle
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
  !
  type(refine_t), public:: refine

CONTAINS
!> =============================================================================
!> Initialize refinement parameters
!> =============================================================================
SUBROUTINE init (self)
  class(refine_t):: self
  integer:: iostat
  integer, save:: levelmin, levelmax
  integer, save:: ratio=2                     ! How does a patch split when it is refined
  logical, save:: first_time = .true.
  namelist /refine_params/ on, levelmin, levelmax, check_interval, ratio, &
    refine_from_t, refine_until_t, min_dx, max_dx, max_jeans, min_jeans, &
    max_shock, min_shock, max_contact, min_contact, max_grad, min_grad, &
    grad_var, region_llc, region_urc, force_cubical
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
  !$omp end critical (refine_cr)
  self%levelmin = levelmin
  self%levelmax = levelmax
  self%ratio = ratio
  call sinkparticles%init
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Check if refinement is desired and needed for the current task
!> =============================================================================
FUNCTION refine_needed (self, tasklist, link)
  class(refine_t):: self
  class(list_t):: tasklist
  class(link_t), pointer:: link
  logical:: refine_needed
  !.............................................................................
  class(patch_t), pointer:: patch
  real:: dxmin, dxmax
  integer, save:: itimer=0
  logical:: new_sink
  !----------------------------------------------------------------------------
  refine_needed = .false.
  if (.not. on) return
  call trace%begin('refine_t%needed', itimer=itimer)
  !
  patch => task2patch(link%task)
  !
  ! filter out patches which do not match the necessary criteria
  if (patch%is_set(bits%virtual)) then
    call trace%end (itimer)
    return
  end if
  if (patch%istep < patch%check_refine_next) then
    call trace%end (itimer)
    return
  end if
  dxmin = minval(patch%ds)
  dxmax = maxval(patch%ds)
  if (dxmin < min_dx .or. dxmax >= max_dx) then
    call trace%end (itimer)
    return
  end if
  if (patch%time < refine_from_t .or. patch%time > refine_until_t) then
    call trace%end (itimer)
    return
  end if
  !
  !------------------------------------------------------------------------
  ! Refinement Criteria
  !------------------------------------------------------------------------
  ! sink particle creation
  new_sink = self%check_sinkparticle (tasklist, link)
  !
  ! Ensure consistent, cubical patch size (e.g. after RAMSES startup)
  if (force_cubical) refine_needed = refine_needed .or. self%check_cubical(patch)>0
  if (io%verbose>1) print *, patch%id, 'refine%needed (patch_n)   =', refine_needed
  !
  ! Jeans criterion
  if (min_jeans > 0.0) refine_needed = refine_needed .or. self%check_Jeans(patch)>0
  if (io%verbose>1) print *, patch%id, 'refine%needed (Jeans)   =', refine_needed
  !
  ! shock and/or contact discontinuity detector
  if (max_shock > 0.0 .or. max_contact > 0.0) &
    refine_needed = refine_needed .or. self%shock_and_cd_detector(patch)>0
  if (io%verbose>1) print *, patch%id, 'refine%needed (shock)   =', refine_needed
  !
  ! gradient detector
  if (any(max_grad > 0.0) ) refine_needed = refine_needed .or. self%gradient_detector(patch)>0
  if (io%verbose>1) print *, patch%id, 'refine%needed (grad)    =', refine_needed
  !
  ! refine on (a) user-specified region(s)
  if (any(region_llc /= 0.0) .and. any(region_urc /= 0.0)) &
    refine_needed = refine_needed .or. self%refine_region(patch)>0
  if (io%verbose>1) print *, patch%id, 'refine%needed (region)  =', refine_needed
  !
  ! Check that neighbors of `self` are not more than one level different.
  ! If they are, mark this patch for splitting so that support can be
  ! provided to the neighbors of `patch`.
  refine_needed = refine_needed .or. self%check_nbor_support (tasklist, link)>0
  if (io%verbose>1) print *, patch%id, 'refine%needed (support) =', refine_needed
  !
  patch%check_refine_next = patch%istep + check_interval
  call trace%end (itimer)
END FUNCTION refine_needed

!> =============================================================================
!> Check if refinement is desired and needed for the current task; if so, push
!> new tasks onto the (ready) queue with the same time as the parent, i.e., to
!> the head of the queue.  If the task is a leaf task, and the level of refinement
!> does not appear to be needed, suggest derefinement.
!> =============================================================================
SUBROUTINE check_current (self, tasklist, link, was_refined)
  class(refine_t):: self
  class(list_t):: tasklist
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch
  logical, intent(out):: was_refined
  real:: dxmin, dxmax
  integer:: refine
  integer, save:: itimer=0
  !.............................................................................
  was_refined = .false.
  if (.not. on) return
  if (io%verbose>1) &
    print *, link%task%id, 'MK: refine_t%check_current, at t =', link%task%time, link%task%level, self%levelmax
  if (link%task%level >= self%levelmax) return
  call trace%begin('refine_t%check_current', itimer=itimer)
  !
  patch => task2patch(link%task)
  if (patch%is_clear(bits%not_leaf)) then
    refine = self%refine_needed (tasklist, link)
    !---------------------------------------------------------------------------
    ! A negative refine value indicates that the patch is not needed
    !---------------------------------------------------------------------------
    if (refine < 0 .and. patch%level > self%levelmin) then
      write (io%output,*) patch%id, 'derefine suggested'
      call tasklist%lock%set
      call self%remove_patch (tasklist, link)
      call check_level_constraint
      call tasklist%lock%unset
    !---------------------------------------------------------------------------
    ! A positive refine value indicates that the patch needs to be refined
    !---------------------------------------------------------------------------
    else if (refine > 0 .and. patch%level < self%levelmax) then
      call tasklist%lock%set
      call self%split_patch (tasklist, link)
      call check_level_constraint
      call tasklist%lock%unset
      was_refined = .true.
    end if
  end if
  !
  call trace%end (itimer)
contains
  !> ===========================================================================
  !> Iteratively check the nbor level constraint, splitting patches until the
  !> there are no more violations.  Since this inserts new patches in the patch
  !> list, and modifies the queue, it must be done in a critical region, blocking
  !> other updates to the patch list.
  !> ===========================================================================
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
            if (io%verbose>2) &
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
      if (io%verbose>2) &
        print *, 'check_level_constraint, iteration =', &
        iteration, still_needed
    end do
  end subroutine check_level_constraint
END SUBROUTINE check_current

!> =============================================================================
!> Create a new patch using an existing "parent" patch.
!> This procedure is used frequently during refinement.
!> =============================================================================
SUBROUTINE make_child_patch (self, tlist, parent_link, parent, pos, child_size)
  class(refine_t):: self
  class(list_t) :: tlist
  class(link_t), pointer :: parent_link, child_link, child_nbor
  class(patch_t), pointer :: parent, child
  real, dimension(3), intent(in) :: pos, child_size
  !------------------------------------------------------------------------
  call trace%begin('refine_mod::make_child_patch')
  !------------------------------------------------------------------------
  ! Allocate a new task, 1/nsplit the size of original one, and adjust its
  ! position in each direction with "pos", in units of the new size.
  !
  ! Using sourced allocation below effectively clones the parent.
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
  child%dt(child%it) = parent%dt(parent%it) / self%ratio ! reset the time-step
  child%t(child%it) = child%time
  if (child%restart <= 0) then
    child%out_next = child%time                        ! force output now
    child%iout = max(0,child%iout - 1)                 ! correct output counter
  end if
  child%parentid = parent%id                         ! for debugging
  child%check_refine_next = 0                        ! allow recursive refine
  !------------------------------------------------------------------------
  ! Append links to the new patch in the patch list, sorted by quality
  ! Neighbour lists are constructed at the end of `split_patch`.
  !------------------------------------------------------------------------
  allocate (child_link)
  child_link%task => child
  child_link%parent => parent_link                   ! link is the link to the new patch
  call tlist%append_link (child_link)                ! which is added to the patch list
  io%ntask = tlist%na                                ! update the task count
  !------------------------------------------------------------------------
  ! Make parent a neighbour of new child patch.
  ! This will be overwritten in `init_all_nbors` but is necessary here to permit
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
  if (parent%time > 0.0 .and. child%restart <= 0) call download%download_link (child_link, all_cells=.true.)
  if (io%verbose>0) then
    1 format(a,3i6,f10.6,2(2x,3f8.4))
    print 1,'refine_t%make_child_patch:  child, level, iout, time, size(3), position(3) =', &
      child%id, child%level, child%iout, child%time, child%size, child%position
    print 1,'                           parent, level, iout, time, size(3), position(3) =', &
      parent%id, parent%level, parent%iout, parent%time, parent%size, parent%position
  end if

  ! increment I/O write counter
  io%nwrite = io%nwrite + 1
  call tlist%check_ready(child_link)

  call trace%end
END SUBROUTINE make_child_patch

!> =============================================================================
!> Split a patch into two (parameter) in each direction, merging the new patches
!> into the list.  Note that the link to the original patch is kept, with a "leaf"
!> flag being cleared, indicating it is not longer a "leaf" patch.
!> =============================================================================
SUBROUTINE split_patch (self, tlist, this)
  class(refine_t):: self
  class(list_t):: tlist
  class(link_t), pointer:: this, link
  class(patch_t), pointer:: this_patch
  integer:: ix, iy, iz
  integer:: nsplit(3)
  real, dimension(3):: pos, rpos, split_size
  !.............................................................................
  if (this%task%is_set(bits%not_leaf)) return
  call trace%begin('refine_t%split_patch', 1)
  !
  this_patch => task2patch(this%task)
  if (this_patch%msplit>0) then
    nsplit = this_patch%msplit
  else
    nsplit = self%ratio
  end if
  !
  where (this_patch%size==0 .or. this_patch%n==1)
    nsplit = 1
  end where
  !
  split_size = this_patch%size / nsplit
  do ix=1,nsplit(1)
    rpos(1) = (ix-(nsplit(1)+1)*0.5) * split_size(1)
    do iy=1,nsplit(2)
      rpos(2) = (iy-(nsplit(2)+1)*0.5) * split_size(2)
      do iz=1,nsplit(3)
        rpos(3) = (iz-(nsplit(3)+1)*0.5) * split_size(3)
        pos = rpos + this_patch%position
        call make_child_patch(self, tlist, this, this_patch, pos, split_size)
      end do
    end do
  end do
  call this_patch%set(bits%not_leaf)
  !
  call tlist%init_all_nbors
  call tlist%reset_status
  link => tlist%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call tlist%check_ready (link)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! FIXME: special care may be needed in multi-rank context
  !-----------------------------------------------------------------------------
  call tlist%init_bdries
  !
  if (io%verbose > 0) call tlist%print_tasks()
  !
  call trace%end
END SUBROUTINE split_patch

!> =============================================================================
!> Remove a patch that has been deemed not necessary.  At this point, the task
!> is not in the ready queue, 
!> =============================================================================
SUBROUTINE remove_patch (self, tlist, this)
  class(refine_t):: self
  class(list_t):: tlist
  class(link_t), pointer:: this, link, parent
  !----------------------------------------------------------------------------
  ! Start processing, blocking (at least for now) all queue-related activity
  !----------------------------------------------------------------------------
  call trace%begin ('refine_t%remove_patch')
  !$omp critical (queue_cr)
  !----------------------------------------------------------------------------
  ! Remove the patch from the task list, and update the nbor list of all nbors
  !----------------------------------------------------------------------------
  call tlist%remove_link (this)
  !----------------------------------------------------------------------------
  ! Send a suicide note to virtual nbors, if this is a boundary patch
  !----------------------------------------------------------------------------
  if (this%task%is_set (bits%boundary)) then
    call this%task%set (bits%remove)
    call this%send_to_vnbors
  end if
  !----------------------------------------------------------------------------
  ! Update the nbor lists of all nbors
  !----------------------------------------------------------------------------
  link => this%nbor
  do while (associated(link))
    call tlist%refresh_nbors (link)
    link => link%next
  end do
  !----------------------------------------------------------------------------
  ! Change the status of the parent patch to leaf
  !----------------------------------------------------------------------------
  associate (parent_task => this%task%parent)
  select type (parent_task)
  class is (patch_t)
  parent => parent_task%link
  end select
  end associate
  call parent%task%clear (bits%not_leaf)
  !----------------------------------------------------------------------------
  ! Deallocate the patch memory, and deallocate the patch instance
  !----------------------------------------------------------------------------
  call this%task%delete()
  deallocate (this%task)
  !----------------------------------------------------------------------------
  ! Decrement the I/O write counter
  !----------------------------------------------------------------------------
  io%nwrite = io%nwrite - 1
  !----------------------------------------------------------------------------
  ! Check the ready status of the parent patch
  !----------------------------------------------------------------------------
  call tlist%check_ready (parent)
  !----------------------------------------------------------------------------
  !$omp end critical (queue_cr)
  call trace%end
END SUBROUTINE remove_patch

!> =============================================================================
!> Check the neighbors of the task associated with `link` to see if they can
!> provide support to newly refined patches made from `link`.
!> The current definition of "support" is a difference in levels of one or less.
!>
!> In the absence of the optional argument, `check_only`, neighbour patches will
!> be split as needed.
!> =============================================================================
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

  refine_it = 0
  this => link%task
  nbor => link%nbor
  do while (associated(nbor))
    npatch => nbor%task
    ! support is defined as having a level difference no greater than one
    if ((npatch%level - this%level) > 1) then
      if (io%verbose > 0) print "('Patch',i5.1,' at level ',i2.1,' cannot provide support to patches at level',i2.1)",this%id,this%level,npatch%level
      if (split_if_necessary) then
        ! mark this patch for splitting to ensure support for neighbor patches
        if (io%verbose > 0) print "('Splitting',i5.1,' in order to provide support')",this%id
        refine_it = 1
        exit
      end if
    end if
    nbor => nbor%next
  end do

  call trace%end
END FUNCTION check_nbor_support

!> =============================================================================
!> Check that the patch dimensions are the desired cubical ones
!> =============================================================================
FUNCTION check_cubical (self, patch) RESULT (refine_it)
  class(refine_t):: self
  class(patch_t), pointer:: patch
  integer:: refine_it
  !.............................................................................
  if (all(patch%n==scene%patch_n)) then
    refine_it = 0
  else
    refine_it = 1
    patch%msplit = 1
  end if
END FUNCTION check_cubical

!> =============================================================================
!> Calculate the Jeans length and compare it to the patch resolution.
!> =============================================================================
FUNCTION check_jeans (self, patch) RESULT(refine_it)
  USE math_mod
  USE scaling_mod
  class(refine_t):: self
  class(patch_t), pointer:: patch
  integer:: refine_it
  real:: Jeans_length, Jeans_num
  real, dimension(:,:,:), allocatable:: c_sound2, pg
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e
  !.............................................................................
  call trace%begin ('refine_t%check_check_jeans')
  associate(l => patch%li, u => patch%ui)
  d => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%it,1)
  e => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%e,patch%it,1)
  allocate(c_sound2(size(d,1),size(d,2),size(d,3)), &
           pg(patch%gn(1),patch%gn(2),patch%gn(3)))

  select type (patch)
  class is (solver_t)
    pg = patch%gas_pressure()
  class default
    call mpi%abort("Unrecognised patch type. Abort!")
  end select
  c_sound2 = patch%gamma * pg(l(1):u(1),l(2):u(2),l(3):u(3)) / d
  Jeans_length = minval(sqrt(math%pi * c_sound2 / (d*scaling%grav)))
  Jeans_num = minval(patch%mesh%d) / Jeans_length
  refine_it = merge(1,0,Jeans_num > max_jeans)
  refine_it = merge(-1,0,Jeans_num < min_jeans)

  if (io%verbose > 0) print*,'Jeans criterion:', patch%id, Jeans_num, min_jeans, refine_it

  deallocate (c_sound2, pg)
  nullify (d,e)
  end associate
  call trace%end
END FUNCTION check_jeans

!> =============================================================================
!> Refine on shocks, with threshold set by tolerance parameter `tol_shock`.
!> =============================================================================
FUNCTION shock_and_cd_detector (self, patch) RESULT (refine)
  USE mesh_mod
  class(refine_t):: self
  class(patch_t), pointer:: patch
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
  !
  nullify (mx, my, mz, d)
  deallocate (pg, v)
  !
  call trace%end
END FUNCTION shock_and_cd_detector

!> =============================================================================
!> Refine on shocks, with threshold set by tolerance parameter `max_shock`.
!> =============================================================================
FUNCTION gradient_detector (self, patch) RESULT (refine)
  USE mesh_mod
  USE solver_mod

  class(refine_t):: self
  class(patch_t), pointer:: patch
  integer:: refine
  !.............................................................................
  class(mesh_t), pointer:: mx, my, mz
  integer:: i, j, k, ione=0, jone=0, kone=0, im1, jm1, km1, ip1, jp1, kp1, n
  real, parameter:: small=1e-10
  real(kind=KindScalarVar):: dqmax
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: q, dq1, dq2, dq3
  real, parameter:: rsmall = 1.0e-10
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

    ! calculate the gradient in each direction
    do k=mz%lo,mz%uo
      km1 = k - kone
      kp1 = k + kone
      do j=my%lo,my%uo
        jm1 = j - jone
        jp1 = j + jone
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
          if (dqmax > max_grad(n)) refine = 1
          if (dqmax < min_grad(n)) refine = -1
        end do
      end do
    end do

  end do

  deallocate (dq1, dq2, dq3, q)
  call trace%end
END FUNCTION gradient_detector

!> =============================================================================
!> Given a box defined by its lower and upper corners, flag a patch for refinement
!> if it lies inside this box, up to a maximum level.
!> This particular criteria is primarily used for testing and debugging.
!> =============================================================================
FUNCTION refine_region (self, patch) RESULT (refine)
  USE mesh_mod
  class(refine_t) :: self
  class(patch_t), pointer :: patch
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
  call trace%end
END FUNCTION refine_region

!> =============================================================================
!> Test for the need to create a new sink particle, and handle the task list
!> aspect of that; creating an nbor list, and adding the link to the task list.
!> =============================================================================
FUNCTION check_sinkparticle (self, task_list, link)
  class(refine_t):: self
  type(list_t):: task_list
  class(link_t), pointer:: link
  logical:: check_sinkparticle
  !.............................................................................
  class(task_t), pointer:: task, patch
  class(patch_t), pointer:: p
  class(link_t), pointer:: new_link, nbor
  class(sinkparticles_t), pointer:: new_sink
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  integer:: m(3)
  logical:: do_trace
  !-----------------------------------------------------------------------------
  if (.not.sinkparticles%on) return
  check_sinkparticle = .false.
  patch => link%task
  select type (patch)
  class is (sinkparticles_t)
    !---------------------------------------------------------------------------
    ! Refresh nbor list
    !---------------------------------------------------------------------------
    call task_list%refresh_nbors(link)
  class is (patch_t)
    p => patch
    check_sinkparticle = sinkparticles%check (p)
    if (check_sinkparticle) then
      do_trace = io%do_trace
      !io%do_trace = .true.
      !-------------------------------------------------------------------------
      ! Allocate a new sinkparticle_t instance, as the task of a new link
      !-------------------------------------------------------------------------
      allocate (new_link)
      allocate (sinkparticles_t:: new_link%task)
      new_link%parent => link
      task => new_link%task
      !-------------------------------------------------------------------------
      ! Create the new sinkparticle, using data from the MHD patch
      !-------------------------------------------------------------------------
      select type (task)
      class is (sinkparticles_t)
        task%link => new_link
        call task%create (p)
        !-------------------------------------------------------------------------
        ! Initialize its neighbor list, and add it to the task list
        !-------------------------------------------------------------------------
        call task_list%init_nbors (new_link)
        call task_list%append_link (new_link)
        nbor => new_link%nbor
        do while (associated(nbor))
          print *, 'nbor:', nbor%task%id, real((nbor%task%position-task%position)/task%ds)
          nbor => nbor%next
        end do
        !------------------------------------------------------------------------
        ! Download variables and accrete
        !------------------------------------------------------------------------
        call task%dnload
        d => task%mem(:,:,:,task%idx%d,task%it,1)
        d = max(d,task%fminval(d))
        if (task%id==io%id_debug) then
          m = nint(task%offset)
          print *, 'sinkparticles_t%create:  task%pos, d',  task%myposition(m), &
            task%mem(m(1),m(2),m(3),task%idx%d,task%it,1)
          print *,'MM', minval(d), task%fminval(task%idx%d), &
                        maxval(d), task%fmaxval(task%idx%d)
        end if
        call task%accrete
      end select
      !-------------------------------------------------------------------------
      ! Refresh the nbor lists of the parent and its nbors; the parent should
      ! be done after the nbors, to avoid putting the child link on its nbor
      ! list until the other ones are updated.  All nbors that happen to have
      ! overlap with the sink particle patch will have it as an nbor
      !-------------------------------------------------------------------------
      nbor => link%nbor
      do while (associated(nbor))
        call task_list%init_nbors (nbor%link)
        nbor => nbor%next
      end do
      call task_list%init_nbors (link)
      io%do_trace = do_trace
    end if
  end select
END FUNCTION check_sinkparticle

END MODULE refine_mod
