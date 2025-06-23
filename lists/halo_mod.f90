!===============================================================================
!> Module that implements halo loading of tasks, using either guard zone
!> loading or boundary conditions.  This relies on going through an init
!> procedure when creating new nbor lists, which allocates either a guard
!> zone loader or a boundary condition data type, depending on the type of
!> boundary.
!>
!> The no-mans-land problem (lack of support for refinement) may be handled
!> by ad-hoc updating the guard zones of the parent patch to the time of
!> child patch.
!>
!> Thread safety: This only requires the normal locking of those nbors, since
!> if another thread is simultaneously updating the parent patch, it will be
!> prevented from doing it at the same time.  However, this requires some
!> checking back and forth, so a simpler approach is to fill in values in a
!> copy of the %mem of the parent, and use that in the final guard zone
!> interpolations.
!===============================================================================
MODULE halo_mod
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE io_mod
  USE os_mod
  USE timer_mod
  USE trace_mod
  USE bits_mod
  USE task_mod
  USE mesh_mod
  USE patch_mod
  USE link_mod
  ! USE guards_mod
  USE region_mod
  USE shared_mod
  ! USE interp_mod  
  implicit none
  private

  type, public, extends(nbor_t):: guards_t
  end type
  !.............................................................................
  type, public, extends(link_t):: halo_t
    class(patch_t), pointer:: coarse => null()
    integer:: n_op(3)=0
    logical:: first_time=.true.
  contains
    procedure:: init
    procedure:: update
    procedure:: initialize
    ! procedure:: download
    ! procedure:: make_coarse
    !procedure:: init_child_nbors
    procedure:: conditionally_add_nbor
    procedure:: info
    procedure:: set_filled
    procedure:: get_filled
    procedure:: not_filled
    procedure, nopass:: alloc_link
    procedure, nopass:: alloc_nbor
    procedure:: check_nbor_list
  end type
  logical, save:: debug=.false.
  logical, save:: check_filled=.true.
  logical, save:: detailed_timer=.false.
  integer, save:: verbose=0
  integer, save:: method=1
  integer, save:: n_tot(4)=0
  real(8), save:: overlap_margin = 1d-10
  !$omp threadprivate (n_tot)

  procedure(update), public, pointer :: halo_update => update

CONTAINS

!===============================================================================
!> Initialize the link, and possibly the nbor nodes
!===============================================================================
SUBROUTINE init (self)
  class(halo_t):: self
  !.............................................................................
  class(nbor_t), pointer:: nbor
  integer:: iostat
  namelist /halo_params/ verbose, debug, detailed_timer, check_filled, method
  logical, save:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%init')
  call self%link_t%init
  !-----------------------------------------------------------------------------
  call io%lock%set ('halo_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, halo_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, halo_params)
    end if
    call io%namelist_warning ('halo_params', iostat)
    if (shared%refine_on) then
      check_filled = .true.
    end if
    write (io%namelist, halo_params)
    call garbage%init_verbose (verbose)
    first_time = .false.
  end if
  self%verbose = verbose
  call self%init_verbose (verbose)
  call io%lock%unset ('halo_t%init')
  !call guards%pre_init
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Walk through the nbor list and update guard zones.  The call will go to
!> the guards_t subclass of nbor_t, if that was used when allocating the nbor
!> list
!===============================================================================
SUBROUTINE update (self, all_tgt, all_src, only)
  implicit none
  class(halo_t)           :: self
  class(nbor_t), pointer  :: nbor
  logical, optional       :: all_tgt, all_src
  integer, optional       :: only
  !.............................................................................
  class(task_t), pointer  :: target, task
  class(nbor_t), pointer  :: nbors
  integer                 :: unit
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  if (self%task%is_set (bits%no_download)) &
    return
  call trace%begin ('halo_t%update',itimer=itimer)
  call self%check_nbor_list ()
  unit = merge (stdout, stdlog, verbose > 2)
  self%n_op = 0
  !-----------------------------------------------------------------------------
  ! If istep is zero the task has just been created, including guard zones
  !-----------------------------------------------------------------------------
  target => self%task
  !.............................................................................
  if (verbose > 1) then
    write (unit,'(a,i6,1p,4e15.6)') &
      'halo begin: id,time,dnltime=', &
      target%id, target%time, target%dnload_time
    flush (unit)
  end if
  !.............................................................................
  if (target%time > target%dnload_time) then
    if (target%is_set (bits%initialize)) then
      call self%initialize (self%nbor, only)
    else
      ! call                                           self%lock%set_read ('halo')
      call self%download (self%nbor, only, all_tgt, all_src)
      ! call                                         self%lock%unset_read ('halo')
    end if
    if (verbose > 1) then
      write (unit,'(a,i6,2e15.6)') &
        'halo done : id,time,dnltime=', &
        target%id, target%time, target%dnload_time
      flush (unit)
    end if
    !...........................................................................
    target%dnload_time = target%time
  else
    !...........................................................................
    if (verbose > 1) then
      write (unit,'(a,i6,2e15.6)') &
        'halo skip', target%id, target%time, target%dnload_time
    end if
    !...........................................................................
    self%n_op = 0
  end if
  !.............................................................................
  if (target%id == io%id_debug) then
    task => self%task
    select type (task)
    class is (patch_t)
      call task%stats ('after halo', stdout)
    end select
  end if
  select type (target)
  class is (patch_t)
    call target%check_gz_density (verbose, label='link%update')
  end select
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Initialize from all sources
!===============================================================================
SUBROUTINE initialize (self, nbors, only)
  class(halo_t)          :: self
  class(task_t), pointer :: tgt, src
  class(nbor_t), pointer :: nbors
  integer,      optional :: only
  class(nbor_t), pointer :: nbor
  integer                :: lt(3), ut(3)
  !-----------------------------------------------------------------------------
  ! tgt => self%task
  ! select type (tgt)
  ! class is (patch_t)
    ! lt = tgt%li
    ! ut = tgt%ui
    ! nbor => nbors
    ! do while (associated(nbor))
      ! src => nbor%task
      ! select type (src)
      ! class is (patch_t)
        ! call interp%from_source (tgt, lt, ut, src, filled=tgt%filled)
      ! end select
      ! nbor => nbor%next
    ! end do
  ! end select
  call tgt%clear (bits%initialize)
END SUBROUTINE initialize

! !===============================================================================
! !> Comprehensive entry point
! !===============================================================================
! SUBROUTINE download (self, nbors, only, all_tgt, all_src)
!   class(halo_t)          :: self
!   class(nbor_t), pointer :: nbors
!   integer,      optional :: only
!   logical,      optional :: all_tgt, all_src
!   !.............................................................................
!   class(nbor_t), pointer :: tnbors
!   integer                :: unit, n_same, tid, i
!   logical                :: use_coarse, all_tgt_L, all_src_L
!   class(nbor_t), pointer :: nbor
!   class(task_t), pointer :: target, source, patch
!   integer, save          :: itimer=0
!   !-----------------------------------------------------------------------------
!   ! Decide whether to use a coarse clone or not
!   !-----------------------------------------------------------------------------
!   call trace%begin ('halo_t%download', itimer=itimer)
!   all_tgt_L = shared%option (all_tgt)
!   all_src_L = shared%option (all_src)
!   unit = merge (stdout, stdlog, verbose > 2)
!   target => self%task
!   tid = target%id
!   !-----------------------------------------------------------------------------
!   ! If the task does not yet have support it will soon have it, and until then
!   ! we make sure to fill all guard zones with values from the parent task, which
!   ! has at least time-extrapolated guard zone values, and possibly already
!   ! authoritative guard zone values
!   !-----------------------------------------------------------------------------
!   call self%set_filled (target)
!   !-----------------------------------------------------------------------------
!   ! Check if we need to use a coarse level clone.
!   ! FIXME: Cases where all level L-1 nbors are covered by level L nbors need
!   ! not use a coarse clone. The case with 26 same level nbors is a special case,
!   ! but there are additional ones to discover.
!   !-----------------------------------------------------------------------------
!   use_coarse = .false.
!   n_same = 0
!   nbor => nbors
!   do while (associated(nbor))
!     if (nbor%task%is_set (bits%remove)) &
!       nbor%download = .false.
!     if (nbor%download) then
!       source => nbor%task
!       if (verbose > 3) then
!         write (unit,'(a8,2x,2(a,i6,":",i2.2,1x,3i4,3x))') 'halo', &
!           'tgt:', target%id, target%key, &
!           'src:', source%id, source%key
!       end if
!       !-------------------------------------------------------------------------
!       ! An nbor ROA that covers the target ROI is sufficient 
!       !-------------------------------------------------------------------------
!       if (source%authority%covers (target%interest, nbor%distance)) then
!         if (verbose > 0) &
!           write (stdlog,*) target%id, 'covered by', source%id
!         exit
!       end if
!       !-------------------------------------------------------------------------
!       ! An nbor at level L-1 is potential reason to use a coarse clone, unless
!       ! another nbor is found above, or 26 nbors on the same level are found
!       !-------------------------------------------------------------------------
!       if (source%level < target%level) then
!         use_coarse = .true.
!       else if (source%level == target%level) then
!         n_same = n_same+1
!       end if
!     end if
!     nbor => nbor%next
!   end do
!   do i=1,3
!     select type(target)
!     class is (patch_t)
!     if (target%mesh(i)%lower_boundary) n_same = n_same+9
!     if (target%mesh(i)%upper_boundary) n_same = n_same+9
!     end select
!   end do
!   !-----------------------------------------------------------------------------
!   ! Tasks with complete cover do not need a coarse clone
!   !-----------------------------------------------------------------------------
!   if (n_same==26) then
!     if (verbose > 0) &
!       write (unit,*) target%id, 'has 26 level L nbors'
!     use_coarse = .false.
!   !-----------------------------------------------------------------------------
!   ! Tasks with extended sorce cover do not need a coarse clone
!   !-----------------------------------------------------------------------------
!   else if (all_src_L) then
!     if (verbose > 0) &
!       write (unit,*) target%id, 'all_src'
!     use_coarse = .false.
!   !-----------------------------------------------------------------------------
!   ! Tasks at or below levelmin do not need a coarse clone
!   !-----------------------------------------------------------------------------
!   else if (target%level <= shared%levelmin) then
!     if (verbose > 0) &
!       write (unit,*) target%id, 'target at levelmin'
!     use_coarse = .false.
!   else if (use_coarse) then
!     if (verbose > 0) &
!       write (unit,*) target%id, 'needs a coarse clone'
!   end if
!   if (verbose > 0) &
!     flush (unit)
!   !-----------------------------------------------------------------------------
!   ! Pass on, or loop over nbors
!   !-----------------------------------------------------------------------------
!   if (use_coarse) then
!     call make_coarse (self, only=only, nbors=nbors)
!   else
!     nbor => nbors
!     do while (associated(nbor))
!       if (nbor%download) then
!         if (nbor%task%is_clear (bits%remove)) then
!           source => nbor%task
!           !---------------------------------------------------------------------
!           ! Only allow prolong if all_src is present
!           !---------------------------------------------------------------------
!           if (n_same /= 26 .or. source%level >= target%level .or. all_src_L) then
!             if (verbose > 2) then
!               write (unit,'(a,2(i6,":",i2.2,3i3,2x))') &
!                 'halo_t%download: tid,sid=', &
!                 target%id, target%key, source%id, source%key
!               flush (unit)
!             end if
!             call guards%source_to_target (source, target, n_op=self%n_op, &
!               only=only, all_src=all_src, all_tgt=all_tgt, &
!               nbstep=nbor%check_step)
!           end if
!         end if
!       end if
!       nbor => nbor%next
!     end do
!   end if
!   if (target%id == io%id_debug) then
!     select type (target)
!     class is (patch_t)
!       block
!       integer:: l(3),u(3)
!       real, pointer, contiguous:: d(:,:,:)
!       l = target%li
!       u = target%ui
!       d => target%mem(:,:,:,target%idx%d,target%it,1)
!       write (stdlog,'(a,i6,"::",i2.2,1p,2g12.3)') &
!         'download: tid,min,max=', &
!         target%id, target%level, &
!         minval(d(l(1):u(1),l(2):u(2),l(3):u(3))), &
!         maxval(d(l(1):u(1),l(2):u(2),l(3):u(3)))
!       flush (stderr)
!       end block
!     end select
!   end if
!   !-----------------------------------------------------------------------------
!   ! Boundary conditions
!   !-----------------------------------------------------------------------------
!   call self%task%bndries()
!   call self%task%periodic_grid()
!   !-----------------------------------------------------------------------------
!   ! Emergency procedure: If some cells are not filled; redo the guard zone
!   ! loading, using full source patches, while printing warnings
!   !-----------------------------------------------------------------------------
!   if (self%not_filled (target)) then
!     write (stderr,'(a,i6,":",i2.2,3i7)') &
!       ' halo_t%not_filled: target=', target%id, target%key
!     write (stdlog,'(a,i6,":",i2.2,3i7)') &
!       ' halo_t%not_filled: target=', target%id, target%key
!     select type (target)
!     class is (patch_t)
!     nbor => nbors
!     do while (associated(nbor))
!       source => nbor%task
!       select type (source)
!       class is (patch_t)
!       write (stdlog,'(a,i6,":",i2.2,3i4,3x,l1,2(2x,3l1),f10.4)') &
!         'not_filled: src,key,download,lb,ub,dphi=', &
!         source%id, source%key, nbor%download, &
!         source%mesh%lower_boundary, source%mesh%upper_boundary, &
!         (source%omega(3)-target%omega(3))*target%time*target%size(2)/target%radius
!       end select
!       nbor => nbor%next
!     end do
!     nbor => nbors
!     do while (associated(nbor))
!       if (nbor%download) then
!         source => nbor%task
!         select type (source)
!         class is (patch_t)
!           call source%extrapolate_guards
!         end select
!         call guards%source_to_target (source, target, n_op=self%n_op, &
!           all_src=.true., all_tgt=all_tgt, only=only, nbstep=nbor%check_step)
!       end if
!       nbor => nbor%next
!     end do
!     end select
!     if (self%not_filled (target)) then
!       call target%log ('repair NOK')
!       write (stderr,*) 'halo_t%update: repair NOK, target%id=', target%id
!       write (io_unit%mpi,*) 'halo_t%update: repair NOK, target%id=', target%id
!     else
!       call target%log ('repair OK')
!       write (stderr,*) 'halo_t%update: repair OK, target%id=', target%id
!       write (io_unit%mpi,*) 'halo_t%update: repair OK, target%id=', target%id
!     end if
!   end if
!   call self%get_filled (target)
!   !.............................................................................
!   if (verbose > 2) then
!     write (unit,'(f12.6,a,3i5)') &
!       wallclock(), pad_trace(1:indent)// &
!       'halo_t%download: n_same,prolong,restrict=', self%n_op
!     self%n_op = 0
!   end if
!   call trace%end (itimer)
! END SUBROUTINE download

! !===============================================================================
! !> Simplify optional argument handling
! !===============================================================================
! FUNCTION my_present(flag) result (yes)
!   logical, optional:: flag
!   logical:: yes
!   !-----------------------------------------------------------------------------
!   if (present(flag)) then
!     yes = flag
!   else
!     yes = .false.
!   end if
! END FUNCTION my_present

! !===============================================================================
! !> Create a clone (C) of the coarse parent patch (P) at level L-1, and fill its
! !> guard zones from patches in the child (self) nbor list at level L-1; these
! !> will all be at most interpolations in time and copies in space.  The completed
! !> coarse patch is able to supply all guard zone values for the level L task
! !> by prolongation.
! !===============================================================================
! SUBROUTINE make_coarse (self, nbors, only, all_tgt, only_lower)
!   class(halo_t)                    :: self
!   integer      , optional          :: only
!   class(nbor_t), optional, pointer :: nbors
!   logical      , optional          :: only_lower, all_tgt
!   !.............................................................................
!   class(nbor_t) , pointer :: nbor, nbor_list
!   class(task_t) , pointer :: target, source, parent, task
!   class(patch_t), pointer :: coarse
!   class(mesh_t) , pointer :: mesh
!   real          , pointer :: target_mem(:,:,:,:)
!   real(8)                 :: wc(2)
!   integer                 :: i, m(6), unit, tid, sid, l(3), u(3), iv, pid
!   logical                 :: use_coarse
!   character(16)           :: type
!   integer, save           :: itimer=0
!   !-----------------------------------------------------------------------------
!   call trace%begin ('halo_t%make_coarse', itimer=itimer)
!   unit = merge (stdout, stdlog, verbose > 2)
!   self%n_op = 0
!   target => self%task
!   if (verbose > 1 .or. target%id == io%id_debug) then
!     if (associated(target%parent)) then
!       pid = target%parent%id
!     else
!       pid = 0
!     end if
!     print *, self%task%id, 'make coarse begin, parent=', pid
!   end if
!   !-----------------------------------------------------------------------------
!   ! Choose the nbor lost
!   !-----------------------------------------------------------------------------
!   if (present(nbors)) then
!     nbor_list => nbors
!   else
!     nbor_list => self%nbor
!   end if
!   !-----------------------------------------------------------------------------
!   ! Detect the case where the parent task is the root grid, and cut directly
!   ! to a prolong, in case the patches above are 2x2x2 periodic
!   !-----------------------------------------------------------------------------
!   if (associated(self%parent)) then
!     source => self%parent%task
!     select type (source)
!     class is (patch_t)
!     if (source%is_set (bits%root_grid)) then
!       select type (target)
!       class is (patch_t)
!         call guards%prolong (target, source, only=only, all_tgt=.true., &
!                              all_src=.true.)
!         !-----------------------------------------------------------------------
!         ! For normal patches, use the nbor list as well
!         !-----------------------------------------------------------------------
!         if (.not. shared%option(only_lower)) then
!           nbor => nbor_list
!           do while (associated(nbor))
!             if (nbor%task%level >= target%level .and. nbor%download) then
!               call guards%source_to_target (nbor%task, target, only=only, &
!                 nbstep=nbor%check_step)
!             end if
!             nbor => nbor%next
!           end do
!         end if
!         call trace%end (itimer)
!         return
!       end select
!     end if
!     end select
!   end if
!   !-----------------------------------------------------------------------------
!   ! Allocate a duplicate of the coarse parent
!   !-----------------------------------------------------------------------------
!   source => self%task%parent
!   tid = target%id             ! for vscode
!   sid = source%id             ! for vscode
!   call io%assert (associated(target), "make_coarse: target not associated")
!   call io%assert (associated(source), "make_coarse: source not associated")
!   !.............................................................................
!   if (verbose > 0 .or. &
!       target%id == io%id_debug .or. &
!       source%id == io%id_debug) then
!     write (unit,'(a,2(i6,":",i2.2))') &
!       'halo_t%make_coarse: target,parent=', &
!       target%id, target%level, source%id, source%level
!   !.............................................................................
!   end if
!   allocate (task, source=source)
!   call os%mem_register (storage_size(task), 1, 'coarse')
!   select type (task)
!   class is (patch_t)
!     coarse => task
!   end select
!   m = shape(coarse%mem)
!   call os%allocate(coarse%mem, [m(1),m(2),m(3),m(4),1,1], 'coarse%mem')
!   call os%allocate(coarse%t,   1, 'coarse%t')
!   call os%allocate(coarse%iit, 1, 'coarse%iit')
!   coarse%mem = tiny(1.)
!   coarse%id = target%id
!   coarse%it = 1
!   coarse%nt = 1
!   coarse%iit = 1
!   coarse%key = self%task%key
!   coarse%key(1) = self%task%key(1) - 1
!   coarse%time = self%task%time
!   coarse%t = coarse%time
!   !-----------------------------------------------------------------------------
!   ! Run through the child nbor list and update the clone from the nbors;
!   ! taking care to only use the nbors that are valid for the child, and also
!   ! taking care not to use any other higher level nbors.  One of the child
!   ! nbors is the parent, and this will take care of the copying of its values,
!   ! in principle only if all_tgt==.true (but this is the default behavior
!   ! of the %same procedure).  The parent task is the only one that can fill
!   ! the interior cells.
!   !-----------------------------------------------------------------------------
!   if (debug) then
!     write (90) self%task%id, self%task%level, self%task%time, &
!       self%task%position, self%task%size
!     write (90) self%task%parent%id, self%task%parent%level, &
!       self%task%parent%time, self%parent%task%position, self%parent%task%size
!     target => self%parent%task
!     select type (target)
!     class is (patch_t)
!       write (90) target%mem(:,:,:,6,coarse%it,1)  ! a
!     end select
!     write (90) coarse%mem(:,:,:,6,coarse%it,1)  ! b
!   end if
!   !.............................................................................
!   type = ''
!   ! nbor => self%parent%nbor
!   ! source => self%parent%task
!   ! select type (source)
!   ! class is (patch_t)
!   !   call guards%same (coarse, source, all_tgt=.true., only=only)
!   ! end select
!   nbor => nbor_list
!   do while (associated(nbor))
!     if (nbor%download) then
!       !-------------------------------------------------------------------------
!       ! Use all nbors at the same level as the clone, including the task%parent
!       !-------------------------------------------------------------------------
!       if (nbor%task%level == self%task%parent%level) then
!         source => nbor%task
!         sid = source%id             ! for vscode
!         select type (source)
!         class is (patch_t)
!           if (source%is_set (bits%root_grid)) then
!             call guards%same (coarse, source, all_tgt=.true., only=only, &
!                               all_src=.true.)
!             type = ', root_grid'
!           else
!             !-------------------------------------------------------------------
!             ! In general we can only trust the interior of these sources, but
!             ! together they should fill up the coarse clone -- one is the parent
!             !----------------------------------------------------.--------------
!             call guards%same (coarse, source, all_tgt=.true., only=only)
!             if (source%id == self%task%parent%id) then
!               type = ', parent used'
!             end if
!           end if
!         end select
!       end if
!     end if
!     nbor => nbor%next
!   end do
!   !.............................................................................
!   if (verbose > 1 .or. target%id == io%id_debug) &
!     print *, self%task%id, 'make coarse: parent level nbors'//type
!   call coarse%check_gz_density (verbose=verbose-1, label='coarse')
!   if (debug) &
!     write (90) coarse%mem(:,:,:,6,coarse%it,1)  ! b
!   !-----------------------------------------------------------------------------
!   ! Restrict the level L target to the coarse clone (FIXME: strictly, only the
!   ! "rim" values are needed)
!   !-----------------------------------------------------------------------------
!   if (verbose > 1 .or. target%id == io%id_debug) &
!     print *, self%task%id, 'make coarse: restrict task to coarse'
!   target => coarse
!   if (.not. shared%option(only_lower)) then
!     call guards%source_to_target (self%task, target, n_op=self%n_op, &
!                                   only=only, all_tgt=.true.)
!     if (debug) &
!       write (90) coarse%mem(:,:,:,6,coarse%it,1)    ! c
!   end if
!   !-----------------------------------------------------------------------------
!   ! Update the child guard-zones from the clone (via prolongation), then
!   ! delete it.  To override the level conditions, mark target unfilled.
!   !-----------------------------------------------------------------------------
!   source => coarse            ! the coarse time is identical to self%time
!   target => self%task         ! the self%task to fill guard zones in
!   select type (target)
!   class is (patch_t)
!     if (check_filled .and. associated(target%filled)) then
!       target%filled = 0
!     end if
!     call dbg (self%task, 'before s2t')
!     if (verbose > 1 .or. target%id == io%id_debug) &
!       print *, self%task%id, 'make coarse: prolong coarse to task'
!     call guards%prolong (target, coarse, only=only, all_src=.true., &
!                          all_tgt=all_tgt)
!     call dbg (self%task, ' after s2t')
!     if (debug) &
!       write (90) target%mem(:,:,:,6,target%it,1)  ! e
!     call os%deallocate(coarse%mem, 'coarse%mem')
!     call os%deallocate(coarse%t, 'coarse%t')
!     call os%deallocate(coarse%iit, 'coarse%iit')
!     call os%mem_register (storage_size(coarse), -1, 'coarse')
!     deallocate (coarse)
!     call target%check_gz_density (verbose, label='prolong coarse')
!     !---------------------------------------------------------------------------
!     ! Add / overwrite the guard zones covered by level L and L+1 nbors
!     !---------------------------------------------------------------------------
!     if (.not. shared%option (only_lower)) then
!       if (verbose > 1 .or. target%id == io%id_debug) &
!         print *, self%task%id, 'make coarse: other nbors'
!       !.........................................................................
!       nbor => nbor_list
!       do while (associated(nbor))
!         if (nbor%task%level >= target%level .and. nbor%download) then
!           call guards%source_to_target (nbor%task, target, only=only, &
!             nbstep=nbor%check_step)
!         end if
!         nbor => nbor%next
!       end do
!       call target%check_gz_density (verbose, label='other nbors')
!       if (debug) then
!         write (90) target%mem(:,:,:,6,target%it,1)  ! e
!       end if
!     end if
!     if (verbose > 1 .or. target%id == io%id_debug) &
!       print *, self%task%id, 'make coarse end'
!   end select
!   call trace%end (itimer)
! END SUBROUTINE make_coarse

! !===============================================================================
! SUBROUTINE dbg (target, label)
!   class(task_t):: target
!   character(*):: label
!   if (target%id == io%id_debug) then
!     select type (target)
!     class is (patch_t)
!       call target%stats (label)
!     end select
!   end if
! END SUBROUTINE dbg

!===============================================================================
!> Check if the nbor should actually be on the nbor list, and if so try to add
!> it, in increasing level order.  If the task is already on the nbor list, the
!> attempt will be quietly ignored.
!===============================================================================
SUBROUTINE conditionally_add_nbor (self, nblink, nbors, added)
  class(halo_t)          :: self                    ! the task nblink
  class(link_t), pointer :: nblink                  ! the nbor task nblink
  class(nbor_t), pointer :: nbors                   ! nbor list
  logical, optional      :: added                   ! added or not?
  !.............................................................................
  type(guards_t), pointer:: guard
  class(nbor_t), pointer :: nbor
  class(task_t), pointer :: task, nbtask, patch
  integer                :: n(4), taskid, nbtaskid, unit
  real(8)                :: distance(3)
  integer, save          :: itimer=0
  logical:: added_loc
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('halo_t%conditionally_add_nbor',itimer=itimer)
  n = 0
  if (self%task%rank == mpi%rank .or. nblink%task%rank == mpi%rank) then
    task   => self%task
    patch  => self%task
    taskid = task%id
    select type (patch)
    class is (patch_t)
      nbtask => nblink%task
      nbtaskid = nbtask%id
      if (nbtaskid /= taskid) then
        select type (nbtask)
        class is (patch_t)
          n(1) = n(1) + 1
          !--- make sure it is not pending remove ---
          if (nbtask%is_clear(bits%remove)) then
            !--- does it have a chance to overlap? ---
            distance = task%distance_to (nbtask)
            n(2) = n(2) + 1
            if (patch%interest%overlaps (nbtask%authority, distance)) then
            !if (patch%overlaps (nbtask)) then
            !if (patch%overlaps2 (nbtask)) then
              n(3) = n(3) + 1
              allocate (guard)
              nbor => guard
              call os%mem_register (storage_size(nbor), 1, 'nbor')
              guard%link => nblink
              guard%task => nblink%task
              call self%add_nbor_by_level (nbors, nbor, added_loc)
              if (.not. added_loc) then
                if (present(added)) added=added_loc
                if (detailed_timer) &
                  call trace%end (itimer)
                return
              end if
              nbor%distance = distance
              !.................................................................
              if (nbtaskid == io%id_debug) write (stdlog,*) &
                nbtaskid, 'overlaps_r'
              unit = merge (stdout, stdlog, task%id==io%id_debug)
              if (verbose > 1 .or. task%id==io%id_debug) then
                write (unit,'(f12.6,2x,a,2i6,":",i2.2,2(1x,3f7.3),2(2x,4i4))') &
                  wallclock(), &
                  'halo_t%cnd_add_nbor: tgt,src,dist,size,n=', &
                  taskid, nbtaskid, nbor%task%level, &
                  nbor%distance/task%size, &
                  nbtask%distance(patch)/task%size, &
                  task%hash_key(), nbtask%hash_key()
              end if
              !.................................................................
            end if
          end if
        end select
      end if
    end select
    n_tot = n_tot + n
  end if
  if (detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE conditionally_add_nbor

!===============================================================================
!> Initialize the link, and possibly the nbor nodes
!===============================================================================
SUBROUTINE info (self, label, unit)
  class(halo_t):: self
  character(*), optional:: label
  integer, optional:: unit
  !.............................................................................
  integer:: unit_L
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%info')
  unit_L = shared%option (unit, stdlog)
  if (verbose < 4) then
    if (verbose > 0) then
      call self%link_t%info (label, unit)
    end if
  else
    if (present(label)) then
      write (unit_L,'(f12.6,a,10i6)') wallclock(), &
        pad_trace(1:indent)//label//': n_tot=', n_tot
      n_tot = 0
    else
      write (unit_L,'(f12.6,a,10i6)') wallclock(), &
        pad_trace(1:indent)//'halo_t%info: n_tot=', n_tot
      n_tot = 0
    end if
    call self%print_nbors (label='halo_t%info', verbose=verbose, unit=unit)
  end if
  call trace%end ()
END SUBROUTINE info

!===============================================================================
!> Initialize the filled(:,:,:) and outside(:,:,:) arrays
!>
!>       HP Zbook 6G, gfortran 8.4
!> --------------------------------
!> method,ns/pt=  1   1.538   3.540
!> method,ns/pt=  2   2.026   3.711
!> method,ns/pt=  3   0.000   3.418
!> method,ns/pt=  4   0.000   3.418
!===============================================================================
SUBROUTINE set_filled (self, task)
  class(halo_t)           :: self
  class(task_t) , pointer :: task
  class(patch_t), pointer :: target
  class(mesh_t) , pointer :: mesh
  real(8)                 :: wc(2)
  integer                 :: i1, i2, i3
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. check_filled .or. task%is_set(bits%no_filled)) return
  call trace%begin ('halo_t%set_filled', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Cast to patch_t
  !-----------------------------------------------------------------------------
  select type (task)
  class is (patch_t)
  target => task
  end select
  !-----------------------------------------------------------------------------
  ! Clear target%filled and initialize internal points to target%level.
  ! Se logical mask 'outside' true in guard zones, false inside.
  ! Test different methods to initialize; bell whether of optimization
  !-----------------------------------------------------------------------------
  call target%support_arrays ()
  ! if (guards%verbose > 0) wc(1) = wallclock()
  if (method == 1) then
    target%filled = 0
  else if (method == 2) then
    do i3=1,target%gn(3)
      do i2=1,target%gn(2)
        !$omp simd
        do i1=1,target%gn(1)
          target%filled(i1,i2,i3) = 0
        end do
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Method 3 and 4 use a single loop nest, with logics inside
  !-----------------------------------------------------------------------------¨
  ! if (guards%verbose > 0) wc(2) = wallclock()
  if (method >= 3) then
    wc(1) = wc(2)
    mesh => target%mesh(3)
    do i3=1,target%gn(3)
      mesh => target%mesh(2)
      do i2=1,target%gn(2)
        !-----------------------------------------------------------------------
        ! Simple logics, easy to read
        !-----------------------------------------------------------------------
        mesh => target%mesh(1)
        if (method == 3) then
          do i1=1,target%gn(1)
            if (target%outside(i1,i2,i3)) then
              target%filled(i1,i2,i3) = 0
            else
              target%filled(i1,i2,i3) = 10*target%level
            end if
          end do
        !-----------------------------------------------------------------------
        ! Attempt to aid vectorization
        !-----------------------------------------------------------------------
        else
          do i1=1,target%gn(1)
            target%filled(i1,i2,i3) = &
              merge (0, 10*target%level, target%outside(i1,i2,i3))
          end do
        end if
      end do
    end do
  else
    !---------------------------------------------------------------------------
    ! Final parts of methods 1 and 2
    !---------------------------------------------------------------------------
    do i3=target%li(3),target%ui(3)
      do i2=target%li(2),target%ui(2)
        !$omp simd
        do i1=target%li(1),target%ui(1)
          target%filled(i1,i2,i3) = 10*target%level
        end do
      end do
    end do
  end if
  ! if (guards%verbose > 0) then
  !   write (stdlog,'(f12.6,a,i3,9p,2f8.3)') &
  !     wallclock(), pad_trace(1:indent)// &
  !     'halo_t%set_filled: filled,outside method,ns/pt=', method, &
  !     (wc(2)-wc(1))/product(target%n), (wallclock()-wc(1))/product(target%n)
  ! end if
  call trace%end (itimer)
END SUBROUTINE set_filled

!===============================================================================
!> Check that target is filled
!===============================================================================
SUBROUTINE get_filled (self, target)
  class(halo_t)          :: self
  class(task_t), pointer :: target
  integer                :: i1, i2, i3
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. check_filled .or. &
      target%is_set (bits%root_grid+bits%no_filled)) return
  call trace%begin ('halo_t%get_filled', itimer=itimer)
  select type (target)
  class is (patch_t)
    if (associated(target%filled)) then
      if (any(target%filled == 0)) then
         call print_filled (self, target)
      end if
      !call target%lock%set ('halo_t%get_filled')
      !deallocate (target%filled, target%outside)
      !nullify    (target%filled, target%outside)
      !call target%lock%unset ('halo_t%get_filled')
    end if
  end select
  call trace%end (itimer)
END SUBROUTINE get_filled

!===============================================================================
!> Check that target is filled
!===============================================================================
SUBROUTINE print_filled (self, target)
  class(halo_t)          :: self
  class(patch_t)         :: target
  class(nbor_t), pointer :: nbor
  class(task_t), pointer :: task
  integer                :: i1, i2, i3, id_debug
  integer, parameter     :: m_max=2000
  integer, save          :: n_max=m_max
  !-----------------------------------------------------------------------------
  ! Log the task + list of nbor tasks and their boundary flags
  !-----------------------------------------------------------------------------
  write (io_unit%mpi,1) &
    'halo_t%print_filled: not filled, task,key,bndries=', target%id, target%key, &
    target%mesh(1)%lower_boundary, target%mesh(1)%upper_boundary, &
    target%mesh(2)%lower_boundary, target%mesh(2)%upper_boundary, &
    target%mesh(3)%lower_boundary, target%mesh(3)%upper_boundary
  1 format(a,i6,":",i2.2,3i4,3(2x,2l1))
  nbor => target%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (patch_t)
      write (io_unit%mpi,1) &
        'halo_t%print_filled: nbtask,key,bdries=', task%id, task%key, &
        task%mesh(1)%lower_boundary, task%mesh(1)%upper_boundary, &
        task%mesh(2)%lower_boundary, task%mesh(2)%upper_boundary, &
        task%mesh(3)%lower_boundary, task%mesh(3)%upper_boundary
    end select
    nbor => nbor%next
  end do
  if (n_max < 0) &
    return
  do i3=1,target%gn(3)
  do i2=1,target%gn(2)
  do i1=1,target%gn(1)
    if (n_max == 0) then
      write (stdlog,*) 'print_filled reached m_max=', m_max
      !$omp atomic write
      n_max = -1
      return
    end if
    if (target%filled(i1,i2,i3)==0) then
      write (io_unit%mpi,'(a,2(2x,4i4))') &
        'halo_t%print_filled: key,ii=', target%key, i1, i2, i3
      !$omp atomic update
      n_max = n_max-1
    end if
  end do
  end do
  end do
END SUBROUTINE print_filled

!===============================================================================
!> Check that target is filled
!===============================================================================
FUNCTION not_filled (self, target)
  class(halo_t)          :: self
  class(task_t), pointer :: target
  logical                :: not_filled
  integer                :: i1, i2, i3
  integer, save          :: n_max=50, itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%not_filled', itimer=itimer)
  not_filled = .false.
  if (check_filled .and. &
      target%is_clear(bits%root_grid+bits%no_filled)) then
    select type (target)
    class is (patch_t)
      if (associated(target%filled)) then
        not_filled = any(target%filled == 0)
        if (not_filled) then
          if (verbose > 0) &
            write (stdout,'(a,i6,i5,2x,3i3)') &
              'id,filled:', target%id, &
              minval(target%filled(:,:,4:17)), &
              minloc(target%filled(:,:,4:17)) + [0,0,3]
          call print_filled (self, target)
        end if
        !.......................................................................
      else
        not_filled = .false.
      end if
    end select
  end if
  call trace%end (itimer)
END FUNCTION not_filled

!===============================================================================
!> Overloaded procedure
!===============================================================================
SUBROUTINE alloc_link (task, link)
  class(task_t), target:: task
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  if (.not. associated(link)) then
    allocate (halo_t:: link)
    call os%mem_register (storage_size(link), 1, 'alloc_link')
    link%task => task
  end if
END SUBROUTINE alloc_link

!===============================================================================
!> Overloaded procedure
!===============================================================================
SUBROUTINE alloc_nbor (link, nbor, label)
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (guards_t:: nbor)
  call os%mem_register (storage_size(nbor), 1, shared%option(label, 'nbor'))
  call nbor%init ()
  nbor%link => link
  nbor%task => link%task
END SUBROUTINE alloc_nbor

SUBROUTINE check_nbor_list (self)
  class(halo_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  logical:: test
  !-----------------------------------------------------------------------------
  if (self%task%id == io%id_debug) then
    task => self%task
    select type (task)
    class is (patch_t)
      nbor => self%nbor
      do while (associated(nbor))
        nbtask => nbor%task
        select type (nbtask)
        class is (patch_t)
          test = task%overlaps2 (nbtask)
        end select
        nbor => nbor%next
      end do
    end select
  end if
END SUBROUTINE check_nbor_list

END MODULE halo_mod
