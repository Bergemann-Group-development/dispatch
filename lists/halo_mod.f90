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
  USE timer_mod
  USE trace_mod
  USE bits_mod
  USE task_mod
  USE mesh_mod
  USE patch_mod
  USE link_mod
  USE guards_mod
  USE region_mod
  USE shared_mod
  implicit none
  private
  !.............................................................................
  type, public, extends(link_t):: halo_t
    class(patch_t), pointer:: coarse => null()
    integer:: n_op(3)=0
    logical:: first_time=.true.
  contains
    procedure:: init
    procedure:: update
    procedure:: download
    procedure:: make_coarse
    !procedure:: init_child_nbors
    procedure:: conditionally_add_nbor
    procedure:: info
    procedure:: set_filled
    procedure:: get_filled
    procedure:: not_filled
    procedure:: has_support
    procedure:: request_support
    procedure:: check_support
    procedure:: find_support
  end type
  logical, save:: check_filled=.true.
  logical, save:: detailed_timer=.false.
  integer, save:: verbose=0
  integer, save:: method=3
  integer, save:: n_tot(4)=0
  real(8), save:: overlap_margin = 1d-10
  !$omp threadprivate (n_tot)
CONTAINS

!===============================================================================
!> Initialize the link, and possibly the nbor nodes
!===============================================================================
SUBROUTINE init (self)
  class(halo_t):: self
  !.............................................................................
  class(nbor_t), pointer:: nbor
  integer:: iostat
  namelist /halo_params/ verbose, detailed_timer, check_filled, method
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
  call io%lock%unset ('halo_t%init')
  call guards%pre_init
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
  call trace%begin ('halo_t%update',itimer=itimer)
  unit = merge (stdlog, stdout, verbose > 2)
  self%n_op = 0
  !-----------------------------------------------------------------------------
  ! If istep is zero the task has just been created, including guard zones
  !-----------------------------------------------------------------------------
  target => self%task
  if (target%time > target%dnload_time) then
    call self%download (self%nbor, only, all_tgt, all_src)
    if (verbose > 0) then
      write (unit,'(a,i6,2e15.6)') &
        'halo done', target%id, target%time, target%dnload_time
    end if
    target%dnload_time = target%time
  else
    if (verbose > 0) then
      write (unit,'(a,i6,2e15.6)') &
        'halo skip', target%id, target%time, target%dnload_time
    end if
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
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Comprehensive entry point
!===============================================================================
SUBROUTINE download (self, nbors, only, all_tgt, all_src)
  class(halo_t)          :: self
  class(nbor_t), pointer :: nbors
  integer,      optional :: only
  logical,      optional :: all_tgt, all_src
  !.............................................................................
  integer                :: unit
  logical                :: use_coarse
  class(nbor_t), pointer :: nbor
  class(task_t), pointer :: target, source
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  ! Decide whether to use a coarse clone or not
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%download', itimer=itimer)
  unit = merge (stdlog, stdout, verbose > 2)
  use_coarse = .false.
  target => self%task
  !-----------------------------------------------------------------------------
  ! If the task does not yet have support it will soon have it, and until then
  ! we make sure to fill all guard zones with values from the parent task, which
  ! has at least time-extrapolated guard zone values, and possibly already
  ! authoritative guard zone values
  !-----------------------------------------------------------------------------
  call self%set_filled (target)
  call self%check_support ()
  if (self%task%is_set (bits%no_support)) then
    source => self%parent%task
    call guards%source_to_target (source, target, n_op=self%n_op, only=only, &
                                  all_src=.true., all_tgt=all_tgt)  
  end if
  !-----------------------------------------------------------------------------
  ! Proceed with the normal guard zone loading
  !-----------------------------------------------------------------------------
  nbor => nbors
  do while (associated(nbor))
    if (nbor%download) then
      source => nbor%task
      if (verbose > 0) then
        write (stdout,'(a8,2x,2(a,i6,":",i2.2,1x,3i4,3x))') 'halo', &
          'tgt:', target%id, target%key, &
          'src:', source%id, source%key
      end if
      !...........................................................................
      if (source%level < target%level) then
        if (source%id /= target%parent%id) then
          use_coarse = .true.
        end if
      end if
    end if
    nbor => nbor%next
  end do
  !-----------------------------------------------------------------------------
  ! Pass on, or loop over nbors
  !-----------------------------------------------------------------------------
  if (use_coarse) then
    call make_coarse (self, only=only, nbors=nbors)
  else
    nbor => nbors
    do while (associated(nbor))
      if (nbor%download) then
        ! if (nbor%task%is_clear (bits%no_support)) then
          source => nbor%task
          call guards%source_to_target (source, target, n_op=self%n_op, &
            only=only, all_src=all_src, all_tgt=all_tgt)
        ! else
        !   print *, target%id, 'skipping', nbor%task%id
        ! end if
      end if
      nbor => nbor%next
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Boundary conditions
  !-----------------------------------------------------------------------------
  call self%task%bndries
  call self%task%periodic_grid()
  !-----------------------------------------------------------------------------
  ! Emergency procedure: If some cells are not filled; redo the guard zone
  ! loading, using full source patches, while printing warnings
  !-----------------------------------------------------------------------------
  if (self%not_filled (target)) then
    write (stderr,'(a,i6,":",i2.2,3i4)') &
      ' halo_t%not_filled: target=', target%id, target%key
    nbor => nbors
    do while (associated(nbor))
      source => nbor%task
      select type (source)
      class is (patch_t)
      write (stderr,'(i6,":",i2.2,3i4,3x,l1,2(2x,3l1))') source%id, source%key, &
          nbor%download, source%mesh%lower_boundary, source%mesh%upper_boundary
      end select
      nbor => nbor%next
    end do
    nbor => nbors
    do while (associated(nbor))
      if (nbor%download) then
        source => nbor%task
        select type (source)
        class is (patch_t)
          call source%extrapolate_guards
        end select
        call guards%source_to_target (source, target, n_op=self%n_op, &
          all_src=.true., all_tgt=all_tgt, only=only)
      end if
      nbor => nbor%next
    end do
    if (self%not_filled (target)) then
      write (stderr,*) 'halo_t%update: repair NOK, target%id=', target%id
    else
      write (stderr,*) 'halo_t%update: repair OK, target%id=', target%id
    end if
  end if
  call self%get_filled (target)
  !.............................................................................
  if (verbose > 0) then
    write (unit,'(f12.6,a,3i5)') &
      wallclock(), pad_trace(1:indent)// &
      'halo_t%download: n_same,prolong,restrict=', self%n_op
    self%n_op = 0
  end if
call trace%end (itimer)
END SUBROUTINE download

!===============================================================================
!> Create a clone of the coarse parent patch a level L-1, and fill its guard
!> zone from patches in the child (self) nbor list a level L-1; these will
!> all be at most interpolations in time and copies in space.  The completed
!> coarse patch is able to supply all guard zone values for the level L task
!> by prolongation.
!===============================================================================
SUBROUTINE make_coarse (self, nbors, only)
  class(halo_t)                    :: self
  integer      , optional          :: only
  class(nbor_t), optional, pointer :: nbors
  !.............................................................................
  class(nbor_t) , pointer :: nbor
  class(task_t) , pointer :: target, source, parent, task
  class(patch_t), pointer :: coarse
  class(mesh_t) , pointer :: mesh
  real(8)                 :: wc(2)
  integer                 :: i, m(6), unit
  logical                 :: use_coarse
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%make_coarse', itimer=itimer)
  unit = merge (stdlog, stdout, verbose > 2)
  if (verbose > 0) &
    print *, self%task%id, 'make coarse'
  self%n_op = 0
  !-----------------------------------------------------------------------------
  ! Allocate a duplicate of the coarse parent
  !-----------------------------------------------------------------------------
  target => self%task
  source => self%task%parent
  !.............................................................................
  if (verbose > 0 .or. &
      target%id == io%id_debug .or. &
      source%id == io%id_debug) then
    write (unit,'(f12.6,a,i6,":",i2.2,2i6)') wallclock(), &
      pad_trace(1:indent)//'halo_t%make_coarse: id,lvl,verb,id_dbg=', &
      source%id, source%level, verbose, io%id_debug
  !.............................................................................
  end if
  allocate (task, source=source)
  select type (task)
  class is (patch_t)
    coarse => task
  end select
  m = shape(coarse%mem)
  nullify (coarse%mem,coarse%t,coarse%iit)
  allocate (coarse%mem(m(1),m(2),m(3),m(4),1,1), coarse%t(1), coarse%iit(1))
  coarse%id = 0
  coarse%it = 1
  coarse%nt = 1
  coarse%iit = 1
  coarse%key = self%task%key
  coarse%key(1) = self%task%key(1) - 1
  coarse%time = self%task%time
  coarse%t = coarse%time
  !-----------------------------------------------------------------------------
  ! Run through the child nbor list and update the clone from the nbors;
  ! taking care to only use the nbors that are valid for the child, and also
  ! taking care not to use any other higher level nbors.  One of the child
  ! nbors is the parent, and this will take care of the copying of its values,
  ! but in principle only if all_tgt==.true (but this is the default behavior
  ! of the %same procedure).  The parent task is the only one that can fill
  ! the interior cells.
  !-----------------------------------------------------------------------------
  if (present(nbors)) then
    nbor => nbors
  else
    nbor => self%nbor
  end if
  do while (associated(nbor))
    if (nbor%download) then
      if (nbor%task%level == self%task%parent%level) then
        if (guards%verbose > 1 .or. self%task%id == io%id_debug) then
          write (unit,'(a,2(i6,":",i2.2),l3)') &
            'halo_t%make_coarse:    using tgt,src,download=', &
            self%task%id, self%task%level, &
            nbor%task%id, nbor%task%level, nbor%download
        end if
        call guards%source_to_target (nbor%task, coarse, all_tgt=.true., &
                                      only=only, n_op=self%n_op)
      else
        if (guards%verbose > 1 .or. self%task%id == io%id_debug) then
          write (unit,'(a,2(i6,":",i2.2),l3)') &
            'halo_t%make_coarse: skipping tgt,src,download=', &
            self%task%id, self%task%level, &
            nbor%task%id, nbor%task%level, nbor%download
        end if
      end if
    end if
    nbor => nbor%next
  end do
  !---------------------------------------------------------------------------
  ! Restrict the level L target to the coarse clone (FIXME: strictly, only the
  ! "rim" values are needed)
  !---------------------------------------------------------------------------
  target => coarse
  call guards%source_to_target (self%task, target, n_op=self%n_op, &
                                only=only, all_tgt=.true.)
  !---------------------------------------------------------------------------
  ! Update the child guard-zones from the clone (via prolongation), then
  ! delete it.  To override the level conditions, mark target unfilled.
  ! All guard zones are to be filled, but not the interior, using all source.
  !---------------------------------------------------------------------------
  source => coarse            ! the coarse time is identical to self%time
  target => self%task         ! the self%task to fill guard zones in
  select type (target)
  class is (patch_t)
    if (check_filled .and. associated(target%filled)) then
      target%filled = 0
    end if
  end select
  call dbg (self%task, 'before s2t')
  call guards%source_to_target (source, target, n_op=self%n_op, &
                                only=only, all_src=.true.)
  call dbg (self%task, ' after s2t')
  deallocate (coarse%mem, coarse%t, coarse%iit)
  deallocate (coarse)
  call trace%end (itimer)
END SUBROUTINE make_coarse

SUBROUTINE dbg (target, label)
  class(task_t):: target
  character(*):: label
  if (target%id == io%id_debug) then
    select type (target)
    class is (patch_t)
      call target%stats (label)
    end select
  end if
END SUBROUTINE

!===============================================================================
!> Initialize the nbor list of a new child task, by adding just the parent
!> task as the only nbor.  Self is the child link
!===============================================================================
SUBROUTINE init_child_nbors (self, child_link)
  class(halo_t):: self
  class(halo_t), pointer:: child_link
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  ! Add the parent as an nbor to child
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%init_child_nbors')
  allocate (guards_t:: nbor)
  nbor%link => self%parent
  nbor%task => self%parent%task
  call nbor%init (self%task)
  call self%add_nbor_by_level (self%nbor, nbor)
  self%task%n_needed = self%task%n_needed + 1
  call trace%end ()
END SUBROUTINE init_child_nbors

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
  class(nbor_t), pointer :: nbor
  class(task_t), pointer :: task, nbtask, patch
  integer                :: n(4)
  real(8)                :: distance(3)
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('halo_t%conditionally_add_nbor',itimer=itimer)
  n = 0
  if (self%task%rank == mpi%rank .or. nblink%task%rank == mpi%rank) then
    task   => self%task
    patch  => self%task
    select type (patch)
    class is (patch_t)
      nbtask => nblink%task
      select type (nbtask)
      class is (patch_t)
        n(1) = n(1) + 1
        !--- make sure it is not pending remove ---
        if (nbtask%is_clear(bits%remove)) then
          !--- does it have a chance to overlap? ---
          distance = task%distance_to (nbtask)
          n(2) = n(2) + 1
          if (patch%interest%overlaps (nbtask%authority, distance)) then
            n(3) = n(3) + 1
            if (nbtask%id == io%id_debug) write (stdlog,*) nbtask%id, 'overlaps_r'
            !--- if yes, allocate nbor associated with task, and add it ---
            allocate (guards_t:: nbor)
            select type(nbor)
            class is (guards_t)
              nbor%link => nblink
              nbor%task => nbtask
              call nbor%init (task)
            class default
              write (stderr,*) 'halo_t%conditionally_add_nbor: WARNING'
              call nbor%init (task)
            end select
            nbor%distance = distance
!            nbor%status = nbor%task%status
            if (verbose > 0) then
              write (stdlog,'(f12.6,a,2i6,":",i2.2,2(1x,3f7.3),2(1x,4i4))') &
                wallclock(), pad_trace(1:indent)// &
                'halo_t%cnd_add_nbor: tgt,src,dist,size,n=', &
                task%id, nbor%task%id, nbor%task%level, &
                nbor%distance/task%size, &
                nbtask%distance(patch)/task%size, &
                patch%hash_key(), nbtask%hash_key()
            end if
            call self%add_nbor_by_level (nbors, nbor, added)
          end if
        end if
      end select
    end select
    n_tot = n_tot + n
  end if
  if (detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE conditionally_add_nbor

!===============================================================================
!> Initialize the link, and possibly the nbor nodes
!===============================================================================
SUBROUTINE info (self, label)
  class(halo_t):: self
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%info')
  if (verbose < 2) then
    call self%link_t%info
  else
    if (present(label)) then
      write (stdlog,'(f12.6,a,10i6)') wallclock(), &
        pad_trace(1:indent)//label//': n_tot=', n_tot
      n_tot = 0
    else
      write (stdlog,'(f12.6,a,10i6)') wallclock(), &
        pad_trace(1:indent)//'halo_t%info: n_tot=', n_tot
      n_tot = 0
    end if
    call self%print_nbors (label='halo_t%info', verbose=verbose)
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
  logical                 :: outside(3)
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. check_filled) return
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
  if (guards%verbose > 0) wc(1) = wallclock()
  if (method == 1) then
    target%filled = 0
    target%outside = .true.
  else if (method == 2) then
    do i3=1,target%gn(3)
      do i2=1,target%gn(2)
        !$omp simd
        do i1=1,target%gn(1)
          target%outside(i1,i2,i3) = .true.
          target%filled(i1,i2,i3) = .false.
        end do
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Method 3 and 4 use a single loop nest, with logics inside
  !-----------------------------------------------------------------------------
  if (guards%verbose > 0) wc(2) = wallclock()
  if (method >= 3) then
    wc(1) = wc(2)
    do i3=1,target%gn(3)
      mesh => target%mesh(3)
      outside(3) = (i3 < mesh%li .or. i3 > mesh%ui)
      do i2=1,target%gn(2)
        mesh => target%mesh(2)
        outside(2) = (i2 < mesh%li .or. i2 > mesh%ui)
        !-----------------------------------------------------------------------
        ! Simple logics, easy to read
        !-----------------------------------------------------------------------
        if (method == 3) then
          !$omp simd
          do i1=1,target%gn(1)
            mesh => target%mesh(1)
            outside(1) = (i1 < mesh%li .or. i1 > mesh%ui)
            if (any(outside)) then
              target%outside(i1,i2,i3) = .true.
              target%filled(i1,i2,i3) = 0
            else
              target%outside(i1,i2,i3) = .false.
              target%filled(i1,i2,i3) = target%level
            end if
          end do
        !-----------------------------------------------------------------------
        ! Attempt to aid vectorization
        !-----------------------------------------------------------------------
        else
          mesh => target%mesh(1)
          !$omp simd
          do i1=1,target%gn(1)
            outside(1) = (i1 < mesh%li .or. i1 > mesh%ui)
            target%outside(i1,i2,i3) = outside(1).or.outside(2).or.outside(3)
            target%filled(i1,i2,i3) = merge(0,target%level,any(outside))
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
          target%outside(i1,i2,i3) = .false.
          target%filled(i1,i2,i3) = target%level
        end do
      end do
    end do
  end if
  if (guards%verbose > 0) then
    write (stdlog,'(f12.6,a,i3,9p,2f8.3)') &
      wallclock(), pad_trace(1:indent)// &
      'guards_t%set_filled: filled,outside method,ns/pt=', method, &
      (wc(2)-wc(1))/product(target%n), (wallclock()-wc(1))/product(target%n)
  end if
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
  if (.not. check_filled .or. target%is_set(bits%root_grid)) return
  call trace%begin ('halo_t%get_filled', itimer=itimer)
  select type (target)
  class is (patch_t)
    if (associated(target%filled)) then
      if (any(target%filled == 0)) then
         call print_filled (target)
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
SUBROUTINE print_filled (target)
  class(patch_t) :: target
  integer        :: i1, i2, i3
  integer, save  :: n_max=500
  !-----------------------------------------------------------------------------
  write (stdlog,'(a,i6)') &
    'guards_t%check_filled: target not filled =', target%id
  do i3=1,target%gn(3)
  do i2=1,target%gn(2)
  do i1=1,target%gn(1)
    if (n_max<=0) &
      return
    if (target%filled(i1,i2,i3)==0) then
      write (stdlog,'(a,i6,2x,3i3)') &
        'guards_t%print_filled: id,ii=', target%id, i1, i2, i3
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
  if (check_filled .and. target%is_clear(bits%root_grid)) then
    select type (target)
    class is (patch_t)
      if (associated(target%filled)) then
        not_filled = any(target%filled == 0)
        if (not_filled .and. verbose > 0) then
          call print_filled (target)
        end if
      else
        not_filled = .false.
      end if
    end select
  end if
  call trace%end (itimer)
END FUNCTION not_filled

!===============================================================================
!> Decide if a link has support.  This uses the %filled(:,:,:) array, which
!> means that this call must not overlap with actual guard zone filling.
!>
!> When this procedure discovers a task that lacks support it is run by the
!> thread that created the task, so no link locking is needed.  Later, when
!> it revisits the task and finds that it has support, if it is run by
!> the same thread that updates it, again no link locking is needed, but if
!> it is run by a thread that is not updating it, link locking is needed.
!===============================================================================
FUNCTION has_support (self)
  class(halo_t):: self
  logical:: has_support
  !.............................................................................
  class(task_t), pointer:: target, source
  class(nbor_t), pointer:: nbor, nbors
  real(8):: distance(3)
  !-----------------------------------------------------------------------------
  call self%lock%set ('download_link')
  call self%copy_nbor_list (self%nbor, nbors)
  call self%lock%unset ('download_link')
  target => self%task
  if (omp%thread_task > 0 .and. target%id /= omp%thread_task) then
    write (stderr,*) 'has_support WARNING: target%id,thread_task=', &
      target%id, omp%thread_task
  end if
  has_support = .true.
  if (target%level > shared%levelmin+1) then
    select type (target)
    class is (patch_t)
      call target%support_arrays ()
      target%filled(:,:,:) = 1
      distance = 0d0
      call target%cover (target, distance, target%filled)
      nbor => nbors
      do while (associated(nbor))
        !-----------------------------------------------------------------------
        ! In case L-2 tasks are included in the nbor list, skip them, since
        ! they are not to be used for downloads
        !-----------------------------------------------------------------------
        if (nbor%task%level >= target%level-1) then
          if (nbor%download) then
            source => nbor%task
            if (verbose > 0) &
            print '("has_support",2(i6,":",i2.2,3i3))', &
            target%id, target%key, source%id, source%key
            select type (source)
            class is (patch_t)
              call target%cover (source, nbor%distance, target%filled)
            end select
          end if
        end if
        nbor => nbor%next
      end do
      has_support = all(target%filled==0)
    end select
  end if
  call self%remove_nbor_list (nbors)
END FUNCTION has_support

!===============================================================================
!> Find the task at level L-2 that should be refined to provide support, and
!> mark the corresponging cells in its patch_t%request(:,:,:) array.  The
!> task must be among the level L-2 nbors of the target parent
!===============================================================================
SUBROUTINE request_support (self, source)
  class(halo_t):: self
  class(task_t), pointer:: source
  !.............................................................................
  class(task_t), pointer:: target
  integer:: i1, i2, i3, key(4), ii(3), n
  real(8):: point(3), distance(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%request_support', itimer=itimer)
  target => self%task
  select type (target)
  class is (patch_t)
    call target%support_arrays ()
    select type (source)
    class is (patch_t)
      call source%request_lock%set ('halo_t%request_support')
      n = 0
      do i3=1,target%gn(3)
      do i2=1,target%gn(2)
      do i1=1,target%gn(1)
        if (target%filled(i1,i2,i3) == 1) then
          !---------------------------------------------------------------------
          ! If the source contains the point, mark it, and clear the cell
          !---------------------------------------------------------------------
          point = target%local_position ([i1,i2,i3]) &
                + source%distance_to (target)
          if (source%authority%contains (point)) then
            ii = source%index_only (point, nearest=.true.)
            source%request(ii(1),ii(2),ii(3)) = 2.0
            target%filled(i1,i2,i3) = 0
            n = n+1
          end if
        end if
      end do
      end do
      end do
      call source%request_lock%unset ('halo_t%request_support')
      !print *, 'cleared', source%id, target%id, n, all(target%filled==0)
    end select
  end select
  call trace%end (itimer)
END SUBROUTINE request_support

!===============================================================================
!> Check if a link task has support, and set status accordingly
!===============================================================================
SUBROUTINE check_support (self)
  class(halo_t):: self
  class(task_t), pointer:: source
  class(nbor_t), pointer:: nbor
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. shared%refine_on) &
    return
  call trace%begin ('halo_t%check_support', itimer=itimer)
  if (self%has_support()) then
    if (self%task%is_set (bits%no_support)) then
      !.........................................................................
!print *, self%task%id, 'clear no_support', omp%thread
      !.........................................................................
      call self%task%clear (bits%no_support)
    end if
  else
    if (self%task%is_clear (bits%no_support)) then
      !.........................................................................
!print *, self%task%id, 'set no_support', omp%thread
      !.........................................................................
      call self%task%set (bits%no_support)
      call self%find_support (source)
      if (associated(source)) then
        call self%request_support (source)
!print *, self%task%id, 'request support from', source%id
      else
        write (stderr,*) 'halo_t%check_support WARNING: did not find support', &
          self%task%id, omp%thread, associated(self%nbor)
        nbor => self%nbor
        do while (associated(nbor))
          print *, nbor%task%id, nbor%task%level, nbor%distance/nbor%task%size
          nbor => nbor%next
        end do
      end if
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE check_support

!===============================================================================
!> Find a level L-2 task that overlaps
!===============================================================================
SUBROUTINE find_support (self, source)
  class(halo_t):: self
  class(task_t), pointer:: source, target
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%find_support')
  target => self%task
  nbor => self%nbor
  do while (associated(nbor))
    source => nbor%task
    if (source%level == self%task%level-2) then
      nbor%distance = target%distance_to (source)
      if (source%authority%overlaps (target%interest, nbor%distance)) then
        call trace%end ()
        return
      end if
    end if
    nbor => nbor%next
  end do
  nullify (source)
  call trace%end ()
END SUBROUTINE find_support

END MODULE halo_mod
