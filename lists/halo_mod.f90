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
    procedure:: make_coarse
    !procedure:: init_child_nbors
    procedure:: purge_nbors
    procedure:: conditionally_add_nbor
    procedure:: no_child_supports
    procedure:: info
    procedure:: set_filled
    procedure:: get_filled
    procedure:: check_support
  end type
  logical, save:: check_filled = .false.
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
    write (io%namelist, halo_params)
    call garbage%init_verbose (verbose)
    first_time = .false.
  end if
  call io%lock%unset ('halo_t%init')
  call face%pre_init
  nbor => self%nbor
  do while (associated(nbor))
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Walk through the nbor list and update guard zones.  The call will go to
!> the guards_t subclass of nbor_t, if that was used when allocating the nbor
!> list
!===============================================================================
SUBROUTINE update (self, all, only)
  implicit none
  class(halo_t)           :: self
  class(nbor_t), pointer  :: nbor
  logical, optional       :: all
  integer, optional       :: only
  !.............................................................................
  class(task_t), pointer  :: target, source, task
  class(patch_t), pointer :: coarse, patch
  class(mesh_t), pointer  :: mesh
  real(8)                 :: wc(2)
  integer                 :: i, m(6)
  logical                 :: use_coarse
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%update',itimer=itimer)
  self%n_op = 0
  !-----------------------------------------------------------------------------
  ! If istep is zero the task has just been created, including guard zones
  !-----------------------------------------------------------------------------
  if (self%task%istep > 0) then
    !---------------------------------------------------------------------------
    ! Check if a coarse update is needed, and make sure nbor%distance is set
    !---------------------------------------------------------------------------
    target => self%task
    use_coarse = .false.
    select type (target)
      class is (patch_t)
      if (target%id == io%id_debug) then
        print *,'id_debug=', io%id_debug
      end if
      call self%set_filled (target)
      !-------------------------------------------------------------------------
      ! Check if we need a coarse level clone; use this for all cases where
      ! there are level L-1 nbor that are NOT the parent.
      !-------------------------------------------------------------------------
      nbor => self%nbor
      do while (associated(nbor))
        if (nbor%task%level < self%task%level) then
          if (nbor%task%id /= target%parent%id) then
            use_coarse = .true.
          end if
        end if
        source => nbor%task
        select type (source)
        class is (patch_t)
          nbor%distance = source%distance (target)
          !.....................................................................
          if (face%verbose > 0) then
            write (stdlog,'(f12.6,a,2(2x,i6,":",i2.2),3x3l1,2x,4i3)') &
              wallclock(), pad_trace(1:indent)// &
              'halo_t%update: id,level,nbor,level,NMD,key=', &
              self%task%id, self%task%level, &
              nbor%task%id, nbor%task%level, &
              nbor%needed, nbor%needs_me, nbor%download, &
              source%hash_key()
          end if
          !.....................................................................
        end select
        nbor => nbor%next
      end do
    end select
    !---------------------------------------------------------------------------
    ! If a coarse patch clone is needed, allocate it, update its guard zones
    ! with a call to guards_r%update, and restrict the child patch to it.
    !---------------------------------------------------------------------------
    if (use_coarse) then
      call self%make_coarse ()
    else
      nbor => self%nbor
      do while (associated(nbor))
        if (nbor%download) then
          call face%source_to_target (nbor%task, self%task, n_op=self%n_op, &
                                      only=only)
        end if
        nbor => nbor%next
      end do
    end if
    call self%get_filled (self%task)
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,3i5)') &
        wallclock(), pad_trace(1:indent)// &
        'halo_t%update: n_same,prolong,restrict=', self%n_op
      self%n_op = 0
    end if
  else
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,i6)') &
        wallclock(), pad_trace(1:indent)// &
        'halo_t%update: no need=', self%task%id
      self%n_op = 0
    end if
  end if
  !-----------------------------------------------------------------------------
  ! With verbose>0, overrule the guards_t%update with the legacy methods (FIXME)
  !-----------------------------------------------------------------------------
  if (verbose > 0) then
    if (self%first_time) then
      self%first_time = .false.
      task => self%task
      select type (task)
      class is (patch_t)
      patch => task
      end select
      call patch%stats ('halo1')
      call self%task%dump (1, ad_hoc=.true., it=self%task%it)
      call self%link_t%update ()
      call patch%stats ('halo2')
      call self%task%dump (2, ad_hoc=.true., it=self%task%it)
    end if
  else
    call self%task%bndries
    call self%task%periodic_grid()
  end if
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Create a clone of the coarse parent patch a level L-1, and fill its guard
!> zone from patches in the child (self) nbor list a level L-1; these will
!> all be at most interpolations in time and copies in space.  The completed
!> coarse patch is able to supply all guard zone values for the level L task
!> by prolongation.
!===============================================================================
SUBROUTINE make_coarse (self)
  class(halo_t)           :: self
  class(nbor_t), pointer  :: nbor
  class(task_t), pointer  :: target, source
  class(patch_t), pointer :: coarse
  class(mesh_t), pointer  :: mesh
  real(8)                 :: wc(2)
  integer                 :: i, m(6)
  logical                 :: use_coarse
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%make_clone', itimer=itimer)
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
    write (stdlog,'(f12.6,a,i6,":",i2.2,2i6)') wallclock(), &
      pad_trace(1:indent)//'halo_t%make_coarse: id,lvl,verb,id_dbg=', &
      source%id, source%level, verbose, io%id_debug
  !.............................................................................
  end if
  select type (source)
  class is (patch_t)
    allocate (coarse, source=source)
  end select
  m = shape(coarse%mem)
  nullify (coarse%mem,coarse%t,coarse%iit)
  allocate (coarse%mem(m(1),m(2),m(3),m(4),1,1), coarse%t(1), coarse%iit(1))
  coarse%it = 1
  coarse%nt = 1
  coarse%iit = 1
  coarse%time = self%task%time
  coarse%t = coarse%time
  !---------------------------------------------------------------------------
  ! Run through the child nbor list and update the clone from the nbors;
  ! taking care to only use the nbors that are valid for the child, and also
  ! taking care not to use any other higher level nbors.  One of the child
  ! nbors is the parent, and this will take care of the copying of its values,
  ! but in principle only if all_tgt==.true (but this is the default behavior
  ! of the %same procedure).  The parent task is the only one that can fill
  ! the interior cells.
  !---------------------------------------------------------------------------
  nbor => self%nbor
  target => coarse
  do while (associated(nbor))
    if (face%verbose > 1) then
      write (stdlog,'(f12.6,a,i6,":",i2.2,l3)') wallclock(), &
        pad_trace(1:indent)// &
        'halo_t%make_coarse: using nbor: id,level,download=', &
        nbor%task%id, nbor%task%level, nbor%download
    end if
    if (nbor%download) then
      if (nbor%task%level == self%task%parent%level) then
        call face%source_to_target (nbor%task, target, all_tgt=.true., &
                                    n_op=self%n_op)
      else
        if (face%verbose > 1) then
          write (stdlog,'(a,i6,":",i2.2,l3)') &
            'skipping nbor: id,level,download=', &
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
  if (face%verbose > 1) then
    write (stdlog,'(a,i6,":",i2.2,l3)') 'restricting self: id,level=', &
      self%task%id, self%task%level
  end if
  target => coarse
  call face%source_to_target (self%task, target, n_op=self%n_op)
  !---------------------------------------------------------------------------
  ! Update the child guard-zones from the clone (via prolongation), then
  ! delete it.  To override the level conditions, mark target unfilled.
  ! All guard zones are to be filled, but not the interior, using all source
  !---------------------------------------------------------------------------
  source => coarse            ! the coarse time is identical to self%time
  target => self%task         ! the self%task to fill guard zones in
  select type (target)
  class is (patch_t)
    target%filled = 0
  end select
  call face%source_to_target (source, target, all_src=.true., n_op=self%n_op)
  deallocate (coarse%mem, coarse%t, coarse%iit)
  deallocate (coarse)
  call trace%end (itimer)
END SUBROUTINE make_coarse

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
!> Set face%download to false on nbors that are not needed, but keep them for
!> other purposes.
!===============================================================================
SUBROUTINE purge_nbors (self)
  class(halo_t)          :: self
  class(nbor_t), pointer :: nbor
  class(task_t), pointer :: task, nbtask
  integer, allocatable   :: filled(:,:,:)
  integer, dimension(3)  :: l, u, ls, us, m
  integer                :: n(2)
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  ! We use the preliminary nbor list, sorted by level, to add values in filled
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%purge_nbors', itimer=itimer)
  task => self%task
  select type (task)
  class is (patch_t)
  m = task%gn
  allocate (filled(m(1),m(2),m(3)))
  filled = task%level
  nbor => self%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (patch_t)
    call task%overlap (l, u, nbtask, ls, us, all_tgt=.true.)
    end select
    filled(l(1):u(1),l(2):u(2),l(3):u(3)) = nbtask%level
    select type (nbor)
    class is (guards_t)
    nbor%lt = l
    nbor%ut = u
    end select
    nbor => nbor%next
  end do
  end select
  !-----------------------------------------------------------------------------
  ! Run through again, and disable download from tasks that do not match filled
  !-----------------------------------------------------------------------------
  n = 0
  nbor => self%nbor
  do while (associated(nbor))
    select type (nbor)
    class is (guards_t)
    l = nbor%lt
    u = nbor%ut
    end select
    if (all(filled(l(1):u(1),l(2):u(2),l(3):u(3)) /= nbor%task%level)) then
      nbor%download = .false.
      nbor%needed = .false.
      n(1) = n(1) + 1
      !.........................................................................
      if (verbose > 1) then
        write (stdlog,'(f12.6,a,i6,":",i2.2)') wallclock(), &
          pad_trace(1:indent)//'update_nbors: disabled', nbtask%id, nbtask%level
      end if
      !.........................................................................
    else
      n(2) = n(2) + 1
    end if
    nbor => nbor%next
  end do
  deallocate (filled)
  !.............................................................................
  if (verbose > 0) then
    write (stdlog,'(f12.6,a,i6,":",i2.2)') wallclock(), &
      pad_trace(1:indent)//'update_nbors: disabled', nbtask%id, nbtask%level
  end if
  !.............................................................................
  call trace%end (itimer)
END SUBROUTINE purge_nbors

!===============================================================================
!> Check if the nbor should actually be on the nbor list, and if so try to add
!> it, in increasing level order.  If the task is already on the nbor list, the
!> attempt will be quietly ignored.
!===============================================================================
SUBROUTINE conditionally_add_nbor (self, nblink)
  class(halo_t)          :: self                    ! the task nblink
  class(link_t), pointer :: nblink                    ! the nbor task nblink
  !.............................................................................
  class(nbor_t), pointer :: nbor
  class(task_t), pointer :: task, nbtask, patch
  integer                :: n(4), added
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('halo_t%conditionally_add_nbor',itimer=itimer)
  n = 0
  task   => self%task
  patch  => self%task
  select type (patch)
  class is (patch_t)
    select type (nblink)
    class is (halo_t)
      nbtask => nblink%task
      select type (nbtask)
      class is (patch_t)
        n(1) = n(1) + 1
        !--- make sure it is not pending remove ---
        if (nbtask%is_clear(bits%remove)) then
          !--- does it have a chance to overlap? ---
          n(2) = n(2) + 1
          if (patch%overlaps_region (nbtask)) then
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
              write (stderr,*) 'WARNING'
              call nbor%init (task)
            end select
            if (verbose > 0) then
              write (stdlog,'(f12.6,a,2i6,":",i2.2,2(1x,3f7.3),2(1x,4i4))') &
                wallclock(), pad_trace(1:indent)// &
                'halo_t%cnd_add_nbor: tgt,src,dist,size,n=', &
                task%id, nbor%task%id, nbor%task%level, &
                nbor%distance/task%size, &
                nbtask%distance(patch)/task%size, &
                patch%hash_key(), nbtask%hash_key()
            end if
            call self%add_nbor_by_level (self%nbor, nbor, added)
          end if
        end if
      end select
    end select
  end select
  n_tot = n_tot + n
  if (detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE conditionally_add_nbor

!===============================================================================
!> Check if any of my children provides support to this nbor
!===============================================================================
FUNCTION no_child_supports (self, this) RESULT (ok)
  class(halo_t):: self, this
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  logical:: ok
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%no_child_supports')
  nbor => self%nbor
  ok = .true.
  do while (associated(nbor))
    if (nbor%task%level == this%task%level+1) then
      if (nbor%task%overlaps (this%task)) then
        ok = .false.
        exit
      end if
    end if
    nbor => nbor%next
  end do
  call trace%end ()
END FUNCTION no_child_supports

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
SUBROUTINE set_filled (self, target)
  class(halo_t)           :: self
  class(patch_t)          :: target
  class(mesh_t), pointer  :: mesh
  real(8)                 :: wc(2)
  integer                 :: i1, i2, i3
  logical                 :: outside(3)
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. check_filled) return
  call trace%begin ('halo_t%set_filled', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Clear target%filled and initialize internal points to target%level.
  ! Se logical mask 'outside' true in guard zones, false inside.
  !-----------------------------------------------------------------------------
  if (.not. associated(target%filled)) then
    allocate(target%filled (target%gn(1),target%gn(2),target%gn(3)))
    allocate(target%outside(target%gn(1),target%gn(2),target%gn(3)))
  end if
  !-----------------------------------------------------------------------------
  ! Test different methods to initialize; bell whether of optimization
  !-----------------------------------------------------------------------------
  if (face%verbose > 0) wc(1) = wallclock()
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
  if (face%verbose > 0) wc(2) = wallclock()
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
  if (face%verbose > 0) then
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
  integer, save          :: n_max=50, itimer=0
  !-----------------------------------------------------------------------------
  if (.not. check_filled) return
  call trace%begin ('halo_t%get_filled', itimer=itimer)
  select type (target)
  class is (patch_t)
    if (any(target%filled == 0)) then
      write (stdlog,'(a,i6)') &
        'guards_t%check_filled: target not filled =', target%id
      do i3=1,target%gn(3)
      do i2=1,target%gn(2)
      do i1=1,target%gn(1)
        if (n_max==0) exit
        if (target%filled(i1,i2,i3)==0) then
          write (stdlog,'(a,i6,2x,3i3)') &
            'guards_t%check_filled: id,ii=', target%id, i1, i2, i3
          !$omp atomic update
          n_max = n_max-1
        end if
      end do
      end do
      end do
      call io%abort ('guards_t%check_filled: some target not filled')
      flush(stdlog)
    end if
    deallocate (target%filled, target%outside)
    nullify    (target%filled, target%outside)
  end select
  call trace%end (itimer)
END SUBROUTINE get_filled

!===============================================================================
!> Check if a halo task has support at level L-1, by finding all L-2 nbors
!> (using the parent nbor list), and checking if any of the L-1 nbors covers
!> the overlap region.  If lack of support is found, trigger refinement.
!===============================================================================
SUBROUTINE check_support (self)
  class(halo_t)          :: self
  class(nbor_t), pointer :: nbors
  !.............................................................................
  type(region_t)         :: region1, region2, task_interest
  class(nbor_t), pointer :: nbor1, nbor2
  class(task_t), pointer :: task, nb1task, nb2task
  integer, dimension(3)  :: lt, ut, l, u
  real(8), dimension(3)  :: d
  logical                :: request_refine
  logical                :: has_support, covers, overlaps1, overlaps2
  integer, allocatable   :: request(:,:,:)
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  if (self%task%level < shared%levelmin + 2) then
    return
  end if
  call trace%begin ('check_support', itimer=itimer)
  task => self%task
  select type (task)
  class is (patch_t)
    !---------------------------------------------------------------------------
    ! Loop over the parent nbors, selecting those at L-2
    !---------------------------------------------------------------------------
    has_support = .true.
    nbor2 => self%parent%nbor
    do while (associated(nbor2))
      if (nbor2%task%level == task%level-2) then
        nb2task => nbor2%task
        select type (nb2task)
        class is (patch_t)
        task_interest = task%ROI (nb2task)
        overlaps2 = nb2task%authority%overlaps (task_interest)
        !-----------------------------------------------------------------------
        ! Select nbor with overlapping regions
        !-----------------------------------------------------------------------
        if (verbose > 1 .or. task%id == io%id_debug) then
          write (stdlog,'(f12.6,a,2i6,l3,2x,4i3)') wallclock(), &
            pad_trace(1:indent)//'check_support: tgt,L-2,ovrl,dist=', &
            task%id, nb2task%id, overlaps2, &
            nb2task%hash_key()
        end if
        !.......................................................................
        if (overlaps2) then
          region2 = nbor2%task%authority%intersection (task_interest)
          request_refine = .true.
          allocate(request(nb2task%gn(1),nb2task%gn(3),nb2task%gn(3)))
          request = 0
          l = ceiling(region2%l/nb2task%ds + nb2task%offset)
          u = floor  (region2%u/nb2task%ds + nb2task%offset)
          l = max(l,nb2task%li)
          u = min(u,nb2task%ui)
          request(l(1):u(1),l(2):u(2),l(3):u(3)) = task%id
          covers = .false.
!write (stdlog,'(a,i6,2x,4i3,2(2x,3i3),i7)') 'request: id,l,u,max=', &
!nb2task%id, nb2task%hash_key(), l, u, maxval(request)
          !---------------------------------------------------------------------
          ! Using the child nbor list here. The one L-1 nbor that does (or
          ! should) cover over the L-2 nbor will also be on the child nbor
          ! list.
          !---------------------------------------------------------------------
          nbor1 => self%nbor
          do while (associated(nbor1))
            nb1task => nbor1%task
            select type (nb1task)
            class is (patch_t)
            if (nb1task%level == task%level-1) then
              region1 = nbor1%task%authority%intersection (task_interest)
              !-----------------------------------------------------------------
              ! If there is some overlap, then chip away at the region
              !-----------------------------------------------------------------
              overlaps1 = region1%overlaps(region2)
              if (overlaps1) then
                !---------------------------------------------------------------
                ! These are the points the L-1 nbor can contribute
                !---------------------------------------------------------------
                l = ceiling(region1%l/nb2task%ds + nb2task%offset)
                u = floor  (region1%u/nb2task%ds + nb2task%offset)
                l = max(l,nb2task%li)
                u = min(u,nb2task%ui)
                request(l(1):u(1),l(2):u(2),l(3):u(3)) = 0
                covers = maxval(request) == 0
!write (stdlog,'(a,i6,2x,4i3,2(2x,3i3),l3)') '  cover: id,l,u=', &
!nb1task%id, nb1task%hash_key(), l, u, covers
              end if
              if (verbose > 1 .or. task%id == io%id_debug) then
                write (stdlog,'(f12.6,a,2i6,2x,4i3,2l3)') wallclock(), &
                  pad_trace(1:indent)//'check_support: tgt,L-1,vol=', &
                  task%id, nb1task%id, nb1task%hash_key(), overlaps1, covers
              end if
              !-----------------------------------------------------------------
              ! If nothing remains there is support
              !-----------------------------------------------------------------
              if (covers) then
                request_refine = .false.
                exit
              end if
            end if
            end select
            nbor1 => nbor1%next
          end do
          !---------------------------------------------------------------------
          ! nbor1 is not covered, so set the no_support bit, and make a
          ! refine_request entry in the nb2task%request array
          !---------------------------------------------------------------------
          if (request_refine) then
            has_support = .false.
            call task%set (bits%no_support)
            !-------------------------------------------------------------------
            ! Compute the lower and uppper index corners in the L-2 nb2task
            ! of the region formed as a intersection of that ones authority and the
            ! level L task's region-of-interest, computed in the nbtask frame
            !-------------------------------------------------------------------
            call nb2task%request_lock%set ('has_support')
            nb2task%request = nb2task%request + request
            call nb2task%request_lock%unset ('has_support')
            if (verbose > 0 .or. task%id == io%id_debug) then
              write (stdlog,'(f12.6,a,2(i6,":",i2.2),2(2x,3i3))') &
                wallclock(), pad_trace(1:indent)//&
                'halo_t%check_support: request id,nbid,l,u=', &
                task%id, task%level, nb2task%id, nb2task%level, l, u
            end if
          end if
          deallocate (request)
        end if
        end select
      end if
      nbor2 => nbor2%next
    end do
  end select
  !-----------------------------------------------------------------------------
  ! If the task now has support, clear the no_support bit
  !-----------------------------------------------------------------------------
  if (has_support) then
    call task%clear (bits%no_support)
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,i6,":",i2.2)') &
        wallclock(), pad_trace(1:indent)//&
        'halo_t%check_support: clear no_support id=', &
        task%id, task%level
    end if
  else if (task%is_clear (bits%no_support)) then
    call task%set (bits%no_support)
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,i6,":",i2.2)') &
        wallclock(), pad_trace(1:indent)//&
        'halo_t%check_support: set no_support id=', &
        task%id, task%level
    end if
  else
    call task%set (bits%no_support)
    if (verbose > 0) then
      write (stdlog,'(f12.6,a,i6,":",i2.2)') &
        wallclock(), pad_trace(1:indent)//&
        'halo_t%check_support: has no_support id=', &
        task%id, task%level
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE check_support

END MODULE halo_mod
