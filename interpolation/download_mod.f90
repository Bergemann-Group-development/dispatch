!===============================================================================
!>  download_link: takes care of downloads to link%task
!>           same: called for patches on the same level
!>         differ: called for patches on different levels
!>
!> Since we keep parent patches, there can be any number of levels overlapping
!> in the guard zones of a patch.  Therefore, to avoid wasting cycles on 
!> computing guard cell values for all levels, we sort the nbor list in order
!> of decreasing levels, and check if a cell is already filled before computing
!> a value from an nbor patch.
!>
!> The "filled" array stores the level of the value stored in the cell, and since
!> staggering can vary for different variables, there needs to be an array for
!> each variable.
!===============================================================================
MODULE download_mod
  USE iso_fortran_env, only: int8
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE bits_mod
  USE kinds_mod
  USE mesh_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE shared_mod
  USE interpolate_mod
  USE lagrange_mod
  USE validate_mod
  USE omp_lock_mod
  USE timer_mod
  USE guard_zones_mod
  USE remesh_mod
  implicit none
  private
  type, public:: download_t
    integer(8):: n_same=0, n_differ=0
    integer:: order_interpolator
    logical:: check_filled=.false.
    logical:: link_locks=.false.
    procedure(interpolator_interface), pointer, nopass:: interpolate => null()
    procedure(same_linear)           , pointer        :: same        => null()
  contains
    procedure:: init
    procedure:: download_link
    procedure, nopass:: dnload_range
    procedure:: prolong
    procedure:: different
    procedure:: different_prolong
    procedure:: test
    procedure, nopass:: set_verbose
  end type
  integer(kind=int8), dimension(:,:,:,:), allocatable:: filled
  integer(kind=4), dimension(:,:,:), allocatable:: src
  !$omp threadprivate (filled, src)
  integer:: verbose=0
  integer:: id_debug=0
  integer:: use_different=0
  logical:: sorted=.false.
  logical:: use_locks=.true.
  logical:: detailed_timer=.false.
  character(len=16):: use_restrict='', use_prolong='', use_same=''
  type(download_t), public:: download
CONTAINS

!===============================================================================
!> Initialise downloading and interpolating options.
!===============================================================================
SUBROUTINE init (self)
  class(download_t):: self
  integer:: iostat
  integer, save:: order_interpolator=-1, order_time=1
  logical, save:: check_filled=.false.
  logical, save:: first_time=.true.
  namelist /download_params/ order_interpolator, order_time, check_filled, &
    use_prolong, use_restrict, use_same, sorted, use_locks, id_debug, &
    detailed_timer, verbose
  character(len=120):: id = &
    '$Id$ interpolation/download_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%begin ('download_t%init')
  call trace%print_id (id)
  !$omp critical (download_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, download_params, iostat=iostat)
    if (io%master) write (*, download_params)
    call guard_zones%init
    call remesh%init
    call lagrange%init
  end if
  if (order_time==1) then
    self%same => same_linear
  else
    self%same => same_lagrange
  end if
  self%check_filled = check_filled
  self%order_interpolator = order_interpolator
  self%interpolate => SelectInterpolator(order_interpolator)
  !$omp end critical (download_cr)
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Select the best method for each case.  So far, this only supports static
!> patches at the same level, or static patches one level up or down
!===============================================================================
SUBROUTINE download_link (self, link, all_cells, only, all_times)
  class(download_t):: self
  class(link_t), pointer:: link
  logical, optional:: all_cells, all_times
  integer, optional:: only
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  !.............................................................................
  class(link_t), pointer:: nbor, sorted_head, nbors
  class(patch_t), pointer:: target, source
  integer:: edges, l(3), u(3), ix, iy, iz, jt(2)
  integer, save:: itimer=0
  real:: pt(2)
  character(len=64):: info
  logical:: all_cells_l, lock_it
  integer:: n_not
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%download_link', itimer=itimer)
  if (.not.associated(link)) then
    call io%abort('download_t%download_link: link not associated')
  else if (.not.associated(link%task)) then
    call io%abort('download_t%download_link: link%task not associated')
  end if
  !-----------------------------------------------------------------------------
  ! Guard against download twice, which could happen if an extras (e.g. RT) task
  ! requests dnload before the parent task
  !-----------------------------------------------------------------------------
  target => task2patch(link%task)
  call target%log ('download', 2)
  if (target%dnload_time >= target%time .or. &
      target%is_set(bits%no_download)) then
    call trace%end (itimer)
    return
  end if
  target%dnload_time = target%time
  if (target%is_clear(bits%no_density)) then
    d => target%mem(:,:,:,target%idx%d,target%it,1)
    call validate%check (target%link, d, 'before download')
    if (target%id == id_debug .or. verbose > 0) then
      write(io_unit%log,*) 'download_t%download_link: id, time =', &
        target%id, target%time
    end if
  end if
  if (.not.associated(link%nbor)) then
    sorted = .false.
  end if
  !-----------------------------------------------------------------------------
  ! Allocate and initialise array to track filling of cells.
  ! A value of `-1` defines unfilled; all other values denote the *level* of the
  ! patch that filled that cell.
  !-----------------------------------------------------------------------------
  if (present(all_cells)) then
    all_cells_l = all_cells
  else if (target%quality == 1024) then
    all_cells_l = .true.
  else
    all_cells_l = .false.
  end if
  if (self%check_filled) then
    allocate (filled(target%gn(1),target%gn(2),target%gn(3),target%nv))
    if (verbose > 0) &
      allocate (src(target%gn(1),target%gn(2),target%gn(3)))
    l = target%li
    u = target%ui
    filled = -1
    if (.not.all_cells_l) then
      filled(l(1):u(1),l(2):u(2),l(3):u(3),:) = target%level
      if (verbose > 0) &
        src(l(1):u(1),l(2):u(2),l(3):u(3)) = target%id
    end if
  end if
  if (target%is_clear(bits%no_density) .and. .not.all_cells_l) &
    call target%check_density ("before download")
  if (target%is_periodic()) call target%periodic_grid
  !-----------------------------------------------------------------------------
  ! Use an nbor list sorted by decreasing level, recording in each cell the
  ! level that supplied values, and skipping nbor contributions from lower levels.
  ! To avoid dead-lock, make a copy of it, and release the lock
  !-----------------------------------------------------------------------------
  nbor => link%nbors_by_level
  !-----------------------------------------------------------------------------
  ! Insert a value only if "filled" indicates that one is missing.  With the 
  ! nbor list sorted in decreasing level order, values are only computed where
  ! actually needed
  !-----------------------------------------------------------------------------
  do while (associated(nbor))
    if (.not.associated(nbor%task)) then
      write (io_unit%log,*) 'WARNING: nbor task not associated'
      flush (io_unit%log)
      nbor => nbor%next
      cycle
    end if
    source => task2patch(nbor%task)
    call source%log ('download', 3)
    if (source%quality == 1024) then
      all_cells_l = .true.
    end if
    if (verbose > 1) then
      write(io_unit%log,'(f12.6,2x,a,2(i6,i3),2x,2l1,i4)') &
        wallclock(), 'target, level, source, level, ND, n_needed =', &
        target%id, target%level, source%id, source%level, &
        nbor%needed, nbor%download, source%n_needed
      flush (io_unit%log)
    end if
    !-------------------------------------------------------------------------
    ! If the source task is being actively updated, and the slot could become
    ! the one after source%new, lock the task mem, so either it cannot change
    ! during the guard zone loading or, if the lock is already set, access to
    ! the source will have to wait until after the source update.
    !-------------------------------------------------------------------------
    call source%time_interval (target%time, jt, pt)
    call source%mem_lock(jt(1))%set ('download_link1')
    if (jt(2) /= jt(1)) &
      call source%mem_lock(jt(2))%set ('download_link2')
    if (verbose > 0) &
      write (io_unit%log,'(a,2i3,1x,2i3,2x,2f6.3,f12.6,2x,10f12.6)') &
        'download_link: it, new, jt, pt, target%time, source%t =', &
        source%it, source%new, jt, pt, target%time, source%t
    if (verbose > 1) then
      write(io_unit%log,'(f12.6,i6,i3,3x,a,i6,i3,2x,2l1)') wallclock(), target%id, target%level, &
        'download_link: source, level =', &
        source%id, source%level, nbor%download, lock_it
      flush (io_unit%log)
    end if
    !---------------------------------------------------------------------------
    ! Skip nbor not marked for download
    !---------------------------------------------------------------------------
    if (nbor%download) then
     if (source%nt == 1) &
       self%same => same_linear
     if (trim(target%kind) /= 'rt_solver') then
       if (source%is_clear(bits%no_download)) then
         !----------------------------------------------------------------------
         ! Default same interpolation when target and source levels agree
         !----------------------------------------------------------------------
         if (source%level==target%level) then
           select case (trim(use_same))
           case ('different')
             call different (self, target, source, jt, pt, all_cells_l, only=only)
           case ('guard_zones')
             call guard_zones%method (target, source, jt, pt)
           case ('guard_same')
             call guard_zones%method (target, source, jt, pt, same=.true.)
           case default
             call self%same (target, source, jt, pt, all_cells_l, only=only)
           end select
         !----------------------------------------------------------------------
         ! Prolong with one level difference and all_cells
         !----------------------------------------------------------------------
         else if (source%level==target%level-1 .and. present(all_cells)) then
           call remesh%prolong (target, source, jt, pt, all_cells)
         !----------------------------------------------------------------------
         ! Prolong cases
         !----------------------------------------------------------------------
         else if ((source%level-target%level) <= -1) then
           select case (trim(use_prolong))
           case ('different')
             call different (self, target, source, jt, pt, all_cells_l, only=only)
           case ('prolong')
             call prolong (self, target, source, jt, pt, all_cells_l, only=only)
           case ('guard_zones')
             call guard_zones%method (target, source, jt, pt)
           case ('remesh')
             call remesh%prolong (target, source, jt, pt, all_cells_l)
           case default
             call different_prolong (self, target, source, jt, pt, all_cells_l, only=only)
           end select
         !----------------------------------------------------------------------
         ! Prolong cases with all_cell
         !----------------------------------------------------------------------
         else if (source%level < target%level .and. present(all_cells)) then
           call different_prolong (self, target, source, jt, pt, all_cells_l, only=only)
         !----------------------------------------------------------------------
         ! Restrict cases
         !----------------------------------------------------------------------
         else if ((source%level-target%level) == 1) then
           select case (trim(use_restrict))
           case ('different')
             call different (self, target, source, jt, pt, all_cells_l, only=only)
           case ('diff_pro')
             call different_prolong (self, target, source, jt, pt, all_cells_l, only=only)
           case ('guard_zones')
             call guard_zones%method (target, source, jt, pt)
           case default
             call different_restrict (self, target, source, jt, pt, only=only)
           end select
         end if
       end if
       if (target%id==id_debug) then
         write(io_unit%log,'(a,i6,1p,4e12.5)') 'source: id,min,max,fmin,fmaxval:', &
                  source%id, &
                  minval(source%mem(:,:,:,source%idx%d,source%it,1)), &
                  maxval(source%mem(:,:,:,source%idx%d,source%it,1)), &
                  source%fminval(source%idx%d), &
                  source%fmaxval(source%idx%d)
         write(io_unit%log,'(a,i6,1p,4e12.5)') 'target: id,min,max,fmin,fmaxval:', &
                  target%id, &
                  minval(target%mem(:,:,:,target%idx%d,target%it,1)), &
                  maxval(target%mem(:,:,:,target%idx%d,target%it,1)), &
                  target%fminval(target%idx%d), &
                  target%fmaxval(target%idx%d)
       end if
     else
       if (trim(source%kind)=='rt_solver') then
         if (abs(source%level-target%level) < 2) then
           if (source%level==target%level) then
             call same_rt (self, target, source, jt, pt)
           else
             call different_rt (self, target, source, jt, pt)
           end if
         end if
       end if
     end if
    end if
    if (jt(2) /= jt(1)) &
      call source%mem_lock(jt(2))%unset ('download_link2')
    call source%mem_lock(jt(1))%unset ('download_link1')
    if (verbose > 1) then
      write(io_unit%log,'(f12.6,i6,i3,3x,a,i6,i3,2x,2l1)') wallclock(), target%id, target%level, &
        'download_link: source, level =', source%id, source%level
      flush (io_unit%log)
    end if
    nbor => nbor%next
  end do
  if (present(only)) then
    if (only==target%idx%d) then
      call target%check_density (" after download")
    end if
  else
    call target%check_density (" after download")
  end if
  if (self%check_filled) then
    if (all_cells_l) then
      l = target%mesh%li
      u = target%mesh%ui
    else
      l = target%mesh%lb
      u = target%mesh%ub
    end if
    n_not = 0
    if (any(filled(l(1):u(1),l(2):u(2),l(3):u(3),:) == -1)) then
      !-------------------------------------------------------------------------
      ! Since only the parent patch is used, this condition is expected, and OK
      !-------------------------------------------------------------------------
      if (verbose > 1) then
        if (all_cells_l) then
          write (io_unit%log,'(a,i6)') &
           'INFO: some internal cells were not filled, iD =', target%id
        else
          write (io_unit%log,'(a,i6)') &
           'INFO: some guard zone cells were not filled, iD =', target%id
        end if
        do iz=l(3),u(3)
        do iy=l(2),u(2)
        do ix=l(1),u(1)
          if (filled(ix,iy,iz,1) == -1) then
            n_not = n_not+1
            write(io_unit%log,'(a,3i3)') 'not filled:', ix, iy, iz
          end if
        end do
        end do
        end do
        if (all_cells_l) then
          write (io%output,'(i6,2x,a,i6,i4)') &
           n_not, 'internal cells were not filled: id, level =', &
           target%id, target%level
        else
          write (io%output,'(i6,2x,a,i6,i4)') &
           n_not, 'guard zone cells were not filled: id, level =', &
           target%id, target%level
        end if
      end if
    end if
    associate (task => link%task)
    select type (task)
    class is (patch_t)
    task%not_filled = n_not
    end select
    end associate
    deallocate (filled)
    if (verbose > 0) &
      deallocate (src)
  end if
  !-----------------------------------------------------------------------------
  call guard_zones%info()
  if (target%is_clear(bits%no_density)) &
    call validate%check (target%link, d, ' after download')
  call trace%end (itimer)
END SUBROUTINE download_link

!===============================================================================
!> Find the target index ranges that correspond to the overlap region, for
!> download from source (source has higher, same, or lower level than target).
!>
!> The first index, l(i), returned is the index of the first target point that
!> is inside the authoritative region; i.e., it must lie to the right of the
!> source point r(li).
!> FIXME: To be validated/improved for moving meshes
!===============================================================================
SUBROUTINE dnload_range (target, source, iv, l, u, offset)
  class(patch_t), pointer :: target, source
  integer                 :: iv, l(3), u(3), offset(3), li, ui
  optional                :: offset
  !.............................................................................
  real(8), dimension(3)   :: h, dist, la, ua
  integer                 :: i
  class(mesh_t), pointer  :: sm, tm
  logical                 :: debug
  !-----------------------------------------------------------------------------
  call trace%begin ('download_t%dnload_range', 2)
  !-----------------------------------------------------------------------------
  ! First take care of mapping the distance between the source and target
  ! centers, periodically, so this can be ignored in the next steps.
  !-----------------------------------------------------------------------------
  dist = source%position - target%position
  if (target%id == id_debug) write (io_unit%log,*) 'dist(1)', dist
  where (source%periodic)
    dist = modulo(dist+0.5_8*source%mesh%b,source%mesh%b) - 0.5_8*source%mesh%b
  end where
  if (target%id == id_debug) write (io_unit%log,*) 'dist(2)', dist
  !-----------------------------------------------------------------------------
  ! Add the distances from the source center to the start and end of the
  ! authoritative domain of each variable.
  !-----------------------------------------------------------------------------
  debug = source%level /= target%level
  la = dist - source%size/2d0
  ua = dist + source%size/2d0
  if (target%id == id_debug) write (io_unit%log,*) 'la,ua(1)', la, ua
  if (source%level /= target%level) then
    do i=1,3
      sm => source%mesh(i)
      tm => target%mesh(i)
      la(i) = la(i) + sm%h(iv)*sm%d - tm%h(iv)*tm%d
      ua(i) = ua(i) + sm%h(iv)*sm%d - tm%h(iv)*tm%d
    end do
  end if
  if (target%id == id_debug) write (io_unit%log,*) 'la,ua(2)', la, ua
  l = ceiling(la/target%mesh%d + target%mesh%o)
  u =   floor(ua/target%mesh%d + target%mesh%o)
  if (target%id == id_debug) write (io_unit%log,*) 'l,u(1)', l, u
  !-----------------------------------------------------------------------------
  ! Make sure the range is legal.  The u index can become smaller than lb, which
  ! means there is no overlap.  Likewise, the l index can become larger than ub
  !-----------------------------------------------------------------------------
  l = max(l,target%mesh%lb)
  u = min(u,target%mesh%ub)
  if (target%id == id_debug) write (io_unit%log,*) 'l,u(2)', l, u
  !-----------------------------------------------------------------------------
  ! Account for collapsed dimensions
  !-----------------------------------------------------------------------------
  where (target%n == 1)
    l = min(l,target%mesh%ub)
    u = max(u,target%mesh%lb)
  end where
  if (target%id == id_debug) write (io_unit%log,*) 'l,u(3)', l, u
  !-----------------------------------------------------------------------------
  ! For same level patch, the integer offset is what should be added to the
  ! target index to get the source index.  For aligned meshes, offset=0.
  !-----------------------------------------------------------------------------
  if (present(offset)) then
    if (source%level==target%level) then
      offset = nint(-dist/target%ds + source%offset - target%offset)
    else
      write(io_unit%log,*) &
        'ERROR: fixed offset not valid for different levels', &
        target%id, target%level, source%id, source%level
      call mpi%abort('dnload_range')
    end if
    if (verbose > 1 .or. io_unit%do_validate .or. target%id==id_debug) &
      write(io_unit%log,1) &
        target%id, target%level, iv, source%id, source%level, l, u, offset
      1 format("download_t%dnload_range: id,lv,iv:",i6,2i3, &
               "  source,lv,l,u,o:",i6,i3,3(2x,3i4))
    if (target%id == id_debug) write (io_unit%log,*) 'l,u,o', l, u, offset
  else
    if (verbose > 1 .or. io_unit%do_validate .or. target%id==id_debug) &
      write(io_unit%log,2) &
        target%id, target%level, iv, source%id, source%level, l, u
      2 format("download_t%dnload_range: id,lv,iv:",i6,2i3, &
               "  source,lv,l,u:",i6,i3,3(2x,3i4))
  end if
  call trace%end()
contains
!===============================================================================
!> Local function that returns the index in target index space, given a position
!> p, in relative target coordinates.
!===============================================================================
  function target_index (p, roundup, nearest) result (ii)
    real(8)           :: p(3)
    logical, optional :: roundup, nearest
    integer           :: ii(3)
    !---------------------------------------------------------------------------
    p = p/target%mesh%d + target%offset                   ! float index
    if (present(roundup)) then
      ii = ceiling(p)                                     ! round up
    else if (present(nearest)) then
      ii = nint(p)                                        ! nearest int
    else
      ii = p + 0.0001                                     ! integer part
    end if
  end function target_index
END SUBROUTINE dnload_range

!===============================================================================
!> Copy guard zone info from source patch at the same level
!>
!> same
!>   dnload_range
!>   time_interval
!===============================================================================
SUBROUTINE same_linear (self, target, source, jt, pt, all_cells, only, all_times)
  class(download_t):: self
  class(patch_t), pointer :: target, source
  logical, optional       :: all_cells, all_times
  integer, optional       :: only
  !.............................................................................
  integer                 :: l(3), u(3), o(3), ii(3), jt(2), i
  real(8)                 :: position(3), xx, dist(3), start
  real                    :: pt(2), fmin
  integer                 :: ix, iy, iz, iv, iv1, iv2, li(3), ui(3), cells
  integer, save           :: itimer=0
  logical                 :: all_cells_l, no
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%same_linear', itimer=itimer)
  if (detailed_timer) then
    cells = 0
    start = wallclock()
  end if
  if (target%id == id_debug) &
    write(io_unit%log,'(a,1p,e14.5,i7,e14.5)') &
      'dbg download_t%same_linear: time, source =', &
      target%time, source%id, source%time
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = target%nv
  end if
  if (present(all_cells)) then
    all_cells_l = all_cells
  else
    all_cells_l = .false.
  end if
  if (verbose > 1) then
    write(io_unit%log,'(a,2i7,2i4)') &
      ' same_linear: target, source, iv1, iv2 =', target%id, source%id, iv1, iv2
  end if
  if (iv2 > source%nv) then
    write(stderr,*) 'iv loop limit:', iv1, iv2
    write(stderr,*) 'source: id, nv =', source%id, source%nv
    write(stderr,*) 'source kind ',source%kind, source%position
    write(stderr,*) 'target kind ',target%kind, target%position
    write(stderr,*) 'target: id, nv =', target%id, target%nv
    call io%abort("download::same_linear::iv2 > source%nv")
  end if
  do iv=iv1, iv2
    if (verbose > 1 .and. source%unsigned(iv)) then
      fmin = min(minval(source%mem(:,:,:,iv,jt(1),1)), &
                 minval(source%mem(:,:,:,iv,jt(2),1)))
      if (fmin <= 0.0) then
        write(io_unit%log,*) 'iv, jt, pt / minval, minloc =', iv, jt, pt
        write(io_unit%log,*) minval(source%mem(:,:,:,iv,jt(1),1)), &
                             minloc(source%mem(:,:,:,iv,jt(1),1))
        write(io_unit%log,*) minval(source%mem(:,:,:,iv,jt(2),1)), &
                             minloc(source%mem(:,:,:,iv,jt(2),1))
        flush(io_unit%log)
      end if
    end if
    if (iv == target%idx%phi) then
      if (verbose > 1) write(io_unit%log,*) &
        target%id, '     using same level potential from', source%id
    end if
    if (source%quality > 1000.) then
      if (iv==target%idx%d .or. (iv >= target%idx%px .and. iv <= target%idx%pz)) then
        all_cells_l = .true.
        write (io_unit%log,*) 'dnload all cells for target, source =', &
          target%id, source%id, source%quality, source%level
      end if
    end if
    ! -- `dphi` is part of `mem` when `time_order` = 0; do not download it.
    if (iv == target%idx%dphi) cycle 
    li = target%mesh%li
    ui = target%mesh%ui
    call dnload_range (target, source, iv, l, u, o)
    if (target%id==id_debug .or. verbose > 1) &
      write(io_unit%log,'(2i7," a l:",3i3,3x,"u:",3i3,3x,"o:",3i3)') &
        target%id,source%id,l,u,o
    l = max(l,source%mesh%lb-o)
    u = min(u,source%mesh%ub-o)
    if (target%id==id_debug .or. verbose > 1) then
      write(io_unit%log,'(2i7," b l:",3i3,3x,"u:",3i3,3x,"o:",3i3)') &
        target%id,source%id,l,u,o
      flush(io_unit%log)
    end if
    if (product(max(0,u-l+1))>0) then
      !-------------------------------------------------------------------------
      ! Extend the pick-up region to fill the guard zones also in the BC regions
      !-------------------------------------------------------------------------
      do i=1,3
        !if (source%mesh(i)%lower_boundary) l(i) = target%mesh(i)%lb
        !if (source%mesh(i)%upper_boundary) u(i) = target%mesh(i)%ub
        !if (i==2) then
          !if (source%boundaries%is_set (bits%yl)) l(i) = target%mesh(i)%lb
          !if (source%boundaries%is_set (bits%yu)) u(i) = target%mesh(i)%ub
        !end if
  !      if (target%mesh(i)%h(iv) /= 0.0) li(i) = li(i)+1
        ! In no-no-mans-land, the `li`-th cell of variables that are staggered in
        ! direction `i` lie just outside the patch edge and should therefore be
        ! interpolated. Currently only implemented for ZEUS solvers; TBD if it
        ! should be applied to STAGGER solvers.
        if (target%kind(1:4) == 'zeus' .and. &
           target%mesh(i)%h(iv) == -0.5 .and. &
           ui(i) > 1) &
            li(i) = li(i)+1
      end do
      ! for staggered quantities, include the ub+1 zone for the appropriate
      ! components note: p* and b* are face-centred in the same direction as the
      ! component, while emfs are face-centred in both directions *orthogonal*
      !  to the component direction.
      if (target%kind(1:4) == 'zeus' .and. &
          iv >= target%idx%px .and. iv <= target%idx%pz) then
        do i=1,3
          if (iv == target%idx%px+i-1 .and. target%n(i) > 1) u(i) = u(i) + 1
        end do
      end if
      if (target%kind(1:4) == 'zeus' .and. &
          iv >= target%idx%bx .and. &
          iv <= target%idx%bz) then
        do i=1,3
          if (iv == target%idx%bx+i-1 .and. target%n(i) > 1) u(i) = u(i) + 1
        end do
      end if
      if (target%kind(1:4) == 'zeus' .and. &
          iv >= target%idx%ex .and. &
          iv <= target%idx%ez) then
        do i=1,3
          if (iv /= target%idx%ex+i-1 .and. target%n(i) > 1) u(i) = u(i) + 1
        end do
      end if
      if (target%id==id_debug) &
        write(io_unit%log,1) source%id, l, u, o, jt, pt
        1 format("download_t%same:", i6,"  l,u,o:",3(3i4,2x)," jt,pt:", &
                 2i4,2x,2f6.2)
      !-------------------------------------------------------------------------
      ! Special case, with nt==1
      !-------------------------------------------------------------------------
      cells = product(u-l+1)
      if (source%nt == 1) then
        if (all_cells_l) then
          do iz=l(3),u(3)
            do iy=l(2),u(2)
              do ix=l(1),u(1)
                target%mem(ix,iy,iz,iv,target%it,1) =  &
                  source%mem(ix+o(1),iy+o(2),iz+o(3),iv,1,1)
              end do
              if (self%check_filled) then
                do ix=l(1),u(1)
                  filled(ix,iy,iz,iv) = source%level
                end do
                if (verbose > 0) then
                  do ix=l(1),u(1)
                    src(ix,iy,iz) = source%id
                  end do
                end if
              end if
            end do
          end do
        else
          do iz=l(3),u(3)
            do iy=l(2),u(2)
              do ix=l(1),u(1)
                ii = [ix,iy,iz]
                if (ix+o(1) > size(source%mem,1)) then
                  print '(5(3i4,2x))', l, u, o, shape(source%mem)
                end if
                target%mem(ix,iy,iz,iv,target%it,1) = &
                  merge(target%mem(ix,iy,iz,iv,target%it,1),  &
                        source%mem(ix+o(1),iy+o(2),iz+o(3),iv,1,1), &
                        all((ii>=li).and.(ii<=ui)))
              end do
              if (self%check_filled) then
                do ix=l(1),u(1)
                  ii = [ix,iy,iz]
                  no = all((ii>=li).and.(ii<=ui))
                  filled(ix,iy,iz,iv) = merge(filled(ix,iy,iz,iv), &
                    int(source%level,kind=int8), no)
                end do
                if (verbose > 0) then
                  do ix=l(1),u(1)
                    ii = [ix,iy,iz]
                    no = all((ii>=li).and.(ii<=ui))
                    src(ix,iy,iz) = merge(src(ix,iy,iz), source%id, no)
                  end do
                end if
              end if
            end do
          end do
        end if      
      !-------------------------------------------------------------------------
      ! Unsigned values; interpolated in the log
      !-------------------------------------------------------------------------
      else if (source%unsigned(iv)) then
        !-----------------------------------------------------------------------
        ! Fill all cells, even interior ones, when all_cells==.true., but only
        ! for cells that are not already filled.  The filling is expensive in
        ! this case, so use an if test.
        !-----------------------------------------------------------------------
        if (all_cells_l) then
          do iz=l(3),u(3)
            do iy=l(2),u(2)
              if (self%check_filled) then
                do ix=l(1),u(1)
                  if (filled(ix,iy,iz,iv) < source%level) then
                    target%mem(ix,iy,iz,iv,target%it,1) = exp( &
                      pt(1)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)) + &
                      pt(2)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)))
                    filled(ix,iy,iz,iv) = source%level
                    if (verbose > 0) &
                      src(ix,iy,iz) = source%id
                  else
                    write (io_unit%log,*)'WARNING: prevented overwrite', target%id, &
                      source%id, ix,iy,iz,o,iv, filled(ix,iy,iz,iv), src(ix,iy,iz)
                  end if
                end do
              else
                do ix=l(1),u(1)
                  target%mem(ix,iy,iz,iv,target%it,1) = exp( &
                    pt(1)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)) + &
                    pt(2)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)))
                end do
              end if
            end do
          end do
        !-----------------------------------------------------------------------
        ! Make sure to change only guard cells, when all_cells==.false..  Here
        ! the calculation is cheap, so we use a merge to allow vectorization.
        !-----------------------------------------------------------------------
        else
          do iz=l(3),u(3)
            do iy=l(2),u(2)
              do ix=l(1),u(1)
                ii = [ix,iy,iz]
!                if (source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) <= 0.0) then
!                  write (io_unit%log,9) wallclock(), &
!                    l,u,ix,iy,iz,o,iv,source%it,source%new,1,jt, &
!                    source%id,source%time,target%id,target%time,source%iit
!                9 format(f12.6,6(3i4,2x),2(i6,f12.6),7i3)
!                end if
!                if (source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1) <= 0.0) then
!                  write (io_unit%log,9) wallclock(), &
!                    l,u,ix,iy,iz,o,iv,source%it,source%new,2,jt, &
!                    source%id,source%time,target%id,target%time,source%iit
!                end if
                target%mem(ix,iy,iz,iv,target%it,1) = &
                  merge(target%mem(ix,iy,iz,iv,target%it,1), &
                    exp(pt(1)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)) + &
                        pt(2)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1))), &
                    all((ii>=li).and.(ii<=ui)))
              end do
              if (self%check_filled) then
                do ix=l(1),u(1)
                  ii = [ix,iy,iz]
                  no = all((ii>=li).and.(ii<=ui))
                  filled(ix,iy,iz,iv) = merge(filled(ix,iy,iz,iv), &
                    int(source%level,kind=int8), no)
                end do
                if (verbose > 0) then
                  do ix=l(1),u(1)
                    ii = [ix,iy,iz]
                    no = all((ii>=li).and.(ii<=ui))
                    src(ix,iy,iz) = merge(src(ix,iy,iz), source%id, no)
                  end do
                end if
              end if
            end do
          end do
        end if
      !-------------------------------------------------------------------------
      ! Signed values
      !-------------------------------------------------------------------------
      else
        !-----------------------------------------------------------------------
        ! All cells
        !-----------------------------------------------------------------------
        if (all_cells_l) then
          do iz=l(3),u(3)
            do iy=l(2),u(2)
              !-----------------------------------------------------------------
              ! Prevent overwriting better values
              !-----------------------------------------------------------------
              if (self%check_filled) then
                do ix=l(1),u(1)
                  no = filled(ix,iy,iz,iv) > source%level
                  target%mem(ix,iy,iz,iv,target%it,1) =  &
                    merge (target%mem(ix,iy,iz,iv,target%it,1), &
                      pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
                      pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1), &
                      no)
                  filled(ix,iy,iz,iv) = merge (filled(ix,iy,iz,iv), &
                    int(source%level, kind=int8), no)
                end do
                if (verbose > 0) then
                  do ix=l(1),u(1)
                    no = filled(ix,iy,iz,iv) > source%level
                    src(ix,iy,iz) = merge(src(ix,iy,iz), source%id, no)
                  end do
                end if
              !-----------------------------------------------------------------
              ! No questions asked
              !-----------------------------------------------------------------
              else
                do ix=l(1),u(1)
                  target%mem(ix,iy,iz,iv,target%it,1) =  &
                    pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
                    pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)
                end do
              end if
            end do
          end do
        !-----------------------------------------------------------------------
        ! Only guard cells.  Prefer same level, so don't check for levels
        !-----------------------------------------------------------------------
        else
          do iz=l(3),u(3)
            do iy=l(2),u(2)
              do ix=l(1),u(1)
                ii = [ix,iy,iz]
                target%mem(ix,iy,iz,iv,target%it,1) = &
                  merge(target%mem(ix,iy,iz,iv,target%it,1), &
                    pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
                    pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1), &
                    all((ii>=li).and.(ii<=ui)))
              end do
              if (self%check_filled) then
                do ix=l(1),u(1)
                  ii = [ix,iy,iz]
                  no = all((ii>=li).and.(ii<=ui))
                  filled(ix,iy,iz,iv) = &
                    merge (filled(ix,iy,iz,iv), int(source%level,kind=int8), no)
                end do
                if (verbose > 0) then
                  do ix=l(1),u(1)
                    ii = [ix,iy,iz]
                    no = all((ii>=li).and.(ii<=ui))
                    src(ix,iy,iz) = merge(src(ix,iy,iz), source%id, no)
                  end do
                end if
              end if
            end do
          end do
        end if
      end if
      if (verbose > 0 .or. target%id==id_debug) then
        if (iv == 1 .or. verbose > 2) &
          write (io_unit%log,2) target%id, target%level, source%id, source%level, &
          'same_line picked up', &
          product(u-l+1), ' values, l, u, mind, pt =', l, u, &
          minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,target%it,1)), pt
        2 format(2(i6,i3),2x,a,i6,a,2(2x,3i4),1p,e14.4,2x,0p2f6.3)
      end if
    end if
  end do
  if (detailed_timer) then
    if (cells > 0) then
      write (stdout,'(a,6p,2(2x,i6,i3),i6," cells,",f9.3," mus/cell")') &
        "download_t%dsame_linear: target, level, source, level,", &
        target%id, target%level, source%id, source%level, cells, &
        (wallclock()-start)/max(1,cells)
    else
      write (stdout,'(a,6p,2(2x,i6,i3),i6," cells,",f9.3," mus")') &
        "download_t%same_linear: target, level, source, level,", &
        target%id, target%level, source%id, source%level, cells, wallclock()-start
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE same_linear

!===============================================================================
!> Copy guard zone info from source patch at the same level
!>
!> same
!>   dnload_range
!>   time_interval
!===============================================================================
SUBROUTINE same_lagrange (self, target, source, jj, pt, all_cells, only, all_times)
  class(download_t):: self
  class(patch_t), pointer :: target, source
  logical, optional       :: all_cells, all_times
  integer, optional       :: only
  integer                 :: jj(2)
  real                    :: pt(2)
  !.............................................................................
  integer                 :: l(3), u(3), o(3), ii(3), jt, i
  real(8)                 :: position(3), xx, dist(3)
  integer                 :: ix, iy, iz, iv, iv1, iv2, li(3), ui(3)
  integer, save           :: itimer=0, nprint=20
  logical                 :: all_cells_l, filter
  integer, parameter      :: order=2
  real(8)                 :: w(order+1)
  integer                 :: j0, j1, iit(source%nt-1)
  real(8)                 :: times(source%nt-1)
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%same_lagrange', itimer=itimer)
  if (target%id == id_debug) &
    write(io_unit%log,*) 'dbg download_t%same_lagrange: time =', target%time
  !-----------------------------------------------------------------------------
  ! Compute weights for Lagrange interpolation
  !-----------------------------------------------------------------------------
  call source%timeslots (iit, times)
  call lagrange%sequence_weights (target%time, times, j0, j1, w, order)
  if (verbose > 2) &
    write(io_unit%log,'(a,2i4,1p,10e14.6)') &
      'same_lagrange: j0, j1, time, times, w:', &
    j0, j1, target%time, times(j0:j1), w
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = target%nv
  end if
  if (present(all_cells)) then
    all_cells_l = all_cells
  else
    all_cells_l = .false.
  end if
  do iv=iv1, iv2
    if (iv == target%idx%phi) then
      if (verbose>1) write(io_unit%log,*) &
        target%id, '     using same level potential from', source%id
    end if
    if (source%quality > 1000.) then
      if (iv==target%idx%d .or. &
         (iv >= target%idx%px .and. iv <= target%idx%pz)) then
          all_cells_l = .true.
          write (io_unit%log,*) '(2)dnload all cells for target, source =', &
            target%id, source%id, source%quality, source%level
      end if
    end if
    ! -- `dphi` is part of `mem` when `time_order` = 0; do not download it.
    if (iv == target%idx%dphi) cycle 
    li = target%mesh%li
    ui = target%mesh%ui
    call dnload_range (target, source, iv, l, u, o)
    if (target%id==id_debug) &
      write(io_unit%log,'(2i7," a l:",3i3,3x,"u:",3i3,3x,"o:",3i3)') &
        target%id,source%id,l,u,o
    l = max(l,source%mesh%lb-o)
    u = min(u,source%mesh%ub-o)
    if (target%id==id_debug) &
      write(io_unit%log,'(2i7," b l:",3i3,3x,"u:",3i3,3x,"o:",3i3)') &
        target%id,source%id,l,u,o
    if (product(max(0,u-l+1))>0) then
      !-------------------------------------------------------------------------
      ! Extend the pick-up region to fill the guard zones also in the BC regions
      !-------------------------------------------------------------------------
      do i=1,3
        if (source%mesh(i)%lower_boundary) l(i) = target%mesh(i)%lb
        if (source%mesh(i)%upper_boundary) u(i) = target%mesh(i)%ub
      end do
      !-------------------------------------------------------------------------
      ! Check if we need to filter away internal points
      !-------------------------------------------------------------------------
      filter = .false.
      if (.not. all_cells_l) then
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              ii = [ix,iy,iz]
              if (all((ii>=li).and.(ii<=ui))) then
                filter = .true.
                exit
              end if
            end do
            if (filter) exit
          end do
          if (filter) exit
        end do
      end if
      if (filter) then
        if (nprint > 0) then
          write(io_unit%log,*) &
'WARNING: attempt to set internal cells in download_t%same_lagrange, using filter!'
          nprint = nprint - 1
          if (nprint==0) then
            write(io_unit%log,*) &
              'target id, position =', target%id, target%position
            write(io_unit%log,*) &
              'source id, position =', source%id, source%position
            write(io_unit%log,*) 'l , u =', l, u
            write(io_unit%log,*) 'gn, o =', target%gn, o
          end if
        end if
      end if
      !---------------------------------------------------------------------------
      ! Unsigned => use interpolation in the log
      !---------------------------------------------------------------------------
      if (source%unsigned(iv)) then
        do iz=l(3),u(3)
        do iy=l(2),u(2)
          if (filter) then
            do ix=l(1),u(1)
              ii = [ix,iy,iz]
              if (all_cells_l .or. any((ii<li).or.(ii>ui))) then
                do jt=j0,j1
                  if (jt==j0) then
                    target%mem(ix,iy,iz,iv,target%it,1) = &
                      w(1+jt-j0)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,iit(jt),1))
                  else
                    target%mem(ix,iy,iz,iv,target%it,1) = &
                    target%mem(ix,iy,iz,iv,target%it,1) + &
                      w(1+jt-j0)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,iit(jt),1))
                  end if
                end do
                target%mem(ix,iy,iz,iv,target%it,1) = &
                  exp(target%mem(ix,iy,iz,iv,target%it,1))
                if (self%check_filled) then
                  filled(ix,iy,iz,iv) = source%level
                  if (verbose > 0) &
                    src(ix,iy,iz) = source%id
                end if
              end if
            end do
          else
            do jt=j0,j1
              if (jt==j0) then
                do ix=l(1),u(1)
                  target%mem(ix,iy,iz,iv,target%it,1) = &
                    w(1+jt-j0)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,iit(jt),1))
                end do
              else
                do ix=l(1),u(1)
                  target%mem(ix,iy,iz,iv,target%it,1) = &
                  target%mem(ix,iy,iz,iv,target%it,1) + &
                    w(1+jt-j0)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,iit(jt),1))
                end do
              end if
              if (self%check_filled) then
                do ix=l(1),u(1)
                  target%mem(ix,iy,iz,iv,target%it,1) = &
                    exp(target%mem(ix,iy,iz,iv,target%it,1))
                  filled(ix,iy,iz,iv) = source%level
                end do
                if (verbose > 0) then
                  do ix=l(1),u(1)
                    src(ix,iy,iz) = source%id
                  end do
                end if
              else
                do ix=l(1),u(1)
                  target%mem(ix,iy,iz,iv,target%it,1) = &
                    exp(target%mem(ix,iy,iz,iv,target%it,1))
                end do
              end if
            end do
          end if
        end do
        end do
      !---------------------------------------------------------------------------
      ! Use linear interpolation
      !---------------------------------------------------------------------------
      else
        do iz=l(3),u(3)
        do iy=l(2),u(2)
          if (filter) then
            do ix=l(1),u(1)
              ii = [ix,iy,iz]
              if (all_cells_l .or. any((ii<li).or.(ii>ui))) then
                do jt=j0,j1
                  if (jt==j0) then
                    target%mem(ix,iy,iz,iv,target%it,1) = &
                      w(1+jt-j0)*(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,iit(jt),1))
                  else
                    target%mem(ix,iy,iz,iv,target%it,1) = &
                    target%mem(ix,iy,iz,iv,target%it,1) + &
                      w(1+jt-j0)*(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,iit(jt),1))
                  end if
                end do
                if (self%check_filled) then
                  filled(ix,iy,iz,iv) = source%level
                  if (verbose > 0) &
                    src(ix,iy,iz) = source%id
                end if
              end if
            end do
          else
            do jt=j0,j1
              if (jt==j0) then
                do ix=l(1),u(1)
                  target%mem(ix,iy,iz,iv,target%it,1) = &
                    w(1+jt-j0)*(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,iit(jt),1))
                end do
              else
                do ix=l(1),u(1)
                  target%mem(ix,iy,iz,iv,target%it,1) = &
                  target%mem(ix,iy,iz,iv,target%it,1) + &
                    w(1+jt-j0)*(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,iit(jt),1))
                end do
              end if
              if (self%check_filled) then
                do ix=l(1),u(1)
                  filled(ix,iy,iz,iv) = source%level
                end do
                if (verbose > 0) then
                  do ix=l(1),u(1)
                    src(ix,iy,iz) = source%id
                  end do
                end if
              end if
            end do
          end if
        end do
        end do
      end if
      if ((verbose>2.or.target%id==id_debug).and.product(max(0,u-l+1))>0) &
        write(io_unit%log,'(i6,4(2x,a,i6),1p,2e12.4)') target%id, &
        'same picked up', product(max(0,u-l+1)), &
        'values of', product(target%n), &
        'from', source%id, &
        ' iv, min, max =', iv, &
        minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,target%it,1)), &
        maxval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,target%it,1))
    end if
  end do
  call trace%end (itimer)
END SUBROUTINE same_lagrange

!===============================================================================
!> Same as above, but for RT_solver type patches
!===============================================================================
SUBROUTINE same_rt (self, target, source, jt, pt)
  class(download_t):: self
  class(patch_t), pointer :: target, source
  !integer(int8), dimension(:,:,:,:), allocatable:: filled
  !.............................................................................
  integer                 :: l(3), u(3), o(3), ii(3), jt(2), li(3), ui(3)
  real(8)                 :: position(3)
  real                    :: pt(2)
  integer                 :: ix, iy, iz, iv
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%same_rt', itimer=itimer)
  call dnload_range (target, source, 1, l, u, o)
  li = target%mesh%li
  ui = target%mesh%ui
  do iv=1,target%nv
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii = [ix,iy,iz]
          if (self%check_filled) then
            if (filled(ix,iy,iz,iv) > source%level) cycle
          end if
          if (all((ii>=li).and.(ii<=ui))) cycle
          target%mem(ix,iy,iz,iv,target%it,1) = &
             pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
             pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)
          if (self%check_filled) then
            filled(ix,iy,iz,iv) = source%level
            if (verbose > 0) &
              src(ix,iy,iz) = source%id
          end if
        end do
      end do
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE same_rt

!===============================================================================
!> Interpolate guard zone info from source patch at a different level
!===============================================================================
SUBROUTINE different (self, target, source, jt, pt, all_cells, only)
  class(download_t):: self
  class(patch_t), pointer :: target, source
  logical, optional:: all_cells
  integer, optional:: only
  !.............................................................................
  integer                 :: ix, iy, iz, iv, jt(2), iv1, iv2, i, cells
  integer, dimension(3)   :: l, u, ii, li, ui, jj
  real(8), dimension(3)   :: pos
  logical                 :: valid
  real                    :: pt(2), pp(3)
  real(kind=KindScalarVar):: qint
  logical                 :: all_cells_l
  real(8)                 :: start
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different', itimer=itimer)
  if (detailed_timer) then
    cells = 0
    start = wallclock()
  end if
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = target%nv
  end if
  if (present(all_cells)) then
    all_cells_l = all_cells
  else
    all_cells_l = .false.
  end if
  do iv=iv1,iv2
    call dnload_range (target, source, iv, l, u)
    li = target%mesh%li
    ui = target%mesh%ui
    all_cells_l = all_cells_l .or. all((l >= li) .and. (u<=ui))
    !---------------------------------------------------------------------------
    ! For no-no-mans-land meshes, use guard zone values for the li index, to
    ! make fluxes up/down symmetric
    !---------------------------------------------------------------------------
!    do i=1,3
!      if (target%mesh(i)%h(iv) /= 0.0) li(i) = li(i)+1
!    end do
    if (verbose > 1 .or. target%id==id_debug) &
      write (io_unit%log,1) &
        source%id, (source%position-target%position)/target%size, l, u, all_cells_l
      1 format("DBG download_t%different:", i6,"  displ:",3f11.6,"  l,u:",2(3i4,2x),l3)
    associate (q1 => source%mem(:,:,:,:,jt(1),1), q2 => source%mem(:,:,:,:,jt(2),1))
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii=[ix,iy,iz]
          if (all_cells_l .or. any((ii-li) < 0 .or. (ii-ui) > 0)) then
            if (iv==1) cells = cells+1
            if (verbose > 1 .and. target%id==id_debug) &
              write (io_unit%log,*) &
                target%id, 'DBG interpolate_t%different:    using', &
                iv,ii, present(all_cells), (ii-li) < 0
            pos = target%myposition (ii, iv)
            call source%index_space (pos, iv, jj, pp)
            if (source%pervolume(iv)) then
              call source%interpolate_specific (target, ii, iv, jt, pt)
            else
              call source%interpolate          (target, ii, iv, jt, pt)
            end if
            if (.not.all_cells_l) then
              if (all((ii>=li .and. ii<=ui))) then
                write (io_unit%log,*)'INTERNAL BREACH!', ii
              end if
            end if
          else
            if (verbose > 1 .or. target%id==id_debug) &
              write (io_unit%log,'(i6,i3,2x,a,i6,i2,i3,2x,3i4,L2,3L1)') &
                target%id, target%level, 'DBG interpolate_t%different: skipping', &
                source%id, source%level, iv,ii, present(all_cells), (ii-li) < 0
            cycle
          end if
        end do
      end do
    end do
    end associate
    if (verbose > 0 .or. target%id==id_debug) then
      if (iv == 1 .or. verbose > 2) &
        write (io_unit%log,2) target%id, target%level, source%id, source%level, &
        'different picked up', &
        product(u-l+1), ' values, l, u, mind, pt =', l, u, &
        minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,target%it,1)), pt
      2 format(2(i6,i3),2x,a,i6,a,2(2x,3i4),1p,e14.4,2x,0p2f6.3)
    end if
  end do
  if (detailed_timer) then
    if (cells > 0) then
      write (stdout,'(a,6p,2(2x,i6,i3),i6," cells,",f9.3," mus/cell")') &
        "download_t%different: target, level, source, level,", &
        target%id, target%level, source%id, source%level, cells, &
        (wallclock()-start)/max(1,cells)
    else
      write (stdout,'(a,6p,2(2x,i6,i3),i6," cells,",f9.3," mus")') &
        "download_t%different: target, level, source, level,", &
        target%id, target%level, source%id, source%level, cells, wallclock()-start
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE different

!===============================================================================
!> Interpolate guard zone info from source patch at a lower level
!===============================================================================
SUBROUTINE different_prolong (self, target, source, jt, pt, all_cells, only)
  class(download_t):: self
  class(patch_t), pointer :: target, source
  logical, optional:: all_cells
  integer, optional:: only
  !.............................................................................
  integer                 :: ntot, nskip, cells, nint, nperv, nexp
  integer                 :: ix, iy, iz, iv, jt(2), iv1, iv2, i, n, srcid, m(3)
  integer, dimension(3)   :: l, u, ii, li, ui, jj
  real(8), dimension(3)   :: pos
  real(4), dimension(3)   :: ot, os, ratio
  real(8)                 :: start
  logical                 :: valid
  real                    :: pt(2), pp(3)
  real(kind=KindScalarVar):: qint, dint
  logical                 :: all_cells_l
  real(kind=KindScalarVar), pointer:: q1(:,:,:), q2(:,:,:), d1(:,:,:), d2(:,:,:)
  real(kind=KindScalarVar), pointer:: q(:,:,:,:)
  real(kind=KindScalarVar), allocatable:: dtmp(:,:,:)
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different_prolong', itimer=itimer)
  start = wallclock()
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = target%nv
  end if
  if (present(all_cells)) then
    all_cells_l = all_cells
  else
    all_cells_l = .false.
  end if
  if (verbose > 1 .or. target%id == id_debug) &
       write(io_unit%log,'(2(a,i7,i4,1p,3e11.3,4x),l3)') &
      'DBG different_prolong; target:', target%id, target%level, target%ds, &
                             'source:', source%id, source%level, source%ds, all_cells_l
  li = target%mesh%li
  ui = target%mesh%ui
  ot = target%mesh%o
  os = target%position - source%position
  os = modulo(os + 0.5_8*target%mesh%b, target%mesh%b) - 0.5_8*target%mesh%b
  os = source%mesh%o + os/source%mesh%d
  ratio = 2.0
  ntot = 0
  nint = 0
  nexp = 0
  nskip = 0
  nperv = 0
  cells = 0
  allocate (dtmp(target%gn(1),target%gn(2),target%gn(3)))
  d1 => source%mem(:,:,:,source%idx%d,jt(1),1)
  d2 => source%mem(:,:,:,source%idx%d,jt(2),1)
  do iv=iv1,iv2
    ! -- `dphi` is part of `mem` when `time_order` = 0; do not download it.
    if (iv == target%idx%dphi) cycle
    call dnload_range (target, source, iv, l, u)
    !---------------------------------------------------------------------------
    ! For no-no-mans-land meshes, use guard zone values for the li index, to
    ! make fluxes up/down symmetric
    !---------------------------------------------------------------------------
!    do i=1,3
!      if (target%mesh(i)%h(iv) /= 0.0) li(i) = li(i)+1
!    end do
    q1 => source%mem(:,:,:,iv,jt(1),1)
    q2 => source%mem(:,:,:,iv,jt(2),1)
    q  => source%mem(:,:,:,iv,:,1)
    all_cells_l = all_cells_l .or. all((l >= li) .and. (u <= ui))
    if (verbose >  0 .or. target%id==id_debug) &
      write (io_unit%log,1) &
        source%id, source%distance(target)/(0.5*(target%size+source%size)), &
        l, u, all_cells_l
      1 format("download_t%different_prolong: source =", i6,"  displ:",3f11.6, &
        "  l,u:",2(3i4,2x)," all_cells:",l2)
    n = 0
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ntot = ntot+1
          ii=[ix,iy,iz]
          if (all_cells_l .or. any((ii-li) < 0 .or. (ii-ui) > 0)) then
            if (self%check_filled) then
              !-----------------------------------------------------------------
              ! If either the cell is already filled with higher level data,
              ! or filled with target level data (which is preferred), cycle.
              !-----------------------------------------------------------------
              if (filled(ix,iy,iz,iv) >  source%level .or. &
                  filled(ix,iy,iz,iv) == target%level) then
                if (allocated(src)) then
                  srcid = src(ix,iy,iz)
                else
                  srcid = -1
                end if
                if (verbose > 2 .or. verbose > 1 .and. target%id == id_debug) &
                  write (io_unit%log,'(a,2(i6,i3),i4,2x,3i4,2x,i4,i6)') &
                    'different_prolong: skipping', target%id, target%level, &
                    source%id, source%level, iv,ii, filled(ix,iy,iz,iv), srcid
                nskip = nskip+1
                cycle
              end if
            end if
            if (iv==1) cells = cells+1
            if (verbose > 1 .and. target%id == id_debug) &
              write (io_unit%log,*) target%id, 'different_prolong:    using', &
                iv,ii, present(all_cells), (ii-li) < 0
            pos = target%myposition (ii, iv)
            call source%index_space (pos, iv, jj, pp)
            !call source%index_space_of (target, ii, iv, jj, pp)
            if (source%unsigned(iv)) then
              nexp = nexp+1
              qint = exp(interpolator%four_d_log (q, jj(1), jj(2), jj(3), jt, [pp,pt(2)]))
              if (iv==target%idx%d) &
                dtmp(ix,iy,iz) = qint
            else if (source%pervolume(iv)) then
              nperv = nperv+1
              qint = interpolator%four_d_pv (q, d1, d2, jj(1), jj(2), jj(3), jt, [pp,pt(2)])
              qint = qint * dtmp(ix,iy,iz)
            else
              qint = interpolator%four_d (q, jj(1), jj(2), jj(3), jt, [pp,pt(2)])
            end if
            n = n + 1
            !if (all((ii >= li .and. ii <= ui)) .and. .not. all_cells_l) then
            !  write (io_unit%log,*) 'INTERNAL BREACH!', ii
            !end if
            target%mem(ix,iy,iz,iv,target%it,1) = qint
            if (self%check_filled) filled(ix,iy,iz,iv) = source%level
          else
            nint = nint+1
            if (verbose > 1 .or. target%id==id_debug) &
              write (io_unit%log,'(i6,i3,2x,a,i6,i2,i3,2x,3i4,L2,2(1x,3L1))') &
                target%id, target%level, 'different_prolong: skipping internal', &
                source%id, source%level, iv, ii, all_cells_l, (ii-li) < 0, (ii-ui) > 0
            cycle
          end if
        end do
      end do
    end do
    if (verbose > 0 .or. target%id==id_debug) then
      if (iv == 1 .or. verbose > 2) &
        write (io_unit%log,2) target%id, target%level, source%id, source%level, &
        'diff_prol picked up', &
        cells, nskip, nint, ' values, l, u, mind, pt =', l, u, &
        minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,target%it,1)), pt
      2 format(2(i6,i3),2x,a,3i6,a,2(2x,3i4),1p,e14.4,2x,0p2f6.3)
    end if
  end do
  deallocate (dtmp)
  if (detailed_timer) then
    if (cells > 0) then
      write (stdout,'(a,6p,2(2x,i6,i3),",",i6," cells,",f9.3," mus/cell")') &
        "download_t%different_prolong: target, level, source, level =", &
        target%id, target%level, source%id, source%level, cells, &
        (wallclock()-start)/cells
    else
      write (stdout,'(a,2(2x,i6,i3),2x,3f8.3,",",i6," cells,",6p,f9.3," mus")') &
        "download_t%different_prolong: target,lvl,source,lvl,rel=", &
        target%id, target%level, source%id, source%level, &
        source%distance(target)/(0.5d0*(target%size+source%size)), &
        cells, wallclock()-start
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE different_prolong

!> =============================================================================
!> Interpolate guard zone info from source patch at a lower level
!> =============================================================================
SUBROUTINE prolong (self, target, source, jt, pt, all_cells, only)
  class(download_t):: self
  class(patch_t), pointer :: target, source
  logical, optional:: all_cells
  integer, optional:: only
  !.............................................................................
  integer                 :: ntot, nskip, cells, nint, nperv, nexp
  integer                 :: ix, iy, iz, iv, jt(2), iv1, iv2, i, m, n, srcid
  integer, dimension(3)   :: l, u, ii, li, ui, jj, jjmin, jjmax, iimin, iimax
  real(8), dimension(3)   :: pos
  real(4), dimension(3)   :: ot, os, ratio
  real(8)                 :: start
  logical                 :: valid
  integer, save           :: itimer=0
  real                    :: pt(2), pp(3)
  real(kind=KindScalarVar):: qint, dint
  logical                 :: all_cells_l
  real(kind=KindScalarVar), allocatable::  sbuf(:,:,:,:,:),  tbuf(:), pbuf(:,:)
  real(kind=KindScalarVar), allocatable:: sdbuf(:,:,:,:,:), tdbuf(:)
  integer                 :: mv(target%nv)
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%prolong', itimer=itimer)
  start = wallclock()
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = target%nv
  end if
  if (present(all_cells)) then
    all_cells_l = all_cells
  else
    all_cells_l = .false.
  end if
  if (verbose > 1 .or. target%id == id_debug) &
    write(io_unit%log,'(2(a,i7,i4,1p,3e11.3,4x),l3)') &
      'DBG prolong; target:', target%id, target%level, target%ds, &
                   'source:', source%id, source%level, source%ds, all_cells_l
  li = target%mesh%li
  ui = target%mesh%ui
  call dnload_range (target, source, target%idx%d, l, u)
  n = product(u-l+2)
  allocate ( sbuf(n,2,2,2,2))
  allocate (sdbuf(n,2,2,2,2))
  allocate ( pbuf(n,4      ))
  allocate ( tbuf(n        ))
  allocate (tdbuf(n        ))
  ot = target%mesh%o
  os = source%mesh%o
  ratio = 2.0
  ntot = 0
  nint = 0
  nexp = 0
  nskip = 0
  nperv = 0
  cells = 0
  do iv=iv1,iv2
    ! `dphi` is part of `mem` when `time_order` = 0; do not download it.
    if (iv == target%idx%dphi) cycle 
    call dnload_range (target, source, iv, l, u)
    !---------------------------------------------------------------------------
    ! For no-no-mans-land meshes, use guard zone values for the li index, to
    ! make fluxes up/down symmetric
    !---------------------------------------------------------------------------
!    do i=1,3
!      if (target%mesh(i)%h(iv) /= 0.0) li(i) = li(i)+1
!    end do
    if (verbose > 1 .or. target%id==id_debug) &
      write (io_unit%log,1) &
        source%id, (source%position-target%position)/target%size, l, u
      1 format("download_t%different_prolong:", i6,"  displ:",3f11.6,"  l,u:", &
        2(3i4,2x))
    all_cells_l = all_cells_l .or. all((l >= li) .and. (u <= ui))
    n = 0
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ntot = ntot+1
          ii=[ix,iy,iz]
          if (all_cells_l .or. any((ii-li) < 0 .or. (ii-ui) > 0)) then
            !if (self%check_filled) then
            !  !-----------------------------------------------------------------
            !  ! If either the cell is already filled with higher level data,
            !  ! or filled with target level data (which is preferred), cycle.
            !  !-----------------------------------------------------------------
            !  if (filled(ix,iy,iz,iv) >  source%level .or. &
            !      filled(ix,iy,iz,iv) == target%level) then
            !    if (allocated(src)) then
            !      srcid = src(ix,iy,iz)
            !    else
            !      srcid = -1
            !    end if
            !    nskip = nskip+1
            !    cycle
            !  end if
            !end if
            if (iv==1) cells = cells+1
            pos = target%myposition (ii, iv)
            call source%index_space (pos, iv, jj, pp)
            !pp = (ii-ot)*ratio + os
            !jj = pp
            !pp = pp - jj
            n = n + 1
            !if (n==1) then
            !  iimin = ii
            !  iimax = ii
            !  jjmin = jj
            !  jjmax = jj
            !else
            !  iimin = min(ii,iimin)
            !  iimax = max(ii,iimax)
            !  jjmin = min(jj,jjmin)
            !  jjmax = max(jj,jjmax)
            !end if
            pbuf(n,1:3) = pp
            pbuf(n,  4) = pt(1)
            sbuf(n,:,:,:,:) = &
              source%mem(jj(1):jj(1)+1,jj(2):jj(2)+1,jj(3):jj(3)+1,iv,jt(:),1)
            !if (self%check_filled) filled(ix,iy,iz,iv) = source%level
            !if (iv==1) &
              !write (io_unit%output,'(a,i4,2(3i4,2x),2x,4f7.3)') &
              !'prol', n, ii, jj, pbuf(n,:)
          else
            nint = nint+1
            !if (verbose > 1 .and. target%id==id_debug) &
              !write (io_unit%output,'(i6,i3,2x,a,i6,i2,i3,2x,3i4,L2,2(1x,3L1))') &
              !  target%id, target%level, 'prolong: skipping internal', &
              !  source%id, source%level, iv, ii, all_cells_l, (ii-li) < 0, 
              !  (ii-ui) > 0
            cycle
          end if
        end do
      end do
    end do
    m = n
    mv(iv) = m
    call mem_interp
    n = 0
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ntot = ntot+1
          ii=[ix,iy,iz]
          if (all_cells_l .or. any((ii-li) < 0 .or. (ii-ui) > 0)) then
            n = n + 1
            target%mem(ix,iy,iz,iv,target%it,1) = tbuf(n)
          end if
        end do
      end do
    end do
  end do
  !print '(a,2(2x,3i4))', 'prolong: iimin, iimax =', iimin, iimax
  !print '(a,2(2x,3i4))', 'prolong: jjmin, jjmax =', jjmin, jjmax
  deallocate (sbuf, pbuf, tbuf, sdbuf, tdbuf)
  if (detailed_timer) then
    if (cells > 0) then
      write (stdout,'(a,6p,2(2x,i6,i3),i6," cells,",f9.3," mus/cell")') &
        "download_t%prolong: target, level, source, level,", &
        target%id, target%level, source%id, source%level, cells, &
        (wallclock()-start)/max(1,cells)
    else
      write (stdout,'(a,6p,2(2x,i6,i3),i6," cells,",f9.3," mus")') &
        "download_t%prolong: target, level, source, level,", &
        target%id, target%level, source%id, source%level, cells, wallclock()-start
    end if
  end if
  call trace%end (itimer)
contains
  subroutine mem_interp 
    real:: q1, q2, q3, q4
    integer:: j
    !print *, iv, m, target%pervolume(iv), target%unsigned(iv)
    !---------------------------------------------------------------------------
    ! Save density corner values
    !---------------------------------------------------------------------------
    if (iv == 1) then
      sdbuf = sbuf
    end if
    !---------------------------------------------------------------------------
    ! Use density corner values for per-volume variables
    !---------------------------------------------------------------------------
    if (target%pervolume(iv)) then
      sbuf(1:m,:,:,:,:) = sbuf(1:m,:,:,:,:)/sdbuf(1:m,:,:,:,:)
    end if
    !---------------------------------------------------------------------------
    ! Take log of unsigned variables
    !---------------------------------------------------------------------------
    if (target%unsigned(iv)) then
      sbuf(1:m,:,:,:,:) = log(sbuf(1:m,:,:,:,:))
    end if
    !---------------------------------------------------------------------------
    ! Interpolate, vectorized
    !---------------------------------------------------------------------------
    do j=1,m
      associate (p1=>pbuf(j,1), p2=>pbuf(j,2), p3=>pbuf(j,3), p4=>pbuf(j,4))
      q1 = 1.0-p1
      q2 = 1.0-p2
      q3 = 1.0-p3
      q4 = 1.0-p4
      tbuf(j) = q4*(q3*(q2*(q1*sbuf(j,1,1,1,1)+p1*sbuf(j,2,1,1,1))   + &
                        p2*(q1*sbuf(j,1,2,1,1)+p1*sbuf(j,2,2,1,1)))  + &
                    p3*(q2*(q1*sbuf(j,1,1,2,1)+p1*sbuf(j,2,1,2,1))   + &
                        p2*(q1*sbuf(j,1,2,2,1)+p1*sbuf(j,2,2,2,1)))) + &
                p4*(q3*(q2*(q1*sbuf(j,1,1,1,2)+p1*sbuf(j,2,1,1,2))   + &
                        p2*(q1*sbuf(j,1,2,1,2)+p1*sbuf(j,2,2,1,2)))  + &
                    p3*(q2*(q1*sbuf(j,1,1,2,2)+p1*sbuf(j,2,1,2,2))   + &
                        p2*(q1*sbuf(j,1,2,2,2)+p1*sbuf(j,2,2,2,2))))
      end associate
    end do
    !---------------------------------------------------------------------------
    ! Take exp of unsigned variables
    !---------------------------------------------------------------------------
    if (target%unsigned(iv)) then
      tbuf(1:m) = exp(tbuf(1:m))
    end if
    !---------------------------------------------------------------------------
    ! Save density values
    !---------------------------------------------------------------------------
    if (iv == 1) then
      tdbuf = tbuf
    end if
    !---------------------------------------------------------------------------
    ! Use density values for per-volume variables
    !---------------------------------------------------------------------------
    if (target%pervolume(iv)) then
      tbuf(1:m) = tbuf(1:m)*tdbuf(1:m)
    end if
  end subroutine mem_interp 
END SUBROUTINE prolong

!===============================================================================
!> Same as above, but for RT_solver type patches
!===============================================================================
SUBROUTINE different_rt (self, target, source, jt, pt)
  class(download_t):: self
  class(patch_t), pointer :: target, source
  !.............................................................................
  integer                 :: ix, iy, iz, iv, jt(2)
  integer, dimension(3)   :: l, u, ii, lb, ub, li, ui
  real(8), dimension(3)   :: pos
  logical                 :: valid
  integer, save           :: itimer=0
  real                    :: pt(2)
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different_rt', itimer=itimer)
  lb = target%mesh%lb
  ub = target%mesh%ub
  call dnload_range (target, source, 1, l, u)
  li = target%mesh%li
  ui = target%mesh%ui
  do iv=1,target%nv
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii=[ix,iy,iz]
          if (self%check_filled) then
            if (filled(ix,iy,iz,iv) > source%level) cycle
          end if
          if (all((ii>=li).and.(ii<=ui))) cycle
          call source%interpolate (target, ii, iv, jt, pt)
          if (self%check_filled) filled(ix,iy,iz,iv) = source%level
        end do
      end do
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE different_rt

!===============================================================================
!> Restrict guard zone info from source patch at a higher level
!===============================================================================
SUBROUTINE different_restrict (self, target, source, jt, pt, only)
  class(download_t):: self
  class(patch_t), pointer :: target, source
  !integer(int8), dimension(:,:,:,:), allocatable:: filled
  integer, optional:: only
  integer                 :: ix, iy, iz, iv, jt(2), iv1, iv2, i, jx, jy, jz
  integer, dimension(3)   :: l, u, ii, li, ui, jj
  real(8), dimension(3)   :: target_pos, source_pos, target_upper_edges, &
                             target_lower_edges, source_upper_edges, &
                             source_lower_edges, le, ue
  logical                 :: valid
  integer, save           :: itimer=0
  real                    :: pt(2), pp(3)
  real(kind=KindScalarVar):: qint
  logical                 :: qmask(3)
  real(8)                 :: start
  integer                 :: cells
  real(8):: rratio(3), dVsource, dVtarget, dWsource, dWtarget, filling_factor, &
            accumulate, filling_total, dWtargeti
  integer:: ilook(3)
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different_restrict', itimer=itimer)
  if (detailed_timer) then
    cells = 0
    start = wallclock()
  end if
  rratio = target%ds / source%ds
  ilook = merge(floor(rratio)/2 + 1, 0, source%n>1)
  dVsource = product(source%ds) ! only valid for Cartesian
  dVtarget = product(target%ds) ! onlt valid for Cartesian
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = target%nv
  end if
  li = target%mesh%li
  ui = target%mesh%ui
  do iv=iv1,iv2
    ! -- `dphi` is part of `mem` when `time_order` = 0; do not download it.
    if (iv == target%idx%dphi) cycle
    ! -- logical mask to track if variable is averaged over a volume or an area
    qmask(:) = .true. 
    if (iv == target%idx%bx) then
      dWsource = source%ds(2) * source%ds(3) ! only valid for Cartesian
      dWtarget = target%ds(2) * target%ds(3) ! only valid for Cartesian
      qmask(1) = .false.
    else if (iv == target%idx%by) then
      dWsource = source%ds(3) * source%ds(1) ! only valid for Cartesian
      dWtarget = target%ds(3) * target%ds(1) ! only valid for Cartesian
      qmask(2) = .false.
    else if (iv == target%idx%bz) then
      dWsource = source%ds(1) * source%ds(2) ! only valid for Cartesian
      dWtarget = target%ds(1) * target%ds(2) ! only valid for Cartesian
      qmask(3) = .false.
    else
      dWsource = dVsource
      dWtarget = dVtarget
    end if
    dWtargeti = 1.0 / dWtarget
    call dnload_range (target, source, iv, l, u)
    associate (q1 => source%mem(:,:,:,iv,jt(1),1), &
               q2 => source%mem(:,:,:,iv,jt(2),1))
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii = [ix,iy,iz]
          if (iv==1) cells = cells+1
          ! prefer same level for guard zone filling
          if (self%check_filled) then
            if (filled(ix,iy,iz,iv) == target%level .and. &
                (any(ii < target%li) .or. any(ii > target%ui))) cycle
          end if
          target_pos = target%myposition (ii, iv)
          target_lower_edges = target_pos - target%ds * 0.499999_8
          target_upper_edges = target_pos + target%ds * 0.499999_8
          ! skip cells that are not *entirely* covered by the finer patch
          if (any(target_lower_edges < source%llc_nat)) cycle
          if (any(target_upper_edges > (source%llc_nat+source%size))) cycle
          call source%index_space (target_pos, iv, jj, pp)
          accumulate = 0.0
          filling_total = 0.0
          do jz=jj(3)-ilook(3),jj(3)+ilook(3)
            do jy=jj(2)-ilook(2),jj(2)+ilook(2)
              do jx=jj(1)-ilook(1),jj(1)+ilook(1)
                ! skip ghost cells of source
                if (any([jx,jy,jz] < source%li .or. [jx,jy,jz] > source%ui)) cycle
                ! find intersection between source and target cells
                source_pos = source%myposition([jx,jy,jz],iv)
                source_lower_edges = source_pos - source%ds * 0.499999_8
                source_upper_edges = source_pos + source%ds * 0.499999_8
                ue = min(source_upper_edges, target_upper_edges)
                le = max(source_lower_edges, target_lower_edges)
                ! filling factor; if zero, then no overlap.
                ! the logical mask is used to produce an area filling factor for
                ! mag. fields (but a volume filling factor otherwise)
                filling_factor = product(max(ue - le, 0.0),mask=qmask) 
                filling_total = filling_total + filling_factor
                ! `accumulate` should be multipled by `dWsource`, but 
                ! `filling_factor` should be divided by `dWsource`, so they cancel.
                if (filling_factor > 0.0) then
                  accumulate = accumulate &
                    + (q1(jx,jy,jz) * pt(1) + q2(jx,jy,jz) * pt(2)) * &
                    filling_factor
                end if
              end do
            end do
          end do
          ! replace a portion of the current cell with values from the finer patch
          if (filling_total > 0.0) then
            target%mem(ix,iy,iz,iv,target%it,1) = accumulate * dWtargeti &
              + target%mem(ix,iy,iz,iv,target%it,1) &
              * (1.0 - filling_total * dWtargeti)
            if (self%check_filled) then
              filled(ix,iy,iz,iv) = source%level
            end if
          end if
        end do
      end do
    end do
    end associate
    if (verbose > 0 .or. target%id==id_debug) then
      if (iv == 1 .or. verbose > 2) &
        write (io_unit%log,2) target%id, target%level, source%id, source%level, &
        'diff_rest picked up', &
        product(u-l+1), ' values, l, u, mind, pt =', l, u, &
        minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,target%it,1)), pt
      2 format(2(i6,i3),2x,a,i6,a,2(2x,3i4),1p,e14.4,2x,0p2f6.3)
    end if
  end do
  if (detailed_timer) then
    write (stdout,'(a,6p,f6.3," mus/cell")') &
      "download_t%different_restrict:", (wallclock()-start)/cells
  end if
  call trace%end (itimer)
END SUBROUTINE different_restrict

!===============================================================================
!> Restrict guard zone info from source patch at a higher level
!===============================================================================
SUBROUTINE restrict_centered (self, target, source, jt, pt, only)
  class(download_t):: self
  class(patch_t), pointer :: target, source
  !integer(int8), dimension(:,:,:,:), allocatable:: filled
  integer, optional:: only
  integer                 :: ix, iy, iz, iv, jt(2), iv1, iv2, i, jx, jy, jz
  integer, dimension(3)   :: l, u, ii, li, ui, jj
  real(8), dimension(3)   :: target_pos, source_pos, target_upper_edges, &
                             target_lower_edges, source_upper_edges, &
                             source_lower_edges, le, ue
  logical                 :: valid
  integer, save           :: itimer=0
  real                    :: pt(2), pp(3)
  real(kind=KindScalarVar):: qint
  logical                 :: qmask(3)
  real(8):: rratio(3), dVsource, dVtarget, dWsource, dWtarget, filling_factor, &
            accumulate(target%nv), filling_total, dWtargeti
  integer:: ilook(3)
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different_restrict', itimer=itimer)
  rratio = target%ds / source%ds
  ilook = merge(floor(rratio)/2 + 1, 0, source%n>1)
  dVsource = product(source%ds) ! only valid for Cartesian
  dVtarget = product(target%ds) ! onlt valid for Cartesian
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = target%nv
  end if
  li = target%mesh%li
  ui = target%mesh%ui
    dWsource = dVsource
    dWtarget = dVtarget
    dWtargeti = 1.0 / dWtarget
    call dnload_range (target, source, 1, l, u)
    associate (q1 => source%mem(:,:,:,:,jt(1),1), q2 => source%mem(:,:,:,:,jt(2),1))
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii = [ix,iy,iz]
          ! prefer same level for guard zone filling
          if (self%check_filled) then
            if (filled(ix,iy,iz,iv) == target%level .and. &
                (any(ii < target%li) .or. any(ii > target%ui))) cycle
          end if
          target_pos = target%myposition (ii, iv)
          target_lower_edges = target_pos - target%ds * 0.499999_8
          target_upper_edges = target_pos + target%ds * 0.499999_8
          ! skip cells that are not *entirely* covered by the finer patch
          if (any(target_lower_edges < source%llc_nat)) cycle
          if (any(target_upper_edges > (source%llc_nat+source%size))) cycle
          call source%index_space (target_pos, iv, jj, pp)
          accumulate = 0.0
          filling_total = 0.0
          do jz=jj(3)-ilook(3),jj(3)+ilook(3)
            do jy=jj(2)-ilook(2),jj(2)+ilook(2)
              do jx=jj(1)-ilook(1),jj(1)+ilook(1)
                ! skip ghost cells of source
                if (any([jx,jy,jz] < source%li .or. [jx,jy,jz] > source%ui)) cycle
                ! find intersection between source and target cells
                source_pos = source%myposition([jx,jy,jz],iv)
                source_lower_edges = source_pos - source%ds * 0.499999_8
                source_upper_edges = source_pos + source%ds * 0.499999_8
                ue = min(source_upper_edges, target_upper_edges)
                le = max(source_lower_edges, target_lower_edges)
                ! filling factor; if zero, then no overlap.
                ! the logical mask is used to produce an area filling factor for
                ! mag. fields (but a volume filling factor otherwise)
                filling_factor = product(max(ue - le, 0.0),mask=qmask)
                filling_total = filling_total + filling_factor
                ! `accumulate` should be multipled by `dWsource`. but
                ! `filling_factor` should be divided by `dWsource`, so they cancel.
                if (filling_factor > 0.0) then
                  do iv=iv1,iv2
                    accumulate(iv) = accumulate(iv) &
                             + (q1(jx,jy,jz,iv) * pt(1) + q2(jx,jy,jz,iv) &
                             * pt(2)) * filling_factor
                  end do
                end if
              end do
            end do
          end do
          ! replace a portion of the current cell with values from the finer patch
          if (filling_total > 0.0) then
            do iv=iv1,iv2
              target%mem(ix,iy,iz,iv,target%it,1) = accumulate(iv) * dWtargeti &
                 + target%mem(ix,iy,iz,iv,target%it,1) &
                 * (1.0 - filling_total * dWtargeti)
            end do
            if (self%check_filled) then
              filled(ix,iy,iz,iv) = source%level
            end if
          end if
        end do
      end do
    end do
    end associate
  call trace%end (itimer)
END SUBROUTINE restrict_centered

!===============================================================================
SUBROUTINE set_verbose (new)
  integer:: new
  verbose = new
END SUBROUTINE set_verbose

!===============================================================================
!> Initialise downloading and interpolating options.
!===============================================================================
SUBROUTINE test (self)
  class(download_t):: self
  !.............................................................................
  class(patch_t), pointer :: source(:), target, src
  integer, dimension(3)  :: l, u, o, ii
  integer                :: i, j, k, iv, dir, id_start, id_end, jt(2)
  real                   :: pt(2)
  real(8), dimension(3)  :: pos
  real(4), dimension(3)  :: pp
  class(mesh_t), pointer :: mesh
  logical                :: ok, allok
  !-----------------------------------------------------------------------------
  return
  call trace%begin ('download_t%test')
  id_start = mpi%id%update(0)
  self%interpolate => SelectInterpolator(-1)
  !-----------------------------------------------------------------------------
  ! Generate meshes with three resolutions, differing with a factor of three,
  ! with and w/o staggering
  !-----------------------------------------------------------------------------
  allocate (target)
  call target%setup (size=0.1_8*[1,1,1], n=[16,16,16], nv=2, nt=2, nw=1)
  filled = -1
  !-----------------------------------------------------------------------------
  ! Loop over the three resolutions for source patches
  !-----------------------------------------------------------------------------
  allocate (source(3))
  do i=1,3
    call source(i)%setup (size=0.1_8*[1,1,1]*(shared%f_zoom)**(i-2), &
      n=[16,16,16], nv=2, nt=2, nw=1)
  end do
  !---------------------------------------------------------------------------
  ! Loop over source size
  !---------------------------------------------------------------------------
  target%level = 2
  do i=1,3
    source(i)%level = i
    !---------------------------------------------------------------------------
    ! Loop over no-stagger / stagger
    !---------------------------------------------------------------------------
    do iv=1,2
      do dir=1,3
        mesh => source(i)%mesh(dir)
        mesh%h(iv)    = -0.5*(iv-1)
        mesh => target%mesh(dir)
        mesh%h(iv)    = -0.5*(iv-1)
      end do
      if (io%verbose>0 .or. io_unit%do_validate) then
        mesh => source(i)%mesh(1)
        write (io_unit%log,*) 'source size, stagger:', source(i)%size(1), mesh%h(iv)
      end if
      !-------------------------------------------------------------------------
      ! Loop over placement of the source mesh: ol, il, c, ir, or
      !-------------------------------------------------------------------------
      src => source(i)
      do j=1,5
        call set_position (src, target, j)
        if (i==2.and.(j==2.or.j==4)) cycle
        if (src%level==target%level) then
          call dnload_range (target, src, iv, l, u, o)
        else
          call dnload_range (target, src, iv, l, u)
        end if
        do k=l(1),u(1)
          do dir=1,3
            mesh => target%mesh(dir)
            pos(dir) = target%position(dir) + mesh%r(k) + mesh%h(iv)*target%ds(dir)
          end do
          call source(i)%index_space (pos, iv, ii, pp)
          if (io%verbose>0 .or. io_unit%do_validate) then
            if (src%level==target%level) then
              print 1, i, iv, j, "k:", k, "o:", o(1), "ii:", ii(1), "pp:", pp(1)
              1 format("download_t%test",3i3,3(2x,a,i4),2x,a,f5.2)
            else
              print 2, i, iv, j, "k:", k, "o: N/A", "ii:", ii(1), "pp:", pp(1)
              2 format("download_t%test",3i3,2x,a,i4,2x,a,2x,a,i4,2x,a,f5.2)
            end if
          end if
        end do
        if (io%verbose>0 .or. io_unit%do_validate) &
          write (io_unit%log,*)' '
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Test actual download values
  !-----------------------------------------------------------------------------
  allok = .true.
  do i=1,3
    src => source(i)
    do j=1,5
      call set_position (src, target, j)
      if (i==2.and.(j==2.or.j==4)) cycle
      if (io%verbose>1) & 
        print "('source size:',i2,3x,'position:',i2,2x,'levels:',2i2)", i, j, &
          src%level, target%level
      call set_values (src)
      call set_values (target)
      if (src%level==target%level) then
        call src%time_interval (target%time, jt, pt)
        call self%same (target, src, jt, pt, all_cells=.false.)
        if (self%check_filled) then
          if (any(filled < 0)) then
            write (io_unit%log,*) target%id, &
              'WARNING: some cells not filled by download_t%same'
          end if
        end if
      else
        !call different_prolong (self, target, src, all_cells=.false.)
        call src%time_interval (target%time, jt, pt)
        call prolong (self, target, src, jt, pt, all_cells=.false.)
      end if
      call check_values (target, ok)
      allok = allok.and.ok
    end do
  end do
  if (io%master) write (io_unit%log,*)io%hl
  if (allok) then
    if (io%master) then
      write (io_unit%log,*)'download_t%test passed'
      write (io_unit%log,*)io%hl
    end if
  else
    write (io_unit%log,*)mpi%rank,'download_t%test failed'
    if (io%master) write (io_unit%log,*)io%hl
    call mpi%abort('download_t%test::check_values')
  end if
  !-----------------------------------------------------------------------------
  ! Clean up
  !-----------------------------------------------------------------------------
  deallocate (filled)
  call target%dealloc
  do i=1,3
    call source(i)%dealloc
  end do
  !-----------------------------------------------------------------------------
  ! Set the MPI-global id counter back to what it was
  !-----------------------------------------------------------------------------
  id_end = mpi%id%update(0)
  id_start = mpi%id%update (id_start-id_end)
  deallocate (target, source)
  self%interpolate => SelectInterpolator(self%order_interpolator)
  call trace%end()
contains
!===============================================================================
!> Set position systematically
!===============================================================================
  subroutine set_position (source, target, j)
    type(patch_t) :: source, target
    integer       :: j
    !---------------------------------------------------------------------------
    select case (j)
    case (1)
      source%position = target%position - 0.5_8*(source%size+target%size) + &
        2.*abs(source%ds-target%ds)
    case (2)
      source%position = target%position - 0.5_8*(source%size+target%size) + &
        2.*abs(source%ds-target%ds)
    case (3)
      source%position = target%position
    case (4)
      source%position = target%position + 0.5_8*(source%size+target%size) - &
        2.*abs(source%ds-target%ds)
    case (5)
      source%position = target%position + 0.5_8*(source%size+target%size) - &
        2.*abs(source%ds-target%ds)
    end select
    source%mesh%p = source%position
    source%t = 0d0; source%time=0d0; source%iit=1
    target%t = 0d0; target%time=0d0; target%iit=1
end subroutine set_position
!===============================================================================
!> Set some values to recover
!===============================================================================
subroutine set_values (patch)
  type(patch_t)          :: patch
  class(mesh_t), pointer :: m1, m2, m3
  integer                :: ix, iy, iz, iv
  real                   :: x, y, z, v
  !-----------------------------------------------------------------------------
  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)
  do iv=1,2
    patch%mem(:,:,:,iv,:,:) = 0.0
    do iz=m3%li,m3%ui
      z = m3%p + m3%r(iz) + m3%h(iv)*m3%d
      do iy=m2%li,m2%ui
        y = m2%p + m2%r(iy) + m2%h(iv)*m2%d
        do ix=m1%li,m1%ui
          x = m1%p + m1%r(ix) + m1%h(iv)*m1%d
          v = x + y*2.0 + z*3.0
          patch%mem(ix,iy,iz,iv,:,:) = v
        end do
      end do
    end do
  end do
end subroutine set_values
!===============================================================================
!> Check values, including in guard zones
!===============================================================================
subroutine check_values (patch, ok)
  type(patch_t)          :: patch
  class(mesh_t), pointer :: m1, m2, m3
  logical                :: ok
  integer                :: ix, iy, iz, iv, n(3)
  real                   :: x, y, z, v, w, eps
  !-----------------------------------------------------------------------------
  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)
  eps=0.1*max(m1%d,m2%d,m3%d)
  n = 0
  do iv=1,2
    do iz=m3%lb,m3%ub
      z = m3%p + m3%r(iz) + m3%h(iv)*m3%d
      do iy=m2%lb,m2%ub
        y = m2%p + m2%r(iy) + m2%h(iv)*m2%d
        do ix=m1%lb,m1%ub
          x = m1%p + m1%r(ix) + m1%h(iv)*m1%d
          v = x + y*2.0 + z*3.0
          w = patch%mem(ix,iy,iz,iv,1,1)
          if (w == 0.0) then
            n(1) = n(1)+1
          else
            if (abs(w-v) > eps) then
              print "('error: ix, iy, iz, iv, value:',4i4,1p,2e15.5)", &
               ix, iy, iz, iv, w, v
              n(3) = n(3)+1
            else
              n(2) = n(2)+1
            end if
          end if
        end do
      end do
    end do
  end do
  if (io%verbose>0) & 
    write (io_unit%log,*)'check_values: n_zero, n_ok, n_nok =', n
  if (n(3) == 0) then
    ok = .true.
  else
    ok = .false.
    write (io_unit%log,*)'check_values: n_zero, n_ok, n_nok =', n
  end if
end subroutine check_values
!===============================================================================
END SUBROUTINE test

END MODULE
