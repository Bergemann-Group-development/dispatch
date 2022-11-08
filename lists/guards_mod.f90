!===============================================================================
!> Handles guard zone interpolations and boundare conditions.  For each
!> nbor node, this gets called by the halo_mod.   Since it inherits the
!> nbor_t data type, it also has links to target and source.
!>
!> The responsibility of this module is to decide when to do BCs and
!> guard zone interpolations, respectively, and to decide the order of
!> interpolations in x, y, and z, such that the perpendicular direction
!> is always interpolated immediately after the time interpolation (or the
!> other way around, possibly)
!===============================================================================
MODULE guards_mod
  USE io_mod
  USE os_mod
  USE mpi_mod
  USE trace_mod
  USE omp_timer_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE boundary_mod
  USE interp_mod
  USE mesh_mod
  USE conserve_mod
  USE validate_mod
  USE math_mod
  USE shared_mod
  USE bits_mod
  implicit none
  private
  type, public, extends(nbor_t):: guards_t
    integer:: verbose = 0
    integer:: n_count(4)=0
    integer:: lt(3)=0, ut(3)=0
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: update
    procedure:: source_to_target
    procedure:: same
    procedure:: prolong
    procedure:: restrict
    procedure:: info
    procedure:: differ
    procedure:: tilted
    procedure:: nbor_info
    procedure, nopass:: alloc_nbor
  end type
  integer:: verbose=0
  logical:: detailed_timer = .false.
  type(guards_t), public:: guards
CONTAINS

!===============================================================================
!> Initialize the nbor node
!===============================================================================
SUBROUTINE pre_init (self)
  class(guards_t) :: self
  !.............................................................................
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /guards_params/ verbose, detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('guards_t%pre_init')
  call interp%init
  !-----------------------------------------------------------------------------
  ! Read namelist parameters
  !-----------------------------------------------------------------------------
  if (first_time) then
    call io%lock%set ('guards_t%pre_init')
    if (first_time) then
      rewind (io%input)
      read (io%input, guards_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, guards_params)
      end if
      call io%namelist_warning ('guards_params', iostat)
      write (io%namelist, guards_params)
      first_time = .false.
    end if
    call io%lock%unset ('guards_t%pre_init')
  end if
  self%verbose = verbose
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the nbor node
!===============================================================================
SUBROUTINE init (self)
  class(guards_t):: self
  !-----------------------------------------------------------------------------
  ! Setup nbor distance and relations, in case these differ for patch_t
  !-----------------------------------------------------------------------------
  call trace%begin ('guards_t%init')
  call self%pre_init ()
  call self%nbor_t%init ()
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> The face update procedure organizes calls to the interpolation routines,
!> with the target and source patch data types served in the argument list
!===============================================================================
SUBROUTINE update (self, all, only)
  class(guards_t)        :: self
  logical, optional      :: all
  integer, optional      :: only
  !-----------------------------------------------------------------------------
  call trace%begin ('guards_t%update')
  call source_to_target (self, self%source, self%target, all)
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!===============================================================================
SUBROUTINE source_to_target (self, source, target, all_tgt, all_src, n_op, &
                             only, nbstep)
  class(guards_t)        :: self
  class(task_t)          :: source, target
  logical, optional      :: all_tgt, all_src
  integer, optional      :: n_op(3)
  integer, optional      :: only, nbstep
  !.............................................................................
  character(10)          :: type
  logical                :: tilted
  integer                :: m(6), tid, sid, unit
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('guards_t%source_to_target', itimer=itimer)
  self%distance = source%distance (target)
  tid = target%id
  sid = source%id
  !.............................................................................
  if (verbose > 2) then
    unit = stdout
  else
    unit = stdlog
  end if
  if (verbose > 1) then
    write (unit,'(a,2(i6,":",i2.2),3f7.3)') &
      'guards_t%source_to_target: src,lvl,tgt,lvl,dist=', &
      source%id, source%level, target%id, target%level, &
      self%distance/target%size
  end if
  !.............................................................................
  select type (target)
  class is (patch_t)
    m = shape(target%mem)
    select type (source)
    class is (patch_t)
      !-------------------------------------------------------------------------
      ! Initialize local test if not done
      !-------------------------------------------------------------------------
      self%distance = target%distance_to (source)
      tilted = associated (source%erot)
      if (source%level < shared%levelmin) then
        type = 'none'
      else if (target%level < shared%levelmin) then
        type = 'none'
        if (source%level == target%level+1) then
          call self%restrict (target, source, all_tgt=.true., nbstep=nbstep)
          type = 'restrict1'
        end if
      else if (tilted) then
        type = 'tilted'
        call self%tilted (target, source, all_tgt, all_src, only)
      else if (source%level == target%level) then
        type = 'same'
        call self%same (target, source, all_tgt, all_src, only, nbstep)
        if (present(n_op)) n_op  (1) = n_op(1) + 1
      else if (source%level == target%level-1) then
        type = 'prolong'
        call self%prolong (target, source, all_tgt, all_src, only, nbstep)
        if (present(n_op)) n_op(2) = n_op(2) + 1
      else if (source%level == target%level+1) then
        type = 'restrict'
        !-----------------------------------------------------------------------
        ! FIXME: For now, the conserve%restrict() can only handle interior cells
        ! Note: Restrict operations are either done only in guard zones, or if
        ! the source is contained in the target then all cells should be used,
        ! so there is no need to test on all_tgt here:
        !-----------------------------------------------------------------------
        if (associated(source%parent)) then
          if (source%parent%id == target%id .and. &
              .not. present(all_src)) then
            call self%restrict (target, source, all_tgt, all_src, only, nbstep)
            !call conserve%restrict (target, source, all_tgt, all_src, only)
            type = 'restrict-c'
          else
            if (source%parent%id == target%id) then
              write (stdlog,'(a,2(i6,":",i2.2,2x),2l1)') &
                'guards_t%restrict: warning, source,target,all=', &
                source%id, source%level, target%id, target%level, &
                present(all_tgt), present(all_src)
            end if
            call self%restrict (target, source, all_tgt, all_src, only, nbstep)
            !call conserve%restrict (target, source, all_tgt, all_src, only)
          end if
        else
          call source%log ('missing parent')
          write (stderr,'(/,2(1x,a,i6,":",i2.2,3i3,2x))') &
            'source:', source%id, source%key, &
            'target:', target%id, target%key
          call io%warning ('guards_t%source_to_target: missing source%parent')
        end if
        if (present(n_op)) n_op(3) = n_op(3) + 1
      end if
      if (verbose > 1 .or. target%id == io%id_debug) then
        if (target%idx%d > 0) then
          write (unit,'(a8,2x,2(a,i6,":",i2.2,3x),a,1p,e12.2)') type, &
            'tgt:', tid, target%level, &
            'src:', sid, source%level, &
            'dmin:', minval(target%mem(:,:,:,target%idx%d,target%it,1))
        end if
      end if
    end select
    !...........................................................................
    if (verbose > 1) then
      write (unit,'(a,2x,i6,":",i2.2,f12.6,2x,i6,":",i2.2,2x,3f8.3,2x,a)') &
        'guards_t%source_to_target: tgt,time,src,dist,type=', &
        target%id, target%level, target%time, source%id, source%level, &
        self%distance/target%size, trim(type)
    end if
    !...........................................................................
  end select
  if (detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE source_to_target

!===============================================================================
!> Interpolate in time, leaving the result in the target.   This procedure
!> ignores the values of all_tgt and all_src, always using as much as possible
!> from the internal parts of the source, corrsponging to all_tgt=T, all_scr=F
!===============================================================================
SUBROUTINE same (self, tgt, src, all_tgt, all_src, only, nbstep)
  class(guards_t)                      :: self
  class(patch_t)                       :: tgt, src
  logical, optional                    :: all_tgt, all_src
  integer, optional                    :: only, nbstep
  !.............................................................................
  class(mesh_t), pointer, dimension(:) :: tm, sm
  type(meshp_t)                        :: om(3)
  real(8), dimension(3)                :: sm_p, tl, tu, sl, su, fl, fu, pl, pu,d
  real, dimension(:,:,:,:), pointer    :: mem
  !contiguous                           :: mem
  real, allocatable                    :: tmp(:,:,:,:)
  integer                              :: i, m(6)
  integer, dimension(3)                :: ls, us, lt, ut
  integer, save                        :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('guards_t%same', itimer=itimer)
  call self%info ('same', tgt, src, all_tgt, all_src, only)
  tm => tgt%mesh                             ! tgt mesh
  sm => src%mesh                             ! src mesh
  do i=1,3
    om(i)%p => sm(i)
  end do
  d = tgt%distance_to (src)                  ! wrapped vector (src-tgt)
  tl =   - tm%s*0.5d0 - tm%ng*tm%d           ! tgt lower corner
  tu =   + tm%s*0.5d0 + tm%ng*tm%d           ! tgt upper corner
  sl = d - sm%s*0.5d0                        ! src lower corner
  su = d + sm%s*0.5d0                        ! src upper corner
  fl = max(tl,sl)                            ! overlap lower corner
  fu = min(tu,su)                            ! overlap upper corner
  call io%assert (all(fu > fl), 'guards_t%differ: corner reversal')
  ls = ceiling((fl-d)/sm%d + sm%o)      ! src lower index
  us = floor  ((fu-d)/sm%d + sm%o)      ! src upper index
  lt = ceiling((fl  )/tm%d + tm%o)      ! tgt lowwe index
  ut = floor  ((fu  )/tm%d + tm%o)      ! tgt upper index
  if (verbose > 2) then
    write (stdlog,'(a,2i6,4(2x,3i7))') 'guards_t%same: tgt,src,ls,us,lt,ut=', &
      tgt%id, src%id, ls, us, lt, ut
    call tgt%stats ('guards_t%same tgt')
    call src%stats ('guards_t%same src')
  end if
  call io%assert (all(lt >= tm%lb), 'guards_t%same: lt breach')
  call io%assert (all(ut <= tm%ub), 'guards_t%same: ut breach')
  call io%assert (all(ls >= sm%lb), 'guards_t%same: ls breach')
  call io%assert (all(us <= sm%ub), 'guards_t%same: us breach')
  !-----------------------------------------------------------------------------
  ! Time interpolate in the target region
  !-----------------------------------------------------------------------------
  if (tgt%nt==1 .and. src%nt==1) then
    call interp%no_t   (tgt, src, ls, us, src%mem, tgt%mem, lt, ut, only=only)
  else
    m = shape(tgt%mem)
    mem => tgt%mem(:,:,:,:,tgt%it,1)
    call interp%in_t   (tgt, src, ls, us, src%mem, mem, lt, ut, only=only, &
                          extrap=tgt%allow_extrapolate, nbstep=nbstep)
  endif
  if (associated(tgt%filled)) &
    tgt%filled (lt(1):ut(1),lt(2):ut(2),lt(3):ut(3)) = 10*src%level
  !.............................................................................
  if (verbose > 1 .or. tgt%id == io%id_debug) then
    if (tgt%idx%d > 0) then
      write (stdlog,'(a,2(i6,"::",i2.2),2(2x,3i3),1p,2g14.3)') &
        '  same: tid,sid,l,u,min=', &
        tgt%id, tgt%level, src%id, src%level, lt, ut, &
        minval(tgt%mem(lt(1):ut(1),lt(2):ut(2),lt(3):ut(3), &
          tgt%idx%d,tgt%it,1)), &
        minval(tgt%mem(:,:,:,tgt%idx%d,tgt%it,1))
    else
      write (stdlog,'(a,2(i6,"::",i2.2),2(2x,3i3))') &
        '  same: tid,sid,l,u,min=', &
        tgt%id, tgt%level, src%id, src%level, lt, ut
    end if
    flush (stdout)  
  end if
  !.............................................................................
  call trace%end (itimer)
END SUBROUTINE same

!===============================================================================
SUBROUTINE prolong (self, tgt, src, all_tgt, all_src, only, nbstep)
  class(guards_t)   :: self
  class(patch_t)    :: tgt, src
  logical, optional :: all_tgt, all_src
  integer, optional :: only, nbstep
  !-----------------------------------------------------------------------------
  call trace%begin ('guards_t%prolong')
  call self%info ('prolong', tgt, src, all_tgt, all_src, only)
  call self%differ (tgt, src, all_tgt, all_src, only, nbstep)
  call trace%end ()
END SUBROUTINE prolong

!===============================================================================
SUBROUTINE restrict (self, tgt, src, all_tgt, all_src, only, nbstep)
  class(guards_t)   :: self
  class(patch_t)    :: tgt, src
  logical, optional :: all_tgt, all_src
  integer, optional :: only, nbstep
  !-----------------------------------------------------------------------------
  call trace%begin ('guards_t%restrict')
  call self%info ('restrict', tgt, src, all_tgt, all_src, only)
  call self%differ (tgt, src, all_tgt, all_src, only, nbstep)
  call trace%end ()
END SUBROUTINE restrict

!===============================================================================
!===============================================================================
SUBROUTINE info (self, label, target, source, all_tgt, all_src, only)
  class(guards_t)   :: self
  class(patch_t)    :: target, source
  character(len=*)  :: label
  logical, optional :: all_tgt, all_src
  integer, optional :: only
  integer           :: o
  !.............................................................................
  if (verbose > 1) then
    if (present(only)) then
      o = only
    else
      o = 0
    end if
    write (stdlog,'(2(f12.6,2x,a,2x,i6,":",i2.2),2(2x,3f7.3),2x,2l1,i3)') &
      wallclock(), 'guards_t%'// &
      trim(label)//': tgt,time,src,dist,all_tgt,all_src=', &
      target%id, target%level, &
      target%time, '  src:', source%id, source%level, &
      self%distance/target%size, source%distance(target)/target%size, &
      present(all_tgt), present(all_src), o
  end if
  !.............................................................................
END SUBROUTINE info

!===============================================================================
!> Interpolate in time, leaving the result in the target.  Unless all_tgt is
!> set the target memory will only be changed in the guard zones, which is of
!> no concern to other tasks, so then there is no need to lock the memory slot.
!===============================================================================
SUBROUTINE differ (self, tgt, src, all_tgt, all_src, only, nbstep)
  class(guards_t)                      :: self
  class(patch_t)                       :: tgt, src
  logical, optional                    :: all_tgt, all_src
  integer, optional                    :: only, nbstep
  !.............................................................................
  class(task_t), pointer               :: nbtask
  class(mesh_t), pointer, dimension(:) :: tm, sm
  type(meshp_t)                        :: om(3)
  real(8), dimension(3)                :: tl, tu, sl, su, fl, fu, pl, pu, sm_p
  real(8), dimension(3)                :: d
  logical, dimension(3)                :: below, above
  integer, dimension(3)                :: ls, us, lt, ut, sm_lb, sm_ub
  real, dimension(:,:,:,:), pointer    :: tgt0, tgt1, tgt2, mem
  contiguous                           :: tgt0, tgt1, tgt2, mem
  logical                              :: src_all, tgt_int
  integer                              :: i, m(5), unit, sid, tid
  integer, save                        :: itimer=0
  !-----------------------------------------------------------------------------
  ! Use all source cells if the all_src is present, or the rarget (child) is
  ! taking the first step.
  !-----------------------------------------------------------------------------
  call trace%begin ('guards_t%differ', itimer=itimer)
  tid = tgt%id
  sid = src%id
  unit = merge (stdout, stdlog, verbose > 2)
  src_all = present(all_src) !.or. tgt%istep==0  ! use full source region
 !src_all = src_all .or. tgt%istep==0            ! first child update
  !-----------------------------------------------------------------------------
  ! Produce values also in internal cells if specifically ordered to, or if
  ! the source is at a higher level, or if a child is taking the first step.
  !-----------------------------------------------------------------------------
  tgt_int = tgt%istep==0                        ! first child update
  tgt_int = tgt_int .or. src%level > tgt%level  ! restrict
  tgt_int = tgt_int .and. tgt%time > 0.0        ! FIXME: restarts!
  tgt_int = tgt_int .or. present(all_tgt)       ! produce full target
  !.............................................................................
  if (verbose > 1) then
    write (unit,'(f12.6,a,2(i6,":",i2.2,2x),4l1)') &
      wallclock(), pad_trace(1:indent)// &
      'guards_t%differ: src,tgt,all_tgt,all_src,tgt_int,src_all=', &
      src%id, src%level, tgt%id, tgt%level, &
      present(all_tgt), present(all_src), tgt_int, src_all
    if (unit == stdout) then
      write (unit,'(f12.6,2x,a,2(i6,":",i2.2,2x),4l1)') &
        wallclock(), &
        'guards_t%differ: src,tgt,all_tgt,all_src,tgt_int,src_all=', &
        src%id, src%level, tgt%id, tgt%level, &
        present(all_tgt), present(all_src), tgt_int, src_all
    end if
  end if
  !.............................................................................
  tm => tgt%mesh                                ! target mesh
  sm => src%mesh                                ! source mesh
  do i=1,3
    om(i)%p => src%mesh(i)                      ! output mesh
  end do
  !-----------------------------------------------------------------------------
  ! Here we need to either work exclusively with nbor%distance, or else rely
  ! exclusively on task %position
  !-----------------------------------------------------------------------------
  tl = - tm%s*0.5d0 - tm%ng*tm%d                ! tgt lower corner
  tu = + tm%s*0.5d0 + tm%ng*tm%d                ! tgt upper corner
  d = tgt%distance_to (src)                     ! vector (src%p-tgt%p)
  if (src_all) then                             ! new child
    sl = d - sm%s*0.5d0 - sm%ng*sm%d            ! src lower corner
    su = d + sm%s*0.5d0 + sm%ng*sm%d            ! src upper corner
    sm_lb = sm%lb
    sm_ub = sm%ub
  else
    sl = d - sm%s*0.5d0                         ! src lower corner
    su = d + sm%s*0.5d0                         ! src upper corner
    sm_lb = sm%li
    sm_ub = sm%ui
  end if
  !-----------------------------------------------------------------------------
  ! The (tl,tu) and (sl,sy) boxes include all cells -- the last cell centers
  ! are half a cell inside the corners, so taking the max() and min() of the
  ! corners gives a box with a half cell margin to cell centers.
  !-----------------------------------------------------------------------------
  fl = max(tl,sl)                               ! overlap lower corner
  fu = min(tu,su)                               ! overlap upper corner
  if (.not. all(fu > fl)) &
    call io%abort ('guards_t%differ: corner reversal')
  !---------------------------------------------------------------------------
  ! Restrict operation:  The source cells are smaller, and are hence closer
  ! to the cell corners.
  !---------------------------------------------------------------------------
  lt = ceiling((fl  )/tm%d + tm%o)              ! tgt lower index
  ut = floor  ((fu  )/tm%d + tm%o)              ! tgt upper index
  fl = (lt-tm%o)*tm%d                           ! floating lower
  fu = (ut-tm%o)*tm%d                           ! floating upper
  ls = floor  ((fl-d)/sm%d + sm%o)              ! src lower index
  us = ceiling((fu-d)/sm%d + sm%o)              ! src upper index
  ls = max(ls,sm_lb)                            ! bounds
  us = min(us,sm_ub)                            ! bounds
  lt = max(lt,tm%lb)                            ! bounds
  ut = min(ut,tm%ub)                            ! bounds
  if (src%is_set (bits%root_grid)) then
    lt = tm%lb
    ut = tm%ub
  end if
  !-----------------------------------------------------------------------------
  ! By construction, all indices should already be within bounds
  !-----------------------------------------------------------------------------
  if (src_all) then
    call io%assert (all(ls >= sm%lb), 'guards_t%differ: ls index error')
    call io%assert (all(us <= sm%ub), 'guards_t%differ: us index error')
  else
    call io%assert (all(ls >= sm%li), 'guards_t%differ: ls index error')
    call io%assert (all(us <= sm%ui), 'guards_t%differ: us index error')
  end if
  call io%assert (all(lt >= tm%lb), 'guards_t%differ: lt index error')
  call io%assert (all(ut <= tm%ub), 'guards_t%differ: ut index error')
  call io%assert (all(us >= ls), 'guards_t%differ: ls/us index reversal')
  call io%assert (all(ut >= lt), 'guards_t%differ: lt/ut index reversal')
  !-----------------------------------------------------------------------------
  ! Decide relative positioning, in terms of axes
  !-----------------------------------------------------------------------------
  below = src%mesh%p < tgt%mesh%p - tgt%mesh%s/2d0
  above = src%mesh%p > tgt%mesh%p + tgt%mesh%s/2d0
  !.............................................................................
  if (verbose > 1) then
    write (stdlog,'(f12.6,a,2(i6,":",i2.2,1x),2(1x,2(1x,3i7)),3(2x,3l1))') &
      wallclock(), pad_trace(1:indent)// &
      'guards_t%differ: src,tgt,ls,us,lt,ut,bel,abv=', &
      src%id, src%level, tgt%id, tgt%level, ls, us, lt, ut,  &
      below, above, tgt_int, src_all
    if (unit == stdout) then
      write (unit,2) wallclock(), &
        'guards_t%differ: src,tgt,ls,us,lt,ut,bel,abv=', &
        src%id, src%key, tgt%id, tgt%key, tgt%distance_to(src)/tgt%size, &
        ls, us, lt, ut, below, above, tgt_int, src_all
    2 format(f12.6,2x,a,2(i6,":",i2.2,3i4,1x),3f10.4,2(1x,2(1x,3i7)),3(2x,3l1))
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays with common size
  !-----------------------------------------------------------------------------
  nullify (tgt0, tgt1)
  m = shape(tgt%mem(:,:,:,:,:,1))
  call os%allocate (tgt0, m, 'tgt0')
  call os%allocate (tgt1, m, 'tgt1')
  call os%allocate (tgt2, m, 'tgt2')
  ! Interpolate in time, unless already done
  call interp%in_t (tgt, src, ls, us, src%mem, tgt0, only=only, &
                    extrap=tgt%allow_extrapolate, nbstep=nbstep)
  !-----------------------------------------------------------------------------
  ! The choice of interpolation order is for optimization only; the results
  ! would not differ if a fixed order was used.  The om(i)%p pointer points
  ! to the mesh that is being interpolated to
  !-----------------------------------------------------------------------------
  mem => tgt%mem(:,:,:,:,tgt%it,1)
  if (below(3) .or. above(3)) then
    call interp%in_z (tgt, lt, ut, src, ls, us, om, tgt0, tgt1,  .true., only=only)
    call interp%in_x (tgt, lt, ut, src, ls, us, om, tgt1, tgt2,  .true., only=only)
    if (tgt_int) then
      call interp%in_y (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, &
                        only=only, lock=tgt%mem_lock(tgt%it))
    else
      call interp%in_y (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, &
                        only=only)
    end if
  else if (below(2) .or. above(2)) then
    call interp%in_y (tgt, lt, ut, src, ls, us, om, tgt0, tgt1,  .true., only=only)
    call interp%in_x (tgt, lt, ut, src, ls, us, om, tgt1, tgt2,  .true., only=only)
    if (tgt_int) then
      call interp%in_z (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, &
                        only=only, lock=tgt%mem_lock(tgt%it))
    else
      call interp%in_z (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, &
                        only=only)
    end if
  else if (below(1) .or. above(1)) then
    call interp%in_x (tgt, lt, ut, src, ls, us, om, tgt0, tgt1,  .true., only=only)
    call interp%in_y (tgt, lt, ut, src, ls, us, om, tgt1, tgt2,  .true., only=only)
    call interp%in_z (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, only=only)
  else
    1 format(a,2(2x,3i3))
    call interp%in_z (tgt, lt, ut, src, ls, us, om, tgt0, tgt1,  .true., only=only)
    call interp%in_y (tgt, lt, ut, src, ls, us, om, tgt1, tgt2,  .true., only=only)
    if (tgt_int) then
      call interp%in_x (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, &
                        only=only, lock=tgt%mem_lock(tgt%it))
    else
      call interp%in_x (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, &
                        only=only)
    end if
  end if
  call os%deallocate (tgt0, 'tgt0')
  call os%deallocate (tgt1, 'tgt1')
  call os%deallocate (tgt2, 'tgt2')
  nullify    (tgt0, tgt1, tgt2)
  !-----------------------------------------------------------------------------
  ! When the source level is larger than the target level (restriction), record
  ! the higher level, to prevent later filling by same level sources.  When the
  ! source level is lower (prolongation), record the target level, since that
  ! is what is deposited
  !-----------------------------------------------------------------------------
  if (verbose > 1 .or. tgt%id == io%id_debug) then
   if (tgt%idx%d > 0) then
    write (stdout,'(a,2(i6,":",i2.2,3i4,1x),3f10.4,2(2x,3i3),1p,2g14.3)') &
      'differ: tgt,src,dist,l,u,max=', tgt%id, tgt%key, src%id, src%key, &
      tgt%distance_to(src)/tgt%size, lt, ut, &
      minval(tgt%mem(lt(1):ut(1),lt(2):ut(2),lt(3):ut(3), &
        tgt%idx%d,tgt%it,1)), &
      minval(tgt%mem(:,:,:,tgt%idx%d,tgt%it,1))
   end if
  end if
  if (associated(tgt%filled)) &
    tgt%filled (lt(1):ut(1),lt(2):ut(2),lt(3):ut(3)) = &
      10*max(src%level, tgt%level)
  call trace%end (itimer)
END SUBROUTINE differ

!===============================================================================
!> Interpolte btw patches of different size
!===============================================================================
SUBROUTINE tilted (self, tgt, src, all_tgt, all_src, only)
  class(guards_t)                      :: self
  class(patch_t)                       :: tgt, src
  logical, optional                    :: all_tgt, all_src
  integer, optional                    :: only
  !.............................................................................
  class(mesh_t), pointer, dimension(:) :: tm, sm
  type(meshp_t)                        :: om(3)
  real(8), dimension(3)                :: tl, tu, sl, su, fl, fu, pl, pu, sm_p
  logical, dimension(3)                :: below, above
  integer, dimension(3)                :: ls, us, lt, ut
  real, dimension(:,:,:,:), pointer    :: tgt0, tgt1, tgt2, mem
  !contiguous                           :: tgt0, tgt1, tgt2, mem
  logical                              :: src_all, tgt_int
  integer                              :: i, j, k, m(5), tid, sid
  real(8)                              :: src_p(3), src_v(3), src_e(3,3), e(3,3)
  real(8)                              :: s_corners(2,2,2,3), ss(2,3)
  integer, save                        :: itimer=0, nprint=10
  real(8), parameter                   :: eps=1e-6
  !-----------------------------------------------------------------------------
  ! Use all source cells if the all_src is present, or the rarget (child) is
  ! taking the first step.
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('guards_t%tilted', itimer=itimer)
  src_all = present(all_src) !.or. tgt%istep==0  ! use full source region
  tid = tgt%id
  sid = src%id
  !-----------------------------------------------------------------------------
  ! Produce values allso in internal cells if specifically ordered to, or if
  ! the source is at a higher level, or if a child is taking the first step.
  !-----------------------------------------------------------------------------
  tgt_int = present(all_tgt)                    ! produce full target
  tgt_int = tgt_int !.or. tgt%istep==0           ! first child update
  tgt_int = tgt_int .or. src%level > tgt%level  ! restrict
  !.............................................................................
  if (verbose > 2) then
    write (stdlog,'(f12.6,a,2(i6,":",i2.2,2x),4l1)') &
      wallclock(), pad_trace(1:indent)// &
      'guards_t%tilted: tgt,src,all_tgt,all_src,tgt_int,src_all=', &
      tid, tgt%level, sid, src%level, &
      present(all_tgt), present(all_src), tgt_int, src_all
  end if
  !.............................................................................
  tm => tgt%mesh                                ! target mesh
  sm => src%mesh                                ! source mesh
  do i=1,3
    om(i)%p => src%mesh(i)                      ! output mesh
  end do
  !-----------------------------------------------------------------------------
  ! Positions are in the target local coordinate system.  matmul(p,tgt%erot)
  ! tranforms from Cartesian to local, matmul(src%erot,p) from local to Cartesian
  !-----------------------------------------------------------------------------
  tl = - tm%s*0.5d0 - tm%ng*tm%d                ! tgt lower corner
  tu = + tm%s*0.5d0 + tm%ng*tm%d                ! tgt upper corner
  !-----------------------------------------------------------------------------
  ! Get src position, velocity, and orientation, in target coordinates
  !-----------------------------------------------------------------------------
  call tgt%get_position (src, src_p, src_v, src_e)
  !-----------------------------------------------------------------------------
  ! Get src corners in src local coordinates
  !-----------------------------------------------------------------------------
  if (src_all) then                             ! new child
    sl = - sm%s*0.5d0 - sm%ng*sm%d              ! src lower corner
    su = + sm%s*0.5d0 + sm%ng*sm%d              ! src upper corner
  else
    !---------------------------------------------------------------------------
    ! we add half the cell size margin as the face corners furthest away from the
    ! experiment centre in outer boundary patches can become no-mans-land. The
    ! shift is actually really really small, especially if tilt between patches
    ! is limited to max one cell, extrapolation is small and thus the error
    ! should be small.
    !--------------------------------------------------------------------------- 
    sl = - sm%s*0.5d0 - 0.5d0*sm%d                       ! src lower corner
    su = + sm%s*0.5d0 + 0.5d0*sm%d                       ! src upper corner
  end if
  !-----------------------------------------------------------------------------
  ! Rotate ALL src corners to tgt coordinates and get the extrema
  !-----------------------------------------------------------------------------
  ss(1,:)=sl
  ss(2,:)=su
  do k=1,2
    do j=1,2
      do i=1,2
        s_corners(i,j,k,:) = matmul (src_e, [ss(i,1),ss(j,2),ss(k,3)])
      end do
    end do
  end do
  do i=1,3
    sl(i) = minval(s_corners(:,:,:,i))
    su(i) = maxval(s_corners(:,:,:,i))
  end do
  if (verbose > 0) then
    write (stdlog,*) sl, 'sl src'
    write (stdlog,*) su, 'su src'
  end if
  !.............................................................................
  if (verbose > 0) then
    write (stdlog,*) sl, 'sl tgt'
    write (stdlog,*) su, 'su tgt'
  end if
  !-----------------------------------------------------------------------------
  ! Add relative position to corner values
  !-----------------------------------------------------------------------------
  sl = sl + src_p
  su = su + src_p
  !.............................................................................
  if (verbose > 1) then
    write (stdout,'(a,2(2(1x,3f6.2)))') &
      'guards_t%tilted: tl,tu,sl,su=', &
      tl/tm%d, tu/tm%d, sl/tm%d, su/tm%d
  end if
  !-----------------------------------------------------------------------------
  ! Use the whole target range?  Not that this should be redundant, since if we
  ! are to use the whole target range it must be included in the source range in
  ! any case
  !-----------------------------------------------------------------------------
  fl = max(tl,sl) - eps                         ! overlap lower corner
  fu = min(tu,su) + eps                         ! overlap upper corner
  lt =   floor(fl/tm%d + tm%o)                  ! integet tgt lower index
  ut = ceiling(fu/tm%d + tm%o)                  ! integer tgt upper index
  lt = max(lt,tm%lb)                            ! inside target halo
  ut = min(ut,tm%ub)                            ! inside target halo
  !.............................................................................
  if (verbose > 0) &
    write (stdout,'(a,2(i6,":",i2.2,3i3),2(2x,3i3),5(2x,3f10.6))') &
      'guards_t%tilted: tgt,key,src,key,lt,ut=', &
      tgt%id, tgt%key, src%id, src%key, lt, ut, &
      src_p, src_v, src_e(:,:), md(tgt%centre_nat(2)), md(src%centre_nat(2))
  !.............................................................................
  if (tgt%id == 1) then
    call validate%check (tgt, src%id, 'source ID')
    call validate%check (tgt, lt, 'lt')
    call validate%check (tgt, ut, 'ut')
    call validate%check (tgt, tgt%mem(:,:,:,:,tgt%it,1), 'before from_source')
  end if
  !.............................................................................
  if (present(only)) then
    call interp%from_only (tgt,lt,ut, src, only, filled=tgt%filled, &
                          all_src=all_src)
  else
    if (interp%use_staggered) then
      call interp%from_source_stag (tgt,lt,ut, src, filled=tgt%filled, all_src=all_src)
    else
      call interp%from_source (tgt,lt,ut, src, filled=tgt%filled, all_src=all_src)
    end if
  end if
  !if (associated(tgt%filled)) &
  !  tgt%filled (lt(1):ut(1),lt(2):ut(2),lt(3):ut(3)) = &
  !    10*max(src%level, tgt%level)
  !.............................................................................
  if (tgt%id == io%id_debug) &
    call validate%check (tgt, tgt%mem(:,:,:,:,tgt%it,1), 'after from_source')
  if (detailed_timer) &
    call trace%end (itimer)
!===============================================================================
contains
  FUNCTION md(r)
    real(8):: r, md
    md = modulo(r+math%pi,math%pi*2d0) - math%pi
  END FUNCTION md
END SUBROUTINE tilted

!===============================================================================
!> Info for one nbor
!===============================================================================
SUBROUTINE nbor_info (self, task, unit)
  class(guards_t):: self
  class(task_t):: task
  class(task_t), pointer:: nbtask
  integer, optional:: unit
  real(8):: distance(3)
  integer:: unit_L
  !-----------------------------------------------------------------------------
  unit_L = shared%option (unit, io_unit%mpi)
  nbtask => self%task
  distance = task%distance_to (nbtask)
  select type (nbtask)
  class is (patch_t)
  write(unit_L,'(3x,2a,i6,":",i2.2,i4,2x,3l1,2x,2l1,2(2x,3f9.3))') &
    'pa_t%nbor_info: kind,id,rank,level,needed,needs_me,dnload,BV,pos= ', &
    self%task%kind(1:6),self%task%id, self%task%level, self%task%rank, &
    self%needed, self%needs_me, self%download, self%task%is_set(bits%boundary), &
    self%task%is_set(bits%virtual), nbtask%centre_nat, &
    distance/(0.5d0*(task%size+self%task%size))
  end select
END SUBROUTINE nbor_info

!===============================================================================
!> Overloaded procedure
!===============================================================================
SUBROUTINE alloc_nbor (link, nbor)
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  allocate (guards_t:: nbor)
  call os%mem_register (storage_size(nbor), 1, 'nbor')
  call nbor%init ()
  nbor%link => link
  nbor%task => link%task
END SUBROUTINE alloc_nbor

END MODULE guards_mod
