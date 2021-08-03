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
  USE trace_mod
  USE omp_timer_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE boundary_mod
  USE interp_mod
  USE mesh_mod
  USE conserve_mod
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
SUBROUTINE init (self, target)
  class(guards_t)        :: self
  class(task_t), pointer :: target
  !-----------------------------------------------------------------------------
  ! Setup nbor distance and relations, in case these differ for patch_t
  !-----------------------------------------------------------------------------
  call trace%begin ('guards_t%init')
  call self%nbor_t%init (target)
  call self%pre_init
  select type (target)
  class is (patch_t)
  self%distance = self%task%distance (target)
  call target%nbor_relations(self%task, self%needed, self%needs_me, self%download)
  end select
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
                             only)
  class(guards_t)        :: self
  class(task_t)          :: source, target
  logical, optional      :: all_tgt, all_src
  integer, optional      :: n_op(3)
  integer, optional      :: only
  !.............................................................................
  character(len=8)       :: type
  integer                :: m(6), tid, sid, unit
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('guards_t%source_to_target', itimer=itimer)
  self%distance = source%distance (target)
  tid = target%id
  sid = source%id
  !.............................................................................
  if (verbose > 2 .or. target%id == io%id_debug &
                  .or. source%id == io%id_debug) then
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
      self%distance = source%distance (target)
      if (source%level == target%level) then
        type = 'same'
        call self%same (target, source, all_tgt, all_src, only)
        n_op  (1) = n_op(1) + 1
      else if (source%level == target%level-1) then
        type = 'prolong'
        call self%prolong (target, source, all_tgt, all_src, only)
        n_op(2) = n_op(2) + 1
      else if (source%level == target%level+1) then
        type = 'restrict'
        !-----------------------------------------------------------------------
        ! FIXME: For now, the conserve%restrict() can only handle interior cells
        !-----------------------------------------------------------------------
        if (source%parent%id == target%id .and. &
            .not. present(all_tgt) .and. .not. present(all_src)) then
          call conserve%restrict (target, source, all_tgt, all_src, only)
        else
          if (source%parent%id == target%id) then
            write (stderr,'(a,2(i6,":",i2.2,2x,2l1))') &
              'guards_t%source_to_target: warning, source,target,all=', &
              source%id, source%level, target%id, target%level, &
              present(all_tgt), present(all_src)
          end if
          call self%restrict (target, source, all_tgt, all_src, only)
        end if
        n_op(3) = n_op(3) + 1
      end if
      if (verbose > 0 .or. target%id == io%id_debug &
                      .or. source%id == io%id_debug) then
        write (unit,'(a8,2x,2(a,i6,":",i2.2,1x,4i3,3x),2x,3l1)') type, &
          'tgt:', tid, target%level, target%key, &
          'src:', sid, source%level, source%key, &
          present(all_tgt), present(all_src), present(only)
      end if
    end select
    !...........................................................................
    if (verbose > 1) then
      write (unit,'(a,2x,i6,";",i2.2,f12.6,2x,i6,":",i2.2,2x,3f8.3,2x,a)') &
        'guards_t%source_to_target: tgt,time,src,dist,type=', &
        target%id, target%level, target%time, source%id, source%level, &
        self%distance/target%size, type
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
SUBROUTINE same (self, tgt, src, all_tgt, all_src, only)
  class(guards_t)                      :: self
  class(patch_t)                       :: tgt, src
  logical, optional                    :: all_tgt, all_src
  integer, optional                    :: only
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
  if (detailed_timer) &
    call trace%begin ('guards_t%same', itimer=itimer)
  call self%info ('same', tgt, src, all_tgt, all_src, only)
  tm => tgt%mesh                             ! tgt mesh
  sm => src%mesh                             ! src mesh
  do i=1,3
    om(i)%p => sm(i)
  end do
  d = src%distance (tgt)
  sm_p    = tm%p + d                         ! src position
  tl = tm%p - tm%s*0.5d0 - tm%ng*tm%d        ! tgt lower corner
  tu = tm%p + tm%s*0.5d0 + tm%ng*tm%d        ! tgt upper corner
  sl = sm_p - sm%s*0.5d0                     ! src lower corner
  su = sm_p + sm%s*0.5d0                     ! src upper corner
  fl = max(tl,sl)                            ! overlap lower corner
  fu = min(tu,su)                            ! overlap upper corner
  ls = src%index_only (fl, roundup=.true.)   ! src lower index
  us = src%index_only (fu, rounddn=.true.)   ! src upper index
  lt = tgt%index_only (fl, roundup=.true.)   ! tgt lowwe index
  ut = tgt%index_only (fu, rounddn=.true.)   ! tgt upper index
  if (verbose > 2) then
    write (stdlog,'(a,2i6,4(2x,3i4))') 'guards_t%same: tgt,src,ls,us,lt,ut=', &
      tgt%id, src%id, ls, us, lt, ut
    call tgt%stats ('guards_t%same tgt')
    call src%stats ('guards_t%same src')
  end if
  !-----------------------------------------------------------------------------
  ! Time interpolate in the target region
  !-----------------------------------------------------------------------------
  m = shape(tgt%mem)
  mem => tgt%mem(:,:,:,:,tgt%it,1)
  call interp%in_t   (tgt, src, ls, us, src%mem, mem, lt, ut, only=only, &
                        extrap=tgt%allow_extrapolate)
  if (associated(tgt%filled)) &
    tgt%filled (lt(1):ut(1),lt(2):ut(2),lt(3):ut(3)) = src%level
  !.............................................................................
  if (detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE same

!===============================================================================
SUBROUTINE prolong (self, tgt, src, all_tgt, all_src, only)
  class(guards_t)   :: self
  class(patch_t)    :: tgt, src
  logical, optional :: all_tgt, all_src
  integer, save     :: itimer=0
  integer, optional :: only
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('guards_t%prolong', itimer=itimer)
  call self%info ('prolong', tgt, src, all_tgt, all_src, only)
  call self%differ (tgt, src, all_tgt, all_src, only)
  if (detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE prolong

!===============================================================================
SUBROUTINE restrict (self, tgt, src, all_tgt, all_src, only)
  class(guards_t)   :: self
  class(patch_t)    :: tgt, src
  logical, optional :: all_tgt, all_src
  integer, save     :: itimer=0
  integer, optional :: only
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('guards_t%restrict',  itimer=itimer)
  call self%info ('restrict', tgt, src, all_tgt, all_src, only)
  call self%differ (tgt, src, all_tgt, all_src, only)
  if (detailed_timer) &
    call trace%end (itimer)
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
  if (verbose > 0) then
    if (present(only)) then
      o = only
    else
      o = 0
    end if
    write (stdlog,'(2(f12.6,a,2x,i6,":",i2.2),2(2x,3f7.3),2x,2l1,i3)') &
      wallclock(), pad_trace(1:indent)//'guards_t%'// &
      trim(label)//': tgt,time,src,dist,all_tgt,all_src=', &
      target%id, target%level, &
      target%time, '  src:', source%id, source%level, &
      self%distance/target%size, source%distance(target)/target%size, &
      present(all_tgt), present(all_src), o
  end if
  !.............................................................................
END SUBROUTINE info

!===============================================================================
!> Interpolate in time, leaving the result in the target
!===============================================================================
SUBROUTINE differ (self, tgt, src, all_tgt, all_src, only)
  class(guards_t)                      :: self
  class(patch_t)                       :: tgt, src
  logical, optional                    :: all_tgt, all_src
  integer, optional                    :: only
  !.............................................................................
  class(task_t), pointer               :: nbtask
  class(mesh_t), pointer, dimension(:) :: tm, sm
  type(meshp_t)                        :: om(3)
  real(8), dimension(3)                :: tl, tu, sl, su, fl, fu, pl, pu, sm_p
  logical, dimension(3)                :: below, above
  integer, dimension(3)                :: ls, us, lt, ut
  real, dimension(:,:,:,:), pointer    :: tgt0, tgt1, tgt2, mem
  !contiguous                           :: tgt0, tgt1, tgt2, mem
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
  unit = merge (stdout, stdlog, verbose > 2 .or. tgt%id==io%id_debug)
  src_all = present(all_src) !.or. tgt%istep==0  ! use full source region
 !src_all = src_all .or. tgt%istep==0            ! first child update
  !-----------------------------------------------------------------------------
  ! Produce values also in internal cells if specifically ordered to, or if
  ! the source is at a higher level, or if a child is taking the first step.
  !-----------------------------------------------------------------------------
  tgt_int = present(all_tgt)                    ! produce full target
  tgt_int = tgt_int .or. tgt%istep==0           ! first child update
  tgt_int = tgt_int .or. src%level > tgt%level  ! restrict
  tgt_int = tgt_int .and. tgt%time > 0.0        ! FIXME: restarts!
  !.............................................................................
  if (verbose > 1) then
    write (unit,'(f12.6,a,2(i6,":",i2.2,2x),4l1)') &
      wallclock(), pad_trace(1:indent)// &
      'guards_t%differ: src,tgt,all_tgt,all_src,tgt_int,src_all=', &
      src%id, src%level, tgt%id, tgt%level, &
      present(all_tgt), present(all_src), tgt_int, src_all
  end if
  if (unit == stdout) then
    write (unit,'(f12.6,2x,a,2(i6,":",i2.2,2x),4l1)') &
      wallclock(), &
      'guards_t%differ: src,tgt,all_tgt,all_src,tgt_int,src_all=', &
      src%id, src%level, tgt%id, tgt%level, &
      present(all_tgt), present(all_src), tgt_int, src_all
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
  tl = tm%p - tm%s*0.5d0 - tm%ng*tm%d           ! tgt lower corner
  tu = tm%p + tm%s*0.5d0 + tm%ng*tm%d           ! tgt upper corner
  sm_p = tm%p + tgt%distance_to (src)           ! position relative target
  if (src_all) then                             ! new child
    sl = sm_p - sm%s*0.5d0 - sm%ng*sm%d         ! src lower corner
    su = sm_p + sm%s*0.5d0 + sm%ng*sm%d         ! src upper corner
  else
    sl = sm_p - sm%s*0.5d0                      ! src lower corner
    su = sm_p + sm%s*0.5d0                      ! src upper corner
  end if
  if (verbose > 2 .or. unit==stdout) then
    write (unit,'(a,i6,2x,3i5,5(2x,3f8.3))') &
      'guards_t%differ src mesh: id,s,p,b=', &
      sm(1)%id, sm%n, sm%s, sm%p, sm%b
    write (unit,'(a,i6,2x,3i5,5(2x,3f8.3))') &
      'guards_t%differ tgt mesh: id,s,p,b=', &
      tm(1)%id, tm%n, tm%s, tm%p, tm%b
  end if
  !-----------------------------------------------------------------------------
  ! Use the whole target range?  Not that this should be redundant, since if we
  ! are to use the whole target range it must be included in the source range in
  ! any case
  !-----------------------------------------------------------------------------
  fl = max(tl,sl)                               ! overlap lower corner
  fu = min(tu,su)                               ! overlap upper corner
  if (tgt%level < src%level) then
    !-------------------------------------------------------------------------
    ! Restrict operation:  Use indices that are strictly inside intersection
    !-------------------------------------------------------------------------
    ls = src%index_only (fl, roundup=.true.)    ! src lower index
    us = src%index_only (fu, rounddn=.true.)    ! src upper index
    tl = src%myposition (ls)                    ! position of lt
    tu = src%myposition (us)                    ! position of ut
    lt = tgt%index_only (tl, roundup=.true.)    ! target lower corner
    ut = tgt%index_only (tu, rounddn=.true.)    ! target upper corner
  else
    !-------------------------------------------------------------------------
    ! Prolongation operation: Use tgt-indices that are strictly inside, and
    ! the correspondingly necessary source indices (which may fall outside
    ! the common intersection).
    !-------------------------------------------------------------------------
    lt = tgt%index_only (fl, roundup=.true.)    ! target lower corner
    ut = tgt%index_only (fu, rounddn=.true.)    ! target upper corner
    tl = tgt%myposition (lt)                    ! position of lt
    tu = tgt%myposition (ut)                    ! position of ut
    ls = src%index_only (fl, rounddn=.true.)    ! src lower index
    us = src%index_only (fu, roundup=.true.)    ! src upper index
  end if
  if (src_all) then
    ls = max(ls,sm%lb)                          ! inside source halo
    us = min(us,sm%ub)                          ! inside source halo
  else
    ls = max(ls,sm%li)                          ! inside source interior
    us = min(us,sm%ui)                          ! inside source interior
  end if
  if (any(lt < tm%lb .or. ut > tm%ub)) then
    call io%warning ('target limits outside range')
  end if
  lt = max(lt,tm%lb)                            ! inside target halo
  ut = min(ut,tm%ub)                            ! inside target halo
  !-----------------------------------------------------------------------------
  ! Decide relative positioning, in terms of axes
  !-----------------------------------------------------------------------------
  below = src%mesh%p < tgt%mesh%p - tgt%mesh%s/2d0
  above = src%mesh%p > tgt%mesh%p + tgt%mesh%s/2d0
  !.............................................................................
  if (verbose > 1) then
    write (stdlog,'(f12.6,a,2(i6,":",i2.2,1x),2(1x,2(1x,3i4)),3(2x,3l1))') &
      wallclock(), pad_trace(1:indent)// &
      'guards_t%differ: src,tgt,ls,us,lt,ut,bel,abv=', &
      src%id, src%level, tgt%id, tgt%level, ls, us, lt, ut,  &
      below, above, tgt_int, src_all
  end if
  if (unit == stdout) then
    write (unit,'(f12.6,2x,a,2(i6,":",i2.2,1x),2(1x,2(1x,3i4)),3(2x,3l1))') &
      wallclock(), &
      'guards_t%differ: src,tgt,ls,us,lt,ut,bel,abv=', &
      src%id, src%level, tgt%id, tgt%level, ls, us, lt, ut,  &
      below, above, tgt_int, src_all
  end if
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays with common size
  !-----------------------------------------------------------------------------
  nullify (tgt0, tgt1)
  m = shape(tgt%mem(:,:,:,:,:,1))
  allocate (tgt0(m(1),m(2),m(3),m(4)))
  allocate (tgt1(m(1),m(2),m(3),m(4)))
  allocate (tgt2(m(1),m(2),m(3),m(4)))
  ! Interpolate in time, unless already done
  call interp%in_t (tgt, src, ls, us, src%mem, tgt0, only=only, &
                    extrap=tgt%allow_extrapolate)
  !-----------------------------------------------------------------------------
  ! The choice of interpolation order is for optimization only; the results
  ! would not differ if a fixed order was used.  The om(i)%p pointer points
  ! to the mesh that is being interpolated to
  !-----------------------------------------------------------------------------
  mem => tgt%mem(:,:,:,:,tgt%it,1)
  if (below(3) .or. above(3)) then
    call interp%in_z (tgt, lt, ut, src, ls, us, om, tgt0, tgt1,  .true., only=only)
    call interp%in_x (tgt, lt, ut, src, ls, us, om, tgt1, tgt2,  .true., only=only)
    call interp%in_y (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, only=only)
  else if (below(2) .or. above(2)) then
    call interp%in_y (tgt, lt, ut, src, ls, us, om, tgt0, tgt1,  .true., only=only)
    call interp%in_x (tgt, lt, ut, src, ls, us, om, tgt1, tgt2,  .true., only=only)
    call interp%in_z (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, only=only)
  else if (below(1) .or. above(1)) then
    call interp%in_x (tgt, lt, ut, src, ls, us, om, tgt0, tgt1,  .true., only=only)
    call interp%in_y (tgt, lt, ut, src, ls, us, om, tgt1, tgt2,  .true., only=only)
    call interp%in_z (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, only=only)
  else
    1 format(a,2(2x,3i3))
    call interp%in_z (tgt, lt, ut, src, ls, us, om, tgt0, tgt1,  .true., only=only)
    call interp%in_y (tgt, lt, ut, src, ls, us, om, tgt1, tgt2,  .true., only=only)
    call interp%in_x (tgt, lt, ut, src, ls, us, om, tgt2, mem , tgt_int, only=only)
  end if
  deallocate (tgt0, tgt1, tgt2)
  nullify    (tgt0, tgt1, tgt2)
  !-----------------------------------------------------------------------------
  ! When the source level is larger than the target level (restriction), record
  ! the higher level, to prevent later filling by same level sources.  When the
  ! source level is lower (prolongation), record the target level, since that
  ! is what is deposited
  !-----------------------------------------------------------------------------
  if (associated(tgt%filled)) &
    tgt%filled (lt(1):ut(1),lt(2):ut(2),lt(3):ut(3)) = max(src%level, tgt%level)
  call trace%end (itimer)
END SUBROUTINE differ

END MODULE guards_mod
