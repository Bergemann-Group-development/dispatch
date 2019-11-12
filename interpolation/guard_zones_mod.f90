!===============================================================================
!> Optimized restrict operations, intended for guard zones (not conservative)
!>
!>    1   2   3   4   5   6   7   8   9   0  11    17  18  19  20  21  22
!>    +---+---+-|-o---o---o---o---o---o---0---0-| --0---o---o---o-|-+---+---+
!>  -o-o-o-o-o-o|+-+-+                          |            +-+-+|o-o-o-o-o-o
!>   4 5 6 7 8 9 0 1 2                                       1 2 3 4 5 6 7 8 9
!>         o-o-o|o-o-o-o-o-o-o-o-o-o-o-o-o-o-o-o|o-o-o
!>         1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2
!>
!> When boundaries are lined up for 16x16x16 patches, restrict operations should
!> involve interpolations over index 14-15, 16-17, and 18-19 on the high source
!> side, and interpolations over index 4-5, 6-7, and 8-9 on the low source side,
!> while of coure 1-3 and 20-22 on the low and high target side.  When the
!> source overlaps with the target it covers target index 4-11 left and 12-19
!> right, with its own range 4-19.
!===============================================================================
MODULE guard_zones_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE patch_mod
  USE mesh_mod
  USE omp_timer_mod
  implicit none
  private
  type, public:: guard_zones_t
    procedure(restrict_simple), pointer:: method => null()
  contains
    procedure:: init
    procedure:: restrict_simple
    procedure:: different_opt
    procedure:: restrict_vect
    procedure, nopass:: info
  end type
  real(8):: update_t(3)=0d0
  integer:: verbose=0
  logical:: detailed_timer=.false.
  type(guard_zones_t), public:: guard_zones
CONTAINS

!===============================================================================
!> Choose method and verbosity
!===============================================================================
SUBROUTINE init (self)
  class(guard_zones_t):: self
  character(len=6):: method='opt'
  namelist /guard_zone_params/ method, detailed_timer, verbose
  integer:: iostat
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_zones_t%init')
  if (first_time) then
    !$omp critical (init_cr)
    if (first_time) then
      rewind (io%input)
      read (io%input, guard_zone_params, iostat=iostat)
      write (io%output, guard_zone_params)
      first_time = .false.
    end if
    !$omp end critical (init_cr)
  end if
  select case (trim(method))
  case ('simple')
    self%method => restrict_simple
  case ('opt')
    self%method => different_opt
  case ('vect')
    self%method => restrict_vect
  case default
    call io%abort ('guard zone method unknown')
  end select
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Compute interpolated guard zone values by setting up interface to the
!> slow but correct patch_t%interpolate() point-by-point interpolation.
!===============================================================================
SUBROUTINE restrict_simple (self, target, source, jt, pt, same)
  class(guard_zones_t):: self
  class(patch_t), pointer:: target, source
  integer:: jt(2)
  real(4):: pt(2)
  logical, optional:: same
  !.............................................................................
  real, dimension(:,:,:), pointer:: tmem, smem, td, sd
  class(mesh_t), pointer:: mesh
  integer:: i1, i2, i3, iv, jv, tl(3), tu(3), sl(3), su(3)
  real(8):: tlc(3), tuc(3), slc(3), suc(3), dist(3)
  real(4):: pl(3), pu(3)
  real(8):: used
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Interpolate in time and space, point-by-point
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_zones_t%restrict_simple', itimer=itimer)
  if (detailed_timer) &
    used = wallclock()
  call source%time_interval (target%time, jt, pt)
  do iv=1,target%nv
    jv = 0
    !---------------------------------------------------------------------------
    ! Target lower and upper corner in real space
    !---------------------------------------------------------------------------
    tlc = target%myposition (target%mesh%lb, jv)
    tuc = target%myposition (target%mesh%ub, jv)
    !---------------------------------------------------------------------------
    ! Map to source index space and limit to source interior
    !---------------------------------------------------------------------------
    sl = source%index_only (tlc, jv, nowrap=.true.)
    su = source%index_only (tuc, jv, nowrap=.true.)
    sl = min(max(sl,source%mesh%li),source%mesh%ui-1)
    su = min(max(su,source%mesh%li),source%mesh%ui-1)
    !---------------------------------------------------------------------------
    ! Map source points to target legal index space
    !---------------------------------------------------------------------------
    slc = source%myposition(sl, jv)
    suc = source%myposition(su, jv)
    tl = target%index_only (slc, jv, nowrap=.true., roundup=.true.)
    tu = target%index_only (suc, jv, nowrap=.true., roundup=.true.)
    tl = min(max(tl,target%mesh%lb),target%mesh%ub)
    tu = min(max(tu,target%mesh%lb),target%mesh%ub)
    if (verbose > 0) &
      write(io_unit%log,'(a,i3,2(2x,3i4),2(2x,3f7.3))') &
        'guard_zones_t%restrict_simple1: iv, tl, tu =', iv, tl, tu
    !call dnload_range (target, source, iv, tl, tu)
    !if (verbose > 0) &
    !  write(io_unit%log,'(a,i3,2(2x,3i4),2(2x,3f7.3))') &
    !    'guard_zones_t%restrict_simple2: iv, tl, tu =', iv, tl, tu
    do i3=tl(3),tu(3); do i2=tl(2),tu(2); do i1=tl(1),tu(1)
      call interpolate (source, target, [i1,i2,i3], iv, jt, pt)
    end do; end do; end do
  end do
  if (detailed_timer) then
    used = wallclock()-used
    write (stdout,'(a,6p,f6.3," mus/cell")') &
      "guard_zone_t%restrict_simple:", used/product(tu-tl+1)
  end if
  call trace%end (itimer)
END SUBROUTINE restrict_simple

!===============================================================================
!> Compute interpolated guard zone values by setting up interface to the
!> slow but correct patch_t%interpolate() point-by-point interpolation.
!===============================================================================
SUBROUTINE different_opt (self, target, source, jt, pt, same)
  class(guard_zones_t):: self
  class(patch_t), pointer:: target, source
  integer:: jt(2)
  real(4):: pt(2)
  logical, optional:: same
  !.............................................................................
  real, parameter:: eps=1e-4
  real, dimension(:,:,:), pointer:: tmem, smem, td, sd
  class(mesh_t), pointer:: mesh
  integer:: i, i1, i2, i3, iv, jv, tl(3), tu(3), sl(3), su(3), cells
  real:: pi(16)
  real(8):: tlc, tuc, slc, suc, dist, box, tc
  class(mesh_t), pointer:: tm, sm
  real(8):: start
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Interpolate in time and space, point-by-point
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_zones_t%different_opt', itimer=itimer)
  if (detailed_timer) then
    start = wallclock()
  end if
  call source%time_interval (target%time, jt, pt)
  do iv=1,target%nv
    jv = 0
    do i=1,3
      !-------------------------------------------------------------------------
      ! Source distance from target center
      !-------------------------------------------------------------------------
      tm => target%mesh(i)
      sm => source%mesh(i)
      box = tm%b
      dist = sm%p - tm%p
      if (target%periodic(i)) &
        dist = modulo(dist+0.5d0*box,box) - 0.5d0*box
      !-------------------------------------------------------------------------
      ! Target corner distances from center
      !-------------------------------------------------------------------------
      tlc = tm%r(tm%lb)
      tuc = tm%r(tm%ub)
      !-------------------------------------------------------------------------
      ! Source lower and upper li and ui points, relative to target center
      !-------------------------------------------------------------------------
      slc = sm%r(sm%li) + dist
      suc = sm%r(sm%ui) + dist
      !-------------------------------------------------------------------------
      ! Lower and upper index distance from target center
      !-------------------------------------------------------------------------
      tlc = max(tlc,slc)
      tuc = min(tuc,suc)
      !-------------------------------------------------------------------------
      ! Target indices to lower and upper guard zone limits
      !-------------------------------------------------------------------------
      tc = 1d0/tm%d
      tl(i) = ceiling(tlc*tc + tm%o - eps)
      tu(i) = ceiling(tuc*tc + tm%o - eps)
      if (verbose > 0) &
        write(io_unit%log,'(a,i3,2x,5f8.3,2x,2i4)') &
        'guard_zones_t%different_opt: i,dist,tlc,tuc,slc,suc,tl,tu =', &
        i, dist, tlc, tuc, slc, suc, tl(i), tu(i)
    end do
    !call dnload_range (target, source, iv, tl, tu)
    !if (verbose > 0) &
    !  write(io_unit%log,'(a,i3,2(2x,3i4),2(2x,3f7.3))') &
    !    'guard_zones_t%restrict_simple2: iv, tl, tu =', iv, tl, tu
    if (present(same)) then
      call interpolate (source, target, [i1,i2,i3], iv, jt, pt, pi=pi, set=.true.)
      do i3=tl(3),tu(3); do i2=tl(2),tu(2); do i1=tl(1),tu(1)
        call interpolate (source, target, [i1,i2,i3], iv, jt, pt, pi=pi)
      end do; end do; end do
    else
      do i3=tl(3),tu(3); do i2=tl(2),tu(2); do i1=tl(1),tu(1)
        call interpolate (source, target, [i1,i2,i3], iv, jt, pt)
      end do; end do; end do
    end if
  end do
  if (detailed_timer) then
    cells = product(tu-tl+1)
    write (stdout,'(a,6p,2(2x,i6,i3),i6," cells,",f9.3," mus/cell")') &
      "guard_zones_t%different_opt: target, level, source, level,", &
      target%id, target%level, source%id, source%level, cells, &
      (wallclock()-start)/max(1,cells)
  end if
  call trace%end (itimer)
END SUBROUTINE different_opt

!===============================================================================
!> Compute interpolated guard zone values
!===============================================================================
SUBROUTINE restrict_vect (self, target, source, jt, pt, same)
  class(guard_zones_t):: self
  class(patch_t), pointer:: target, source
  integer:: jt(2)
  real(4):: pt(2)
  logical, optional:: same
  !.............................................................................
  real, dimension(:,:,:), pointer:: tmem, smem, td, sd
  class(mesh_t), pointer:: mesh
  integer:: i, i1, i2, i3, j, k, iv, tl(3), tu(3), sl(3), su(3), j0(3)
  real(8):: tlc(3), tuc(3), slc(3), suc(3)
  real(4):: p(3), p1, p2, p3, w(8)
  real(8):: used
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Map the source corners into target index space
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_zones_t%restrict_vect', itimer=itimer)
  if (detailed_timer) &
    used = wallclock()
  !-----------------------------------------------------------------------------
  ! Interpolate in time and space
  !-----------------------------------------------------------------------------
  call source%time_interval (target%time, jt, pt)
  if (verbose > 0) &
    write(io_unit%log,'(a,2(2x,i6,i3,2x,3f7.3),2x,3f7.3)') &
      'target, level, lc, source, level, lc, relative =', &
      target%id, target%level, target%llc_cart, &
      source%id, source%level, source%llc_cart, &
      source%distance(target)/(0.5d0*(target%size+source%size))
  do iv=1,target%nv
    !---------------------------------------------------------------------------
    ! Target lower and upper corner in real space; these are the positions of
    ! the first and last target points.
    !---------------------------------------------------------------------------
    tlc = target%myposition (target%mesh%lb, iv)
    tuc = target%myposition (target%mesh%ub, iv)
    !---------------------------------------------------------------------------
    ! Map to source index space and limit to source interior. These are the last
    ! source points to the left of a target point, with the next source point to
    ! the right of the target point, and still inside the source interior.
    !---------------------------------------------------------------------------
    sl = source%index_only (tlc, 1, nowrap=.true.)
    su = source%index_only (tuc, 1, nowrap=.true.)
    sl = min(max(sl,source%mesh%li),source%mesh%ui-1)
    su = min(max(su,source%mesh%li),source%mesh%ui-1)
    !---------------------------------------------------------------------------
    ! Map those source points back to the target index space, rounding up to
    ! find the target points between the source points.
    !---------------------------------------------------------------------------
    slc = source%myposition(sl, iv)
    suc = source%myposition(su, iv)
    !---------------------------------------------------------------------------
    ! Compensate for periodic wrapping of corner points
    !---------------------------------------------------------------------------
    !where (target%position < source%position .and. slc > source%position)
      !slc = slc - target%mesh%b
    !end where
    !where (target%position > source%position .and. suc < source%position)
      !suc = suc + target%mesh%b
    !end where
    tl = target%index_only (slc, 1, roundup=.true., nowrap=.true.)
    tu = target%index_only (suc, 1, roundup=.true., nowrap=.true.)
    !---------------------------------------------------------------------------
    ! The brackets below should not really be needed, since the points should
    ! all be inside the target legal range
    !---------------------------------------------------------------------------
    tl = min(max(tl,target%mesh%lb),target%mesh%ub)
    tu = min(max(tu,target%mesh%lb),target%mesh%ub)
    if (verbose > 0) then
      write(io_unit%log,1) &
        'guard_zones_t%restrict_vect: iv, tlc, sl, slc, tl =', &
        iv, tlc, sl, slc, tl
      1 format(a,i3,2(2(2x,3f8.4,2x,3i4)))
      write(io_unit%log,1) &
        'guard_zones_t%restrict_vect: iv, tuc, su, suc, tu =', &
        iv, tuc, su, suc, tu
    end if
    !---------------------------------------------------------------------------
    ! Compute indices and interpolation weights for the first target point, and
    ! save the weights for the 8 nearest neighbor points in source space
    !---------------------------------------------------------------------------
    call source%index_space (target%myposition(tl), 1, j0, p)
    j = 0
    do i3=0,1
      p3 = merge (p(3), 1.0-p(3), i3==1)
      do i2=0,1
        p2 = merge (p(2), 1.0-p(2), i2==1)
        do i1=0,1
          p1 = merge (p(1), 1.0-p(1), i1==1)
          j = j+1
          w(j) = p1*p2*p3
        end do
      end do
    end do
    !---------------------------------------------------------------------------
    if (verbose > 0) &
      write(io_unit%log,'(a,i3,4(2x,3i4))') &
        'guard_zones_t%restrict_vect: iv, tl, tu, j0 =', &
        iv, tl, tu, j0
    !---------------------------------------------------------------------------
    ! Map target and source memory to contiguous 3-D arrays
    !---------------------------------------------------------------------------
    tmem => target%mem(:,:,:,iv,target%it,1)
    td => target%mem(:,:,:,target%idx%d,target%it,1)
    do k=1,2
      smem => source%mem(:,:,:,iv,jt(k),1)
      sd => source%mem(:,:,:,source%idx%d,jt(k),1)
      call update (target, tmem, td, tl, tu, iv, source, smem, sd, j0, k, &
                   pt(k), w, verbose)
    end do
  end do
  if (detailed_timer) then
    used = wallclock()-used
    write (stdout,'(a,6p,f6.3," mus/cell")') &
      "guard_zone_t%restrict_vect:", used/product(tu-tl+1)
  end if
  call trace%end (itimer)
END SUBROUTINE restrict_vect

!===============================================================================
!> Update interpolated guard zone values by looping over target indices,
!> incrementing source indices by two (FIXME: this assumes factor two AMR)
!> for each step in target indices.
!>
!> For a 16x16x16 test case, the expected source index ranges are:
!>  case  source   target   j0    l   u
!>  -----------------------------------
!>  left:   4-9      1-3     4    1   3
!> right:  14-19    19-22   14   19  22
!>  full:   4-19     4-11    4    4  11
!>  full:   4-19    12-19    4   12  19
!===============================================================================
SUBROUTINE update (target, tmem, td, l, u, iv, source, smem, sd, j0, k, &
                   pt, w, verbose)
  class(patch_t), pointer:: target, source
  real, dimension(:,:,:), contiguous, intent(out):: tmem
  real, dimension(:,:,:), contiguous, intent(in):: smem
  real, dimension(:,:,:), contiguous, intent(inout):: td, sd
  integer:: l(3), u(3), iv, j0(3), k, verbose
  real:: pt, w(8)
  !.............................................................................
  real(8):: wc(4)
  real, allocatable:: f(:,:), tot(:)
  integer:: i, j, n, i1 , i2, i3, j1, j2, j3, jj(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_zones_t%update')
  if (detailed_timer) then
    wc(1) = wallclock()
  end if
  n = product(u-l+1)
  allocate (f(n,8), tot(n))
  !-----------------------------------------------------------------------------
  ! Gather
  !-----------------------------------------------------------------------------
  i = 0
  j3 = j0(3)
  do i3=l(3),u(3)
    j2 = j0(2)
    do i2=l(2),u(2)
      j1 = j0(1)
      if (target%pervolume(iv)) then
        do i1=l(1),u(1)
          i = i+1
!define VALIDATE
#ifdef VALIDATE
          jj = [j1,j2,j3]
          if (any(jj < source%li .or. jj > source%ui)) then
            write(stderr     ,1) iv, jj, j0, l, u
            write(io_unit%log,1) iv, jj, j0, l, u
            1 format("guard_zones_t%update: ERROR iv,jj,j0,l,u =",i3,4(2x,3i4))
          end if
#endif
          f(i,1) = smem(j1  ,j2  ,j3  )/sd(j1  ,j2  ,j3  )
          f(i,2) = smem(j1+1,j2  ,j3  )/sd(j1+1,j2  ,j3  )
          f(i,3) = smem(j1  ,j2+1,j3  )/sd(j1  ,j2+1,j3  )
          f(i,4) = smem(j1+1,j2+1,j3  )/sd(j1+1,j2+1,j3  )
          f(i,5) = smem(j1  ,j2  ,j3+1)/sd(j1  ,j2  ,j3+1)
          f(i,6) = smem(j1+1,j2  ,j3+1)/sd(j1+1,j2  ,j3+1)
          f(i,7) = smem(j1  ,j2+1,j3+1)/sd(j1  ,j2+1,j3+1)
          f(i,8) = smem(j1+1,j2+1,j3+1)/sd(j1+1,j2+1,j3+1)
          j1 = j1+2
        end do
      else
        do i1=l(1),u(1)
          i = i+1
#ifdef VALIDATE
          jj = [j1,j2,j3]
          if (any(jj < source%li .or. jj > source%ui)) then
            write(stderr     ,1) iv, jj, j0, l, u
            write(io_unit%log,1) iv, jj, j0, l, u
          end if
#endif
          f(i,1) = smem(j1  ,j2  ,j3  )
          f(i,2) = smem(j1+1,j2  ,j3  )
          f(i,3) = smem(j1  ,j2+1,j3  )
          f(i,4) = smem(j1+1,j2+1,j3  )
          f(i,5) = smem(j1  ,j2  ,j3+1)
          f(i,6) = smem(j1+1,j2  ,j3+1)
          f(i,7) = smem(j1  ,j2+1,j3+1)
          f(i,8) = smem(j1+1,j2+1,j3+1)
          j1 = j1+2
        end do
      end if
      j2 = j2+2
    end do
    j3 = j3+2
  end do
  if (detailed_timer) then
    wc(2) = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Compute weighted results, vectorized
  !-----------------------------------------------------------------------------
  tot = 0.0
  if (target%unsigned(iv)) then
    do j=1,8
    do i=1,n
      tot(i) = tot(i) + pt*w(j)*log(f(i,j))
    end do
    end do
    tot = exp(tot)
  else
    do j=1,8
    do i=1,n
      tot(i) = tot(i) + pt*w(j)*f(i,j)
    end do
    end do
  end if
  if (detailed_timer) then
    wc(3) = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Scatter, vectorized
  !-----------------------------------------------------------------------------
  i = 0
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      if (k==1) then
        do i1=l(1),u(1)
          i = i+1
          tmem(i1,i2,i3) = tot(i)
        end do
      else
        do i1=l(1),u(1)
          i = i+1
          tmem(i1,i2,i3) = tmem(i1,i2,i3) + tot(i)
        end do
      end if
    end do
  end do
  !-----------------------------------------------------------------------------
  ! For per-volume quantities, multiply with mass density
  !-----------------------------------------------------------------------------
  if (target%pervolume(iv) .and. k==2) then
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        do i1=l(1),u(1)
          tmem(i1,i2,i3) = tmem(i1,i2,i3)*td(i1,i2,i3)
        end do
      end do
    end do
  end if
  deallocate (f, tot)
  if (detailed_timer) then
    wc(4) = wallclock()
    do i=1,3
      !$omp atomic update
      update_t(i) = update_t(i) + (wc(i+1)-wc(i))
    end do
  end if
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Interpolate in space and time, for a case where it is safe = the source point
!> is inside the valid domain (inside a distance of half a cell from the source
!> patch boundary for non-staggered points)
!===============================================================================
SUBROUTINE interpolate (source, target, ii, iv, jt, pt, pi, set)
  class(patch_t)           :: source
  class(patch_t), pointer  :: target
  integer                  :: ii(3), iv, jt(2)
  real                     :: pt(2)
  real, optional           :: pi(:)
  logical, optional        :: set
  !.............................................................................
  real(8)                  :: pos(3)
  integer                  :: j(3)
  real                     :: p(3)
  integer, save            :: itimer=0
  !-----------------------------------------------------------------------------
  if (present(set)) then
    pi(01) = pt(1)*(1.-p(3))*(1.-p(2))*(1.0-p(1))
    pi(02) = pt(1)*(1.-p(3))*(1.-p(2))*(    p(1))
    pi(03) = pt(1)*(1.-p(3))*(   p(2))*(1.0-p(1))
    pi(04) = pt(1)*(1.-p(3))*(   p(2))*(    p(1))
    pi(05) = pt(1)*(   p(3))*(1.-p(2))*(1.0-p(1))
    pi(06) = pt(1)*(   p(3))*(1.-p(2))*(    p(1))
    pi(07) = pt(1)*(   p(3))*(   p(2))*(1.0-p(1))
    pi(08) = pt(1)*(   p(3))*(   p(2))*(    p(1))
    pi(09) = pt(2)*(1.-p(3))*(1.-p(2))*(1.0-p(1))
    pi(10) = pt(2)*(1.-p(3))*(1.-p(2))*(    p(1))
    pi(11) = pt(2)*(1.-p(3))*(   p(2))*(1.0-p(1))
    pi(12) = pt(2)*(1.-p(3))*(   p(2))*(    p(1))
    pi(13) = pt(2)*(   p(3))*(1.-p(2))*(1.0-p(1))
    pi(14) = pt(2)*(   p(3))*(1.-p(2))*(    p(1))
    pi(15) = pt(2)*(   p(3))*(   p(2))*(1.0-p(1))
    pi(16) = pt(2)*(   p(3))*(   p(2))*(    p(1))
  else if (source%pervolume(iv)) then
    call interpolate_specific (source, target, ii, iv, jt, pt)
  else
    pos = target%myposition (ii, iv)
    call source%index_space (pos, iv, j, p)
    !if (io%verbose > 1) &
    !  write(io_unit%log,'(a,2i6,2x,3i4,2x,3f7.3)') &
    !    'interpolate: target, source, j, p =', target%id, source%id, j, p
    !call trace%begin ('patch_t%interpolate', itimer=itimer)
    if (source%unsigned(iv)) then
     if (present(pi)) then
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = exp( &
        pi(01)*log(source%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(1),1)) + &
        pi(02)*log(source%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(1),1)) + &
        pi(03)*log(source%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(1),1)) + &
        pi(04)*log(source%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(1),1)) + &
        pi(05)*log(source%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(1),1)) + &
        pi(06)*log(source%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(1),1)) + &
        pi(07)*log(source%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(1),1)) + &
        pi(08)*log(source%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(1),1)) + &
        pi(09)*log(source%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(2),1)) + &
        pi(10)*log(source%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(2),1)) + &
        pi(11)*log(source%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(2),1)) + &
        pi(12)*log(source%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(2),1)) + &
        pi(13)*log(source%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(2),1)) + &
        pi(14)*log(source%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(2),1)) + &
        pi(15)*log(source%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(2),1)) + &
        pi(16)*log(source%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(2),1)))
     else
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = exp( &
       pt(1)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*log(source%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(1),1))    + &
                                    (    p(1))*log(source%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(1),1)))   + &
                         (   p(2))*((1.0-p(1))*log(source%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(1),1))    + &
                                    (    p(1))*log(source%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(1),1))))  + &
              (   p(3))*((1.-p(2))*((1.0-p(1))*log(source%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(1),1))    + &
                                    (    p(1))*log(source%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(1),1)))   + &
                         (   p(2))*((1.0-p(1))*log(source%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(1),1))    + &
                                    (    p(1))*log(source%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(1),1))))) + &
       pt(2)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*log(source%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(2),1))    + &
                                    (    p(1))*log(source%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(2),1)))   + &
                         (   p(2))*((1.0-p(1))*log(source%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(2),1))    + &
                                    (    p(1))*log(source%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(2),1))))  + &
              (   p(3))*((1.-p(2))*((1.0-p(1))*log(source%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(2),1))    + &
                                    (    p(1))*log(source%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(2),1)))   + &
                         (   p(2))*((1.0-p(1))*log(source%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(2),1))    + &
                                    (    p(1))*log(source%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(2),1))))))
     end if
    else
     if (present(pi)) then
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = &
        pi(01)*source%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(1),1) + &
        pi(02)*source%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(1),1) + &
        pi(03)*source%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(1),1) + &
        pi(04)*source%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(1),1) + &
        pi(05)*source%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(1),1) + &
        pi(06)*source%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(1),1) + &
        pi(07)*source%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(1),1) + &
        pi(08)*source%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(1),1) + &
        pi(09)*source%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(2),1) + &
        pi(10)*source%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(2),1) + &
        pi(11)*source%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(2),1) + &
        pi(12)*source%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(2),1) + &
        pi(13)*source%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(2),1) + &
        pi(14)*source%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(2),1) + &
        pi(15)*source%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(2),1) + &
        pi(16)*source%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(2),1)
     else
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = &
       pt(1)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*source%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(1),1)    + &
                                    (    p(1))*source%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(1),1))   + &
                         (   p(2))*((1.0-p(1))*source%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(1),1)    + &
                                    (    p(1))*source%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(1),1)))  + &
              (   p(3))*((1.-p(2))*((1.0-p(1))*source%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(1),1)    + &
                                    (    p(1))*source%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(1),1))   + &
                         (   p(2))*((1.0-p(1))*source%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(1),1)    + &
                                    (    p(1))*source%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(1),1)))) + &
       pt(2)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*source%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(2),1)    + &
                                    (    p(1))*source%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(2),1))   + &
                         (   p(2))*((1.0-p(1))*source%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(2),1)    + &
                                    (    p(1))*source%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(2),1)))  + &
              (   p(3))*((1.-p(2))*((1.0-p(1))*source%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(2),1)    + &
                                    (    p(1))*source%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(2),1))   + &
                         (   p(2))*((1.0-p(1))*source%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(2),1)    + &
                                    (    p(1))*source%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(2),1))))
     end if
    end if
  end if
  if (target%id==io%id_debug .and. io%verbose>2) &
    write(io%output,1) target%id, ii, iv, source%id, j, jt, p, pt, &
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1)
    1 format(i6,1x,3i3,i4,2x,"DBG patch_t%interpolate src:",i6, &
      "   j,p;",3i4,1x,2i4,3f6.2,1x,2f6.2,"   out:",1p,e15.6)
  !call trace%end (itimer)
END SUBROUTINE interpolate

!===============================================================================
!> Interpolate in space and time, for a case where it is safe = the source point
!> is inside the valid domain (inside a distance of half a cell from the source
!> patch boundary for non-staggered points)
!===============================================================================
SUBROUTINE interpolate_specific (source, target, ii, iv, jt, pt)
  class(patch_t)           :: source
  class(patch_t), pointer  :: target
  integer                  :: ii(3), iv, jt(2)
  real                     :: pt(2)
  real                     :: tmp(2,2,2,2)
  !.............................................................................
  real(8)                  :: pos(3)
  integer                  :: j(3)
  real                     :: p(3)
  integer, save            :: itimer=0
  !-----------------------------------------------------------------------------
  pos = target%myposition (ii, iv)
  call source%index_space (pos, iv, j, p)
  !call trace%begin ('patch_t%interpolate', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Pick up 16 values for tri-linear interpolation, divided by density
  !-----------------------------------------------------------------------------
  associate (id=>source%idx%d)
  tmp(:,:,:,:) = source%mem(j(1):j(1)+1,j(2):j(2)+1,j(3):j(3)+1,iv,jt(1:2),1) &
               / source%mem(j(1):j(1)+1,j(2):j(2)+1,j(3):j(3)+1,id,jt(1:2),1)
  target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = &
   pt(1)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*tmp(1,1,1,1)    + &
                                (    p(1))*tmp(2,1,1,1))   + &
                     (   p(2))*((1.0-p(1))*tmp(1,2,1,1)    + &
                                (    p(1))*tmp(2,2,1,1)))  + &
          (   p(3))*((1.-p(2))*((1.0-p(1))*tmp(1,1,2,1)    + &
                                (    p(1))*tmp(2,1,2,1))   + &
                     (   p(2))*((1.0-p(1))*tmp(1,2,2,1)    + &
                                (    p(1))*tmp(2,2,2,1)))) + &
   pt(2)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*tmp(1,1,1,2)    + &
                                (    p(1))*tmp(2,1,1,2))   + &
                     (   p(2))*((1.0-p(1))*tmp(1,2,1,2)    + &
                                (    p(1))*tmp(2,2,1,2)))  + &
          (   p(3))*((1.-p(2))*((1.0-p(1))*tmp(1,1,2,2)    + &
                                (    p(1))*tmp(2,1,2,2))   + &
                     (   p(2))*((1.0-p(1))*tmp(1,2,2,2)    + &
                                (    p(1))*tmp(2,2,2,2))))
  !-----------------------------------------------------------------------------
  ! Multiply by target denssity
  !-----------------------------------------------------------------------------
  target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = &
    target%mem(ii(1),ii(2),ii(3),iv,target%it,1) * &
    target%mem(ii(1),ii(2),ii(3),id,target%it,1)
  end associate
  if (target%id==io%id_debug .and. io%verbose>2) &
    write(io%output,1) target%id, ii, iv, source%id, j, jt, p, pt, &
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1)
    1 format(i6,1x,3i3,i4,2x,"DBG patch_t%interpolate_specific src:", &
      i6,"   j,p;",3i4,1x,2i4,3f6.2,1x,2f6.2,"   out:",1p,e15.6)
  !call trace%end (itimer)
END SUBROUTINE interpolate_specific

SUBROUTINE info
  if (verbose > 0) then
    write (stdout,*) 'guard_zones_t%info: update-t =', update_t
  end if
END SUBROUTINE info

END MODULE guard_zones_mod
