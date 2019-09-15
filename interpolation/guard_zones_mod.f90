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
  implicit none
  private
  type, public:: guard_zones_t
    procedure(restrict_simple), pointer:: restrict => null()
  contains
    procedure:: init
    procedure:: restrict_simple
    procedure:: restrict_opt
  end type
  integer:: verbose=0
  type(guard_zones_t), public:: guard_zones
CONTAINS

!===============================================================================
!> Choose method and verbosity
!===============================================================================
SUBROUTINE init (self)
  class(guard_zones_t):: self
  character(len=6):: method='opt'
  namelist /guard_zone_params/ method, verbose
  integer:: iostat
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (init_cr)
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, guard_zone_params, iostat=iostat)
      write (io%output, guard_zone_params)
    end if
    !$omp end critical (init_cr)
  end if
  select case (method)
  case ('simple')
    self%restrict => restrict_simple
  case ('opt')
    self%restrict => restrict_opt
  case default
    call io%abort ('guard zone method unknown')
  end select
END SUBROUTINE init

!===============================================================================
!> Compute interpolated guard zone values by setting up interface to the
!> slow but correct patch_t%interpolate() point-by-point interpolation.
!===============================================================================
SUBROUTINE restrict_simple (self, target, source, jt, pt)
  class(guard_zones_t):: self
  class(patch_t), pointer:: target, source
  integer:: jt(2)
  real(4):: pt(2)
  !.............................................................................
  real, dimension(:,:,:), pointer:: tmem, smem, td, sd
  class(mesh_t), pointer:: mesh
  integer:: i1, i2, i3, iv, tl(3), tu(3), sl(3), su(3)
  real(8):: tlc(3), tuc(3), slc(3), suc(3), dist(3)
  real(4):: pl(3), pu(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Interpolate in time and space, point-by-point
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_zones_t%restrict_simple', itimer=itimer)
  call source%time_interval (target%time, jt, pt)
  do iv=1,target%nv
    !---------------------------------------------------------------------------
    ! Target lower and upper corner in real space
    !---------------------------------------------------------------------------
    tlc = target%myposition (target%mesh%lb, iv)
    tuc = target%myposition (target%mesh%ub, iv)
    !---------------------------------------------------------------------------
    ! Map to source index space and limit to source interior
    !---------------------------------------------------------------------------
    sl = source%index_only (tlc)
    su = source%index_only (tuc,roundup=.true.)
    sl = max(sl,source%mesh%li)
    su = min(su,source%mesh%ui)
    !---------------------------------------------------------------------------
    ! Map source points to target legal index space
    !---------------------------------------------------------------------------
    slc = source%myposition(sl, iv) + source%ds
    suc = source%myposition(su, iv)
    tl = target%index_only (slc, iv)
    tu = target%index_only (suc, iv)
    tl = max(tl,target%mesh%lb)
    tu = min(tu,target%mesh%ub)
    if (verbose > 0) &
      write(io_unit%log,'(a,i3,2(2x,3i4),2(2x,3f7.3))') &
        'guard_zones_t%restrict_simple1: iv, tl, tu =', iv, tl, tu
    !call dnload_range (target, source, iv, tl, tu)
    !if (verbose > 0) &
    !  write(io_unit%log,'(a,i3,2(2x,3i4),2(2x,3f7.3))') &
    !    'guard_zones_t%restrict_simple2: iv, tl, tu =', iv, tl, tu
    do i3=tl(3),tu(3); do i2=tl(2),tu(2); do i1=tl(1),tu(1)
      call source%interpolate (target, [i1,i2,i3], iv, jt, pt)
    end do; end do; end do
  end do
  call trace%end (itimer)
END SUBROUTINE restrict_simple

!===============================================================================
!> Compute interpolated guard zone values
!===============================================================================
SUBROUTINE restrict_opt (self, target, source, jt, pt)
  class(guard_zones_t):: self
  class(patch_t), pointer:: target, source
  integer:: jt(2)
  real(4):: pt(2)
  !.............................................................................
  real, dimension(:,:,:), pointer:: tmem, smem, td, sd
  class(mesh_t), pointer:: mesh
  integer:: i, i1, i2, i3, j, k, iv, tl(3), tu(3), sl(3), su(3), j0(3)
  real(8):: tlc(3), tuc(3), slc(3), suc(3)
  real(4):: p(3), p1, p2, p3, w(8)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Map the source corners into target index space
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_zones_t%restrict_opt', itimer=itimer)
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
    sl = source%index_only (tlc, iv)
    su = source%index_only (tuc, iv)
    sl = max(sl,source%mesh%li)
    su = min(su,source%mesh%ui-1)
    !---------------------------------------------------------------------------
    ! Map those source points back to the target index space, rounding up to
    ! find the target points between the source points.
    !---------------------------------------------------------------------------
    slc = source%myposition(sl, iv)
    suc = source%myposition(su, iv)
    tl = target%index_only (slc, iv, roundup=.true.)
    tu = target%index_only (suc, iv, roundup=.true.)
    !---------------------------------------------------------------------------
    ! The brackets below should not really be needed, since the points should
    ! all be inside the target legal range
    !---------------------------------------------------------------------------
    tl = max(tl,target%mesh%lb)
    tu = min(tu,target%mesh%ub)
    if (verbose > 0) then
      write(io_unit%log,1) &
        'guard_zones_t%restrict_opt: iv, tlc, sl, slc, tl =', &
        iv, tlc, sl, slc, tl
      1 format(a,i3,2(2(2x,3f8.4,2x,3i4)))
      write(io_unit%log,1) &
        'guard_zones_t%restrict_opt: iv, tuc, su, suc, tu =', &
        iv, tuc, su, suc, tu
    end if
    !---------------------------------------------------------------------------
    ! Compute indices and interpolation weights for the first target point, and
    ! save the weights for the 8 nearest neighbor points in source space
    !---------------------------------------------------------------------------
    call source%index_space (target%myposition(tl), iv, j0, p)
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
        'guard_zones_t%restrict_opt: iv, tl, tu, j0 =', &
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
  call trace%end (itimer)
END SUBROUTINE restrict_opt

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
  real, allocatable:: f(:,:), tot(:)
  integer:: i, j, n, i1 , i2, i3, j1, j2, j3, jj(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('guard_zones_t%update')
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
  call trace%end
END SUBROUTINE update

!===============================================================================
!> OBSOLETE
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
  where (source%periodic)
    dist = modulo(dist+0.5_8*source%mesh%b,source%mesh%b) - 0.5_8*source%mesh%b
  end where
  !-----------------------------------------------------------------------------
  ! Add the distances from the source center to the start and end of the
  ! authoritative domain of each variable.
  !-----------------------------------------------------------------------------
  debug = source%level /= target%level
  la = dist - source%size/2d0
  ua = dist + source%size/2d0
  if (source%level /= target%level) then
    do i=1,3
      sm => source%mesh(i)
      tm => target%mesh(i)
      la(i) = la(i) + sm%h(iv)*sm%d - tm%h(iv)*tm%d
      ua(i) = ua(i) + sm%h(iv)*sm%d - tm%h(iv)*tm%d
    end do
  end if
  l = ceiling(la/target%mesh%d + target%mesh%o)
  u =   floor(ua/target%mesh%d + target%mesh%o)
  !-----------------------------------------------------------------------------
  ! Make sure the range is legal.  The u index can become smaller than lb, which
  ! means there is no overlap.  Likewise, the l index can become larger than ub
  !-----------------------------------------------------------------------------
  l = max(l,target%mesh%lb)
  u = min(u,target%mesh%ub)
  !-----------------------------------------------------------------------------
  ! Account for collapsed dimensions
  !-----------------------------------------------------------------------------
  where (target%n == 1)
    l = min(l,target%mesh%ub)
    u = max(u,target%mesh%lb)
  end where
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
      !call mpi%abort('dnload_range')
    end if
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

END MODULE guard_zones_mod
