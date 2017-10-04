!> =============================================================================
!>  download_link: takes care of downloads to link%task
!>           same: called for patches on the same level
!>         differ: called for patches on different levels
!> =============================================================================
MODULE download_mod
  USE io_mod
  USE omp_mod
  USE bits_mod
  USE kinds_mod
  USE mesh_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE interpolate_mod
  implicit none
  private
  type, public:: download_t
    integer(8):: n_same=0, n_differ=0
    procedure(interpolator_interface), pointer, nopass:: interpolate => null()
  contains
    procedure, nopass:: download_link
    procedure:: init
  end type
  type(download_t), public:: download
CONTAINS

!> =============================================================================
!> Initialise downloading and interpolating options.
!> =============================================================================
SUBROUTINE init (self)
  class(download_t):: self
  integer:: iostat
  integer, save:: order_interpolator=-1
  logical, save:: first_time=.true.
  namelist /download_params/ order_interpolator
  !-----------------------------------------------------------------------------
  call trace%begin ('download_mod%init')
  !$omp critical (download_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, download_params, iostat=iostat)
    if (io%master) write (*, download_params)
    call self%test
  end if
  self%interpolate => SelectInterpolator(order_interpolator)
  !$omp end critical (download_cr)
  call trace%end
END SUBROUTINE init

  call trace%end
END SUBROUTINE

!> =============================================================================
!> Select the best method for each case.  So far, this only supports static
!> patches at the same level, or static patches one level up or down
!> =============================================================================
SUBROUTINE download_link (link, all_cells, only)
  class(link_t), pointer:: link
  logical, optional:: all_cells
  integer, optional:: only
  !.............................................................................
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: target, source
  integer:: edges
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%download_link', itimer=itimer)
  target => task2patch(link%task)
  if (trim(target%kind)/='rt_solver') then
    if (io%verbose>1) &
      call target%check_density ("before download")
  end if
  if (target%is_set(bits%root_grid)) then
    call target%periodic_grid
  else
    do edges=3,1,-1
      nbor => link%nbor
      do while (associated(nbor))
        source => task2patch(nbor%task)
        if (target%count_edges(source) == edges) then
          if (target%id == io%id_debug) then
            print *, source%id, ' EDGE COUNT:', edges
          end if
          if (trim(target%kind)/='rt_solver') then
            if (source%is_clear(bits%no_download)) then
              if (abs(source%level-target%level) < 2) then
                if (source%level==target%level) then
                  call same (target, source, all_cells, only=only)
                else
                  call different (target, source, all_cells, only=only)
                end if
              end if
            end if
          else
            if (trim(source%kind)=='rt_solver') then
              if (abs(source%level-target%level) < 2) then
                if (source%level==target%level) then
                  call same_rt (target, source)
                else
                  call different_rt (target, source)
                end if
              end if
            end if
          end if
        end if
        nbor => nbor%next
      end do
    end do
  end if
  call target%check_density (" after download")
  call trace%end (itimer)
END SUBROUTINE download_link

!> =============================================================================
!> Find the target index ranges that correspond to the overlap region, for
!> download from source (source has higher, same, or lower level than target).
!>
!> The first index, l(i), returned is the index of the first target point that
!> is inside the authoritative region; i.e., it must lie to the right of the
!> source point r(li).
!> FIXME: To be validated/improved for moving meshes
!> =============================================================================
SUBROUTINE dnload_range (target, source, iv, l, u, offset)
  class(patch_t), pointer :: target, source
  integer                 :: iv, l(3), u(3), offset(3), one(3)
  optional                :: offset
  !.............................................................................
  real(8), dimension(3):: h, dist, rl, ru
  integer:: i, li, ui
  !-----------------------------------------------------------------------------
  call trace%begin ('download_mod::dnload_range', 2)
  ! for collapsed dimensions
  where (source%n > 1)
    one = 1
  elsewhere
    one = 0
  endwhere
  h = source%index_stagger(iv)
  !-----------------------------------------------------------------------------
  ! Map to the first and last point inside the source domain, computing first
  ! the distance btw patch midpoints, then adding the offsets to the li and ui
  ! points. By using rl and ru the procedure becomes independent of no_mans_land.
  !-----------------------------------------------------------------------------
  dist = source%position - target%position + h*(source%mesh%d-target%mesh%d)
  where (source%periodic)
    dist = modulo(dist+0.5_8*source%mesh%b,source%mesh%b) -0.5_8*source%mesh%b
  end where
  !-----------------------------------------------------------------------------
  ! The lf and uf values are floating point relative coordinates corresponding
  ! to li and ui.  These are the boundaries of the 'authoritative domain'
  !-----------------------------------------------------------------------------
  do i=1,3
    rl(i) = dist(i) + source%mesh(i)%r(source%mesh(i)%li+one(i))
    ru(i) = dist(i) + source%mesh(i)%r(source%mesh(i)%ui  )
  end do
  l = target_index (rl)
  u = target_index (ru)
  !-----------------------------------------------------------------------------
  ! Make sure the range is legal.  The u index can become smaller than lb, which
  ! means there is no overlap.  Likewise, the l index can become larger than ub
  !-----------------------------------------------------------------------------
  l = max(target%mesh%lb,l)
  if (target%kind(1:14) == 'zeus_mhd_patch') then
    where (target%mesh%n > 1 .and. h /= 0.0)
      u = min(target%mesh%ub+1,u)
    elsewhere
      u = min(target%mesh%ub,u)
    end where
  else
    u = min(target%mesh%ub,u)
  end if
  l = max(min(l,target%mesh%uo),target%mesh%lb)
  u = max(min(u,target%mesh%ub),target%mesh%lo)
  if (target%id==io%id_debug) print *, 'mk3', l(1), u(1)
  !-----------------------------------------------------------------------------
  ! For same level patch, the integer offset is what should be added to the
  ! target index to get the source index.  For aligned meshes, offset=0, while
  ! for adjacent meshes, the offset is equal to the number of authoritative points.
  !-----------------------------------------------------------------------------
  if (present(offset)) then
    if (source%level==target%level) then
      where (source%mesh%ui > u)
        offset = source%mesh%ui     - u
      else where
        offset = source%mesh%li+one - l
      end where
      dist = source%distance(target)  ! positive if source to the right of target
    else
      print *, 'ERROR: fixed offset not valid for different levels', &
        target%id, target%level, source%id, source%level
      error stop "dnload_range"
    end if
    if (io%verbose>3.or.io_unit%do_validate.or.target%id==io%id_debug) &
      print 1, target%id, target%level, iv, source%id, source%level, l, u, offset
      1 format("patch_t%dnload_range: id,lv,iv:",i6,2i3, &
               "  source,lv,l,u,o:",i6,i3,3(2x,3i4))
  else
    if (io%verbose>3.or.io_unit%do_validate.or.target%id==io%id_debug) &
      print 1, target%id, target%level, iv, source%id, source%level, l, u
  end if
  call trace%end
contains
!> =============================================================================
  function target_index (p, roundup) result (ii)
    real(8)           :: p(3)
    logical, optional :: roundup
    integer           :: ii(3)
    !---------------------------------------------------------------------------
    p = p/target%mesh%d + target%offset                   ! float index
    if (present(roundup)) then
      ii = ceiling(p)                                     ! round up
    else
      ii = p + 0.0001                                     ! integer part
    end if
  end function target_index
END SUBROUTINE dnload_range

!> =============================================================================
!> Copy guard zone info from source patch at the same level
!>
!> same
!>   dnload_range
!>   time_interval
!> =============================================================================
SUBROUTINE same (target, source, all_cells, only)
  class(patch_t), pointer :: target, source
  logical, optional:: all_cells
  integer, optional:: only
  !.............................................................................
  integer                 :: l(3), u(3), o(3), ii(3), jj(3), jt(2)
  real(8)                 :: position(3), xx
  real                    :: pt(2)
  integer                 :: ix, iy, iz, iv, iv1, iv2
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%same', itimer=itimer)
  call source%time_interval (target%time, jt, pt)
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = target%nv
  end if
  do iv=iv1, iv2
    call dnload_range (target, source, iv, l, u, o)
    if (target%id==io%id_debug) &
      print 1, source%id, l, u, o, jt, pt
      1 format("DBG download_t%same:", i6,"  l,u,o:",3(3i4,2x)," jt,pt:",2i4,2x,2f6.2)
    if (source%unsigned(iv)) then
      if (present(all_cells)) then
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              target%mem(ix,iy,iz,iv,target%it,1) = exp( &
                pt(1)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)) + &
                pt(2)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)))
            end do
          end do
        end do
      else
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              ii = [ix,iy,iz]
              jj = ii + o
              target%mem(ix,iy,iz,iv,target%it,1) = merge( &
                exp(pt(1)*log(source%mem(jj(1),jj(2),jj(3),iv,jt(1),1)) + &
                    pt(2)*log(source%mem(jj(1),jj(2),jj(3),iv,jt(2),1))), &
                target%mem(ix,iy,iz,iv,target%it,1), &
                any((ii-target%mesh%li-1) < 0 .or. (ii-target%mesh%ui) > 0))
            end do
          end do
        end do
      end if
    else
      if (present(all_cells)) then
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              target%mem(ix,iy,iz,iv,target%it,1) =  &
                pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
                pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)
            end do
          end do
        end do
      else
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              ii = [ix,iy,iz]
              jj = ii + o
              target%mem(ix,iy,iz,iv,target%it,1) = merge( &
                pt(1)*source%mem(jj(1),jj(2),jj(3),iv,jt(1),1) + &
                pt(2)*source%mem(jj(1),jj(2),jj(3),iv,jt(2),1), &
                target%mem(ix,iy,iz,iv,target%it,1), &
                any((ii-target%mesh%li-1) < 0 .or. (ii-target%mesh%ui) > 0))
            end do
          end do
        end do
      end if
    end if
  end do
  if (io%verbose>2.or.target%id==io%id_debug) &
    print *, target%id, source%id, 'same picked up', product(u-l+1), ' values, mind =', &
    minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,target%it,1))
  call trace%end (itimer)
END SUBROUTINE same

!> =============================================================================
!> Same as above, but for RT_solver type patches
!> =============================================================================
SUBROUTINE same_rt (target, source)
  class(patch_t), pointer :: target, source
  !.............................................................................
  integer                 :: l(3), u(3), o(3), ii(3), jt(2)
  real(8)                 :: position(3)
  real                    :: pt(2)
  integer                 :: ix, iy, iz, iv
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%same_rt', itimer=itimer)
  call source%time_interval (target%time, jt, pt)
  call dnload_range (target, source, 1, l, u, o)
  do iv=1,target%nv
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii = [ix,iy,iz]
          target%mem(ix,iy,iz,iv,target%it,1) = &
             pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
             pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)
          !if (all(ii>target%mesh%li.and.ii<=target%mesh%ui)) print *, &
            !target%id, 'ERROR: internal point modified', source%id, ii, iv
        end do
      end do
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE same_rt

!> =============================================================================
!> Interpolate guard zone info from source patch at a different level
SUBROUTINE different (target, source, all_cells, only)
  class(patch_t), pointer :: target, source
  logical, optional:: all_cells
  integer, optional:: only
  !.............................................................................
  integer                 :: ix, iy, iz, iv, jt(2), iv1, iv2
  integer, dimension(3)   :: l, u, ii, li, jj
  real(8), dimension(3)   :: pos, h
  logical                 :: valid
  integer, save           :: itimer=0
  real                    :: pt(2), pp(3)
  real(kind=KindScalarVar):: qint
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different', itimer=itimer)
  call source%time_interval (target%time, jt, pt)
  if (present(only)) then
    iv1 = only
    iv2 = only
  else
    iv1 = 1
    iv2 = target%nv
  end if
  do iv=iv1,iv2
    call dnload_range (target, source, iv, l, u)
    !---------------------------------------------------------------------------
    ! For no-no-mans-land meshes, use guard zone values for the li index, to
    ! make fluxes up/down symmetric
    !---------------------------------------------------------------------------
    li = target%mesh%li
    if (.not.target%no_mans_land) then
      h = target%index_stagger(iv)
      where (h /= 0.0)
        li = target%mesh%li + 1
      end where
      if (target%id==io%id_debug) print '(a,i3,3(2x,3i4))', 'XX: iv,li', iv,li,l,u
    end if
    if (target%id==io%id_debug) &
      print 1, source%id, (source%position-target%position)/target%size, l, u
      1 format("DBG download_t%different:", i6,"  displ:",3f11.6,"  l,u:",2(3i4,2x))
    associate (q1 => source%mem(:,:,:,:,jt(1),1), q2 => source%mem(:,:,:,:,jt(2),1))
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii=[ix,iy,iz]
          if (present(all_cells) .or. &
              any((ii-li) < 0 .or. (ii-target%mesh%ui) > 0)) then
            if (target%id==io%id_debug .and. io%verbose>1) &
              print *, target%id, 'DBG interpolate_t%different:    using', iv,ii, present(all_cells), (ii-li) < 0
            pos = target%myposition (ii, iv)
            call source%index_space (pos, iv, jj, pp)
            if (source%pervolume(iv)) then
              qint = interpolator%interpolator_specific (q1, jj(1), jj(2), jj(3), iv, pp) * pt(1) &
                   + interpolator%interpolator_specific (q2, jj(1), jj(2), jj(3), iv, pp) * pt(2)
            else
              qint = interpolator%interpolator (q1, jj(1), jj(2), jj(3), iv, pp) * pt(1) &
                   + interpolator%interpolator (q2, jj(1), jj(2), jj(3), iv, pp) * pt(2)
            end if
            target%mem(ix,iy,iz,iv,target%it,1) = qint
          else
            if (target%id==io%id_debug .and. io%verbose>1) &
              print *, target%id, 'DBG interpolate_t%different: skipping', iv,ii, present(all_cells), (ii-li) < 0
            cycle
          end if
        end do
      end do
    end do
    end associate
  end do
  if (io%verbose>2.or.target%id==io%id_debug) &
    print 2, target%id, target%level, source%id, source%level, 'different picked up', &
    product(u-l+1), ' values, mind =', &
    minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,target%it,1)), l, u
  2 format(2(i6,i3),2x,a,i6,a,1p,e12.4,2(2x,3i4))
  call trace%end (itimer)
END SUBROUTINE different

!> =============================================================================
!> Same as above, but for RT_solver type patches
!> =============================================================================
SUBROUTINE different_rt (target, source)
  class(patch_t), pointer :: target, source
  !.............................................................................
  integer                 :: ix, iy, iz, iv, jt(2)
  integer, dimension(3)   :: l, u, ii, lb, ub
  real(8), dimension(3)   :: pos
  logical                 :: valid
  integer, save           :: itimer=0
  real                    :: pt(2)
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different_rt', itimer=itimer)
  lb = target%mesh%lb
  ub = target%mesh%ub
  call source%time_interval (target%time, jt, pt)
  call dnload_range (target, source, 1, l, u)
  do iv=1,target%nv
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii=[ix,iy,iz]
          call source%interpolate (target, ii, iv, jt, pt)
        end do
      end do
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE different_rt

END MODULE
