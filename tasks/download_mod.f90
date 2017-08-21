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
  use interpolate_mod
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
  integer:: order_interpolator, iostat
  logical:: first_time=.true.
  namelist /download_params/ order_interpolator
  !-----------------------------------------------------------------------------
  call trace%begin ('download_mod%init')

  if (first_time) then
    order_interpolator = -1
    first_time = .false.
    rewind (io%input)
    read (io%input, download_params, iostat=iostat)
    if (io%master) write (*, download_params)
    self%interpolate => SelectInterpolator(order_interpolator)
  end if

  call trace%end
END SUBROUTINE

!> =============================================================================
!> Select the best method for each case.  So far, this only supports static
!> patches at the same level, or static patches one level up or down
!> =============================================================================
SUBROUTINE download_link (link, all_cells)
  class(link_t), pointer:: link
  logical, optional:: all_cells
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: target, source
  integer:: edges
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%download_link')
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
                  call same (target, source, all_cells)
                else
                  call different (target, source, all_cells)
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
  call trace%end
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
  integer                 :: iv, l(3), u(3), offset(3)
  optional                :: offset
  !.............................................................................
  real(8), dimension(3):: h, dist, rl, ru
  integer:: i, li, ui
  !-----------------------------------------------------------------------------
  call trace%begin ('download_mod::dnload_range')
  h = source%index_stagger(iv)
  !-----------------------------------------------------------------------------
  ! Map to the first and last point inside the source domain, computing first
  ! the distance btw patch midpoints, then adding the offsets to the li and ui
  ! points. By using rl and ru the procedure becomes independent of no_mans_land.
  !-----------------------------------------------------------------------------
  dist = source%position - target%position + h*(source%mesh%d-target%mesh%d)
  !-----------------------------------------------------------------------------
  ! The lf and uf values are floating point relative coordinates corresponding
  ! to li and ui.  These are the boundaries of the 'authoritative domain'
  !-----------------------------------------------------------------------------
  l = target%index_only (dist + source%mesh%lf)
  u = target%index_only (dist + source%mesh%uf)
  !-----------------------------------------------------------------------------
  ! Make sure the range is legal.  The u index can become smaller than lb, which
  ! means there is no overlap.  Likewise, the l index can become larger than ub
  !-----------------------------------------------------------------------------
  l = max(target%mesh%lb,l)
  u = min(target%mesh%ub,u)
  !-----------------------------------------------------------------------------
  ! For same level patch, the integer offset is what should be added to the
  ! target index to get the source index.  For aligned meshes, offset=0, while
  ! for adjacent meshes, the offset is equal to the number of authoritative points.
  !-----------------------------------------------------------------------------
  if (present(offset)) then
    if (source%level==target%level) then
      where (source%mesh%ui > u)
        offset = source%mesh%ui - u
      else where
        offset = source%mesh%li - l
      end where
      dist = source%distance(target)  ! positive if source to the right of target
    else
      print *, 'ERROR: fixed offset not valid for different levels', &
        target%id, target%level, source%id, source%level
      error stop "dnload_range"
    end if
    if (io%verbose>1.or.target%id==io%id_debug) &
      print 1, target%id, target%level, iv, source%id, source%level, l, u, offset
      1 format("patch_t%dnload_range: id,lv,iv:",i6,2i3, &
               "  source,lv,l,u,o:",i6,i3,3(2x,3i4))
  else
    if (io%verbose>1.or.target%id==io%id_debug) &
      print 1, target%id, target%level, iv, source%id, source%level, l, u
  end if
  call trace%end
END SUBROUTINE dnload_range

!> =============================================================================
!> Copy guard zone info from source patch at the same level
!>
!> same
!>   dnload_range
!>   time_interval
!> =============================================================================
SUBROUTINE same (target, source, all_cells)
  class(patch_t), pointer :: target, source
  logical, optional:: all_cells
  !.............................................................................
  integer                 :: l(3), u(3), o(3), ii(3), jt(2)
  real(8)                 :: position(3), xx
  real                    :: pt(2)
  integer                 :: ix, iy, iz, iv
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%same', itimer=itimer)
  call source%time_interval (target%time, jt, pt)
  do iv=1,target%nv
    call dnload_range (target, source, iv, l, u, o)
    if (target%id==io%id_debug) &
      print 1, source%id, l, u, o, jt, pt
      1 format("DBG download_t%same:", i6,"  l,u,o:",3(3i4,2x)," jt,pt:",2i4,2x,2f6.2)
    if (source%unsigned(iv)) then
      do iz=l(3),u(3)
        do iy=l(2),u(2)
          do ix=l(1),u(1)
            ii = [ix,iy,iz]
            if (present(all_cells) .or. &
                any((ii-target%mesh%li) < 0 .or. (ii-target%mesh%ui) > 0)) then
              target%mem(ix,iy,iz,iv,target%it,1) = exp( &
                 pt(1)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)) + &
                 pt(2)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)))
            else
              if (target%id==io%id_debug .and. io%verbose>1) &
                print *, target%id, 'DBG interpolate_t%same: skipping', ii
              cycle
            end if
          end do
        end do
      end do
    else
      do iz=l(3),u(3)
        do iy=l(2),u(2)
          do ix=l(1),u(1)
            ii = [ix,iy,iz]
            if (present(all_cells) .or. &
                any((ii-target%mesh%li) < 0 .or. (ii-target%mesh%ui) > 0)) then
              target%mem(ix,iy,iz,iv,target%it,1) = &
                pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
                pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)
            else
              if (target%id==io%id_debug .and. io%verbose>1) &
                print *, target%id, 'DBG interpolate_t%same: skipping', ii
              cycle
            end if
          end do
        end do
      end do
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
  call trace%begin('download_t%same', itimer=itimer)
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
SUBROUTINE different (target, source, all_cells)
  class(patch_t), pointer :: target, source
  logical, optional:: all_cells
  !.............................................................................
  integer                 :: ix, iy, iz, iv, jt(2)
  integer, dimension(3)   :: l, u, ii, lb, ub
  real(8), dimension(3)   :: pos
  logical                 :: valid
  integer, save           :: itimer=0
  real                    :: pt(2)
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different', itimer=itimer)
  lb = target%mesh%lb
  ub = target%mesh%ub
  call source%time_interval (target%time, jt, pt)
  do iv=1,target%nv
    call dnload_range (target, source, iv, l, u)
    if (target%id==io%id_debug) &
      print 1, source%id, (source%position-target%position)/target%size, l, u
      1 format("DBG download_t%different:", i6,"  displ:",3f11.6,"  l,u:",2(3i4,2x))
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii=[ix,iy,iz]
          if (present(all_cells) .or. &
              any((ii-target%mesh%li) < 0 .or. (ii-target%mesh%ui) > 0)) then
            if (source%pervolume(iv)) then
              call source%interpolate_specific (target, ii, iv, jt, pt)
            else
              call source%interpolate (target, ii, iv, jt, pt)
            end if
          else
            if (target%id==io%id_debug .and. io%verbose>1) &
              print *, target%id, 'DBG interpolate_t%different: skipping', ii
            cycle
          end if
        end do
      end do
    end do
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
  call trace%begin('download_t%different', itimer=itimer)
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
