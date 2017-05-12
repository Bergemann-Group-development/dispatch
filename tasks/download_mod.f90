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
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%download_link')
  target => task2patch(link%task)
  if (io%verbose>1) &
    call target%check_density ("before download")
  if (target%is_set(bits%root_grid)) then
    call target%periodic_grid
  else
    nbor => link%nbor
    associate (d=>target%mem(:,:,:,target%i_d,target%it,1))
    do while (associated(nbor))
      source => task2patch(nbor%task)
      if (source%is_clear(bits%no_download)) then
        if (abs(source%level-target%level) < 2) then
          if (source%level==target%level) then
            call same (target, source)
          else
            call different (target, source)
          end if
        end if
      end if
      nbor => nbor%next
    end do
    end associate
  end if
  call target%check_density (" after download")
  call trace%end
END SUBROUTINE download_link

!> =============================================================================
!> Find the target index ranges that correspond to the overlap region, for
!> download from source (source has higher, same, or lower level than target)
!> FIXME: To be validated/improved for moving meshes
!> =============================================================================
SUBROUTINE dnload_range (target, source, iv, l, u, offset)
  class(patch_t), pointer :: target, source
  integer                 :: iv, l(3), u(3), offset(3)
  optional                :: offset
  !.............................................................................
  real, dimension(3):: f, pos1, pos2, dist
  real(8):: eps(3), h(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('download_mod::dnload_range')
  eps = 0.001_8*target%ds
  h = source%index_stagger (iv)
  !-----------------------------------------------------------------------------
  ! Map to the first and last point inside the source domain
  !-----------------------------------------------------------------------------
  l = target%index_only (source%mesh%p+source%mesh%lf-eps, iv, roundup=.true.)
  u = target%index_only (source%mesh%p+source%mesh%uf+eps, iv)
  !-----------------------------------------------------------------------------
  ! Make sure the range is legal, and includes all guard cell layers
  !-----------------------------------------------------------------------------
  l = max(target%mesh%lb,min(target%mesh%uo,l))
  u = max(target%mesh%lo,min(target%mesh%ub,u))
  !-----------------------------------------------------------------------------
  ! If the source is to the left or right of the target in a direction, limit
  ! the region to the ghost zones on that side.
  !-----------------------------------------------------------------------------
  if (present(offset)) then
    if (source%level==target%level) then
      offset = 0
      dist = source%distance(target)  ! positive if source to the right of target
      where (dist > source%ds)
        offset = -target%mesh%n+1
        l = target%mesh%uo
        u = target%mesh%ub
      else where (dist < -source%ds)
        offset = target%mesh%n-1
        l = target%mesh%lb
        u = target%mesh%lo
      else where
        l = target%mesh%li
        u = target%mesh%ui
      end where
    else
      print *, 'WARNING: fixed offset not valid for different levels', &
        target%id, target%level, source%id, source%level
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
SUBROUTINE same (target, source)
  class(patch_t), pointer :: target, source
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
            target%mem(ix,iy,iz,iv,target%it,1) = exp( &
               pt(1)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)) + &
               pt(2)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)))
            !if (all(ii>target%mesh%li.and.ii<=target%mesh%ui)) print *, &
              !target%id, 'ERROR: internal point modified', source%id, ii, iv
          end do
        end do
      end do
    else
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
    end if
  end do
  if (io%verbose>2.or.target%id==io%id_debug) &
    print *, target%id, source%id, 'same picked up', product(u-l+1), ' values, mind =', &
    minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,target%it,1))
  call trace%end (itimer)
END SUBROUTINE same

!> =============================================================================
!> Interpolate guard zone info from source patch at a different level
SUBROUTINE different (target, source)
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
  do iv=1,target%nv
    call dnload_range (target, source, iv, l, u)
    if (target%id==io%id_debug) &
      print 1, source%id, (source%position-target%position)/target%size, l, u
      1 format("DBG download_t%different:", i6,"  displ:",3f11.6,"  l,u:",2(3i4,2x))
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii=[ix,iy,iz]
          call source%interpolate (target, ii, iv, jt, pt)
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

END MODULE
