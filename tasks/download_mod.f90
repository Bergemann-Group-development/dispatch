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
  ! If the source and target have the same size, then map point to point
  !-----------------------------------------------------------------------------
  if (source%level == target%level) then
    l = target%index_only (source%myposition (source%mesh%li, iv), iv)
    u = target%index_only (source%myposition (source%mesh%ui, iv), iv)
  !-----------------------------------------------------------------------------
  ! If the source is smaller, map the first possible source point to the next
  ! available target point to the right (or coinciding)
  !-----------------------------------------------------------------------------
!  else if (source%level > target%level) then
  else
    l = target%index_only (source%mesh%p+source%mesh%lf-eps, iv, roundup=.true.)
    u = target%index_only (source%mesh%p+source%mesh%uf-eps , iv)
  !-----------------------------------------------------------------------------
  ! If the source is larger, use the entire range as a first guess
  !-----------------------------------------------------------------------------
!  else if (source%level < target%level) then
!    l = target%index_only (source%mesh%p + source%mesh%lf + 0.5_8*target%ds + eps, iv)
!    u = target%index_only (source%mesh%p + source%mesh%uf - 0.5_8*target%ds + eps, iv)
  end if
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
      where (dist > 0.0)
        offset = -target%mesh%n
      else where (dist < 0)
        offset = target%mesh%n
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
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii = [ix,iy,iz]
          target%mem(ix,iy,iz,iv,target%it,1) = &
             pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
             pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)
          if (all(ii>target%mesh%li.and.ii<=target%mesh%ui)) print *, &
            target%id, 'ERROR: internal point modified', source%id, ii, iv
        end do
      end do
    end do
  end do
  if (io%verbose>2.or.target%id==io%id_debug) &
    print *, target%id, source%id, 'same picked up', product(u-l+1), ' values, mind =', &
    minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,target%it,1))
  call trace%end (itimer)
END SUBROUTINE same

SUBROUTINE edge_fill (target, ix, iy, iz, iv)
  class(patch_t), pointer :: target
  integer        :: ix, iy, iz, iv, xb, yb, zb, l(3), u(3)
  logical        :: xedge, yedge, zedge, xlo, xhi, ylo, yhi, zlo, zhi
  !-----------------------------------------------------------------------------
  l = target%mesh%lb
  u = target%mesh%ub
  associate (mem => target%mem(:,:,:,iv,target%it,1))
  if (ix==l(1)+1) then
    xlo = .true.
    xb = l(1)
  end if
  if (ix==u(1)-1) then
    xhi = .true.
    xb = u(1)
  end if
  if (iy==l(2)+1) then
    ylo = .true.
    yb = l(2)
  end if
  if (iy==u(2)-1) then
    yhi = .true.
    yb = u(2)
  end if
  if (iz==l(3)+1) then
    zlo = .true.
    zb = l(3)
  end if
  if (iz==u(3)-1) then
    zhi = .true.
    zb = u(3)
  end if
  xedge = xlo.or.xhi
  yedge = ylo.or.yhi
  zedge = zlo.or.zhi
  if (xedge.and.yedge.and.zedge) then
    mem(xb,iy,iz) = 0.5*(mem(xb,iy-1,iz)+mem(xb,iy+1,iz))
    mem(ix,yb,iz) = 0.5*(mem(ix-1,yb,iz)+mem(ix+1,yb,iz))
    mem(ix,iy,zb) = 0.5*(mem(ix-1,iy,zb)+mem(ix+1,iy,zb))
  else if (xedge.and.yedge) then
    mem(xb,iy,iz) = 0.5*(mem(xb,iy-1,iz)+mem(xb,iy+1,iz))
    mem(ix,yb,iz) = 0.5*(mem(ix-1,yb,iz)+mem(ix+1,yb,iz))
  else if (yedge.and.zedge) then
    mem(ix,yb,iz) = 0.5*(mem(ix,yb,iz-1)+mem(ix,yb,iz+1))
    mem(ix,iy,zb) = 0.5*(mem(ix,iy-1,zb)+mem(ix,iy+1,zb))
  else if (zedge.and.xedge) then
    mem(ix,iy,zb) = 0.5*(mem(ix-1,iy,zb)+mem(ix+1,iy,zb))
    mem(xb,iy,iz) = 0.5*(mem(xb,iy,iz-1)+mem(xb,iy,iz+1))
  end if
  if (xedge) then
    mem(ix,iy,iz) = 0.5*(mem(ix-1,iy,iz)+mem(ix+1,iy,iz))
  end if
  if (yedge) then
    mem(ix,iy,iz) = 0.5*(mem(ix,iy-1,iz)+mem(ix,iy+1,iz))
  end if
  if (zedge) then
    mem(ix,iy,iz) = 0.5*(mem(ix,iy,iz-1)+mem(ix,iy,iz+1))
  end if
  end associate
END SUBROUTINE edge_fill

!> =============================================================================
!> Fix the problematic layer just inside the boundary of the ghost zones
!> =============================================================================
SUBROUTINE edge_fill2 (self)
  class(patch_t):: self
  integer:: ix, iy, iz, it0, it1, it2
  integer, dimension(3):: lb, lo, li, ub, uo, ui
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin ('patch_t%edge_fill2', itimer=itimer)
  lb = self%mesh%lb
  li = self%mesh%li
  lo = self%mesh%lo
  ub = self%mesh%ub
  ui = self%mesh%ui
  uo = self%mesh%uo
  !-----------------------------------------------------------------------------
  ! XY-layers
  !-----------------------------------------------------------------------------
  self%mem(li(1):ui(1),li(2):ui(2),lb(3)+1,:,it0,1) = 0.5*( &
  self%mem(li(1):ui(1),li(2):ui(2),lb(3)  ,:,it1,1) + &
  self%mem(li(1):ui(1),li(2):ui(2),lb(3)+2,:,it2,1) )
  self%mem(li(1):ui(1),li(2):ui(2),ub(3)-1,:,it0,1) = 0.5*( &
  self%mem(li(1):ui(1),li(2):ui(2),ub(3)  ,:,it1,1) + &
  self%mem(li(1):ui(1),li(2):ui(2),ub(3)-2,:,it2,1) )
  !-----------------------------------------------------------------------------
  ! YZ-layers
  !-----------------------------------------------------------------------------
  self%mem(lb(1)+1,li(2):ui(2),li(3):ui(3),:,it0,1) = 0.5*( &
  self%mem(lb(1)  ,li(2):ui(2),li(3):ui(3),:,it1,1) + &
  self%mem(lb(1)+2,li(2):ui(2),li(3):ui(3),:,it2,1) )
  self%mem(ub(1)-1,li(2):ui(2),li(3):ui(3),:,it0,1) = 0.5*( &
  self%mem(ub(1)  ,li(2):ui(2),li(3):ui(3),:,it1,1) + &
  self%mem(ub(1)-2,li(2):ui(2),li(3):ui(3),:,it2,1) )
  !-----------------------------------------------------------------------------
  ! ZX-layers
  !-----------------------------------------------------------------------------
  self%mem(li(1):ui(1),lb(2)+1,li(3):ui(3),:,it0,1) = 0.5*( &
  self%mem(li(1):ui(1),lb(2)  ,li(3):ui(3),:,it1,1) + &
  self%mem(li(1):ui(1),lb(2)+2,li(3):ui(3),:,it2,1) )
  self%mem(li(1):ui(1),ub(2)-1,li(3):ui(3),:,it0,1) = 0.5*( &
  self%mem(li(1):ui(1),ub(2)  ,li(3):ui(3),:,it1,1) + &
  self%mem(li(1):ui(1),ub(2)-2,li(3):ui(3),:,it2,1) )
  call trace_end (itimer)
END SUBROUTINE edge_fill2

!> =============================================================================
!> Interpolate guard zone info from source patch at a different level
!>
!> different
!>   dnload_range
!>   interpolate_different
!>     myposition
!>     is_valid_at
!>       index_stagger
!>     interpolate
!>       index_space
!>         time_interval
!> =============================================================================
SUBROUTINE different (target, source)
  class(patch_t), pointer :: target, source
  !.............................................................................
  integer                 :: ix, iy, iz, iv
  integer, dimension(3)   :: l, u, ii, lb, ub
  real(8), dimension(3)   :: pos
  logical                 :: valid
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different', itimer=itimer)
  lb = target%mesh%lb
  ub = target%mesh%ub
  do iv=1,target%nv
    call dnload_range (target, source, iv, l, u)
    if (target%id==io%id_debug) &
      print 1, source%id, (source%position-target%position)/target%size, l, u
      1 format("DBG download_t%different:", i6,"  displ:",3f11.6,"  l,u:",2(3i4,2x))
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii=[ix,iy,iz]
          call source%interpolate (target, ii, iv)
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
