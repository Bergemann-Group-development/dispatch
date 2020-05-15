!===============================================================================
!>  download_link: takes care of downloads to link%task
!>           same: called for patches on the same level
!>         differ: called for patches on different levels
!===============================================================================
MODULE download_mod
  USE iso_fortran_env, only: int8
  USE io_mod
  USE omp_mod
  USE mpi_mod
  USE bits_mod
  USE kinds_mod
  USE mesh_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE scene_mod
  USE interpolate_mod
  implicit none
  private
  type, public:: download_t
    integer(8):: n_same=0, n_differ=0
    integer:: order_interpolator
    procedure(interpolator_interface), pointer, nopass:: interpolate => null()
  contains
    procedure:: download_link
    procedure:: init
    procedure:: test
  end type
  type(download_t), public:: download
CONTAINS

!===============================================================================
!> Initialise downloading and interpolating options.
!===============================================================================
SUBROUTINE init (self)
  class(download_t):: self
  integer:: iostat
  integer, save:: order_interpolator=-1
  logical, save:: first_time=.true.
  namelist /download_params/ order_interpolator
  !-----------------------------------------------------------------------------
  call trace%begin ('download_t%init')
  !$omp critical (download_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, download_params, iostat=iostat)
    call io%namelist_warning ('download_params', iostat)
    write (io%namelist, download_params)
  end if
  self%order_interpolator = order_interpolator
  self%interpolate => SelectInterpolator(order_interpolator)
  !$omp end critical (download_cr)
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialise downloading and interpolating options.
!===============================================================================
SUBROUTINE test (self)
  class(download_t):: self
  !.............................................................................
  class(patch_t), pointer :: source(:), target, src
  integer, dimension(3)  :: l, u, o, ii
  integer                :: i, j, k, iv, dir
  real(8), dimension(3)  :: pos
  real(4), dimension(3)  :: pp
  class(mesh_t), pointer :: mesh
  logical                :: ok, allok
  integer(int8), dimension(:,:,:,:), allocatable:: filled
  !-----------------------------------------------------------------------------
  call trace%begin ('download_t%test')
  self%interpolate => SelectInterpolator(-1)
  !-----------------------------------------------------------------------------
  ! Generate meshes with three resolutions, differing with a factor of three,
  ! with and w/o staggering
  !-----------------------------------------------------------------------------
  allocate (target)
  call target%setup (size=0.1_8*[1,1,1], n=[16,16,16], nv=2, nt=2, nw=1)
  ! Array to track filling of cells; `-1` means unfilled.
  allocate (filled(target%gn(1),target%gn(2),target%gn(3),target%nv))
  filled = -1
  !-----------------------------------------------------------------------------
  ! Loop over the three resolutions for source patches
  !-----------------------------------------------------------------------------
  allocate (source(3))
  do i=1,3
    call source(i)%setup (size=0.1_8*[1,1,1]*(scene%f_zoom)**(i-2), n=[16,16,16], nv=2, nt=2, nw=1)
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
        print *, 'source size, stagger:', source(i)%size(1), mesh%h(iv)
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
          print *,' '
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
        print "('source size:',i2,3x,'position:',i2,2x,'levels:',2i2)", i, j, src%level, target%level
      call set_values (src)
      call set_values (target)
      if (src%level==target%level) then
        call same (target, src, filled, all_cells=.false.)
      else
        call different_prolong (target, src, filled, all_cells=.false.)
      end if
      call check_values (target, ok)
      allok = allok.and.ok
    end do
  end do
  if (io%master) print *,io%hl
  if (allok) then
    if (io%master) then
      print *,'download_t%test passed'
      print *,io%hl
    end if
  else
    print *,mpi%rank,'download_t%test failed'
    if (io%master) print *,io%hl
    call mpi%abort('download_t%test::check_values')
  end if
  !-----------------------------------------------------------------------------
  ! Clean up
  !-----------------------------------------------------------------------------
  deallocate (filled)
  call target%delete
  do i=1,3
    call source(i)%delete
  end do
  call target%init_id_reset
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
      source%position = target%position - 0.5_8*(source%size+target%size) + 2.*abs(source%ds-target%ds)
    case (2)
      source%position = target%position - 0.5_8*(source%size+target%size) + 2.*abs(source%ds-target%ds)
    case (3)
      source%position = target%position
    case (4)
      source%position = target%position + 0.5_8*(source%size+target%size) - 2.*abs(source%ds-target%ds)
    case (5)
      source%position = target%position + 0.5_8*(source%size+target%size) - 2.*abs(source%ds-target%ds)
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
    print *,'check_values: n_zero, n_ok, n_nok =', n
  if (n(3) == 0) then
    ok = .true.
  else
    ok = .false.
    print *,'check_values: n_zero, n_ok, n_nok =', n
  end if
end subroutine check_values
!===============================================================================
END SUBROUTINE test

!===============================================================================
!> Select the best method for each case.  So far, this only supports static
!> patches at the same level, or static patches one level up or down
!===============================================================================
SUBROUTINE download_link (self, link, all_cells, only)
  class(download_t):: self
  class(link_t), pointer:: link
  logical, optional:: all_cells
  integer, optional:: only
  !.............................................................................
  class(link_t), pointer:: nbor, sorted_head
  class(patch_t), pointer:: target, source
  integer:: edges
  integer, save:: itimer=0
  integer(kind=int8), dimension(:,:,:,:), allocatable:: filled
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%download_link', itimer=itimer)
  target => task2patch(link%task)
  !-----------------------------------------------------------------------------
  ! Allocate and initialise array to track filling of cells.
  ! A value of `-1` defines unfilled; all other values denote the *level* of the
  ! patch that filled that cell.
  !-----------------------------------------------------------------------------
  allocate (filled(target%gn(1),target%gn(2),target%gn(3),target%nv))
  filled = -1

  if (trim(target%kind)/='rt_solver') then
    if (io%verbose>1) &
      call target%check_density ("before download")
  end if
  if (target%is_periodic()) call target%periodic_grid

  ! Take the nbor list and produce a list sorted by decreasing distance from `target`.
  !!!call sort_nbors_by_distance (link,sorted_head)

  ! Walk through the sorted neighbour list
  !!!nbor => sorted_head
  nbor => link%nbor
  do while (associated(nbor))
    source => task2patch(nbor%task)
    if (trim(target%kind)/='rt_solver') then
      if (source%is_clear(bits%no_download)) then
        if (source%level==target%level) then
          call same (target, source, filled, all_cells, only=only)
        else if ((source%level-target%level) == -1) then
          call different_prolong (target, source, filled, all_cells, only=only)
        else if ((source%level-target%level) == 1) then
          call different_restrict (target, source, filled, only=only)
        end if
      end if
    else
      if (trim(source%kind)=='rt_solver') then
        if (abs(source%level-target%level) < 2) then
          if (source%level==target%level) then
            call same_rt (target, source, filled)
          else
            call different_rt (target, source, filled)
          end if
        end if
      end if
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
  deallocate (filled)
  !!!nbor => sorted_head
  !!!if (associated(nbor)) then
  !!!  call clean_sorted_nbor_list (sorted_head)
  !!!else
  !!!  print *,' sorted_head not associated'
  !!!end if

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
  do i=1,3
    sm => source%mesh(i)
    tm => target%mesh(i)
    !la(i) = dist(i) + sm%r(sm%li) + sm%h(iv)*sm%d - tm%h(iv)*tm%d
    !ua(i) = dist(i) + sm%r(sm%ui) + sm%h(iv)*sm%d - tm%h(iv)*tm%d
    la(i) = dist(i) - 0.5*source%size(i) + sm%h(iv)*sm%d - tm%h(iv)*tm%d
    ua(i) = dist(i) + 0.5*source%size(i) + sm%h(iv)*sm%d - tm%h(iv)*tm%d
  end do
  l = target_index (la)
  u = target_index (ua)
  !-----------------------------------------------------------------------------
  ! Make sure the range is legal.  The u index can become smaller than lb, which
  ! means there is no overlap.  Likewise, the l index can become larger than ub
  !-----------------------------------------------------------------------------
  l = max(target%mesh%lb,l)
  u = min(target%mesh%ub,u)
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
      where      (source%mesh%ui > u)
        offset = source%mesh%ui - u
      else where (source%mesh%li < l)
        offset = source%mesh%li - l
      else where
        offset = 0
      end where
    else
      print *, 'ERROR: fixed offset not valid for different levels', &
        target%id, target%level, source%id, source%level
      call mpi%abort('dnload_range')
    end if
    if (io%verbose>3.or.io_unit%do_validate.or.target%id==io%id_debug) &
      print 1, target%id, target%level, iv, source%id, source%level, l, u, offset
      1 format("download_t%dnload_range: id,lv,iv:",i6,2i3, &
               "  source,lv,l,u,o:",i6,i3,3(2x,3i4))
  else
    if (io%verbose>3.or.io_unit%do_validate.or.target%id==io%id_debug) &
      print 1, target%id, target%level, iv, source%id, source%level, l, u
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
SUBROUTINE same (target, source, filled, all_cells, only)
  class(patch_t), pointer :: target, source
  integer(int8), dimension(:,:,:,:):: filled
  logical, optional       :: all_cells
  integer, optional       :: only
  !.............................................................................
  integer                 :: l(3), u(3), o(3), ii(3), jt(2), i
  real(8)                 :: position(3), xx
  real                    :: pt(2)
  integer                 :: ix, iy, iz, iv, iv1, iv2, li(3), ui(3)
  integer, save           :: itimer=0
  logical                 :: all_cells_l
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
  if (present(all_cells)) then
    all_cells_l = all_cells
  else
    all_cells_l = .false.
  end if
  do iv=iv1, iv2
    if (iv == target%idx%phi) then
      if (target%poisson_time > source%poisson_time + source%dtime*source%grace) then
        if (io%verbose>1) print *, target%id, ' not using same level potential from', source%id
        cycle
      else
        if (io%verbose>1) print *, target%id, '     using same level potential from', source%id
      end if
    end if
    if (iv == target%idx%dphi) cycle ! `dphi` is part of `mem` when `time_order` = 0; do not download it.
    li = target%mesh%li
    ui = target%mesh%ui
    call dnload_range (target, source, iv, l, u, o)
    !---------------------------------------------------------------------------
    ! Extend the pick-up region to fill the guard zones also in the BC regions
    !---------------------------------------------------------------------------
    do i=1,3
      if (source%mesh(i)%lower_boundary) l(i) = target%mesh(i)%lb
      if (source%mesh(i)%upper_boundary) u(i) = target%mesh(i)%ub
      !if (i==2) then
        !if (source%boundaries%is_set (bits%yl)) l(i) = target%mesh(i)%lb
        !if (source%boundaries%is_set (bits%yu)) u(i) = target%mesh(i)%ub
      !end if
!      if (target%mesh(i)%h(iv) /= 0.0) li(i) = li(i)+1
    end do
    if (target%id==io%id_debug) &
      print 1, source%id, l, u, o, jt, pt
      1 format("DBG download_t%same:", i6,"  l,u,o:",3(3i4,2x)," jt,pt:",2i4,2x,2f6.2)
    if (source%unsigned(iv)) then
      if (all_cells_l) then
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              !!!if (filled(ix,iy,iz,iv) > source%level) cycle
              target%mem(ix,iy,iz,iv,target%it,1) = exp( &
                pt(1)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)) + &
                pt(2)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)))
              filled(ix,iy,iz,iv) = source%level
            end do
          end do
        end do
      else
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              ii = [ix,iy,iz]
              !!!if (filled(ix,iy,iz,iv) > source%level) cycle
              target%mem(ix,iy,iz,iv,target%it,1) = &
                merge(target%mem(ix,iy,iz,iv,target%it,1), &
                  exp(pt(1)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)) + &
                      pt(2)*log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1))), &
                  all((ii>=li).and.(ii<=ui)))
                filled(ix,iy,iz,iv) = source%level
            end do
          end do
        end do
      end if
    else
      if (all_cells_l) then
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              !!!if (filled(ix,iy,iz,iv) > source%level) cycle
              target%mem(ix,iy,iz,iv,target%it,1) =  &
                pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
                pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)
              filled(ix,iy,iz,iv) = source%level
            end do
          end do
        end do
      else
        do iz=l(3),u(3)
          do iy=l(2),u(2)
            do ix=l(1),u(1)
              ii = [ix,iy,iz]
              ! prefer same level for filling guard zones
              !!!if (filled(ix,iy,iz,iv) > source%level) cycle
              target%mem(ix,iy,iz,iv,target%it,1) = &
                merge(target%mem(ix,iy,iz,iv,target%it,1), &
                  pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
                  pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1), &
                  all((ii>=li).and.(ii<=ui)))
              filled(ix,iy,iz,iv) = source%level
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

!===============================================================================
!> Same as above, but for RT_solver type patches
!===============================================================================
SUBROUTINE same_rt (target, source, filled)
  class(patch_t), pointer :: target, source
  integer(int8), dimension(:,:,:,:):: filled
  !.............................................................................
  integer                 :: l(3), u(3), o(3), ii(3), jt(2), li(3), ui(3)
  real(8)                 :: position(3)
  real                    :: pt(2)
  integer                 :: ix, iy, iz, iv
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%same_rt', itimer=itimer)
  call source%time_interval (target%time, jt, pt)
  call dnload_range (target, source, 1, l, u, o)
  li = target%mesh%li
  ui = target%mesh%ui
  do iv=1,target%nv
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii = [ix,iy,iz]
          if (filled(ix,iy,iz,iv) > source%level) cycle
          if (all((ii>=li).and.(ii<=ui))) cycle
          target%mem(ix,iy,iz,iv,target%it,1) = &
             pt(1)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1) + &
             pt(2)*source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)
          filled(ix,iy,iz,iv) = source%level
        end do
      end do
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE same_rt

!===============================================================================
!> Interpolate guard zone info from source patch at a lower level
!===============================================================================
SUBROUTINE different_prolong (target, source, filled, all_cells, only)
  class(patch_t), pointer :: target, source
  integer(int8), dimension(:,:,:,:):: filled
  logical, optional:: all_cells
  integer, optional:: only
  !.............................................................................
  integer                 :: ix, iy, iz, iv, jt(2), iv1, iv2, i
  integer, dimension(3)   :: l, u, ii, li, ui, jj
  real(8), dimension(3)   :: pos
  logical                 :: valid
  integer, save           :: itimer=0
  real                    :: pt(2), pp(3)
  real(kind=KindScalarVar):: qint
  logical                 :: all_cells_l
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different_prolong', itimer=itimer)
  if (io%verbose > 1) &
    print *,'MK different_prolong; target:', target%id, target%level, &
                                  'source:', source%id, source%level
  call source%time_interval (target%time, jt, pt)
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
  li = target%mesh%li
  ui = target%mesh%ui
  do iv=iv1,iv2
    if (iv == target%idx%dphi) cycle ! `dphi` is part of `mem` when `time_order` = 0; do not download it.
    call dnload_range (target, source, target%idx%d, l, u)
    !---------------------------------------------------------------------------
    ! For no-no-mans-land meshes, use guard zone values for the li index, to
    ! make fluxes up/down symmetric
    !---------------------------------------------------------------------------
!    do i=1,3
!      if (target%mesh(i)%h(iv) /= 0.0) li(i) = li(i)+1
!    end do
    if (target%id==io%id_debug) &
      print 1, source%id, (source%position-target%position)/target%size, l, u
      1 format("DBG download_t%different_prolong:", i6,"  displ:",3f11.6,"  l,u:",2(3i4,2x))
    associate (q1 => source%mem(:,:,:,:,jt(1),1), q2 => source%mem(:,:,:,:,jt(2),1))
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii=[ix,iy,iz]
          if (all_cells_l .or. any((ii-li) < 0 .or. (ii-ui) > 0)) then
            if (filled(ix,iy,iz,iv) > source%level) then
              if (target%id==io%id_debug .and. io%verbose>1) &
                print *, target%id, 'DBG different_prolong: skipping', iv,ii, filled(ix,iy,iz,iv)
              cycle
            end if
            if (target%id==io%id_debug .and. io%verbose>1) &
              print *, target%id, 'DBG different_prolong:    using', iv,ii, present(all_cells), (ii-li) < 0
            pos = target%myposition (ii, iv)
            call source%index_space (pos, iv, jj, pp)
            if (source%pervolume(iv)) then
              qint = interpolator%interpolator_specific (q1, jj(1), jj(2), jj(3), iv, pp) * pt(1) &
                   + interpolator%interpolator_specific (q2, jj(1), jj(2), jj(3), iv, pp) * pt(2)
            else
              qint = interpolator%interpolator (q1, jj(1), jj(2), jj(3), iv, pp) * pt(1) &
                   + interpolator%interpolator (q2, jj(1), jj(2), jj(3), iv, pp) * pt(2)
            end if
            if (all((ii>li.and.ii<ui)) .and. .not. all_cells_l) then
              print *,'INTERNAL BREACH!', ii
            end if
            target%mem(ix,iy,iz,iv,target%it,1) = qint
            filled(ix,iy,iz,iv) = source%level
          else
            if (target%id==io%id_debug .and. io%verbose>1) &
              print *, target%id, 'DBG different_prolong: skipping', iv,ii, present(all_cells), (ii-li) < 0
            cycle
          end if
        end do
      end do
    end do
    end associate
  end do
  if (io%verbose>2.or.target%id==io%id_debug) &
    print 2, target%id, target%level, source%id, source%level, 'different_prolong picked up', &
    product(u-l+1), ' values, mind =', &
    minval(target%mem(l(1):u(1),l(2):u(2),l(3):u(3),1,target%it,1)), l, u
  2 format(2(i6,i3),2x,a,i6,a,1p,e12.4,2(2x,3i4))
  call trace%end (itimer)
END SUBROUTINE different_prolong

!===============================================================================
!> Same as above, but for RT_solver type patches
!===============================================================================
SUBROUTINE different_rt (target, source, filled)
  class(patch_t), pointer :: target, source
  integer(int8), dimension(:,:,:,:):: filled
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
  call source%time_interval (target%time, jt, pt)
  call dnload_range (target, source, 1, l, u)
  li = target%mesh%li
  ui = target%mesh%ui
  do iv=1,target%nv
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii=[ix,iy,iz]
          if (filled(ix,iy,iz,iv) > source%level) cycle
          if (all((ii>=li).and.(ii<=ui))) cycle
          call source%interpolate (target, ii, iv, jt, pt)
          filled(ix,iy,iz,iv) = source%level
        end do
      end do
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE different_rt

!===============================================================================
!> Restrict guard zone info from source patch at a higher level
!===============================================================================
SUBROUTINE different_restrict (target, source, filled, only)
  class(patch_t), pointer :: target, source
  integer(int8), dimension(:,:,:,:):: filled
  integer, optional:: only
  integer                 :: ix, iy, iz, iv, jt(2), iv1, iv2, i, jx, jy, jz
  integer, dimension(3)   :: l, u, ii, li, ui, jj
  real(8), dimension(3)   :: target_pos, source_pos, target_upper_edges, target_lower_edges, &
                             source_upper_edges, source_lower_edges, le, ue
  logical                 :: valid
  integer, save           :: itimer=0
  real                    :: pt(2), pp(3)
  real(kind=KindScalarVar):: qint
  logical                 :: qmask(3)
  real(8):: rratio(3), dVsource, dVtarget, dWsource, dWtarget, filling_factor, accumulate, &
            filling_total, dWtargeti
  integer:: ilook(3)
  !-----------------------------------------------------------------------------
  call trace%begin('download_t%different_restrict', itimer=itimer)

  call source%time_interval (target%time, jt, pt)
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
    if (iv == target%idx%phi) then
      if (target%poisson_time > source%poisson_time + source%dtime*source%grace) then
        cycle
      end if
    end if
    if (iv == target%idx%dphi) cycle ! `dphi` is part of `mem` when `time_order` = 0; do not download it.
    qmask(:) = .true. ! logical mask to track if variable is averaged over a volume or an area
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
    associate (q1 => source%mem(:,:,:,iv,jt(1),1), q2 => source%mem(:,:,:,iv,jt(2),1))
    do iz=l(3),u(3)
      do iy=l(2),u(2)
        do ix=l(1),u(1)
          ii = [ix,iy,iz]
          ! prefer same level for guard zone filling
          if (filled(ix,iy,iz,iv) == target%level .and. &
              (any(ii < target%li) .or. any(ii > target%ui))) cycle
          target_pos = target%myposition (ii, iv)
          target_lower_edges = target_pos - target%ds * 0.499999_8 ! to 1e-6 tolerance
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
                source_lower_edges = source_pos - source%ds * 0.499999_8 ! to 1e-6 tolerance
                source_upper_edges = source_pos + source%ds * 0.499999_8
                ue = min(source_upper_edges, target_upper_edges)
                le = max(source_lower_edges, target_lower_edges)
                ! filling factor; if zero, then no overlap.
                ! the logical mask is used to produce an area filling factor for
                ! mag. fields (but a volume filling factor otherwise)
                filling_factor = product(max(ue - le, 0.0),mask=qmask) ! to ~1.2e-5 tolerance
                filling_total = filling_total + filling_factor
                ! one might note that `accumulate` should be multipled by `dWsource`.
                ! but `filling_factor` should be divided by `dWsource`, so they cancel.
                if (filling_factor > 0.0) then
                  accumulate = accumulate &
                             + (q1(jx,jy,jz) * pt(1) + q2(jx,jy,jz) * pt(2)) * filling_factor
                end if
              end do
            end do
          end do
          ! replace a portion of the current cell with values from the finer patch
          if (filling_total > 0.0) then
            target%mem(ix,iy,iz,iv,target%it,1) = accumulate * dWtargeti &
                                                + target%mem(ix,iy,iz,iv,target%it,1) &
                                                  * (1.0 - filling_total * dWtargeti)
            filled(ix,iy,iz,iv) = source%level
          end if
        end do
      end do
    end do
    end associate
  end do
  call trace%end (itimer)
END SUBROUTINE different_restrict

!===============================================================================
!> Create a nbors list sorted by decreasing distance.
!===============================================================================
SUBROUTINE sort_nbors_by_distance (this, head)
  class(link_t), pointer:: this, head
  class(link_t), pointer:: nbor, xx, pp, nn
  real(8):: r2n, r2x
  !.............................................................................
  allocate (head)
  if (.not. associated(this%nbor)) then
    nullify(head)
    return
  end if
  head%task => this%nbor%task
  nbor => this%nbor%next

  do while (associated(nbor))
    allocate (nn)
    nn%task => nbor%task
    r2n = sum(this%task%distance(nbor%task)**2)
    xx => head
    nullify(pp)
    do while (associated(xx))
      if (.not.associated(xx%task,nbor%task)) then
        r2x = sum(this%task%distance(xx%task)**2)
        if (r2n > r2x) then
          ! insert before
          nn%next => xx
          if (.not. associated(pp)) then
            head => nn
          else
            pp%next => nn
          end if
          exit
        end if
        if (.not. associated(xx%next)) xx%next => nn
      end if
      pp => xx
      xx => xx%next
    end do
    nbor => nbor%next
  end do

  if (io%verbose > 1) then
    print*, "target: ",this%task%id
    print*, "sorted nbors:"
    xx => head
    do while(associated(xx))
      print*,xx%task%id,sum(this%task%distance(xx%task)**2)
      xx => xx%next
    end do
    print*,'----'
  end if
END SUBROUTINE sort_nbors_by_distance

!===============================================================================
SUBROUTINE clean_sorted_nbor_list (head)
  class(link_t), pointer:: head
  class(link_t), pointer:: next, prev
  !.............................................................................
  next => head%next
  deallocate (head)
  do while(associated(next))
    prev => next
    next => next%next
    deallocate (prev)
  end do
END SUBROUTINE clean_sorted_nbor_list

END MODULE
