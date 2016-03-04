!> =============================================================================
!> =============================================================================
MODULE download_mod
  USE io_mod
  USE bits_mod
  USE mesh_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE list_mod
  implicit none
  private
  public download_link
CONTAINS

!> =============================================================================
!> Select the best method for each case.  So far, this only supports static
!> patches at the same level, or static patches one level up or down
!> =============================================================================
SUBROUTINE download_link (link, all_cells)
  class(link_t), pointer:: link
  logical, optional:: all_cells
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: patch
  logical:: static, same
  integer:: nt, nv
  !-----------------------------------------------------------------------------
  call trace_begin('download_t%download')
  patch => task2patch(link%task)
  if (patch%is_set(bits%root_grid)) then
    call patch%periodic_grid
  else
    nt = patch%nt
    nv = patch%nv
    static = .true.; same = .true.
    if (patch%is_clear(bits%static)) static = .false.
    nbor => link%nbor
    do while (associated(nbor))
      if (io%verbose>4) &
        print*,patch%id,nbor%task%id,nbor%task%is_clear(bits%static),nbor%task%level
      if (nbor%task%is_clear(bits%static)) static = .false.
      if (nbor%task%level /= patch%level) same = .false.
      if (io%do_debug) print*,patch%id,patch%level,nbor%task%id,nbor%task%level,same
      nbor => nbor%next
    end do
    if (static) then
      if (same) then
        call download_same_static (link)
      else
        call download_different_static (link)
      end if
    else
      print *, 'ERROR: moving patches found'
    end if
  end if
  call trace_end
END SUBROUTINE

!> =============================================================================
!> We speed up the simple case when nbor patches have the same resolution, and
!> are static, in several ways:
!>   1) linear interpolation in time, with fixed weights
!>   2) vectorizing over the x-direction, using as long vectors as possible
!>   3) looping instead over iv when the x-loop is too short
!>   4) using the "end"-patches also for edge- and corner- patches -- this is
!>      ok, since we only need these values to be sensible. This way, we only
!>      need 6 of the 26 neighbor patches in 3D!
!> =============================================================================
SUBROUTINE download_same_static (link)
  class(patch_t), pointer:: patch, other
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task
  integer:: ix, iy, iz, o, jt1, jt2, it, iv, nt, nv
  real:: pt, qt
  type(mesh_t), pointer:: m(:)
  !-----------------------------------------------------------------------------
  call trace_begin('download::download_same_static')
  patch => task2patch (link%task)
  m => patch%mesh
  it = patch%it
  nv = patch%nv
  !-----------------------------------------------------------------------------
  ! XY-planes in the two Z guard zones
  !-----------------------------------------------------------------------------
  if (m(3)%n > 1) then
    other => find_nbor (link, -3, o, patch%time, jt1, jt2, pt, qt)
    if (associated(other)) then
      do iz=m(3)%lb,m(3)%lo
        do iy=m(2)%lb,m(2)%ub
          do iv=1,nv
            do ix=m(1)%lb,m(1)%ub
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy,iz+o,iv,jt1,1)*qt &
                                          + other%mem(ix,iy,iz+o,iv,jt2,1)*pt
            end do
          end do
        end do
      end do
    end if
    other => find_nbor (link, +3, o, patch%time, jt1, jt2, pt, qt)
    if (associated(other)) then
      do iz=m(3)%uo,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            do iv=1,nv
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy,iz+o,iv,jt1,1)*qt &
                                          + other%mem(ix,iy,iz+o,iv,jt2,1)*pt
            end do
          end do
        end do
      end do
    end if
  end if
  !-----------------------------------------------------------------------------
  ! XZ-planes in the two Y guard zones
  !-----------------------------------------------------------------------------
  if (m(2)%n > 1) then
    other => find_nbor (link, -2, o, patch%time, jt1, jt2, pt, qt)
    if (associated(other)) then
      do iz=m(3)%li,m(3)%ui
        do iy=m(2)%lb,m(2)%lo
          do iv=1,nv
            do ix=m(1)%lb,m(1)%ub
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy+o,iz,iv,jt1,1)*qt &
                                          + other%mem(ix,iy+o,iz,iv,jt2,1)*pt
            end do
          end do
        end do
      end do
    end if
    other => find_nbor (link, +2, o, patch%time, jt1, jt2, pt, qt)
    if (associated(other)) then
      do iz=m(3)%li,m(3)%ui
        do iy=m(2)%uo,m(2)%ub
          do iv=1,nv
            do ix=m(1)%lb,m(1)%ub
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy+o,iz,iv,jt1,1)*qt &
                                          + other%mem(ix,iy+o,iz,iv,jt2,1)*pt
            end do
          end do
        end do
      end do
    end if
  end if
  !-----------------------------------------------------------------------------
  ! YZ-planes in the two X guard zones
  !-----------------------------------------------------------------------------
  if (m(1)%n > 1) then
    other => find_nbor (link, -1, o, patch%time, jt1, jt2, pt, qt)
    if (associated(other)) then
      do iz=m(3)%li,m(3)%ui
        do iy=m(2)%li,m(2)%ui
          do ix=m(1)%lb,m(1)%lo
            do iv=1,nv
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix+o,iy,iz,iv,jt1,1)*qt &
                                          + other%mem(ix+o,iy,iz,iv,jt2,1)*pt
            end do
          end do
        end do
      end do
    end if
    other => find_nbor (link, +1, o, patch%time, jt1, jt2, pt, qt)
    if (associated(other)) then
      do iz=m(3)%li,m(3)%ui
        do iy=m(2)%li,m(2)%ui
          do ix=m(1)%uo,m(1)%ub
            do iv=1,nv
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix+o,iy,iz,iv,jt1,1)*qt &
                                          + other%mem(ix+o,iy,iz,iv,jt2,1)*pt
            end do
          end do
        end do
      end do
    end if
  end if
  call trace_end
END SUBROUTINE download_same_static

!> =============================================================================
!> Find the nbor task in direction idir (plus/minos)
!> =============================================================================
FUNCTION find_nbor (link, idir, o, time, jt1, jt2, pt, qt) RESULT (other)
  class(link_t), pointer:: link, nbor
  class(patch_t), pointer:: patch, other
  type(mesh_t), pointer:: m(:)
  real(8):: time, position(3), dist
  integer:: idir, adir, o, jt1, jt2, it, iv, nv, nt
  real:: pt, qt, dt
  !-----------------------------------------------------------------------------
  !$omp critical (list_cr)
  call trace_begin('download::find_nbor')
  patch => task2patch (link%task)
  m => patch%mesh
  nt = patch%nt
  nv = patch%nv
  adir = abs(idir)
  position = patch%position
  position(adir) = position(adir) + patch%size(adir)*idir/adir
  nbor => link%nbor
  do while (associated(nbor))
    other => task2patch(nbor%task)
    dist = patch%point_distance(other%position, position)
    if (io%verbose>5) &
      print'(a,2i5,f8.3,2(2x,3f8.3))', 'find_nbor:', &
        patch%id, nbor%task%id, dist, (position-other%position)/patch%size
    if (dist < patch%size(abs(idir))*0.5) exit
    nbor => nbor%next
  end do
  if (.not.associated(other)) then
    print*,'find_nbors: ERROR'
    stop
  end if
  !-----------------------------------------------------------------------------
  ! The offset maps, when added to a target index uo, the source lb
  !-----------------------------------------------------------------------------
  o = (other%mesh(adir)%li-patch%mesh(adir)%uo)*idir/abs(idir)
  !-----------------------------------------------------------------------------
  ! Find the right time interval in the other patch
  !-----------------------------------------------------------------------------
  jt1=1
  do it=1,nt-2
    if (other%t(other%iit(it))<time) jt1=it
    if (io%verbose>4) print*,it,time, other%t(other%iit(it))
  end do
  jt2 = jt1+1
  jt1 = other%iit(jt1)
  jt2 = other%iit(jt2)
  dt = other%t(jt2)-other%t(jt1)
  pt = (time-other%t(jt1))/max(dt,1e-30)
  qt = 1.-pt
  if (abs(pt)< 1e-4) then
    if (io%verbose>1) print*,'EXACT'
    pt=0.0; qt=1.0
  else if (abs(qt) < 1e-4) then
    if (io%verbose>1) print*,'EXACT'
    pt=1.0; qt=0.0
  end if
  if (io%verbose>2) &
    print '(a,i3,i7,i4,i4,i3,f6.2,2x,3i3)', 'download::find_nbor: idir, id, o, jt, pt', &
      idir, other%id, o, jt1, jt2, pt, nint((other%position-patch%position)/patch%size)
  !$omp end critical (list_cr)
  call trace_end
END FUNCTION find_nbor

!> =============================================================================
!> =============================================================================
SUBROUTINE download_different_static (link)
  class(patch_t), pointer:: patch, other
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task
  type:: target_t
    integer:: ix        ! index in x
    integer:: iy        ! index in y
    integer:: iz        ! index in z
    real:: time         ! target time
  end type
  type:: source_t
    integer:: jx(2)     ! indices in x
    integer:: jy(2)     ! indices in y
    integer:: jz(2)     ! indices in z
    integer:: jt(2)     ! indices in t
    real:: px(2)        ! weights in x
    real:: py(2)        ! weights in y
    real:: pz(2)        ! weights in z
    real:: pt(2)        ! weights in t
  end type
  integer:: ix, iy, iz, jx(2), jy(2), jz(2), ii(3), jt(2)
  real(8):: pt, time
  type(mesh_t), pointer:: m(:)
  !-----------------------------------------------------------------------------
  ! XY-planes in the two Z guard zones
  !-----------------------------------------------------------------------------
  call trace_begin('download::download_different_static')
  patch => task2patch (link%task)
  m => patch%mesh
  time = patch%time
  do iz=m(3)%lb,m(3)%lo
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ii = [ix,iy,iz]
        call find_best_nbor  (link, ii, other, jt, pt, -3, time)
        call interpolate_to (patch, ii, other, jt, pt, -3)
      end do
    end do
  end do
  do iz=m(3)%uo,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ii = [ix,iy,iz]
        call find_best_nbor  (link, ii, other, jt, pt, +3, time)
        call interpolate_to (patch, ii, other, jt, pt, +3)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! XZ-planes in the two Y guard zones
  !-----------------------------------------------------------------------------
  do iz=m(3)%li,m(3)%ui
    do iy=m(2)%lb,m(2)%lo
      do ix=m(1)%lb,m(1)%ub
        ii = [ix,iy,iz]
        call find_best_nbor  (link, ii, other, jt, pt, -2, time)
        call interpolate_to (patch, ii, other, jt, pt, -2)
      end do
    end do
  end do
  do iz=m(3)%li,m(3)%ui
    do iy=m(2)%uo,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        call find_best_nbor  (link, ii, other, jt, pt, +2, time)
        call interpolate_to (patch, ii, other, jt, pt, +2)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! YZ-planes in the two X guard zones
  !-----------------------------------------------------------------------------
  do iz=m(3)%li,m(3)%ui
    do ix=m(1)%lb,m(1)%lo
      do iy=m(2)%li,m(2)%ui
        ii = [ix,iy,iz]
        call find_best_nbor  (link, ii, other, jt, pt, -1, time)
        call interpolate_to (patch, ii, other, jt, pt, -1)
      end do
    end do
  end do
  do iz=m(3)%li,m(3)%ui
    do ix=m(1)%uo,m(1)%ub
      do iy=m(2)%li,m(2)%ui
        ii = [ix,iy,iz]
        call find_best_nbor  (link, ii, other, jt, pt, +1, time)
        call interpolate_to (patch, ii, other, jt, pt, +1)
      end do
    end do
  end do
  call trace_end
END SUBROUTINE download_different_static

!> =============================================================================
!> =============================================================================
SUBROUTINE interpolate_to (target, ii, source, jt, pt, idir)
  class (patch_t):: target
  class (patch_t), pointer:: source
  integer:: ii(3), jj(3), jt(2), iv, i, it, j, nmax, idir
  real(8):: pp(3), pos(3), pt
  real:: p(4), q(4)
  !-----------------------------------------------------------------------------
  if (.not.associated(source)) return
  call trace_begin ('download_mod::interpolate_to',2)
  if (io%do_debug) print*,'interpolate_to',target%it,ii,source%id,jt,pt
  it = target%it
  jj = 0
  do iv=1,source%nv
    do i=1,3
      pos(i) = target%position(i) + (ii(i)-target%offset(i)+target%mesh(i)%h(iv))*target%ds(i)      ! position of target ii
      pp(i) = (pos(i)-source%position(i))/source%ds(i) + source%offset(i) - source%mesh(i)%h(iv)    ! float index in source
      jj(i) = pp(i)                                                                                 ! integer part
      jj(i) = max(1,min(source%mesh(i)%gn-1,jj(i)))
      pp(i) = pp(i) - jj(i)                                                                         ! fractions
      if (pp(i)<-0.35d0 .or. pp(i)>1.35d0) then
        pp(i) = (pos(i)-source%position(i))/source%ds(i) + source%offset(i) - source%mesh(i)%h(iv)  ! float index in source
        print '(a,2i3,2f6.1,2(2x,a,i6,i3," i:",i4,f13.8))', 'INDEX ERROR:', idir, i, &
        target%mesh(i)%h(iv), source%mesh(i)%h(iv), &
        'target: id,lv=',target%id, target%level, ii(i), (source%position(i)-target%position(i))/source%mesh(i)%d, &
        'source: id,lv=',source%id, source%level, jj(i), pp(i)
        print*,'   source shift:', (target%position-source%position)/source%ds + source%offset
        print*,'    target size:', target%size
        print*,'   target index:', (pos(i)-target%position(i))/target%ds(i) + target%offset(i)
        print*,'    source size:', source%size
        print*,'   source index:', (pos(i)-source%position(i))/source%ds(i) + source%offset(i)
        print*,'staggered index:', (pos(i)-source%position(i))/source%ds(i) + source%offset(i) - source%mesh(i)%h(iv)
        stop
      end if
    end do
    p(1:3)=pp; p(4)=pt
    q = 1.-p
    target%mem(ii(1),ii(2),ii(3),iv,it,1) = &
      q(4)*(q(3)*(q(2)*(q(1)*source%mem(jj(1)  ,jj(2)  ,jj(3)  ,iv,jt(1),1)    + &
                        p(1)*source%mem(jj(1)+1,jj(2)  ,jj(3)  ,iv,jt(1),1))   + &
                  p(2)*(q(1)*source%mem(jj(1)  ,jj(2)+1,jj(3)  ,iv,jt(1),1)    + &
                        p(1)*source%mem(jj(1)+1,jj(2)+1,jj(3)  ,iv,jt(1),1)))  + &
            p(3)*(q(2)*(q(1)*source%mem(jj(1)  ,jj(2)  ,jj(3)+1,iv,jt(1),1)    + &
                        p(1)*source%mem(jj(1)+1,jj(2)  ,jj(3)+1,iv,jt(1),1))   + &
                  p(2)*(q(1)*source%mem(jj(1)  ,jj(2)+1,jj(3)+1,iv,jt(1),1)    + &
                        p(1)*source%mem(jj(1)+1,jj(2)+1,jj(3)+1,iv,jt(1),1)))) + &
      p(4)*(q(3)*(q(2)*(q(1)*source%mem(jj(1)  ,jj(2)  ,jj(3)  ,iv,jt(2),1)    + &
                        p(1)*source%mem(jj(1)+1,jj(2)  ,jj(3)  ,iv,jt(2),1))   + &
                  p(2)*(q(1)*source%mem(jj(1)  ,jj(2)+1,jj(3)  ,iv,jt(2),1)    + &
                        p(1)*source%mem(jj(1)+1,jj(2)+1,jj(3)  ,iv,jt(2),1)))  + &
            p(3)*(q(2)*(q(1)*source%mem(jj(1)  ,jj(2)  ,jj(3)+1,iv,jt(2),1)    + &
                        p(1)*source%mem(jj(1)+1,jj(2)  ,jj(3)+1,iv,jt(2),1))   + &
                  p(2)*(q(1)*source%mem(jj(1)  ,jj(2)+1,jj(3)+1,iv,jt(2),1)    + &
                        p(1)*source%mem(jj(1)+1,jj(2)+1,jj(3)+1,iv,jt(2),1))))
  end do
  call trace_end
END SUBROUTINE interpolate_to

!> =============================================================================
!> Find the best nbor patch overlapping the point with indices ii in link%patch
!> =============================================================================
SUBROUTINE find_best_nbor (link, ii, best, jt, pt, idir, time)
  class(link_t), pointer:: link
  class(patch_t), pointer:: best
  integer:: ii(3), jt(2), idir, adir
  real(8):: pt, dt
  !.............................................................................
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: other, patch
  type(mesh_t), pointer:: m(:)
  real(8):: time, position(3), dist, best_dist, pp(3)
  integer:: it, iv, nv, nt
  !-----------------------------------------------------------------------------
  call trace_begin('download_mod::find_best_nbor',2)
  if (io%do_debug) print*,'find_best_nbor:', link%task%id, ii, time
  patch => task2patch (link%task)
  m => patch%mesh
  nt = patch%nt
  nv = patch%nv
  !
  position = patch%myposition (ii)
  !if (associated(best)) then
  !  if (all(abs(position-best%position) < best%size*0.5d0)) then
  !    call trace_end
  !    return
  !  end if
  !end if
  nullify(best)
  best_dist = 1e30
  nbor => link%nbor
  do while (associated(nbor))
    !---------------------------------------------------------------------------
    ! Skíp patches more than one level up/down, and the parent patch
    !---------------------------------------------------------------------------
    other => task2patch(nbor%task)
    if (abs(patch%level-other%level)<2) then
      pp = (position-other%position)
      if (all(abs(pp) < other%size*0.5d0)) then
        !-----------------------------------------------------------------------
        ! Choose the "best" of the overlapping ones -- the smallest distance
        ! among those with the highest level
        !-----------------------------------------------------------------------
        if (.not.associated(best)) then
          best => other
        else if (other%level > best%level) then
          best => other
        end if
        if (io%verbose>1) &
          print 1,'patch ',patch%id, patch%level, &
            '   best nbor:', other%id, other%level, &
            (other%position-patch%position)/other%ds, &
            patch%distance(other)/other%size
          1 format(a,i6,i4,a,i6,i4,2(2x,3f10.4))
      else
          if (io%verbose>2) &
            print 1,'patch ',patch%id, patch%level, &
              ' outside nbor:', other%id, other%level, &
              (other%position-patch%position)/other%ds, &
              patch%distance(other)/other%size
      end if
    end if
    nbor => nbor%next
  end do
  !-----------------------------------------------------------------------------
  ! Find the right time interval in the other patch
  !-----------------------------------------------------------------------------
  if (associated(best)) then
    jt(1)=1
    do it=1,nt-2
      if (best%t(best%iit(it))<time) jt(1)=it
      if (io%verbose>4) print*,it,time, best%t(best%iit(it))
    end do
    jt(2) = jt(1)+1
    jt(1) = best%iit(jt(1))
    jt(2) = best%iit(jt(2))
    dt = best%t(jt(2))-best%t(jt(1))
    pt = (time-best%t(jt(1)))/max(dt,1e-30)
    if (abs(pt)< 1e-4) then
      if (io%verbose>2) print*,'EXACT'
      pt=0.0
    else if (abs(1.-pt) < 1d-4) then
      if (io%verbose>2) print*,'EXACT'
      pt=1.0
    end if
    if (io%do_debug) print*,'find_best_nbor:', patch%id, jt, pt
  else if (io%verbose>1) then
    print*,'find_best_nbor: ERROR', patch%id, ii
    print*,'   level:', patch%level, patch%size/patch%ds
    nullify(best)
    nbor => link%nbor
    do while (associated(nbor))
      other => task2patch(nbor%task)
      if (abs(patch%level-other%level)<2) then
        print*,'other%lv:', other%level
        print*,'   index:', (position-other%position)/other%ds + other%offset
        print'(2(i5,i3,3f12.4),l4)', &
          patch%id, patch%level, (patch%position-other%position)/other%size, &
          other%id, other%level, &
          other%offset+(position-other%position)/other%ds, &
          other%overlaps_point(idir,position)
        if (other%overlaps_point(idir,position)) then
          if (.not.associated(best)) then
            best => other
          else if (other%level > best%level) then
            best => other
          end if
        end if
      end if
      nbor => nbor%next
    end do
  end if
  call trace_end
END SUBROUTINE find_best_nbor

END MODULE
