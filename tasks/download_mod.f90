!> =============================================================================
!> =============================================================================
MODULE download_mod
  USE io_mod
  USE omp_mod
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
END SUBROUTINE download_link

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
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: patch, other
  class(task_t), pointer:: task
  integer:: ix, iy, iz, o, jt(2), it, iv, nt, nv, neg
  real:: pt(2)
  type(mesh_t), pointer:: m(:)
  real, dimension(:,:,:,:), pointer:: guard
  logical:: restore(-3:3)
  !-----------------------------------------------------------------------------
  call trace_begin('download::download_same_static')
  patch => task2patch (link%task)
  m => patch%mesh
  it = patch%it
  nv = patch%nv
  !-----------------------------------------------------------------------------
  ! Save guard zone values and restore the planes that have no neighbors; they
  ! serve as static boundary conditions, but may be overwritten in the corners,
  ! since we only download from 6 neighbors.
  !-----------------------------------------------------------------------------
  allocate (guard(m(1)%gn,m(2)%gn,m(3)%gn,nv))
  restore = .false.
  guard = patch%mem(:,:,:,:,it,1)
  !-----------------------------------------------------------------------------
  ! XY-planes in the two Z guard zones
  !-----------------------------------------------------------------------------
  if (m(3)%n > 1) then
    other => find_nbor (link, -3, o, patch%time, jt, pt)
    if (associated(other)) then
      neg = 0
      call other%lock%set ('download_same -3')
      do iz=m(3)%lb,m(3)%lo
        do iy=m(2)%lb,m(2)%ub
          do iv=1,nv
            do ix=m(1)%lb,m(1)%ub
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy,iz+o,iv,jt(1),1)*pt(1) &
                                          + other%mem(ix,iy,iz+o,iv,jt(2),1)*pt(2)
              neg = neg+merge(1,0,iv==1.and.patch%mem(ix,iy,iz,iv,it,1)<=0.0)
            end do
          end do
        end do
      end do
      call other%lock%unset ('download_same -3')
      if (neg>0) print*, other%id, 'MK: -3, neg, o =', neg, o, jt, pt
      call patch%check_density ('after -3')
    else
      restore(-3) = .true.
          guard(1:m(1)%gn,1:m(2)%gn,m(3)%lb:m(3)%lo,:) = &
      patch%mem(1:m(1)%gn,1:m(2)%gn,m(3)%lb:m(3)%lo,:,it,1)
    end if
    other => find_nbor (link, +3, o, patch%time, jt, pt)
    if (associated(other)) then
      call other%lock%set ('download_same +3')
      neg = 0
      do iz=m(3)%uo,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%ub
            do iv=1,nv
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy,iz+o,iv,jt(1),1)*pt(1) &
                                          + other%mem(ix,iy,iz+o,iv,jt(2),1)*pt(2)
              neg = neg+merge(1,0,iv==1.and.patch%mem(ix,iy,iz,iv,it,1)<=0.0)
            end do
          end do
        end do
      end do
      if (neg>0) print*, other%id, 'MK: +3, neg, o =', neg, o, jt, pt
      call other%lock%unset ('download_same +3')
      call patch%check_density ('after +3')
    else
      restore(+3) = .true.
          guard(1:m(1)%gn,1:m(2)%gn,m(3)%uo:m(3)%ub,:) = &
      patch%mem(1:m(1)%gn,1:m(2)%gn,m(3)%uo:m(3)%ub,:,it,1)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! XZ-planes in the two Y guard zones
  !-----------------------------------------------------------------------------
  if (m(2)%n > 1) then
    other => find_nbor (link, -2, o, patch%time, jt, pt)
    if (associated(other)) then
      call other%lock%set ('download_same -2')
      neg = 0
      do iz=m(3)%lb,m(3)%ub
        do iy=m(2)%lb,m(2)%lo
          do iv=1,nv
            do ix=m(1)%lb,m(1)%ub
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy+o,iz,iv,jt(1),1)*pt(1) &
                                          + other%mem(ix,iy+o,iz,iv,jt(2),1)*pt(2)
              neg = neg+merge(1,0,iv==1.and.patch%mem(ix,iy,iz,iv,it,1)<=0.0)
            end do
          end do
        end do
      end do
      call other%lock%unset ('download_same -2')
      if (neg>0) print*, other%id, 'MK: -2, neg, o =', neg, o, jt, pt
      call patch%check_density ('after -2')
    else
      restore(-2) = .true.
          guard(1:m(1)%gn,m(2)%lb:m(2)%lo,1:m(3)%gn,:) = &
      patch%mem(1:m(1)%gn,m(2)%lb:m(2)%lo,1:m(3)%gn,:,it,1)
    end if
    other => find_nbor (link, +2, o, patch%time, jt, pt)
    if (associated(other)) then
      call other%lock%set ('download_same +2')
      neg = 0
      do iz=m(3)%lb,m(3)%ub
        do iy=m(2)%uo,m(2)%ub
          do iv=1,nv
            do ix=m(1)%lb,m(1)%ub
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy+o,iz,iv,jt(1),1)*pt(1) &
                                          + other%mem(ix,iy+o,iz,iv,jt(2),1)*pt(2)
              neg = neg+merge(1,0,iv==1.and.patch%mem(ix,iy,iz,iv,it,1)<=0.0)
            end do
          end do
        end do
      end do
      call other%lock%unset ('download_same +2')
      if (neg>0) print*, other%id, 'MK: +2, neg, o =', neg, o, jt, pt
      call patch%check_density ('after +2')
    else
      restore(+2) = .true.
          guard(1:m(1)%gn,m(2)%uo:m(2)%ub,1:m(3)%gn,:) = &
      patch%mem(1:m(1)%gn,m(2)%uo:m(2)%ub,1:m(3)%gn,:,it,1)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! YZ-planes in the two X guard zones
  !-----------------------------------------------------------------------------
  if (m(1)%n > 1) then
    other => find_nbor (link, -1, o, patch%time, jt, pt)
    if (associated(other)) then
      call other%lock%set ('download_same -1')
      neg = 0
      do iz=m(3)%lb,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%lb,m(1)%lo
            do iv=1,nv
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix+o,iy,iz,iv,jt(1),1)*pt(1) &
                                          + other%mem(ix+o,iy,iz,iv,jt(2),1)*pt(2)
              neg = neg+merge(1,0,iv==1.and.patch%mem(ix,iy,iz,iv,it,1)<=0.0)
            end do
          end do
        end do
      end do
      call other%lock%unset ('download_same -1')
      if (neg>0) print*, other%id, 'MK: -1, neg, o =', neg, o, jt, pt
      call patch%check_density ('after -1')
    else
      restore(-1) = .true.
          guard(m(1)%lb:m(1)%lo,1:m(2)%gn,1:m(3)%gn,:) = &
      patch%mem(m(1)%lb:m(1)%lo,1:m(2)%gn,1:m(3)%gn,:,it,1)
    end if
    other => find_nbor (link, +1, o, patch%time, jt, pt)
    if (associated(other)) then
      call other%lock%set ('download_same +1')
      neg = 0
      do iz=m(3)%lb,m(3)%ub
        do iy=m(2)%lb,m(2)%ub
          do ix=m(1)%uo,m(1)%ub
            do iv=1,nv
              patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix+o,iy,iz,iv,jt(1),1)*pt(1) &
                                          + other%mem(ix+o,iy,iz,iv,jt(2),1)*pt(2)
              neg = neg+merge(1,0,iv==1.and.patch%mem(ix,iy,iz,iv,it,1)<=0.0)
            end do
          end do
        end do
      end do
      call other%lock%unset ('download_same +1')
      if (neg>0) print*, other%id, 'MK: +1, neg, o =', neg, o, jt, pt
      call patch%check_density ('after +1')
    else
      restore(+1) = .true.
          guard(m(1)%uo:m(1)%ub,1:m(2)%gn,1:m(3)%gn,:) = &
      patch%mem(m(1)%uo:m(1)%ub,1:m(2)%gn,1:m(3)%gn,:,it,1)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Restore the guard zones that may have been affected by corner downloads
  !-----------------------------------------------------------------------------
  if (restore(-3)) patch%mem(1:m(1)%gn,1:m(2)%gn,m(3)%lb:m(3)%lo,:,it,1) = &
                       guard(1:m(1)%gn,1:m(2)%gn,m(3)%lb:m(3)%lo,:)
  if (restore(+3)) patch%mem(1:m(1)%gn,1:m(2)%gn,m(3)%uo:m(3)%ub,:,it,1) = &
                       guard(1:m(1)%gn,1:m(2)%gn,m(3)%uo:m(3)%ub,:)
  if (restore(-2)) patch%mem(1:m(1)%gn,m(2)%lb:m(2)%lo,1:m(3)%gn,:,it,1) = &
                       guard(1:m(1)%gn,m(2)%lb:m(2)%lo,1:m(3)%gn,:)
  if (restore(+2)) patch%mem(1:m(1)%gn,m(2)%uo:m(2)%ub,1:m(3)%gn,:,it,1) = &
                       guard(1:m(1)%gn,m(2)%uo:m(2)%ub,1:m(3)%gn,:)
  if (restore(-1)) patch%mem(m(1)%lb:m(1)%lo,1:m(2)%gn,1:m(3)%gn,:,it,1) = &
                       guard(m(1)%lb:m(1)%lo,1:m(2)%gn,1:m(3)%gn,:)
  if (restore(+1)) patch%mem(m(1)%uo:m(1)%ub,1:m(2)%gn,1:m(3)%gn,:,it,1) = &
                       guard(m(1)%uo:m(1)%ub,1:m(2)%gn,1:m(3)%gn,:)
  deallocate (guard)
  call trace_end
END SUBROUTINE download_same_static

!> =============================================================================
!> Scan the available times, and set jt(1),jt(2) to the pair that gives the
!> smallest interpolation / extrapolation weight pt, without using the iit()
!> array, and taking care to access the other%t times in atomic read statements
!> =============================================================================
SUBROUTINE find_time (other, time, jt, pt)
  class(patch_t), pointer:: other
  real(8):: time
  integer:: jt(2)
  real:: pt(2)
  !.............................................................................
  integer:: it1, it2, it
  real(8):: t1, t2
  real:: dt, pt2
  logical:: found
  !-----------------------------------------------------------------------------
  call other%lock%set ('find_time')
  jt = 1
  found = .false.
  if (time==0.0) then
    found = .true.
    pt(2) = 0.0
  else
    pt(2) = 1+other%grace
    !---------------------------------------------------------------------------
    ! Check all nt-2 time intervals, and settle for the one that gives the
    ! smallest extrapolation factor pt(2).  Initially, when all time slots
    ! have time=0.0,
    !---------------------------------------------------------------------------
    do it=1,other%nt-2
      it1 = other%iit(it)
      it2 = other%iit(it+1)
      t1 = other%t(it1)
      t2 = other%t(it2)
      dt = t2 - t1
      if (dt > 0.0) then
        pt2 = (time-other%t(it1))/dt
        if (pt2 < pt(2)) then
          found = .true.
          pt(2) = pt2
          jt(1) = it1
          jt(2) = it2
          !print*,omp_mythread(),' MK1: id, it', other%id, it, jt, dt, pt2
        else
          !print*,omp_mythread(),' MK2: id, it', other%id, it, jt, dt, pt2
        end if
      end if
    end do
  end if
  pt(1) = 1.0-pt(2)
  if (.not.found) then
    print '(i6,i4,2x,a,2i4,2f9.4,3x,"time",f9.4,3x,"other%t",7f9.5)', other%id, &
      omp_mythread(), &
      'WARNING: using fall-back', jt, pt, time, other%t(other%iit(1:other%nt-1))
  end if
  if (abs(pt(1))< 1e-4) then
    if (io%verbose>2) print*,'EXACT'
    pt(1)=0.0; pt(2)=1.0
  else if (abs(pt(2)) < 1e-4) then
    if (io%verbose>2) print*,'EXACT'
    pt(1)=1.0; pt(2)=0.0
  end if
  if (io%verbose>1) &
    print'(a,i6,f9.4,2i4,2f7.3,10f9.4)', 'find_time: ', other%id, time, jt, pt, other%t
  call other%lock%unset ('find_time')
END SUBROUTINE find_time

!> =============================================================================
!> Find the nbor task in direction idir (plus/minos)
!> =============================================================================
FUNCTION find_nbor (link, idir, o, time, jt, pt) RESULT (other)
  class(link_t), pointer:: link, nbor
  class(patch_t), pointer:: patch, other, try
  type(mesh_t), pointer:: m(:)
  real(8):: time, position(3), dist
  integer:: idir, adir, o, jt(2), it, iv, nv, nt
  real:: pt(2), dt
  !-----------------------------------------------------------------------------
  call trace_begin('download::find_nbor')
  patch => task2patch (link%task)
  !print*, patch%id, 'looking at nbors'
  m => patch%mesh
  nt = patch%nt
  nv = patch%nv
  adir = abs(idir)
  position = patch%position
  position(adir) = position(adir) + patch%size(adir)*idir/adir
  nbor => link%nbor
  nullify (other)
  do while (associated(nbor))
    try => task2patch(nbor%task)
    dist = patch%point_distance(try%position, position)
    if (io%verbose>2) &
      print'(a,i3,2i5,f8.2,2(2x,3f8.3))', 'find_nbor:', &
        idir, patch%id, nbor%task%id, dist/patch%size(adir), (position-try%position)/patch%size
    if (dist < patch%size(adir)*0.1) then
      other => try
      !print '(i6,a,f9.5,3x,"other%t:",7f9.5)', patch%id, ' found nbor, times =', patch%time, other%t(1:nt-1)
      exit
    end if
    nbor => nbor%next
  end do
  if (.not.associated(other)) then
    if (io%verbose>1) &
      print*, patch%id, 'find_nbors: not found in direction', idir
    call trace_end
    return
  end if
  !-----------------------------------------------------------------------------
  ! The offset maps, when added to a target index uo, the source lb
  !-----------------------------------------------------------------------------
  o = (other%mesh(adir)%li-patch%mesh(adir)%uo)*idir/abs(idir)
  !-----------------------------------------------------------------------------
  ! Find the time interval to interpolate in / extrapolate from
  !-----------------------------------------------------------------------------
  call find_time (other, time, jt, pt)
  if (io%verbose>3) &
    print '(a,i3,i7,i4,i4,i3,2f6.2,2x,3i3)', 'download::find_nbor: idir, id, o, jt, pt', &
      idir, other%id, o, jt, pt, nint((other%position-patch%position)/patch%size)
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
  real(8):: time
  real:: pt(2)
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
        if (patch%mem(ix,iy,iz,1,patch%it,1) <= 0.0) then
          call find_best_nbor  (link, ii, other, jt, pt, -3, time, 'XY1')
          print 2, patch%id, ii, other%id, patch%mem(ix,iy,iz,1,patch%it,1), 'XY1'
        end if
        2 format(40x,i8,3i4,i8,1p,e11.2,3x,a)
      end do
    end do
  end do
  do iz=m(3)%uo,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ii = [ix,iy,iz]
        call find_best_nbor  (link, ii, other, jt, pt, +3, time)
        call interpolate_to (patch, ii, other, jt, pt, +3)
        if (patch%mem(ix,iy,iz,1,patch%it,1) <= 0.0) then
          call find_best_nbor  (link, ii, other, jt, pt, +3, time, 'XY2')
          print 2, patch%id, ii, other%id, patch%mem(ix,iy,iz,1,patch%it,1), 'XY2'
        end if
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
        if (patch%mem(ix,iy,iz,1,patch%it,1) <= 0.0) then
          call find_best_nbor  (link, ii, other, jt, pt, -2, time, 'XZ1')
          print 2, patch%id, ii, other%id, patch%mem(ix,iy,iz,1,patch%it,1), 'XZ1'
        end if
      end do
    end do
  end do
  do iz=m(3)%li,m(3)%ui
    do iy=m(2)%uo,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ii = [ix,iy,iz]
        call find_best_nbor  (link, ii, other, jt, pt, +2, time)
        call interpolate_to (patch, ii, other, jt, pt, +2)
        if (patch%mem(ix,iy,iz,1,patch%it,1) <= 0.0) then
          call find_best_nbor  (link, ii, other, jt, pt, +2, time, 'XZ2')
          print 2, patch%id, ii, other%id, patch%mem(ix,iy,iz,1,patch%it,1), 'XZ2'
        end if
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
        if (patch%mem(ix,iy,iz,1,patch%it,1) <= 0.0) then
          call find_best_nbor  (link, ii, other, jt, pt, -1, time, 'YZ1')
          print 2, patch%id, ii, other%id, patch%mem(ix,iy,iz,1,patch%it,1)
        end if
      end do
    end do
  end do
  do iz=m(3)%li,m(3)%ui
    do ix=m(1)%uo,m(1)%ub
      do iy=m(2)%li,m(2)%ui
        ii = [ix,iy,iz]
        call find_best_nbor  (link, ii, other, jt, pt, +1, time)
        call interpolate_to (patch, ii, other, jt, pt, +1)
        if (patch%mem(ix,iy,iz,1,patch%it,1) <= 0.0) then
          call find_best_nbor  (link, ii, other, jt, pt, +1, time, 'YZ2')
          print 2, patch%id, ii, other%id, patch%mem(ix,iy,iz,1,patch%it,1)
        end if
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
  integer:: ii(3), jj(3), jt(2), iv, i, it, j, nmax, idir, adir
  real(8):: pp(3), pos(3)
  real:: pt(2)
  real:: p(4), q(4)
  !-----------------------------------------------------------------------------
  if (.not.associated(source)) then
    if (io%verbose>1) print*,'WARNING, no source in interpolate_to!'
    return
  end if
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
        if (io%verbose > 1) &
        print '(a,2i3,2f6.1,2(2x,a,i6,i3," i:",i4,f13.8))', 'INDEX ERROR:', idir, i, &
          target%mesh(i)%h(iv), source%mesh(i)%h(iv), &
          'target: id,lv=',target%id, target%level, ii(i), (source%position(i)-target%position(i))/source%mesh(i)%d, &
          'source: id,lv=',source%id, source%level, jj(i), pp(i)
        if (io%verbose > 2) then
          print*,'   source shift:', (target%position-source%position)/source%ds + source%offset
          print*,'    target size:', target%size
          print*,'   target index:', (pos(i)-target%position(i))/target%ds(i) + target%offset(i)
          print*,'    source size:', source%size
          print*,'   source index:', (pos(i)-source%position(i))/source%ds(i) + source%offset(i)
          print*,'staggered index:', (pos(i)-source%position(i))/source%ds(i) + source%offset(i) - source%mesh(i)%h(iv)
        end if
        pp = max(-0.35d0,min(pp,1.35d0))
      end if
    end do
    adir = abs(idir)
    p(1:3)=pp; p(4)=pt(2)
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
SUBROUTINE find_best_nbor (link, ii, best, jt, pt, idir, time, debug)
  class(link_t), pointer:: link
  class(patch_t), pointer:: best
  integer:: ii(3), jt(2), idir
  real:: pt(2)
  real(8):: time
  character(len=*), optional:: debug
  !.............................................................................
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: other, patch
  type(mesh_t), pointer:: m(:)
  real(8):: dt, position(3), dist, pp(3)
  real:: score, best_score
  integer:: i, it, iv, nv, nt, adir, best_id
  !-----------------------------------------------------------------------------
  call trace_begin('download_mod::find_best_nbor',2)
  if (io%do_debug) print*,'find_best_nbor:', link%task%id, ii, time
  patch => task2patch (link%task)
  m => patch%mesh
  nt = patch%nt
  nv = patch%nv
  adir = abs(idir)
  !
  position = patch%myposition (ii)
  !if (associated(best)) then
  !  if (all(abs(position-best%position) < best%size*0.5d0)) then
  !    call trace_end
  !    return
  !  end if
  !end if
  nullify(best)
  best_id = 0
  best_score = 0.0
  nbor => link%nbor
  do while (associated(nbor))
    !---------------------------------------------------------------------------
    ! Skíp patches more than one level up/down, and the parent patch
    !---------------------------------------------------------------------------
    other => task2patch(nbor%task)
    if (abs(patch%level-other%level)<2) then
      pp = (position-other%position)
      !-----------------------------------------------------------------------
      ! Design a score that is negative if any coordinate is outside the guard
      ! zone, and which counts +3 for an inside primary coordinate, and +1 for
      ! an inside secondary coordinate.  Levels count less, but enough to favor
      ! the higher level against a lower level
      !-----------------------------------------------------------------------
      score = other%level*0.25
      do i=1,3
        score = score + merge(3.0,1.0,i==adir)* &
          merge(1.0,0.0,abs(pp(i)) < other%size(i)*0.5d0) - &
          merge(99.,0.0,abs(pp(i)) > other%gsize(i)*0.5d0)
      end do
      if (score > best_score) then
          best => other
          best_score = score
          best_id = best%id
      end if
      if (io%verbose>1.or.present(debug)) then
        block
          logical dbg
          dbg = (abs(idir)==1 .and. ii(2)==m(2)%gn/2 .and. ii(3)==m(3)%gn/2) &
           .or. (abs(idir)==2 .and. ii(1)==m(1)%gn/2 .and. ii(3)==m(3)%gn/2) &
           .or. (abs(idir)==3 .and. ii(1)==m(1)%gn/2 .and. ii(2)==m(2)%gn/2)
          if (dbg) then
            if (associated(best)) then
              print 1,'patch ',patch%id, patch%level, idir, &
                      'nbor:', other%id, other%level, score, &
                      'best:',  best%id,  best%level, best_score
            1 format(a,i6,i4,i3,2(4x,a,i6,i4,f10.4))
            else
              print 1,'patch ',patch%id, patch%level, idir, &
                '   nbor:', other%id, other%level, score
            end if
          end if
        end block
      end if
    end if
    nbor => nbor%next
  end do
  !-----------------------------------------------------------------------------
  ! Find the right time interval in the other patch
  !-----------------------------------------------------------------------------
  if (associated(best)) &
    call find_time (best, time, jt, pt)
  !
  call trace_end
END SUBROUTINE find_best_nbor

END MODULE
