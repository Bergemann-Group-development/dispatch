!> =============================================================================
!>  download_link: takes care of downloads to link%task
!>    same_static: called for static patches on the same level
!>  differ_static: called for static patches on different levels
!> interpolate_to: called for nbors on different levels
!>      find_time: computes time interpolation parameters
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
!> At the moment, this is limited to the order of the interpolator used.
!> =============================================================================
SUBROUTINE init (self)
  class(download_t):: self
  integer:: order_interpolator, iostat
  logical:: first_time=.true.
  namelist /download_params/ order_interpolator
  !-----------------------------------------------------------------------------
  call trace_begin ('download_mod%init')

  if (first_time) then
    order_interpolator = -1
    first_time = .false.
    rewind (io%input)
    read (io%input, download_params, iostat=iostat)
    if (io%master) write (*, download_params)

    self%interpolate => SelectInterpolator(order_interpolator)
  end if

  call trace_end
END SUBROUTINE

!> =============================================================================
!> Select the best method for each case.  So far, this only supports static
!> patches at the same level, or static patches one level up or down
!> =============================================================================
SUBROUTINE download_link (link, all_cells)
  class(link_t), pointer:: link
  logical, optional:: all_cells
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: patch
  logical:: static
  integer:: nt, nv
  !-----------------------------------------------------------------------------
  call trace_begin('download_t%download')
  patch => task2patch(link%task)
  if (io%verbose>2) &
    call patch%check_density ("before download")
  if (patch%is_set(bits%root_grid)) then
    call patch%periodic_grid
  else
    nt = patch%nt
    nv = patch%nv
    static = .true.
    if (patch%is_clear(bits%static)) static = .false.
    nbor => link%nbor
    do while (associated(nbor))
      if (nbor%task%is_clear(bits%no_download)) then
        if (io%master.and.io%verbose>4) &
          write (io_unit%log,*) patch%id,nbor%task%id,nbor%task%is_clear(bits%static),nbor%task%level
        if (nbor%task%is_clear(bits%static)) static = .false.
      end if
      nbor => nbor%next
    end do

    if (static) then
      call same_static (link)
      !call different_static (link)
    else
      write (io_unit%log,*) 'ERROR: moving patches found'
    end if
  end if
  call patch%check_density (" after download")
  call trace_end
END SUBROUTINE download_link

!> =============================================================================
!> Run through neighbor list and pick up guard cells. A crucial question is if
!> it is necessary to use OMP locks to pick up the guard zone values. As long as
!> the iit array is reconstucted from the single value source%it, which is atomic,
!> there should be no need for locking.  This is in patch_mod::time_interpolate.
!> =============================================================================
SUBROUTINE same_static (self)
  class (link_t):: self
  class (link_t),  pointer:: nbor
  class (patch_t), pointer:: source, target
  integer:: i, it, iv, l1(3), u1(3), o(3), jt(2), ntot, ix, iy, iz, n, verbose
  real(8):: pt(2)
  real:: dmin
  integer, save:: itimer=0
  logical:: ahead
  !.............................................................................
  call trace_begin ('download_mod::same_static', itimer=itimer)
  target => task2patch(self%task)
  n = 0
  nbor => self%nbor
  do while (associated (nbor))
    if (nbor%needed .and. nbor%task%level /= self%task%level) then
      nbor => nbor%next
      cycle
    end if
    if (nbor%task%is_clear(bits%no_download)) then
      n = n + 1
      associate (this => nbor%task)
      select type (this)
      class is (patch_t)
        source => this
        it = target%it
        call target%time_interpolate (source, jt, pt)
        call target%dnload_range (source, l1, u1, o)
        !-----------------------------------------------------------------------
        ! Ensure guard zones are filled in physical boundary regions
        !-----------------------------------------------------------------------
        where (target%mesh%lower_boundary)
          l1 = target%mesh%lb
        end where
        where (target%mesh%upper_boundary)
          u1 = target%mesh%ub
        end where
        ntot = product(u1-l1+1)
        if (io%master .and. io%verbose>3) then
          write (io_unit%log,'(a,3(3i4,3x),i6)') 'target ranges: l1, u1, o, n =', l1, u1, o, ntot
        end if
        do iv=1,source%nv
          dmin = 1e30
          do iz=l1(3),u1(3)
          do iy=l1(2),u1(2)
          do ix=l1(1),u1(1)
            target%mem(ix,iy,iz,iv,it,1) = source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)*pt(1) &
                                         + source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)*pt(2)
            dmin = min(dmin,target%mem(ix,iy,iz,iv,it,1))
          end do
          end do
          end do
          if (iv==1) then
            if (io%verbose==3 .or. dmin<0.0 .or. target%id==io%id_debug) write (io_unit%log,1) &
            'patch', target%id, ' level', target%level, ' picked up', ntot, &
            ' points from patch', source%id, ' at level', source%level, &
            'position', (source%position-target%position)/target%size, pt, dmin
          1 format(5(1x,a,i7,1x),1x,a,5f6.2,1p,e12.3)
            !-----------------------------------------------------------------------
            ! Detailed analysis of rare negative density cases, ending with a patch-up
            ! that interpolates in the log to replace the negative values.  The cases
            ! looked at generally are caused by time extrapolation, and even if they
            ! are extremely rare they need to be handled, lest the solver chokes on
            ! the unphysical values that result.
            !-----------------------------------------------------------------------
            if (dmin < 0.0) then
              write (io_unit%log,'(a,i9,3x,a,i9,3x,a,2i4,3x,a,2f7.3)') &
                'download_mod::same_static: ERROR, negative density contribution to', &
                target%id, 'from', source%id,'jt:',jt,'pt:',pt
              verbose = io_unit%verbose
              io_unit%verbose = 5
              ahead = source%is_ahead_of (target)
              call target%time_interpolate (source, jt, pt)
              io_unit%verbose = verbose
              write (io_unit%log,'(a,l2,4(i9,2f12.6),f10.3)') &
                'same_static: target,source id,time,dtime,jt1,t1,jt2,t2 =', &
                ahead, &
                target%id, target%time, target%dtime, &
                source%id, source%time, source%dtime, &
                jt(1), source%t(jt(1)), source%dt(jt(1)), &
                jt(2), source%t(jt(2)), source%dt(jt(1)), &
                (target%time-source%time)/source%dtime
              !---------------------------------------------------------------------
              ! If there were negative values produced, redo this guard zone, with
              ! debug printout and log-interpolation
              !---------------------------------------------------------------------
              do iz=l1(3),u1(3)
              do iy=l1(2),u1(2)
              do ix=l1(1),u1(1)
                target%mem(ix,iy,iz,iv,it,1) = source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1)*pt(1) &
                                             + source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1)*pt(2)
                if (target%mem(ix,iy,iz,iv,it,1) < 0.0) then
                  write (io_unit%log,'(a,1p,7g16.6)') 'same_static: source t =', &
                    (source%t(source%iit(i)), i=1,source%nt-1)
                  write (io_unit%log,'(a,1p,7g16.6)') 'same_static: source d =', &
                    (source%mem(ix+o(1),iy+o(2),iz+o(3),iv,source%iit(i),1), i=1,source%nt-1)
                  !-----------------------------------------------------------------
                  ! Extrapolate in the log, to ensure positive results
                  !-----------------------------------------------------------------
                  target%mem(ix,iy,iz,iv,it,1) = exp( &
                    log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1))*pt(1) + &
                    log(source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1))*pt(2))
                  write (io_unit%log,'(a,2(3i4,1x),1p,3g16.6)') 'same_static: ix,iy,iz,src1,src2,trg =', &
                    ix     , iy     , iz     , &
                    ix+o(1), iy+o(2), iz+o(3), &
                    source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(1),1), &
                    source%mem(ix+o(1),iy+o(2),iz+o(3),iv,jt(2),1), &
                    target%mem(ix,iy,iz,iv,it,1)
                end if
              end do
              end do
              end do
            end if
          end if
          if (io_unit%verbose>3) write (io_unit%log,2) &
            'patch', target%id, ' level', target%level, 'var', iv, ' picked up', &
            ntot, ' points from patch', source%id, ' at level', source%level, &
            'position', (source%position-target%position)/target%size, pt, dmin
          2 format(6(1x,a,i7,1x),1x,a,5f6.2,1p,e12.3)
        end do
      end select
      end associate
    end if
    nbor => nbor%next
  end do
  if (io%verbose > 2 .and. n /= 26) then
    write (io_unit%log,'(a,i4,a,i9)') &
      'download_mod::same_static: WARNING,',n,' nbors contributed to', target%id
  end if
  if (n > 0) then
    !$omp atomic
    download%n_same = download%n_same+1
  end if

  call trace_end(itimer)
END SUBROUTINE same_static

!> =============================================================================
!> Scan the available times, and set jt(1),jt(2) to the pair that gives the
!> smallest interpolation / extrapolation weight pt
!> =============================================================================
SUBROUTINE find_time (other, time, jt, pt)
  class(patch_t), pointer:: other
  real(8):: time
  integer:: jt(2)
  real:: pt(2)
  !.............................................................................
  integer:: it1, it2, it, iit(other%nt)
  real(8):: t1, t2
  real:: dt, pt2
  logical:: found
  !-----------------------------------------------------------------------------
  jt = 1
  found = .false.
  if (time==0.0) then
    found = .true.
    pt(2) = 0.0
  else
    pt(2) = 1.0001+other%grace
    !---------------------------------------------------------------------------
    ! Check all nt-2 time intervals, and settle for the one that gives the
    ! smallest extrapolation factor pt(2).
    !---------------------------------------------------------------------------
    iit = other%local_iit()
    do it=1,other%nt-2
      it1 = iit(it)
      it2 = iit(it+1)
      t1 = other%t(it1)
      t2 = other%t(it2)
      dt = t2 - t1
      if (dt > 0.0) then
        pt2 = (time-other%t(it1))/dt
        if (pt2 <= pt(2)) then
          found = .true.
          pt(2) = pt2
          jt(1) = it1
          jt(2) = it2
          !write (io_unit%log,*) omp_mythread,' MK1: id, it', other%id, it, jt, dt, pt2
        else
          !write (io_unit%log,*) omp_mythread,' MK2: id, it', other%id, it, jt, dt, pt2
        end if
      end if
    end do
  end if
  if (.not.found) then
    pt(2) = 0.0
    if (io%master .and. other%is_clear(bits%virtual) .and. io%verbose>2) &
      write (io_unit%log,'(i6,i4,2x,a,2i4,2f9.4,3x,"time",f9.4,3x,"other%t",7f9.5)') &
      other%id, omp_mythread, 'WARNING: using fall-back', &
      jt, pt, time, other%t(iit(1:other%nt-1))
  end if
  pt(1) = 1.0-pt(2)
  if (abs(pt(1))< 1e-4) then
    if (io%master .and. io%verbose>2) write (io_unit%log,*) 'EXACT'
    pt(1)=0.0; pt(2)=1.0
  else if (abs(pt(2)) < 1e-4) then
    if (io%master .and. io%verbose>2) write (io_unit%log,*) 'EXACT'
    pt(1)=1.0; pt(2)=0.0
  end if
  if (io%verbose>1) &
    write (io_unit%log,'(a,i6,f9.4,2i4,2f7.3,10f9.4)') &
    'find_time: ', other%id, time, jt, pt, other%t
END SUBROUTINE find_time

!> =============================================================================
!> =============================================================================
SUBROUTINE different_static (link)
  class(link_t), pointer:: link
  !.............................................................................
  class(patch_t), pointer:: target, source
  integer:: ix, iy, iz, jx(2), jy(2), jz(2), ii(3), jt(2), n_interp
  real:: pt(2), old_score
  class(mesh_t), pointer:: m(:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! XY-planes in the two Z guard zones
  !-----------------------------------------------------------------------------
  call trace_begin('download::different_static', itimer=itimer)
  target => task2patch (link%task)
  n_interp = 0
  m => target%mesh
  nullify(source)
  old_score = -1.
  do iz=m(3)%lb,m(3)%lo
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ii = [ix,iy,iz]
        call interpolate_to (link, target, ii, source, jt, pt, -3, old_score, n_interp)
      end do
    end do
  end do
  call reset
  do iz=m(3)%uo,m(3)%ub
    do iy=m(2)%lb,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ii = [ix,iy,iz]
        call interpolate_to (link, target, ii, source, jt, pt, +3, old_score, n_interp)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! XZ-planes in the two Y guard zones
  !-----------------------------------------------------------------------------
  call reset
  do iz=m(3)%li,m(3)%ui
    do iy=m(2)%lb,m(2)%lo
      do ix=m(1)%lb,m(1)%ub
        ii = [ix,iy,iz]
        call interpolate_to (link, target, ii, source, jt, pt, -2, old_score, n_interp)
      end do
    end do
  end do
  call reset
  do iz=m(3)%li,m(3)%ui
    do iy=m(2)%uo,m(2)%ub
      do ix=m(1)%lb,m(1)%ub
        ii = [ix,iy,iz]
        call interpolate_to (link, target, ii, source, jt, pt, +2, old_score, n_interp)
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! YZ-planes in the two X guard zones
  !-----------------------------------------------------------------------------
  call reset
  do iz=m(3)%li,m(3)%ui
    do ix=m(1)%lb,m(1)%lo
      do iy=m(2)%li,m(2)%ui
        ii = [ix,iy,iz]
        call interpolate_to (link, target, ii, source, jt, pt, -1, old_score, n_interp)
      end do
    end do
  end do
  call reset
  do iz=m(3)%li,m(3)%ui
    do ix=m(1)%uo,m(1)%ub
      do iy=m(2)%li,m(2)%ui
        ii = [ix,iy,iz]
        call interpolate_to (link, target, ii, source, jt, pt, +1, old_score, n_interp)
      end do
    end do
  end do
  call reset

  if (n_interp > 0) then
    !$omp atomic
    download%n_differ = download%n_differ+1
  end if

  call trace_end (itimer)
contains
  subroutine reset
    if (associated(source)) then
      nullify (source)
    end if
    old_score = -1.
  end subroutine
END SUBROUTINE different_static

!> =============================================================================
!> =============================================================================
SUBROUTINE interpolate_to (link, target, ii, source, jt, pt, idir, old_score, n_interp)
  class(link_t), pointer:: link
  class(patch_t), pointer:: target, source
  integer:: ii(3), jt(2), idir, n_interp
  real:: pt(2), old_score
  integer:: jj(3), iv, i, it, j
  real(8):: pp(3), cellpos(3)
  real(kind=KindScalarVar):: p(4), q(4), qint1, qint2, result(target%nv)
  !-----------------------------------------------------------------------------
  ! Find the best nbor patch.
  !-----------------------------------------------------------------------------
  cellpos = target%myposition (ii)
  call find_best_nbor  (link, target, cellpos, source, jt, pt, idir, old_score)
  if (.not.associated(source)) then
    !if (target%id==121) write (io_unit%log,*) 'WARNING: no best_nbor for 121', ii, idir, old_score
    return
  else
    !if (target%id==121) write (io_unit%log,*) 'best_nbor for 121', source%id, ii, idir, old_score
  end if
  
  ! if the two patches are at the same level, then we should be using `same_static`.
  if (source%level == target%level) return

  !call trace_begin ('download_mod::interpolate_to',2)
  if (io%do_debug) write (io_unit%log,*) 'interpolate_to',target%it,ii,source%id,jt,pt
  it = target%it
  jj = 0
  do iv=1,source%nv
    do i=1,3
      cellpos(i) = target%position(i) + (ii(i)-target%offset(i)+target%mesh(i)%h(iv))*target%ds(i)      ! position of target ii
      pp(i) = (cellpos(i)-source%position(i))/source%ds(i) + source%offset(i) - source%mesh(i)%h(iv)    ! float index in source
      jj(i) = pp(i)                                                                                 ! integer part
      jj(i) = max(1,min(source%mesh(i)%gn-1,jj(i)))
      !---------------------------------------------------------------------------------------------
      ! Note that interpolating in the ghost zones is safe, since they are constantly maintained
      ! and up-to-date; first via extrapolation in time immediately after update, and then via the
      ! next download.  This is the best way to handle transitions in resolution when time slot
      ! history is available-- cf. the way it is done in RAMSES, where one insteads goes one level
      ! down in resolution
      !---------------------------------------------------------------------------------------------
      pp(i) = pp(i) - jj(i)                                                                         ! fractions
      if (pp(i)<-0.35d0 .or. pp(i)>1.35d0) then
        pp(i) = (cellpos(i)-source%position(i))/source%ds(i) + source%offset(i) - source%mesh(i)%h(iv)  ! float index in source
        if (io%master .and. io%verbose > 1) &
        write (io_unit%log,'(a,2i3,2f6.1,2(2x,a,i6,i3," i:",i4,f13.8))') 'INDEX ERROR:', idir, i, &
          target%mesh(i)%h(iv), source%mesh(i)%h(iv), &
          'target: id,lv=',target%id, target%level, ii(i), (source%position(i)-target%position(i))/source%mesh(i)%d, &
          'source: id,lv=',source%id, source%level, jj(i), pp(i)
        if (io%master .and. io%verbose > 2) then
          write (io_unit%log,*) '   source shift:', (target%position-source%position)/source%ds + source%offset
          write (io_unit%log,*) '    target size:', target%size
          write (io_unit%log,*) '   target index:', (cellpos(i)-target%position(i))/target%ds(i) + target%offset(i)
          write (io_unit%log,*) '    source size:', source%size
          write (io_unit%log,*) '   source index:', (cellpos(i)-source%position(i))/source%ds(i) + source%offset(i)
          write (io_unit%log,*) 'staggered index:', (cellpos(i)-source%position(i))/source%ds(i) + source%offset(i) - source%mesh(i)%h(iv)
        end if
        pp(i) = max(-0.35d0,min(pp(i),1.35d0))
      end if
    end do
    p(1:3)=pp; p(4)=pt(2)
    q = 1.-p
    if (target%id==io%id_debug) &
      write(io_unit%log,'(2(2i5,2x),2(3i4,2x),3f6.2)') target%id, idir, source%id, source%level, ii, jj, pp

    associate (qt1 => source%mem(:,:,:,:,jt(1),1), qt2 => source%mem(:,:,:,:,jt(2),1))
      qint1 = download%interpolate(qt1,jj(1),jj(2),jj(3),iv,p(1),p(2),p(3))
      qint2 = download%interpolate(qt2,jj(1),jj(2),jj(3),iv,p(1),p(2),p(3))
      result(iv) = q(4) * qint1 + p(4) * qint2
    end associate
  end do

  ! we now need to account for `target` and `source` potentially having different
  ! coord. systems, and therefore not a one-to-one mapping between vector components.
  call source%MapVectors(result, target, ii, jj, target%mem(ii(1),ii(2),ii(3),1:source%nv,it,1))
  n_interp = n_interp + 1

  !call trace_end
END SUBROUTINE interpolate_to

!> =============================================================================
!> Find the best nbor patch overlapping the point with position in link%patch
!> =============================================================================
SUBROUTINE find_best_nbor (link, target, position, best, jt, pt, idir, old_score, debug)
  class(link_t), pointer:: link
  class(patch_t), pointer:: target, best
  integer:: jt(2), idir
  real:: pt(2), old_score
  character(len=*), optional:: debug
  !.............................................................................
  class(link_t), pointer:: nbor
  class(patch_t), pointer:: other
  real(8):: dt, position(3), dist, pp(3)
  real:: score, best_score
  integer:: i, it, iv, adir, best_id
  !-----------------------------------------------------------------------------
  ! Check if the new point is fully inside the previous best patch.  If not,
  ! try finding a new one
  !-----------------------------------------------------------------------------
  adir = abs(idir)
  if (associated(best)) then
    if (fscore(best) == old_score) then
      if (target%id == io%id_debug) &
      print *,target%id,'keeping best nbor:', best%id, ' level:', best%level, ' direction:', idir
      return
    else
      nullify (best)
    end if
  end if
  call trace_begin('download_mod::find_best_nbor',2)
  best_id = 0
  best_score = 0.0
  nbor => link%nbor
  do while (associated(nbor))
    !---------------------------------------------------------------------------
    ! Skíp patches more than one level up/down, and the parent patch
    !---------------------------------------------------------------------------
    other => task2patch(nbor%task)
    if (abs(target%level-other%level)<2) then
      !-----------------------------------------------------------------------
      ! Design a score that is negative if any coordinate is outside the guard
      ! zone, and which counts +3 for an inside primary coordinate, and +1 for
      ! an inside secondary coordinate.  Levels count less, but enough to favor
      ! the higher level against a lower level
      !-----------------------------------------------------------------------
      score = fscore(other)
      if (score > best_score) then
          best => other
          best_score = score
          best_id = best%id
      end if
    end if
    if (target%id==io%id_debug.and.adir==3) write (*,'(i4,a,i6,2i4,3f8.2,5(3x,3f8.2))') &
      target%id,' search best_nbor', other%id, target%level, other%level, &
      score, best_score, old_score, abs(pp)/other%size, abs(pp)/other%gsize, pp, &
      (other%position-target%position)/target%size
    nbor => nbor%next
  end do
  !-----------------------------------------------------------------------------
  ! Find the right time interval in the other patch, and save score
  !-----------------------------------------------------------------------------
  if (associated(best)) then
    call best%time_interval (target%time, jt, pt)
    old_score = best_score
      if (target%id == io%id_debug) &
      print *,target%id,'picking best nbor:', best%id, ' level:', best%level, ' direction:', idir
  else
      if (target%id == io%id_debug) &
      print *,target%id,'failed'
  end if
  call trace_end
contains
  function fscore(test) result(score)
    class(patch_t):: test
    real:: score
    adir = abs(idir)
    pp = (position-test%position)
    !pp = target%point_distances (position, test%position)
    score = test%level*2.5
    do i=1,3
      score = score + merge(3.0,1.0,i==adir)* &
        merge(1.0,-9.0,abs(pp(i)) < test%size(i)*0.5d0)
    end do
  end function
END SUBROUTINE find_best_nbor

END MODULE
