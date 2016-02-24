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
!> patches at the same level, with no other types even allowed in the nbor list
!> =============================================================================
SUBROUTINE download_link (link, all_cells)
  type(link_t), pointer:: link
  logical, optional:: all_cells
  type(link_t), pointer:: nbor
  class(patch_t), pointer:: patch
  logical:: static, same
  integer:: nt, nv
  !-----------------------------------------------------------------------------
  call trace_begin('download_t%download')
  patch => task2patch(link%task)
  nt = patch%nt
  nv = patch%nv
  static = .true.; same = .true.
  if (patch%is_clear(bits%static)) static = .false.
  nbor => link%nbor
  do while (associated(nbor))
    print*,patch%id,nbor%task%id,nbor%task%is_clear(bits%static),nbor%task%level
    if (nbor%task%is_clear(bits%static)) static = .false.
    if (nbor%task%level /= patch%level) same = .false.
    nbor => nbor%next
  end do
  if (static) then
    if (same) then
      call download_same_static (link)
    else
      print *, 'ERROR: different level patches'
    end if
  else
    print *, 'ERROR: moving patches found'
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
  type(link_t), pointer:: link, nbor
  class(task_t), pointer:: task
  integer:: ix, iy, iz, o, jt1, jt2, it, iv, nt, nv
  real:: pt, qt
  type(mesh_t), pointer:: m(:)
  !-----------------------------------------------------------------------------
  ! XY-planes in the two Z guard zones
  !-----------------------------------------------------------------------------
  call trace_begin('download::download_same_static')
  patch => task2patch (link%task)
  m => patch%mesh
  it = patch%it
  nv = patch%nv
  print*,nv
  other => find_nbor (link, -3, o, patch%time, jt1, jt2, pt, qt)
  do iz=m(3)%lo,m(3)%li-1
    do iy=m(2)%lo,m(2)%uo
      do iv=1,nv
        do ix=m(1)%lo,m(1)%uo
          patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy,iz+o,iv,jt1,1)*qt &
                                      + other%mem(ix,iy,iz+o,iv,jt2,1)*pt
        end do
      end do
    end do
  end do
  other => find_nbor (link, +3, o, patch%time, jt1, jt2, pt, qt)
  do iz=m(3)%ui+1,m(3)%uo
    do iy=m(3)%lo,m(3)%uo
      do ix=m(1)%lo,m(1)%uo
        do iv=1,nv
          patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy,iz+o,iv,jt1,1)*qt &
                                      + other%mem(ix,iy,iz+o,iv,jt2,1)*pt
        end do
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! XZ-planes in the two Y guard zones
  !-----------------------------------------------------------------------------
  other => find_nbor (link, -2, o, patch%time, jt1, jt2, pt, qt)
  do iz=m(3)%li,m(3)%ui
    do iy=m(2)%lo,m(2)%li-1
      do iv=1,nv
        do ix=m(1)%lo,m(1)%uo
          patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy+o,iz,iv,jt1,1)*qt &
                                      + other%mem(ix,iy+o,iz,iv,jt2,1)*pt
        end do
      end do
    end do
  end do
  other => find_nbor (link, +2, o, patch%time, jt1, jt2, pt, qt)
  do iz=m(3)%li,m(3)%ui
    do iy=m(2)%ui+1,m(2)%uo
      do iv=1,nv
        do ix=m(1)%lo,m(1)%uo
          patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix,iy+o,iz,iv,jt1,1)*qt &
                                      + other%mem(ix,iy+o,iz,iv,jt2,1)*pt
        end do
      end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! YZ-planes in the two X guard zones
  !-----------------------------------------------------------------------------
  other => find_nbor (link, -1, o, patch%time, jt1, jt2, pt, qt)
  do iz=m(3)%li,m(3)%ui
    do iy=m(2)%lo,m(2)%li-1
      do ix=m(1)%lo,m(1)%li-1
        do iv=1,nv
          patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix+o,iy,iz,iv,jt1,1)*qt &
                                      + other%mem(ix+o,iy,iz,iv,jt2,1)*pt
        end do
      end do
    end do
  end do
  other => find_nbor (link, +1, o, patch%time, jt1, jt2, pt, qt)
  do iz=m(3)%li,m(3)%ui
    do iy=m(2)%li,m(2)%ui
      do ix=m(1)%ui+1,m(1)%uo
        do iv=1,nv
          patch%mem(ix,iy,iz,iv,it,1) = other%mem(ix+o,iy,iz,iv,jt1,1)*qt &
                                      + other%mem(ix+o,iy,iz,iv,jt2,1)*pt
        end do
      end do
    end do
  end do
  call trace_end
END SUBROUTINE download_same_static

!> =============================================================================
!> Find the nbor task in direction idir (plus/minos)
!> =============================================================================
FUNCTION find_nbor (link, idir, o, time, jt1, jt2, pt, qt) RESULT (other)
  type(link_t), pointer:: link, nbor
  class(patch_t), pointer:: patch, other
  type(mesh_t), pointer:: m(:)
  real(8):: time, position(3), dist
  integer:: idir, adir, o, jt1, jt2, it, iv, nv, nt
  real:: pt, qt, dt
  !-----------------------------------------------------------------------------
  call trace_begin('download::find_nbor',1)
  patch => task2patch (link%task)
  m => patch%mesh
  nt = patch%nt
  nv = patch%nv
  nbor => link%nbor
  nullify(other)
  adir = abs(idir)
  position = patch%position + patch%size(abs(idir))*idir/adir
  do while (associated(nbor))
    other => task2patch (nbor%task)
    dist = point_distance(patch, other%position, position)
    print*,patch%it,nbor%task%id, dist
    if (dist < patch%size(abs(idir))*0.5) then
      exit
    end if
    nbor => nbor%next
  end do
  o = -(other%mesh(adir)%ui+1-patch%mesh(adir)%lo)*idir/abs(idir)
  print*,associated(other)
  !-----------------------------------------------------------------------------
  ! Find the right time interval in the other patch
  !-----------------------------------------------------------------------------
  jt1=1
  do it=1,nt-2
    if (other%t(other%iit(it))<time) jt1=it
    print*,it,time, other%t(other%iit(it))
  end do
  jt2 = jt1+1
  jt1 = other%iit(jt1)
  jt2 = other%iit(jt2)
  dt = other%t(jt2)-other%t(jt1)
  pt = (time-other%t(jt1))/max(dt,1e-30)
  qt = 1.-pt
  if (io%verbose>0) &
    print '(a,i3,i7,i4,i4,i3,2f6.2)', 'download::find_nbor: idir, id, o, jt, pt', &
      idir, other%id, o, jt1, jt2, pt
  call trace_end
END FUNCTION find_nbor

END MODULE
