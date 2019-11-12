!===============================================================================
!> This module handles remeshing, eg from lower resolution to higher resolution,
!> with special attention to maximizing speed, by identifying cases where cells
!> are aligned (static meshes).  FIXME: initially assumed to always be the case.
!===============================================================================
MODULE remesh_mod
  USE io_mod
  USE trace_mod
  USE omp_timer_mod
  USE patch_mod
  implicit none
  private
  type, public:: remesh_t
  contains
    procedure:: init
    procedure:: prolong
    procedure:: scatter_0th
    procedure:: scatter_1st
    procedure, nopass:: weights
  end type
  integer:: verbose=0
  integer:: order=0
  type(remesh_t), public:: remesh
CONTAINS

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE init (self)
  class(remesh_t):: self
  !.............................................................................
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /remesh_params/ verbose, order
  !-----------------------------------------------------------------------------
  call trace%begin ('remesh_t%init')
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      rewind (io_unit%input)
      read (io_unit%input, remesh_params, iostat=iostat)
      write (io_unit%output, remesh_params)
      first_time = .false.
    end if
    !$omp end critical (input_cr)
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Interface btw prolong and scatter.  In the basic case with no_mans_land=T,
!> the lower left corner in source space is the index point to the left of the
!> target index point; i.e., the source%index_only point corresponding to
!> target%myposition (target%mesh%li).
!===============================================================================
SUBROUTINE prolong (self, target, source, jt, pt, all_cells)
  class(remesh_t):: self
  class(patch_t), pointer:: target, source
  integer:: jt(2)
  real:: pt(2)
  logical:: all_cells
  !.............................................................................
  real(8):: pnt(3)
  integer:: ll(3), nn(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('remesh_t%prolong')
  call self%init
  if (all_cells) then
    pnt = target%myposition (target%mesh%li)
    if (order==0) then
      ll = source%index_only (pnt, roundup=.true.)
      nn = source%mesh%n/2
      call self%scatter_0th (source, target, ll, nn)
    else
      ll = source%index_only (pnt)
      nn = source%mesh%n/2+1
      call self%scatter_1st (source, target, ll, nn)
    end if
    if (verbose > 0) then
      write (io_unit%output,'(a,2i6,3f10.6,2(2x,3i4))') &
        'remesh_t%prolong: source, target, pnt, ll, nn =', &
        source%id, target%id, pnt, ll, nn
    end if
  else
    call io%abort('remesh_t%prolong: cannot handle all_cells=.false.')
  end if
  call trace%end ()
END SUBROUTINE prolong

!===============================================================================
!> Scatter = remesh from coarse to fine mesh.  The index corner ll gives the
!> lower left corner in index space.  If this is, for example, mesh%li, the
!> region covered will start at li1 for the coarse mesh, and at li for the
!> fine mesh.
!===============================================================================
SUBROUTINE scatter_0th (self, coarse, fine, ll, nn)
  class(remesh_t):: self
  class(patch_t), pointer:: coarse, fine
  integer:: ll(3), nn(3)
  optional:: nn
  !.............................................................................
  integer:: i1, i2, i3, j1, j2, j3, iv, lc(3), uc(3), lf(3), uf(3)
  real, dimension(:,:,:), pointer:: c, f
  real:: s(2)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('remesh_t%scatter_0th', itimer=itimer)
  lc = ll
  if (present(nn)) then
    uc = ll + nn - 1
  else
    uc = ll + coarse%mesh%n/2 - 1
  end if
  lf = fine%mesh%li
  uf = fine%mesh%ui
  !-----------------------------------------------------------------------------
  ! Loop over variables, duplicating the coarse value in all three directions
  !-----------------------------------------------------------------------------
  do iv=1,fine%nv
    c => coarse%mem(:,:,:,iv,coarse%it,1)
    f =>   fine%mem(:,:,:,iv,  fine%it,1)
    j3 = lf(3)
    do i3=lc(3),uc(3)
      j2 = lf(2)
      do i2=lc(2),uc(2)
        !-----------------------------------------------------------------------
        !-----------------------------------------------------------------------
        j1 = lf(1)
        do i1=lc(1),uc(1)
          f(j1  ,j2  ,j3  ) = c(i1,i2,i3)
          f(j1+1,j2  ,j3  ) = c(i1,i2,i3)
          f(j1  ,j2+1,j3  ) = c(i1,i2,i3)
          f(j1+1,j2+1,j3  ) = c(i1,i2,i3)
          f(j1  ,j2  ,j3+1) = c(i1,i2,i3)
          f(j1+1,j2  ,j3+1) = c(i1,i2,i3)
          f(j1  ,j2+1,j3+1) = c(i1,i2,i3)
          f(j1+1,j2+1,j3+1) = c(i1,i2,i3)
          j1 = j1+2
        end do
        j2 = j2+2
      end do
      j3 = j3+2
    end do
    if (verbose > 1) then
      s(1) = sum(c(lc(1):uc(1),lc(2):uc(2),lc(3):uc(3)))*product(coarse%ds)
      s(2) = sum(f(lf(1):uf(1),lf(2):uf(2),lf(3):uf(3)))*product(  fine%ds)
      write (io_unit%output,'(a,i4,1p,3e12.3,4(2x,3i4))') &
        'remesh_t%scatter_0th: iv, sum_c, sum_f =', &
        iv, s, fine%fmaxval(f), lc, uc, lf, uf
    end if
  end do
  call trace%end (itimer)
END SUBROUTINE scatter_0th

!===============================================================================
!> Scatter = remesh from coarse to fine mesh.  The index corner ll gives the
!> lower left corner in index space.  If this is, for example, mesh%li, the
!> region covered will start at li-1 for the coarse mesh, and at li-2 for the
!> fine mesh.  The li-2 point will lack contributions, so the the fine mesh will
!> only be valid from lo=li-1 to uo=ui+1.
!===============================================================================
SUBROUTINE scatter_1st (self, coarse, fine, ll, nn)
  class(remesh_t):: self
  class(patch_t), pointer:: coarse, fine
  integer:: ll(3), nn(3)
  optional:: nn
  !.............................................................................
  integer:: i1, i2, i3, j1, j2, j3, iv, j, lc(3), uc(3), lf(3), ic(3)
  real, dimension(:,:,:), pointer:: c, f
  real, allocatable:: out(:,:), w(:,:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('remesh_t%scatter_1st', itimer=itimer)
  ic = coarse%mesh%li
  lc = ll
  if (present(nn)) then
    uc = ll + nn
  else
    uc = ll + coarse%mesh%n/2 - 1
  end if
  lf = fine%mesh%li
  allocate (out(fine%gn(1),8), w(8,fine%gn(1)))
  !-----------------------------------------------------------------------------
  ! Loop over variables, computing remeshing weights for each type of centering
  !-----------------------------------------------------------------------------
  do iv=1,fine%nv
    c => coarse%mem(:,:,:,iv,coarse%it,1)
    f =>   fine%mem(:,:,:,iv,  fine%it,1)
    w = weights(iv)
    f = 0.0
    do i3=lc(3)-1,uc(3)+1
    do i2=lc(2)-1,uc(2)+1
      !-------------------------------------------------------------------------
      ! Loop over the 8 fine points for each coarse point, and over coarse points
      ! in the x-direction (this loop vectorizes)
      !-------------------------------------------------------------------------
      do j=1,8
        do i1=lc(1)-1,uc(1)+1
          out(i1,j) = w(1,j)*c(i1,i2  ,i3  ) + w(2,j)*c(i1+1,i2  ,i3  ) &
                    + w(3,j)*c(i1,i2+1,i3  ) + w(4,j)*c(i1+1,i2+1,i3  ) &
                    + w(5,j)*c(i1,i2  ,i3+1) + w(6,j)*c(i1+1,i2  ,i3+1) &
                    + w(7,j)*c(i1,i2+1,i3+1) + w(8,j)*c(i1+1,i2+1,i3+1)
        end do
      end do
      !-------------------------------------------------------------------------
      ! Add up the contributions to the 8 fine points for each coarse point.
      ! This loop does not vectorize, but has a significant load per iteration
      !-------------------------------------------------------------------------
      i1 = lc(1)-1
      j1 = lf(1) + (i1-ic(1))*2
      j2 = lf(2) + (i2-ic(2))*2
      j3 = lf(3) + (i3-ic(3))*2
      do i1=lc(1)-1,uc(1)+2
        !if (i1 >= lc(1) .and. i2 >= lc(2) .and. i3 >= lc(3)) &
          f(j1  ,j2  ,j3  ) = f(j1  ,j2  ,j3  ) + out(i1,1)
        !if (i1 <= uc(2) .and. i2 >= lc(2) .and. i3 >= lc(3)) &
          f(j1+1,j2  ,j3  ) = f(j1+1,j2  ,j3  ) + out(i1,2)
        !if (i1 >= lc(1) .and. i2 <= uc(2) .and. i3 >= lc(3)) &
          f(j1  ,j2+1,j3  ) = f(j1  ,j2+1,j3  ) + out(i1,3)
        !if (i1 <= uc(1) .and. i2 <= uc(2) .and. i3 >= lc(3)) &
          f(j1+1,j2+1,j3  ) = f(j1+1,j2+1,j3  ) + out(i1,4)
        !if (i1 >= lc(1) .and. i2 >= lc(2) .and. i3 <= uc(3)) &
          f(j1  ,j2  ,j3+1) = f(j1  ,j2  ,j3+1) + out(i1,5)
        !if (i1 <= uc(1) .and. i2 >= lc(2) .and. i3 <= uc(3)) &
          f(j1+1,j2  ,j3+1) = f(j1+1,j2  ,j3+1) + out(i1,6)
        !if (i1 >= lc(1) .and. i2 <= uc(2) .and. i3 <= uc(3)) &
          f(j1  ,j2+1,j3+1) = f(j1  ,j2+1,j3+1) + out(i1,7)
        !if (i1 <= uc(1) .and. i2 <= uc(2) .and. i3 <= uc(3)) &
          f(j1+1,j2+1,j3+1) = f(j1+1,j2+1,j3+1) + out(i1,8)
        j1 = j1+2
      end do
    end do
    end do
  end do
  deallocate (out)
  call trace%end (itimer)
END SUBROUTINE scatter_1st

!===============================================================================
!> Weight factors for RAMSES type arrangement of coarse / fine (no_mans_land=t).
!> The weight w(i,j) is the contribution of coarse point i to the value at fine
!> point j.  Both i and j are numbered in increasing order of first x, then y,
!> then z.
!>
!> FIXME: There should be a solver-specific weight function defined for each
!> solver -- this one is only intended for testing.
!===============================================================================
FUNCTION weights (iv) RESULT (w)
  integer:: iv
  real:: w(8,8), w1, w2, w3
  integer:: i1, i2, i3, i, j1, j2, j3, j
  !-----------------------------------------------------------------------------
  w = 0.0
  i = 1
  do i3=0,1
    do i2=0,1
      do i1=0,1
        j = 1
        do j3=0,1
          w3 = merge(0.75,0.25,i3==j3)
          do j2=0,1
            w2 = merge(0.75,0.25,i2==j2)
            !-------------------------------------------------------------------
            ! i1 is the coarse point location in x: 0 for 0.0, 1 for 1.0
            ! j1 is the   fine point location in x: 0 for .25, 1 for .75
            !-------------------------------------------------------------------
            do j1=0,1
              w1 = merge(0.75,0.25,i1==j1)
              w(i,j) = w1*w2*w3
              j = j+1
            end do
          end do
        end do
        i = i+1
      end do
    end do
  end do
END FUNCTION weights

END MODULE remesh_mod
