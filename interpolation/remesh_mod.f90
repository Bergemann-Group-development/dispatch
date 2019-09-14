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
    procedure, nopass:: weights
  end type
  type(remesh_t), public:: remesh
CONTAINS

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE init (self)
  class(remesh_t):: self
END SUBROUTINE init

!===============================================================================
!> Prolong = remesh from coarse to fine mesh.  The index corner ll gives the
!> lower left corner in index space.  If this is, for example, mesh%li, the
!> region covered will start at li-1 for the coarse mesh, and at li-2 for the
!> fine mesh.  The li-2 point will lack contributions, so the the fine mesh will
!> only be valid from lo=li-1 to uo=ui+1.
!===============================================================================
SUBROUTINE prolong (self, coarse, fine, ll, nn)
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
  call trace%begin ('remesh_t%prolong', itimer=itimer)
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
END SUBROUTINE prolong

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
