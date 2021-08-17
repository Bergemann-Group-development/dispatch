MODULE disk_selfgravity_mod
  USE selfgravity_mod
  USE poisson_mod
  USE link_mod
  USE task_mod
  USE kinds_mod
  USE patch_mod
  USE mesh_mod
  implicit none
  private
  type, public, extends(selfgravity_t):: disk_selfgravity_t
  contains
    procedure:: downloads
    procedure:: uploads
    !procedure:: restrict
    !procedure:: prolong
  end type
CONTAINS

!===============================================================================
!> Download residuals (only) from higher level patches that overlap with the
!> internal cells of the patch
!> Verbatim copy form selfgravity_mod, with correct calculation of distance
!> between patches
!===============================================================================
SUBROUTINE downloads (self)
  class(disk_selfgravity_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: patch, nbpatch
  !-----------------------------------------------------------------------------
  ! Download values from the L+1 tasks contained in this task
  !-----------------------------------------------------------------------------
  patch => self%patch
  select type (patch)
    class is (patch_t)
    self%d   => patch%mem(:,:,:,patch%idx%d   ,patch%it,1)
    self%phi => patch%mem(:,:,:,patch%idx%phi ,patch%it,1)
    self%res => patch%mem(:,:,:,patch%idx%rphi,patch%it,1)
    self%src = self%fourPiG*self%d
    call poisson%laplace_2nd (patch%ds, self%phi, self%src, self%res)
    nbor => patch%link%nbor
    do while (associated(nbor))
      nbpatch => nbor%task
      select type (nbpatch)
      class is (patch_t)
        if (nbpatch%level == patch%level+1 ) then
          if (all(abs(patch%distance_to(nbpatch)) < 0.5d0*patch%size)) then
            call self%restrict (nbpatch, patch, patch%idx%rphi)
          end if
        end if
      end select
      nbor => nbor%next
    end do
  end select
END SUBROUTINE downloads

!===============================================================================
!> Upload potentials corrections (only) to higher level patches that overlap
!> with the internal cells of the patch
!> Verbatim copy form selfgravity_mod, with correct calculation of distance
!> between patches
!===============================================================================
SUBROUTINE uploads (self)
  class(disk_selfgravity_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: patch, nbpatch
  !-----------------------------------------------------------------------------
  ! Upload values to the L-1 tasks on the nbor list
  !-----------------------------------------------------------------------------
  patch => self%patch
  select type (patch)
    class is (patch_t)
    self%d   => patch%mem(:,:,:,patch%idx%d  ,patch%it,1)
    self%phi => patch%mem(:,:,:,patch%idx%phi,patch%it,1)
    nbor => patch%link%nbor
    do while (associated(nbor))
      nbpatch => nbor%task
      select type (nbpatch)
        class is (patch_t)
        if (nbpatch%level == patch%level-1) then
          if (all(abs(nbpatch%distance_to(patch)) < 0.5d0*nbpatch%size)) then
            call self%prolong (nbpatch, patch, patch%idx%phi)
          end if
        end if
      end select
      nbor => nbor%next
    end do
  end select
END SUBROUTINE uploads

!===============================================================================
!> Restrict
!> Verbatim copy form selfgravity_mod, with correct calculation of distance
!> between patches as well as conversion from source to target coordinate
!> system
!===============================================================================
SUBROUTINE restrict (self, source, target, iv)
  class(disk_selfgravity_t):: self
  class(patch_t)          :: source, target
  integer                 :: iv
  real, pointer           :: src(:,:,:,:), tgt(:,:,:)
  integer                 :: i1, i2, i3, j1, j2, j3, k(2), l(3), u(3)
  real                    :: pt(2), p1, p2, p3, offset(3)
  real(8)                 :: d(3), fl(3), fu(3), size(3), ds(3)
  logical                 :: outside(3)
  class(mesh_t), pointer  :: sm, tm
  !-----------------------------------------------------------------------------
  call source%time_interval (target%time, k, pt, &
    extrap=target%allow_extrapolate)
  src => source%mem(:,:,:,iv,:,1)
  tgt => target%mem(:,:,:,iv,target%it,1)
  d = target%distance_to (source)
  !-----------------------------------------------------------------------------
  ! Convert source%size, source%mesh%d, and source%mesh%o to coordinate system
  ! of target
  !-----------------------------------------------------------------------------
  size = matmul (matmul (source%e, source%size), target%e)
  ds = matmul (matmul (source%e, source%ds), target%e)
  offset = matmul (matmul (source%e, source%offset), target%e)
  fl = d - size/2d0
  fu = d + size/2d0
  l = ceiling(fl/target%mesh%d + target%mesh%o) + 1
  u =   floor(fu/target%mesh%d + target%mesh%o) - 1
  do i3=l(3),u(3)
    tm => target%mesh(3)
    outside(3) = i3 < tm%li .or. i3 > tm%ui
    p3 = (tm%r(i3) + d(3))/ds(3) + offset(3)
    j3 = p3
    p3 = p3 - j3
    do i2=l(2),u(2)
      tm => target%mesh(2)
      outside(2) = i2 < tm%li .or. i2 > tm%ui
      p2 = (tm%r(i2) + d(2))/ds(2) + offset(2)
      j2 = p2
      p2 = p2 - j2
      do i1=l(1),u(1)
        outside(1) = i1 < tm%li .or. i1 > tm%ui
        if (any(outside)) then
          tm => target%mesh(1)
          p1 = (tm%r(i1) + d(1))/ds(1) + offset(1)
          j1 = p1
          p1 = p1 - j1
          tgt(i1,i2,i3) = &
              pt(1)*((1.-p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3  ,k(1))    + &
                                       (    p1)*src(j1+1,j2  ,j3  ,k(1)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3  ,k(1))    + &
                                       (    p1)*src(j1+1,j2+1,j3  ,k(1))))  + &
                     (   p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3+1,k(1))    + &
                                       (    p1)*src(j1+1,j2  ,j3+1,k(1)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3+1,k(1))    + &
                                       (    p1)*src(j1+1,j2+1,j3+1,k(1))))) + &
              pt(2)*((1.-p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3  ,k(2))    + &
                                       (    p1)*src(j1+1,j2  ,j3  ,k(2)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3  ,k(2))    + &
                                       (    p1)*src(j1+1,j2+1,j3  ,k(2))))  + &
                     (   p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3+1,k(2))    + &
                                       (    p1)*src(j1+1,j2  ,j3+1,k(2)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3+1,k(2))    + &
                                      (    p1)*src(j1+1,j2+1,j3+1,k(2)))))
        end if
      end do
    end do
  end do
END SUBROUTINE restrict

!===============================================================================
!> Prolong
!> Verbatim copy form selfgravity_mod, with correct calculation of distance
!> between patches as well as conversion from source to target coordinate
!> system
!===============================================================================
SUBROUTINE prolong (self, source, target, iv)
  class(disk_selfgravity_t):: self
  class(patch_t)          :: source, target
  integer                 :: iv
  real, pointer           :: src(:,:,:,:), tgt(:,:,:)
  integer                 :: i1, i2, i3, j1, j2, j3, k(2)
  real                    :: pt(2), p1, p2, p3, offset(3)
  real(8)                 :: d(3), size(3), ds(3)
  logical                 :: outside(3)
  class(mesh_t), pointer  :: sm, tm
  integer                 :: tid, sid
  !-----------------------------------------------------------------------------
  call source%time_interval (target%time, k, pt, &
    extrap=target%allow_extrapolate)
  src => source%mem(:,:,:,iv,:,1)
  tgt => target%mem(:,:,:,iv,target%it,1)
  tid = target%id
  sid = source%id
  d = target%distance_to (source)
  !-----------------------------------------------------------------------------
  ! Convert source%mesh%d, and source%mesh%o to coordinate system
  ! of target
  !-----------------------------------------------------------------------------
  ds = matmul (matmul (source%e, source%ds), target%e)
  offset = matmul (matmul (source%e, source%offset), target%e)
  do i3=target%mesh(3)%lb+1,target%mesh(3)%ub-1
    tm => target%mesh(3)
    sm => source%mesh(3)
    outside(3) = i3 < tm%li .or. i3 > tm%ui
    p3 = (tm%r(i3) + d(3))/ds(3) + offset(3)
    j3 = p3
    if (j3 > sm%gn-1) then
      j3 = sm%gn - 1
    end if
    p3 = p3 - j3
    do i2=target%mesh(2)%lb+1,target%mesh(2)%ub-1
      tm => target%mesh(2)
      sm => source%mesh(2)
      outside(2) = i2 < tm%li .or. i2 > tm%ui
      p2 = (tm%r(i2) + d(2))/ds(2) + offset(2)
      j2 = p2
      if (j2 > sm%gn-1) then
        j2 = sm%gn - 1
      end if
      p2 = p2 - j2
      do i1=target%mesh(1)%lb+1,target%mesh(1)%ub-1
        outside(1) = i1 < tm%li .or. i1 > tm%ui
        if (any(outside)) then
          tm => target%mesh(1)
          sm => source%mesh(1)
          p1 = (tm%r(i1) + d(1))/ds(1) + offset(1)
          j1 = p1
          if (j1 > sm%gn-1) then
            j1 = sm%gn - 1
          end if
          if (j1 < 1) then
            j1 = 1
          end if
          p1 = p1 - j1
          tgt(i1,i2,i3) = &
              pt(1)*((1.-p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3  ,k(1))    + &
                                       (    p1)*src(j1+1,j2  ,j3  ,k(1)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3  ,k(1))    + &
                                       (    p1)*src(j1+1,j2+1,j3  ,k(1))))  + &
                     (   p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3+1,k(1))    + &
                                       (    p1)*src(j1+1,j2  ,j3+1,k(1)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3+1,k(1))    + &
                                       (    p1)*src(j1+1,j2+1,j3+1,k(1))))) + &
              pt(2)*((1.-p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3  ,k(2))    + &
                                       (    p1)*src(j1+1,j2  ,j3  ,k(2)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3  ,k(2))    + &
                                       (    p1)*src(j1+1,j2+1,j3  ,k(2))))  + &
                     (   p3)*((1.-p2)*((1.0-p1)*src(j1  ,j2  ,j3+1,k(2))    + &
                                       (    p1)*src(j1+1,j2  ,j3+1,k(2)))   + &
                              (   p2)*((1.0-p1)*src(j1  ,j2+1,j3+1,k(2))    + &
                                      (    p1)*src(j1+1,j2+1,j3+1,k(2)))))
        end if
      end do
    end do
  end do
END SUBROUTINE prolong

END MODULE disk_selfgravity_mod
