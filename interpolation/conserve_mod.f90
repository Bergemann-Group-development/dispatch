!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> Module for conservative interpolation operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE conserve_mod
  USE trace_mod
  USE io_mod
  USE patch_mod
  USE mesh_mod
  USE shared_mod
  USE interp_mod
  implicit none
  private
  type, public:: conserve_t
    class(patch_t), pointer:: patch => null()
  contains
    procedure:: prolong
    procedure:: restrict
    procedure:: test
  end type
  type(conserve_t), public:: conserve
CONTAINS

!===============================================================================
!> Conservative restrict operation, in the special case with aligned source and
!> target cells
!===============================================================================
SUBROUTINE prolong (self, target, source, all_src, only, count)
  class(conserve_t)      :: self
  class(patch_t)         :: target, source
  class(mesh_t), pointer :: m, m1, m2, m3
  logical,      optional :: all_src, count
  integer,      optional :: only
  !.............................................................................
  integer                :: i1, i2, i3, i4, iv0, iv1, j1, j2, j3, jt(2), l, n
  integer,  dimension(3) :: lt, ut, ls, us, k
  real(8)                :: p(3), d(3)
  real                   :: pt(2), c(2)
  integer,          save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('conserve_t%prolong', itimer=itimer)
  d = target%distance_to (source)
  !.............................................................................
  if (target%id == io%id_debug) &
    write (stdout,'(a,2(i6,":",i2.2,3i5),3x,3g12.3)') &
      'conserve_t%prolong: target,source,distance/size=', &
      target%id, target%key, source%id, source%key, d/target%size
  !.............................................................................
  call source%time_interval (target, target%time, jt, pt)
  if (present(only)) then
    iv0 = only
    iv1 = only
  else
    iv0 = 1
    iv1 = target%nv
  end if
  if (option(all_src)) then
    ls = source%mesh%lb
    us = source%mesh%ub
  else
    ls = source%mesh%li
    us = source%mesh%ui
  end if
  lt = target%mesh%li
  ut = target%mesh%ui
  do i4=iv0,iv1
    !---------------------------------------------------------------------------
    ! [i1,i2,i3] is the target index coordinates of the cell to give a value
    !---------------------------------------------------------------------------
    !$omp parallel do collapse(2) private (i1,i2,i3,p,k,j1,j2,j3)
    do i3=lt(3),ut(3)
      do i2=lt(2),ut(2)
        do i1=lt(1),ut(1)
          !---------------------------------------------------------------------
          ! p is the target cell center location in target index coordinates,
          ! to which we add 27 contributions from source cells with rounded
          ! integer indices k(:)
          !---------------------------------------------------------------------
          if (target%filled(i1,i2,i3) < target%level) then
            target%mem(i1,i2,i3,i4,target%it,1) = 0.0
            c = pt/27.
            do j3=-1,1
              do j2=-1,1
                do j1=-1,1
                  p = ([i1+j1,i2+j2,i3+j3]-target%offset)*target%ds
                  k = (p - d)/source%ds + source%offset
                  k = max(ls,min(k,us))
                    target%mem(i1,i2,i3,i4,target%it,1) = &
                    target%mem(i1,i2,i3,i4,target%it,1) + &
                    source%mem(k(1),k(2),k(3),i4,jt(1),1)*c(1) + &
                    source%mem(k(1),k(2),k(3),i4,jt(2),1)*c(2)
                end do
              end do
            end do
          end if
        end do
      end do
    end do
    !$omp end parallel do
  end do
  !-----------------------------------------------------------------------------
  ! Mark filled
  !-----------------------------------------------------------------------------
  n = 0
  do i3=lt(3),ut(3)
    do i2=lt(2),ut(2)
      do i1=lt(1),ut(1)
        if (target%filled(i1,i2,i3) < target%level) then
          target%filled(i1,i2,i3) =  target%level
          n = n + 1
        end if
      end do
    end do
  end do
  if (option(count)) &
    print '("ancestor",i6,":",i2.2," filled",i6," out of",i6," cells")', &
    target%id, target%level, n, product(ut-lt+1)
  call trace%end (itimer)
END SUBROUTINE prolong

!===============================================================================
!> Convinience function
!===============================================================================
LOGICAL FUNCTION option (flag)
  logical, optional:: flag
  if (present(flag)) then
    option = flag
  else
    option = .false.
  end if
END FUNCTION option

!===============================================================================
!> Conservative restrict operation, in the special case with aligned source and
!> target cells
!===============================================================================
SUBROUTINE restrict (self, target, source, all_tgt, all_src, only)
  class(conserve_t):: self
  class(patch_t):: target, source
  class(mesh_t), pointer:: m, m1, m2, m3
  logical, optional:: all_tgt, all_src
  integer, optional:: only
  !.............................................................................
  real, allocatable:: src(:,:,:,:)
  integer:: i, i1, i2, i3, iv, iv0, iv1, j1, j2, j3, tid
  integer:: j1m, j2m, j3m, j1p, j2p, j3p
  integer, dimension(3):: j, l, u, ls, us
  real:: c
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('conserve_t%restrict', itimer=itimer)
  if (target%id == io%id_debug) then
    print *,'DEBUG'
  end if
  if (present(only)) then
    iv0 = only
    iv1 = only
  else
    iv0 = 1
    iv1 = target%nv
  end if
  do i=1,3
    !---------------------------------------------------------------------------
    ! l(:) is the index of the nearest cell center in the coarse target mesh
    ! that corresponds to the first internal point of the fine source mesh
    !---------------------------------------------------------------------------
    m => source%mesh(i)
    if (present(all_src)) then
      l(i) = index (target%mesh(i), m%p + m%r(m%lb))
      u(i) = index (target%mesh(i), m%p + m%r(m%ub))
    else
      l(i) = index (target%mesh(i), m%p + m%r(m%li))
      u(i) = index (target%mesh(i), m%p + m%r(m%ui))
    end if
    l(i) = max(l(i),m%lb)
    u(i) = min(u(i),m%ub)
    !---------------------------------------------------------------------------
    ! j(:) is the first mid-point in the source, which corresponds in position
    ! to the first point in the coarse (target) mesh
    !---------------------------------------------------------------------------
    m => target%mesh(i)
    j(i) = index (source%mesh(i), m%p + m%r(l(i)))
  end do
  allocate (src(source%gn(1),source%gn(2),source%gn(3),source%mv))
  if (present(all_src)) then
    ls = source%mesh%lb
    us = source%mesh%ub
  else
    ls = source%mesh%li
    us = source%mesh%ui
  end if
  call interp%in_t (target, source, ls, us, source%mem, src, only=only)
  if (shared%refine_ratio == 2.0) then
    c = 1./8.
    do iv=iv0,iv1
      j3 = j(3)
      do i3=l(3),u(3)
        j2 = j(2)
        do i2=l(2),u(2)
          j1 = j(1)
          do i1=l(1),u(1)
            target%mem(i1,i2,i3,iv,target%it,1) = &
             (src(j1  ,j2  ,j3  ,iv)  + &
              src(j1-1,j2  ,j3  ,iv) + &
              src(j1  ,j2-1,j3  ,iv) + &
              src(j1-1,j2-1,j3  ,iv) + &
              src(j1  ,j2  ,j3-1,iv) + &
              src(j1-1,j2  ,j3-1,iv) + &
              src(j1  ,j2-1,j3-1,iv) + &
              src(j1-1,j2-1,j3-1,iv))*c
              j1 = j1 + 2
            end do
            j2 = j2 + 2
          end do
          j3 = j3 + 2
      end do
    end do
  else if (shared%refine_ratio == 3.0) then
    c = 1./27.
    do iv=iv0,iv1
      j3 = j(3)
      do i3=l(3),u(3)
        j2 = j(2)
        do i2=l(2),u(2)
          j1 = j(1)
          if (j(1) > 1) then
            do i1=l(1),u(1)
              target%mem(i1,i2,i3,iv,target%it,1) = &
               (src(j1-1,j2-1,j3-1,iv) + &
                src(j1  ,j2-1,j3-1,iv) + &
                src(j1+1,j2-1,j3-1,iv) + &
                src(j1-1,j2  ,j3-1,iv) + &
                src(j1  ,j2  ,j3-1,iv) + &
                src(j1+1,j2  ,j3-1,iv) + &
                src(j1-1,j2+1,j3-1,iv) + &
                src(j1  ,j2+1,j3-1,iv) + &
                src(j1+1,j2+1,j3-1,iv) + &
                src(j1-1,j2-1,j3  ,iv) + &
                src(j1  ,j2-1,j3  ,iv) + &
                src(j1+1,j2-1,j3  ,iv) + &
                src(j1-1,j2  ,j3  ,iv) + &
                src(j1  ,j2  ,j3  ,iv) + &
                src(j1+1,j2  ,j3  ,iv) + &
                src(j1-1,j2+1,j3  ,iv) + &
                src(j1  ,j2+1,j3  ,iv) + &
                src(j1+1,j2+1,j3  ,iv) + &
                src(j1-1,j2-1,j3+1,iv) + &
                src(j1  ,j2-1,j3+1,iv) + &
                src(j1+1,j2-1,j3+1,iv) + &
                src(j1-1,j2  ,j3+1,iv) + &
                src(j1  ,j2  ,j3+1,iv) + &
                src(j1+1,j2  ,j3+1,iv) + &
                src(j1-1,j2+1,j3+1,iv) + &
                src(j1  ,j2+1,j3+1,iv) + &
                src(j1+1,j2+1,j3+1,iv))*c
              j1 = j1 + 3
            end do
          else
            do i1=l(1),u(1)
              j1m = max(j1-1,1)
              j2m = max(j2-1,1)
              j3m = max(j3-1,1)
              j1p = min(j1+1,source%gn(1))
              j2p = min(j2+1,source%gn(2))
              j3p = min(j3+1,source%gn(3))
              target%mem(i1,i2,i3,iv,target%it,1) = &
               (src(j1m,j2m,j3m,iv) + &
                src(j1 ,j2m,j3m,iv) + &
                src(j1p,j2m,j3m,iv) + &
                src(j1m,j2 ,j3m,iv) + &
                src(j1 ,j2 ,j3m,iv) + &
                src(j1p,j2 ,j3m,iv) + &
                src(j1m,j2p,j3m,iv) + &
                src(j1 ,j2p,j3m,iv) + &
                src(j1p,j2p,j3m,iv) + &
                src(j1m,j2m,j3 ,iv) + &
                src(j1 ,j2m,j3 ,iv) + &
                src(j1p,j2m,j3 ,iv) + &
                src(j1m,j2 ,j3 ,iv) + &
                src(j1 ,j2 ,j3 ,iv) + &
                src(j1p,j2 ,j3 ,iv) + &
                src(j1m,j2p,j3 ,iv) + &
                src(j1 ,j2p,j3 ,iv) + &
                src(j1p,j2p,j3 ,iv) + &
                src(j1m,j2m,j3p,iv) + &
                src(j1 ,j2m,j3p,iv) + &
                src(j1p,j2m,j3p,iv) + &
                src(j1m,j2 ,j3p,iv) + &
                src(j1 ,j2 ,j3p,iv) + &
                src(j1p,j2 ,j3p,iv) + &
                src(j1m,j2p,j3p,iv) + &
                src(j1 ,j2p,j3p,iv) + &
                src(j1p,j2p,j3p,iv))*c
              j1 = j1 + 3
            end do
          end if
          j2 = j2 + 3
        end do
        j3 = j3 + 3
      end do
    end do
  else
    write (stderr,*) 'refine_ratio =', shared%refine_ratio
    call io%abort ('conserve_t%restrict: refine_ratio not supported')
  end if
  deallocate (src)
  call trace%end (itimer)
contains
  !=============================================================================
  !> Return the index of the point in the give mesh
  !=============================================================================
  function index (m, p)
    integer:: index
    class(mesh_t):: m
    real(8):: p
    !...........................................................................
    index = nint ((p-m%p)/m%d + m%o)
  end function index
END SUBROUTINE restrict

SUBROUTINE test (self)
  class(conserve_t):: self
  call test1 (3.0, 15)
  call test1 (2.0, 16)
END SUBROUTINE test

SUBROUTINE test1 (refine_ratio, n)
  real:: refine_ratio
  integer:: n
  !.............................................................................
  type(patch_t):: source, target
  !-----------------------------------------------------------------------------
  ! Note that the number of points must be divisible by the refinement factor
  !-----------------------------------------------------------------------------
  call source%pre_init
  call target%pre_init
  target%size = target%box
  source%n = n
  target%n = n
  shared%refine_ratio = refine_ratio
  source%size = target%size/shared%refine_ratio
  call source%init
  call target%init
  source%idx%d = 1
  target%idx%d = 1
  call conserve%restrict (target, source, only=source%idx%d)
  call target%dealloc
  call source%dealloc
END SUBROUTINE test1

END MODULE conserve_mod
