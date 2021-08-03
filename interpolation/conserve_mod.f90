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
    procedure:: restrict
    procedure:: test
  end type
  type(conserve_t), public:: conserve
CONTAINS

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
    l(i) = index (target%mesh(i), m%p + m%r(m%li))
    u(i) = index (target%mesh(i), m%p + m%r(m%ui))
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