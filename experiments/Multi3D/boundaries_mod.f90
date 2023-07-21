!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE standard_boundaries_mod
  USE eos_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE patch_mod
  USE units_mod
  USE op_math_mod, only: pchip_mono
  ! ----------------------------------------------------------------------------
  ! Experiment - specific USE modules
  ! ----------------------------------------------------------------------------
  !USE lower_damp_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
  contains
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
    procedure:: mirror
    procedure:: extrap
  end type
  integer, save:: verbose=1
CONTAINS

!===============================================================================
SUBROUTINE mirror(self)
  class(boundaries_t):: self
  integer :: i, b, it
  ! ----------------------------------------------------------------------------
  if (self%patch%ng(3) == 0) return
  it = self%patch%it

  if (self%is_set(bits%zu)) then
      b = ubound(self%patch%mem,3)
      self%patch%mem(:,:,b,:,it,:) = self%patch%mem(:,:,b-1,:,it,:)
  endif

  if (self%is_set(bits%zl)) then
      b = lbound(self%patch%mem,3)
      self%patch%mem(:,:,b,:,it,:) = self%patch%mem(:,:,b+1,:,it,:)
  endif
END SUBROUTINE mirror

!===============================================================================
SUBROUTINE extrap(self)
  class(boundaries_t):: self
  integer :: i, j, k, iv, iw, b, it, ng
  real(4) :: xxold(3), xxnew(self%patch%ng(3))
  real(4), dimension(:,:,:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  if (self%patch%ng(3) == 0) return
  it = self%patch%it
  ng = self%patch%ng(3)
  mem => self%patch%mem

  if (self%is_set(bits%zu)) then
      b = self%patch%ui(3)
      xxold = [(i*self%patch%ds(3), i=1, 3)]
      xxnew = [(i*self%patch%ds(3), i=4, 3+ng)]
      do iw=1, self%patch%nw
      do iv=1, self%patch%nv
      do j=1,  self%patch%gn(2)
      do i=1,  self%patch%gn(1)
        call pchip_mono(xxold, mem(i,j,b-2:b,   iv,it,iw), &
                        xxnew, mem(i,j,b+1:b+ng,iv,it,iw), sorted=.true.)
        mem(i,j,b+1:b+ng,iv,it,iw) = max(mem(i,j,b+1:b+ng,iv,it,iw), 0.)
      end do
      end do
      end do
      end do
  endif

  if (self%is_set(bits%zl)) then
      b = self%patch%li(3)
      xxold = [(i*self%patch%ds(3), i=1, 3)]
      xxnew = [(i*self%patch%ds(3), i=1-ng, 0)]
      do iw=1, self%patch%nw
      do iv=1, self%patch%nv
      do j=1,  self%patch%gn(2)
      do i=1,  self%patch%gn(1)
        call pchip_mono(xxold, mem(i,j,b:b+2,   iv,it,iw), &
                        xxnew, mem(i,j,b-ng:b-1,iv,it,iw), sorted=.true.)
        mem(i,j,b-ng:b-1,iv,it,iw) = max(mem(i,j,b-ng:b-1,iv,it,iw), 0.)
      end do
      end do
      end do
      end do
  endif
END SUBROUTINE extrap

END MODULE boundaries_mod
