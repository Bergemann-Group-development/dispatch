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
  USE patch_mod
  !
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    real::   e_l = 10.33333
    real::   s_l = 10.33333
    real:: rho_l =  3.857143
    real::  px_l =  10.141852233 !2.629369 * 3.857143
    real:: vy_l, vz_l, Bx_l, By_l, Bz_l, &
           P_r, rho_r, e1_r, px_r, py_r, pz_r, Bx_r, By_r, Bz_r
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: init
    procedure:: primary
    procedure:: derivative
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
  end type
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self, patch)
  class(boundaries_t)     :: self
  class(patch_t),target          :: patch
  !.............................................................................
  call self%standard_boundaries_t%init (patch)
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  integer:: i
  ! ............................................................................
   call trace_begin('boundaries_t%primary')
  associate(m => self%mesh)
  if (self%mesh(1)%lower_boundary) then
    do i=m(1)%lb,m(1)%lo
      mem(i,:,:,:) = mem(m(1)%li,:,:,:)
    end do
  end if
  !
  if (self%mesh(1)%upper_boundary) then
    do i=m(1)%uo,m(1)%gn
      mem(i,:,:,:) = mem(m(1)%ui,:,:,:)
    end do
  end if
  end associate
  call trace%end()
END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  !
  !if (self%mesh(3)%lower_boundary) call self%lower_damp(mem)

END SUBROUTINE derivative


END MODULE boundaries_mod
