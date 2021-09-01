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
  USE patch_mod
  USE standard_boundaries_mod
  !
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
    real::   e_l = 10.33333

    real:: Bx_l, By_l, Bz_l, &
           rho_l, e_r, px_l, py_l, pz_l, s_l
     real:: Bx_r, By_r, Bz_r, &
           rho_r, px_r, py_r, pz_r,  s_r
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
SUBROUTINE init (self)
  class(boundaries_t)              :: self
  !.............................................................................

  call self%standard_boundaries_t%init()
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  real(8):: dtime
  ! ----------------------------------------------------------------------------
  integer:: i
  ! ............................................................................
   call trace_begin('boundaries_t%primary')
  associate(m => self%mesh)
  if (self%mesh(1)%lower_boundary) then
    do i=m(1)%lb,m(1)%lo
      mem(i,:,:,idx%d) = self%rho_l
      mem(i,:,:,idx%px) = self%px_l
      mem(i,:,:,idx%py) = self%py_l
      mem(i,:,:,idx%pz) = self%pz_l
      mem(i,:,:,idx%bx) = self%bx_l
      mem(i,:,:,idx%by) = self%by_l
      mem(i,:,:,idx%bz) = self%bz_l
          
      if (idx%e>0) then
        mem(i,:,:,idx%e) = self%e_l
      else if (idx%s > 0) then
        mem(i,:,:,idx%s) = self%s_l
      end if
      if (associated(self%patch%filled)) then
        self%patch%filled(i,:,:) = 1
      end if
    end do  
    
  end if
  !
  if (self%mesh(1)%upper_boundary) then
    do i=m(1)%uo,m(1)%gn
      mem(i,:,:,idx%d) = self%rho_r
      mem(i,:,:,idx%px) = self%px_r
      mem(i,:,:,idx%py) = self%py_r
      mem(i,:,:,idx%pz) = self%pz_r
      mem(i,:,:,idx%bx) = self%bx_r
      mem(i,:,:,idx%by) = self%by_r
      mem(i,:,:,idx%bz) = self%bz_r
          
      if (idx%e>0) then
        mem(i,:,:,idx%e) = self%e_r
      else if (idx%s > 0) then
        mem(i,:,:,idx%s) = self%s_r
      end if
      if (associated(self%patch%filled)) then
        self%patch%filled(i,:,:) = 1
      end if
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
