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
           P_r, rho_r, e1_r, px_r, py_r, pz_r, Bx_r, By_r, Bz_r, e_r, s_r
    character(32):: solver
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
  class(boundaries_t)     :: self
  !.............................................................................
  call self%standard_boundaries_t%init
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
      mem(i,:,:,idx%d) = self%rho_l
      mem(i,:,:,idx%px) = self%px_l
      mem(i,:,:,idx%py) = 0.0
      mem(i,:,:,idx%pz) = 0.0
          
      if (idx%e>0) then
        mem(i,:,:,idx%e) = self%e_l
      else if (idx%s>0) then
        mem(i,:,:,idx%s) = self%s_l
      end if
      if (associated(self%patch%filled)) then
      self%patch%filled(i,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%lo) = 1
    end if 
    end do  
    if (associated(self%patch%filled)) then
      self%patch%filled(m(1)%lb:m(1)%lo,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%lo) = 1
    end if 
  end if
  !
  if (self%mesh(1)%upper_boundary) then
    do i=m(1)%uo,m(1)%gn
      mem(i,:,:,idx%d) = self%rho_r
      mem(i,:,:,idx%px) = self%px_r
      mem(i,:,:,idx%py) = 0.0
      mem(i,:,:,idx%pz) = 0.0
          
      if (idx%e>0) then
        mem(i,:,:,idx%e) = self%e_r
      else if (idx%s>0) then
        mem(i,:,:,idx%s) = self%s_r
      end if
      if (associated(self%patch%filled)) then
      self%patch%filled(i,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%lo) = 1
    end if
    end do
  end if
  !
  if (self%mesh(1)%lower_boundary) then
    do i=m(1)%lb,m(1)%lo
      mem(i,:,:,idx%d) = self%rho_l
      mem(i,:,:,idx%px) = self%px_l
      mem(i,:,:,idx%py) = 0.0
      mem(i,:,:,idx%pz) = 0.0
          
      if (idx%e>0) then
        mem(i,:,:,idx%e) = self%e_l
      else if (idx%s>0) then
        mem(i,:,:,idx%s) = self%s_l
      end if
      if (associated(self%patch%filled)) then
      self%patch%filled(i,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%lo) = 1
    end if 
    end do  
  end if
  !
  if (self%mesh(1)%upper_boundary) then
    do i=m(1)%uo,m(1)%gn
      mem(i,:,:,idx%d) = self%rho_r
      mem(i,:,:,idx%px) = self%px_r
      mem(i,:,:,idx%py) = 0.0
      mem(i,:,:,idx%pz) = 0.0
          
      if (idx%e>0) then
        mem(i,:,:,idx%e) = self%e_r
      else if (idx%s>0) then
        mem(i,:,:,idx%s) = self%s_r
      end if
      if (associated(self%patch%filled)) then
      self%patch%filled(i,m(2)%lb:m(2)%ub,m(3)%lb:m(3)%ub) = 1
    end if
    end do
  end if
  !
  if (self%mesh(2)%lower_boundary) then
    do i=m(2)%lb,m(2)%lo
      mem(:,i,:,idx%d) = self%rho_l
      mem(:,i,:,idx%px) = self%px_l
      mem(:,i,:,idx%py) = 0.0
      mem(:,i,:,idx%pz) = 0.0
          
      if (idx%e>0) then
        mem(:,i,:,idx%e) = self%e_l
      else if (idx%s>0) then
        mem(:,i,:,idx%s) = self%s_l
      end if
      if (associated(self%patch%filled)) then
      self%patch%filled(m(1)%lb:m(1)%ub,i,m(3)%lb:m(3)%ub) = 1
    end if 
    end do  
  end if
  !
  if (self%mesh(2)%upper_boundary) then
    do i=m(2)%uo,m(2)%gn
      mem(:,i,:,idx%d) = self%rho_r
      mem(:,i,:,idx%px) = self%px_r
      mem(:,i,:,idx%py) = 0.0
      mem(:,i,:,idx%pz) = 0.0
          
      if (idx%e>0) then
        mem(:,i,:,idx%e) = self%e_r
      else if (idx%s>0) then
        mem(:,i,:,idx%s) = self%s_r
      end if
      if (associated(self%patch%filled)) then
      self%patch%filled(m(1)%lb:m(1)%ub,i,m(3)%lb:m(3)%ub) = 1
    end if
    end do
  end if
  !
  if (self%mesh(3)%lower_boundary) then
    do i=m(3)%lb,m(3)%lo
      mem(:,:,i,idx%d) = self%rho_l
      mem(:,:,i,idx%px) = self%px_l
      mem(:,:,i,idx%py) = 0.0
      mem(:,:,i,idx%pz) = 0.0
          
      if (idx%e>0) then
        mem(:,:,i,idx%e) = self%e_l
      else if (idx%s>0) then
        mem(:,:,i,idx%s) = self%s_l
      end if
      if (associated(self%patch%filled)) then
      self%patch%filled(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,i) = 1
    end if 
    end do  
  end if
  !
  if (self%mesh(3)%upper_boundary) then
    do i=m(3)%uo,m(3)%gn
      mem(:,:,i,idx%d) = self%rho_r
      mem(:,:,i,idx%px) = self%px_r
      mem(:,:,i,idx%py) = 0.0
      mem(:,:,i,idx%pz) = 0.0
          
      if (idx%e>0) then
        mem(:,:,i,idx%e) = self%e_r
      else if (idx%s>0) then
        mem(:,:,i,idx%s) = self%s_r
      end if
      if (associated(self%patch%filled)) then
      self%patch%filled(m(1)%lb:m(1)%ub,m(2)%lb:m(2)%ub,i) = 1
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
