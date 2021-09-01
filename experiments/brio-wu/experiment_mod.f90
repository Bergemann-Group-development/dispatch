!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!=======================================================================
!> Shu-Osher shock tube problem.
!> Eample 8 in https://doi.org/10.1016/0021-9991(89)90222-2
!=======================================================================
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE bits_mod
  USE mesh_mod
  USE kinds_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
    procedure:: bndries
  end type

CONTAINS

!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  logical:: first_time = .true.
  class(mesh_t), pointer:: m(:)
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  self%mhd = .true. 
  call self%solver_t%init
  
  self%unsigned(self%idx%d) = .true.
  ! assign boundary conditions.
  m=>self%mesh
  self%boundaries%rho_l = self%mem(m(1)%lo,m(2)%li,m(3)%li,self%idx%d,1,1)
  self%boundaries%rho_r = self%mem(m(1)%uo,m(2)%li,m(3)%li,self%idx%d,1,1)
  if (self%idx%e>0) then
  self%boundaries%e_l   = self%mem(m(1)%lo,m(2)%li,m(3)%li,self%idx%e,1,1)
  self%boundaries%e_r   = self%mem(m(1)%uo,m(2)%li,m(3)%li,self%idx%e,1,1)
  end if
  if (self%idx%s>0) then
  self%boundaries%s_l   = self%mem(m(1)%lo,m(2)%li,m(3)%li,self%idx%s,1,1)
  self%boundaries%s_r   = self%mem(m(1)%uo,m(2)%li,m(3)%li,self%idx%s,1,1)
  end if
  self%boundaries%px_l  = self%mem(m(1)%lo,m(2)%li,m(3)%li,self%idx%px,1,1)
  self%boundaries%px_r  = self%mem(m(1)%uo,m(2)%li,m(3)%li,self%idx%px,1,1)
  self%boundaries%py_l  = self%mem(m(1)%lo,m(2)%li,m(3)%li,self%idx%py,1,1)
  self%boundaries%py_r  = self%mem(m(1)%uo,m(2)%li,m(3)%li,self%idx%py,1,1)
  self%boundaries%pz_l  = self%mem(m(1)%lo,m(2)%li,m(3)%li,self%idx%pz,1,1)
  self%boundaries%pz_r  = self%mem(m(1)%uo,m(2)%li,m(3)%li,self%idx%pz,1,1)

  self%boundaries%bx_l  = self%mem(m(1)%lo,m(2)%li,m(3)%li,self%idx%bx,1,1)
  self%boundaries%bx_r  = self%mem(m(1)%uo,m(2)%li,m(3)%li,self%idx%bx,1,1)
  self%boundaries%by_l  = self%mem(m(1)%lo,m(2)%li,m(3)%li,self%idx%by,1,1)
  self%boundaries%by_r  = self%mem(m(1)%uo,m(2)%li,m(3)%li,self%idx%by,1,1)
  self%boundaries%bz_l  = self%mem(m(1)%lo,m(2)%li,m(3)%li,self%idx%bz,1,1)
  self%boundaries%bz_r  = self%mem(m(1)%uo,m(2)%li,m(3)%li,self%idx%bz,1,1)

  call trace_end
END SUBROUTINE init

!===============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).  The next step
!> is to call a timestep%udpate, whose actual implementation is in a file
!> $(TOP)/timestep/$(UPD)/timestep_mod.f90 -- currently an Adams-Bashforth
!> method.
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: iv,i1,i2
  logical :: first_time=.true.
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  
  call self%solver_t%update

  call trace_end
END SUBROUTINE update

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:), pointer, contiguous:: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%bndries')
  mem => self%mem(:,:,:,:,self%it,1)
  call self%boundaries%primary (mem)
  call trace%end()
END SUBROUTINE bndries
!===============================================================================
END MODULE experiment_mod
