!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Alfven wave problem to test the hlld solver implementation
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE pboundary_mod
  USE bits_mod
  USE index_mod
  USE solver_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
  end type

CONTAINS

!> =============================================================================
!> Setup a patch with an Alfven wave (see Masson et al. 2012)
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  logical:: first_time = .true.
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Initial conditions
  !----------------------------------------------------------------------------
  self%mhd = .false.
  call self%solver_t%init

  call self%set (bits%static)
  if (self%initial%time > 0d0) then
    self%time = self%initial%time
    self%t(self%it) = self%initial%time
  end if
  self%gamma=1.4
  self%unsigned(self%idx%d) = .true.
  self%periodic = [.false., .true., .true.]
  ! assign boundary conditions.
  call pboundary%condition(self%mem(:,:,:,:,self%it,1), self%mesh, self%idx, .true.)
  call pboundary%condition(self%mem(:,:,:,:,self%it,1), self%mesh, self%idx, .false.)
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).  The next step
!> is to call a timestep%udpate, whose actual implementation is in a file
!> $(TOP)/timestep/$(UPD)/timestep_mod.f90 -- currently an Adams-Bashforth
!> method.
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: iv,i1,i2
  logical :: first_time=.true.
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  if ((self%llc_nat(1)-self%ds(1)) .lt. self%mesh(1)%origin) then
    call pboundary%condition (self%mem(:,:,:,:,self%it,1), self%mesh, self%idx,.true.)
  endif
  if ((self%llc_nat(1)+self%size(1)+self%ds(1)) .gt. (self%mesh(1)%origin + self%box(1))) then
    call pboundary%condition (self%mem(:,:,:,:,self%it,1), self%mesh, self%idx, .false.)
  endif
  
  call self%mhd_t%update

  call trace_end
END SUBROUTINE update
!> =============================================================================
END MODULE experiment_mod
