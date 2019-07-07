!===============================================================================
!> This module contains all experiment specific information.  The experiment_t
!> object extends the generic MHD object mhd_t, which is again an extension of
!> the basic patch_t object for tasks (object task_t) with meshes.
!>
!> The MHD object is chosen with the SOLVER macro in the Makefile in the curren
!> directory, while the generic task_t and patch_t objects are defined in the
!> $(TOP)/tasks/ directory, in task_mod.f90 and patch_mod.f90.
!>
!> The init procedure here (experiment_t::init) is called from the cartesian_mod
!> module, which is specified as the patch distributor in the main program
!> dispatch.f90.   The update procedure here is called by the task_list_mod
!> update procedure in $(TOP)/lists/task_list_mod.f90.
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE initial_mod
  USE boundaries_mod
  USE kinds_mod
  USE eos_mod
  USE scaling_mod
  USE scalar_mod
  USE timestep_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure init
    procedure update
    !procedure q_diffuse
  end type
  logical, save:: first_time=.true.
  real, save:: temp0=0.0
  real, save:: me_per_year=0.0
  real, save:: epsilon=1.0
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real, dimension(:,:,:,:), pointer :: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%mhd = .false.
  call self%solver_t%init
  self%unsigned = .false.
  !----------------------------------------------------------------------------
  ! Boundary setup
  !----------------------------------------------------------------------------
  call self%init_bdries
  call self%boundaries%init(mesh=self%mesh)
  self%boundaries%id = self%id
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
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: mem
  integer:: iv
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  ! -- apply BCs only to the parent MHD task
  do iv=1, self%nv
    mem => self%mem(:,:,:,iv,self%it,1)
    call self%boundaries%extrapolate(mem, iv)
  end do
  call self%solver_t%update
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> If RT is on, add the radiative cooling / heating to the entropy time derivative
!===============================================================================
!SUBROUTINE q_diffuse (self, q)
!  class (experiment_t) :: self
!  real(kind=KindScalarVar), dimension(:,:,:), target  :: q
!  real(kind=KindScalarVar), dimension(:,:,:), allocatable :: fmax
!  real(kind=KindScalarVar), dimension(:,:,:), pointer :: dedt, d, dsdt
!  real:: u_max, rk_min
!  integer:: ix
!  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
!  if (eos%do_rt.and.associated(self%rt)) then
!    call trace%begin ('experiment_t%q_diffuse', itimer=itimer)
!    self%rt%epsilon = epsilon
!    call self%rt%update
!    associate(l => self%li, u => self%ui)
!    dedt => q(l(1):u(1),l(2):u(2),l(3):u(3))
!    call self%rt%get_q (self%time, dedt)
!    end associate
!    !---------------------------------------------------------------------------
!    ! Time step limit
    !---------------------------------------------------------------------------
!    rk_min=minval(self%rt%rk)
!    allocate(fmax(self%gn(1), self%gn(2), self%gn(3)))
!    print *, 'pg, d, pg/d max', maxval(self%pg), maxval(self%d), minval(self%d), maxval(self%pg/self%d)
!    print *, 'pg/d 4', maxval((self%pg/self%d)**4)
!    fmax = (self%pg/self%d)**4/self%pg*(rk_min*minval(self%ds))
!    fmax = fmax/(1.+(rk_min*minval(self%ds))**2)
!    u_max = scaling%stefan*16./3.*self%fmaxval(fmax)*self%courant/self%cdtd*epsilon
!    deallocate(fmax)
!    if (io%debug(2)) print 1, self%id, 'u_max(RT)', u_max, self%u_max
!    1 format(i6,3x,a,1p,2g12.3)
!    self%u_max = max(self%u_max,u_max)
!    call trace%end (itimer)
!  end if
!END SUBROUTINE q_diffuse

END MODULE experiment_mod
