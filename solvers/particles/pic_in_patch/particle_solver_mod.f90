!==============================================================================
!> Particle data type for PIC solvers
!===============================================================================
!define DBG
MODULE particle_solver_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE species_mod
  USE fields_mod
  USE kinds_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Everything needed by a particle solver
  !-----------------------------------------------------------------------------
  type, public, extends(task_t):: particle_solver_t
    class(patch_t),    pointer :: patch         => null() ! pointer to patch
    type(fields_t)             :: fields                  ! particles to move
    type(species_t)            :: species                 ! task link
    logical                    :: on = .false.            ! on/off
    logical                    :: initialized = .false.   ! ICs done
    integer, dimension(3)      :: lb = 0,  ub = 0         ! lower & upper bounds
    real(var_k), pointer, contiguous, dimension(:,:,:,:):: B, E
    !---------------------------------------------------------------------------
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: initial_condition
    procedure:: pre_update
    procedure:: post_update
    procedure:: update
  end type
  integer:: verbose = 0
  logical:: on = .false.
  logical:: first_time = .true.
CONTAINS

!===============================================================================
!> Initialize particle solver, reading parameters and adding task link
!===============================================================================
SUBROUTINE pre_init (self)
  class(particle_solver_t), target:: self
  class(patch_t), pointer:: patch
  integer:: iostat
  namelist /particle_solver_params/ on, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pre_init')
  patch => self%patch
  patch%nv = 11
  patch%mv = 11
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read   (io%input, particle_solver_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read   (io%input, particle_solver_params)
    end if
    call io%namelist_warning ('particle_solver_params', iostat)
    write (io%namelist, particle_solver_params)
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
 call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize particle solver, after patch is initialized
!===============================================================================
SUBROUTINE init (self)
  class(particle_solver_t)    :: self
  class(patch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%init')
  !-----------------------------------------------------------------------------
  patch => self%patch
  self%lb = patch%li
  self%ub = patch%ui
  self%fields%lb = self%lb
  self%fields%ub = self%ub
  self%species%lb = self%lb
  self%species%ub = self%ub
  call self%fields%pre_init (patch)
  call self%species%pre_init (patch)
  call self%fields%init (patch)
  call self%species%init (patch)
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize particle solver, after patch is initialized
!===============================================================================
SUBROUTINE initial_condition (self)
  class(particle_solver_t)    :: self
  !-----------------------------------------------------------------------------
  if (.not. on .or. self%initialized) return
  self%initialized = .true.
  call trace%begin ('particle_solver_t%initial_condition')
  !-----------------------------------------------------------------------------
  ! call self%fields%initial_condition
  ! call self%species%initial_condition
  call trace%end ()
END SUBROUTINE initial_condition

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE pre_update (self)
  class(particle_solver_t), target:: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%pre_update', itimer=itimer)
  ! call self%fields%pre_update
  ! call self%species%pre_update
  call self%update ()
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE post_update (self)
  class(particle_solver_t):: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%post_update', itimer=itimer)
  ! call self%fields%post_update
  ! call self%species%post_update
  call trace%end (itimer)
END SUBROUTINE post_update

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE update (self)
  class(particle_solver_t):: self
  class(patch_t), pointer:: patch
  integer:: m(6)
  real(var_k), pointer, contiguous:: J(:,:,:,:)
  real(var_k), pointer, contiguous, dimension(:,:,:,:):: B0, B1, E0, E1
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  patch => self%patch
  m = shape(patch%mem)
  B0 => patch%mem(:,:,:,patch%idx%bx:patch%idx%bz,patch%it ,1)
  E0 => patch%mem(:,:,:,patch%idx%ex:patch%idx%ez,patch%it ,1)
  B1 => patch%mem(:,:,:,patch%idx%bx:patch%idx%bz,patch%new,1)
  E1 => patch%mem(:,:,:,patch%idx%ex:patch%idx%ez,patch%new,1)
  m = shape(self%patch%mem)
  allocate (J(3,m(1),m(2),m(3)))
  call self%species%update (B0, E0, J)
  call self%fields%update  (J)
  deallocate (J)
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE update

END MODULE particle_solver_mod
