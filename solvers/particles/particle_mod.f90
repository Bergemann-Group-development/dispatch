!===============================================================================
!> Particle data type, extends a dll_node data type, so it can become part of
!> a particles_t data type
!===============================================================================
MODULE particle_mod
  USE iso_fortran_env, only: int8
  USE io_unit_mod
  USE trace_mod
  USE dll_mod
  USE patch_forces_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: particle_t
    integer:: id=0
    integer:: it=1
    real:: ds=1.0
    real:: mass=1.0
    real(8):: time=0d0, dtime=0d0
    integer(kind=int8), allocatable:: iit(:)
    real, allocatable:: v(:,:)
    real(8), allocatable:: r(:,:), t(:), m(:)
    real(8):: f(3)
    type(patch_forces_t), pointer:: patch_forces
  contains
    procedure:: init
    procedure:: update
  end type
  integer, save:: id=0
  integer, save:: verbose=0
  integer, save:: nt=4
  type(particle_t), public:: particle
  !-----------------------------------------------------------------------------
CONTAINS

!===============================================================================
!> Read parameters, and initialize a particle with a unique ID
!===============================================================================
SUBROUTINE init (self)
  class (particle_t):: self
  logical, save:: first_time=.true.
  integer:: iostat
  namelist /particle_params/ verbose, nt
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_t%init')
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io_unit%input)
      read (io_unit%input, particle_params,iostat=iostat)
      write (io_unit%output, particle_params)
    end if
    !$omp end critical (input_cr)
  end if
  !$omp atomic
  id = id+1
  if (verbose > 1) &
    print *, 'particle_t%init: id =', id
  self%id = id
  allocate (self%r(3,nt), self%v(3,nt), self%m(nt), self%t(nt), self%iit(nt))
  self%r  = 0.0_8
  self%v  = 0.0_8
  self%f  = 0.0_8
  self%t  = 0.0_8
  self%it = 1
  self%iit = 1
  allocate (self%patch_forces)
  call self%patch_forces%init
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Read parameters, and initialize a particle with a unique ID
!===============================================================================
SUBROUTINE update (self)
  class (particle_t):: self
END SUBROUTINE update

END MODULE particle_mod
