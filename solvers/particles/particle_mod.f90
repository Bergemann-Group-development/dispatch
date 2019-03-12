!===============================================================================
!> Particle data type, extends a dll_node data type, so it can become part of
!> a particle_list_t data type
!===============================================================================
MODULE particle_mod
  USE iso_fortran_env, only: int8
  USE io_unit_mod
  USE dll_mod
  implicit none
  private
  integer, save:: id=0
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: particle_t
    integer(8):: id=0
    integer:: it=1
    real:: ds=1.0
    real:: mass=1.0
    real(8):: time=0d0, dtime=0d0
    integer(kind=int8), allocatable:: iit(:)
    real, allocatable:: v(:,:)
    real(8), allocatable:: r(:,:), t(:)
  contains
    procedure:: init
  end type
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
  namelist /particle_params/ verbose, nt
  !-----------------------------------------------------------------------------
  if (first_time) then
    !$omp critical (input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io_unit%input)
      read (io_unit%input, particle_params)
      write (io_unit%output, particle_params)
    end if
    !$omp end critical (input_cr)
  end if
  !$omp atomic
  id = id+1
  self%id = id
  allocate (self%r(3,nt), self%v(3,nt), self%t(nt), self%iit(nt))
  self%r  = 0.0_8
  self%v  = 0.0_8
  self%t  = 0.0_8
  self%iit = 1
END SUBROUTINE init

END MODULE particle_mod
