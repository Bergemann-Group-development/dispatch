!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the C-shock problem to test ambipolar diffusion
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE mpi_mod
  USE kinds_mod
  USE mesh_mod
  USE index_mod
  USE solver_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
  end type
  logical, save:: first_time = .true.
CONTAINS

!===============================================================================
!> A series of single point perturbation tests
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:)  , pointer:: d, e, px, py, pz, bx, by, bz
  real(var_k), dimension(:,:,:,:), pointer:: p, b
  integer:: i, j, m
  !-----------------------------------------------------------------------------
  ! Initialize solver and define variables
  !-----------------------------------------------------------------------------
  call trace_begin ('experiment_t%init')
  call self%solver_t%init
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  e  => self%mem(:,:,:,self%idx%e ,self%it,1)
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  py => self%mem(:,:,:,self%idx%py,self%it,1)
  pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  bx => self%mem(:,:,:,self%idx%bx,self%it,1)
  by => self%mem(:,:,:,self%idx%by,self%it,1)
  bz => self%mem(:,:,:,self%idx%bz,self%it,1)
  p  => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  b  => self%mem(:,:,:,self%idx%bx:self%idx%bz,self%it,1)
  !-----------------------------------------------------------------------------
  ! Test 1: Constant B test, with a single velocity point
  !-----------------------------------------------------------------------------
  d = 1.0
  e = 1.0
  self%gamma = 1.4
  self%dt_fixed = 0.01
  m = (self%gn(1)+1)/2
  do j=0,1
    do i=1,3
      print *,''
      print *,'====================== ',i,j,' ============================='
      p = 0.0
      p(m,m,m,mod(i+j,3)+1) = 1.0
      b = 0.0
      b(:,:,:,i) = 1.0
      !-------------------------------------------------------------------------
      ! The update() call causes output of the ICs, the %output writes results
      !-------------------------------------------------------------------------
      self%out_next = 0d0
      call self%update
      call self%rotate
      self%out_next = self%dt_fixed*0.99
      call self%output
    end do
  end do
  !-----------------------------------------------------------------------------
  stop
  call trace_end ()
END SUBROUTINE init

!===============================================================================
END MODULE experiment_mod
