!> =============================================================================
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE mhd_mod
  USE initial_mod
  USE timestep_mod
  USE atmosphere_mod
  USE units_mod
  USE scaling_mod
  USE planet_mod
  USE mesh_mod
  USE boundary_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
    type(initial_t):: initial
  contains
    procedure init
    procedure update
  end type
  integer(8), public:: n_update=0
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  call scaling%init
  !----------------------------------------------------------------------------
  ! Allocate memory and set IC values
  !----------------------------------------------------------------------------
  self%nv = 8
  self%msplit = 3
  call self%initial%init
  call self%patch_t%init
  call self%set(bits%static)
  call self%initial%selected(self%mem, self%mesh)
  call self%mhd_t%init
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Impose boundary conditions, taking staggering into account
!> =============================================================================
SUBROUTINE boundary_condition (self, mem, h)
  class(experiment_t)    :: self
  real, dimension(:,:,:) :: mem
  real(8)                :: h(3)
  !.............................................................................
  type(mesh_t) :: m(3)
  integer      :: ix, iy, iz, iv, nr
  real(8)      :: x, y, z, r
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%boundary_conditions')
  nr = atmosphere%nr
  do iv=1,self%nv
    do iz=m(3)%lb,m(3)%ub
      z = m(3)%p + m(3)%r(iz) + h(3) - planet%position(3)
      do iy=m(2)%lb,m(2)%ub
        y = m(2)%p + m(2)%r(iy) + h(2) - planet%position(2)
        do ix=m(1)%lb,m(1)%ub
          x = m(1)%p + m(1)%r(ix) + h(1) - planet%position(1)
          r = sqrt(x**2+y**2+z**2)
          if (r < planet%radius) then
            mem(ix,iy,iz) = 0.0
          end if
        end do
      end do
    end do
  end do
  call trace_end
END SUBROUTINE boundary_condition

!> =============================================================================
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%mhd_t%update
  call boundary%condition (self%mem(:,:,:,1,self%it,2), self%mesh, [-0.5d0,0d0,0d0], self%nv)
  call boundary%condition (self%mem(:,:,:,2,self%it,2), self%mesh, [0d0,-0.5d0,0d0], self%nv)
  call boundary%condition (self%mem(:,:,:,3,self%it,2), self%mesh, [0d0,0d0,-0.5d0], self%nv)
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%initial%selected(self%mem, self%mesh)
  !$omp atomic
  n_update = n_update + product(self%n)
  call trace_end
END SUBROUTINE update

END MODULE
