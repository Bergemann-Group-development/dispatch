!> =============================================================================
!> Setup the initial conditions and experiment update
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE pebbles_mod
  USE initial_mod
  USE timestep_mod
  USE atmosphere_mod
  USE units_mod
  USE scaling_mod
  USE planet_mod
  USE mesh_mod
  USE boundary_mod
  USE force_mod
  implicit none
  private
  type, public, extends(pebble_patch_t):: experiment_t
    type(initial_t):: initial
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')
  call scaling%init
  !----------------------------------------------------------------------------
  ! Allocate memory and set IC values
  !----------------------------------------------------------------------------
  self%nv = 8
  self%msplit = 3
  call self%initial%init
  call self%mhd_t%init
  call self%patch_t%init
  call self%set(bits%static)
  call self%initial%selected(self%mem, self%mesh, self%omega)
  call force%init
  call boundary%condition (self%mem(:,:,:,3,self%it,1), self%mesh, [-0.5d0,0d0,0d0], self%nv)
  call boundary%condition (self%mem(:,:,:,4,self%it,1), self%mesh, [0d0,-0.5d0,0d0], self%nv)
  call boundary%condition (self%mem(:,:,:,5,self%it,1), self%mesh, [0d0,0d0,-0.5d0], self%nv)
  call self%pebble_patch_t%init
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Choice of update steps for this experiment
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !.............................................................................
  type(timestep_t):: timestep
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%check_density ('before update')
  !-----------------------------------------------------------------------------
  ! Option for transporting pebbles with frozen MHD values, in this case, make
  ! sure to wrap around in the time buffer before freezing, so the value needed
  ! by the pebbles have been computed and stored everywhere.
  !-----------------------------------------------------------------------------
  if (self%only_pebbles) then
    if (self%istep<self%nt+1) then
      call self%mhd_t%update
    end if
    self%mem(:,:,:,:,self%it,2) = 0.0
  else
    call self%mhd_t%update
    call boundary%condition (self%mem(:,:,:,1,self%it,2), self%mesh, [0d0,0d0,0d0], self%nv)
    call boundary%condition (self%mem(:,:,:,2,self%it,2), self%mesh, [0d0,0d0,0d0], self%nv)
    call boundary%condition (self%mem(:,:,:,3,self%it,2), self%mesh, [-0.5d0,0d0,0d0], self%nv)
    call boundary%condition (self%mem(:,:,:,4,self%it,2), self%mesh, [0d0,-0.5d0,0d0], self%nv)
    call boundary%condition (self%mem(:,:,:,5,self%it,2), self%mesh, [0d0,0d0,-0.5d0], self%nv)
  end if
  block
  !.............................................................................
  integer      :: ix, iy, iz, iv, nv
  real(8)      :: x, y, z, r, dr, damp
  type(mesh_t), pointer:: m(:)
  !-----------------------------------------------------------------------------
  if (.not.self%only_pebbles) then
    associate(m => self%mesh)
    do iv=3,5
      do iz=m(3)%lb,m(3)%ub
        z = m(3)%p + m(3)%r(iz) + m(3)%h(3)*m(3)%d - planet%position(3)
        do iy=m(2)%lb,m(2)%ub
          y = m(2)%p + m(2)%r(iy) + m(2)%h(2)*m(2)%d - planet%position(2)
          do ix=m(1)%lb,m(1)%ub
            x = m(1)%p + m(1)%r(ix) + m(1)%h(1)*m(1)%d - planet%position(1)
            r = sqrt(x**2+y**2+z**2)
            dr = r-planet%radius
            damp = exp(-dr/(0.1*planet%radius))
            if (dr>0.0) then
              self%mem(ix,iy,iz,iv,self%it,2) = &
              self%mem(ix,iy,iz,iv,self%it,2) - 300*damp* &
              self%mem(ix,iy,iz,iv,self%it,1)
            end if
            if (r-planet%radius < self%ds(1)*2) then
              self%mem(ix,iy,iz,iv,self%it,1) = 0.0
            end if
          end do
        end do
      end do
    end do
    end associate
  end if
  end block
  call self%courant_condition
  call self%pebble_patch_t%update_pebbles
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density ('after timestep')
  call trace_end
END SUBROUTINE update

END MODULE
