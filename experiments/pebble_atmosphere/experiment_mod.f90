!> =============================================================================
!> Setup the initial conditions and experiment update
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE pebble_patch_mod
  USE initial_mod
  !USE timestep_mod
  USE atmosphere_mod
  USE units_mod
  USE scaling_mod
  USE planet_mod
  USE mesh_mod
  USE boundary_mod
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
!> The initialization must be done in the right order.  The rotation speed
!> omega must be known before the initial%init runs, for example.
!>
!> dispatch
!>   setup::init
!>     scaling::init                            ! code units
!>     planet::init                             ! planet size, orbit
!>     atmosphere::init                         ! construct hydrostatic state
!>     rubiks_mod::init                         ! makes all patches
!>       experiment_mod::init                   ! calls other inits
!>         mhd_mod::init                        ! reads gamma
!>           force_mod::init                    ! initialize forcing
!>         pebbles_mod::init                    ! reads omega
!>         initial_mod::init                    ! selects ICs
!>         initial_mod::selected                ! computes ICs
!>         boundary_mod::condition              ! apply BCs to ICs
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real, dimension(:,:,:,:), pointer:: ff
  integer:: iv
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Allocate memory and set IC values
  !----------------------------------------------------------------------------
  self%nv = 5
  self%msplit = 3
  call scaling%init
  self%nv = 5
  self%mhd = .false.
  call self%pebble_patch_t%init
  self%mhd = .false.
  self%nv = 5
  call self%set (bits%static)
  call self%initial%init
  ff => self%mem(:,:,:,:,1,1)
  call self%initial%selected (self%mesh, ff)
  do iv=self%i_px,self%i_pz
    call boundary%fixed (self%mem(:,:,:,iv,self%it,1), 0.0, self%mesh, iv)
  end do
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Choice of update steps for this experiment
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: mem
  real, dimension(:,:,:), pointer:: d, e
  real(8), dimension(3):: h
  integer:: it, iv
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%check_density ('before update')
  !-----------------------------------------------------------------------------
  ! Option for transporting pebbles with frozen MHD values, in this case, make
  ! sure to wrap around in the time buffer before freezing, so the value needed
  ! by the pebbles have been computed and stored everywhere.
  !-----------------------------------------------------------------------------
  if (self%only_pebbles .or. self%is_set(bits%not_leaf)) then
    if (self%istep<self%nt+1) then
      call self%mhd_t%update
    end if
    if (self%is_set (bits%has_BC)) then
      do iv=self%i_px,self%i_pz
        call boundary%fixed (self%mem(:,:,:,iv,it,2), 0.0, self%mesh, iv)
      end do
    end if
    self%mem(:,:,:,:,self%it,2) = 0.0
  else
    if (self%is_set (bits%has_BC)) then
      allocate (mem(size(self%mem,1),size(self%mem,2),size(self%mem,3),2))
      mem(:,:,:,1) = self%mem(:,:,:,self%i_d,self%it,1)
      mem(:,:,:,2) = self%mem(:,:,:,self%i_e,self%it,1)
      !print *,self%id,' fixed BC 1', mem(12,12,12,:), self%it, self%nw
    end if
    call self%mhd_t%update
    if (self%is_set (bits%has_BC)) then
      if (self%nw > 1) then
        it = self%it
        do iv=1,self%nv
          call boundary%fixed (self%mem(:,:,:,iv,it,2), 0.0, self%mesh, iv)
        end do
      else
        it = self%new
      end if
      d => mem(:,:,:,1)
      e => mem(:,:,:,2)
      call boundary%fixed (self%mem(:,:,:,self%i_d,it,1), d, self%mesh, self%i_d)
      call boundary%fixed (self%mem(:,:,:,self%i_e,it,1), e, self%mesh, self%i_e)
      do iv=self%i_px,self%i_pz
        call boundary%fixed (self%mem(:,:,:,iv,it,1), 0.0, self%mesh, iv)
      end do
      deallocate (mem)
    end if
  end if
  call self%courant_condition
  call self%pebble_patch_t%update_pebbles
  !if (self%is_set (bits%has_BC)) &
    !print *,self%id,' fixed BC 2', self%mem(12,12,12,1:2,it,1), it
  !call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  !if (self%is_set (bits%has_BC)) &
    !print *,self%id,' fixed BC 3', self%mem(12,12,12,1:2,it,1), it
  call self%check_density ('after timestep')
  call trace_end
END SUBROUTINE update

END MODULE
