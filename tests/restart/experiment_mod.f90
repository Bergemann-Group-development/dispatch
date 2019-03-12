!===============================================================================
!> Setup the initial conditions and experiment update
!===============================================================================
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
  USE boundaries_mod
  implicit none
  private
  type, public, extends(pebble_patch_t):: experiment_t
  contains
    procedure init
    procedure update
    procedure output
    procedure edge_boundary_conditions
  end type
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
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
!>         initial_mod::condition              ! computes ICs
!>         boundaries_mod::condition            ! apply BCs to ICs
!===============================================================================
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
  call self%initial%init (self%kind)
  ff => self%mem(:,:,:,:,1,1)
  call self%initial%condition (self%mesh, ff, self%idx, self%idx, self%level)
  call self%boundaries%init (self%mesh, position=planet%position, radius=planet%radius)
  do iv=self%idx%px,self%idx%pz
    call self%boundaries%spherical (self%mem(:,:,:,iv,self%it,1), 0.0, self%mesh, iv)
  end do
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Choice of update steps for this experiment
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: mem
  real, dimension(:,:,:), pointer:: d, e
  real(8), dimension(3):: h
  integer:: it, iv, nb
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%check_density ('before update')
  !-----------------------------------------------------------------------------
  ! Option for transporting pebbles with frozen MHD values, in this case, make
  ! sure to wrap around in the time buffer before freezing, so the value needed
  ! by the pebbles have been computed and stored everywhere.
  !-----------------------------------------------------------------------------
  !call self%boundaries%init (self%mesh, position=planet%position, radius=planet%radius)
  if (self%only_pebbles .or. self%is_set(bits%not_leaf)) then
    if (self%istep<self%nt+1) then
      call self%mhd_t%update
    end if
    if (self%is_set (bits%spherical)) then
      do iv=self%idx%px,self%idx%pz
        call self%boundaries%spherical (self%mem(:,:,:,iv,it,2), 0.0, self%mesh, iv)
      end do
    end if
    self%mem(:,:,:,:,self%it,2) = 0.0
  !-----------------------------------------------------------------------------
  ! Normal time stepping, with boundary conditions
  !-----------------------------------------------------------------------------
  else
    !------------------------------------------------------------------------
    ! If spherical inner BCs, save the incoming values
    !------------------------------------------------------------------------
    if (self%boundaries%is_set (bits%spherical)) then
      allocate (mem(size(self%mem,1),size(self%mem,2),size(self%mem,3),2))
      mem(:,:,:,1) = self%mem(:,:,:,self%idx%d,self%it,1)
      mem(:,:,:,2) = self%mem(:,:,:,self%idx%e,self%it,1)
    end if
    !--------------------------------------------------------------------------
    ! Update the MHD state
    !--------------------------------------------------------------------------
    call self%mhd_t%update
    if (self%boundaries%status /= 0) then
      print *, self%id, 'BOUNDARY CONDITION BITS:', self%boundaries%status
      !------------------------------------------------------------------------
      ! If spherical inner BCs, save the incoming values
      !------------------------------------------------------------------------
      if (self%boundaries%is_set (bits%spherical)) then
        print *, self%id, 'BOUNDARY CONDITION: spherical'
        if (self%nw > 1) then
          it = self%it
          do iv=1,self%nv
            call self%boundaries%spherical (self%mem(:,:,:,iv,it,2), 0.0, self%mesh, iv)
          end do
        else
          it = self%new
        end if
        d => mem(:,:,:,self%idx%d)
        e => mem(:,:,:,self%idx%e)
        call self%boundaries%init (radius=planet%radius-2.*self%ds(1))
        call self%boundaries%spherical_extrapolate (self%mem(:,:,:,self%idx%d,it,1), d, self%mesh, self%idx%d)
        call self%boundaries%spherical_extrapolate (self%mem(:,:,:,self%idx%e,it,1), e, self%mesh, self%idx%e)
        call self%boundaries%init (radius=planet%radius)
        do iv=self%idx%px,self%idx%pz
          call self%boundaries%spherical (self%mem(:,:,:,iv,it,1), 0.0, self%mesh, iv)
        end do
        deallocate (mem)
      end if
      !------------------------------------------------------------------------
      ! Edge boundary conditions
      !------------------------------------------------------------------------
      if (self%boundaries%is_set (bits%xdir + bits%ydir + bits%zdir)) then
        nb = 0
        if (self%boundaries%is_set (bits%xdir)) nb = nb+1
        if (self%boundaries%is_set (bits%ydir)) nb = nb+1
        if (self%boundaries%is_set (bits%zdir)) nb = nb+1
        if (nb == 1) print *, self%id, 'BOUNDARY CONDITION: face'
        if (nb == 2) print *, self%id, 'BOUNDARY CONDITION: edge'
        if (nb == 3) print *, self%id, 'BOUNDARY CONDITION: corner'
        call self%edge_boundary_conditions
      end if
    end if
  end if
  call self%courant_condition
  call self%pebble_patch_t%update_pebbles
  call self%check_density ('after timestep')
  call trace_end
END SUBROUTINE update

SUBROUTINE output (self, experiment_name)
  class(experiment_t):: self
  character(len=64), optional:: experiment_name
  integer:: iv, it
  if (self%id==5) then
    it = self%iit(self%nt-2)
    do iv=2,4
      print *,'OUT', real(self%t(it)), iv, self%mem(18,18,1:4,iv,it,1)/self%mem(18,18,1:4,1,it,1)
    end do
  end if
  call self%patch_t%output
END SUBROUTINE

!===============================================================================
!> Boundary conditions at the edges of the outermost rubiks cubes; extrapolate
!> the density and pressure, and enforce vanishing derivative on the perpendicular
!> velocities. To what extend the boundary conditions are applied in the x,y,z
!> directions depends on the self%boundary%status bits {xdir,ydir,zdir}.
!===============================================================================
SUBROUTINE edge_boundary_conditions (self)
  class(experiment_t):: self
  real, dimension(:,:,:), pointer:: d, e, p, u
  integer:: iv
  !-----------------------------------------------------------------------------
  ! Edge boundary conditions
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%edge_boundary_conditions')
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  e => self%mem(:,:,:,self%idx%e,self%it,1)
  allocate (u(size(d,1),size(d,2),size(d,3)))
  call self%boundaries%extrapolate (d)
  u = e/d
  call self%boundaries%symmetric (u)
  e = u*d
  do iv=self%idx%px,self%idx%pz
    p => self%mem(:,:,:,iv,self%it,1)
    u = p/d
    call self%boundaries%symmetric (u)
    p = u*d
  end do
  deallocate (u)
  call trace%end
END SUBROUTINE edge_boundary_conditions

END MODULE
