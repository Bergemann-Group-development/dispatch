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
  USE boundaries_mod
  implicit none
  private
  type, public, extends(pebble_patch_t):: experiment_t
    type(initial_t):: initial
  contains
    procedure :: init
    procedure :: update
    procedure :: output
    procedure :: edge_boundary_conditions
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
!>         boundaries_mod::condition            ! apply BCs to ICs
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
  call self%initial%selected (self%mesh, ff, self%level)
  call self%boundaries%init (self%mesh, position=planet%position, radius=planet%radius)
  do iv=self%i_px,self%i_pz
    call self%boundaries%spherical (self%mem(:,:,:,iv,self%it,1), 0.0, self%mesh, iv)
  end do
  call self%output
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
  integer:: it, iv, nb,i
  !-----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%check_density ('before update')

  !-----------------------------------------------------------------------------
  ! Option for transporting pebbles with frozen MHD values, in this case, make
  ! sure to wrap around in the time buffer before freezing, so the value needed
  ! by the pebbles have been computed and stored everywhere.
  !-----------------------------------------------------------------------------
  call self%boundaries%init (self%mesh, position=planet%position, radius=planet%radius)
  if (self%only_pebbles .or. self%is_set(bits%not_leaf)) then
    if (self%istep<self%nt+1) then
      call self%mhd_t%update
    end if
    if (self%is_set (bits%spherical)) then
      do iv=self%i_px,self%i_pz
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
      mem(:,:,:,1) = self%mem(:,:,:,self%i_d,self%it,1)
      mem(:,:,:,2) = self%mem(:,:,:,self%i_e,self%it,1)
    end if
    !--------------------------------------------------------------------------
    ! Update the MHD state
    !--------------------------------------------------------------------------
    call self%mhd_t%update
    if (self%boundaries%status /= 0) then
    !  print *, self%id, 'BOUNDARY CONDITION BITS:', self%boundaries%status
      !------------------------------------------------------------------------
      ! If spherical inner BCs, save the incoming values
      !------------------------------------------------------------------------
      if (self%boundaries%is_set (bits%spherical)) then
      !  print *, self%id, 'BOUNDARY CONDITION: spherical'
        if (self%nw > 1) then
          it = self%it
          do iv=1,self%nv
            call self%boundaries%spherical (self%mem(:,:,:,iv,it,2), 0.0, self%mesh, iv)
          end do
        else
          it = self%new
        end if
        d => mem(:,:,:,1)
        e => mem(:,:,:,2)
        call self%boundaries%init (radius=planet%radius-2.*self%ds(1))
        call self%boundaries%spherical_extrapolate (self%mem(:,:,:,self%i_d,it,1), d, self%mesh, self%i_d)
        call self%boundaries%spherical_extrapolate (self%mem(:,:,:,self%i_e,it,1), e, self%mesh, self%i_e)
        call self%boundaries%init (radius=planet%radius)
        do iv=self%i_px,self%i_pz
          call self%boundaries%spherical (self%mem(:,:,:,iv,it,1), 0.0, self%mesh, iv)
        end do
        deallocate (mem)
      end if
      !------------------------------------------------------------------------
      ! Edge boundary conditions
      !------------------------------------------------------------------------
      if (self%boundaries%status /= 0) then
        nb = 0
        if (self%boundaries%is_set (bits%xl+bits%xu)) then
          nb = nb+1
          call self%edge_boundary_conditions(1)
        end if
        if (self%boundaries%is_set (bits%yl+bits%yu)) then
          nb = nb+1
          call self%edge_boundary_conditions(2)
        end if
        if (self%boundaries%is_set (bits%zl+bits%zu)) then
          nb = nb+1
          call self%edge_boundary_conditions(3)
        end if
    !    if (nb == 1) print *, self%id, 'BOUNDARY CONDITION: face'
    !    if (nb == 2) print *, self%id, 'BOUNDARY CONDITION: edge'
    !    if (nb == 3) print *, self%id, 'BOUNDARY CONDITION: corner'
      end if
    end if
  end if
  !call self%courant_condition
   if (self%id==42) then
   associate (r1=>self%mesh(1)%r, r2=>self%mesh(2)%r, r3=>self%mesh(3)%r)
 
    print *, '=== VERTICAL forces in position ===', self%mesh%p
    print *, 'time, dtime', self%time, self%dtime
    print *, 'iz, pos, p_x, p_y, p_z, ee, d'
    do i = self%mesh(3)%lb, self%mesh(3)%ub
      print *, i, r3(i), self%mem(14,14,i,self%i_px,self%it,1), self%mem(14,14,i,self%i_py,self%it,1), &
      self%mem(14,14,i,self%i_pz,self%it,1),self%mem(14,14,i,self%i_e,self%it,1),self%mem(14,14,i,self%i_d,self%it,1)
    end do
  end associate
  end if
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

!> =============================================================================
!> Boundary conditions at the edges of the outermost rubiks cubes; extrapolate
!> the density and pressure, and enforce vanishing derivative on the perpendicular
!> velocities. To what extend the boundary conditions are applied in the x,y,z
!> directions depends on the self%boundary%status bits
!> =============================================================================
SUBROUTINE edge_boundary_conditions (self, ii)
  class(experiment_t):: self
  integer, optional:: ii
  real, dimension(:,:,:), pointer:: d, e, p, u
  integer:: iv
  !-----------------------------------------------------------------------------
  ! Edge boundary conditions
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%edge_boundary_conditions')
  d => self%mem(:,:,:,self%i_d,self%it,1)
  e => self%mem(:,:,:,self%i_e,self%it,1)
  allocate (u(size(d,1),size(d,2),size(d,3)))
  call self%boundaries%extrapolate (d)
  u = e/d
  call self%boundaries%symmetric (u)
  e = u*d
  do iv=self%i_px,self%i_pz
    p => self%mem(:,:,:,iv,self%it,1)
    u = p/d
    call self%boundaries%symmetric (u)
    p = u*d
  end do
  deallocate (u)
  call trace%end
END SUBROUTINE edge_boundary_conditions

SUBROUTINE constant_boundary_conditions (self, wl)
  class(experiment_t):: self
  integer, optional :: wl
  class(mesh_t), pointer, dimension(:) :: m
  integer:: iv, iz, iy, ix, lz, uz, ly, uy, lx, ux
  real(8) :: x, y, z, r, p, d, px, py, pz
  !-----------------------------------------------------------------------------
  ! Edge boundary conditions
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%constant_boundary_conditions')

  m => self%mesh
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
    
  select case(wl)
  case(1)
    if(self%boundaries%is_set(bits%xl)) then
      lx = m(1)%lb; ux = m(1)%lo
      ly = m(2)%lb; uy = m(2)%ub
      lz = m(3)%lb; uz = m(3)%ub
    
    else if (self%boundaries%is_set(bits%xu)) then
      lx = m(1)%uo; ux = m(1)%ub
      ly = m(2)%lb; uy = m(2)%ub
      lz = m(3)%lb; uz = m(3)%ub
    else
      error stop "cbc (1)"
    end if
  case(2)
    if(self%boundaries%is_set(bits%yl)) then
      lx = m(1)%lb; ux = m(1)%ub
      ly = m(2)%lb; uy = m(2)%lo
      lz = m(3)%lb; uz = m(3)%ub
    
    else if (self%boundaries%is_set(bits%yu)) then
      lx = m(1)%lb; ux = m(1)%ub
      ly = m(2)%uo; uy = m(2)%ub
      lz = m(3)%lb; uz = m(3)%ub
    else
      error stop "cbc (2)"
    end if
  case(3)
    if(self%boundaries%is_set(bits%zl)) then
      lx = m(1)%lb; ux = m(1)%ub
      ly = m(2)%lb; uy = m(2)%ub
      lz = m(3)%lb; uz = m(3)%lo
    
    else if (self%boundaries%is_set(bits%zu)) then
      lx = m(1)%lb; ux = m(1)%ub
      ly = m(2)%lb; uy = m(2)%ub
      lz = m(3)%uo; uz = m(3)%ub
    else
      error stop "cbc (3)"
    end if
  end select
  
  
  do iz=lz,uz
    z = m(3)%p + r3(iz) - planet%position(3)
    do iy=ly,uy
      y = m(2)%p + r2(iy) - planet%position(2)
      do ix=lx,ux
        x = m(1)%p + r1(ix) - planet%position(1)
        r = sqrt(x**2+y**2+z**2)
        d = atmosphere%external%density
        p = atmosphere%external%pressure
        self%mem(ix,iy,iz,self%i_d,self%it,1) = d
        self%mem(ix,iy,iz,self%i_e,self%it,1) = p/self%gamma
        !-----------------------------------------------------------------------
        ! A rotation profile, where the rotation center is at y=+1 AU, and
        ! the direction of rotation is +x, so the head wind is in the -x
        ! direction.  The slower rotation at larger radii means that Ux is
        ! proportional to +y.  This gives a speed difference of only 130 cm/s
        ! across the Earth, and plus/minus 0.15 km/s across two Hill radii.
        ! The effect is not enough to significantly flatten the outer part of
        ! the Hill sphere, so we will probably only be seeing effects of the
        ! Coriolis force beyond the Hill radius.
        !-----------------------------------------------------------------------
        px = -self%initial%v_accrete
        py = -1.5*planet%omega(3)*x - self%initial%v_head
        pz = 0.0
        
        self%mem(ix,iy,iz,self%i_px, self%it, 1) = d*px ! radial speed
        self%mem(ix,iy,iz,self%i_py, self%it, 1) = d*py ! azimuthal speed
        self%mem(ix,iy,iz,self%i_pz, self%it, 1) = d*pz ! vertical speed

      end do
    end do
  end do
  end associate
  
  call trace%end
END SUBROUTINE constant_boundary_conditions

END MODULE
