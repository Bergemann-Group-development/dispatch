!> =============================================================================
!> Setup the initial conditions and experiment update
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE pebble_patch_mod
  USE initial_mod
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
    procedure :: init
    procedure :: update
    procedure :: output
    procedure :: boundary_conditions
    procedure :: edge_boundary_conditions
    procedure :: constant_boundary_conditions
    procedure :: density_boundary_condition
  end type
  real, save:: temp0=0.0
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
  integer:: iv, iostat
  logical, save:: first_time=.true.
  character(len=32):: type='planet_atmosphere'
  namelist /experiment_params/ type
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')
  !$omp critical (io_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input,experiment_params) !, iostat=iostat)
    if (io%master) write (*,experiment_params)
  end if
  !$omp end critical (io_cr)
  !----------------------------------------------------------------------------
  ! Allocate memory and set IC values.  The initial%init (type) specifies the
  ! type of IC to be applied when the call to gpatch_t%init is made.
  !----------------------------------------------------------------------------
  self%nv = 5
  self%msplit = 3
  call scaling%init
  call self%initial%init (type)
  call self%force%init (type)
  self%nv = 5
  self%mhd = .false.
  !----------------------------------------------------------------------------
  ! Through the hierarchy of class(task_t)%init calls, the call to %init below
  ! sets the IC values, either from reading in snapshot values, or by calling
  ! the IC condition procedure self%initial%selected() from gpatch_t%init.
  !----------------------------------------------------------------------------
  call self%pebble_patch_t%init
  self%mhd = .false.
  self%nv = 5
  call self%set (bits%static)
  call self%boundaries%init (self%mesh, position=planet%position, radius=planet%radius)
  if (trim(type)/='disk' .and. all(abs(self%position) < self%ds)) then
    print *, 'experiment_t%init: setting spherical BC bits on patch', self%id
    call self%set (bits%has_BC)
    call self%boundaries%set (bits%spherical)
    do iv=self%i_px,self%i_pz
      call self%boundaries%spherical (self%mem(:,:,:,iv,self%it,1), 0.0, self%mesh, iv)
    end do
  end if
  if (self%restart<0) &
    call self%output
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Choice of update steps for this experiment
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !.............................................................................
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
  if (self%only_pebbles .or. self%is_set(bits%not_leaf)) then
    if (self%istep<self%nt+1) then
      call self%mhd_t%update
    end if
    if (self%is_set (bits%spherical)) then
      do iv=self%i_px,self%i_pz
        call self%boundaries%spherical (self%mem(:,:,:,iv,self%it,2), 0.0, self%mesh, iv)
      end do
    end if
  !-----------------------------------------------------------------------------
  ! Normal time stepping, with boundary conditions applied after the update
  !-----------------------------------------------------------------------------
  else
    !--------------------------------------------------------------------------
    ! Update the MHD state
    !--------------------------------------------------------------------------
    call self%boundary_conditions (self%it)
    call self%mhd_t%update
  end if
  !-----------------------------------------------------------------------------
  if (self%id==io%id_debug) then
   associate (r1=>self%mesh(1)%r, r2=>self%mesh(2)%r, r3=>self%mesh(3)%r)
    print *, '=== VERTICAL forces in position ===', self%mesh%p
    print *, 'time, dtime', self%time, self%dtime
    print *, 'iz, pos, p_x, p_y, p_z, ee, d'
    do i = self%mesh(3)%lb, self%mesh(3)%ub
      print *, i, r3(i), self%mem(14,14,i,self%i_px,self%it,1), self%mem(14,14,i,self%i_py,self%it,1), &
      self%mem(14,14,i,self%i_pz,self%it,1),self%mem(14,14,i,self%i_e,self%it,1),self%mem(14,14,i,self%i_d,self%it,1)
    end do
    print *, 'min T', minval((self%mem(:,:,:,self%i_e,self%it,1)/self%mem(:,:,:,self%i_d,self%it,1)*(self%gamma-1)*scaling%temp))
    print *, 'max T', maxval((self%mem(:,:,:,self%i_e,self%it,1)/self%mem(:,:,:,self%i_d,self%it,1)*(self%gamma-1)*scaling%temp))
    print *, 'min max d', minval(self%mem(:,:,:,self%i_d,self%it,1)), maxval (self%mem(:,:,:,self%i_d,self%it,1))
    print *, 'min max e', minval(self%mem(:,:,:,self%i_e,self%it,1)), maxval (self%mem(:,:,:,self%i_e,self%it,1))
   end associate
  end if
  call self%pebble_patch_t%update_pebbles
  call self%check_density ('after timestep')
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Interface to patch_t%output, for debug printout
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(experiment_t):: self
  character(len=64), optional:: experiment_name
  integer:: iv, it, m(3)
  if (self%id==io%id_debug) then
    m = self%mesh%gn/2
    it = self%iit(self%nt-2)
    do iv=2,4
      print *,'OUT: t, iv, u(z) =', real(self%t(it)), iv, &
        self%mem(m(1),m(2),m(3)-2:m(3)+2,iv,it,1)/ &
        self%mem(m(1),m(2),m(3)-2:m(3)+2, 1,it,1)
    end do
  end if
  call self%boundary_conditions (self%it)
  call self%patch_t%output
END SUBROUTINE

!> =============================================================================
!> Boundary conditions for this experiment
!> =============================================================================
SUBROUTINE boundary_conditions (self, it)
  class(experiment_t):: self
  integer:: it
  !.............................................................................
  real(8), dimension(3):: h
  integer:: iv, nb,i
  !-----------------------------------------------------------------------------
  if (self%boundaries%status /= 0) then
  !  print *, self%id, 'BOUNDARY CONDITION BITS:', self%boundaries%status
    !---------------------------------------------------------------------------
    ! If spherical inner BCs, save the incoming values
    !---------------------------------------------------------------------------
    if (self%boundaries%is_set (bits%spherical)) then
      if (io%verbose>1) print *, self%id, 'BOUNDARY CONDITION: spherical'
      print *, self%id, 'BOUNDARY CONDITION: spherical, t =', self%t(self%it)
      call self%boundaries%spherical_extrapolate (self%mem(:,:,:,self%i_d,it,1), self%mesh, self%i_d)
      call self%boundaries%spherical_extrapolate (self%mem(:,:,:,self%i_e,it,1), self%mesh, self%i_e)
      do iv=self%i_px,self%i_pz
        call self%boundaries%spherical (self%mem(:,:,:,iv,it,1), 0.0, self%mesh, iv)
        if (self%nw>1) &
          call self%boundaries%spherical (self%mem(:,:,:,iv,it,2), 0.0, self%mesh, iv)
      end do
    end if
    !------------------------------------------------------------------------
    ! Edge boundary conditions
    !------------------------------------------------------------------------
    if (self%boundaries%status /= 0) then
      nb = 0
      if (self%boundaries%is_set (bits%xl+bits%xu)) then
        nb = nb+1
        call self%edge_boundary_conditions(it, 1)
      end if
      if (self%boundaries%is_set (bits%yl+bits%yu)) then
        nb = nb+1
        call self%edge_boundary_conditions(it, 2)
      end if
      if (self%boundaries%is_set (bits%zl+bits%zu)) then
        nb = nb+1
        call self%edge_boundary_conditions(it, 3)
      end if
  !    if (nb == 1) print *, self%id, 'BOUNDARY CONDITION: face'
  !    if (nb == 2) print *, self%id, 'BOUNDARY CONDITION: edge'
  !    if (nb == 3) print *, self%id, 'BOUNDARY CONDITION: corner'
    end if
  end if
END SUBROUTINE boundary_conditions

!> =============================================================================
!> Boundary conditions at the edges of the outermost rubiks cubes; extrapolate
!> the density and pressure, and enforce vanishing derivative on the perpendicular
!> velocities. To what extend the boundary conditions are applied in the x,y,z
!> directions depends on the self%boundary%status bits
!> =============================================================================
SUBROUTINE edge_boundary_conditions (self, it, ii)
  class(experiment_t):: self
  integer:: it
  integer, optional:: ii
  real, dimension(:,:,:), pointer:: d, e, p, u, k
  integer:: iv
  real:: ee
  !-----------------------------------------------------------------------------
  ! Edge boundary conditions
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%edge_boundary_conditions')
  self%boundaries%id = self%id
  d => self%mem(:,:,:,self%i_d,it,1)
  e => self%mem(:,:,:,self%i_e,it,1)
  allocate (u(size(d,1),size(d,2),size(d,3)))
  allocate (k(size(d,1),size(d,2),size(d,3)))
  call self%density_boundary_condition (it,ii)  ! vertical density BC
  k = 0.0                                       ! initialize kinetic energy
  do iv=self%i_px,self%i_pz
    p => self%mem(:,:,:,iv,it,1)                ! momentum
    u = p/d                                     ! velocity
    call self%boundaries%symmetric (u)          ! zero gradient velocity
    p = u*d                                     ! recontruct momentum
    k = k+p**2                                  ! accumulate kinetic energy
  end do
  k = 0.5*k/d                                   ! normalize kinetic energy
  u = (e-k)/d                                   ! internal energy
  ee = atmosphere%external%temperature &        ! external disk e/d
    /(self%gamma-1.0)
  call self%boundaries%constant (u, ee)         ! constant bdry value
  e = u*d+k                                     ! reconstruct total energy
  deallocate (u,k)
  call trace%end
END SUBROUTINE edge_boundary_conditions

!> =====================================================================
!> fixed boundary conditions to test that normal BCs work
!> =====================================================================
SUBROUTINE constant_boundary_conditions (self, it, wl)
  class(experiment_t):: self
  integer:: it
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
        self%mem(ix,iy,iz,self%i_d,it,1) = d
        self%mem(ix,iy,iz,self%i_e,it,1) = p/self%gamma
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
        px = -planet%v_accrete
        py = -1.5*planet%omega(3)*x - planet%v_head
        pz = 0.0
        self%mem(ix,iy,iz,self%i_px,it,1) = d*px ! radial speed
        self%mem(ix,iy,iz,self%i_py,it,1) = d*py ! azimuthal speed
        self%mem(ix,iy,iz,self%i_pz,it,1) = d*pz ! vertical speed
      end do
    end do
  end do
  end associate
  call trace%end
END SUBROUTINE constant_boundary_conditions

!> =====================================================================
!> Constrain the density to follow an achored Gaussian profile
!> =====================================================================
SUBROUTINE density_boundary_condition (self, it, ii)
  class(experiment_t):: self
  integer:: it, ii
  class(mesh_t), pointer, dimension(:) :: m
  integer:: iv, iz, iy, ix, lz, uz, ly, uy, lx, ux, ix0, iz0
  real(8) :: x, x0, z, z0, r, d1, d0, tt, c2
  real, dimension(:,:,:), pointer:: d
  !-----------------------------------------------------------------------------
  ! Vertical edge boundary conditions
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%density_boundary_conditions')
  m => self%mesh
  d => self%mem(:,:,:,1,it,1)
  associate (r1=>m(1)%r, r3=>m(3)%r, mem=>self%mem)
  if (self%id==io%id_debug) &
    print '(a,2i5,2l4,3(2x,2i3))', 'MK', self%id, it, &
      self%boundaries%is_set(bits%zl), &
      self%boundaries%is_set(bits%zu), lx,ux, ly,uy, lz,uz
  tt = atmosphere%external%temperature
  c2 = 0.5*planet%omega(3)**2/tt
  if (ii==1) then
    if(self%boundaries%is_set(bits%xl)) then
      lx = m(1)%lb; ux = m(1)%lo
      ly = m(2)%lb; uy = m(2)%ub
      lz = m(3)%lb; uz = m(3)%ub
      ix0 = ux+1
      x0 = m(1)%p + r1(ix0) - planet%position(1)
    else if (self%boundaries%is_set(bits%xu)) then
      lx = m(1)%uo; ux = m(1)%ub
      ly = m(2)%lb; uy = m(2)%ub
      lz = m(3)%lb; uz = m(3)%ub
      ix0 = lx-1
      x0 = m(1)%p + r1(ix0) - planet%position(1)
    else
      error stop "dbc"
    end if
    c2 = c2*2.0*planet%trap
    do iz=lz,uz
      do iy=ly,uy
        do ix=lx,ux
          x = m(1)%p + r1(ix) - planet%position(1)
          d0 = mem(ix0,iy,iz,1,it,1)
          self%mem(ix,iy,iz,1,it,1) = d0*exp(-c2*(x**2-x0**2))
        end do
      end do
    end do
  else if (ii==2) then
    call self%boundaries%symmetric (d)
  else if (ii==3) then
    if(self%boundaries%is_set(bits%zl)) then
      lx = m(1)%lb; ux = m(1)%ub
      ly = m(2)%lb; uy = m(2)%ub
      lz = m(3)%lb; uz = m(3)%lo
      iz0 = uz+1
      z0 = m(3)%p + r3(iz0) - planet%position(3)
    else if (self%boundaries%is_set(bits%zu)) then
      lx = m(1)%lb; ux = m(1)%ub
      ly = m(2)%lb; uy = m(2)%ub
      lz = m(3)%uo; uz = m(3)%ub
      iz0 = lz-1
      z0 = m(3)%p + r3(iz0) - planet%position(3)
    else
      error stop "dbc"
    end if
    do iz=lz,uz
      z = m(3)%p + r3(iz) - planet%position(3)
      do iy=ly,uy
        do ix=lx,ux
          d0 = mem(ix,iy,iz0,1,it,1)
          d1 = d0*exp(-c2*(z**2-z0**2))
          if (self%id==io%id_debug.and.iy==uy/2.and.iz==uz) &
            print '(a,i4,1p,8e14.6)', &
            'BCd: ix,d0,p0/do,exp,do,d =',ix, d0, tt, (z-z0), &
            exp(-c2*(z**2-z0**2)), &
            self%mem(ix,iy,iz,1,it,1), d1
          self%mem(ix,iy,iz,1,it,1) = d1
        end do
      end do
    end do
  end if
  end associate
  call trace%end
END SUBROUTINE density_boundary_condition

END MODULE
