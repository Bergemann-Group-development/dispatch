!===============================================================================
!> Setup the initial conditions and experiment update
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE os_mod
  USE omp_mod
  USE bits_mod
  USE trace_mod
  USE kinds_mod
  USE pebble_patch_mod
  USE initial_mod
  USE atmosphere_mod
  USE units_mod
  USE scaling_mod
  USE planet_mod
  USE mesh_mod
  USE boundaries_mod
  USE scalar_mod
  USE eos_mod
  USE lagrange_mod
  USE timestep_mod
  !USE rt_boundaries_mod
  USE eos_mod
  USE force_mod
  USE scratch_mod
  implicit none
  private
  type, public, extends(pebble_patch_t):: experiment_t
    type(force_t):: force
    type(initial_t):: initial
    real:: mdot           ! accretion rate in code units
    real(8):: x_sponge    ! distance from the physical boundary where
                          ! "wave-killing" boundaries are turned on.
  contains
    procedure :: init
    procedure :: update
    procedure :: rotate
    procedure :: output
    procedure :: boundary_conditions
    procedure :: edge_boundary_conditions
    procedure :: symmetric_z_u
    procedure :: density_boundary_condition
    procedure :: is_relevant
  end type
  integer, save:: verbose=0
  integer, save:: track_level=99
  logical, save:: unsigned=.true.
  logical, save:: pervolume=.false.
  logical, save:: first_time=.true.
  character(len=32), save:: type='planet_atmosphere'
CONTAINS

!===============================================================================
!> Initializations
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  ! ............................................................................
  integer:: iv, iostat
  logical, save:: first_time=.true.
  real(8):: x_sponge=-1.0e30
  real, pointer:: mem(:,:,:,:)
  namelist /experiment_params/ type, x_sponge, verbose, track_level, unsigned, &
    pervolume
  ! namelist /experiment_params/ type, x_sponge, me_per_year, spherical_bc_corners, &
  !   use_point_source, use_attenuation, epsilon, diffusion_radius, diffuse_radius, &
  !   disk_heat_cgs, track_level, verbose
  !-----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace%begin('experiment_t%init')
  call io%lock%set ()
  if (first_time) then
    rewind (io%input)
    read (io%input,experiment_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input,experiment_params)
    end if
    write (io%namelist,experiment_params)
  end if
  call io%lock%unset ()
  self%x_sponge = x_sponge
  if (self%level >= track_level) self%track = .true.
  !-----------------------------------------------------------------------------
  ! Allocate memory and set IC values.  The initial%init (type) specifies the
  ! type of IC to be applied when the call to gpatch_t%init is made.
  !-----------------------------------------------------------------------------
  self%msplit = 3
  call scaling%init
  self%nv = 5
  self%periodic = .false.
  self%periodic(2) = .true.
  self%initial%initial_condition = type
  !-----------------------------------------------------------------------------
  ! Through the hierarchy of class(task_t)%init calls, the call to %init below
  ! sets the IC values, either from reading in snapshot values, or by calling
  ! the IC condition procedure self%initial%condition() from gpatch_t%init.
  !-----------------------------------------------------------------------------
  call self%solver_t%init
  self%idx%e = max(self%idx%s, self%idx%e)
  self%idx%s = max(self%idx%s, self%idx%e)
  call self%idx%copy
  if (first_time) then
    call self%idx%print (stdout)
    first_time = .false.
  end if
  !-----------------------------------------------------------------------------
  ! Initial conditions
  !-----------------------------------------------------------------------------
  call self%initial%init (self%kind, self%etype)
  self%initial%id = self%id
  call self%force%init (self%staggered, self%id, self%mesh)
  mem => self%mem(:,:,:,:,self%it,1)
  call self%initial%condition (self%mesh, mem, self%idx)
  !-----------------------------------------------------------------------------
  ! Pebbles
  !-----------------------------------------------------------------------------
  call self%pebble_patch_t%init
  call self%set (bits%static)
  !-----------------------------------------------------------------------------
  ! Boundary conditions
  !-----------------------------------------------------------------------------
  call self%boundaries%init ()
  call self%boundaries%set_attribute (radius = planet%radius)
  if (self%boundaries%is_set (bits%spherical)) then
    !---------------------------------------------------------------------------
    ! Only the innermost central patch needs spherica BCs
    !---------------------------------------------------------------------------
    if (all(self%size > planet%radius*4.0)) then
      call self%boundaries%clear (bits%spherical)
      if (verbose > 0) &
        write (stdout,'(3x,a,i6)') &
          'experiment_t%init: clearing spherical BC bits on patch', self%id
    else
      if (verbose > 0) &
        write (stdout,'(3x,a,i6)') &
          'experiment_t%init: setting spherical BC bits on patch', self%id
      call self%set (bits%has_BC)
      do iv=self%idx%px,self%idx%pz
        call self%boundaries%spherical (mem(:,:,:,iv), 0.0, self%mesh, iv)
      end do
    end if
  end if
  call self%boundary_conditions (self%it, self%it)
  !-----------------------------------------------------------------------------
  ! Choose interpolation types. The smallest initial glitch (positive) is with
  !
  ! self%unsigned(self%idx%d) = .true.
  ! self%pervolume(self%idx%s) = .false.
  !
  ! All other combinations give a stronger negative initial glitch in p1 at
  ! mesh%li. Both types of glitches come from the different precision in the
  ! initial hydrostatic state, and go away with time. The glitch is visible
  ! in u1 at index li when the pressure on one side comes from lower resolution.
  !
  ! Note that the current default in stagger2/mhd_mod.f90 is
  !
  ! self%unsigned(self%idx%d) = .true.
  ! self%pervolume(self%idx%s) = .true.
  !
  ! Tests show that these should either both be true, or both be false.
  !-----------------------------------------------------------------------------
  self%unsigned(self%idx%d) = unsigned
  self%pervolume(self%idx%px) = pervolume
  self%pervolume(self%idx%py) = pervolume
  self%pervolume(self%idx%pz) = pervolume
  self%pervolume(self%idx%s) = unsigned
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Choice of update steps for this experiment
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: d, e
  real(var_k), dimension(:,:,:,:), pointer :: p
  real(var_k), dimension(:,:,:),   pointer, contiguous :: Ux, Uy, Uz, pg
  real(var_k), dimension(:,:,:,:), pointer, contiguous :: U
  real(8), dimension(3):: h
  integer:: it, iv, nb, i
  logical:: panic
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%update', itimer=itimer)
  !.............................................................................
  call self%stats ('start update', it=self%it, verbose=verbose)
  call self%output
  !-----------------------------------------------------------------------------
  ! Add forcing -- this could be moved to extras_mod
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  call os%allocate (U, [self%gn,3])
  call self%p2U (U, self%it, linear=.true.)
  Ux => U(:,:,:,1)
  Uy => U(:,:,:,2)
  Uz => U(:,:,:,3)
  !-----------------------------------------------------------------------------
  ! Thread-private scratch%force needs to be allocated here; %init() runs single
  !-----------------------------------------------------------------------------
  call scratch%allocate (scratch%force_per_unit_mass, [self%gn,3])
  scratch%force_per_unit_mass = &
    self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)
  !-----------------------------------------------------------------------------
  ! Update the HD patch, extrapolate guard cells (needed in BC)
  !-----------------------------------------------------------------------------
  self%pre_updated = .true.   ! prevent force resetting!
  call self%solver_t%update
  call self%stats ('after rmht_t%update', it=self%new, verbose=verbose)
  call self%extrapolate_guards ()
  call self%boundary_conditions (self%it, self%new)
  call self%stats ('after BCs', it=self%new, verbose=verbose)
  call os%deallocate (U)
  !-----------------------------------------------------------------------------
  ! Update the pebbles, if used
  !-----------------------------------------------------------------------------
  if (self%do_pebbles) then
    pg => self%pressure(:,:,:,self%it)
  end if
  ! call self%boundary_conditions (self%it, self%new)
  call self%stats ('after BCs(2)', it=self%new, verbose=verbose)
  if (self%do_pebbles) then
    call self%pebble_patch_t%update_pebbles
  end if
  !-----------------------------------------------------------------------------
  call self%check_nan (panic, label='after timestep')
  call self%counter_update
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Make sure the RT task rotates slots correctly, in perfect sync with this task
!===============================================================================
SUBROUTINE rotate (self)
  class(experiment_t) :: self
  !-----------------------------------------------------------------------------
  ! Check (already here -- to cover for rt%rotate) that we don't do it twice,
  ! and that we actually have the same information in the rt%task_t data type
  !-----------------------------------------------------------------------------
  if (self%rotated) return
  call self%pebble_patch_t%rotate()
END SUBROUTINE rotate

!===============================================================================
!> Interface to patch_t%output, for debug printout
!===============================================================================
SUBROUTINE output (self)
  class(experiment_t):: self
  integer:: iv, it, m(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%output', itimer=itimer)
  if (self%id==io%id_debug) then
    m = self%mesh%gn/2
    it = self%iit(self%nt-2)
    do iv=2,4
      print *,'OUT: t, iv, u(z) =', real(self%t(it)), iv, &
        self%mem(m(1),m(2),m(3)-2:m(3)+2,iv,it,1)/ &
        self%mem(m(1),m(2),m(3)-2:m(3)+2, 1,it,1)
    end do
  end if
  !if (self%time > self%out_next) &
  !  call self%boundary_conditions (self%it)
  call self%pebble_patch_t%output
  call trace%end (itimer)
END SUBROUTINE

!===============================================================================
!> Boundary conditions for this experiment
!===============================================================================
SUBROUTINE boundary_conditions (self, it, new)
  class(experiment_t) :: self
  integer             :: it, new
  !.............................................................................
  real(8)             :: h(3)
  integer             :: iv, nb, i
  real, pointer       :: d(:,:,:), e(:,:,:), p(:,:,:,:), mem(:,:,:,:,:)
  real, pointer       :: ss(:,:,:), ss2(:,:,:)
  contiguous          :: ss, ss2
  real                :: s0
  integer, save       :: itimer=0
  !-----------------------------------------------------------------------------
  if (self%boundaries%bits /= 0) then
    call trace%begin ('experiment_t%boundary_conditions', itimer=itimer)
    call self%check_density ('before boundary_conditions')
    if (self%boundaries%is_set (bits%spherical)) then
      if (verbose > 0) print *, self%id, 'BOUNDARY CONDITION: spherical, t =', &
         self%t(self%it)
      !-------------------------------------------------------------------------
      ! Spherical inner BCs: vanishing mass flux & constant density and thermal
      ! energy / entropy
      !-------------------------------------------------------------------------
      mem => self%mem(:,:,:,:,:,1)
      p => mem(:,:,:,self%idx%px:self%idx%pz,new)
      call self%boundaries%spherical_p_ns (p, &
        self%mesh, planet%position, planet%radius)
      call self%boundaries%spherical_keep (mem, self%idx%d, it, new, &
        self%mesh, planet%position, planet%radius)
      call self%boundaries%spherical_keep (mem, self%idx%e, it, new, &
        self%mesh, planet%position, planet%radius)
    !---------------------------------------------------------------------------
    ! Edge boundary conditions
    !---------------------------------------------------------------------------
    else if (self%boundaries%bits /= 0) then
      if (verbose > 2) then
        print *, self%id, 'MK x-BC:', self%boundaries%is_set (bits%xl), self%boundaries%is_set (bits%xu)
        print *, self%id, 'MK y-BC:', self%boundaries%is_set (bits%yl), self%boundaries%is_set (bits%yu)
        print *, self%id, 'MK z-BC:', self%boundaries%is_set (bits%zl), self%boundaries%is_set (bits%zu)
      end if
      call self%edge_boundary_conditions(new)
    end if
    call trace%end (itimer)
  end if
END SUBROUTINE boundary_conditions

!===============================================================================
!> Boundary conditions at the edges of the outermost rubiks cubes; extrapolate
!> the density and pressure, and enforce vanishing derivative on the perpendicular
!> velocities. To what extend the boundary conditions are applied in the x,y,z
!> directions depends on the self%boundaries%bits bits. Note that the same
!> condition is applied in all directions that have boundary bits set!
!===============================================================================
SUBROUTINE edge_boundary_conditions (self, it)
  class(experiment_t):: self
  integer:: it
  real(var_k), dimension(:,:,:), pointer, contiguous:: d, s, ss
  real(var_k), dimension(:,:,:), pointer, contiguous:: px, py, pz
  real(var_k), dimension(:,:,:,:), pointer:: u
  integer:: iv, select, xbits, zbits
  logical:: radial, vertical
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%edge_boundary_conditions', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Allocate and load values into velocity, consistent with current density
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,it,1)
  !-----------------------------------------------------------------------------
  ! Density boundary conditions
  !-----------------------------------------------------------------------------
  call self%stats ('before density BC', it=it, verbose=verbose)
  call self%density_boundary_condition (it)
  !call self%boundaries%log_extrapolate (d, self%idx%d)
  call self%stats ('after density BC', it=it, verbose=verbose)
  !-----------------------------------------------------------------------------
  ! Velocity boundary conditions
  !-----------------------------------------------------------------------------
  px => self%mem(:,:,:,self%idx%px,it,1)
  py => self%mem(:,:,:,self%idx%py,it,1)
  pz => self%mem(:,:,:,self%idx%pz,it,1)
  !-----------------------------------------------------------------------------
  ! The azimuthal velocity should be safe to extrapolate in directions
  ! where it it parallel to the boundary, while the horizontal one should
  ! vanish and the vertical one should have vanishing derivative
  !-----------------------------------------------------------------------------
  radial   = self%boundaries%is_set (bits%xl) .or. &
             self%boundaries%is_set (bits%xu)
  vertical = self%boundaries%is_set (bits%zl) .or. &
             self%boundaries%is_set (bits%zu)
  !-----------------------------------------------------------------------------
  ! Radial BCs: zero radial, extrapolate azim, zero-gradient vertical
  !-----------------------------------------------------------------------------
  xbits = bits%xl + bits%xu
  if (radial) then
    call self%boundaries%antisymmetric (px, self%idx%px, select=xbits)
    call self%boundaries%extrapolate   (py, self%idx%py, select=xbits)
    call self%boundaries%symmetric     (pz, self%idx%pz, select=xbits)
  end if
  !-----------------------------------------------------------------------------
  ! Vertical BCs: zero-gradient of velocity components
  !-----------------------------------------------------------------------------
  if (vertical) then
    call self%symmetric_z_u   (px, self%idx%px)
    call self%symmetric_z_u   (py, self%idx%py)
    call self%symmetric_z_u   (pz, self%idx%pz)
  end if
  call self%stats ('after veocity BC', it=it, verbose=verbose)
  !-----------------------------------------------------------------------------
  ! Energy BCs: zero-gradient on energy or entropy per unit mass
  !-----------------------------------------------------------------------------
  allocate (ss(size(d,1),size(d,2),size(d,3)))
  s => self%mem(:,:,:,self%idx%s,it,1)
  ss = s/d
  zbits = bits%zl + bits%zu
  call self%boundaries%symmetric (ss, self%idx%s, select=xbits+zbits)
  s = ss*d
  deallocate (ss)
  !-----------------------------------------------------------------------------
  ! Apply the "wave killing" boundaries to all variables near the x-boundaries;
  !-----------------------------------------------------------------------------
  if (self%x_sponge > -1.0e-32) then
    call apply_sponge_boundaries (self)
  end if
  call trace%end (itimer)
END SUBROUTINE edge_boundary_conditions

!===============================================================================
!> Special case procedures to enforce zero-gradient on velocity components
!===============================================================================
SUBROUTINE symmetric_z_u (self, v, iv)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:), pointer, contiguous:: v
  integer:: iv
  !.............................................................................
  real(var_k), dimension(:,:,:), pointer, contiguous:: d
  integer:: ix, iy, iz, jz, iz1, iz2, jz1, jz2
  class(mesh_t), pointer:: m
  !-----------------------------------------------------------------------------
  m => self%mesh(3)
  d  => self%mem(:,:,:,self%idx%d,self%it,1)
  !-----------------------------------------------------------------------------
  ! Cell-centered case
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%zl)) then
    iz1 = 1   ; iz2 = m%lo; jz1 = iz2+1; jz2 = iz2
  else if (self%boundaries%is_set (bits%zu)) then
    iz1 = m%uo; iz2 = m%gn; jz1 = iz1-1; jz2 = iz1
  end if
  if (all(m%h(self%idx%px:self%idx%pz) == 0d0)) then
    do iz=iz1,iz2
      jz = jz1 + (jz2-iz)
      do iy=1,self%gn(2)
      do ix=1,self%gn(1)
        v(ix,iy,iz) = v(ix,iy,jz)*d(ix,iy,iz)/d(ix,iy,jz)
      end do
      end do
    end do
  !-----------------------------------------------------------------------------
  ! Staggered horizontally
  !-----------------------------------------------------------------------------
  else if (m%h(iv) == 0.0 .and. iv==self%idx%px) then
    do iz=iz1,iz2
      jz = jz1 + (jz2-iz)
      do iy=1,self%gn(2)
      do ix=1,self%gn(1)-1
        v(ix,iy,iz) = &
        v(ix,iy,jz)*(d(ix,iy,iz)+d(ix+1,iy,iz))/(d(ix,iy,jz)+d(ix+1,iy,jz))
      end do
      end do
    end do
  else if (m%h(iv) == 0.0 .and. iv==self%idx%py) then
    do iz=iz1,iz2
      jz = jz1 + (jz2-iz)
      do iy=1,self%gn(2)-1
      do ix=1,self%gn(1)
        v(ix,iy,iz) = &
        v(ix,iy,jz)*(d(ix,iy,iz)+d(ix,iy+1,iz))/(d(ix,iy,jz)+d(ix,iy+1,jz))
      end do
      end do
    end do
  !-----------------------------------------------------------------------------
  ! Staggered vertically, with one more guard zone index on the low side
  !-----------------------------------------------------------------------------
  else if (m%h(iv) < 0.0) then
    if (self%boundaries%is_set (bits%zl)) then
      iz1 = 1   ; iz2 = m%li; jz1 = iz2+1; jz2 = iz2
    else if (self%boundaries%is_set (bits%zu)) then
      iz1 = m%uo; iz2 = m%gn-1; jz1 = iz1-1; jz2 = iz1
    end if
    do iz=iz1,iz2
      jz = jz1 + (jz2-iz)
      do iy=1,self%gn(2)-1
      do ix=1,self%gn(1)-1
        v(ix,iy,iz) = &
        v(ix,iy,jz)*(d(ix,iy,iz)+d(ix,iy,iz+1))/(d(ix,iy,jz)+d(ix,iy,jz+1))
      end do
      end do
    end do
  end if
END SUBROUTINE symmetric_z_U

!===============================================================================
!> Constrain the density to follow an anchored Gaussian profile
!===============================================================================
SUBROUTINE density_boundary_condition (self, it)
  class(experiment_t):: self
  integer:: it
  class(mesh_t), pointer, dimension(:) :: m
  integer:: iv, iz, iy, ix, lz, uz, ly, uy, lx, ux, ix0, iy0, iz0
  real(8) :: x, x0, z, z0, r, d1, d0, tt, c2
  real(var_k), dimension(:,:,:), pointer:: d, dd
  !-----------------------------------------------------------------------------
  ! Vertical edge boundary conditions
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%density_boundary_conditions')
  m => self%mesh
  d => self%mem(:,:,:,1,it,1)
  if (self%nw>1) dd => self%mem(:,:,:,self%idx%d,it,2)
  associate (r1=>m(1)%r, r3=>m(3)%r)
  tt = atmosphere%external%temperature
  if(self%boundaries%is_set(bits%xl)) then
    lx = m(1)%lb; ux = m(1)%lo
    ly = m(2)%lb; uy = m(2)%ub
    lz = m(3)%lb; uz = m(3)%ub
    ix0 = m(1)%li
    x0 = m(1)%p + r1(ix0) - planet%position(1)
    call apply_x
  end if
  if (self%boundaries%is_set(bits%xu)) then
    lx = m(1)%uo; ux = m(1)%ub
    ly = m(2)%lb; uy = m(2)%ub
    lz = m(3)%lb; uz = m(3)%ub
    ix0 = m(1)%ui
    x0 = m(1)%p + r1(ix0) - planet%position(1)
    call apply_x
  end if
  if(self%boundaries%is_set(bits%yl)) then
    lx = m(1)%lb; ux = m(1)%ub
    ly = m(2)%lb; uy = m(2)%lo
    lz = m(3)%lb; uz = m(3)%ub
    iy0 = m(2)%li
    call apply_y
  end if
  if (self%boundaries%is_set(bits%yu)) then
    lx = m(1)%lb; ux = m(1)%ub
    ly = m(2)%uo; uy = m(2)%ub
    lz = m(3)%lb; uz = m(3)%ub
    iy0 = m(2)%ui
    call apply_y
  end if
  if(self%boundaries%is_set(bits%zl)) then
    lx = m(1)%lb; ux = m(1)%ub
    ly = m(2)%lb; uy = m(2)%ub
    lz = m(3)%lb; uz = m(3)%lo
    iz0 = m(3)%li
    z0 = m(3)%p + r3(iz0) - planet%position(3)
    call apply_z
  end if
  if (self%boundaries%is_set(bits%zu)) then
    lx = m(1)%lb; ux = m(1)%ub
    ly = m(2)%lb; uy = m(2)%ub
    lz = m(3)%uo; uz = m(3)%ub
    iz0 = m(3)%ui
    z0 = m(3)%p + r3(iz0) - planet%position(3)
    call apply_z
  end if
  !if (self%id==io%id_debug) &
  !  print '(a,2i5,2l4,3(2x,2i5))', 'MK', self%id, it, &
  !    self%boundaries%is_set(bits%zl), &
  !    self%boundaries%is_set(bits%zu), lx,ux, ly,uy, lz,uz
  call self%check_density ('after density BC')
  end associate
  call trace%end()
  return
contains
  subroutine apply_x
    if (self%boundaries%check_filled) call trace%begin ('apply_x')
    if (verbose > 2) &
      print *, self%id, 'MK: x-density BC', lx, ux, ly, uy, lz, uz
    c2 = planet%omega(3)**2*planet%trap/tt
    associate (r1=>m(1)%r)
    do iz=lz,uz
      do iy=ly,uy
        do ix=lx,ux
          x = m(1)%p + r1(ix) - planet%position(1)
          d0 = self%mem(ix0,iy,iz,self%idx%d,it,1)
          self%mem(ix,iy,iz,self%idx%d,it,1) = d0*exp(-c2*(x**2-x0**2))
        end do
        if (self%boundaries%check_filled) then
          do ix=lx,ux
            self%boundaries%filled(ix,iy,iz,self%idx%d) = .true.
          end do
        end if
      end do
    end do
    end associate
    if (self%nw>1) then
      do iz=lz,uz
        do iy=ly,uy
          do ix=lx,ux
            self%mem(ix,iy,iz,self%idx%d,it,2) = 0.0
          end do
        end do
      end do
    end if
    if (self%boundaries%check_filled) call trace%end()
  end subroutine apply_x
  subroutine apply_y
    if (self%boundaries%check_filled) call trace%begin ('apply_y')
    if (verbose > 2) &
      print *, self%id, 'MK: y-density BC', lx, ux, ly, uy, lz, uz
    do iz=lz,uz
      do iy=ly,uy
        do ix=lx,ux
          self%mem(ix,iy,iz,self%idx%d,it,1) = self%mem(ix,iy0,iz,self%idx%d,it,1)
        end do
        if (self%boundaries%check_filled) then
          do ix=lx,ux
            self%boundaries%filled(ix,iy,iz,self%idx%d) = .true.
          end do
        end if
      end do
    end do
    if (self%nw>1) then
      do iz=lz,uz
        do iy=ly,uy
          do ix=lx,ux
            self%mem(ix,iy,iz,self%idx%d,it,2) = 0.0
          end do
        end do
      end do
    end if
    if (self%boundaries%check_filled) call trace%end()
  end subroutine apply_y
  subroutine apply_z
    if (self%boundaries%check_filled) call trace%begin ('apply_z')
    if (verbose > 2) &
    print *, self%id, 'MK: z-density BC', lx, ux, ly, uy, lz, uz
    c2 = 0.5*planet%omega(3)**2/tt
    associate (r3=>m(3)%r)
    do iz=lz,uz
      z = m(3)%p + r3(iz) - planet%position(3)
      do iy=ly,uy
        do ix=lx,ux
          d0 = self%mem(ix,iy,iz0,self%idx%d,it,1)
          d1 = d0*exp(-c2*(z**2-z0**2))
          if (self%id==io%id_debug.and.iy==uy/2.and.iz==uz) &
            print '(a,i4,1p,8e14.6)', &
            'BCd: ix,d0,p0/do,exp,do,d =',ix, d0, tt, (z-z0), &
            exp(-c2*(z**2-z0**2)), &
            self%mem(ix,iy,iz,self%idx%d,it,1), d1
          self%mem(ix,iy,iz,self%idx%d,it,1) = d1
        end do
        if (self%boundaries%check_filled) then
          do ix=lx,ux
            self%boundaries%filled(ix,iy,iz,self%idx%d) = .true.
          end do
        end if
      end do
    end do
    end associate
    if (self%nw>1) then
      do iz=lz,uz
        do iy=ly,uy
          do ix=lx,ux
            self%mem(ix,iy,iz,self%idx%d,it,2) = 0.0
          end do
        end do
     end do
    end if
    if (self%boundaries%check_filled) call trace%end()
  end subroutine apply_z
END SUBROUTINE density_boundary_condition

!===============================================================================
!> Apply the "wave killing" boundaries of de Val-Borro et al. (2006) to all
!> variables.
!> Currently implemented for x-boundaries only.
!===============================================================================
SUBROUTINE apply_sponge_boundaries (self)
  USE units_mod, only: cgs
  class(experiment_t):: self
  integer:: iv, i, j, k, eye
  real(var_k):: x, z, xkill_lo, xkill_hi, xbndy_lo, xbndy_hi, tau_kill
  real(var_k):: d0, d, p0, p, T0, entropy0, px0, py0, pz0, &
                             asymptote
  class(mesh_t), pointer:: m1, m2, m3
  real(var_k), dimension(:,:,:,:), pointer:: mem
  !-----------------------------------------------------------------------------
  if (.not. (self%boundaries%is_set(bits%xl) .or. &
      self%boundaries%is_set(bits%xu))) &
    return
  call trace%begin ("experiment_t::apply_sponge_boundaries")

  if (io%debug(2) .or. self%track .and. verbose > 0) &
    print*, self%id, ": Wave killing ON"
  mem => self%mem(:,:,:,:,self%it,1)
  m1 => self%mesh(1)
  m2 => self%mesh(2)
  m3 => self%mesh(3)
  tau_kill = 2.0 * cgs%pi / planet%omega(3)
  d0 = atmosphere%external%density
  p0 = atmosphere%external%pressure
  T0 = atmosphere%external%temperature
  px0 = 0.0
  pz0 = 0.0
  !
  xbndy_lo = m1%p + m1%r(m1%li) - 0.5 * m1%d
  xkill_lo = xbndy_lo + self%x_sponge
  xbndy_hi = m1%p + m1%r(m1%ui) + 0.5 * m1%d
  xkill_hi = xbndy_hi - self%x_sponge
  !
  do iv=1,self%nv
    do k=m3%lb,m3%ub
      z = m3%p + m3%r(k)
      do j=m2%lb,m2%ub
        do i=m1%lb,m1%ub
          x = m1%p + m1%r(i)
          d = d0 * exp(-(planet%omega(3)**2 * 0.5 * (z**2 + 2.0 * &
            planet%trap * (x - planet%position(1))**2) / T0))
          p = p0*(d/d0)
          if (iv == self%idx%d) then
            asymptote = d
          else if (iv == self%idx%py) then
            asymptote = -d * ((1.5+planet%trap)*planet%omega(3)*x - planet%v_head)
          else if (iv == self%idx%e) then
            if (self%etype == 'total') then
              py0 = -d * ((1.5+planet%trap)*planet%omega(3)*x - planet%v_head)
              asymptote = p / (eos%gamma - 1.0) + 0.5 * d * py0**2
            else if (self%etype == 'thermal') then
              asymptote = p / (eos%gamma - 1.0)
            else if (self%etype == 'entropy') then
              asymptote = d * (log(p) - eos%gamma * log(d)) / (eos%gamma - 1.0)
            else
              call io%abort ('unknown %etype in apply_sponge')
            end if
          else if (iv == self%idx%px) then
            asymptote = px0
          else if (iv == self%idx%pz) then
            asymptote = pz0
          end if
          !
          ! apply damping procedure to each variable.
          if (self%boundaries%is_set(bits%xl) .and. (x <= xkill_lo)) then
            mem(i,j,k,iv) = mem(i,j,k,iv) + self%dtime* &
              self%boundaries%sponge(mem(i,j,k,iv), asymptote, tau_kill, &
              xkill_lo, xbndy_lo, x)
          else if (self%boundaries%is_set(bits%xu) .and. (x >= xkill_hi)) then
            mem(i,j,k,iv) = mem(i,j,k,iv) + self%dtime* &
              self%boundaries%sponge(mem(i,j,k,iv), asymptote, tau_kill, &
              xkill_hi, xbndy_hi, x)
          else if (self%boundaries%is_set(bits%xl) .and. (i <= m1%lo)) then
            ! reset ghost cell values; in the case of the x-velocity, make it
            ! anti-symmetric lower x-boundary
            if (iv == self%idx%px) then
              eye = 2*m1%lo - i + 1
              mem(i,j,k,iv) = -mem(eye,j,k,iv)
            else
              mem(i,j,k,iv) = asymptote
            end if
          else if (self%boundaries%is_set(bits%xu) .and. (i >= m1%uo)) then
            ! upper x-boundary
            if (iv == self%idx%px) then
              eye = 2*m1%uo - i - 1
              mem(i,j,k,iv) = -mem(eye,j,k,iv)
            else
              mem(i,j,k,iv) = asymptote
            end if
          end if
        end do
      end do
    end do
  end do
  call trace%end()
END SUBROUTINE apply_sponge_boundaries

!===============================================================================
!===============================================================================
FUNCTION is_relevant (self) RESULT (ok)
  class(experiment_t):: self
  logical:: ok
  real(8):: r(3), rmax
  integer:: ix, iy, iz
  !-----------------------------------------------------------------------------
  ok = .false.
  rmax = 0.0
  do iz=-1,1,2
    r(3) = self%position(3) -planet%position(3) + iz*self%size(3)*0.5_8
    do iy=-1,1,2
      r(2) = self%position(2) -planet%position(2) + iy*self%size(2)*0.5_8
      do ix=-1,1,2
        r(1) = self%position(1) -planet%position(1) + ix*self%size(1)*0.5_8
        rmax = max(rmax,sqrt(sum(r**2)))
        if (rmax > planet%radius) ok=.true.
      end do
    end do
  end do
  if (verbose > 1) &
    print '("is_relevant:",i6,1p,2e11.3,l3)', self%level, rmax, planet%radius, ok
END FUNCTION is_relevant

END MODULE experiment_mod
