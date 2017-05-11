!> =============================================================================
!> Setup the initial conditions and experiment update
!> =============================================================================
MODULE experiment_mod
  USE io_mod
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
  implicit none
  private
  type, public, extends(pebble_patch_t):: experiment_t
  contains
    procedure :: init
    procedure :: debug
    procedure :: update
    procedure :: output
    procedure :: boundary_conditions
    procedure :: edge_boundary_conditions
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
!>         initial_mod::condition               ! computes ICs
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
  self%mhd = .false.
  self%msplit = 3
  call scaling%init
  call self%initial%init (self%kind, type)
  call self%force%init (self%kind, type)
  self%nv = 5
  self%mhd = .false.
  !----------------------------------------------------------------------------
  ! Through the hierarchy of class(task_t)%init calls, the call to %init below
  ! sets the IC values, either from reading in snapshot values, or by calling
  ! the IC condition procedure self%initial%condition() from gpatch_t%init.
  !----------------------------------------------------------------------------
  self%periodic = .false.
  self%periodic(2) = .true.
  call self%pebble_patch_t%init
  ! In the case that we're using ZEUS, `self%nv` should now be 10.
  !self%mhd = .false.
  !self%nv = 5
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
  call self%debug ('end init')
  !----------------------------------------------------------------------------
  ! Export the patch ids for debugging purposes
  !----------------------------------------------------------------------------
  self%force%id = self%id
  self%initial%id = self%id
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE debug (self, label, it)
  class(experiment_t):: self
  integer, optional:: it
  character(len=*):: label
  integer:: it_l, n
  !.............................................................................
  if (present(it)) then
    it_l = it
  else
    it_l = self%it
  end if
  if (self%id==io%id_debug) then
    if (self%position(1) > 0.0) then
      n = self%mesh(1)%gn
      print '(i4,1x,a,f6.1,1p,3(1x,3e14.6),1x,a)', &
        self%id, 'DBG:', self%time, &
        self%mem(n-3:n-1,18,18,self%i_d ,it_l,1), &
        self%mem(n-3:n-1,18,18,self%i_s ,it_l,1)/self%mem(n-3:n-1,18,18,self%i_d ,it_l,1), &
        self%mem(n-3:n-1,18,18,self%i_px,it_l,1), trim(label)
    else
      print '(i4,1x,a,f6.1,1p,2(1x,3e14.6),1x,a)', &
        self%id, 'DBG:', self%time, &
        self%mem(2:4,18,18,self%i_d ,it_l,1), &
        self%mem(2:4,18,18,self%i_s ,it_l,1)/self%mem(2:4,18,18,self%i_d ,it_l,1), &
        self%mem(2:4,18,18,self%i_px,it_l,1), trim(label)
    end if
  end if
END SUBROUTINE debug

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
  call self%debug('start update')
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
    call self%debug ('after BCs')
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
    call self%debug ('after update', self%new)
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
  if (self%time > self%out_next) &
    call self%boundary_conditions (self%it)
  call self%pebble_patch_t%output
END SUBROUTINE

!> =============================================================================
!> Boundary conditions for this experiment
!> =============================================================================
SUBROUTINE boundary_conditions (self, it)
  class(experiment_t) :: self
  integer             :: it
  !.............................................................................
  real(8)             :: h(3)
  integer             :: iv, nb,i
  integer, save       :: itimer=0
  real, pointer       :: d(:,:,:), e(:,:,:), s(:,:,:)
  real                :: s0
  !-----------------------------------------------------------------------------
  if (self%boundaries%bits /= 0) then
    call trace%begin ('experiment_mod::boundary_conditions', itimer=itimer)
  !  print *, self%id, 'BOUNDARY CONDITION BITS:', self%boundaries%bits
    call self%check_density ('before boundary_conditions')
    !---------------------------------------------------------------------------
    ! If spherical inner BCs, save the incoming values
    !---------------------------------------------------------------------------
    if (self%boundaries%is_set (bits%spherical)) then
      if (io%verbose>1) print *, self%id, 'BOUNDARY CONDITION: spherical'
      print *, self%id, 'BOUNDARY CONDITION: spherical, t =', self%t(self%it)
      call self%boundaries%spherical_extrapolate (self%mem(:,:,:,self%i_d,it,1), self%mesh, self%i_d)
      d => self%mem(:,:,:,self%i_d,it,1)
      s0 = (log(atmosphere%surface%pressure) - &
        log(atmosphere%surface%density)*self%gamma)/(self%gamma-1.0)
      if (self%kind(1:9)=='stagger2_') then
        s => self%mem(:,:,:,self%i_s,it,1)
        s = s/d
        call self%boundaries%spherical (s, s0, self%mesh, self%i_s)
        s = s*d
      else
        allocate (s(self%gn(1),self%gn(2),self%gn(3)))
        call self%e2S (s, it)
        call self%boundaries%spherical (s, s0, self%mesh, self%i_s)
        call self%S2e (s, it)
        deallocate (s)
      end if
      if (self%nw>1) then
        call self%boundaries%spherical (self%mem(:,:,:,self%i_d,it,2), 0.0, self%mesh, self%i_d)
        call self%boundaries%spherical (self%mem(:,:,:,self%i_e,it,2), 0.0, self%mesh, self%i_e)
      end if
      do iv=self%i_px,self%i_pz
        call self%boundaries%spherical (self%mem(:,:,:,iv,it,1), 0.0, self%mesh, iv)
        if (self%nw>1) &
          call self%boundaries%spherical (self%mem(:,:,:,iv,it,2), 0.0, self%mesh, iv)
      end do
    end if
    !------------------------------------------------------------------------
    ! Edge boundary conditions
    !------------------------------------------------------------------------
    if (self%boundaries%bits /= 0) then
      if (io%verbose>2) then
        print *, self%id, 'MK x-BC:', self%boundaries%is_set (bits%xl), self%boundaries%is_set (bits%xu)
        print *, self%id, 'MK y-BC:', self%boundaries%is_set (bits%yl), self%boundaries%is_set (bits%yu)
        print *, self%id, 'MK z-BC:', self%boundaries%is_set (bits%zl), self%boundaries%is_set (bits%zu)
      end if
      call self%edge_boundary_conditions(it)
    end if
    call self%check_density ('after boundary_conditions')
    call trace%end
  end if
END SUBROUTINE boundary_conditions

!> =============================================================================
!> Boundary conditions at the edges of the outermost rubiks cubes; extrapolate
!> the density and pressure, and enforce vanishing derivative on the perpendicular
!> velocities. To what extend the boundary conditions are applied in the x,y,z
!> directions depends on the self%boundaries%bits bits. Note that the same
!> condition is applied in all directions that have boundary bits set!
!> =============================================================================
SUBROUTINE edge_boundary_conditions (self, it)
  class(experiment_t):: self
  integer:: it
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e_th, uu, s, ss
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: u
  integer:: iv, select
  real:: ee
  !-----------------------------------------------------------------------------
  ! Edge boundary conditions
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%edge_boundary_conditions')
  self%boundaries%id = self%id
  d => self%mem(:,:,:,self%i_d,it,1)
  allocate (u   (size(d,1),size(d,2),size(d,3),3))
  allocate (e_th(size(d,1),size(d,2),size(d,3)  ))
                            call self%check_density ('before density_BC')
  call self%density_boundary_condition (it)                             ! vertical density BC
  !call self%boundaries%log_extrapolate (d, self%i_d)
  if (self%position(1)<0.0 .and. self%position(2)>0.0 .and. self%boundaries%is_set(bits%yu)) then
    select = bits%yu
    if (io%verbose>2) &
      print *, self%id, 'EXTRAP rho at yu', &
      self%boundaries%is_set(bits%xl), self%boundaries%is_set(bits%xu), &
      self%boundaries%is_set(bits%yl), self%boundaries%is_set(bits%yu), &
      self%boundaries%is_set(bits%zl), self%boundaries%is_set(bits%zu)
    call self%boundaries%log_extrapolate (d, self%i_d, select=select)
    if (self%id==io%id_debug) &
      print *, self%id, 'EXTRAP D', d(6,33:37,31)
  else if (self%position(1)>0.0 .and. self%position(2)<0.0 .and. self%boundaries%is_set(bits%yl)) then
    select = bits%yl
    if (io%verbose>2) &
      print *, self%id, 'EXTRAP rho at yl', &
      self%boundaries%is_set(bits%xl), self%boundaries%is_set(bits%xu), &
      self%boundaries%is_set(bits%yl), self%boundaries%is_set(bits%yu), &
      self%boundaries%is_set(bits%zl), self%boundaries%is_set(bits%zu)
    call self%boundaries%log_extrapolate (d, self%i_d, select=select)
    if (self%id==io%id_debug) &
      print *, self%id, 'EXTRAP D', d(6,1:4,31)
  end if
                            call self%check_density ('after density_BC')
  !---------------------------------------------------------------------
  ! In the vertical and radial directions, where density varies
  ! systematically, apply BCs to velocity
  !---------------------------------------------------------------------
  call self%p2U(u,it)                                                      ! velocities
  do iv=1,3
    uu => u(:,:,:,iv)
!    call self%boundaries%symmetric (uu, self%i_px+iv-1, select=select) ! apply symmetric boundaries to velocities
    !-------------------------------------------------------------------
    ! Extrapolate at boundaries with strong supersonic outflows
    !-------------------------------------------------------------------
    if      (self%position(1)<0.0 .and. self%position(2)>0.0 .and. self%boundaries%is_set(bits%yu)) then
      select = iand(self%boundaries%bits,not(bits%yu))
      call self%boundaries%symmetric (uu, self%i_px+iv-1, select=select)
      select = bits%yu
      call self%boundaries%extrapolate (uu, self%i_px+iv-1, select=select)
      if (io%verbose>2) &
        print *, self%id, 'EXTRAP U1', iv, u(6,33:37,31,2)
    else if (self%position(1)>0.0 .and. self%position(2)<0.0 .and. self%boundaries%is_set(bits%yl)) then
      select = iand(self%boundaries%bits,not(bits%yl))
      call self%boundaries%symmetric (uu, self%i_px+iv-1, select=select)
      select = bits%yl
      call self%boundaries%extrapolate (uu, self%i_px+iv-1, select=select)
      if (io%verbose>2) &
        print *, self%id, 'EXTRAP U2', iv, u(6,33:37,31,2)
    !-------------------------------------------------------------------
    ! The azimuthal velocity should be safe to extrapolate in directions
    ! where it it paralle to the boundary
    !-------------------------------------------------------------------
    !-------------------------------------------------------------------
    else if (iv==2) then
      call self%boundaries%extrapolate (uu, self%i_px+iv-1)
    !-------------------------------------------------------------------
    ! Use vanishing perpendicular derivative in other cases
    !-------------------------------------------------------------------
    else
      call self%boundaries%symmetric (uu, self%i_px+iv-1)
    end if
  end do
                            call self%check_density ('after u BC')
  call self%U2p(u,it)                                                      ! convert velocities to momentum
  !---------------------------------------------------------------------
  ! Extrapolate in per unit mass variables
  !---------------------------------------------------------------------
  if (self%kind(1:9) == 'stagger2_') then
    s => self%mem(:,:,:,self%i_s,it,1)
    ss => e_th                                                          ! borrow local allocation
    ss = s/d
    call self%boundaries%extrapolate(ss, self%i_s)                      ! extrapolate in entropy per unit mass
    s = ss*d
  else
    call self%e2E_th(e_th, it)                                          ! thermal energy per unit mass
    call self%boundaries%log_extrapolate(e_th, self%i_e)                ! extrapolate in the log
    call self%E_th2e(e_th, it)                                          ! convert to energy per unit volume
  end if
                            call self%check_nan ('after th2e')
  deallocate (u,e_th)
  call trace%end
END SUBROUTINE edge_boundary_conditions

!> =====================================================================
!> Constrain the density to follow an anchored Gaussian profile
!> =====================================================================
SUBROUTINE density_boundary_condition (self, it)
  class(experiment_t):: self
  integer:: it
  class(mesh_t), pointer, dimension(:) :: m
  integer:: iv, iz, iy, ix, lz, uz, ly, uy, lx, ux, ix0, iy0, iz0
  real(8) :: x, x0, z, z0, r, d1, d0, tt, c2
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, dd
  !-----------------------------------------------------------------------------
  ! Vertical edge boundary conditions
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%density_boundary_conditions')
  m => self%mesh
  d => self%mem(:,:,:,1,it,1)
  if (self%nw>1) dd => self%mem(:,:,:,self%i_d,it,2)
  associate (r1=>m(1)%r, r3=>m(3)%r)
  if (self%id==io%id_debug) &
    print '(a,2i5,2l4,3(2x,2i3))', 'MK', self%id, it, &
      self%boundaries%is_set(bits%zl), &
      self%boundaries%is_set(bits%zu), lx,ux, ly,uy, lz,uz
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
  call self%check_density ('after density BC')
  end associate
  call trace%end
  return
contains
  subroutine apply_x
    if (io%verbose>2) &
    print *, self%id, 'MK: x-density BC', lx, ux, ly, uy, lz, uz
    c2 = planet%omega(3)**2*planet%trap/tt
    associate (r1=>m(1)%r)
    do iz=lz,uz
      do iy=ly,uy
        do ix=lx,ux
          x = m(1)%p + r1(ix) - planet%position(1)
          d0 = self%mem(ix0,iy,iz,self%i_d,it,1)
          self%mem(ix,iy,iz,self%i_d,it,1) = d0*exp(-c2*(x**2-x0**2))
        end do
      end do
    end do
    end associate
    if (self%nw>1) then
      do iz=lz,uz
        do iy=ly,uy
          do ix=lx,ux
            self%mem(ix,iy,iz,self%i_d,it,2) = 0.0
          end do
        end do
      end do
    end if
  end subroutine apply_x
  subroutine apply_y
    if (io%verbose>2) &
      print *, self%id, 'MK: y-density BC', lx, ux, ly, uy, lz, uz
    do iz=lz,uz
      do iy=ly,uy
        do ix=lx,ux
          self%mem(ix,iy,iz,self%i_d,it,1) = self%mem(ix,iy0,iz,self%i_d,it,1)
        end do
      end do
    end do
    if (self%nw>1) then
      do iz=lz,uz
        do iy=ly,uy
          do ix=lx,ux
            self%mem(ix,iy,iz,self%i_d,it,2) = 0.0
          end do
        end do
      end do
    end if
  end subroutine apply_y
  subroutine apply_z
    if (io%verbose>2) &
    print *, self%id, 'MK: z-density BC', lx, ux, ly, uy, lz, uz
    c2 = 0.5*planet%omega(3)**2/tt
    associate (r3=>m(3)%r)
    do iz=lz,uz
      z = m(3)%p + r3(iz) - planet%position(3)
      do iy=ly,uy
        do ix=lx,ux
          d0 = self%mem(ix,iy,iz0,self%i_d,it,1)
          d1 = d0*exp(-c2*(z**2-z0**2))
          if (self%id==io%id_debug.and.iy==uy/2.and.iz==uz) &
            print '(a,i4,1p,8e14.6)', &
            'BCd: ix,d0,p0/do,exp,do,d =',ix, d0, tt, (z-z0), &
            exp(-c2*(z**2-z0**2)), &
            self%mem(ix,iy,iz,self%i_d,it,1), d1
          self%mem(ix,iy,iz,self%i_d,it,1) = d1
        end do
      end do
    end do
    end associate
    if (self%nw>1) then
      do iz=lz,uz
        do iy=ly,uy
          do ix=lx,ux
            self%mem(ix,iy,iz,self%i_d,it,2) = 0.0
          end do
        end do
     end do
    end if
  end subroutine apply_z
END SUBROUTINE density_boundary_condition

END MODULE
